/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "../extern/pdqsort/pdqsort.h"
#include "mip/HighsMipSolverData.h"
#include "mip/MipTimer.h"
#include "presolve/HPresolve.h"
#include "presolve/HPresolveUtils.h"

namespace presolve {

HPresolve::Result HPresolve::implAwareConstrPropagation(
    HighsPostsolveStack& postsolve_stack) {
  if (mipsolver == nullptr) return Result::kOk;

  HighsProfilingScope clock(*mipsolver->profiling_,
                            kMipClockImplAwareConstrPropagation);

  HighsCliqueTable& cliquetable = mipsolver->mipdata_->cliquetable;
  HighsImplications& implications = mipsolver->mipdata_->implications;

  struct VariableData {
    double val;
    HighsCDouble weightLower;
    HighsCDouble weightUpper;
    // position of a row binary in blockMembers, -1 if it is not in a clique
    // block with other row binaries
    HighsInt member;

    VariableData(double myVal = 0.0)
        : val(myVal), weightLower(0.0), weightUpper(0.0), member(-1) {}

    void update(bool upper, const HighsCDouble& delta) {
      if (upper)
        weightUpper += delta;
      else
        weightLower += delta;
    }
  };

  // map from column indices to variable data with direct access; insertion,
  // lookup and clearing in the number of entries. the entries are iterated
  // in the order of insertion
  struct IndexMap {
    struct Entry {
      HighsInt index;
      VariableData data;
      HighsInt key() const { return index; }
      VariableData& value() { return data; }
      const VariableData& value() const { return data; }
    };
    std::vector<Entry> entries;
    // position of an index in entries, -1 if it is not in the map
    std::vector<HighsInt> position;

    // prepares the map for the indices [0, size)
    void resize(HighsInt size) { position.resize(size, -1); }

    VariableData* find(HighsInt index) {
      return position[index] < 0 ? nullptr : &entries[position[index]].data;
    }

    VariableData& operator[](HighsInt index) {
      if (position[index] < 0) {
        position[index] = static_cast<HighsInt>(entries.size());
        entries.push_back({index, VariableData()});
      }
      return entries[position[index]].data;
    }

    void clear() {
      for (const Entry& entry : entries) position[entry.index] = -1;
      entries.clear();
    }

    size_t size() const { return entries.size(); }
    std::vector<Entry>::iterator begin() { return entries.begin(); }
    std::vector<Entry>::iterator end() { return entries.end(); }
  };

  struct Breakpoint {
    double value;
    // row binary forced away from its min-contribution value
    HighsInt binCol;
    bool isLowerType;
  };

  struct ImpliedBoundKey {
    HighsInt binCol;
    HighsInt binVal;
    HighsInt nonBinCol;
    bool operator==(const ImpliedBoundKey& other) const {
      return binCol == other.binCol && binVal == other.binVal &&
             nonBinCol == other.nonBinCol;
    }
  };

  // positive lifts of binaries due to implied bounds on row non-binaries
  struct Lift {
    HighsInt binCol;
    HighsInt binVal;
    HighsCDouble value;
  };

  // binary in the row with its sort key for the clique propagation
  struct SortedBin {
    HighsInt col;
    double sortKey;
  };

  // member of a clique block: row binaries whose min-contribution literals
  // form a clique, so at most one of them takes its min-contribution value
  struct BlockMember {
    HighsInt col;
    HighsInt block;
    // number of reasons for the member to leave its min-contribution value
    HighsInt count;
  };

  // literal considered for two-column clique extraction with its weight
  struct CliqueCandidate {
    HighsCliqueTable::CliqueVar var;
    HighsCDouble weight;
  };

  // set of indices with insertion and clearing in the number of entries
  struct IndexSet {
    std::vector<HighsInt> entries;
    std::vector<HighsBool> flags;

    // prepares the set for the indices [0, size)
    void resize(HighsInt size) {
      flags.resize(size);
      entries.reserve(size);
    }

    // adds an index unless it is already in the set
    bool add(HighsInt index) {
      if (flags[index]) return false;
      flags[index] = true;
      entries.push_back(index);
      return true;
    }

    void clear() {
      for (HighsInt index : entries) flags[index] = false;
      entries.clear();
    }

    // moves the entries to target and clears the set
    void moveEntriesTo(std::vector<HighsInt>& target) {
      for (HighsInt index : entries) flags[index] = false;
      std::swap(entries, target);
      entries.clear();
    }
  };

  // bound on a non-binary implied by a binary
  struct CachedBound {
    HighsInt binCol;
    HighsInt binVal;
    double lb;
    double ub;
  };

  // bounds on a non-binary implied by binaries. they do not depend on the
  // row, and the variable bounds and implications do not change during the
  // call, so they are collected once
  struct CachedColumn {
    std::vector<CachedBound> bounds;
    // variable bounds and reverse implications visited to collect the
    // bounds, -1 if not collected yet
    HighsInt numVisits = -1;
  };

  // data structures
  IndexMap binNonZeros;
  IndexMap nonBinNonZeros;
  std::vector<SortedBin> sortedBins;
  std::vector<CliqueCandidate> candidates;
  std::vector<Breakpoint> breakpoints;
  // bounds implied by binaries on the non-binaries of the current row, as
  // read by the propagation (filtered and merged), keyed by literal and
  // non-binary
  HighsHashTable<ImpliedBoundKey, HighsImplications::Implication> impliedBounds;
  // implied bounds on non-binaries, indexed by column
  std::vector<CachedColumn> cachedColumns;
  // position of a literal (2 * col + val) in the cached bounds of the column
  // being collected, -1 if it is not there
  std::vector<HighsInt> literalPos;
  std::vector<Lift> lifts;
  std::vector<HighsInt> liftedBins;
  // clique blocks of the row binaries: the members of block b are
  // blockMembers[blockStarts[b], blockStarts[b + 1]), sorted by decreasing
  // absolute coefficient
  std::vector<HighsInt> blockStarts;
  std::vector<BlockMember> blockMembers;
  IndexSet touchedBlocks;
  // literals and the block members they force away from their
  // min-contribution value
  std::vector<std::pair<HighsCliqueTable::CliqueVar, HighsInt>> blockPairs;
  std::vector<HighsCliqueTable::CliqueVar> partitionVars;
  std::vector<HighsInt> partitionStarts;
  std::vector<HighsCliqueTable::CliqueVar> impliedFixings;
  std::vector<
      std::pair<HighsCliqueTable::CliqueVar, HighsCliqueTable::CliqueVar>>
      conflictingLiterals;
  std::vector<HighsCliqueTable::CliqueVar> infeasibleLiterals;
  IndexSet modifiedRows;
  IndexSet changedNonBinaries;
  std::vector<HighsInt> objInds;
  bool checkObjective = false;
  bool objectiveAffected = false;

  // counters
  HighsInt numVarsFixed = 0;
  HighsInt numCliquesAdded = 0;
  HighsInt numBoundsTightened = 0;
  // calls of checkRow, for the rows and the objective (numObjectiveChecks of
  // them for the objective)
  HighsInt numCheckRowCalls = 0;
  HighsInt numObjectiveChecks = 0;
  // work in elementary operations: clique queries, look-ups of implied bounds
  // and visits of implications, variable bounds and clique neighbours
  int64_t numWork = 0;

  // maximum number of non-zeros in rows or objective
  const HighsInt maxRowSize = 300;
  const HighsInt maxObjSize = 20000;
  // maximum work per row to discover non-binaries: implications visited and
  // the visits to collect the bounds of discovered non-binaries
  const HighsInt maxDiscoveryWork = 10 * maxRowSize;
  // maximum number of cliques added
  const HighsInt maxNumCliquesAdded =
      std::max(HighsInt{1000}, std::min(HighsInt{10000}, numNonzeros() / 10));
  // maximum number of clique queries to partition a row into clique blocks
  const int64_t maxPartitionQueries = int64_t{maxRowSize} * maxRowSize / 2;
  // check for timeout according to this frequency (number of checkRow calls)
  const HighsInt checkTimeFrequency = 100;
  // work limit per call: as in probing, it scales with the number of
  // non-zeros, is smaller in sub-MIPs and grows with the reductions found in
  // the main presolve
  const int64_t workPerNonzero = mipsolver->submip ? 10 : 30;
  const int64_t minWork = mipsolver->submip ? 0 : 100000;
  const int64_t workBonusPerBoundChange = 100;
  const int64_t workBonusPerClique = 1000;
  const int64_t workLimit = minWork + workPerNonzero * numNonzeros();

  auto isBinary = [&](HighsInt col) {
    // implicit integers are binary as well, consistent with the clique table
    // and the implications
    return model->integrality_[col] != HighsVarType::kContinuous &&
           model->col_lower_[col] == 0.0 && model->col_upper_[col] == 1.0;
  };

  // value of a binary with row coefficient val that minimises its
  // contribution to the row
  auto minContributionVal = [](double val) {
    return val < 0 ? HighsInt{1} : HighsInt{0};
  };

  auto trackColChange = [&](HighsInt col) {
    for (const auto& nz : getColumnVector(col))
      if (rowsize[nz.index()] > 1 && rowsize[nz.index()] <= maxRowSize)
        modifiedRows.add(nz.index());
    if (checkObjective && model->col_cost_[col] != 0.0)
      objectiveAffected = true;
    // non-binaries with changed bounds are checked for conflicting
    // implications in the next round
    if (!isBinary(col)) changedNonBinaries.add(col);
  };

  auto fixCol = [&](HighsInt col, bool toUpper) {
    numVarsFixed++;
    trackColChange(col);
    HPRESOLVE_CHECKED_CALL(toUpper ? fixColToUpper(postsolve_stack, col)
                                   : fixColToLower(postsolve_stack, col));
    return checkLimits(postsolve_stack);
  };

  auto changeNonBinaryLower = [&](HighsInt col, double newLowerBnd) {
    if (colDeleted[col] ||
        newLowerBnd <= model->col_lower_[col] + primal_feastol)
      return Result::kOk;
    bool accept = true;
    if (model->integrality_[col] != HighsVarType::kContinuous)
      newLowerBnd = std::ceil(newLowerBnd - primal_feastol);
    else
      // continuous columns: same acceptance criteria as domain propagation
      newLowerBnd = mipsolver->adjustedLb(
          static_cast<HighsCDouble>(newLowerBnd), model->col_lower_[col],
          model->col_upper_[col], false, accept);
    assert(newLowerBnd <= model->col_upper_[col] + primal_feastol);
    if (newLowerBnd == model->col_upper_[col]) {
      HPRESOLVE_CHECKED_CALL(fixCol(col, true));
    } else if (accept) {
      numBoundsTightened++;
      trackColChange(col);
      HPRESOLVE_CHECKED_CALL(changeColLower(col, newLowerBnd));
    }
    return Result::kOk;
  };

  auto changeNonBinaryUpper = [&](HighsInt col, double newUpperBnd) {
    if (colDeleted[col] ||
        newUpperBnd >= model->col_upper_[col] - primal_feastol)
      return Result::kOk;
    bool accept = true;
    if (model->integrality_[col] != HighsVarType::kContinuous)
      newUpperBnd = std::floor(newUpperBnd + primal_feastol);
    else
      // continuous columns: same acceptance criteria as domain propagation
      newUpperBnd = mipsolver->adjustedUb(
          static_cast<HighsCDouble>(newUpperBnd), model->col_lower_[col],
          model->col_upper_[col], false, accept);
    assert(newUpperBnd >= model->col_lower_[col] - primal_feastol);
    if (newUpperBnd == model->col_lower_[col]) {
      HPRESOLVE_CHECKED_CALL(fixCol(col, false));
    } else if (accept) {
      numBoundsTightened++;
      trackColChange(col);
      HPRESOLVE_CHECKED_CALL(changeColUpper(col, newUpperBnd));
    }
    return Result::kOk;
  };

  auto addClique = [&](HighsCliqueTable::CliqueVar v1,
                       HighsCliqueTable::CliqueVar v2) {
    std::array<HighsCliqueTable::CliqueVar, 2> clique = {v1, v2};
    StatusResult added = updateCliqueTableAddClique(
        clique.data(), static_cast<HighsInt>(clique.size()), impliedFixings);
    HPRESOLVE_CHECKED_CALL(static_cast<Result>(added));
    if (!added) return Result::kOk;
    numCliquesAdded++;
    // remember affected columns
    trackColChange(v1.col);
    trackColChange(v2.col);
    // apply fixings derived by the clique table
    for (const HighsCliqueTable::CliqueVar& fixing : impliedFixings)
      HPRESOLVE_CHECKED_CALL(fixCol(fixing.col, fixing.val == 1));
    return Result::kOk;
  };

  auto addNonZero = [&](HighsInt col, double val) {
    // skip fixed columns that have not been removed yet; their contribution is
    // already accounted for in the threshold
    if (model->col_lower_[col] == model->col_upper_[col]) return true;
    if (isBinary(col)) {
      VariableData data(val);
      data.update(val > 0, std::abs(val));
      binNonZeros[col] = data;
    } else {
      if ((val > 0 && model->col_lower_[col] <= -kHighsInf) ||
          (val < 0 && model->col_upper_[col] >= kHighsInf))
        return false;
      VariableData data(val);
      if (model->col_lower_[col] > -kHighsInf &&
          model->col_upper_[col] < kHighsInf) {
        HighsCDouble range =
            std::abs(val) * (static_cast<HighsCDouble>(model->col_upper_[col]) -
                             model->col_lower_[col]);
        data.update(val > 0, range);
      }
      nonBinNonZeros[col] = data;
    }
    return true;
  };

  auto loadModelRow = [&](HighsInt row, HighsInt direction) {
    binNonZeros.clear();
    nonBinNonZeros.clear();

    for (const auto& nz : getRowVector(row))
      // add non-zero
      if (!addNonZero(nz.index(), direction * nz.value())) return false;
    return true;
  };

  auto loadObjective = [&](HighsInt direction, HighsCDouble& minActivity) {
    binNonZeros.clear();
    nonBinNonZeros.clear();
    minActivity = 0;

    for (HighsInt col : objInds) {
      // costs of active columns do not change in this method
      if (colDeleted[col]) continue;
      double cost = direction * model->col_cost_[col];
      // add non-zero; rejects infinite bounds on the min-activity side
      if (!addNonZero(col, cost)) return false;
      // accumulate standard minimum
      double bound = cost > 0 ? model->col_lower_[col] : model->col_upper_[col];
      assert(std::abs(bound) < kHighsInf);
      minActivity += cost * static_cast<HighsCDouble>(bound);
    }
    return true;
  };

  auto addImpliedBound = [&](HighsInt binCol, HighsInt binVal,
                             HighsInt nonBinCol, double lb, double ub) {
    // only the bounds that are read are stored: both sides for row binaries
    // at their min-contribution value (looked up by collectBreakpoints),
    // otherwise only the side that lifts the non-binary (lower bound for a
    // positive coefficient, upper bound for a negative one). discovered
    // non-binaries (zero coefficient) do not lift. the binaries in
    // binNonZeros are the row binaries at this point, since lifted binaries
    // and clique neighbours are only added in updateBinaryWeights
    const VariableData* bin = binNonZeros.find(binCol);
    bool isBreakpointKey = bin != nullptr && bin->val != 0.0 &&
                           binVal == minContributionVal(bin->val);
    if (!isBreakpointKey) {
      // implied bounds are only collected for loaded non-binaries
      const VariableData* nonBinVar = nonBinNonZeros.find(nonBinCol);
      assert(nonBinVar != nullptr);
      double val = nonBinVar->val;
      if (val <= 0) lb = -kHighsInf;
      if (val >= 0) ub = kHighsInf;
      if (lb == -kHighsInf && ub == kHighsInf) return;
    }
    // keep the tightest bounds. each key is currently inserted at most once
    // per row (the cached bounds are merged per literal), so this only guards
    // against future sources of implied bounds
    HighsImplications::Implication& bound =
        impliedBounds[ImpliedBoundKey{binCol, binVal, nonBinCol}];
    bound.lb = std::max(bound.lb, lb);
    bound.ub = std::min(bound.ub, ub);
  };

  // value of a variable bound coef * x_bin + constant at x_bin = binVal
  auto varBoundValue = [](const HighsImplications::VarBound& vb,
                          HighsInt binVal) {
    return binVal == 1 ? vb.coef + vb.constant : vb.constant;
  };

  auto cacheImpliedBounds = [&](HighsInt col, CachedColumn& cached) {
    // bounds on a non-binary implied by variable bound constraints and
    // implications, collected unless they are cached already
    if (cached.numVisits >= 0) return;
    cached.numVisits = 0;
    auto addCachedBound = [&](HighsInt binCol, HighsInt binVal, double lb,
                              double ub) {
      // round bounds of integer non-binaries; variable bounds may imply
      // fractional values
      if (model->integrality_[col] != HighsVarType::kContinuous) {
        lb = std::ceil(lb - primal_feastol);
        ub = std::floor(ub + primal_feastol);
      }
      // bounds implied by the same literal are merged
      HighsInt& pos = literalPos[2 * binCol + binVal];
      if (pos < 0) {
        pos = static_cast<HighsInt>(cached.bounds.size());
        cached.bounds.push_back({binCol, binVal, lb, ub});
      } else {
        cached.bounds[pos].lb = std::max(cached.bounds[pos].lb, lb);
        cached.bounds[pos].ub = std::min(cached.bounds[pos].ub, ub);
      }
    };
    implications.getVlbs(col).for_each(
        [&](HighsInt binCol, const HighsImplications::VarBound& vlb) {
          cached.numVisits++;
          if (colDeleted[binCol] || !isBinary(binCol)) return;
          // x_bin = val --> x_j >= coef * val + constant
          for (HighsInt val = 1; val >= 0; --val)
            addCachedBound(binCol, val, varBoundValue(vlb, val), kHighsInf);
        });
    implications.getVubs(col).for_each(
        [&](HighsInt binCol, const HighsImplications::VarBound& vub) {
          cached.numVisits++;
          if (colDeleted[binCol] || !isBinary(binCol)) return;
          // x_bin = val --> x_j <= coef * val + constant
          for (HighsInt val = 1; val >= 0; --val)
            addCachedBound(binCol, val, -kHighsInf, varBoundValue(vub, val));
        });
    // implications of all binaries on the non-binary. the reverse
    // implications are never cleaned up, so look up each one
    implications.getReverseImplications(col).for_each([&](HighsInt binCol,
                                                          bool) {
      cached.numVisits++;
      if (colDeleted[binCol] || !isBinary(binCol)) return;
      for (HighsInt val = 0; val <= 1; ++val) {
        const auto* impl = implications.getImplications(binCol, val).find(col);
        if (impl != nullptr) addCachedBound(binCol, val, impl->lb, impl->ub);
      }
    });
    // reset the positions for the next column
    for (const CachedBound& bound : cached.bounds)
      literalPos[2 * bound.binCol + bound.binVal] = -1;
  };

  auto addBoundsFromCache = [&](HighsInt col) {
    for (const CachedBound& bound : cachedColumns[col].bounds) {
      // binaries may have been fixed since the bounds were cached
      if (colDeleted[bound.binCol] || !isBinary(bound.binCol)) continue;
      addImpliedBound(bound.binCol, bound.binVal, col, bound.lb, bound.ub);
    }
  };

  auto addDiscoveredBounds = [&](HighsInt col, HighsInt& numDiscoveryWork) {
    // bounds on a discovered non-binary implied by binaries. as for row
    // non-binaries, the variable bounds and reverse implications visited to
    // collect them count as work; here, they count towards the discovery
    // budget, which is checked before each visit in discoverNonBinaries and
    // added to the work at the end of collectImpliedBounds
    CachedColumn& cached = cachedColumns[col];
    cacheImpliedBounds(col, cached);
    numDiscoveryWork += cached.numVisits;
    addBoundsFromCache(col);
  };

  auto discoverNonBinaries = [&](HighsInt binCol, HighsInt& numDiscoveryWork) {
    // non-binaries outside the row implied by a binary are discovered; their
    // bounds can be tightened via the piecewise walk even though their row
    // coefficient is zero. non-binaries linked to a row binary only by a
    // variable bound are not discovered (no reverse index; no reductions on
    // MIPLIB3 when tested). returns false when the row is full or the
    // discovery budget is used up, i.e. when discovery should stop
    for (HighsInt val = 0; val <= 1; ++val) {
      if (implications.getImplications(binCol, val)
              .for_each([&](HighsInt targetCol,
                            const HighsImplications::Implication&) {
                if (numDiscoveryWork >= maxDiscoveryWork ||
                    static_cast<HighsInt>(binNonZeros.size() +
                                          nonBinNonZeros.size()) >= maxRowSize)
                  return true;
                numDiscoveryWork++;
                if (colDeleted[targetCol] || isBinary(targetCol)) return false;
                // the bounds of known non-binaries are collected already;
                // they include this implication (reverse implications)
                if (nonBinNonZeros.find(targetCol) != nullptr) {
                  assert(std::any_of(cachedColumns[targetCol].bounds.begin(),
                                     cachedColumns[targetCol].bounds.end(),
                                     [&](const CachedBound& bound) {
                                       return bound.binCol == binCol &&
                                              bound.binVal == val;
                                     }));
                  return false;
                }
                if (model->col_lower_[targetCol] <= -kHighsInf ||
                    model->col_upper_[targetCol] >= kHighsInf)
                  return false;
                nonBinNonZeros[targetCol] = VariableData(0.0);
                addDiscoveredBounds(targetCol, numDiscoveryWork);
                return false;
              }))
        return false;
    }
    return true;
  };

  auto computeLift = [&](HighsInt nonBinCol,
                         const HighsImplications::Implication& impl) {
    // increase of a non-binary's minimum contribution to the row due to an
    // implied bound: a_j > 0 uses the implied lower bound, a_j < 0 the
    // implied upper bound
    const VariableData* nonBinVar = nonBinNonZeros.find(nonBinCol);
    // implied bounds are only collected for loaded non-binaries
    assert(nonBinVar != nullptr);
    double val = nonBinVar->val;
    // addNonZero rejects columns with an infinite bound on the
    // min-activity side
    assert(val <= 0 || model->col_lower_[nonBinCol] > -kHighsInf);
    assert(val >= 0 || model->col_upper_[nonBinCol] < kHighsInf);
    if (val > 0 && impl.lb > -kHighsInf)
      return val * (static_cast<HighsCDouble>(impl.lb) -
                    model->col_lower_[nonBinCol]);
    if (val < 0 && impl.ub < kHighsInf)
      return val * (static_cast<HighsCDouble>(impl.ub) -
                    model->col_upper_[nonBinCol]);
    return HighsCDouble(0.0);
  };

  auto collectImpliedBounds = [&]() {
    // collect the bounds on non-binaries implied by binaries from the cached
    // bounds of the row non-binaries and by discovering non-binaries
    impliedBounds.clear();
    // work to discover non-binaries, limited by maxDiscoveryWork
    HighsInt numDiscoveryWork = 0;
    for (const auto& nz : nonBinNonZeros) {
      HighsInt nonBinCol = nz.key();
      CachedColumn& cached = cachedColumns[nonBinCol];
      cacheImpliedBounds(nonBinCol, cached);
      // the visits count as work, also when the bounds are cached
      numWork += cached.numVisits;
      addBoundsFromCache(nonBinCol);
    }
    // discover non-binaries outside the row
    for (const auto& binVar : binNonZeros)
      if (!discoverNonBinaries(binVar.key(), numDiscoveryWork)) break;

    // compute the lifts and collect the binaries outside the row whose implied
    // bounds lift them. all implied bounds on row non-binaries have been
    // collected, so this is complete. the implied bounds added by discovering
    // non-binaries below are on columns outside the row and do not lift
    lifts.clear();
    liftedBins.clear();
    for (const auto& entry : impliedBounds) {
      const ImpliedBoundKey& key = entry.key();
      HighsCDouble lift = computeLift(key.nonBinCol, entry.value());
      if (lift <= 0) continue;
      lifts.push_back({key.binCol, key.binVal, lift});
      if (binNonZeros.find(key.binCol) == nullptr)
        liftedBins.push_back(key.binCol);
    }
    // a binary may lift several non-binaries
    pdqsort(liftedBins.begin(), liftedBins.end());
    liftedBins.erase(std::unique(liftedBins.begin(), liftedBins.end()),
                     liftedBins.end());
    for (HighsInt binCol : liftedBins)
      if (!discoverNonBinaries(binCol, numDiscoveryWork)) break;
    numWork += numDiscoveryWork;
  };

  auto memberVal = [&](HighsInt member) {
    return binNonZeros.find(blockMembers[member].col)->val;
  };

  auto blockIncrease = [&](HighsInt block) {
    // increase of the row's minimum activity when the members with a positive
    // count leave their min-contribution value: the largest remaining member
    // takes over from the largest one, or none remains
    double maxAbsVal = std::abs(memberVal(blockStarts[block]));
    for (HighsInt k = blockStarts[block]; k < blockStarts[block + 1]; ++k)
      if (blockMembers[k].count == 0) return maxAbsVal - std::abs(memberVal(k));
    return maxAbsVal;
  };

  auto partitionRowBinaries = [&]() {
    // partition the min-contribution literals of the row binaries into
    // cliques. at most one literal of a block holds, so all but one member
    // leave their min-contribution value. returns the resulting increase of
    // the row's minimum activity
    blockStarts.assign(1, 0);
    blockMembers.clear();
    partitionVars.clear();
    for (const auto& bin : binNonZeros) {
      HighsInt col = bin.key();
      double val = bin.value().val;
      if (colDeleted[col] || val == 0.0) continue;
      partitionVars.emplace_back(col, minContributionVal(val));
    }

    // largest coefficients first, ties by column index, so the members of
    // each block are sorted by decreasing absolute coefficient
    auto byDecreasingAbsVal = [&](const HighsCliqueTable::CliqueVar& a,
                                  const HighsCliqueTable::CliqueVar& b) {
      return std::make_pair(-std::abs(binNonZeros.find(a.col)->val), a.col) <
             std::make_pair(-std::abs(binNonZeros.find(b.col)->val), b.col);
    };
    // the clique separation budgets its work with the number of neighbourhood
    // queries, so the queries of the partition are not counted there; they
    // count as work of this method instead
    int64_t oldNumQueries = cliquetable.numNeighbourhoodQueries;
    cliquetable.cliquePartition(partitionVars, partitionStarts,
                                byDecreasingAbsVal, maxPartitionQueries);
    numWork += cliquetable.numNeighbourhoodQueries - oldNumQueries;
    cliquetable.numNeighbourhoodQueries = oldNumQueries;

    // a single block of all binaries with equal coefficients and no
    // non-binaries, e.g. the >= side of a set partitioning row, yields the
    // same reductions without blocks, which is cheaper. with blocks, its
    // clique-aware threshold would be zero and skip the pair test
    HighsInt numVars = static_cast<HighsInt>(partitionVars.size());
    if (partitionStarts.size() == 2 && numVars >= 2 &&
        nonBinNonZeros.size() == 0 &&
        std::abs(binNonZeros.find(partitionVars[0].col)->val) ==
            std::abs(binNonZeros.find(partitionVars[numVars - 1].col)->val))
      return HighsCDouble(0.0);

    HighsCDouble increase = 0.0;
    for (size_t p = 0; p + 1 < partitionStarts.size(); ++p) {
      HighsInt blockStart = partitionStarts[p];
      HighsInt blockEnd = partitionStarts[p + 1];
      if (blockEnd - blockStart < 2) continue;

      HighsInt block = static_cast<HighsInt>(blockStarts.size()) - 1;
      for (HighsInt k = blockStart; k < blockEnd; ++k) {
        HighsInt col = partitionVars[k].col;
        VariableData& data = binNonZeros[col];
        data.member = static_cast<HighsInt>(blockMembers.size());
        blockMembers.push_back({col, block, 0});
        // all members but the largest one leave their min-contribution value
        if (k > blockStart) increase += std::abs(data.val);
      }
      blockStarts.push_back(static_cast<HighsInt>(blockMembers.size()));
    }
    touchedBlocks.resize(static_cast<HighsInt>(blockStarts.size()) - 1);
    return increase;
  };

  auto updateBinaryWeights = [&](double threshold) {
    // lift binaries: their implied bounds increase the minimum contribution
    // of non-binaries to the row
    for (const Lift& lift : lifts)
      binNonZeros[lift.binCol].update(lift.binVal == 1, lift.value);

    // consider cliques: sort binaries for early termination.
    // few clique connections first (cheap), large coefficients last.
    // remainingWeight tracks unprocessed clique propagation potential,
    // maxWeight tracks the largest binary weight seen so far.
    sortedBins.clear();
    HighsCDouble remainingWeight = 0;
    HighsCDouble maxWeight = 0;
    // the threshold may be zero or slightly negative (within the tolerance),
    // which would make the sort keys infinite or flip their order
    double sortThreshold = std::max(threshold, primal_feastol);
    for (const auto& binVar : binNonZeros) {
      HighsInt col = binVar.key();
      double val = binVar.value().val;
      const HighsCDouble& weightLower = binVar.value().weightLower;
      const HighsCDouble& weightUpper = binVar.value().weightUpper;
      maxWeight = max(maxWeight, max(weightUpper, weightLower));
      // binaries outside the row (zero coefficient) only receive weight; they
      // pass nothing on to their clique neighbours
      if (val == 0.0) continue;
      sortedBins.push_back({col, cliquetable.numCliques(col, val < 0) -
                                     10.0 * std::abs(val) / sortThreshold});
      remainingWeight += std::abs(val);
    }

    // equal sort keys are ordered by column index
    pdqsort(sortedBins.begin(), sortedBins.end(),
            [](const SortedBin& a, const SortedBin& b) {
              return std::make_pair(a.sortKey, a.col) <
                     std::make_pair(b.sortKey, b.col);
            });

    // the weights are first accumulated additively, i.e. |a| for every binary
    // forced away from its min-contribution value. for members of clique
    // blocks this overestimates the increase relative to the clique-aware
    // minimum activity, so it is corrected below. the additive weights are
    // upper bounds on the exact ones, so the early termination stays valid
    blockPairs.clear();
    for (const SortedBin& sortedBin : sortedBins) {
      HighsInt col = sortedBin.col;
      // copies, since the propagation below inserts into binNonZeros
      const VariableData* bin = binNonZeros.find(col);
      double val = bin->val;
      HighsInt member = bin->member;
      double absval = std::abs(val);

      // no binary's weight can exceed threshold after processing the
      // current and remaining binaries, so skip clique propagation
      if (maxWeight + remainingWeight <= threshold + primal_feastol) break;
      remainingWeight -= absval;

      // collect neighbors and update weights
      cliquetable.forEachUniqueNeighbor(
          HighsCliqueTable::CliqueVar(col, minContributionVal(val)),
          [&](HighsCliqueTable::CliqueVar neighbor) {
            numWork++;
            // substitutions in the clique table can leave cliques {x, ~x}
            // until the next rebuild, so the column can be its own
            // neighbour; skip it, otherwise its coefficient is counted twice
            if (colDeleted[neighbor.col] || neighbor.col == col) return;
#ifndef NDEBUG
            // binaries outside the row that are not lifted: their implied
            // bounds on row non-binaries were collected in
            // collectImpliedBounds and do not lift them
            if (binNonZeros.find(neighbor.col) == nullptr) {
              HighsInt neighborCol = neighbor.col;
              for (HighsInt val = 0; val <= 1; ++val) {
                implications.getImplications(neighborCol, val)
                    .for_each([&](HighsInt targetCol,
                                  const HighsImplications::Implication&) {
                      const auto* impl = impliedBounds.find(
                          ImpliedBoundKey{neighborCol, val, targetCol});
                      assert(impl == nullptr ||
                             computeLift(targetCol, *impl) <= 0);
                    });
              }
            }
#endif
            VariableData& data = binNonZeros[neighbor.col];
            data.update(neighbor.val == 1, absval);
            maxWeight = max(maxWeight, neighbor.val == 1 ? data.weightUpper
                                                         : data.weightLower);
            if (member >= 0) blockPairs.emplace_back(neighbor, member);
          });
    }

    if (blockMembers.empty()) return;

    // a block member set to the value opposite to its min-contribution value
    // leaves it as well; addNonZero has charged |a| for this
    for (HighsInt k = 0; k < static_cast<HighsInt>(blockMembers.size()); ++k)
      blockPairs.emplace_back(
          HighsCliqueTable::CliqueVar(blockMembers[k].col,
                                      1 - minContributionVal(memberVal(k))),
          k);

    // per literal and block, replace the additive contributions of the
    // members forced away by the exact increase of the block. the members of
    // a block are contiguous in blockMembers, so sorting by member groups the
    // pairs by block
    pdqsort(blockPairs.begin(), blockPairs.end(),
            [](const std::pair<HighsCliqueTable::CliqueVar, HighsInt>& a,
               const std::pair<HighsCliqueTable::CliqueVar, HighsInt>& b) {
              return std::make_pair(a.first.index(), a.second) <
                     std::make_pair(b.first.index(), b.second);
            });
    size_t groupStart = 0;
    while (groupStart < blockPairs.size()) {
      HighsCliqueTable::CliqueVar literal = blockPairs[groupStart].first;
      HighsInt block = blockMembers[blockPairs[groupStart].second].block;
      size_t groupEnd = groupStart;
      HighsCDouble additive = 0.0;
      while (groupEnd < blockPairs.size() &&
             blockPairs[groupEnd].first.index() == literal.index() &&
             blockMembers[blockPairs[groupEnd].second].block == block) {
        HighsInt member = blockPairs[groupEnd].second;
        blockMembers[member].count++;
        additive += std::abs(memberVal(member));
        groupEnd++;
      }
      binNonZeros[literal.col].update(
          literal.val == 1,
          static_cast<HighsCDouble>(blockIncrease(block)) - additive);
      for (size_t k = groupStart; k < groupEnd; ++k)
        blockMembers[blockPairs[k].second].count--;
      groupStart = groupEnd;
    }
  };

  auto findCliques = [&](double threshold) {
    // two-column clique extraction
    if (threshold <= primal_feastol) return Result::kOk;

    // no more cliques can be added
    if (numCliquesAdded >= maxNumCliquesAdded) return Result::kOk;

    // stop if no cliques were found after a number of tries
    HighsInt consecutiveNoClique = 0;
    const HighsInt maxConsecutiveNoClique = 10;

    // collect candidates
    candidates.clear();
    for (const auto& bin : binNonZeros) {
      HighsInt col = bin.key();
      if (colDeleted[col]) continue;
      const HighsCDouble& weightUpper = bin.value().weightUpper;
      const HighsCDouble& weightLower = bin.value().weightLower;
      if (weightUpper > primal_feastol &&
          weightUpper <= threshold + primal_feastol)
        candidates.push_back(
            {HighsCliqueTable::CliqueVar{col, 1}, weightUpper});
      if (weightLower > primal_feastol &&
          weightLower <= threshold + primal_feastol)
        candidates.push_back(
            {HighsCliqueTable::CliqueVar{col, 0}, weightLower});
    }

    // sort candidates by decreasing weight; equal weights are ordered by
    // literal index
    pdqsort(candidates.begin(), candidates.end(),
            [](const CliqueCandidate& a, const CliqueCandidate& b) {
              if (a.weight > b.weight) return true;
              if (b.weight > a.weight) return false;
              return a.var.index() < b.var.index();
            });

    for (size_t i = 0;
         i < candidates.size() && numCliquesAdded < maxNumCliquesAdded; ++i) {
      const auto& v1 = candidates[i].var;
      // fixed by the clique table
      if (colDeleted[v1.col]) continue;

      for (size_t j = i + 1; j < candidates.size(); ++j) {
        const auto& v2 = candidates[j].var;
        if (colDeleted[v2.col]) continue;

        // skip pair of variables if they are already in a clique together
        // (this includes complementary literals of the same column)
        numWork++;
        if (presolveCliqueTable.haveCommonClique(v1, v2)) continue;

        // combined weight cannot exceed threshold; since candidates are
        // sorted by decreasing weight, no later pair with v1 can either
        if (candidates[i].weight + candidates[j].weight <=
            threshold + primal_feastol)
          break;

        // the pair test looks up the implied bounds of both literals on the
        // non-binaries and queries the clique table for each row binary
        numWork +=
            2 * static_cast<int64_t>(nonBinNonZeros.size() + sortedBins.size());

        // "simulate" setting v1 and v2 active: what is the effect on the
        // constraint's minimum activity?
        //
        // non-binary contributions: when v1 or v2 is active, implied bounds
        // may tighten a non-binary's feasible range, increasing its minimum
        // contribution to the row. intersect the bounds implied by v1 and v2
        // and add the excess over the standard bound. contradicting implied
        // bounds are handled by findConflictingImplications at the start of
        // each round.
        HighsCDouble sum = 0;
        for (const auto& nz : nonBinNonZeros) {
          // discovered non-binaries outside the row (zero coefficient) do not
          // contribute to the row's activity
          if (nz.value().val == 0.0) continue;
          HighsInt col = nz.key();
          // start from the column's bounds; implied bounds may be looser
          HighsImplications::Implication intersection(model->col_lower_[col],
                                                      model->col_upper_[col]);
          for (const auto& v : {v1, v2}) {
            const auto* impl = impliedBounds.find(
                ImpliedBoundKey{static_cast<HighsInt>(v.col),
                                static_cast<HighsInt>(v.val), col});
            if (impl == nullptr) continue;
            intersection.lb = std::max(intersection.lb, impl->lb);
            intersection.ub = std::min(intersection.ub, impl->ub);
          }
          sum += computeLift(col, intersection);
        }

        // binary contributions: add |a_j| for each binary forced away from
        // its min-contribution value when v1 or v2 is active (via clique).
        // only binaries in the row contribute; sortedBins (built in
        // updateBinaryWeights) holds exactly these, while binNonZeros also
        // contains lifted binaries and clique neighbours outside the row
        for (const SortedBin& sortedBin : sortedBins) {
          HighsInt bcol = sortedBin.col;
          const VariableData& bin = *binNonZeros.find(bcol);
          double bval = bin.val;
          double absbval = std::abs(bval);

          // (col, val) setting that minimises this binary's contribution
          HighsCliqueTable::CliqueVar minContribCliqueVar(
              bcol, minContributionVal(bval));

          // binary is forced away from min-contribution if v1 or v2
          // is the complement, or shares a clique with it
          if (presolveCliqueTable.haveCommonClique(v1, minContribCliqueVar) ||
              presolveCliqueTable.haveCommonClique(v2, minContribCliqueVar)) {
            if (bin.member < 0) {
              sum += absbval;
            } else {
              // members of clique blocks are accounted for per block below
              blockMembers[bin.member].count++;
              touchedBlocks.add(blockMembers[bin.member].block);
            }
          }
        }
        for (HighsInt block : touchedBlocks.entries) {
          sum += blockIncrease(block);
          for (HighsInt k = blockStarts[block]; k < blockStarts[block + 1]; ++k)
            blockMembers[k].count = 0;
        }
        touchedBlocks.clear();

        // combined weight exceeds threshold: v1 and v2 form a clique
        if (sum > threshold + primal_feastol) {
          HPRESOLVE_CHECKED_CALL(addClique(v1, v2));
          consecutiveNoClique = 0;
          // fixed by the clique table
          if (colDeleted[v1.col]) break;
        } else {
          if (++consecutiveNoClique >= maxConsecutiveNoClique)
            return Result::kOk;
        }
      }
    }
    return Result::kOk;
  };

  auto collectBreakpoints = [&](HighsInt nonBinCol) {
    // the implied activity w(d) for non-binary x_r at value d is piecewise
    // linear: between breakpoints, it changes linearly with slope a_r (the
    // row coefficient). at breakpoints, discrete jumps occur because binary
    // implications activate or deactivate.
    //
    // a breakpoint arises from a binary->non-binary implication read in
    // reverse (contrapositive). for binary x_i at its min-contribution value v:
    //   (x_i = v -> x_r >= lambda)  implies  (x_r < lambda -> x_i != v)
    //   (x_i = v -> x_r <= mu)      implies  (x_r > mu     -> x_i != v)
    // forcing x_i away from its min-contribution value costs |a_i| in
    // activity.
    //
    // lower-type breakpoint at lambda: active when x_r < lambda (cost |a_i|)
    // upper-type breakpoint at mu:     active when x_r > mu     (cost |a_i|)
    breakpoints.clear();
    numWork += static_cast<int64_t>(binNonZeros.size());
    for (const auto& binVar : binNonZeros) {
      HighsInt binCol = binVar.key();
      double binVal = binVar.value().val;
      // binaries outside the row cost nothing when forced away from a value
      if (binVal == 0.0) continue;

      // look up implications from the binary's min-contribution value:
      // a_i > 0: min at x_i = 0, a_i < 0: min at x_i = 1
      HighsInt val = minContributionVal(binVal);
      const auto* impl =
          impliedBounds.find(ImpliedBoundKey{binCol, val, nonBinCol});
      if (impl == nullptr) continue;

      if (impl->lb > model->col_lower_[nonBinCol] + primal_feastol)
        breakpoints.push_back({impl->lb, binCol, true});
      if (impl->ub < model->col_upper_[nonBinCol] - primal_feastol)
        breakpoints.push_back({impl->ub, binCol, false});
    }
    // sort by increasing value so computeBound can walk from either end.
    // at equal values, lower-type breakpoints come first: neither type is
    // active at the breakpoint itself, so the walk has to deactivate before
    // it activates to detect a crossing at that point
    pdqsort(breakpoints.begin(), breakpoints.end(),
            [](const Breakpoint& a, const Breakpoint& b) {
              if (a.value != b.value) return a.value < b.value;
              return a.isLowerType && !b.isLowerType;
            });
  };

  auto activateBlockMembers = [&](bool lowerType) {
    // members of clique blocks with breakpoints of the given type leave their
    // min-contribution value; returns the increase of the minimum activity
    HighsCDouble increase = 0.0;
    for (const Breakpoint& bp : breakpoints) {
      HighsInt member = binNonZeros.find(bp.binCol)->member;
      if (member < 0 || bp.isLowerType != lowerType) continue;
      blockMembers[member].count++;
      touchedBlocks.add(blockMembers[member].block);
    }
    for (HighsInt block : touchedBlocks.entries)
      increase += blockIncrease(block);
    touchedBlocks.clear();
    return increase;
  };

  auto resetBlockMembers = [&]() {
    for (const Breakpoint& bp : breakpoints) {
      HighsInt member = binNonZeros.find(bp.binCol)->member;
      if (member >= 0) blockMembers[member].count = 0;
    }
  };

  auto tightenNonBinaryBounds = [&](double threshold) {
    // non-binary piecewise bound tightening
    for (auto& nz : nonBinNonZeros) {
      HighsInt nonBinCol = nz.key();
      double nonBinVal = nz.value().val;

      // weights incomplete when a bound is infinite (addNonZero does not add
      // the column's own range)
      if (model->col_lower_[nonBinCol] <= -kHighsInf ||
          model->col_upper_[nonBinCol] >= kHighsInf)
        continue;

      // collect breakpoints and update endpoint weights: at x_r = lb all
      // lower-type breakpoints are active (x_r < lambda holds), at x_r = ub
      // all upper-type breakpoints are active (x_r > mu holds)
      collectBreakpoints(nonBinCol);
      for (const auto& bp : breakpoints) {
        const VariableData& bin = *binNonZeros.find(bp.binCol);
        if (bin.member < 0)
          nz.value().update(!bp.isLowerType, std::abs(bin.val));
      }
      // members of clique blocks are accounted for per block
      if (!blockMembers.empty()) {
        for (bool lowerType : {true, false}) {
          nz.value().update(!lowerType, activateBlockMembers(lowerType));
          resetBlockMembers();
        }
      }

      // w(lb) > threshold means lb is infeasible; walk right to find new lb
      // w(ub) > threshold means ub is infeasible; walk left to find new ub
      const HighsCDouble& weightLower = nz.value().weightLower;
      const HighsCDouble& weightUpper = nz.value().weightUpper;
      bool tightenLower = weightLower > threshold + primal_feastol;
      bool tightenUpper = weightUpper > threshold + primal_feastol;
      if (!tightenLower && !tightenUpper) continue;

      // walk the piecewise linear implied activity w(d) from colBound toward
      // otherColBound, looking for where w drops from above threshold to
      // below. between breakpoints w changes with slope val (the non-binary
      // row coefficient). at each breakpoint a discrete jump occurs:
      // 1. walking right (+1): lower-type deactivates,
      //                        upper-type activates
      // 2. walking left  (-1): upper-type deactivates,
      //                        lower-type activates
      // the new bound is found by interpolation within a segment, or
      // snapped to the breakpoint where the jump crosses the threshold.
      auto computeBound = [&](double val, const HighsCDouble& inputWeight,
                              double colBound, double otherColBound,
                              HighsInt direction, double& newColBound) {
        HighsCDouble weight = inputWeight;
        double d = colBound;
        newColBound = colBound;

        // point where w crosses the threshold on the segment starting at
        // segStart with weight segWeight
        auto interpolate = [&](const HighsCDouble& segWeight, double segStart) {
          assert(val != 0.0);
          return static_cast<double>(segStart + (threshold - segWeight) / val);
        };

        HighsInt start = 0;
        HighsInt end = static_cast<HighsInt>(breakpoints.size()) - 1;
        if (direction < 0) std::swap(start, end);

        for (HighsInt i = start; direction * i <= direction * end;
             i += direction) {
          double bp = breakpoints[i].value;
          // skip breakpoints more than primal_feastol outside [d,
          // otherColBound] (consistent with the infeasible-literal check in
          // findConflictingImplications). breakpoints at otherColBound are not
          // skipped: the jump there can make otherColBound feasible
          if (direction * bp < direction * d - primal_feastol ||
              direction * bp > direction * otherColBound + primal_feastol)
            continue;

          // linear change from d to bp: w changes by val * (bp - d)
          HighsCDouble weightAtBreakpoint =
              weight + val * (bp - static_cast<HighsCDouble>(d));

          // threshold crossed mid-segment: interpolate exact crossing point
          if (weight > threshold + primal_feastol &&
              weightAtBreakpoint <= threshold + primal_feastol) {
            newColBound = interpolate(weight, d);
            return true;
          }

          weight = weightAtBreakpoint;

          // discrete jump at the breakpoint; direction multiplier handles
          // both cases:
          // 1. going right -> lower-type deactivates / upper activates,
          // 2. going left  -> upper-type deactivates / lower activates
          const VariableData& bin = *binNonZeros.find(breakpoints[i].binCol);
          if (bin.member >= 0) {
            // the jump of a block member depends on the other members of its
            // block that are forced away at the same time
            HighsInt block = blockMembers[bin.member].block;
            double increaseBefore = blockIncrease(block);
            bool activate = (direction > 0) != breakpoints[i].isLowerType;
            blockMembers[bin.member].count += activate ? 1 : -1;
            weight += blockIncrease(block) - increaseBefore;
          } else if (breakpoints[i].isLowerType)
            weight -= direction * std::abs(bin.val);
          else
            weight += direction * std::abs(bin.val);
          d = bp;

          // threshold crossed at the discrete jump: snap bound to bp.
          // must be checked here: the next iteration's interpolation check
          // requires weight > threshold at the segment start, which no longer
          // holds after the jump
          if (weight <= threshold + primal_feastol) {
            newColBound = bp;
            return true;
          }
        }

        // check the final segment from last breakpoint to otherColBound
        if (weight > threshold + primal_feastol) {
          HighsCDouble weightAtBound =
              weight + val * (otherColBound - static_cast<HighsCDouble>(d));
          if (weightAtBound <= threshold + primal_feastol) {
            newColBound = interpolate(weight, d);
            return true;
          }
        }
        return false;
      };

      // tighten lower bound: walk from lb toward ub. at lb, the block members
      // with lower-type breakpoints are forced away
      if (!colDeleted[nonBinCol] && tightenLower) {
        double newLowerBnd;
        if (!blockMembers.empty()) activateBlockMembers(true);
        bool feasible = computeBound(
            nonBinVal, weightLower, model->col_lower_[nonBinCol],
            model->col_upper_[nonBinCol], HighsInt{1}, newLowerBnd);
        if (!blockMembers.empty()) resetBlockMembers();
        // w(d) > threshold for all d in [lb, ub]: no feasible value
        if (!feasible) return Result::kPrimalInfeasible;
        HPRESOLVE_CHECKED_CALL(changeNonBinaryLower(nonBinCol, newLowerBnd));
      }

      // tighten upper bound: walk from ub toward lb. at ub, the block members
      // with upper-type breakpoints are forced away
      if (!colDeleted[nonBinCol] && tightenUpper) {
        double newUpperBnd;
        if (!blockMembers.empty()) activateBlockMembers(false);
        bool feasible = computeBound(
            nonBinVal, weightUpper, model->col_upper_[nonBinCol],
            model->col_lower_[nonBinCol], HighsInt{-1}, newUpperBnd);
        if (!blockMembers.empty()) resetBlockMembers();
        // w(d) > threshold for all d in [lb, ub]: no feasible value
        if (!feasible) return Result::kPrimalInfeasible;
        HPRESOLVE_CHECKED_CALL(changeNonBinaryUpper(nonBinCol, newUpperBnd));
      }
    }
    return Result::kOk;
  };

  auto checkRow = [&](double threshold) {
    numCheckRowCalls++;
    if (numCheckRowCalls % checkTimeFrequency == 0)
      HPRESOLVE_CHECKED_CALL(checkTimeLimit());
    // the min-contribution literals of the row binaries in a clique block
    // exclude each other, which increases the row's minimum activity
    threshold -= static_cast<double>(partitionRowBinaries());
    if (threshold < -primal_feastol) return Result::kPrimalInfeasible;

    // collect implied bounds on non-binaries
    collectImpliedBounds();

    // lift / update binary weights
    updateBinaryWeights(threshold);

    // binary fixing
    for (const auto& binVar : binNonZeros) {
      HighsInt col = binVar.key();
      if (colDeleted[col]) continue;

      bool lowerInfeasible =
          binVar.value().weightLower > threshold + primal_feastol;
      bool upperInfeasible =
          binVar.value().weightUpper > threshold + primal_feastol;
      // neither value of the binary satisfies the row
      if (lowerInfeasible && upperInfeasible) return Result::kPrimalInfeasible;
      if (lowerInfeasible || upperInfeasible)
        HPRESOLVE_CHECKED_CALL(fixCol(col, lowerInfeasible));
    }

    // two-column clique extraction
    HPRESOLVE_CHECKED_CALL(findCliques(threshold));

    // non-binary piecewise bound tightening
    HPRESOLVE_CHECKED_CALL(tightenNonBinaryBounds(threshold));

    return Result::kOk;
  };

  auto currentWorkLimit = [&]() {
    // reductions add to the work limit in the main presolve
    int64_t limit = workLimit;
    if (!mipsolver->submip)
      limit += workBonusPerBoundChange * (numVarsFixed + numBoundsTightened) +
               workBonusPerClique * numCliquesAdded;
    return limit;
  };

  auto isWorkLimitReached = [&]() { return numWork > currentWorkLimit(); };

  auto sortRowsByProbingScore = [&](std::vector<HighsInt>& rows) {
    // rows whose binaries have high probing scores (many implications) come
    // first, so that they are checked before the work limit is reached. the
    // row score is the average probing score of the row binaries, as in
    // enumeration; ties are broken by row index
    std::vector<std::tuple<double, double, HighsInt>> rowScores;
    rowScores.reserve(rows.size());
    for (HighsInt row : rows) {
      int64_t score = 0;
      HighsInt score2 = 0;
      HighsInt numBins = 0;
      for (const auto& nz : getRowVector(row)) {
        if (!isBinary(nz.index())) continue;
        auto probingScore = computeProbingScore(nz.index());
        score += probingScore.first;
        score2 += probingScore.second;
        numBins++;
      }
      double divisor = std::max(HighsInt{1}, numBins);
      rowScores.emplace_back(-score / divisor, -score2 / divisor, row);
    }
    pdqsort(rowScores.begin(), rowScores.end());
    rows.clear();
    for (const auto& rowScore : rowScores)
      rows.push_back(std::get<2>(rowScore));
  };

  // prepare vectors
  modifiedRows.resize(model->num_row_);
  changedNonBinaries.resize(model->num_col_);
  binNonZeros.resize(model->num_col_);
  nonBinNonZeros.resize(model->num_col_);
  cachedColumns.resize(model->num_col_);
  literalPos.resize(2 * model->num_col_, -1);

  // packed storage for objective
  if (mipsolver->mipdata_->upper_limit < kHighsInf ||
      mipsolver->mipdata_->lower_bound > -kHighsInf) {
    objInds.reserve(model->num_col_);
    for (HighsInt col = 0; col < model->num_col_; col++) {
      if (colDeleted[col] || model->col_cost_[col] == 0.0) continue;
      objInds.push_back(col);
    }
    // the objective is only checked if it is not too large
    checkObjective = static_cast<HighsInt>(objInds.size()) <= maxObjSize;
  }

  // initialise
  bool myObjectiveAffected = checkObjective;
  std::vector<HighsInt> myModifiedRows;
  myModifiedRows.reserve(model->num_row_);
  for (HighsInt row = 0; row < model->num_row_; row++)
    if (!rowDeleted[row] && rowsize[row] > 1 && rowsize[row] <= maxRowSize)
      myModifiedRows.push_back(row);
  std::vector<HighsInt> myChangedNonBinaries;
  myChangedNonBinaries.reserve(model->num_col_);
  for (HighsInt col = 0; col < model->num_col_; col++)
    if (!colDeleted[col] && !isBinary(col)) myChangedNonBinaries.push_back(col);

  // check rows
  while (myObjectiveAffected || !myModifiedRows.empty() ||
         !myChangedNonBinaries.empty()) {
    // stop when the work limit is reached; it is checked again before each
    // row. the checks of conflicting implications below do not count as work
    if (isWorkLimitReached()) break;

    // literals implying contradicting bounds on a non-binary do not depend on
    // a row: check all non-binaries in the first round, afterwards those with
    // changed bounds
    for (HighsInt col : myChangedNonBinaries) {
      if (colDeleted[col] || isBinary(col)) continue;
      conflictingLiterals.clear();
      infeasibleLiterals.clear();
      double impliedLower;
      double impliedUpper;
      implications.findConflictingImplications(
          col, model->col_lower_[col], model->col_upper_[col],
          model->integrality_[col] != HighsVarType::kContinuous,
          conflictingLiterals, infeasibleLiterals, impliedLower, impliedUpper,
          maxNumCliquesAdded - numCliquesAdded);
      for (const HighsCliqueTable::CliqueVar& literal : infeasibleLiterals) {
        if (colDeleted[literal.col] || !isBinary(literal.col)) continue;
        HPRESOLVE_CHECKED_CALL(fixCol(literal.col, literal.val == 0));
      }
      for (const auto& conflict : conflictingLiterals) {
        const HighsCliqueTable::CliqueVar& v1 = conflict.first;
        const HighsCliqueTable::CliqueVar& v2 = conflict.second;
        if (colDeleted[v1.col] || colDeleted[v2.col] || !isBinary(v1.col) ||
            !isBinary(v2.col))
          continue;
        if (presolveCliqueTable.haveCommonClique(v1, v2)) continue;
        if (numCliquesAdded >= maxNumCliquesAdded) break;
        HPRESOLVE_CHECKED_CALL(addClique(v1, v2));
      }
      // bounds implied by the union of the ranges of both literals of a binary
      if (colDeleted[col]) continue;
      if (impliedLower > impliedUpper + primal_feastol)
        return Result::kPrimalInfeasible;
      HPRESOLVE_CHECKED_CALL(changeNonBinaryLower(col, impliedLower));
      HPRESOLVE_CHECKED_CALL(changeNonBinaryUpper(col, impliedUpper));
    }

    // check objective function
    if (myObjectiveAffected) {
      HighsCDouble minActivity;
      if (mipsolver->mipdata_->upper_limit < kHighsInf &&
          loadObjective(HighsInt{1}, minActivity)) {
        numObjectiveChecks++;
        HPRESOLVE_CHECKED_CALL(checkRow(static_cast<double>(
            static_cast<HighsCDouble>(mipsolver->mipdata_->upper_limit) -
            model->offset_ - minActivity)));
      }
      if (mipsolver->mipdata_->lower_bound > -kHighsInf &&
          loadObjective(HighsInt{-1}, minActivity)) {
        numObjectiveChecks++;
        HPRESOLVE_CHECKED_CALL(checkRow(static_cast<double>(
            static_cast<HighsCDouble>(model->offset_) -
            mipsolver->mipdata_->lower_bound - minActivity)));
      }
    }

    sortRowsByProbingScore(myModifiedRows);
    for (HighsInt row : myModifiedRows) {
      if (rowDeleted[row]) continue;
      if (isWorkLimitReached()) break;
      if (model->row_upper_[row] < kHighsInf && loadModelRow(row, HighsInt{1}))
        HPRESOLVE_CHECKED_CALL(checkRow(
            -impliedRowBounds.getSumLowerOrig(row, -model->row_upper_[row])));
      if (model->row_lower_[row] > -kHighsInf &&
          loadModelRow(row, HighsInt{-1}))
        HPRESOLVE_CHECKED_CALL(checkRow(
            impliedRowBounds.getSumUpperOrig(row, -model->row_lower_[row])));
    }

    // re-check changed rows and non-binaries
    myObjectiveAffected = objectiveAffected;
    objectiveAffected = false;
    modifiedRows.moveEntriesTo(myModifiedRows);
    changedNonBinaries.moveEntriesTo(myChangedNonBinaries);
  }

  highsLogDev(options->log_options, HighsLogType::kInfo,
              "Implication-aware constraint propagation: %" HIGHSINT_FORMAT
              " checks (%" HIGHSINT_FORMAT
              " of the objective), fixed %" HIGHSINT_FORMAT
              " columns, added %" HIGHSINT_FORMAT
              " cliques, tightened %" HIGHSINT_FORMAT " bounds, used %" PRId64
              " of %" PRId64 " work units%s\n",
              numCheckRowCalls, numObjectiveChecks, numVarsFixed,
              numCliquesAdded, numBoundsTightened, numWork, currentWorkLimit(),
              isWorkLimitReached() ? " (limit reached)" : "");

  return checkLimits(postsolve_stack);
}

}  // namespace presolve
