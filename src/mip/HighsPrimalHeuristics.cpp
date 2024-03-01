/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsPrimalHeuristics.h"

#include <numeric>
#include <unordered_set>

//#include "io/HighsIO.h"
//#include "lp_data/HConst.h"
//#include "lp_data/HighsLpUtils.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsDomainChange.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "pdqsort/pdqsort.h"
//#include "util/HighsHash.h"
#include "util/HighsIntegers.h"

// GCC floating point errors are well-known for 32-bit architectures;
// see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=323.
// An easy workaround is to add the "volatile" keyword to avoid
// problematic GCC optimizations that impact precision.
#ifdef __i386__
#define FP_32BIT_VOLATILE volatile
#else
#define FP_32BIT_VOLATILE
#endif

HighsPrimalHeuristics::HighsPrimalHeuristics(HighsMipSolver& mipsolver)
    : mipsolver(mipsolver),
      lp_iterations(0),
      randgen(mipsolver.options_mip_->random_seed) {
  successObservations = 0;
  numSuccessObservations = 0;
  infeasObservations = 0;
  numInfeasObservations = 0;
}

void HighsPrimalHeuristics::setupIntCols() {
  intcols = mipsolver.mipdata_->integer_cols;

  pdqsort(intcols.begin(), intcols.end(), [&](HighsInt c1, HighsInt c2) {
    const FP_32BIT_VOLATILE double lockScore1 =
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->uplocks[c1]) *
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->downlocks[c1]);

    const FP_32BIT_VOLATILE double lockScore2 =
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->uplocks[c2]) *
        (mipsolver.mipdata_->feastol + mipsolver.mipdata_->downlocks[c2]);

    if (lockScore1 > lockScore2) return true;
    if (lockScore2 > lockScore1) return false;

    const FP_32BIT_VOLATILE double cliqueScore1 =
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c1, 1)) *
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c1, 0));

    const FP_32BIT_VOLATILE double cliqueScore2 =
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c2, 1)) *
        (mipsolver.mipdata_->feastol +
         mipsolver.mipdata_->cliquetable.getNumImplications(c2, 0));

    return std::make_tuple(cliqueScore1, HighsHashHelpers::hash(uint64_t(c1)),
                           c1) >
           std::make_tuple(cliqueScore2, HighsHashHelpers::hash(uint64_t(c2)),
                           c2);
  });
}

bool HighsPrimalHeuristics::solveSubMip(
    const HighsLp& lp, const HighsBasis& basis, double fixingRate,
    std::vector<double> colLower, std::vector<double> colUpper,
    HighsInt maxleaves, HighsInt maxnodes, HighsInt stallnodes) {
  HighsOptions submipoptions = *mipsolver.options_mip_;
  HighsLp submip = lp;

  // set bounds and restore integrality of the lp relaxation copy
  submip.col_lower_ = std::move(colLower);
  submip.col_upper_ = std::move(colUpper);
  submip.integrality_ = mipsolver.model_->integrality_;
  submip.offset_ = 0;

  // set limits
  submipoptions.mip_max_leaves = maxleaves;
  submipoptions.output_flag = false;

  const bool allow_submip_log = true;
  if (allow_submip_log && lp.num_col_ == -54 && lp.num_row_ == -172) {
    submipoptions.output_flag = true;
    printf(
        "HighsPrimalHeuristics::solveSubMip (%d, %d) with output_flag = %s\n",
        int(lp.num_col_), int(lp.num_row_),
        highsBoolToString(submipoptions.output_flag).c_str());
  }

  submipoptions.mip_max_nodes = maxnodes;
  submipoptions.mip_max_stall_nodes = stallnodes;
  submipoptions.mip_pscost_minreliable = 0;
  submipoptions.time_limit -=
      mipsolver.timer_.read(mipsolver.timer_.solve_clock);
  submipoptions.objective_bound = mipsolver.mipdata_->upper_limit;

  if (!mipsolver.submip) {
    double curr_abs_gap =
        mipsolver.mipdata_->upper_limit - mipsolver.mipdata_->lower_bound;

    if (curr_abs_gap == kHighsInf) {
      curr_abs_gap = fabs(mipsolver.mipdata_->lower_bound);
      if (curr_abs_gap == kHighsInf) curr_abs_gap = 0.0;
    }

    submipoptions.mip_rel_gap = 0.0;
    submipoptions.mip_abs_gap =
        mipsolver.mipdata_->feastol * std::max(curr_abs_gap, 1000.0);
  }

  submipoptions.presolve = "on";
  submipoptions.mip_detect_symmetry = false;
  submipoptions.mip_heuristic_effort = 0.8;
  // setup solver and run it

  HighsSolution solution;
  solution.value_valid = false;
  solution.dual_valid = false;
  // Create HighsMipSolver instance for sub-MIP
  HighsMipSolver submipsolver(*mipsolver.callback_, submipoptions, submip,
                              solution, true, mipsolver.submip_level + 1);

  submipsolver.rootbasis = &basis;
  HighsPseudocostInitialization pscostinit(mipsolver.mipdata_->pseudocost, 1);
  submipsolver.pscostinit = &pscostinit;
  submipsolver.clqtableinit = &mipsolver.mipdata_->cliquetable;
  submipsolver.implicinit = &mipsolver.mipdata_->implications;
  // Solve the sub-MIP
  submipsolver.run();
  mipsolver.max_submip_level =
      std::max(submipsolver.max_submip_level + 1, mipsolver.max_submip_level);
  if (submipsolver.mipdata_) {
    double numUnfixed = mipsolver.mipdata_->integral_cols.size() +
                        mipsolver.mipdata_->continuous_cols.size();
    double adjustmentfactor = submipsolver.numCol() / std::max(1.0, numUnfixed);
    // (double)mipsolver.orig_model_->a_matrix_.value_.size();
    int64_t adjusted_lp_iterations =
        (size_t)(adjustmentfactor * submipsolver.mipdata_->total_lp_iterations);
    lp_iterations += adjusted_lp_iterations;

    if (mipsolver.submip)
      mipsolver.mipdata_->num_nodes += std::max(
          int64_t{1}, int64_t(adjustmentfactor * submipsolver.node_count_));
    TrivialHeuristicData& mipsolver_submip_statistics =
        mipsolver.mipdata_->submip_trivial_heuristics_statistics_;
    upCopyLocalTrivialHeuristicsStatistics(mipsolver_submip_statistics);
  }

  if (submipsolver.modelstatus_ == HighsModelStatus::kInfeasible) {
    infeasObservations += fixingRate;
    ++numInfeasObservations;
  }
  if (submipsolver.node_count_ <= 1 &&
      submipsolver.modelstatus_ == HighsModelStatus::kInfeasible)
    return false;
  HighsInt oldNumImprovingSols = mipsolver.mipdata_->numImprovingSols;
  if (submipsolver.modelstatus_ != HighsModelStatus::kInfeasible &&
      !submipsolver.solution_.empty()) {
    mipsolver.mipdata_->trySolution(submipsolver.solution_, 'L');
  }

  if (mipsolver.mipdata_->numImprovingSols != oldNumImprovingSols) {
    // remember fixing rate as good
    successObservations += fixingRate;
    ++numSuccessObservations;
  }

  return true;
}

double HighsPrimalHeuristics::determineTargetFixingRate() {
  double lowFixingRate = 0.6;
  double highFixingRate = 0.6;

  if (numInfeasObservations != 0) {
    double infeasRate = infeasObservations / numInfeasObservations;
    highFixingRate = 0.9 * infeasRate;
    lowFixingRate = std::min(lowFixingRate, highFixingRate);
  }

  if (numSuccessObservations != 0) {
    double successFixingRate = successObservations / numSuccessObservations;
    lowFixingRate = std::min(lowFixingRate, 0.9 * successFixingRate);
    highFixingRate = std::max(successFixingRate * 1.1, highFixingRate);
  }

  double fixingRate = randgen.real(lowFixingRate, highFixingRate);
  // if (!mipsolver.submip) printf("fixing rate: %.2f\n", 100.0 * fixingRate);
  return fixingRate;
}

class HeuristicNeighbourhood {
  HighsDomain& localdom;
  HighsInt numFixed;
  HighsHashTable<HighsInt> fixedCols;
  size_t startCheckedChanges;
  size_t nCheckedChanges;
  HighsInt numTotal;

 public:
  HeuristicNeighbourhood(HighsMipSolver& mipsolver, HighsDomain& localdom)
      : localdom(localdom),
        numFixed(0),
        startCheckedChanges(localdom.getDomainChangeStack().size()),
        nCheckedChanges(startCheckedChanges) {
    for (HighsInt i : mipsolver.mipdata_->integral_cols)
      if (localdom.col_lower_[i] == localdom.col_upper_[i]) ++numFixed;

    numTotal = mipsolver.mipdata_->integral_cols.size() - numFixed;
  }

  double getFixingRate() {
    while (nCheckedChanges < localdom.getDomainChangeStack().size()) {
      HighsInt col = localdom.getDomainChangeStack()[nCheckedChanges++].column;
      if (localdom.variableType(col) == HighsVarType::kContinuous) continue;
      if (localdom.isFixed(col)) fixedCols.insert(col);
    }

    return numTotal ? static_cast<double>(fixedCols.size()) /
                          static_cast<double>(numTotal)
                    : 0.0;
  }

  void backtracked() {
    nCheckedChanges = startCheckedChanges;
    if (fixedCols.size()) fixedCols.clear();
  }
};

void HighsPrimalHeuristics::rootReducedCost() {
  std::vector<std::pair<double, HighsDomainChange>> lurkingBounds =
      mipsolver.mipdata_->redcostfixing.getLurkingBounds(mipsolver);
  if (10 * lurkingBounds.size() < mipsolver.mipdata_->integral_cols.size())
    return;
  pdqsort(lurkingBounds.begin(), lurkingBounds.end(),
          [](const std::pair<double, HighsDomainChange>& a,
             const std::pair<double, HighsDomainChange>& b) {
            return a.first > b.first;
          });

  auto localdom = mipsolver.mipdata_->domain;

  HeuristicNeighbourhood neighbourhood(mipsolver, localdom);

  double currCutoff = kHighsInf;
  double lower_bound;

  lower_bound = mipsolver.mipdata_->lower_bound + mipsolver.mipdata_->feastol;

  for (const std::pair<double, HighsDomainChange>& domchg : lurkingBounds) {
    currCutoff = domchg.first;

    if (currCutoff <= lower_bound) break;

    if (localdom.isActive(domchg.second)) continue;
    localdom.changeBound(domchg.second);

    while (true) {
      localdom.propagate();
      if (localdom.infeasible()) {
        localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
        mipsolver.mipdata_->lower_bound =
            std::max(mipsolver.mipdata_->lower_bound, currCutoff);
        localdom.backtrack();
        if (localdom.getBranchDepth() == 0) break;
        neighbourhood.backtracked();
        continue;
      }
      break;
    }
    double fixingRate = neighbourhood.getFixingRate();
    if (fixingRate >= 0.5) break;
    // double gap = (currCutoff - mipsolver.mipdata_->lower_bound) /
    //             std::max(std::abs(mipsolver.mipdata_->lower_bound), 1.0);
    // if (gap < 0.001) break;
  }

  double fixingRate = neighbourhood.getFixingRate();
  if (fixingRate < 0.3) return;

  solveSubMip(*mipsolver.model_, mipsolver.mipdata_->firstrootbasis, fixingRate,
              localdom.col_lower_, localdom.col_upper_,
              500,  // std::max(50, int(0.05 *
                    // (mipsolver.mipdata_->num_leaves))),
              200 + mipsolver.mipdata_->num_nodes / 20, 12);
}

void HighsPrimalHeuristics::RENS(const std::vector<double>& tmp) {
  HighsPseudocost pscost(mipsolver.mipdata_->pseudocost);
  HighsSearch heur(mipsolver, pscost);
  HighsDomain& localdom = heur.getLocalDomain();
  heur.setHeuristic(true);

  intcols.erase(std::remove_if(intcols.begin(), intcols.end(),
                               [&](HighsInt i) {
                                 return mipsolver.mipdata_->domain.isFixed(i);
                               }),
                intcols.end());

  HighsLpRelaxation heurlp(mipsolver.mipdata_->lp);
  // only use the global upper limit as LP limit so that dual proofs are valid
  heurlp.setObjectiveLimit(mipsolver.mipdata_->upper_limit);
  heurlp.setAdjustSymmetricBranchingCol(false);
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                        localdom.col_lower_.data(),
                                        localdom.col_upper_.data());
  localdom.clearChangedCols();
  heur.createNewNode();

  // determine the initial number of unfixed variables fixing rate to decide if
  // the problem is restricted enough to be considered for solving a submip
  double maxfixingrate = determineTargetFixingRate();
  double fixingrate = 0.0;
  bool stop = false;
  // heurlp.setIterationLimit(2 * mipsolver.mipdata_->maxrootlpiters);
  // printf("iterlimit: %" HIGHSINT_FORMAT "\n",
  //       heurlp.getLpSolver().getOptions().simplex_iteration_limit);
  HighsInt targetdepth = 1;
  HighsInt nbacktracks = -1;
  HeuristicNeighbourhood neighbourhood(mipsolver, localdom);
retry:
  ++nbacktracks;
  neighbourhood.backtracked();
  // printf("current depth : %" HIGHSINT_FORMAT
  //        "   target depth : %" HIGHSINT_FORMAT "\n",
  //        heur.getCurrentDepth(), targetdepth);
  if (heur.getCurrentDepth() > targetdepth) {
    if (!heur.backtrackUntilDepth(targetdepth)) {
      lp_iterations += heur.getLocalLpIterations();
      return;
    }
  }

  // printf("fixingrate before loop is %g\n", fixingrate);
  assert(heur.hasNode());
  while (true) {
    // printf("evaluating node\n");
    heur.evaluateNode();
    // printf("done evaluating node\n");
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      if (mipsolver.mipdata_->domain.infeasible()) {
        lp_iterations += heur.getLocalLpIterations();
        return;
      }

      if (!heur.backtrack()) break;
      neighbourhood.backtracked();
      continue;
    }

    fixingrate = neighbourhood.getFixingRate();
    // printf("after evaluating node current fixingrate is %g\n", fixingrate);
    if (fixingrate >= maxfixingrate) break;
    if (stop) break;
    if (nbacktracks >= 10) break;

    HighsInt numBranched = 0;
    double stopFixingRate = std::min(
        1.0 - (1.0 - neighbourhood.getFixingRate()) * 0.9, maxfixingrate);
    const auto& relaxationsol = heurlp.getSolution().col_value;
    for (HighsInt i : intcols) {
      if (localdom.col_lower_[i] == localdom.col_upper_[i]) continue;

      double downval =
          std::floor(relaxationsol[i] + mipsolver.mipdata_->feastol);
      double upval = std::ceil(relaxationsol[i] - mipsolver.mipdata_->feastol);

      downval = std::min(downval, localdom.col_upper_[i]);
      upval = std::max(upval, localdom.col_lower_[i]);
      if (localdom.col_lower_[i] < downval) {
        ++numBranched;
        heur.branchUpwards(i, downval, downval - 0.5);
        localdom.propagate();
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
          break;
        }
      }
      if (localdom.col_upper_[i] > upval) {
        ++numBranched;
        heur.branchDownwards(i, upval, upval + 0.5);
        localdom.propagate();
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
          break;
        }
      }

      if (neighbourhood.getFixingRate() >= stopFixingRate) break;
    }

    if (numBranched == 0) {
      auto getFixVal = [&](HighsInt col, double fracval) {
        double fixval;

        // reinforce direction of this solution away from root
        // solution if the change is at least 0.4
        // otherwise take the direction where the objective gets worse
        // if objective is zero round to nearest integer
        double rootchange = mipsolver.mipdata_->rootlpsol.empty()
                                ? 0.0
                                : fracval - mipsolver.mipdata_->rootlpsol[col];
        if (rootchange >= 0.4)
          fixval = std::ceil(fracval);
        else if (rootchange <= -0.4)
          fixval = std::floor(fracval);
        if (mipsolver.model_->col_cost_[col] > 0.0)
          fixval = std::ceil(fracval);
        else if (mipsolver.model_->col_cost_[col] < 0.0)
          fixval = std::floor(fracval);
        else
          fixval = std::floor(fracval + 0.5);
        // make sure we do not set an infeasible domain
        fixval = std::min(localdom.col_upper_[col], fixval);
        fixval = std::max(localdom.col_lower_[col], fixval);
        return fixval;
      };

      pdqsort(heurlp.getFractionalIntegers().begin(),
              heurlp.getFractionalIntegers().end(),
              [&](const std::pair<HighsInt, double>& a,
                  const std::pair<HighsInt, double>& b) {
                return std::make_pair(
                           std::abs(getFixVal(a.first, a.second) - a.second),
                           HighsHashHelpers::hash(
                               (uint64_t(a.first) << 32) +
                               heurlp.getFractionalIntegers().size())) <
                       std::make_pair(
                           std::abs(getFixVal(b.first, b.second) - b.second),
                           HighsHashHelpers::hash(
                               (uint64_t(b.first) << 32) +
                               heurlp.getFractionalIntegers().size()));
              });

      double change = 0.0;
      // select a set of fractional variables to fix
      for (auto fracint : heurlp.getFractionalIntegers()) {
        double fixval = getFixVal(fracint.first, fracint.second);

        if (localdom.col_lower_[fracint.first] < fixval) {
          ++numBranched;
          heur.branchUpwards(fracint.first, fixval, fracint.second);
          localdom.propagate();
          if (localdom.infeasible()) {
            localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
            break;
          }

          fixingrate = neighbourhood.getFixingRate();
        }

        if (localdom.col_upper_[fracint.first] > fixval) {
          ++numBranched;
          heur.branchDownwards(fracint.first, fixval, fracint.second);
          localdom.propagate();
          if (localdom.infeasible()) {
            localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
            break;
          }

          fixingrate = neighbourhood.getFixingRate();
        }

        if (fixingrate >= maxfixingrate) break;

        change += std::abs(fixval - fracint.second);
        if (change >= 0.5) break;
      }
    }

    if (numBranched == 0) break;
    heurlp.flushDomain(localdom);
  }

  // printf("stopped heur dive with fixing rate %g\n", fixingrate);
  // if there is no node left it means we backtracked to the global domain and
  // the subproblem was solved with the dive
  if (!heur.hasNode()) {
    lp_iterations += heur.getLocalLpIterations();
    return;
  }
  // determine the fixing rate to decide if the problem is restricted enough to
  // be considered for solving a submip

  fixingrate = neighbourhood.getFixingRate();
  // printf("fixing rate is %g\n", fixingrate);
  if (fixingrate < 0.1 ||
      (mipsolver.submip && mipsolver.mipdata_->numImprovingSols != 0)) {
    // heur.childselrule = ChildSelectionRule::kBestCost;
    heur.setMinReliable(0);
    heur.solveDepthFirst(10);
    lp_iterations += heur.getLocalLpIterations();
    if (mipsolver.submip) mipsolver.mipdata_->num_nodes += heur.getLocalNodes();
    // lpiterations += heur.lpiterations;
    // pseudocost = heur.pseudocost;
    return;
  }

  heurlp.removeObsoleteRows(false);
  if (!solveSubMip(heurlp.getLp(), heurlp.getLpSolver().getBasis(), fixingrate,
                   localdom.col_lower_, localdom.col_upper_,
                   500,  // std::max(50, int(0.05 *
                         // (mipsolver.mipdata_->num_leaves))),
                   200 + mipsolver.mipdata_->num_nodes / 20, 12)) {
    int64_t new_lp_iterations = lp_iterations + heur.getLocalLpIterations();
    if (new_lp_iterations + mipsolver.mipdata_->heuristic_lp_iterations >
        100000 + ((mipsolver.mipdata_->total_lp_iterations -
                   mipsolver.mipdata_->heuristic_lp_iterations -
                   mipsolver.mipdata_->sb_lp_iterations) >>
                  1)) {
      lp_iterations = new_lp_iterations;
      return;
    }

    targetdepth = heur.getCurrentDepth() / 2;
    if (targetdepth <= 1 || mipsolver.mipdata_->checkLimits()) {
      lp_iterations = new_lp_iterations;
      return;
    }
    maxfixingrate = fixingrate * 0.5;
    // printf("infeasible in root node, trying with lower fixing rate %g\n",
    //        maxfixingrate);
    goto retry;
  }

  lp_iterations += heur.getLocalLpIterations();
}

void HighsPrimalHeuristics::RINS(const std::vector<double>& relaxationsol) {
  if (int(relaxationsol.size()) != mipsolver.numCol()) return;

  intcols.erase(std::remove_if(intcols.begin(), intcols.end(),
                               [&](HighsInt i) {
                                 return mipsolver.mipdata_->domain.isFixed(i);
                               }),
                intcols.end());

  HighsPseudocost pscost(mipsolver.mipdata_->pseudocost);
  HighsSearch heur(mipsolver, pscost);
  HighsDomain& localdom = heur.getLocalDomain();
  heur.setHeuristic(true);

  HighsLpRelaxation heurlp(mipsolver.mipdata_->lp);
  // only use the global upper limit as LP limit so that dual proofs are valid
  heurlp.setObjectiveLimit(mipsolver.mipdata_->upper_limit);
  heurlp.setAdjustSymmetricBranchingCol(false);
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                        localdom.col_lower_.data(),
                                        localdom.col_upper_.data());
  localdom.clearChangedCols();
  heur.createNewNode();

  // determine the initial number of unfixed variables fixing rate to decide if
  // the problem is restricted enough to be considered for solving a submip
  double maxfixingrate = determineTargetFixingRate();
  double minfixingrate = 0.25;
  double fixingrate = 0.0;
  bool stop = false;
  HighsInt nbacktracks = -1;
  HighsInt targetdepth = 1;
  HeuristicNeighbourhood neighbourhood(mipsolver, localdom);
retry:
  ++nbacktracks;
  neighbourhood.backtracked();
  // printf("current depth : %" HIGHSINT_FORMAT "   target depth : %"
  // HIGHSINT_FORMAT "\n", heur.getCurrentDepth(),
  //       targetdepth);
  if (heur.getCurrentDepth() > targetdepth) {
    if (!heur.backtrackUntilDepth(targetdepth)) {
      lp_iterations += heur.getLocalLpIterations();
      return;
    }
  }

  assert(heur.hasNode());

  while (true) {
    heur.evaluateNode();
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      // printf("backtrack1\n");
      if (mipsolver.mipdata_->domain.infeasible()) {
        lp_iterations += heur.getLocalLpIterations();
        return;
      }

      if (!heur.backtrack()) break;
      neighbourhood.backtracked();
      continue;
    }

    fixingrate = neighbourhood.getFixingRate();

    if (stop) break;
    if (fixingrate >= maxfixingrate) break;
    if (nbacktracks >= 10) break;

    std::vector<std::pair<HighsInt, double>>::iterator fixcandend;

    // partition the fractional variables to consider which ones should we fix
    // in this dive first if there is an incumbent, we dive towards the RINS
    // neighbourhood
    fixcandend = std::partition(
        heurlp.getFractionalIntegers().begin(),
        heurlp.getFractionalIntegers().end(),
        [&](const std::pair<HighsInt, double>& fracvar) {
          return std::abs(relaxationsol[fracvar.first] -
                          mipsolver.mipdata_->incumbent[fracvar.first]) <=
                 mipsolver.mipdata_->feastol;
        });

    bool fixtolpsol = true;

    auto getFixVal = [&](HighsInt col, double fracval) {
      double fixval;
      if (fixtolpsol) {
        // RINS neighbourhood (with extension)
        fixval = std::floor(relaxationsol[col] + 0.5);
      } else {
        // reinforce direction of this solution away from root
        // solution if the change is at least 0.4
        // otherwise take the direction where the objective gets worse
        // if objective is zero round to nearest integer
        double rootchange = fracval - mipsolver.mipdata_->rootlpsol[col];
        if (rootchange >= 0.4)
          fixval = std::ceil(fracval);
        else if (rootchange <= -0.4)
          fixval = std::floor(fracval);
        if (mipsolver.model_->col_cost_[col] > 0.0)
          fixval = std::ceil(fracval);
        else if (mipsolver.model_->col_cost_[col] < 0.0)
          fixval = std::floor(fracval);
        else
          fixval = std::floor(fracval + 0.5);
      }
      // make sure we do not set an infeasible domain
      fixval = std::min(localdom.col_upper_[col], fixval);
      fixval = std::max(localdom.col_lower_[col], fixval);
      return fixval;
    };

    // no candidates left to fix for getting to the neighbourhood, therefore we
    // switch to a different diving strategy until the minimal fixing rate is
    // reached
    HighsInt numBranched = 0;
    if (heurlp.getFractionalIntegers().begin() == fixcandend) {
      fixingrate = neighbourhood.getFixingRate();
      double stopFixingRate =
          std::min(maxfixingrate, 1.0 - (1.0 - fixingrate) * 0.9);
      const auto& currlpsol = heurlp.getSolution().col_value;
      for (HighsInt i : intcols) {
        if (localdom.col_lower_[i] == localdom.col_upper_[i]) continue;

        if (std::abs(currlpsol[i] - mipsolver.mipdata_->incumbent[i]) <=
            mipsolver.mipdata_->feastol) {
          double fixval = HighsIntegers::nearestInteger(currlpsol[i]);
          if (localdom.col_lower_[i] < fixval) {
            ++numBranched;
            heur.branchUpwards(i, fixval, fixval - 0.5);
            localdom.propagate();
            if (localdom.infeasible()) {
              localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
              break;
            }

            fixingrate = neighbourhood.getFixingRate();
          }
          if (localdom.col_upper_[i] > fixval) {
            ++numBranched;
            heur.branchDownwards(i, fixval, fixval + 0.5);
            localdom.propagate();
            if (localdom.infeasible()) {
              localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
              break;
            }

            fixingrate = neighbourhood.getFixingRate();
          }

          if (fixingrate >= stopFixingRate) break;
        }
      }

      if (numBranched != 0) {
        // printf(
        //    "fixed %" HIGHSINT_FORMAT " additional cols, old fixing rate:
        //    %.2f%%, new fixing " "rate: %.2f%%\n", numBranched, fixingrate,
        //    getFixingRate());
        heurlp.flushDomain(localdom);
        continue;
      }

      if (fixingrate >= minfixingrate)
        break;  // if the RINS neighbourhood achieved a high enough fixing rate
                // by itself we stop here
      fixcandend = heurlp.getFractionalIntegers().end();
      // now sort the variables by their distance towards the value they will
      // be fixed to
      fixtolpsol = false;
    }

    // now sort the variables by their distance towards the value they will be
    // fixed to
    pdqsort(heurlp.getFractionalIntegers().begin(), fixcandend,
            [&](const std::pair<HighsInt, double>& a,
                const std::pair<HighsInt, double>& b) {
              return std::make_pair(
                         std::abs(getFixVal(a.first, a.second) - a.second),
                         HighsHashHelpers::hash(
                             (uint64_t(a.first) << 32) +
                             heurlp.getFractionalIntegers().size())) <
                     std::make_pair(
                         std::abs(getFixVal(b.first, b.second) - b.second),
                         HighsHashHelpers::hash(
                             (uint64_t(b.first) << 32) +
                             heurlp.getFractionalIntegers().size()));
            });

    double change = 0.0;
    // select a set of fractional variables to fix
    for (auto fracint = heurlp.getFractionalIntegers().begin();
         fracint != fixcandend; ++fracint) {
      double fixval = getFixVal(fracint->first, fracint->second);

      if (localdom.col_lower_[fracint->first] < fixval) {
        ++numBranched;
        heur.branchUpwards(fracint->first, fixval, fracint->second);
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
          break;
        }

        fixingrate = neighbourhood.getFixingRate();
      }

      if (localdom.col_upper_[fracint->first] > fixval) {
        ++numBranched;
        heur.branchDownwards(fracint->first, fixval, fracint->second);
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
          break;
        }

        fixingrate = neighbourhood.getFixingRate();
      }

      if (fixingrate >= maxfixingrate) break;

      change += std::abs(fixval - fracint->second);
      if (change >= 0.5) break;
    }

    if (numBranched == 0) break;

    heurlp.flushDomain(localdom);

    // printf("%" HIGHSINT_FORMAT "/%" HIGHSINT_FORMAT " fixed, fixingrate is
    // %g\n", nfixed, ntotal, fixingrate);
  }

  // if there is no node left it means we backtracked to the global domain and
  // the subproblem was solved with the dive
  if (!heur.hasNode()) {
    lp_iterations += heur.getLocalLpIterations();
    return;
  }
  // determine the fixing rate to decide if the problem is restricted enough
  // to be considered for solving a submip

  // printf("fixing rate is %g\n", fixingrate);
  fixingrate = neighbourhood.getFixingRate();
  if (fixingrate < 0.1 ||
      (mipsolver.submip && mipsolver.mipdata_->numImprovingSols != 0)) {
    // heur.childselrule = ChildSelectionRule::kBestCost;
    heur.setMinReliable(0);
    heur.solveDepthFirst(10);
    lp_iterations += heur.getLocalLpIterations();
    if (mipsolver.submip) mipsolver.mipdata_->num_nodes += heur.getLocalNodes();
    // lpiterations += heur.lpiterations;
    // pseudocost = heur.pseudocost;
    return;
  }

  heurlp.removeObsoleteRows(false);
  if (!solveSubMip(heurlp.getLp(), heurlp.getLpSolver().getBasis(), fixingrate,
                   localdom.col_lower_, localdom.col_upper_,
                   500,  // std::max(50, int(0.05 *
                         // (mipsolver.mipdata_->num_leaves))),
                   200 + mipsolver.mipdata_->num_nodes / 20, 12)) {
    int64_t new_lp_iterations = lp_iterations + heur.getLocalLpIterations();
    if (new_lp_iterations + mipsolver.mipdata_->heuristic_lp_iterations >
        100000 + ((mipsolver.mipdata_->total_lp_iterations -
                   mipsolver.mipdata_->heuristic_lp_iterations -
                   mipsolver.mipdata_->sb_lp_iterations) >>
                  1)) {
      lp_iterations = new_lp_iterations;
      return;
    }

    targetdepth = heur.getCurrentDepth() / 2;
    if (targetdepth <= 1 || mipsolver.mipdata_->checkLimits()) {
      lp_iterations = new_lp_iterations;
      return;
    }
    // printf("infeasible in root node, trying with lower fixing rate\n");
    maxfixingrate = fixingrate * 0.5;
    goto retry;
  }

  lp_iterations += heur.getLocalLpIterations();
}

bool HighsPrimalHeuristics::tryRoundedPoint(const std::vector<double>& point,
                                            char source) {
  auto localdom = mipsolver.mipdata_->domain;

  HighsInt numintcols = intcols.size();
  for (HighsInt i = 0; i != numintcols; ++i) {
    HighsInt col = intcols[i];
    double intval = point[col];
    intval = std::min(localdom.col_upper_[col], intval);
    intval = std::max(localdom.col_lower_[col], intval);

    localdom.fixCol(col, intval, HighsDomain::Reason::branching());
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
      return false;
    }
    localdom.propagate();
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
      return false;
    }
  }

  if (numintcols != mipsolver.numCol()) {
    HighsLpRelaxation lprelax(mipsolver);
    lprelax.loadModel();
    lprelax.setIterationLimit(
        std::max(int64_t{10000}, 2 * mipsolver.mipdata_->firstrootlpiters));
    lprelax.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                           localdom.col_lower_.data(),
                                           localdom.col_upper_.data());

    if (numintcols / (double)mipsolver.numCol() >= 0.2)
      lprelax.getLpSolver().setOptionValue("presolve", "on");
    else
      lprelax.getLpSolver().setBasis(mipsolver.mipdata_->firstrootbasis,
                                     "HighsPrimalHeuristics::tryRoundedPoint");

    HighsLpRelaxation::Status st = lprelax.resolveLp();

    if (st == HighsLpRelaxation::Status::kInfeasible) {
      std::vector<HighsInt> inds;
      std::vector<double> vals;
      double rhs;
      if (lprelax.computeDualInfProof(mipsolver.mipdata_->domain, inds, vals,
                                      rhs)) {
        HighsCutGeneration cutGen(lprelax, mipsolver.mipdata_->cutpool);
        cutGen.generateConflict(localdom, inds, vals, rhs);
      }
      return false;
    } else if (lprelax.unscaledPrimalFeasible(st)) {
      mipsolver.mipdata_->addIncumbent(
          lprelax.getLpSolver().getSolution().col_value, lprelax.getObjective(),
          source);
      return true;
    }
  }

  return mipsolver.mipdata_->trySolution(localdom.col_lower_, source);
}

bool HighsPrimalHeuristics::linesearchRounding(
    const std::vector<double>& point1, const std::vector<double>& point2,
    char source) {
  std::vector<double> roundedpoint;

  HighsInt numintcols = intcols.size();
  roundedpoint.resize(mipsolver.numCol());

  double alpha = 0.0;
  assert(int(mipsolver.mipdata_->uplocks.size()) == mipsolver.numCol());
  assert(int(point1.size()) == mipsolver.numCol());
  assert(int(point2.size()) == mipsolver.numCol());

  while (alpha < 1.0) {
    double nextalpha = 1.0;
    bool reachedpoint2 = true;
    // printf("trying alpha = %g\n", alpha);
    for (HighsInt i = 0; i != numintcols; ++i) {
      HighsInt col = intcols[i];
      assert(col >= 0);
      assert(col < mipsolver.numCol());
      if (mipsolver.mipdata_->uplocks[col] == 0) {
        roundedpoint[col] = std::ceil(std::max(point1[col], point2[col]) -
                                      mipsolver.mipdata_->feastol);
        continue;
      }

      if (mipsolver.mipdata_->downlocks[col] == 0) {
        roundedpoint[col] = std::floor(std::min(point1[col], point2[col]) +
                                       mipsolver.mipdata_->feastol);
        continue;
      }

      double convexcomb = (1.0 - alpha) * point1[col] + alpha * point2[col];
      double intpoint2 = std::floor(point2[col] + 0.5);
      roundedpoint[col] = std::floor(convexcomb + 0.5);

      if (roundedpoint[col] == intpoint2) continue;

      reachedpoint2 = false;
      double tmpalpha = (roundedpoint[col] + 0.5 + mipsolver.mipdata_->feastol -
                         point1[col]) /
                        std::abs(point2[col] - point1[col]);
      if (tmpalpha < nextalpha && tmpalpha > alpha + 1e-2) nextalpha = tmpalpha;
    }

    if (tryRoundedPoint(roundedpoint, source)) return true;

    if (reachedpoint2) return false;

    alpha = nextalpha;
  }

  return false;
}

void HighsPrimalHeuristics::randomizedRounding(
    const std::vector<double>& relaxationsol) {
  if (int(relaxationsol.size()) != mipsolver.numCol()) return;

  auto localdom = mipsolver.mipdata_->domain;

  for (HighsInt i : intcols) {
    double intval;
    if (mipsolver.mipdata_->uplocks[i] == 0)
      intval = std::ceil(relaxationsol[i] - mipsolver.mipdata_->feastol);
    else if (mipsolver.mipdata_->downlocks[i] == 0)
      intval = std::floor(relaxationsol[i] + mipsolver.mipdata_->feastol);
    else
      intval = std::floor(relaxationsol[i] + randgen.real(0.1, 0.9));

    intval = std::min(localdom.col_upper_[i], intval);
    intval = std::max(localdom.col_lower_[i], intval);

    localdom.fixCol(i, intval, HighsDomain::Reason::branching());
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
      return;
    }
    localdom.propagate();
    if (localdom.infeasible()) {
      localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
      return;
    }
  }

  if (int(mipsolver.mipdata_->integer_cols.size()) != mipsolver.numCol()) {
    HighsLpRelaxation lprelax(mipsolver);
    lprelax.loadModel();
    lprelax.setIterationLimit(
        std::max(int64_t{10000}, 2 * mipsolver.mipdata_->firstrootlpiters));
    lprelax.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                           localdom.col_lower_.data(),
                                           localdom.col_upper_.data());
    if ((5 * intcols.size()) / mipsolver.numCol() >= 1)
      lprelax.getLpSolver().setOptionValue("presolve", "on");
    else
      lprelax.getLpSolver().setBasis(
          mipsolver.mipdata_->firstrootbasis,
          "HighsPrimalHeuristics::randomizedRounding");
    HighsLpRelaxation::Status st = lprelax.resolveLp();

    if (st == HighsLpRelaxation::Status::kInfeasible) {
      std::vector<HighsInt> inds;
      std::vector<double> vals;
      double rhs;
      if (lprelax.computeDualInfProof(mipsolver.mipdata_->domain, inds, vals,
                                      rhs)) {
        HighsCutGeneration cutGen(lprelax, mipsolver.mipdata_->cutpool);
        cutGen.generateConflict(localdom, inds, vals, rhs);
      }

    } else if (lprelax.unscaledPrimalFeasible(st))
      mipsolver.mipdata_->addIncumbent(
          lprelax.getLpSolver().getSolution().col_value, lprelax.getObjective(),
          'R');
  } else {
    mipsolver.mipdata_->trySolution(localdom.col_lower_, 'R');
  }
}

void HighsPrimalHeuristics::feasibilityPump() {
  HighsLpRelaxation lprelax(mipsolver.mipdata_->lp);
  std::unordered_set<std::vector<HighsInt>, HighsVectorHasher, HighsVectorEqual>
      referencepoints;
  std::vector<double> roundedsol;
  HighsLpRelaxation::Status status = lprelax.resolveLp();
  lp_iterations += lprelax.getNumLpIterations();

  std::vector<double> fracintcost;
  std::vector<HighsInt> fracintset;

  std::vector<HighsInt> mask(mipsolver.model_->num_col_, 1);
  std::vector<double> cost(mipsolver.model_->num_col_, 0.0);

  lprelax.getLpSolver().setOptionValue("simplex_strategy",
                                       kSimplexStrategyPrimal);
  lprelax.setObjectiveLimit();
  lprelax.getLpSolver().setOptionValue(
      "primal_simplex_bound_perturbation_multiplier", 0.0);

  lprelax.setIterationLimit(5 * mipsolver.mipdata_->avgrootlpiters);

  while (!lprelax.getFractionalIntegers().empty()) {
    const auto& lpsol = lprelax.getLpSolver().getSolution().col_value;
    roundedsol = lprelax.getLpSolver().getSolution().col_value;

    std::vector<HighsInt> referencepoint;
    referencepoint.reserve(mipsolver.mipdata_->integer_cols.size());

    auto localdom = mipsolver.mipdata_->domain;
    for (HighsInt i : mipsolver.mipdata_->integer_cols) {
      assert(mipsolver.variableType(i) == HighsVarType::kInteger);
      double intval = std::floor(roundedsol[i] + randgen.real(0.4, 0.6));
      intval = std::max(intval, localdom.col_lower_[i]);
      intval = std::min(intval, localdom.col_upper_[i]);
      roundedsol[i] = intval;
      referencepoint.push_back((HighsInt)intval);
      if (!localdom.infeasible()) {
        localdom.fixCol(i, intval, HighsDomain::Reason::branching());
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
          continue;
        }
        localdom.propagate();
        if (localdom.infeasible()) {
          localdom.conflictAnalysis(mipsolver.mipdata_->conflictPool);
          continue;
        }
      }
    }

    bool havecycle = !referencepoints.emplace(referencepoint).second;
    for (HighsInt k = 0; havecycle && k < 2; ++k) {
      for (HighsInt i = 0; i != 10; ++i) {
        HighsInt flippos =
            randgen.integer(mipsolver.mipdata_->integer_cols.size());
        HighsInt col = mipsolver.mipdata_->integer_cols[flippos];
        if (roundedsol[col] > lpsol[col])
          roundedsol[col] = (HighsInt)std::floor(lpsol[col]);
        else if (roundedsol[col] < lpsol[col])
          roundedsol[col] = (HighsInt)std::ceil(lpsol[col]);
        else if (roundedsol[col] < mipsolver.mipdata_->domain.col_upper_[col])
          roundedsol[col] = mipsolver.mipdata_->domain.col_upper_[col];
        else
          roundedsol[col] = mipsolver.mipdata_->domain.col_lower_[col];

        referencepoint[flippos] = (HighsInt)roundedsol[col];
      }
      havecycle = !referencepoints.emplace(referencepoint).second;
    }

    if (havecycle) return;

    if (linesearchRounding(lpsol, roundedsol, 'F')) return;

    if (lprelax.getNumLpIterations() >=
        1000 + mipsolver.mipdata_->avgrootlpiters * 5)
      break;

    for (HighsInt i : mipsolver.mipdata_->integer_cols) {
      assert(mipsolver.variableType(i) == HighsVarType::kInteger);

      if (mipsolver.mipdata_->uplocks[i] == 0 ||
          mipsolver.mipdata_->downlocks[i] == 0)
        cost[i] = 0.0;
      else if (lpsol[i] > roundedsol[i] - mipsolver.mipdata_->feastol)
        cost[i] = -1.0 + randgen.real(-1e-4, 1e-4);
      else
        cost[i] = 1.0 + randgen.real(-1e-4, 1e-4);
    }

    lprelax.getLpSolver().changeColsCost(mask.data(), cost.data());
    int64_t niters = -lprelax.getNumLpIterations();
    status = lprelax.resolveLp();
    niters += lprelax.getNumLpIterations();
    if (niters == 0) break;
    lp_iterations += niters;
  }

  if (lprelax.getFractionalIntegers().empty() &&
      lprelax.unscaledPrimalFeasible(status))
    mipsolver.mipdata_->addIncumbent(
        lprelax.getLpSolver().getSolution().col_value, lprelax.getObjective(),
        'F');
}

void HighsPrimalHeuristics::centralRounding() {
  if (HighsInt(mipsolver.mipdata_->analyticCenter.size()) != mipsolver.numCol())
    return;

  if (!mipsolver.mipdata_->firstlpsol.empty())
    linesearchRounding(mipsolver.mipdata_->firstlpsol,
                       mipsolver.mipdata_->analyticCenter, 'C');
  else if (!mipsolver.mipdata_->rootlpsol.empty())
    linesearchRounding(mipsolver.mipdata_->rootlpsol,
                       mipsolver.mipdata_->analyticCenter, 'C');
  else
    linesearchRounding(mipsolver.mipdata_->analyticCenter,
                       mipsolver.mipdata_->analyticCenter, 'C');
}

#if 0
void HighsPrimalHeuristics::clique() {
  HighsHashTable<HighsInt, double> entries;
  double offset = 0.0;

  HighsDomain& globaldom = mipsolver.mipdata_->domain;
  for (HighsInt j = 0; j != mipsolver.numCol(); ++j) {
    HighsInt col = j;
    double val = mipsolver.colCost(col);
    if (val == 0.0) continue;

    if (!globaldom.isBinary(col)) {
      offset += val * globaldom.col_lower_[col];
      continue;
    }

    mipsolver.mipdata_->cliquetable.resolveSubstitution(col, val, offset);
    entries[col] += val;
  }

  std::vector<double> profits;
  std::vector<HighsCliqueTable::CliqueVar> objvars;

  for (const auto& entry : entries) {
    double objprofit = -entry.value();
    if (objprofit < 0) {
      offset += objprofit;
      profits.push_back(-objprofit);
      objvars.emplace_back(entry.key(), 0);
    } else {
      profits.push_back(objprofit);
      objvars.emplace_back(entry.key(), 1);
    }
  }

  std::vector<double> solution(mipsolver.numCol());

  HighsInt nobjvars = profits.size();
  for (HighsInt i = 0; i != nobjvars; ++i) solution[objvars[i].col] = objvars[i].val;

  std::vector<std::vector<HighsCliqueTable::CliqueVar>> cliques;
  double bestviol;
  HighsInt bestviolpos;
  HighsInt numcliques;

  cliques = mipsolver.mipdata_->cliquetable.separateCliques(
      solution, mipsolver.mipdata_->domain, mipsolver.mipdata_->feastol);
  numcliques = cliques.size();
  while (numcliques != 0) {
    bestviol = 0.5;
    bestviolpos = -1;

    for (HighsInt c = 0; c != numcliques; ++c) {
      double viol = -1.0;
      for (HighsCliqueTable::CliqueVar clqvar : cliques[c])
        viol += clqvar.weight(solution);

      if (viol > bestviolpos) {
        bestviolpos = c;
        bestviol = viol;
      }
    }

    cliques = mipsolver.mipdata_->cliquetable.separateCliques(
        solution, mipsolver.mipdata_->domain, mipsolver.mipdata_->feastol);
    numcliques = cliques.size();
  }
}
#endif

void HighsPrimalHeuristics::flushStatistics() {
  mipsolver.mipdata_->heuristic_lp_iterations += lp_iterations;
  mipsolver.mipdata_->total_lp_iterations += lp_iterations;
  lp_iterations = 0;
}

void HighsPrimalHeuristics::trivial() {
  if (mipsolver.options_mip_->mip_trivial_heuristics == kHighsOffString) return;
  /*
  printf("HighsPrimalHeuristics::trivial() submip = %d(%d); numRestarts = %d;
  dimensions (%d, %d); solution_objective = %g\n", mipsolver.submip,
  int(mipsolver.submip_level), int(mipsolver.mipdata_->numRestarts),
         int(mipsolver.model_->num_col_),
         int(mipsolver.model_->num_row_),
         mipsolver.solution_objective_);
  */

  const bool try_heuristics = true;
  if (!try_heuristics) {
    for (HighsInt heuristic = 0; heuristic < kTrivialHeuristicCount;
         heuristic++)
      trivial_heuristics_statistics_.record[heuristic].not_run++;
    return;
  }

  const bool is_ip = mipsolver.mipdata_->continuous_cols.size() == 0;
  if (is_ip) {
    // Pure IPs are easy, as an LP doesn't have to be solved to assess
    // feasibility
    runTrivial();
  } else {
    Highs lc_highs;
    HighsLp lc_lp = *mipsolver.model_;
    lc_lp.integrality_.clear();
    lc_highs.setOptionValue("output_flag", false);
    lc_highs.passModel(std::move(lc_lp));
    runTrivial(&lc_highs);
  }
}

void HighsPrimalHeuristics::runTrivial(Highs* highs) {
  const std::vector<HighsInt>& integer_cols = mipsolver.mipdata_->integer_cols;
  const std::vector<HighsInt>& continuous_cols =
      mipsolver.mipdata_->continuous_cols;
  HighsInt num_integer_col = integer_cols.size();
  HighsInt num_continuous_col = continuous_cols.size();
  const bool is_ip = num_continuous_col == 0;
  if (!is_ip) assert(highs);
  const HighsInt use_num_heuristic = kTrivialHeuristicCount;

  // Now try trivial heuristics
  const std::vector<char> heuristic_source = {'z', 'u', 'l', 'p'};
  //  printf("Number of continuous columns is %d\n", int(num_continuous_col));
  //  assert(num_continuous_col == 0);

  const std::vector<double>& col_lower = mipsolver.model_->col_lower_;
  const std::vector<double>& col_upper = mipsolver.model_->col_upper_;
  const std::vector<double>& row_lower = mipsolver.model_->row_lower_;
  const std::vector<double>& row_upper = mipsolver.model_->row_upper_;
  const HighsSparseMatrix& matrix = mipsolver.model_->a_matrix_;
  // Determine the following properties, according to which some
  // trivial heuristics are duplicated or fail immediately
  bool all_integer_lower_non_positive = true;
  bool all_integer_lower_zero = true;
  bool all_integer_lower_finite = true;
  bool all_integer_upper_finite = true;
  for (HighsInt integerCol = 0; integerCol < num_integer_col; integerCol++) {
    HighsInt iCol = integer_cols[integerCol];
    if (col_lower[iCol] > 0) all_integer_lower_non_positive = false;
    if (col_lower[iCol]) all_integer_lower_zero = false;
    if (col_lower[iCol] <= -kHighsInf) all_integer_lower_finite = false;
    if (col_upper[iCol] >= kHighsInf) all_integer_upper_finite = false;
    // Only continue if one of the properties still holds
    if (!(all_integer_lower_non_positive || all_integer_lower_zero ||
          all_integer_upper_finite))
      break;
  }
  const bool all_integer_boxed =
      all_integer_lower_finite && all_integer_upper_finite;
  /*
  printf(
      "\nTrying trivial heuristics\n"
      "   all_integer_lower_non_positive = %d\n"
      "   all_integer_lower_zero = %d\n"
      "   all_integer_upper_finite = %d\n"
      "   all_integer_boxed = %d\n",
      all_integer_lower_non_positive, all_integer_lower_zero,
      all_integer_upper_finite, all_integer_boxed);
  */
  const double feasibility_tolerance =
      mipsolver.options_mip_->mip_feasibility_tolerance;
  // Loop through the trivial heuristics
  std::vector<double> solution(mipsolver.model_->num_col_);
  std::vector<HighsInt> prev_num_try;
  // Value of use_presolve is used to decide whether to use presolve
  // or the current optimal basis when an LP is solved for a
  // heuristic, and it depends on the outcome of previous heuristic LP
  // solves
  bool use_presolve = true;
  for (HighsInt heuristic = 0; heuristic < use_num_heuristic; heuristic++) {
    TrivialHeuristicRecord& record =
        trivial_heuristics_statistics_.record[heuristic];
    prev_num_try.push_back(record.not_run + record.cannot_run + record.fail +
                           record.feasible + record.improvement);
    // Heuristics are tried locally for IPs. For non-IPs, the values
    // of the integer variables are set in the solution vector, and
    // then the same code is used to test all the heuristics
    if (heuristic == 0) {
      // First heuristic is to see whether all-zero for integer
      // variables is feasible
      //
      // If there is a positive lower bound then the heuristic fails
      bool heuristic_failed = false;
      if (!all_integer_lower_non_positive) {
        record.cannot_run++;
        continue;
      }
      if (is_ip) {
        // Determine whether a zero row activity is feasible
        for (HighsInt iRow = 0; iRow < mipsolver.model_->num_row_; iRow++) {
          if (row_lower[iRow] > feasibility_tolerance ||
              row_upper[iRow] < -feasibility_tolerance) {
            heuristic_failed = true;
            break;
          }
        }
        if (heuristic_failed) {
          record.fail++;
          continue;
        }
        solution.assign(mipsolver.model_->num_col_, 0);
      } else {
        // Accumulate the integer variable assignments in the solution
        // vector
        solution.assign(num_integer_col, 0);
      }

    } else if (heuristic == 1) {
      // Second heuristic is to see whether all-upper for integer
      // variables is feasible
      //
      // If there is an infinite upper bound then the heuristic fails
      if (!all_integer_upper_finite) {
        record.cannot_run++;
        continue;
      }
      if (is_ip) {
        // Trivially feasible for columns, but what about rows?
        if (!mipsolver.mipdata_->solutionRowFeasible(col_upper)) {
          record.fail++;
          continue;
        }
        solution = col_upper;
      } else {
        // Accumulate the integer variable assignments in the solution
        // vector
        for (HighsInt integerCol = 0; integerCol < num_integer_col;
             integerCol++) {
          HighsInt iCol = integer_cols[integerCol];
          solution[integerCol] = col_upper[iCol];
        }
      }
    } else if (heuristic == 2) {
      // Third heuristic is to see whether all-lower for integer
      // variables (if distinct from all-zero) is feasible
      if (all_integer_lower_zero) {
        record.cannot_run++;
        continue;
      }
      if (is_ip) {
        // Trivially feasible for columns, but what about rows?
        if (!mipsolver.mipdata_->solutionRowFeasible(col_lower)) {
          record.fail++;
          continue;
        }
        solution = col_lower;
      } else {
        // Accumulate the integer variable assignments in the solution
        // vector
        for (HighsInt integerCol = 0; integerCol < num_integer_col;
             integerCol++) {
          HighsInt iCol = integer_cols[integerCol];
          solution[integerCol] = col_lower[iCol];
        }
      }
    } else if (heuristic == 3) {
      // Fourth heuristic is to see whether the "lock point" is feasible
      if (!all_integer_boxed) {
        record.cannot_run++;
        continue;
      }
      for (HighsInt integerCol = 0; integerCol < num_integer_col;
           integerCol++) {
        HighsInt iCol = integer_cols[integerCol];
        HighsInt num_positive_values = 0;
        HighsInt num_negative_values = 0;
        for (HighsInt iEl = matrix.start_[iCol]; iEl < matrix.start_[iCol + 1];
             iEl++) {
          if (matrix.value_[iEl] > 0)
            num_positive_values++;
          else
            num_negative_values++;
        }
        solution[integerCol] = num_positive_values > num_negative_values
                                   ? col_lower[iCol]
                                   : col_upper[iCol];
      }
      if (is_ip) {
        // Trivially feasible for columns, but what about rows?
        if (!mipsolver.mipdata_->solutionRowFeasible(solution)) {
          record.fail++;
          continue;
        }
      }
      // For non-IPs, the values of the integer variables are already
      // set in the solution vector
    }
    double check_objective = 0;
    if (!is_ip) {
      const HighsLp& presolved_lp = highs->getPresolvedLp();
      const HighsLp& lp = highs->getLp();
      // Integrality of LP should already have been cleared
      assert(!lp.isMip());
      bool heuristic_failed = false;
      highs->changeColsBounds(num_integer_col, integer_cols.data(),
                              solution.data(), solution.data());
      highs->run();
      use_presolve = 10 * (presolved_lp.num_col_ + presolved_lp.num_row_) <
                     (lp.num_col_ + lp.num_row_);
      HighsModelStatus status = highs->getModelStatus();
      heuristic_failed = status == HighsModelStatus::kInfeasible;
      if (!heuristic_failed) assert(status == HighsModelStatus::kOptimal);
      if (heuristic_failed) {
        record.fail++;
        continue;
      }
      // Heuristic has yielded an integer feasible solution, so copy
      // it into the solution vector
      solution = highs->getSolution().col_value;
      check_objective = highs->getObjectiveValue() - highs->getLp().offset_;
    }
    HighsCDouble cdouble_obj = 0.0;
    for (HighsInt iCol = 0; iCol < mipsolver.model_->num_col_; iCol++)
      cdouble_obj += mipsolver.colCost(iCol) * solution[iCol];
    double obj = double(cdouble_obj);
    if (!is_ip) {
      double dl_obj =
          std::abs(check_objective - obj) / std::max(1.0, std::abs(obj));
      const bool small_dl_obj = dl_obj < 1e-5;
      assert(small_dl_obj);
    }
    const double save_upper_bound = mipsolver.mipdata_->upper_bound;
    const bool new_incumbent = obj < mipsolver.mipdata_->upper_bound &&
                               mipsolver.mipdata_->addIncumbent(
                                   solution, obj, heuristic_source[heuristic]);
    if (new_incumbent) {
      printf(
          "Trivial heuristic %d has succeeded: objective = %g; upper bound "
          "reduced by %g (%g to %g)\n",
          int(heuristic), obj,
          mipsolver.mipdata_->upper_bound - save_upper_bound,
          mipsolver.mipdata_->upper_bound, save_upper_bound);
      record.improvement++;
    } else {
      record.feasible++;
    }
    if (use_presolve) highs->clearSolver();
  }
  for (HighsInt heuristic = 0; heuristic < use_num_heuristic; heuristic++) {
    TrivialHeuristicRecord& record =
        trivial_heuristics_statistics_.record[heuristic];
    const HighsInt num_try = record.not_run + record.cannot_run + record.fail +
                             record.feasible + record.improvement;
    const bool num_try_ok = num_try == prev_num_try[heuristic] + 1;
    if (!num_try_ok)
      printf(
          "HighsPrimalHeuristics::trivial() Heuristic %d: %d = num_try != "
          "prev_num_try[heuristic]+1 = %d\n",
          int(heuristic), int(num_try), int(prev_num_try[heuristic] + 1));
    assert(num_try_ok);
  }
}

void HighsPrimalHeuristics::initialiseLocalTrivialHeuristicsStatistics() {
  initialiseTrivialHeuristicsStatistics(this->trivial_heuristics_statistics_);
}

void HighsPrimalHeuristics::downCopyLocalTrivialHeuristicsStatistics(
    const TrivialHeuristicData& from_statistics) {
  copyTrivialHeuristicsStatistics(from_statistics,
                                  this->trivial_heuristics_statistics_);
}

void HighsPrimalHeuristics::upCopyLocalTrivialHeuristicsStatistics(
    TrivialHeuristicData& to_statistics) {
  copyTrivialHeuristicsStatistics(this->trivial_heuristics_statistics_,
                                  to_statistics);
}

// Not in class since they are called from methods in
// HighsMipSolverData::runSetup()
void initialiseTrivialHeuristicsStatistics(TrivialHeuristicData& statistics) {
  statistics.record.clear();
  TrivialHeuristicRecord record_;
  for (HighsInt heuristic = kTrivialHeuristicZero;
       heuristic < kTrivialHeuristicCount; heuristic++)
    statistics.record.push_back(record_);
}

void copyTrivialHeuristicsStatistics(
    const TrivialHeuristicData& from_statistics,
    TrivialHeuristicData& to_statistics) {
  if (!(from_statistics.record.size() == kTrivialHeuristicCount)) {
    printf("Cannot copy from statistics\n");
    return;
  }
  assert(to_statistics.record.size() == kTrivialHeuristicCount);
  /*
  HighsInt data = from_statistics.record[0].not_run +
    from_statistics.record[0].cannot_run +
    from_statistics.record[0].fail +
    from_statistics.record[0].feasible +
    from_statistics.record[0].improvement;
    printf("Can    copy from statistics: %d\n", int(data));
  */
  for (HighsInt heuristic = kTrivialHeuristicZero;
       heuristic < kTrivialHeuristicCount; heuristic++) {
    const TrivialHeuristicRecord& from_record =
        from_statistics.record[heuristic];
    TrivialHeuristicRecord& to_record = to_statistics.record[heuristic];
    to_record.not_run = from_record.not_run;
    to_record.cannot_run = from_record.cannot_run;
    to_record.fail = from_record.fail;
    to_record.feasible = from_record.feasible;
    to_record.improvement = from_record.improvement;
  }
}

void reportTrivialHeuristicsStatistics(const HighsLogOptions& log_options,
                                       const TrivialHeuristicData& statistics) {
  // If presolve is sufficient to determine the status of a sub-MIP,
  // then there is no record of trivial heuristic statistics
  if (!(statistics.record.size() == kTrivialHeuristicCount)) return;
  for (HighsInt heuristic = kTrivialHeuristicZero;
       heuristic < kTrivialHeuristicCount; heuristic++) {
    const TrivialHeuristicRecord& record = statistics.record[heuristic];
    HighsInt num_try = record.not_run + record.cannot_run + record.fail +
                       record.feasible + record.improvement;
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Heuristic %d tried %d times\n", int(heuristic), int(num_try));
  }
}
