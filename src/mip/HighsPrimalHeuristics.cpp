/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsPrimalHeuristics.h"

#include <numeric>
#include <unordered_set>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLpUtils.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsDomainChange.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolverData.h"
#include "util/HighsHash.h"

HighsPrimalHeuristics::HighsPrimalHeuristics(HighsMipSolver& mipsolver)
    : mipsolver(mipsolver),
      lp_iterations(0),
      randgen(mipsolver.options_mip_->highs_random_seed) {
  successObservations = 0;
  numSuccessObservations = 0;
  infeasObservations = 0;
  numInfeasObservations = 0;
}

void HighsPrimalHeuristics::setupIntCols() {
  intcols = mipsolver.mipdata_->integer_cols;

  std::sort(intcols.begin(), intcols.end(), [&](HighsInt c1, HighsInt c2) {
    HighsInt uplocks1 = mipsolver.mipdata_->uplocks[c1];
    HighsInt downlocks1 = mipsolver.mipdata_->downlocks[c1];

    HighsInt cliqueImplicsUp1 =
        mipsolver.mipdata_->cliquetable.getNumImplications(c1, 1);
    HighsInt cliqueImplicsDown1 =
        mipsolver.mipdata_->cliquetable.getNumImplications(c1, 0);

    HighsInt uplocks2 = mipsolver.mipdata_->uplocks[c2];
    HighsInt downlocks2 = mipsolver.mipdata_->downlocks[c2];

    HighsInt cliqueImplicsUp2 =
        mipsolver.mipdata_->cliquetable.getNumImplications(c2, 1);
    HighsInt cliqueImplicsDown2 =
        mipsolver.mipdata_->cliquetable.getNumImplications(c2, 0);

    return std::make_tuple(uplocks1 * downlocks1,
                           cliqueImplicsUp1 * cliqueImplicsDown1,
                           HighsHashHelpers::hash(uint64_t(c1)), c1) >
           std::make_tuple(uplocks2 * downlocks2,
                           cliqueImplicsUp2 * cliqueImplicsDown2,
                           HighsHashHelpers::hash(uint64_t(c2)), c2);
  });
}

bool HighsPrimalHeuristics::solveSubMip(
    const HighsLp& lp, const HighsBasis& basis, double fixingRate,
    std::vector<double> colLower, std::vector<double> colUpper,
    HighsInt maxleaves, HighsInt maxnodes, HighsInt stallnodes) {
  HighsOptions submipoptions = *mipsolver.options_mip_;
  HighsLp submip = lp;

  // set bounds and restore integrality of the lp relaxation copy
  submip.colLower_ = std::move(colLower);
  submip.colUpper_ = std::move(colUpper);
  submip.integrality_ = mipsolver.model_->integrality_;
  submip.offset_ = 0;

  // set limits
  submipoptions.mip_max_leaves = maxleaves;
  submipoptions.output_flag = false;
  submipoptions.mip_max_nodes = maxnodes;
  submipoptions.mip_max_stall_nodes = stallnodes;
  submipoptions.mip_pscost_minreliable = 0;
  submipoptions.time_limit -=
      mipsolver.timer_.read(mipsolver.timer_.solve_clock);
  submipoptions.dual_objective_value_upper_bound =
      mipsolver.mipdata_->upper_limit;
  submipoptions.presolve = "on";
  // setup solver and run it

  HighsMipSolver submipsolver(submipoptions, submip, true);
  submipsolver.rootbasis = &basis;
  HighsPseudocostInitialization pscostinit(mipsolver.mipdata_->pseudocost, 1);
  submipsolver.pscostinit = &pscostinit;
  submipsolver.clqtableinit = &mipsolver.mipdata_->cliquetable;
  submipsolver.implicinit = &mipsolver.mipdata_->implications;
  submipsolver.run();
  if (submipsolver.mipdata_) {
    double adjustmentfactor =
        submipsolver.numNonzero() / (double)mipsolver.numNonzero();
    size_t adjusted_lp_iterations =
        (size_t)(adjustmentfactor * adjustmentfactor *
                 submipsolver.mipdata_->total_lp_iterations);
    lp_iterations += adjusted_lp_iterations;

    if (mipsolver.submip)
      mipsolver.mipdata_->num_nodes += std::max(
          size_t{1}, size_t(adjustmentfactor * submipsolver.node_count_));
  }

  if (submipsolver.modelstatus_ == HighsModelStatus::kPrimalInfeasible) {
    infeasObservations += fixingRate;
    ++numInfeasObservations;
  }
  if (submipsolver.node_count_ <= 1 &&
      submipsolver.modelstatus_ == HighsModelStatus::kPrimalInfeasible)
    return false;
  HighsInt oldNumImprovingSols = mipsolver.mipdata_->numImprovingSols;
  if (submipsolver.modelstatus_ != HighsModelStatus::kPrimalInfeasible &&
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

void HighsPrimalHeuristics::RENS(const std::vector<double>& tmp) {
  HighsSearch heur(mipsolver, mipsolver.mipdata_->pseudocost);
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
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                        localdom.colLower_.data(),
                                        localdom.colUpper_.data());
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
  std::shared_ptr<const HighsBasis> basis;
retry:
  ++nbacktracks;
  // printf("current depth : %" HIGHSINT_FORMAT
  //        "   target depth : %" HIGHSINT_FORMAT "\n",
  //        heur.getCurrentDepth(), targetdepth);
  if (heur.getCurrentDepth() > targetdepth) {
    if (!heur.backtrackUntilDepth(targetdepth)) return;
  }
  HighsHashTable<HighsInt> fixedCols;
  HighsInt numGlobalFixed = 0;
  for (HighsInt i : mipsolver.mipdata_->integral_cols) {
    // skip fixed and continuous variables
    if (mipsolver.mipdata_->domain.colLower_[i] ==
        mipsolver.mipdata_->domain.colUpper_[i])
      ++numGlobalFixed;

    // count locally fixed variable
    if (localdom.isFixed(i)) fixedCols.insert(i);
  }
  HighsInt ntotal = mipsolver.mipdata_->integral_cols.size() - numGlobalFixed;
  size_t nCheckedChanges = 0;
  auto getFixingRate = [&]() {
    while (nCheckedChanges < localdom.getDomainChangeStack().size()) {
      HighsInt col = localdom.getDomainChangeStack()[nCheckedChanges++].column;
      if (mipsolver.variableType(col) == HighsVarType::kContinuous) continue;

      if (localdom.isFixed(col)) fixedCols.insert(col);
    }

    return fixedCols.size() / (double)ntotal;
  };
  // printf("fixingrate before loop is %g\n", fixingrate);
  assert(heur.hasNode());
  if (basis) {
    heurlp.setStoredBasis(basis);
    heurlp.recoverBasis();
  }
  while (true) {
    // printf("evaluating node\n");
    heur.evaluateNode();
    if (!basis) {
      heurlp.storeBasis();
      basis = heurlp.getStoredBasis();
    }
    // printf("done evaluating node\n");
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      if (mipsolver.mipdata_->domain.infeasible()) {
        lp_iterations += heur.getLocalLpIterations();
        return;
      }

      if (!heur.backtrack()) break;
      continue;
    }

    // printf("num backtracks: %" HIGHSINT_FORMAT "\n", nbacktracks);
    // if we estimate that there is no improving solution in this subtree, we
    // stop fixing variables but still backtrack to a node that has a good
    // estimate and is not pruned as the stop flag is checked after the
    // acktracking
    if (heur.getCurrentEstimate() > mipsolver.mipdata_->upper_limit) {
      heur.cutoffNode();
      ++nbacktracks;
      // printf("node cutoff due to bad estimate\n");
      if (!heur.backtrack()) break;
      stop = true;
      continue;
    }

    fixingrate = getFixingRate();
    // printf("after evaluating node current fixingrate is %g\n", fixingrate);
    if (fixingrate >= maxfixingrate) break;
    if (stop) break;
    if (nbacktracks >= 10) break;

    HighsInt numBranched = 0;
    double stopFixingRate =
        std::min(1.0 - (1.0 - getFixingRate()) * 0.9, maxfixingrate);
    const auto& relaxationsol = heurlp.getSolution().col_value;
    for (HighsInt i : intcols) {
      if (localdom.colLower_[i] == localdom.colUpper_[i]) continue;

      double downval =
          std::floor(relaxationsol[i] + mipsolver.mipdata_->feastol);
      double upval = std::ceil(relaxationsol[i] - mipsolver.mipdata_->feastol);

      downval = std::min(downval, localdom.colUpper_[i]);
      upval = std::max(upval, localdom.colLower_[i]);
      if (localdom.colLower_[i] < downval) {
        heur.branchUpwards(i, downval, downval - 0.5);
        ++numBranched;
      }
      if (localdom.colUpper_[i] > upval) {
        heur.branchDownwards(i, upval, upval + 0.5);
        ++numBranched;
      }

      localdom.propagate();
      if (localdom.infeasible()) break;

      if (getFixingRate() >= stopFixingRate) break;
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
        if (mipsolver.model_->colCost_[col] > 0.0)
          fixval = std::ceil(fracval);
        else if (mipsolver.model_->colCost_[col] < 0.0)
          fixval = std::floor(fracval);
        else
          fixval = std::floor(fracval + 0.5);
        // make sure we do not set an infeasible domain
        fixval = std::min(localdom.colUpper_[col], fixval);
        fixval = std::max(localdom.colLower_[col], fixval);
        return fixval;
      };

      std::sort(heurlp.getFractionalIntegers().begin(),
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

        if (localdom.colLower_[fracint.first] < fixval) {
          heur.branchUpwards(fracint.first, fixval, fracint.second);
          if (localdom.infeasible()) break;
        }

        if (localdom.colUpper_[fracint.first] > fixval) {
          heur.branchDownwards(fracint.first, fixval, fracint.second);
          if (localdom.infeasible()) break;
        }

        localdom.propagate();
        if (localdom.infeasible()) break;

        fixingrate = getFixingRate();
        if (fixingrate >= maxfixingrate) break;

        change += std::abs(fixval - fracint.second);
        if (change >= 0.5) break;
      }
    }
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

  fixingrate = getFixingRate();
  // printf("fixing rate is %g\n", fixingrate);
  if (fixingrate < 0.1 ||
      (mipsolver.submip && mipsolver.mipdata_->numImprovingSols != 0)) {
    // heur.childselrule = ChildSelectionRule::BestCost;
    heur.setMinReliable(0);
    heur.solveDepthFirst(10);
    lp_iterations += heur.getLocalLpIterations();
    // lpiterations += heur.lpiterations;
    // pseudocost = heur.pseudocost;
    return;
  }

  heurlp.removeObsoleteRows(false);
  if (!solveSubMip(heurlp.getLp(), heurlp.getLpSolver().getBasis(),
                   getFixingRate(), localdom.colLower_, localdom.colUpper_,
                   500,  // std::max(50, int(0.05 *
                         // (mipsolver.mipdata_->num_leaves))),
                   200 + int(0.05 * (mipsolver.mipdata_->num_nodes)), 12)) {
    targetdepth = heur.getCurrentDepth() / 2;
    if (targetdepth <= 1 || mipsolver.mipdata_->checkLimits()) return;
    maxfixingrate = fixingrate * 0.5;
    // printf("infeasible in in root node, trying with lower fixing rate %g\n",
    //        maxfixingrate);
    goto retry;
  }
}

void HighsPrimalHeuristics::RINS(const std::vector<double>& relaxationsol) {
  if (int(relaxationsol.size()) != mipsolver.numCol()) return;

  intcols.erase(std::remove_if(intcols.begin(), intcols.end(),
                               [&](HighsInt i) {
                                 return mipsolver.mipdata_->domain.isFixed(i);
                               }),
                intcols.end());

  HighsSearch heur(mipsolver, mipsolver.mipdata_->pseudocost);
  HighsDomain& localdom = heur.getLocalDomain();
  heur.setHeuristic(true);

  HighsLpRelaxation heurlp(mipsolver.mipdata_->lp);
  // only use the global upper limit as LP limit so that dual proofs are valid
  heurlp.setObjectiveLimit(mipsolver.mipdata_->upper_limit);
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                        localdom.colLower_.data(),
                                        localdom.colUpper_.data());
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
  std::shared_ptr<const HighsBasis> basis;
retry:
  ++nbacktracks;
  // printf("current depth : %" HIGHSINT_FORMAT "   target depth : %"
  // HIGHSINT_FORMAT "\n", heur.getCurrentDepth(),
  //       targetdepth);
  if (heur.getCurrentDepth() > targetdepth) {
    if (!heur.backtrackUntilDepth(targetdepth)) return;
  }
  HighsHashTable<HighsInt> fixedCols;
  HighsInt numGlobalFixed = 0;
  for (HighsInt i : mipsolver.mipdata_->integral_cols) {
    // skip fixed and continuous variables
    if (mipsolver.mipdata_->domain.colLower_[i] ==
        mipsolver.mipdata_->domain.colUpper_[i])
      ++numGlobalFixed;

    // count locally fixed variable
    if (localdom.isFixed(i)) fixedCols.insert(i);
  }
  HighsInt ntotal = mipsolver.mipdata_->integral_cols.size() - numGlobalFixed;
  size_t nCheckedChanges = 0;
  auto getFixingRate = [&]() {
    while (nCheckedChanges < localdom.getDomainChangeStack().size()) {
      HighsInt col = localdom.getDomainChangeStack()[nCheckedChanges++].column;
      if (mipsolver.variableType(col) == HighsVarType::kContinuous) continue;

      if (localdom.isFixed(col)) fixedCols.insert(col);
    }

    return fixedCols.size() / (double)ntotal;
  };
  assert(heur.hasNode());
  if (basis) {
    heurlp.setStoredBasis(basis);
    heurlp.recoverBasis();
  }

  while (true) {
    heur.evaluateNode();
    if (!basis) {
      heurlp.storeBasis();
      basis = heurlp.getStoredBasis();
    }
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      // printf("backtrack1\n");
      if (mipsolver.mipdata_->domain.infeasible()) {
        lp_iterations += heur.getLocalLpIterations();
        return;
      }

      if (!heur.backtrack()) break;
      continue;
    }

    // if we estimate that there is no improving solution in this subtree, we
    // stop fixing variables but still backtrack to a node that has a good
    // estimate and is not pruned as the stop flag is checked after the
    // acktracking
    if (heur.getCurrentEstimate() > mipsolver.mipdata_->upper_limit) {
      heur.cutoffNode();
      ++nbacktracks;
      // printf("backtrack2\n");
      if (!heur.backtrack()) break;
      stop = true;
      continue;
    }

    fixingrate = getFixingRate();

    if (stop) break;
    if (fixingrate >= maxfixingrate) break;
    if (nbacktracks >= 10) break;

    std::vector<std::pair<HighsInt, double>>::iterator fixcandend;

    // partition the fractional variables to consider which ones should we fix
    // in this dive first if there is an incumbent, we dive towards the RINS
    // neighborhood
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
        // RINS neighborhood (with extension)
        fixval = std::floor(relaxationsol[col] + 0.5);
      } else {
        // reinforce direction of this solution away from root
        // solution if the change is at least 0.4
        // otherwise take the direction where the objective gets worse
        // if objcetive is zero round to nearest integer
        double rootchange = fracval - mipsolver.mipdata_->rootlpsol[col];
        if (rootchange >= 0.4)
          fixval = std::ceil(fracval);
        else if (rootchange <= -0.4)
          fixval = std::floor(fracval);
        if (mipsolver.model_->colCost_[col] > 0.0)
          fixval = std::ceil(fracval);
        else if (mipsolver.model_->colCost_[col] < 0.0)
          fixval = std::floor(fracval);
        else
          fixval = std::floor(fracval + 0.5);
      }
      // make sure we do not set an infeasible domain
      fixval = std::min(localdom.colUpper_[col], fixval);
      fixval = std::max(localdom.colLower_[col], fixval);
      return fixval;
    };

    // no candidates left to fix for getting to the neighborhood, therefore we
    // switch to a different diving strategy until the minimal fixing rate is
    // reached
    if (heurlp.getFractionalIntegers().begin() == fixcandend) {
      HighsInt numBranched = 0;

      fixingrate = getFixingRate();
      double stopFixingRate =
          std::min(maxfixingrate, 1.0 - (1.0 - getFixingRate()) * 0.9);
      const auto& currlpsol = heurlp.getSolution().col_value;
      for (HighsInt i : intcols) {
        if (localdom.colLower_[i] == localdom.colUpper_[i]) continue;

        if (std::abs(currlpsol[i] - mipsolver.mipdata_->incumbent[i]) <=
            mipsolver.mipdata_->feastol) {
          double fixval = std::round(currlpsol[i]);
          if (localdom.colLower_[i] < fixval) {
            heur.branchUpwards(i, fixval, fixval - 0.5);
            ++numBranched;
          }
          if (localdom.colUpper_[i] > fixval) {
            heur.branchDownwards(i, fixval, fixval + 0.5);
            ++numBranched;
          }

          localdom.propagate();
          if (localdom.infeasible()) break;
          if (getFixingRate() >= stopFixingRate) break;
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
        break;  // if the RINS neigborhood achieved a high enough fixing rate
                // by itself we stop here
      fixcandend = heurlp.getFractionalIntegers().end();
      // now sort the variables by their distance towards the value they will
      // be fixed to
      fixtolpsol = false;
    }

    // now sort the variables by their distance towards the value they will be
    // fixed to
    std::sort(heurlp.getFractionalIntegers().begin(), fixcandend,
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

      if (localdom.colLower_[fracint.first] < fixval) {
        heur.branchUpwards(fracint.first, fixval, fracint.second);
        if (localdom.infeasible()) break;
      }

      if (localdom.colUpper_[fracint.first] > fixval) {
        heur.branchDownwards(fracint.first, fixval, fracint.second);
        if (localdom.infeasible()) break;
      }

      localdom.propagate();
      if (localdom.infeasible()) break;

      fixingrate = getFixingRate();
      if (fixingrate >= maxfixingrate) break;

      change += std::abs(fixval - fracint.second);
      if (change >= 0.5) break;
    }

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
  fixingrate = getFixingRate();
  if (fixingrate < 0.1) {
    // heur.childselrule = ChildSelectionRule::BestCost;
    heur.setMinReliable(0);
    heur.solveDepthFirst(10);
    lp_iterations += heur.getLocalLpIterations();
    // lpiterations += heur.lpiterations;
    // pseudocost = heur.pseudocost;
    return;
  }

  heurlp.removeObsoleteRows(false);
  if (!solveSubMip(heurlp.getLp(), heurlp.getLpSolver().getBasis(),
                   getFixingRate(), localdom.colLower_, localdom.colUpper_,
                   500,  // std::max(50, int(0.05 *
                         // (mipsolver.mipdata_->num_leaves))),
                   200 + int(0.05 * (mipsolver.mipdata_->num_nodes)), 12)) {
    targetdepth = heur.getCurrentDepth() / 2;
    if (targetdepth <= 1 || mipsolver.mipdata_->checkLimits()) return;
    // printf("infeasible in in root node, trying with lower fixing rate\n");
    maxfixingrate = fixingrate * 0.5;
    goto retry;
  }
}

bool HighsPrimalHeuristics::tryRoundedPoint(const std::vector<double>& point,
                                            char source) {
  auto localdom = mipsolver.mipdata_->domain;

  HighsInt numintcols = intcols.size();
  for (HighsInt i = 0; i != numintcols; ++i) {
    HighsInt col = intcols[i];
    double intval = point[col];
    intval = std::min(localdom.colUpper_[col], intval);
    intval = std::max(localdom.colLower_[col], intval);

    localdom.fixCol(col, intval);
    if (localdom.infeasible()) return false;
    localdom.propagate();
    if (localdom.infeasible()) return false;
  }

  if (numintcols != mipsolver.numCol()) {
    HighsLpRelaxation lprelax(mipsolver);
    lprelax.loadModel();
    lprelax.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                           localdom.colLower_.data(),
                                           localdom.colUpper_.data());

    if (numintcols / (double)mipsolver.numCol() >= 0.2)
      lprelax.getLpSolver().setOptionValue("presolve", "on");
    else
      lprelax.getLpSolver().setBasis(mipsolver.mipdata_->firstrootbasis);

    HighsLpRelaxation::Status st = lprelax.resolveLp();

    if (st == HighsLpRelaxation::Status::Infeasible) {
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

  return mipsolver.mipdata_->trySolution(localdom.colLower_, source);
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

    intval = std::min(localdom.colUpper_[i], intval);
    intval = std::max(localdom.colLower_[i], intval);

    localdom.fixCol(i, intval, HighsDomain::Reason::branching());
    if (localdom.infeasible()) return;
    localdom.propagate();
    if (localdom.infeasible()) return;
  }

  if (int(mipsolver.mipdata_->integer_cols.size()) != mipsolver.numCol()) {
    HighsLpRelaxation lprelax(mipsolver.mipdata_->lp);
    lprelax.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                           localdom.colLower_.data(),
                                           localdom.colUpper_.data());
    HighsLpRelaxation::Status st = lprelax.resolveLp();

    if (st == HighsLpRelaxation::Status::Infeasible) {
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
    mipsolver.mipdata_->trySolution(localdom.colLower_, 'R');
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

  std::vector<HighsInt> mask(mipsolver.model_->numCol_, 1);
  std::vector<double> cost(mipsolver.model_->numCol_, 0.0);
  if (mipsolver.mipdata_->upper_limit != kHighsInf) {
    std::vector<HighsInt> objinds;
    std::vector<double> objval;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.colCost(i) != 0) {
        objinds.push_back(i);
        objval.push_back(mipsolver.colCost(i));
      }
    }

    lprelax.getLpSolver().addRow(-kHighsInf, mipsolver.mipdata_->upper_limit,
                                 objinds.size(), objinds.data(), objval.data());
  }

  lprelax.getLpSolver().setOptionValue("simplex_strategy",
                                       SIMPLEX_STRATEGY_PRIMAL);
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
      intval = std::max(intval, localdom.colLower_[i]);
      intval = std::min(intval, localdom.colUpper_[i]);
      roundedsol[i] = intval;
      referencepoint.push_back((HighsInt)intval);
      if (!localdom.infeasible()) {
        localdom.fixCol(i, intval);
        if (localdom.infeasible()) continue;
        localdom.propagate();
        if (localdom.infeasible()) continue;
      }
    }

    bool havecycle = !referencepoints.emplace(referencepoint).second;

    while (havecycle) {
      for (HighsInt i = 0; i != 10; ++i) {
        HighsInt flippos =
            randgen.integer(mipsolver.mipdata_->integer_cols.size());
        HighsInt col = mipsolver.mipdata_->integer_cols[flippos];
        if (roundedsol[col] > lpsol[col])
          roundedsol[col] = (HighsInt)std::floor(lpsol[col]);
        else if (roundedsol[col] < lpsol[col])
          roundedsol[col] = (HighsInt)std::ceil(lpsol[col]);
        else if (roundedsol[col] < mipsolver.mipdata_->domain.colUpper_[col])
          roundedsol[col] = mipsolver.mipdata_->domain.colUpper_[col];
        else
          roundedsol[col] = mipsolver.mipdata_->domain.colLower_[col];

        referencepoint[flippos] = (HighsInt)roundedsol[col];
      }
      havecycle = !referencepoints.emplace(referencepoint).second;
    }

    if (linesearchRounding(lpsol, roundedsol, 'F')) return;

    if (lprelax.getNumLpIterations() >=
        1000 + mipsolver.mipdata_->avgrootlpiters * 5)
      break;

    for (HighsInt i : mipsolver.mipdata_->integer_cols) {
      assert(mipsolver.variableType(i) == HighsVarType::kInteger);

      if (lpsol[i] > roundedsol[i] - mipsolver.mipdata_->feastol)
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
  Highs ipm;
  ipm.setOptionValue("solver", "ipm");
  ipm.setOptionValue("run_crossover", false);
  ipm.setOptionValue("presolve", "off");
  ipm.setOptionValue("output_flag", false);
  HighsLp lpmodel(
      *mipsolver.model_);  // mipsolver.mipdata_->lp.getLpSolver().getLp());
  // lpmodel.colLower_ = mipsolver.mipdata_->domain.colLower_;
  // lpmodel.colUpper_ = mipsolver.mipdata_->domain.colUpper_;
  lpmodel.colCost_.assign(lpmodel.numCol_, 0.0);
  ipm.passModel(std::move(lpmodel));

  if (mipsolver.mipdata_->upper_limit != kHighsInf) {
    std::vector<HighsInt> objinds;
    std::vector<double> objval;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.colCost(i) != 0) {
        objinds.push_back(i);
        objval.push_back(mipsolver.colCost(i));
      }
    }

    ipm.addRow(-kHighsInf, mipsolver.mipdata_->upper_limit, objinds.size(),
               objinds.data(), objval.data());
  }
  ipm.run();
  const std::vector<double>& sol = ipm.getSolution().col_value;
  if (int(sol.size()) != mipsolver.numCol()) return;
  if (ipm.getModelStatus() == HighsModelStatus::kOptimal) {
    HighsInt nfixed = 0;
    HighsInt nintfixed = 0;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.mipdata_->domain.colLower_[i] ==
          mipsolver.mipdata_->domain.colUpper_[i])
        continue;
      if (sol[i] <=
          mipsolver.model_->colLower_[i] + mipsolver.mipdata_->feastol) {
        mipsolver.mipdata_->domain.changeBound(HighsBoundType::kUpper, i,
                                               mipsolver.model_->colLower_[i]);
        if (mipsolver.mipdata_->domain.infeasible()) return;
        ++nfixed;
        if (mipsolver.variableType(i) == HighsVarType::kInteger) ++nintfixed;
      } else if (sol[i] >=
                 mipsolver.model_->colUpper_[i] - mipsolver.mipdata_->feastol) {
        mipsolver.mipdata_->domain.changeBound(HighsBoundType::kLower, i,
                                               mipsolver.model_->colUpper_[i]);
        if (mipsolver.mipdata_->domain.infeasible()) return;
        ++nfixed;
        if (mipsolver.variableType(i) == HighsVarType::kInteger) ++nintfixed;
      }
    }
    if (nfixed > 0)
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::kInfo,
                   "Fixing %" HIGHSINT_FORMAT " columns (%" HIGHSINT_FORMAT
                   " integers) sitting at bound at "
                   "analytic center\n",
                   nfixed, nintfixed);
    mipsolver.mipdata_->domain.propagate();
    if (mipsolver.mipdata_->domain.infeasible()) return;
  }

  if (!mipsolver.mipdata_->firstlpsol.empty())
    linesearchRounding(mipsolver.mipdata_->firstlpsol, sol, 'C');
  else if (!mipsolver.mipdata_->rootlpsol.empty())
    linesearchRounding(mipsolver.mipdata_->rootlpsol, sol, 'C');
  else
    linesearchRounding(sol, sol, 'C');
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
      offset += val * globaldom.colLower_[col];
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
