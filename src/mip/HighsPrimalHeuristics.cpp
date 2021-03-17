/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
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
    : mipsolver(mipsolver), lp_iterations(0), randgen(mipsolver.numNonzero()) {
  successObservations = 0;
  numSuccessObservations = 0;
  infeasObservations = 0;
  numInfeasObservations = 0;
}

bool HighsPrimalHeuristics::solveSubMip(
    const HighsLp& lp, const HighsBasis& basis, double fixingRate,
    std::vector<double> colLower, std::vector<double> colUpper, int maxleaves,
    int maxnodes, int stallnodes) {
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
  submipoptions.time_limit -=
      mipsolver.timer_.read(mipsolver.timer_.solve_clock);
  submipoptions.dual_objective_value_upper_bound =
      mipsolver.mipdata_->upper_limit;
  submipoptions.presolve = "on";
  // setup solver and run it

  HighsMipSolver submipsolver(submipoptions, submip, true);
  submipsolver.rootbasis = &basis;
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

  if (submipsolver.modelstatus_ == HighsModelStatus::PRIMAL_INFEASIBLE) {
    infeasObservations += fixingRate;
    ++numInfeasObservations;
  }
  if (submipsolver.node_count_ <= 1 &&
      submipsolver.modelstatus_ == HighsModelStatus::PRIMAL_INFEASIBLE)
    return false;
  int oldNumImprovingSols = mipsolver.mipdata_->numImprovingSols;
  if (submipsolver.modelstatus_ != HighsModelStatus::PRIMAL_INFEASIBLE &&
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

  std::uniform_real_distribution<double> dist(lowFixingRate, highFixingRate);
  double fixingRate = dist(randgen);
  // if (!mipsolver.submip) printf("fixing rate: %.2f\n", 100.0 * fixingRate);
  return fixingRate;
}

void HighsPrimalHeuristics::RENS(const std::vector<double>& tmp) {
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
  double fixingrate = 0.0;
  bool stop = false;
  // heurlp.setIterationLimit(2 * mipsolver.mipdata_->maxrootlpiters);
  // printf("iterlimit: %d\n",
  //       heurlp.getLpSolver().getHighsOptions().simplex_iteration_limit);
  int targetdepth = 1;
  int nbacktracks = -1;
retry:
  HighsHashTable<int> fixedCols;
  int numGlobalFixed = 0;
  for (int i : mipsolver.mipdata_->integral_cols) {
    // skip fixed and continuous variables
    if (mipsolver.mipdata_->domain.colLower_[i] ==
        mipsolver.mipdata_->domain.colUpper_[i])
      ++numGlobalFixed;

    // count locally fixed variable
    if (localdom.isFixed(i)) fixedCols.insert(i);
  }
  int ntotal = mipsolver.mipdata_->integral_cols.size() - numGlobalFixed;
  size_t nCheckedChanges = 0;
  auto getFixingRate = [&]() {
    while (nCheckedChanges < localdom.getDomainChangeStack().size()) {
      int col = localdom.getDomainChangeStack()[nCheckedChanges++].column;
      if (mipsolver.variableType(col) == HighsVarType::CONTINUOUS) continue;

      if (localdom.isFixed(col)) fixedCols.insert(col);
    }

    return fixedCols.size() / (double)ntotal;
  };
  ++nbacktracks;
  // printf("current depth : %d   target depth : %d\n", heur.getCurrentDepth(),
  //        targetdepth);
  while (heur.getCurrentDepth() > targetdepth) {
    heur.cutoffNode();
    if (!heur.backtrack()) break;
  }
  if (!heur.hasNode()) return;

  while (true) {
    heur.evaluateNode();
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      if (mipsolver.mipdata_->domain.infeasible()) {
        lp_iterations += heur.getLocalLpIterations();
        return;
      }

      if (!heur.backtrack()) break;
      continue;
    }

    // printf("num backtracks: %d\n", nbacktracks);
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
    // printf("%d/%d fixed, fixingrate is %g\n", nfixed, ntotal, fixingrate);
    if (fixingrate >= maxfixingrate) break;
    if (stop) break;
    if (nbacktracks >= 10) break;

    int numBranched = 0;
    double stopFixingRate =
        std::min(1.0 - (1.0 - getFixingRate()) * 0.9, maxfixingrate);
    const auto& relaxationsol = heurlp.getSolution().col_value;
    for (int i : mipsolver.mipdata_->integer_cols) {
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
      auto getFixVal = [&](int col, double fracval) {
        double fixval;

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
        // make sure we do not set an infeasible domain
        fixval = std::min(localdom.colUpper_[col], fixval);
        fixval = std::max(localdom.colLower_[col], fixval);
        return fixval;
      };

      std::sort(heurlp.getFractionalIntegers().begin(),
                heurlp.getFractionalIntegers().end(),
                [&](const std::pair<int, double>& a,
                    const std::pair<int, double>& b) {
                  return std::abs(getFixVal(a.first, a.second) - a.second) <
                         std::abs(getFixVal(b.first, b.second) - b.second);
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

  // printf("fixing rate is %g\n", fixingrate);
  fixingrate = getFixingRate();
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
                   200 + int(0.05 * (mipsolver.mipdata_->num_nodes)), 5)) {
    targetdepth = heur.getCurrentDepth() / 2;
    if (targetdepth <= 1 || mipsolver.mipdata_->checkLimits()) return;
    // printf("infeasible in in root node, trying with lower fixing rate\n");
    maxfixingrate = fixingrate * 0.5;
    goto retry;
  }
}

void HighsPrimalHeuristics::RINS(const std::vector<double>& relaxationsol) {
  if (int(relaxationsol.size()) != mipsolver.numCol()) return;

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
  int nbacktracks = -1;
  int targetdepth = 1;
retry:
  HighsHashTable<int> fixedCols;
  int numGlobalFixed = 0;
  for (int i : mipsolver.mipdata_->integral_cols) {
    // skip fixed and continuous variables
    if (mipsolver.mipdata_->domain.colLower_[i] ==
        mipsolver.mipdata_->domain.colUpper_[i])
      ++numGlobalFixed;

    // count locally fixed variable
    if (localdom.isFixed(i)) fixedCols.insert(i);
  }
  int ntotal = mipsolver.mipdata_->integral_cols.size() - numGlobalFixed;
  size_t nCheckedChanges = 0;
  auto getFixingRate = [&]() {
    while (nCheckedChanges < localdom.getDomainChangeStack().size()) {
      int col = localdom.getDomainChangeStack()[nCheckedChanges++].column;
      if (mipsolver.variableType(col) == HighsVarType::CONTINUOUS) continue;

      if (localdom.isFixed(col)) fixedCols.insert(col);
    }

    return fixedCols.size() / (double)ntotal;
  };
  ++nbacktracks;
  // printf("current depth : %d   target depth : %d\n", heur.getCurrentDepth(),
  //       targetdepth);
  while (heur.getCurrentDepth() > targetdepth) {
    heur.cutoffNode();
    if (!heur.backtrack()) break;
  }
  if (!heur.hasNode()) return;

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

    std::vector<std::pair<int, double>>::iterator fixcandend;

    // partition the fractional variables to consider which ones should we fix
    // in this dive first if there is an incumbent, we dive towards the RINS
    // neighborhood
    fixcandend = std::partition(
        heurlp.getFractionalIntegers().begin(),
        heurlp.getFractionalIntegers().end(),
        [&](const std::pair<int, double>& fracvar) {
          return std::abs(relaxationsol[fracvar.first] -
                          mipsolver.mipdata_->incumbent[fracvar.first]) <=
                 mipsolver.mipdata_->feastol;
        });

    bool fixtolpsol = true;

    auto getFixVal = [&](int col, double fracval) {
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
      int numBranched = 0;

      fixingrate = getFixingRate();
      double stopFixingRate =
          std::min(maxfixingrate, 1.0 - (1.0 - getFixingRate()) * 0.9);
      const auto& currlpsol = heurlp.getSolution().col_value;
      for (int i : mipsolver.mipdata_->integer_cols) {
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
        //    "fixed %d additional cols, old fixing rate: %.2f%%, new fixing "
        //    "rate: %.2f%%\n",
        //    numBranched, fixingrate, getFixingRate());
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
    std::sort(
        heurlp.getFractionalIntegers().begin(), fixcandend,
        [&](const std::pair<int, double>& a, const std::pair<int, double>& b) {
          return std::abs(getFixVal(a.first, a.second) - a.second) <
                 std::abs(getFixVal(b.first, b.second) - b.second);
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

    // printf("%d/%d fixed, fixingrate is %g\n", nfixed, ntotal, fixingrate);
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
                   200 + int(0.05 * (mipsolver.mipdata_->num_nodes)), 5)) {
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

  int numintcols = mipsolver.mipdata_->integer_cols.size();
  for (int i = 0; i != numintcols; ++i) {
    int col = mipsolver.mipdata_->integer_cols[i];
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
    lprelax.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                           localdom.colLower_.data(),
                                           localdom.colUpper_.data());

    if (numintcols / (double)mipsolver.numCol() >= 0.2)
      lprelax.getLpSolver().setHighsOptionValue("presolve", "on");
    else
      lprelax.getLpSolver().setBasis(mipsolver.mipdata_->firstrootbasis);

    HighsLpRelaxation::Status st = lprelax.resolveLp();

    if (st == HighsLpRelaxation::Status::Infeasible) {
      std::vector<int> inds;
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

  int numintcols = mipsolver.mipdata_->integer_cols.size();
  roundedpoint.resize(mipsolver.numCol());

  double alpha = 0.0;
  assert(int(mipsolver.mipdata_->uplocks.size()) == mipsolver.numCol());
  assert(int(point1.size()) == mipsolver.numCol());
  assert(int(point2.size()) == mipsolver.numCol());

  while (alpha < 1.0) {
    double nextalpha = 1.0;
    bool reachedpoint2 = true;
    // printf("trying alpha = %g\n", alpha);
    for (int i = 0; i != numintcols; ++i) {
      int col = mipsolver.mipdata_->integer_cols[i];
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

  std::uniform_real_distribution<double> dist(0.1, 0.9);

  for (int i : mipsolver.mipdata_->integer_cols) {
    double intval;
    if (mipsolver.mipdata_->uplocks[i] == 0)
      intval = std::ceil(relaxationsol[i] - mipsolver.mipdata_->feastol);
    else if (mipsolver.mipdata_->downlocks[i] == 0)
      intval = std::floor(relaxationsol[i] + mipsolver.mipdata_->feastol);
    else
      intval = std::floor(relaxationsol[i] + dist(randgen));

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
      std::vector<int> inds;
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
  std::unordered_set<std::vector<int>, HighsVectorHasher, HighsVectorEqual>
      referencepoints;
  std::vector<double> roundedsol;
  HighsLpRelaxation::Status status = lprelax.resolveLp();
  lp_iterations += lprelax.getNumLpIterations();

  std::uniform_real_distribution<double> roundingdist(0.4, 0.6);
  std::uniform_real_distribution<double> perturbdist(-1e-4, 1e-4);
  std::vector<double> fracintcost;
  std::vector<int> fracintset;

  std::vector<int> mask(mipsolver.model_->numCol_, 1);
  std::vector<double> cost(mipsolver.model_->numCol_, 0.0);
  if (mipsolver.mipdata_->upper_limit != HIGHS_CONST_INF) {
    std::vector<int> objinds;
    std::vector<double> objval;
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.colCost(i) != 0) {
        objinds.push_back(i);
        objval.push_back(mipsolver.colCost(i));
      }
    }

    lprelax.getLpSolver().addRow(-HIGHS_CONST_INF,
                                 mipsolver.mipdata_->upper_limit,
                                 objinds.size(), objinds.data(), objval.data());
  }

  lprelax.getLpSolver().setHighsOptionValue("simplex_strategy",
                                            SIMPLEX_STRATEGY_PRIMAL);
  lprelax.getLpSolver().setHighsOptionValue(
      "primal_simplex_bound_perturbation_multiplier", 0.0);

  lprelax.setIterationLimit(mipsolver.mipdata_->maxrootlpiters);

  while (!lprelax.getFractionalIntegers().empty()) {
    const auto& lpsol = lprelax.getLpSolver().getSolution().col_value;
    roundedsol = lprelax.getLpSolver().getSolution().col_value;

    std::vector<int> referencepoint;
    referencepoint.reserve(mipsolver.mipdata_->integer_cols.size());

    auto localdom = mipsolver.mipdata_->domain;
    for (int i : mipsolver.mipdata_->integer_cols) {
      assert(mipsolver.variableType(i) == HighsVarType::INTEGER);
      double intval = std::floor(roundedsol[i] + roundingdist(randgen));
      intval = std::max(intval, localdom.colLower_[i]);
      intval = std::min(intval, localdom.colUpper_[i]);
      roundedsol[i] = intval;
      referencepoint.push_back((int)intval);
      if (!localdom.infeasible()) {
        localdom.fixCol(i, intval);
        if (localdom.infeasible()) continue;
        localdom.propagate();
        if (localdom.infeasible()) continue;
      }
    }

    bool havecycle = !referencepoints.emplace(referencepoint).second;

    while (havecycle) {
      std::uniform_int_distribution<int> dist(
          0, mipsolver.mipdata_->integer_cols.size() - 1);
      for (int i = 0; i != 10; ++i) {
        int flippos = dist(randgen);
        int col = mipsolver.mipdata_->integer_cols[flippos];
        if (roundedsol[col] > lpsol[col])
          roundedsol[col] = (int)std::floor(lpsol[col]);
        else if (roundedsol[col] < lpsol[col])
          roundedsol[col] = (int)std::ceil(lpsol[col]);
        else if (roundedsol[col] < mipsolver.mipdata_->domain.colUpper_[col])
          roundedsol[col] = mipsolver.mipdata_->domain.colUpper_[col];
        else
          roundedsol[col] = mipsolver.mipdata_->domain.colLower_[col];

        referencepoint[flippos] = (int)roundedsol[col];
      }
      havecycle = !referencepoints.emplace(referencepoint).second;
    }

    if (linesearchRounding(lpsol, roundedsol, 'F')) return;

    if (lprelax.getNumLpIterations() >=
        1000 + mipsolver.mipdata_->maxrootlpiters * 3)
      break;

    for (int i : mipsolver.mipdata_->integer_cols) {
      assert(mipsolver.variableType(i) == HighsVarType::INTEGER);

      if (lpsol[i] > roundedsol[i] - mipsolver.mipdata_->feastol)
        cost[i] = -1.0 + perturbdist(randgen);
      else
        cost[i] = 1.0 + perturbdist(randgen);
    }

    lprelax.getLpSolver().changeColsCost(mask.data(), cost.data());
    size_t oldnumiters = lprelax.getNumLpIterations();
    status = lprelax.resolveLp();
    lp_iterations += lprelax.getNumLpIterations() - oldnumiters;
  }

  if (lprelax.getFractionalIntegers().empty() &&
      lprelax.unscaledPrimalFeasible(status))
    mipsolver.mipdata_->addIncumbent(
        lprelax.getLpSolver().getSolution().col_value, lprelax.getObjective(),
        'F');
}

void HighsPrimalHeuristics::centralRounding() {
  Highs ipm;
  ipm.setHighsOptionValue("solver", "ipm");
  ipm.setHighsOptionValue("run_crossover", false);
  ipm.setHighsOptionValue("presolve", "off");
  ipm.setHighsOptionValue("output_flag", false);
  HighsLp lpmodel(
      *mipsolver.model_);  // mipsolver.mipdata_->lp.getLpSolver().getLp());
  // lpmodel.colLower_ = mipsolver.mipdata_->domain.colLower_;
  // lpmodel.colUpper_ = mipsolver.mipdata_->domain.colUpper_;
  lpmodel.colCost_.assign(lpmodel.numCol_, 0.0);
  ipm.passModel(std::move(lpmodel));

  if (mipsolver.mipdata_->upper_limit != HIGHS_CONST_INF) {
    std::vector<int> objinds;
    std::vector<double> objval;
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.colCost(i) != 0) {
        objinds.push_back(i);
        objval.push_back(mipsolver.colCost(i));
      }
    }

    ipm.addRow(-HIGHS_CONST_INF, mipsolver.mipdata_->upper_limit,
               objinds.size(), objinds.data(), objval.data());
  }
  ipm.run();
  const std::vector<double>& sol = ipm.getSolution().col_value;
  if (int(sol.size()) != mipsolver.numCol()) return;
  if (ipm.getModelStatus() == HighsModelStatus::OPTIMAL) {
    int nfixed = 0;
    int nintfixed = 0;
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.mipdata_->domain.colLower_[i] ==
          mipsolver.mipdata_->domain.colUpper_[i])
        continue;
      if (sol[i] <=
          mipsolver.model_->colLower_[i] + mipsolver.mipdata_->feastol) {
        mipsolver.mipdata_->domain.changeBound(HighsBoundType::Upper, i,
                                               mipsolver.model_->colLower_[i]);
        if (mipsolver.mipdata_->domain.infeasible()) return;
        ++nfixed;
        if (mipsolver.variableType(i) == HighsVarType::INTEGER) ++nintfixed;
      } else if (sol[i] >=
                 mipsolver.model_->colUpper_[i] - mipsolver.mipdata_->feastol) {
        mipsolver.mipdata_->domain.changeBound(HighsBoundType::Lower, i,
                                               mipsolver.model_->colUpper_[i]);
        if (mipsolver.mipdata_->domain.infeasible()) return;
        ++nfixed;
        if (mipsolver.variableType(i) == HighsVarType::INTEGER) ++nintfixed;
      }
    }
    if (nfixed > 0)
      highsLogUser(mipsolver.options_mip_->log_options, HighsLogType::INFO,
                   "Fixing %d columns (%d integers) sitting at bound at "
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

void HighsPrimalHeuristics::clique() {
  HighsHashTable<int, double> entries;
  double offset = 0.0;

  HighsDomain& globaldom = mipsolver.mipdata_->domain;
  for (int j = 0; j != mipsolver.numCol(); ++j) {
    int col = j;
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

  int nobjvars = profits.size();
  for (int i = 0; i != nobjvars; ++i) solution[objvars[i].col] = objvars[i].val;

  std::vector<std::vector<HighsCliqueTable::CliqueVar>> cliques;
  double bestviol;
  int bestviolpos;
  int numcliques;

  cliques = mipsolver.mipdata_->cliquetable.separateCliques(
      solution, mipsolver.mipdata_->domain, mipsolver.mipdata_->feastol);
  numcliques = cliques.size();
  while (numcliques != 0) {
    bestviol = 0.5;
    bestviolpos = -1;

    for (int c = 0; c != numcliques; ++c) {
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

void HighsPrimalHeuristics::flushStatistics() {
  mipsolver.mipdata_->heuristic_lp_iterations += lp_iterations;
  mipsolver.mipdata_->total_lp_iterations += lp_iterations;
  lp_iterations = 0;
}
