/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsSearch.h"

#include <numeric>

#include "mip/HighsMipSolverData.h"

HighsSearch::HighsSearch(HighsMipSolver& mipsolver,
                         const HighsPseudocost& pseudocost)
    : mipsolver(mipsolver),
      lp(nullptr),
      localdom(mipsolver.mipdata_->domain.createChildDomain()),
      pseudocost(pseudocost) {
  nnodes = 0;
  treeweight = 0.0;
  depthoffset = 0;
  lpiterations = 0;
  heurlpiterations = 0;
  sblpiterations = 0;
  upper_limit = HIGHS_CONST_INF;
  inheuristic = false;
  inbranching = false;
  childselrule = ChildSelectionRule::Obj;
  this->localdom.setDomainChangeStack(std::vector<HighsDomainChange>());
}

double HighsSearch::checkSol(const std::vector<double>& sol,
                             bool& integerfeasible) const {
  HighsCDouble objval = 0.0;
  integerfeasible = true;
  for (int i = 0; i != mipsolver.numCol(); ++i) {
    objval += sol[i] * mipsolver.colCost(i);
    assert(std::isfinite(sol[i]));

    if (!integerfeasible || mipsolver.variableType(i) != HighsVarType::INTEGER)
      continue;

    double intval = std::floor(sol[i] + 0.5);
    if (std::abs(sol[i] - intval) > mipsolver.mipdata_->feastol) {
      integerfeasible = false;
    }
  }

  return double(objval);
}

#if 0
void HighsSearch::roundingHeuristic() {
  assert(localdom.getChangedCols().empty());
  assert(lp->unscaledPrimalFeasible(lp->getStatus()) &&
         !lp->getFractionalIntegers().empty());
  const std::vector<double>& lpsol = lp->getLpSolver().getSolution().col_value;

  HighsDomain& domain = localdom;

#ifdef HIGHS_EXTRA_DEBUG
  auto oldlower = domain.colLower_;
  auto oldupper = domain.colUpper_;
  auto oldlowerlp = lp->getLpSolver().getLp().colLower_;
  auto oldupperlp = lp->getLpSolver().getLp().colUpper_;

  for (int i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;
    assert(oldlowerlp[i] == oldlower[i]);
    assert(oldupperlp[i] == oldupper[i]);
  }
#endif

  for (int i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;

    double intval = std::floor(lpsol[i] + 0.5);
    double fixval;
    if (std::abs(lpsol[i] - intval) <= mipsolver.mipdata_->feastol) {
      fixval = intval;
    } else {
      double frac = lpsol[i] - std::floor(lpsol[i]);
      if (random.fraction() > frac)
        fixval = std::floor(lpsol[i]);
      else
        fixval = std::ceil(lpsol[i]);
    }
    bool changed = false;

    if (domain.colLower_[i] < fixval - mipsolver.mipdata_->feastol) {
      domain.changeBound(HighsBoundType::Lower, i, fixval,
                         domain.getChangedCols().empty() ? -1 : -2);
      changed = true;
    }

    if (domain.colUpper_[i] > fixval + mipsolver.mipdata_->feastol) {
      domain.changeBound(HighsBoundType::Upper, i, fixval,
                         domain.getChangedCols().empty() ? -1 : -2);
      changed = true;
    }
    if (changed) {
      domain.propagate();
      if (domain.infeasible()) {
        domain.backtrack();
        domain.clearChangedCols();
        return;
      }
    }
  }

  assert(!domain.getChangedCols().empty());

  // printf("solving LP for rounded solution");
  lp->storeBasis();
  bool success = false;
  lp->flushDomain(domain);
  size_t oldnumiter = lp->getNumLpIterations();
  lp->resolveLp();
  lpiterations += lp->getNumLpIterations() - oldnumiter;
  if (lp->integerFeasible() && lp->getObjective() < upper_bound) {
    success = true;
    addIncumbent(lp->getLpSolver().getSolution().col_value);
    updateUpperBound(lp->getObjective());
  }

  domain.backtrack();
#ifdef HIGHS_EXTRA_DEBUG
  for (int i = 0; i != mipsolver.numCol(); ++i) {
    assert(oldlower[i] == domain.colLower_[i]);
    assert(oldupper[i] == domain.colUpper_[i]);
  }
#endif
  lp->flushDomain(domain);

#ifdef HIGHS_EXTRA_DEBUG
  for (int i = 0; i != mipsolver.numCol(); ++i) {
    assert(oldlowerlp[i] == lp->getLpSolver().getLp().colLower_[i]);
    assert(oldupperlp[i] == lp->getLpSolver().getLp().colUpper_[i]);
  }
#endif

  lp->recoverBasis();
  int lim;
  lp->getLpSolver().getHighsOptionValue("simplex_iteration_limit", lim);
  lp->getLpSolver().setHighsOptionValue("simplex_iteration_limit",
                                        HIGHS_CONST_I_INF);
  oldnumiter = lp->getNumLpIterations();
  lp->resolveLp();
  lpiterations += lp->getNumLpIterations() - oldnumiter;
  lp->getLpSolver().setHighsOptionValue("simplex_iteration_limit", lim);
  // printf("restoring LP after rounding took %d iterations\n",
  // lp->getLpSolver().getSimplexIterationCount());
  if (success) printDisplayLine('R');
}
#endif

double HighsSearch::getCutoffBound() const {
  return std::min(mipsolver.mipdata_->upper_limit, upper_limit);
}

void HighsSearch::setRINSNeighbourhood(const std::vector<double>& basesol,
                                       const std::vector<double>& relaxsol) {
  for (int i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;
    if (localdom.colLower_[i] == localdom.colUpper_[i]) continue;

    double intval = std::floor(basesol[i] + 0.5);
    if (std::abs(relaxsol[i] - intval) < mipsolver.mipdata_->feastol) {
      if (localdom.colLower_[i] < intval)
        localdom.changeBound(HighsBoundType::Lower, i,
                             std::min(intval, localdom.colUpper_[i]), -2);
      if (localdom.colUpper_[i] > intval)
        localdom.changeBound(HighsBoundType::Upper, i,
                             std::max(intval, localdom.colLower_[i]), -2);
    }
  }
}

void HighsSearch::setRENSNeighbourhood(const std::vector<double>& lpsol) {
  for (int i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;
    if (localdom.colLower_[i] == localdom.colUpper_[i]) continue;

    double downval = std::floor(lpsol[i] + mipsolver.mipdata_->feastol);
    double upval = std::ceil(lpsol[i] - mipsolver.mipdata_->feastol);

    if (localdom.colLower_[i] < downval) {
      localdom.changeBound(HighsBoundType::Lower, i,
                           std::min(downval, localdom.colUpper_[i]), -2);
      if (localdom.infeasible()) return;
    }
    if (localdom.colUpper_[i] > upval) {
      localdom.changeBound(HighsBoundType::Upper, i,
                           std::max(upval, localdom.colLower_[i]), -2);
      if (localdom.infeasible()) return;
    }
  }
}

void HighsSearch::heuristicSearchNew() {
  const auto& lpsol = lp->getLpSolver().getSolution().col_value;
  HighsSearch heur(mipsolver, pseudocost);
  heur.localdom.setParentDomain(&localdom);
  heur.inheuristic = true;
  if (mipsolver.mipdata_->incumbent.empty()) {
    heur.setRENSNeighbourhood(lpsol);
    heur.localdom.propagate();
    if (heur.localdom.infeasible() || heur.localdom.getChangedCols().empty())
      return;

    int nfixed = 0;
    int ntotal = 0;
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      // skip fixed and continuous variables
      if (mipsolver.variableType(i) == HighsVarType::CONTINUOUS) continue;
      if (mipsolver.mipdata_->domain.colLower_[i] ==
          mipsolver.mipdata_->domain.colUpper_[i])
        continue;

      ++ntotal;
      // count locally fixed variable
      if (heur.localdom.colLower_[i] == heur.localdom.colUpper_[i]) ++nfixed;
    }
    double fixingrate = nfixed / (double)ntotal;
    // printf("RENS fixing rate: %.2f\n", fixingrate);
    if (fixingrate < 0.1) {
      // if the fixing rate is too small we just perform a dive with some
      // backtracks in this neighborhood and do not create a submip
      HighsLpRelaxation heurlp(*lp);
      // only use the global upper limit as LP limit so that dual proofs are
      // valid
      heurlp.setObjectiveLimit(mipsolver.mipdata_->upper_limit);
      heur.setLpRelaxation(&heurlp);

      heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                            heur.localdom.colLower_.data(),
                                            heur.localdom.colUpper_.data());
      heur.localdom.clearChangedCols();
      heur.nodestack.emplace_back();
      heur.pseudocost.setMinReliable(0);
      heur.solveDepthFirst(10);
      heurlpiterations += heur.lpiterations;
      lpiterations += heur.lpiterations;
      pseudocost = heur.pseudocost;
      return;
    }
    solveSubMip(std::move(heur.localdom.colLower_),
                std::move(heur.localdom.colUpper_),
                500,  // std::max(50, int(0.05 *
                      // (mipsolver.mipdata_->num_leaves))),
                200 + int(0.05 * (mipsolver.mipdata_->num_nodes + nnodes)));
    return;
  }

  HighsLpRelaxation heurlp(*lp);
  // only use the global upper limit as LP limit so that dual proofs are valid
  heurlp.setObjectiveLimit(mipsolver.mipdata_->upper_limit);
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                        heur.localdom.colLower_.data(),
                                        heur.localdom.colUpper_.data());
  heur.localdom.clearChangedCols();
  heur.nodestack.emplace_back();
  int nbacktracks = 0;

  // determine the initial number of unfixed variables fixing rate to decide if
  // the problem is restricted enough to be considered for solving a submip
  static constexpr double minfixingrate = 0.25;
  int ntotal = 0;
  int nfixed = 0;
  double fixingrate = 0.0;
  bool stop = false;

  while (true) {
    heur.evaluateNode();
    if (heur.currentNodePruned()) {
      ++nbacktracks;
      // printf("backtrack1\n");
      if (!heur.backtrack()) break;
      continue;
    }

    // if we estimate that there is no improving solution in this subtree, we
    // stop fixing variables but still backtrack to a node that has a good
    // estimate and is not pruned as the stop flag is checked after the
    // acktracking
    if ((heur.getCurrentEstimate() > getCutoffBound())) {
      heur.nodestack.back().opensubtrees = 0;
      ++nbacktracks;
      // printf("backtrack2\n");
      if (!heur.backtrack()) break;
      stop = true;
      continue;
    }

    nfixed = 0;
    ntotal = 0;
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      // skip fixed and continuous variables
      if (mipsolver.variableType(i) == HighsVarType::CONTINUOUS) continue;
      if (mipsolver.mipdata_->domain.colLower_[i] ==
          mipsolver.mipdata_->domain.colUpper_[i])
        continue;

      ++ntotal;
      // count locally fixed variable
      if (heur.localdom.colLower_[i] == heur.localdom.colUpper_[i]) ++nfixed;
    }
    fixingrate = nfixed / (double)ntotal;

    if (stop) break;
    if (nbacktracks >= 10) break;

    std::vector<std::pair<int, double>>::iterator fixcandend;

    // partition the fractional variables to consider which ones should we fix
    // in this dive first if there is an incumbent, we dive towards the RINS
    // neighborhood
    fixcandend = std::partition(
        heurlp.getFractionalIntegers().begin(),
        heurlp.getFractionalIntegers().end(),
        [&](const std::pair<int, double>& fracvar) {
          return std::abs(lpsol[fracvar.first] -
                          mipsolver.mipdata_->incumbent[fracvar.first]) <=
                 mipsolver.mipdata_->feastol;
        });

    bool fixtolpsol = true;

    auto getFixVal = [&](int col, double fracval) {
      double fixval;
      if (fixtolpsol) {
        // RINS neighborhood (with extension)
        fixval = std::floor(lpsol[col] + 0.5);
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
      fixval = std::min(heur.localdom.colUpper_[col], fixval);
      fixval = std::max(heur.localdom.colLower_[col], fixval);
      return fixval;
    };

    // no candidates left to fix for getting to the neighborhood, therefore we
    // switch to a different diving strategy until the minimal fixing rate is
    // reached
    if (heurlp.getFractionalIntegers().begin() == fixcandend) {
      if (fixingrate >= minfixingrate)
        break;  // if the RINS neigborhood achieved a high enough fixing rate by
                // itself we stop here
      fixcandend = heurlp.getFractionalIntegers().end();
      // now sort the variables by their distance towards the value they will be
      // fixed to
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

      if (heur.localdom.colLower_[fracint.first] < fixval) {
        heur.nodestack.back().branching_point = fracint.second;
        heur.nodestack.back().branchingdecision.column = fracint.first;
        heur.nodestack.back().branchingdecision.boundval = fixval;
        heur.nodestack.back().branchingdecision.boundtype =
            HighsBoundType::Lower;
        heur.nodestack.back().opensubtrees = 1;
        heur.localdom.changeBound(heur.nodestack.back().branchingdecision);
        heur.nodestack.emplace_back(heur.nodestack.back().lower_bound,
                                    heur.nodestack.back().estimate);
        if (heur.localdom.infeasible()) break;
      }

      if (heur.localdom.colUpper_[fracint.first] > fixval) {
        heur.nodestack.back().branching_point = fracint.second;
        heur.nodestack.back().branchingdecision.column = fracint.first;
        heur.nodestack.back().branchingdecision.boundval = fixval;
        heur.nodestack.back().branchingdecision.boundtype =
            HighsBoundType::Upper;
        heur.nodestack.back().opensubtrees = 1;
        heur.localdom.changeBound(heur.nodestack.back().branchingdecision);
        heur.nodestack.emplace_back(heur.nodestack.back().lower_bound,
                                    heur.nodestack.back().estimate);
        if (heur.localdom.infeasible()) break;
      }

      ++nfixed;
      fixingrate = nfixed / (double)ntotal;
      //  if (fixingrate >= maxfixingrate) break;
      change += std::abs(fixval - fracint.second);
      if (change >= 1.0) break;
    }

    heurlp.flushDomain(heur.localdom);

    // printf("%d/%d fixed, fixingrate is %g\n", nfixed, ntotal, fixingrate);
  }

  heurlpiterations += heur.lpiterations;
  lpiterations += heur.lpiterations;
  heur.lpiterations = 0;

  // if there is no node left it means we backtracked to the global domain and
  // the subproblem was solved with the dive
  if (!heur.hasNode()) return;

  // determine the fixing rate to decide if the problem is restricted enough to
  // be considered for solving a submip

  // printf("fixing rate is %g\n", fixingrate);
  if (fixingrate < 0.1) {
    // heur.childselrule = ChildSelectionRule::BestCost;
    heur.pseudocost.setMinReliable(0);
    heur.solveDepthFirst(10);
    heurlpiterations += heur.lpiterations;
    lpiterations += heur.lpiterations;
    pseudocost = heur.pseudocost;
    return;
  }

  solveSubMip(std::move(heur.localdom.colLower_),
              std::move(heur.localdom.colUpper_),
              500,  // std::max(50, int(0.05 *
                    // (mipsolver.mipdata_->num_leaves))),
              200 + int(0.05 * (mipsolver.mipdata_->num_nodes + nnodes)));
}

void HighsSearch::heuristicSearch() {
  const auto& lpsol = lp->getLpSolver().getSolution().col_value;

  HighsSearch heur(mipsolver, pseudocost);
  heur.localdom.setParentDomain(&localdom);
  heur.inheuristic = true;
  if (!mipsolver.mipdata_->incumbent.empty()) {
    heur.setRINSNeighbourhood(mipsolver.mipdata_->incumbent, lpsol);
  } else {
    // in case we have no incumbent we use the rens neighbourhood
    heur.setRENSNeighbourhood(lpsol);
  }

  if (heur.localdom.infeasible()) return;

  heur.localdom.propagate();
  if (heur.localdom.infeasible()) return;

  double objlim = mipsolver.mipdata_->upper_bound != HIGHS_CONST_INF
                      ? mipsolver.mipdata_->upper_bound -
                            std::abs(0.01 * (mipsolver.mipdata_->upper_bound -
                                             nodestack.back().lower_bound))
                      : HIGHS_CONST_INF;
  if (mipsolver.mipdata_->objintscale != 0.0 && objlim != HIGHS_CONST_INF) {
    objlim = std::floor(mipsolver.mipdata_->objintscale * objlim) /
                 mipsolver.mipdata_->objintscale +
             mipsolver.mipdata_->feastol;
  }
  heur.upper_limit = std::min(objlim, getCutoffBound());

  HighsLpRelaxation heurlp(*lp);
  // only use the global upper limit as LP limit so that dual proofs are valid
  heurlp.setObjectiveLimit(mipsolver.mipdata_->upper_limit);
  heur.setLpRelaxation(&heurlp);

  heurlp.getLpSolver().changeColsBounds(
      0, mipsolver.numCol() - 1, mipsolver.mipdata_->domain.colLower_.data(),
      mipsolver.mipdata_->domain.colUpper_.data());
  heurlp.flushDomain(heur.localdom);

  size_t lpitertarget =
      mipsolver.mipdata_->heuristic_effort * getTotalLpIterations() -
      getHeuristicLpIterations();

  size_t backtracklim =
      std::ceil(lpitertarget * mipsolver.mipdata_->num_leaves /
                (double)mipsolver.mipdata_->total_lp_iterations);
  if (backtracklim < 10) backtracklim = 10;

  heur.pseudocost.setMinReliable(
      0);  // std::min(1, pseudocost.getMinReliable()));
  heur.nodestack.emplace_back();
  double cutoffbnd = getCutoffBound();
  heur.solveDepthFirst(backtracklim);

  if (mipsolver.mipdata_->upper_limit < cutoffbnd)
    lp->setObjectiveLimit(mipsolver.mipdata_->upper_limit);

  pseudocost = heur.pseudocost;

  heurlpiterations += heur.lpiterations;
  lpiterations += heur.lpiterations;
}

void HighsSearch::addBoundExceedingConflict() {
  if (mipsolver.mipdata_->upper_limit != HIGHS_CONST_INF) {
    double rhs;
    if (lp->computeDualProof(mipsolver.mipdata_->domain,
                             mipsolver.mipdata_->upper_limit, inds, vals,
                             rhs)) {
      HighsSeparation::computeAndAddConflictCut(mipsolver, localdom, inds, vals,
                                                rhs);
      // int cutind = cutpool.addCut(inds.data(), vals.data(), inds.size(),
      // rhs); localdom.cutAdded(cutind);
    }
  }
}

void HighsSearch::addInfeasibleConflict() {
  double rhs;
  if (lp->computeDualInfProof(mipsolver.mipdata_->domain, inds, vals, rhs)) {
    // double minactlocal = 0.0;
    // double minactglobal = 0.0;
    // for (int i = 0; i < int(inds.size()); ++i) {
    //  if (vals[i] > 0.0) {
    //    minactlocal += localdom.colLower_[inds[i]] * vals[i];
    //    minactglobal += globaldom.colLower_[inds[i]] * vals[i];
    //  } else {
    //    minactlocal += localdom.colUpper_[inds[i]] * vals[i];
    //    minactglobal += globaldom.colUpper_[inds[i]] * vals[i];
    //  }
    //}
    // int oldnumcuts = cutpool.getNumCuts();
    HighsSeparation::computeAndAddConflictCut(mipsolver, localdom, inds, vals,
                                              rhs);

    // if (cutpool.getNumCuts() > oldnumcuts) {
    //  printf(
    //      "added cut from infeasibility proof with local min activity %g, "
    //      "global min activity %g, and rhs %g\n",
    //      minactlocal, minactglobal, rhs);
    //} else {
    //  printf(
    //      "no cut found for infeasibility proof with local min activity %g, "
    //      "global min "
    //      " activity %g, and rhs % g\n ",
    //      minactlocal, minactglobal, rhs);
    //}
    // int cutind = cutpool.addCut(inds.data(), vals.data(), inds.size(), rhs);
    // localdom.cutAdded(cutind);
  }
}

int HighsSearch::selectBranchingCandidate() {
  assert(!lp->getFractionalIntegers().empty());

  static constexpr int basisstart_threshold = 20;
  std::vector<double> upscore;
  std::vector<double> downscore;
  std::vector<uint8_t> upscorereliable;
  std::vector<uint8_t> downscorereliable;

  int numfrac = lp->getFractionalIntegers().size();
  const auto& fracints = lp->getFractionalIntegers();

  upscore.resize(numfrac, HIGHS_CONST_INF);
  downscore.resize(numfrac, HIGHS_CONST_INF);

  upscorereliable.resize(numfrac, 0);
  downscorereliable.resize(numfrac, 0);

  // initialize up and down scores of variables that have a
  // reliable pseudocost so that they do not get evaluated
  for (int k = 0; k != numfrac; ++k) {
    int col = fracints[k].first;

    if (pseudocost.isReliable(col) || branchingVarReliableAtNode(col)) {
      upscore[k] = pseudocost.getPseudocostUp(col, fracints[k].second);
      downscore[k] = pseudocost.getPseudocostDown(col, fracints[k].second);
      upscorereliable[k] = true;
      downscorereliable[k] = true;
    }
  }

  std::vector<int> evalqueue;
  evalqueue.resize(numfrac);
  std::iota(evalqueue.begin(), evalqueue.end(), 0);
  auto selectBestScore = [&]() {
    int best = -1;
    double bestscore = -1.0;
    for (int k : evalqueue) {
      double score;
      if ((upscore[k] == 0.0 && upscorereliable[k]) ||
          (downscore[k] == 0.0 && downscorereliable[k]))
        score = pseudocost.getScore(fracints[k].first, 0.0, 0.0);
      else
        score = upscore[k] == HIGHS_CONST_INF || downscore[k] == HIGHS_CONST_INF
                    ? HIGHS_CONST_INF
                    : pseudocost.getScore(fracints[k].first, upscore[k],
                                          downscore[k]);

      if (score > bestscore) {
        bestscore = score;
        best = k;
      }
    }

    return best;
  };

  lp->storeBasis();
  auto basis = lp->getStoredBasis();

  while (true) {
    int candidate = selectBestScore();

    if (upscorereliable[candidate] && downscorereliable[candidate]) {
      lp->setStoredBasis(std::move(basis));
      lp->recoverBasis();
      lp->run();
      return candidate;
    }

    int col = fracints[candidate].first;
    double fracval = fracints[candidate].second;
    double upval = std::ceil(fracval);
    double downval = std::floor(fracval);

    if (!downscorereliable[candidate]) {
      // evaluate down branch
      localdom.changeBound(HighsBoundType::Upper, col, downval);
      localdom.propagate();
      if (localdom.infeasible()) {
        localdom.backtrack();
        localdom.clearChangedCols();
        localdom.changeBound(HighsBoundType::Lower, col, upval, -2);
        lp->setStoredBasis(std::move(basis));
        return -1;
      }

      lp->flushDomain(localdom);

      size_t numiters = lp->getNumLpIterations();
      HighsLpRelaxation::Status status = lp->run(false);
      numiters = lp->getNumLpIterations() - numiters;
      lpiterations += numiters;
      sblpiterations += numiters;

      if (lp->scaledOptimal(status)) {
        double delta = downval - fracval;
        bool integerfeasible;
        const std::vector<double>& sol =
            lp->getLpSolver().getSolution().col_value;
        double solobj = checkSol(sol, integerfeasible);

        double objdelta = std::max(solobj - lp->getObjective(), 0.0);
        if (objdelta < mipsolver.mipdata_->epsilon) objdelta = 0.0;

        downscore[candidate] = objdelta;
        downscorereliable[candidate] = 1;
        markBranchingVarDownReliableAtNode(col);
        pseudocost.addObservation(col, delta, objdelta);

        for (int k = 0; k != numfrac; ++k) {
          double otherfracval = fracints[k].second;
          double otherdownval = std::floor(fracints[k].second);
          double otherupval = std::ceil(fracints[k].second);
          if (sol[fracints[k].first] <=
              otherdownval + mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && downscore[k] != 0.0) {
              downscorereliable[k] = 1;
              markBranchingVarDownReliableAtNode(fracints[k].first);
              pseudocost.addObservation(fracints[k].first,
                                        otherdownval - otherfracval, objdelta);
            }
            downscore[k] = std::min(downscore[k], objdelta);
          } else if (sol[fracints[k].first] >=
                     otherupval - mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && upscore[k] != 0.0) {
              upscorereliable[k] = 1;
              markBranchingVarUpReliableAtNode(fracints[k].first);
              pseudocost.addObservation(fracints[k].first,
                                        otherupval - otherfracval, objdelta);
            }
            upscore[k] = std::min(upscore[k], objdelta);
          }
        }

        if (lp->unscaledPrimalFeasible(status) && integerfeasible) {
          double cutoffbnd = getCutoffBound();
          mipsolver.mipdata_->addIncumbent(
              lp->getLpSolver().getSolution().col_value, solobj,
              inheuristic ? 'H' : 'B');

          if (mipsolver.mipdata_->upper_limit < cutoffbnd)
            lp->setObjectiveLimit(mipsolver.mipdata_->upper_limit);
        }

        if (lp->unscaledDualFeasible(status)) {
          if (solobj > getCutoffBound()) {
            addBoundExceedingConflict();
            localdom.backtrack();
            lp->flushDomain(localdom);
            localdom.changeBound(HighsBoundType::Lower, col, upval, -2);
            lp->setStoredBasis(std::move(basis));
            if (numiters > basisstart_threshold) lp->recoverBasis();
            return -1;
          }
        } else if (solobj > getCutoffBound()) {
          addBoundExceedingConflict();
          localdom.propagate();
          bool infeas = localdom.infeasible();
          if (infeas) {
            localdom.backtrack();
            lp->flushDomain(localdom);
            localdom.changeBound(HighsBoundType::Lower, col, upval, -2);
            lp->setStoredBasis(std::move(basis));
            if (numiters > basisstart_threshold) lp->recoverBasis();
            return -1;
          }
        }
      } else if (status == HighsLpRelaxation::Status::Infeasible) {
        addInfeasibleConflict();
        localdom.backtrack();
        lp->flushDomain(localdom);
        localdom.changeBound(HighsBoundType::Lower, col, upval, -2);
        lp->setStoredBasis(std::move(basis));
        if (numiters > basisstart_threshold) lp->recoverBasis();
        return -1;
      } else {
        // printf("todo2\n");
        // in case of an LP error we set the score of this variable to zero to
        // avoid choosing it as branching candidate if possible
        downscore[candidate] = 0.0;
        upscore[candidate] = 0.0;
        downscorereliable[candidate] = 1;
        upscorereliable[candidate] = 1;
        markBranchingVarUpReliableAtNode(col);
        markBranchingVarDownReliableAtNode(col);
      }

      localdom.backtrack();
      lp->flushDomain(localdom);
      if (numiters > basisstart_threshold) lp->recoverBasis();
    } else {
      // evaluate up branch
      localdom.changeBound(HighsBoundType::Lower, col, upval);
      localdom.propagate();
      if (localdom.infeasible()) {
        localdom.backtrack();
        localdom.clearChangedCols();
        localdom.changeBound(HighsBoundType::Upper, col, downval, -2);
        lp->setStoredBasis(std::move(basis));
        return -1;
      }

      lp->flushDomain(localdom);

      size_t numiters = lp->getNumLpIterations();
      HighsLpRelaxation::Status status = lp->run(false);
      numiters = lp->getNumLpIterations() - numiters;
      lpiterations += numiters;
      sblpiterations += numiters;

      if (lp->scaledOptimal(status)) {
        double delta = upval - fracval;
        bool integerfeasible;

        const std::vector<double>& sol =
            lp->getLpSolver().getSolution().col_value;
        double solobj = checkSol(sol, integerfeasible);

        double objdelta = std::max(solobj - lp->getObjective(), 0.0);
        if (objdelta < mipsolver.mipdata_->epsilon) objdelta = 0.0;

        upscore[candidate] = objdelta;
        upscorereliable[candidate] = 1;
        markBranchingVarUpReliableAtNode(col);
        pseudocost.addObservation(col, delta, objdelta);

        for (int k = 0; k != numfrac; ++k) {
          double otherfracval = fracints[k].second;
          double otherdownval = std::floor(fracints[k].second);
          double otherupval = std::ceil(fracints[k].second);
          if (sol[fracints[k].first] <=
              otherdownval + mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && downscore[k] != 0.0) {
              downscorereliable[k] = 1;
              markBranchingVarDownReliableAtNode(fracints[k].first);
              pseudocost.addObservation(fracints[k].first,
                                        otherdownval - otherfracval, objdelta);
            }
            downscore[k] = std::min(downscore[k], objdelta);

          } else if (sol[fracints[k].first] >=
                     otherupval - mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && upscore[k] != 0.0) {
              upscorereliable[k] = 1;
              markBranchingVarUpReliableAtNode(fracints[k].first);
              pseudocost.addObservation(fracints[k].first,
                                        otherupval - otherfracval, objdelta);
            }
            upscore[k] = std::min(upscore[k], objdelta);
          }
        }

        if (lp->unscaledPrimalFeasible(status) && integerfeasible) {
          double cutoffbnd = getCutoffBound();
          mipsolver.mipdata_->addIncumbent(
              lp->getLpSolver().getSolution().col_value, solobj,
              inheuristic ? 'H' : 'B');

          if (mipsolver.mipdata_->upper_limit < cutoffbnd)
            lp->setObjectiveLimit(mipsolver.mipdata_->upper_limit);
        }

        if (lp->unscaledDualFeasible(status)) {
          if (solobj > getCutoffBound()) {
            addBoundExceedingConflict();
            localdom.backtrack();
            lp->flushDomain(localdom);
            localdom.changeBound(HighsBoundType::Upper, col, downval, -2);
            lp->setStoredBasis(std::move(basis));
            if (numiters > basisstart_threshold) lp->recoverBasis();
            return -1;
          }
        } else if (solobj > getCutoffBound()) {
          addBoundExceedingConflict();
          localdom.propagate();
          bool infeas = localdom.infeasible();
          if (infeas) {
            localdom.backtrack();
            lp->flushDomain(localdom);
            localdom.changeBound(HighsBoundType::Upper, col, downval, -2);
            lp->setStoredBasis(std::move(basis));
            if (numiters > basisstart_threshold) lp->recoverBasis();
            return -1;
          }
        }
      } else if (status == HighsLpRelaxation::Status::Infeasible) {
        addInfeasibleConflict();
        localdom.backtrack();
        lp->flushDomain(localdom);
        localdom.changeBound(HighsBoundType::Upper, col, downval, -2);
        lp->setStoredBasis(std::move(basis));
        if (numiters > basisstart_threshold) lp->recoverBasis();
        return -1;
      } else {
        // printf("todo2\n");
        // in case of an LP error we set the score of this variable to zero to
        // avoid choosing it as branching candidate if possible
        downscore[candidate] = 0.0;
        upscore[candidate] = 0.0;
        downscorereliable[candidate] = 1;
        upscorereliable[candidate] = 1;
        markBranchingVarUpReliableAtNode(col);
        markBranchingVarDownReliableAtNode(col);
      }

      localdom.backtrack();
      lp->flushDomain(localdom);
      if (numiters > basisstart_threshold) lp->recoverBasis();
    }
  }
}

const HighsSearch::NodeData* HighsSearch::getParentNodeData() const {
  if (nodestack.size() <= 1) return nullptr;

  return &nodestack[nodestack.size() - 2];
}

void HighsSearch::currentNodeToQueue(HighsNodeQueue& nodequeue) {
  localdom.propagate();
  if (!localdom.infeasible()) {
    nodequeue.emplaceNode(localdom.getReducedDomainChangeStack(),
                          nodestack.back().lower_bound,
                          nodestack.back().estimate, getCurrentDepth());
  } else
    treeweight += std::pow(0.5, getCurrentDepth() - 1);
  nodestack.back().opensubtrees = 0;

  backtrack();
  lp->flushDomain(localdom);
}

void HighsSearch::openNodesToQueue(HighsNodeQueue& nodequeue) {
  if (nodestack.empty()) return;
  if (nodestack.back().opensubtrees == 0) backtrack();

  while (!nodestack.empty()) {
    localdom.propagate();
    if (!localdom.infeasible()) {
      nodequeue.emplaceNode(localdom.getReducedDomainChangeStack(),
                            nodestack.back().lower_bound,
                            nodestack.back().estimate, getCurrentDepth());
    } else
      treeweight += std::pow(0.5, getCurrentDepth() - 1);
    nodestack.back().opensubtrees = 0;

    backtrack();
  }

  lp->flushDomain(localdom);
}

void HighsSearch::solveSubMip(std::vector<double> colLower,
                              std::vector<double> colUpper, int maxleaves,
                              int maxnodes) {
  HighsOptions submipoptions = *mipsolver.options_mip_;
  HighsLp submip = *mipsolver.model_;

  // set bounds and restore integrality of the lp relaxation copy
  submip.colLower_ = std::move(colLower);
  submip.colUpper_ = std::move(colUpper);
  submip.integrality_ = mipsolver.model_->integrality_;
  submip.offset_ = 0;

  // set limits
  submipoptions.mip_max_leaves = maxleaves;
  submipoptions.logfile = nullptr;
  submipoptions.output = nullptr;
  submipoptions.mip_max_nodes = maxnodes;
  submipoptions.time_limit -=
      mipsolver.timer_.read(mipsolver.timer_.solve_clock);
  submipoptions.dual_objective_value_upper_bound =
      mipsolver.mipdata_->upper_limit;
  // setup solver and run it
  HighsMipSolver submipsolver(submipoptions, submip, true);
  submipsolver.rootbasis = &mipsolver.mipdata_->firstrootbasis;
  submipsolver.run();

  if (submipsolver.modelstatus_ != HighsModelStatus::PRIMAL_INFEASIBLE &&
      !submipsolver.presolve_.data_.recovered_solution_.col_value.empty()) {
    bool integerfeasible;
    double solobj =
        checkSol(submipsolver.presolve_.data_.recovered_solution_.col_value,
                 integerfeasible);
    assert(std::isfinite(solobj));
    if (integerfeasible)
      mipsolver.mipdata_->addIncumbent(
          submipsolver.presolve_.data_.recovered_solution_.col_value, solobj,
          'L');
  }

  if (submipsolver.mipdata_) {
    double adjustmentfactor =
        submipsolver.numNonzero() / (double)mipsolver.numNonzero();
    size_t adjusted_lp_iterations =
        (size_t)(adjustmentfactor * adjustmentfactor *
                 submipsolver.mipdata_->total_lp_iterations);
    heurlpiterations += adjusted_lp_iterations;
    lpiterations += adjusted_lp_iterations;
  }
}

void HighsSearch::flushStatistics() {
  mipsolver.mipdata_->num_nodes += nnodes;
  nnodes = 0;

  mipsolver.mipdata_->pruned_treeweight += treeweight;
  treeweight = 0;

  mipsolver.mipdata_->total_lp_iterations += lpiterations;
  lpiterations = 0;

  mipsolver.mipdata_->heuristic_lp_iterations += heurlpiterations;
  heurlpiterations = 0;

  mipsolver.mipdata_->sb_lp_iterations += sblpiterations;
  sblpiterations = 0;
}

size_t HighsSearch::getHeuristicLpIterations() const {
  return heurlpiterations + mipsolver.mipdata_->heuristic_lp_iterations;
}

size_t HighsSearch::getTotalLpIterations() const {
  return lpiterations + mipsolver.mipdata_->total_lp_iterations;
}

size_t HighsSearch::getStrongBranchingLpIterations() const {
  return sblpiterations + mipsolver.mipdata_->sb_lp_iterations;
}

void HighsSearch::resetLocalDomain() {
  this->lp->getLpSolver().changeColsBounds(
      0, mipsolver.numCol() - 1, mipsolver.mipdata_->domain.colLower_.data(),
      mipsolver.mipdata_->domain.colUpper_.data());
  localdom = mipsolver.mipdata_->domain.createChildDomain();

#ifndef NDEBUG
  for (int i = 0; i != mipsolver.numCol(); ++i) {
    assert(lp->getLpSolver().getLp().colLower_[i] == localdom.colLower_[i] ||
           mipsolver.variableType(i) == HighsVarType::CONTINUOUS);
    assert(lp->getLpSolver().getLp().colUpper_[i] == localdom.colUpper_[i] ||
           mipsolver.variableType(i) == HighsVarType::CONTINUOUS);
  }
#endif
}

void HighsSearch::installNode(HighsNodeQueue::OpenNode&& node) {
  localdom.setDomainChangeStack(std::move(node.domchgstack));
  nodestack.emplace_back(node.lower_bound, node.estimate);
  depthoffset = node.depth - 1;
}

void HighsSearch::evaluateNode() {
#ifdef HIGHS_DEBUGSOL
  bool debugsolactive = true;
  HighsCDouble debugsolobj = 0;
  for (int i = 0; i != mipsolver.numCol(); ++i) {
    if (highsDebugSolution[i] + mipsolver.mipdata_->epsilon <
            localdom.colLower_[i] ||
        highsDebugSolution[i] - mipsolver.mipdata_->epsilon >
            localdom.colUpper_[i]) {
      debugsolactive = false;
    }

    debugsolobj += highsDebugSolution[i] * mipsolver.colCost(i);
  }
#endif
  localdom.propagate();
#ifdef HIGHS_DEBUGSOL
  if (debugsolactive && mipsolver.mipdata_->upper_bound >
                            debugsolobj + mipsolver.mipdata_->epsilon) {
    bool debugsolstillactive = true;

    for (int i = 0; i != mipsolver.numCol(); ++i) {
      if (highsDebugSolution[i] + mipsolver.mipdata_->epsilon <
              localdom.colLower_[i] ||
          highsDebugSolution[i] - mipsolver.mipdata_->epsilon >
              localdom.colUpper_[i]) {
        debugsolstillactive = false;
        break;
      }
    }

    assert(debugsolstillactive);
  }
#endif
  assert(!nodestack.empty());
  NodeData& currnode = nodestack.back();

  bool prune = false;

  if (localdom.infeasible()) {
#ifdef HIGHS_DEBUGSOL
    assert(!debugsolactive || mipsolver.mipdata_->upper_bound <=
                                  debugsolobj + mipsolver.mipdata_->feastol);
#endif
    localdom.clearChangedCols();
    prune = true;
  } else {
    lp->flushDomain(localdom);

#ifndef NDEBUG
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      assert(lp->getLpSolver().getLp().colLower_[i] == localdom.colLower_[i] ||
             mipsolver.variableType(i) == HighsVarType::CONTINUOUS);
      assert(lp->getLpSolver().getLp().colUpper_[i] == localdom.colUpper_[i] ||
             mipsolver.variableType(i) == HighsVarType::CONTINUOUS);
    }
#endif
    size_t oldnumiters = lp->getNumLpIterations();
    HighsLpRelaxation::Status status = lp->resolveLp();
    lpiterations += lp->getNumLpIterations() - oldnumiters;
    if (lp->scaledOptimal(status)) {
      lp->storeBasis();

      currnode.estimate = lp->computeBestEstimate(pseudocost);

      if (lp->unscaledPrimalFeasible(status)) {
        if (lp->getFractionalIntegers().empty()) {
          double cutoffbnd = getCutoffBound();
          mipsolver.mipdata_->addIncumbent(
              lp->getLpSolver().getSolution().col_value, lp->getObjective(),
              inheuristic ? 'H' : 'T');
          if (mipsolver.mipdata_->upper_limit < cutoffbnd)
            lp->setObjectiveLimit(mipsolver.mipdata_->upper_limit);
        }
      }

      if (lp->unscaledDualFeasible(status)) {
        currnode.lpsolved = true;
        currnode.lower_bound =
            std::max(lp->getObjective(), currnode.lower_bound);

#ifdef HIGHS_DEBUGSOL
        assert(!debugsolactive ||
               currnode.lower_bound <=
                   debugsolobj + mipsolver.mipdata_->epsilon);
#endif

        const NodeData* parent = getParentNodeData();

        if (parent != nullptr && parent->lpsolved &&
            parent->branching_point != parent->branchingdecision.boundval) {
          int col = parent->branchingdecision.column;
          double delta =
              parent->branchingdecision.boundval - parent->branching_point;
          double objdelta =
              std::max(0.0, currnode.lower_bound - parent->lower_bound);

          pseudocost.addObservation(col, delta, objdelta);
        }

        if (currnode.lower_bound > getCutoffBound()) {
#ifdef HIGHS_DEBUGSOL
          if (debugsolactive && mipsolver.mipdata_->upper_bound >
                                    debugsolobj + mipsolver.mipdata_->feastol) {
            lp->getLpSolver().writeModel("wronglp->mps");
            lp->getLpSolver().writeBasis("wronglp->bas");
            assert(false);
          }
#endif
          addBoundExceedingConflict();
          prune = true;
        }
      } else if (lp->getObjective() > getCutoffBound()) {
        // the LP is not solved to dual feasibilty due to scaling/numerics
        // therefore we compute a conflict constraint as if the LP was bound
        // exceeding and propagate the local domain again. The lp relaxation
        // class will take care to consider the dual multipliers with an
        // increased zero tolerance due to the dual infeasibility when computing
        // the proof constraint.
        addBoundExceedingConflict();
        localdom.propagate();
        prune = localdom.infeasible();
      }

    } else if (status == HighsLpRelaxation::Status::Infeasible) {
      addInfeasibleConflict();
#ifdef HIGHS_DEBUGSOL
      if (debugsolactive && mipsolver.mipdata_->upper_bound >
                                debugsolobj + mipsolver.mipdata_->feastol) {
        lp->getLpSolver().writeModel("wronglp->mps");
        lp->getLpSolver().writeBasis("wronglp->bas");
        assert(false);
      }
#endif
      prune = true;
    }
  }

#if 0 
      if( status == HighsModelStatus::PRIMAL_INFEASIBLE )
      {
        double rhs;
        if( lp->computeDualProof(globaldom, currnode.lower_bound, inds, vals, rhs) )
        {
          double glbminact = 0.0;
          double localminact = 0.0;

          for( size_t i = 0; i != inds.size(); ++i )
          {
            if( vals[i] < 0 )
            {
              glbminact += globaldom.colUpper_[inds[i]] * vals[i];
              localminact += localdom.colUpper_[inds[i]] * vals[i];
            }
            else
            {
              glbminact += globaldom.colLower_[inds[i]] * vals[i];
              localminact += localdom.colLower_[inds[i]] * vals[i];
            }
            
          }

          printf("glbminact: %g   localminact: %g   rhs: %g  len: %d\n", glbminact, localminact, rhs, (int)inds.size());
        }

      }
#endif

  if (prune) {
#ifdef HIGHS_DEBUGSOL
    if (debugsolactive && mipsolver.mipdata_->upper_bound >
                              debugsolobj + mipsolver.mipdata_->feastol) {
      bool debugsolstillactive = true;

      for (int i = 0; i != mipsolver.numCol(); ++i) {
        if (highsDebugSolution[i] + mipsolver.mipdata_->feastol <
                localdom.colLower_[i] ||
            highsDebugSolution[i] - mipsolver.mipdata_->feastol >
                localdom.colUpper_[i]) {
          debugsolstillactive = false;
          break;
        }
      }
      assert(debugsolstillactive);
    }
#endif
    treeweight += std::pow(0.5, getCurrentDepth() - 1);
    currnode.opensubtrees = 0;
  }
}

bool HighsSearch::branch() {
  assert(localdom.getChangedCols().empty());

  NodeData& currnode = nodestack.back();
  assert(currnode.opensubtrees == 2);
  currnode.branchingdecision.column = -1;
  inbranching = true;

  int minrel = pseudocost.getMinReliable();
  pseudocost.setSeed(random.integer());

  while (currnode.opensubtrees == 2 && lp->scaledOptimal(lp->getStatus()) &&
         !lp->getFractionalIntegers().empty()) {
    // evalUnreliableBranchCands();
    if (minrel > 0) {
      size_t sbiters = getStrongBranchingLpIterations();
      size_t sbmaxiters =
          100000 + (getTotalLpIterations() - getHeuristicLpIterations() -
                    getStrongBranchingLpIterations()) /
                       2;
      if (sbiters > sbmaxiters) {
        pseudocost.setMinReliable(0);
      } else if (sbiters > sbmaxiters / 2) {
        double reductionratio =
            (sbiters - sbmaxiters / 2) / (double)(sbmaxiters - sbmaxiters / 2);

        int minrelreduced = int(8.0 - reductionratio * 7.0);
        pseudocost.setMinReliable(std::min(minrel, minrelreduced));
      }
    }

    int branchcand = selectBranchingCandidate();

    if (branchcand != -1) {
      auto branching = lp->getFractionalIntegers()[branchcand];
      currnode.branchingdecision.column = branching.first;
      currnode.branching_point = branching.second;

      int col = branching.first;
      switch (childselrule) {
        case ChildSelectionRule::Up:
          currnode.branchingdecision.boundtype = HighsBoundType::Lower;
          currnode.branchingdecision.boundval =
              std::ceil(currnode.branching_point);
          break;
        case ChildSelectionRule::Down:
          currnode.branchingdecision.boundtype = HighsBoundType::Upper;
          currnode.branchingdecision.boundval =
              std::floor(currnode.branching_point);
          break;
        case ChildSelectionRule::RootSol:
          if (currnode.branching_point >= mipsolver.mipdata_->rootlpsol[col]) {
            currnode.branchingdecision.boundtype = HighsBoundType::Lower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::Upper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
          }
          break;
        case ChildSelectionRule::Obj:
          if (mipsolver.colCost(col) >= 0) {
            currnode.branchingdecision.boundtype = HighsBoundType::Lower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::Upper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
          }
          break;
        case ChildSelectionRule::Random:
          if (random.integer() % 2 == 0) {
            currnode.branchingdecision.boundtype = HighsBoundType::Lower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::Upper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
          }
          break;
        case ChildSelectionRule::BestCost: {
          if (pseudocost.getPseudocostUp(col, currnode.branching_point) >
              pseudocost.getPseudocostDown(col, currnode.branching_point)) {
            currnode.branchingdecision.boundtype = HighsBoundType::Upper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::Lower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
          }
          break;
        }
        case ChildSelectionRule::WorstCost:
          if (pseudocost.getPseudocostUp(col, currnode.branching_point) >=
              pseudocost.getPseudocostDown(col, currnode.branching_point)) {
            currnode.branchingdecision.boundtype = HighsBoundType::Lower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
          } else {
            currnode.branchingdecision.boundtype = HighsBoundType::Upper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
          }
          break;
      }
      break;
    }

    assert(!localdom.getChangedCols().empty());
    evaluateNode();
  }
  inbranching = false;
  pseudocost.setMinReliable(minrel);

  assert(currnode.opensubtrees == 2 || currnode.opensubtrees == 0);

  if (currnode.opensubtrees != 2) return false;

  if (currnode.branchingdecision.column == -1) {
    double bestscore = -1.0;
    // solution branching failed, so choose any integer variable to branch
    // on in case we have a different solution status could happen due to a
    // fail in the LP solution process
    pseudocost.setSeed(random.integer());

    for (int i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.variableType(i) == HighsVarType::CONTINUOUS) continue;

      if (localdom.colUpper_[i] - localdom.colLower_[i] < 0.5) continue;

      double fracval;
      if (localdom.colLower_[i] != -HIGHS_CONST_INF)
        fracval = localdom.colLower_[i] + 0.5;
      else if (localdom.colUpper_[i] != HIGHS_CONST_INF)
        fracval = localdom.colUpper_[i] - 0.5;
      else
        fracval = 0.5;

      double score = pseudocost.getScore(i, fracval);
      assert(score >= 0.0);
      if (score > bestscore) {
        bestscore = score;
        double upval = std::ceil(fracval);
        currnode.branching_point = upval;
        currnode.branchingdecision.boundtype = HighsBoundType::Lower;
        currnode.branchingdecision.column = i;
        currnode.branchingdecision.boundval = upval;
      }
    }
  }

  if (currnode.branchingdecision.column == -1) {
    lp->getLpSolver().clearSolver();
    mipsolver.mipdata_->cutpool.removeAllRows(*lp);
    lp->setIterationLimit();
    lp->getLpSolver().setHighsOptionValue("presolve", "on");
    evaluateNode();
    lp->getLpSolver().setHighsOptionValue("presolve", "off");

    if (currnode.opensubtrees != 0) {
      printf(
          "WARNING: all integers colls are fixed, LP may be unstable, possibly "
          "pruning optimal solution, lp status: scaled=%d unscaled=%d\n",
          (int)lp->getLpSolver().getModelStatus(true),
          (int)lp->getLpSolver().getModelStatus(false));
      currnode.opensubtrees = 0;
    }
    return false;
  }

  // finally open a new node with the branching decision added
  // and remember that we have one open subtree left
  localdom.changeBound(currnode.branchingdecision);
  currnode.opensubtrees = 1;
  nodestack.emplace_back(currnode.lower_bound, currnode.estimate);

  return true;
}

bool HighsSearch::backtrack() {
  assert(!nodestack.empty());
  assert(nodestack.back().opensubtrees == 0);

  while (nodestack.back().opensubtrees == 0) {
    nodestack.pop_back();

#ifndef NDEBUG
    HighsDomainChange branchchg =
#endif
        localdom.backtrack();
    if (nodestack.empty()) {
      lp->flushDomain(localdom);
      return false;
    }
    assert(branchchg.boundval == nodestack.back().branchingdecision.boundval);
    assert(branchchg.boundtype == nodestack.back().branchingdecision.boundtype);
    assert(branchchg.column == nodestack.back().branchingdecision.column);
  }

  NodeData& currnode = nodestack.back();

  assert(currnode.opensubtrees == 1);
  currnode.opensubtrees = 0;
  bool fallbackbranch =
      currnode.branchingdecision.boundval == currnode.branching_point;

  if (currnode.branchingdecision.boundtype == HighsBoundType::Lower) {
    currnode.branchingdecision.boundtype = HighsBoundType::Upper;
    currnode.branchingdecision.boundval =
        std::floor(currnode.branchingdecision.boundval - 0.5);
  } else {
    currnode.branchingdecision.boundtype = HighsBoundType::Lower;
    currnode.branchingdecision.boundval =
        std::ceil(currnode.branchingdecision.boundval + 0.5);
  }

  if (fallbackbranch)
    currnode.branching_point = currnode.branchingdecision.boundval;

  localdom.changeBound(currnode.branchingdecision);
  nodestack.emplace_back(currnode.lower_bound, currnode.estimate);
  lp->flushDomain(localdom);

  return true;
}

void HighsSearch::dive() {
  reliableatnode.clear();

  do {
    ++nnodes;
    evaluateNode();

    if (nodestack.back().opensubtrees == 0) return;

    if (!branch()) return;

  } while (true);
}

void HighsSearch::solveDepthFirst(size_t maxbacktracks) {
  do {
    if (maxbacktracks == 0) break;

    dive();

    --maxbacktracks;

  } while (backtrack());
}