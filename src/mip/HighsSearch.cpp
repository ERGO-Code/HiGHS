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

    if (localdom.colLower_[i] < downval)
      localdom.changeBound(HighsBoundType::Lower, i,
                           std::min(downval, localdom.colUpper_[i]), -2);
    if (localdom.colUpper_[i] > upval)
      localdom.changeBound(HighsBoundType::Upper, i,
                           std::max(upval, localdom.colLower_[i]), -2);
  }
}

void HighsSearch::heuristicSearch() {
  HighsSearch heur(mipsolver, pseudocost);

  heur.childselrule = ChildSelectionRule::BestCost;
  heur.localdom.setParentDomain(&localdom);
  heur.inheuristic = true;
  const auto& lpsol = lp->getLpSolver().getSolution().col_value;
  if (!mipsolver.mipdata_->incumbent.empty()) {
    heur.setRINSNeighbourhood(mipsolver.mipdata_->incumbent, lpsol);
  } else {
    // in case we have no incumbent we use the rens neighbourhood
    heur.setRENSNeighbourhood(lpsol);
  }

  heur.localdom.propagate();
  if (heur.localdom.infeasible()) return;

  double objlim = mipsolver.mipdata_->upper_bound != HIGHS_CONST_INF
                      ? mipsolver.mipdata_->upper_bound -
                            std::abs(0.01 * (mipsolver.mipdata_->upper_bound -
                                             nodestack.back().lower_bound))
                      : HIGHS_CONST_INF;

  heur.upper_limit = objlim;

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

  heur.pseudocost.setMinReliable(1);
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

    if (pseudocost.isReliable(col)) {
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

      size_t oldnumiters = lp->getNumLpIterations();
      HighsLpRelaxation::Status status = lp->run();
      lpiterations += lp->getNumLpIterations() - oldnumiters;

      if (lp->scaledOptimal(status)) {
        double delta = downval - fracval;
        bool integerfeasible;
        const std::vector<double>& sol =
            lp->getLpSolver().getSolution().col_value;
        double solobj = checkSol(sol, integerfeasible);

        double objdelta = std::max(solobj - lp->getObjective(), 0.0);
        if (objdelta < 1e-9) objdelta = 0.0;

        downscore[candidate] = objdelta;
        downscorereliable[candidate] = 1;
        pseudocost.addObservation(col, delta, objdelta);

        for (int k = 0; k != numfrac; ++k) {
          double otherfracval = fracints[k].second;
          double otherdownval = std::floor(fracints[k].second);
          double otherupval = std::ceil(fracints[k].second);
          if (sol[fracints[k].first] <=
              otherdownval + mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && downscore[k] != 0.0) {
              downscorereliable[k] = 1;
              pseudocost.addObservation(fracints[k].first,
                                        otherdownval - otherfracval, objdelta);
            }
            downscore[k] = std::min(downscore[k], objdelta);
          } else if (sol[fracints[k].first] >=
                     otherupval - mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && upscore[k] != 0.0) {
              upscorereliable[k] = 1;
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
          if (solobj >= getCutoffBound()) {
            addBoundExceedingConflict();
            localdom.backtrack();
            lp->flushDomain(localdom);
            localdom.changeBound(HighsBoundType::Lower, col, upval, -2);
            lp->setStoredBasis(std::move(basis));
            lp->recoverBasis();
            return -1;
          }
        }
      } else if (status == HighsLpRelaxation::Status::Infeasible) {
        addInfeasibleConflict();
        localdom.backtrack();
        lp->flushDomain(localdom);
        localdom.changeBound(HighsBoundType::Lower, col, upval, -2);
        lp->setStoredBasis(std::move(basis));
        lp->recoverBasis();
        return -1;
      } else {
        // printf("todo2\n");
        // in case of an LP error we set the score of this variable to zero to
        // avoid choosing it as branching candidate if possible
        downscore[candidate] = 0.0;
        upscore[candidate] = 0.0;
        downscorereliable[candidate] = 1;
        upscorereliable[candidate] = 1;
      }

      localdom.backtrack();
      lp->flushDomain(localdom);
      lp->recoverBasis();
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

      size_t oldnumiters = lp->getNumLpIterations();
      HighsLpRelaxation::Status status = lp->run();
      lpiterations += lp->getNumLpIterations() - oldnumiters;

      if (lp->scaledOptimal(status)) {
        double delta = upval - fracval;
        bool integerfeasible;

        const std::vector<double>& sol =
            lp->getLpSolver().getSolution().col_value;
        double solobj = checkSol(sol, integerfeasible);

        double objdelta = std::max(solobj - lp->getObjective(), 0.0);
        if (objdelta < 1e-9) objdelta = 0.0;

        upscore[candidate] = objdelta;
        upscorereliable[candidate] = 1;
        pseudocost.addObservation(col, delta, objdelta);

        for (int k = 0; k != numfrac; ++k) {
          double otherfracval = fracints[k].second;
          double otherdownval = std::floor(fracints[k].second);
          double otherupval = std::ceil(fracints[k].second);
          if (sol[fracints[k].first] <=
              otherdownval + mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && downscore[k] != 0.0) {
              downscorereliable[k] = 1;
              pseudocost.addObservation(fracints[k].first,
                                        otherdownval - otherfracval, objdelta);
            }
            downscore[k] = std::min(downscore[k], objdelta);

          } else if (sol[fracints[k].first] >=
                     otherupval - mipsolver.mipdata_->feastol) {
            if (objdelta == 0.0 && upscore[k] != 0.0) {
              upscorereliable[k] = 1;
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
          if (solobj >= getCutoffBound()) {
            addBoundExceedingConflict();
            localdom.backtrack();
            lp->flushDomain(localdom);
            localdom.changeBound(HighsBoundType::Upper, col, downval, -2);
            lp->setStoredBasis(std::move(basis));
            lp->recoverBasis();
            return -1;
          }
        }
      } else if (status == HighsLpRelaxation::Status::Infeasible) {
        addInfeasibleConflict();
        localdom.backtrack();
        lp->flushDomain(localdom);
        localdom.changeBound(HighsBoundType::Upper, col, downval, -2);
        lp->setStoredBasis(std::move(basis));
        lp->recoverBasis();
        return -1;
      } else {
        // printf("todo2\n");
        // in case of an LP error we set the score of this variable to zero to
        // avoid choosing it as branching candidate if possible
        downscore[candidate] = 0.0;
        upscore[candidate] = 0.0;
        downscorereliable[candidate] = 1;
        upscorereliable[candidate] = 1;
      }

      localdom.backtrack();
      lp->flushDomain(localdom);
      lp->recoverBasis();
    }
  }
}

const HighsSearch::NodeData* HighsSearch::getParentNodeData() const {
  if (nodestack.size() <= 1) return nullptr;

  return &nodestack[nodestack.size() - 2];
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

void HighsSearch::flushStatistics() {
  mipsolver.mipdata_->num_nodes += nnodes;
  nnodes = 0;

  mipsolver.mipdata_->pruned_treeweight += treeweight;
  treeweight = 0;

  mipsolver.mipdata_->total_lp_iterations += lpiterations;
  lpiterations = 0;

  mipsolver.mipdata_->heuristic_lp_iterations += heurlpiterations;
  heurlpiterations = 0;
}

size_t HighsSearch::getHeuristicLpIterations() const {
  return heurlpiterations + mipsolver.mipdata_->heuristic_lp_iterations;
}

size_t HighsSearch::getTotalLpIterations() const {
  return lpiterations + mipsolver.mipdata_->total_lp_iterations;
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
    if (lp->getMip().debugSolution_[i] + mipsolver.mipdata_->feastol <
            localdom.colLower_[i] ||
        lp->getMip().debugSolution_[i] - mipsolver.mipdata_->feastol >
            localdom.colUpper_[i]) {
      debugsolactive = false;
    }

    debugsolobj += lp->getMip().debugSolution_[i] * lp->getMip().colCost_[i];
  }
#endif
  localdom.propagate();
#ifdef HIGHS_DEBUGSOL
  if (debugsolactive &&
      upper_bound > debugsolobj + mipsolver.mipdata_->feastol) {
    bool debugsolstillactive = true;

    for (int i = 0; i != mipsolver.numCol(); ++i) {
      if (lp->getMip().debugSolution_[i] + mipsolver.mipdata_->feastol <
              localdom.colLower_[i] ||
          lp->getMip().debugSolution_[i] - mipsolver.mipdata_->feastol >
              localdom.colUpper_[i]) {
        debugsolstillactive = false;
        break;
      }
    }

    assert(debugsolstillactive);
  }
#endif

  NodeData& currnode = nodestack.back();

  bool prune = false;

  if (localdom.infeasible()) {
#ifdef HIGHS_DEBUGSOL
    assert(!debugsolactive ||
           upper_bound <= debugsolobj + mipsolver.mipdata_->feastol);
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
        } else if (!inbranching && !inheuristic) {
          if (getHeuristicLpIterations() <
                  getTotalLpIterations() *
                      mipsolver.mipdata_->heuristic_effort &&
              currnode.estimate < getCutoffBound()) {
            heuristicSearch();
          }
        }
      }

      if (lp->unscaledDualFeasible(status)) {
        currnode.lower_bound =
            std::max(lp->getObjective(), currnode.lower_bound);

        const NodeData* parent = getParentNodeData();

        if (parent != nullptr &&
            parent->branching_point != parent->branchingdecision.boundval) {
          int col = parent->branchingdecision.column;
          double delta =
              parent->branchingdecision.boundval - parent->branching_point;
          double objdelta =
              std::max(0.0, currnode.lower_bound - parent->lower_bound);

          pseudocost.addObservation(col, delta, objdelta);
        }

        if (currnode.lower_bound >= getCutoffBound()) {
#ifdef HIGHS_DEBUGSOL
          if (debugsolactive &&
              upper_bound > debugsolobj + mipsolver.mipdata_->feastol) {
            lp->getLpSolver().writeModel("wronglp->mps");
            lp->getLpSolver().writeBasis("wronglp->bas");
            assert(false);
          }
#endif
          addBoundExceedingConflict();
          prune = true;
        }
      }

    } else if (status == HighsLpRelaxation::Status::Infeasible) {
      addInfeasibleConflict();
#ifdef HIGHS_DEBUGSOL
      if (debugsolactive &&
          upper_bound > debugsolobj + mipsolver.mipdata_->feastol) {
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

          printf("glbminact: %g   localminact: %g   rhs: %g  len: %lu\n", glbminact, localminact, rhs, inds.size());
        }

      }
#endif

  if (prune) {
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
  pseudocost.setSeed(random.integer());

  while (currnode.opensubtrees == 2 && lp->scaledOptimal(lp->getStatus()) &&
         !lp->getFractionalIntegers().empty()) {
    // evalUnreliableBranchCands();
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

  assert(currnode.opensubtrees == 2 || currnode.opensubtrees == 0);

  if (currnode.opensubtrees != 2) return false;

  if (currnode.branchingdecision.column == -1) {
    double bestscore = -1.0;
    // solution branching failed, so choose any integer variable to branch
    // on in case we have a different solution status could happen due to a
    // fail in the LP solution process
    pseudocost.setSeed(random.integer());

    for (int i = 0; i != mipsolver.numCol(); ++i) {
      if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;

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