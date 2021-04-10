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

#include "lp_data/HConst.h"
#include "mip/HighsCutGeneration.h"
#include "mip/HighsDomainChange.h"
#include "mip/HighsMipSolverData.h"

HighsSearch::HighsSearch(HighsMipSolver& mipsolver,
                         const HighsPseudocost& pseudocost)
    : mipsolver(mipsolver),
      lp(nullptr),
      localdom(mipsolver.mipdata_->domain),
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
  childselrule = ChildSelectionRule::Disjunction;
  this->localdom.setDomainChangeStack(std::vector<HighsDomainChange>());
}

double HighsSearch::checkSol(const std::vector<double>& sol,
                             bool& integerfeasible) const {
  HighsCDouble objval = 0.0;
  integerfeasible = true;
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
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

double HighsSearch::getCutoffBound() const {
  return std::min(mipsolver.mipdata_->upper_limit, upper_limit);
}

void HighsSearch::setRINSNeighbourhood(const std::vector<double>& basesol,
                                       const std::vector<double>& relaxsol) {
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;
    if (localdom.colLower_[i] == localdom.colUpper_[i]) continue;

    double intval = std::floor(basesol[i] + 0.5);
    if (std::abs(relaxsol[i] - intval) < mipsolver.mipdata_->feastol) {
      if (localdom.colLower_[i] < intval)
        localdom.changeBound(HighsBoundType::Lower, i,
                             std::min(intval, localdom.colUpper_[i]),
                             HighsDomain::Reason::unspecified());
      if (localdom.colUpper_[i] > intval)
        localdom.changeBound(HighsBoundType::Upper, i,
                             std::max(intval, localdom.colLower_[i]),
                             HighsDomain::Reason::unspecified());
    }
  }
}

void HighsSearch::setRENSNeighbourhood(const std::vector<double>& lpsol) {
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
    if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;
    if (localdom.colLower_[i] == localdom.colUpper_[i]) continue;

    double downval = std::floor(lpsol[i] + mipsolver.mipdata_->feastol);
    double upval = std::ceil(lpsol[i] - mipsolver.mipdata_->feastol);

    if (localdom.colLower_[i] < downval) {
      localdom.changeBound(HighsBoundType::Lower, i,
                           std::min(downval, localdom.colUpper_[i]),
                           HighsDomain::Reason::unspecified());
      if (localdom.infeasible()) return;
    }
    if (localdom.colUpper_[i] > upval) {
      localdom.changeBound(HighsBoundType::Upper, i,
                           std::max(upval, localdom.colLower_[i]),
                           HighsDomain::Reason::unspecified());
      if (localdom.infeasible()) return;
    }
  }
}

void HighsSearch::createNewNode() {
  nodestack.emplace_back();
  nodestack.back().domgchgStackPos = localdom.getDomainChangeStack().size();
}

void HighsSearch::cutoffNode() { nodestack.back().opensubtrees = 0; }

void HighsSearch::setMinReliable(HighsInt minreliable) {
  pseudocost.setMinReliable(minreliable);
}

void HighsSearch::branchDownwards(HighsInt col, double newub,
                                  double branchpoint) {
  NodeData& currnode = nodestack.back();

  assert(currnode.opensubtrees == 2);
  assert(mipsolver.variableType(col) != HighsVarType::CONTINUOUS);

  currnode.opensubtrees = 1;
  currnode.branching_point = branchpoint;
  currnode.branchingdecision.column = col;
  currnode.branchingdecision.boundval = newub;
  currnode.branchingdecision.boundtype = HighsBoundType::Upper;

  HighsInt domchgPos = localdom.getDomainChangeStack().size();
  localdom.changeBound(currnode.branchingdecision);
  nodestack.emplace_back(currnode.lower_bound, currnode.estimate);
  nodestack.back().domgchgStackPos = domchgPos;
}

void HighsSearch::branchUpwards(HighsInt col, double newlb,
                                double branchpoint) {
  NodeData& currnode = nodestack.back();

  assert(currnode.opensubtrees == 2);
  assert(mipsolver.variableType(col) != HighsVarType::CONTINUOUS);

  currnode.opensubtrees = 1;
  currnode.branching_point = branchpoint;
  currnode.branchingdecision.column = col;
  currnode.branchingdecision.boundval = newlb;
  currnode.branchingdecision.boundtype = HighsBoundType::Lower;

  HighsInt domchgPos = localdom.getDomainChangeStack().size();
  localdom.changeBound(currnode.branchingdecision);
  nodestack.emplace_back(currnode.lower_bound, currnode.estimate);
  nodestack.back().domgchgStackPos = domchgPos;
}

void HighsSearch::addBoundExceedingConflict() {
  if (mipsolver.mipdata_->upper_limit != HIGHS_CONST_INF) {
    double rhs;
    if (lp->computeDualProof(mipsolver.mipdata_->domain,
                             mipsolver.mipdata_->upper_limit, inds, vals,
                             rhs)) {
      HighsCutGeneration cutGen(*lp, mipsolver.mipdata_->cutpool);
      mipsolver.mipdata_->debugSolution.checkCut(inds.data(), vals.data(),
                                                 inds.size(), rhs);
      cutGen.generateConflict(localdom, inds, vals, rhs);
    }
  }
}

void HighsSearch::addInfeasibleConflict() {
  double rhs;
  if (lp->computeDualInfProof(mipsolver.mipdata_->domain, inds, vals, rhs)) {
    // double minactlocal = 0.0;
    // double minactglobal = 0.0;
    // for (HighsInt i = 0; i < int(inds.size()); ++i) {
    //  if (vals[i] > 0.0) {
    //    minactlocal += localdom.colLower_[inds[i]] * vals[i];
    //    minactglobal += globaldom.colLower_[inds[i]] * vals[i];
    //  } else {
    //    minactlocal += localdom.colUpper_[inds[i]] * vals[i];
    //    minactglobal += globaldom.colUpper_[inds[i]] * vals[i];
    //  }
    //}
    // HighsInt oldnumcuts = cutpool.getNumCuts();
    HighsCutGeneration cutGen(*lp, mipsolver.mipdata_->cutpool);
    mipsolver.mipdata_->debugSolution.checkCut(inds.data(), vals.data(),
                                               inds.size(), rhs);
    cutGen.generateConflict(localdom, inds, vals, rhs);

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
    // HighsInt cutind = cutpool.addCut(inds.data(), vals.data(), inds.size(),
    // rhs); localdom.cutAdded(cutind);
  }
}

HighsInt HighsSearch::selectBranchingCandidate(int64_t maxSbIters) {
  assert(!lp->getFractionalIntegers().empty());

  static constexpr HighsInt basisstart_threshold = 20;
  std::vector<double> upscore;
  std::vector<double> downscore;
  std::vector<uint8_t> upscorereliable;
  std::vector<uint8_t> downscorereliable;

  HighsInt numfrac = lp->getFractionalIntegers().size();
  const auto& fracints = lp->getFractionalIntegers();

  upscore.resize(numfrac, HIGHS_CONST_INF);
  downscore.resize(numfrac, HIGHS_CONST_INF);

  upscorereliable.resize(numfrac, 0);
  downscorereliable.resize(numfrac, 0);

  // initialize up and down scores of variables that have a
  // reliable pseudocost so that they do not get evaluated
  for (HighsInt k = 0; k != numfrac; ++k) {
    HighsInt col = fracints[k].first;
    double fracval = fracints[k].second;

    assert(fracval > localdom.colLower_[col] + mipsolver.mipdata_->feastol);
    assert(fracval < localdom.colUpper_[col] - mipsolver.mipdata_->feastol);

    if (pseudocost.isReliable(col) || branchingVarReliableAtNode(col)) {
      upscore[k] = pseudocost.getPseudocostUp(col, fracval);
      downscore[k] = pseudocost.getPseudocostDown(col, fracval);
      upscorereliable[k] = true;
      downscorereliable[k] = true;
    }
  }

  std::vector<HighsInt> evalqueue;
  evalqueue.resize(numfrac);
  std::iota(evalqueue.begin(), evalqueue.end(), 0);

  auto numNodesUp = [&](HighsInt k) {
    if (mipsolver.mipdata_->domain.isBinary(fracints[k].first))
      return mipsolver.mipdata_->nodequeue.numNodesUp(fracints[k].first);

    return mipsolver.mipdata_->nodequeue.numNodesUp(fracints[k].first,
                                                    fracints[k].second);
  };

  auto numNodesDown = [&](HighsInt k) {
    if (mipsolver.mipdata_->domain.isBinary(fracints[k].first))
      return mipsolver.mipdata_->nodequeue.numNodesDown(fracints[k].first);

    return mipsolver.mipdata_->nodequeue.numNodesDown(fracints[k].first,
                                                      fracints[k].second);
  };

  double minScore = mipsolver.mipdata_->feastol;

  auto selectBestScore = [&](bool finalSelection) {
    HighsInt best = -1;
    double bestscore = -1.0;
    double bestnodes = -1.0;
    int64_t bestnumnodes = 0;

    for (HighsInt k : evalqueue) {
      double score;

      double s = 0.01 * std::min(upscorereliable[k] ? upscore[k] : 0,
                                 downscorereliable[k] ? downscore[k] : 0);
      minScore = std::max(s, minScore);

      if ((upscore[k] == 0.0 && upscorereliable[k]) ||
          (downscore[k] == 0.0 && downscorereliable[k]))
        score = pseudocost.getScore(fracints[k].first, 0.0, 0.0);
      else {
        score = upscore[k] == HIGHS_CONST_INF || downscore[k] == HIGHS_CONST_INF
                    ? finalSelection ? pseudocost.getScore(fracints[k].first,
                                                           fracints[k].second)
                                     : HIGHS_CONST_INF
                    : pseudocost.getScore(fracints[k].first, upscore[k],
                                          downscore[k]);
      }

      int64_t upnodes = numNodesUp(k);
      int64_t downnodes = numNodesDown(k);
      double nodes = 0;
      int64_t numnodes = upnodes + downnodes;
      if (upnodes != 0 || downnodes != 0)
        nodes =
            (downnodes / (double)(numnodes)) * (upnodes / (double)(numnodes));
      if (score > bestscore ||
          (score > bestscore - mipsolver.mipdata_->feastol &&
           std::make_pair(nodes, numnodes) >
               std::make_pair(bestnodes, bestnumnodes))) {
        bestscore = score;
        best = k;
        bestnodes = nodes;
        bestnumnodes = numnodes;
      }
    }

    return best;
  };

  lp->storeBasis();
  auto basis = lp->getStoredBasis();

  while (true) {
    bool mustStop = getStrongBranchingLpIterations() >= maxSbIters ||
                    mipsolver.mipdata_->checkLimits();

    HighsInt candidate = selectBestScore(mustStop);

    if ((upscorereliable[candidate] && downscorereliable[candidate]) ||
        mustStop) {
      lp->setStoredBasis(std::move(basis));
      lp->recoverBasis();
      lp->run();
      return candidate;
    }

    HighsInt col = fracints[candidate].first;
    double fracval = fracints[candidate].second;
    double upval = std::ceil(fracval);
    double downval = std::floor(fracval);

    if (!downscorereliable[candidate]) {
      // evaluate down branch
      int64_t inferences = -(int64_t)localdom.getDomainChangeStack().size() - 1;
      localdom.changeBound(HighsBoundType::Upper, col, downval);
      localdom.propagate();
      inferences += localdom.getDomainChangeStack().size();
      if (localdom.infeasible()) {
        pseudocost.addCutoffObservation(col, false);
        localdom.backtrack();
        localdom.clearChangedCols();
        localdom.changeBound(HighsBoundType::Lower, col, upval,
                             HighsDomain::Reason::unspecified());
        lp->setStoredBasis(std::move(basis));
        return -1;
      }

      pseudocost.addInferenceObservation(col, inferences, false);

      lp->flushDomain(localdom);

      int64_t numiters = lp->getNumLpIterations();
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
        if (objdelta < mipsolver.mipdata_->feastol) objdelta = 0.0;

        downscore[candidate] = objdelta;
        downscorereliable[candidate] = 1;
        markBranchingVarDownReliableAtNode(col);
        pseudocost.addObservation(col, delta, objdelta);

        for (HighsInt k = 0; k != numfrac; ++k) {
          double otherfracval = fracints[k].second;
          double otherdownval = std::floor(fracints[k].second);
          double otherupval = std::ceil(fracints[k].second);
          if (sol[fracints[k].first] <=
              otherdownval + mipsolver.mipdata_->feastol) {
            if (objdelta <= minScore && !downscorereliable[k])
              downscorereliable[k] = 1;

            downscore[k] = std::min(downscore[k], objdelta);
          } else if (sol[fracints[k].first] >=
                     otherupval - mipsolver.mipdata_->feastol) {
            if (objdelta <= minScore && !upscorereliable[k])
              upscorereliable[k] = 1;

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
            localdom.changeBound(HighsBoundType::Lower, col, upval,
                                 HighsDomain::Reason::unspecified());
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
            localdom.changeBound(HighsBoundType::Lower, col, upval,
                                 HighsDomain::Reason::unspecified());
            lp->setStoredBasis(std::move(basis));
            if (numiters > basisstart_threshold) lp->recoverBasis();
            return -1;
          }
        }
      } else if (status == HighsLpRelaxation::Status::Infeasible) {
        addInfeasibleConflict();
        pseudocost.addCutoffObservation(col, false);
        localdom.backtrack();
        lp->flushDomain(localdom);
        localdom.changeBound(HighsBoundType::Lower, col, upval,
                             HighsDomain::Reason::unspecified());
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
      int64_t inferences = -(int64_t)localdom.getDomainChangeStack().size() - 1;
      localdom.changeBound(HighsBoundType::Lower, col, upval);
      localdom.propagate();
      inferences += localdom.getDomainChangeStack().size();
      if (localdom.infeasible()) {
        pseudocost.addCutoffObservation(col, true);
        localdom.backtrack();
        localdom.clearChangedCols();
        localdom.changeBound(HighsBoundType::Upper, col, downval,
                             HighsDomain::Reason::unspecified());
        lp->setStoredBasis(std::move(basis));
        return -1;
      }

      pseudocost.addInferenceObservation(col, inferences, true);
      lp->flushDomain(localdom);

      int64_t numiters = lp->getNumLpIterations();
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

        for (HighsInt k = 0; k != numfrac; ++k) {
          double otherfracval = fracints[k].second;
          double otherdownval = std::floor(fracints[k].second);
          double otherupval = std::ceil(fracints[k].second);
          if (sol[fracints[k].first] <=
              otherdownval + mipsolver.mipdata_->feastol) {
            if (objdelta <= minScore && !downscorereliable[k])
              downscorereliable[k] = 1;

            downscore[k] = std::min(downscore[k], objdelta);

          } else if (sol[fracints[k].first] >=
                     otherupval - mipsolver.mipdata_->feastol) {
            if (objdelta <= minScore && !upscorereliable[k])
              upscorereliable[k] = 1;

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
            localdom.changeBound(HighsBoundType::Upper, col, downval,
                                 HighsDomain::Reason::unspecified());
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
            localdom.changeBound(HighsBoundType::Upper, col, downval,
                                 HighsDomain::Reason::unspecified());
            lp->setStoredBasis(std::move(basis));
            if (numiters > basisstart_threshold) lp->recoverBasis();
            return -1;
          }
        }
      } else if (status == HighsLpRelaxation::Status::Infeasible) {
        addInfeasibleConflict();
        pseudocost.addCutoffObservation(col, true);
        localdom.backtrack();
        lp->flushDomain(localdom);
        localdom.changeBound(HighsBoundType::Upper, col, downval,
                             HighsDomain::Reason::unspecified());
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
    } else {
      mipsolver.mipdata_->debugSolution.nodePruned(localdom);
      treeweight += std::pow(0.5, getCurrentDepth() - 1);
    }
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

  mipsolver.mipdata_->sb_lp_iterations += sblpiterations;
  sblpiterations = 0;
}

int64_t HighsSearch::getHeuristicLpIterations() const {
  return heurlpiterations + mipsolver.mipdata_->heuristic_lp_iterations;
}

int64_t HighsSearch::getTotalLpIterations() const {
  return lpiterations + mipsolver.mipdata_->total_lp_iterations;
}

int64_t HighsSearch::getLocalLpIterations() const { return lpiterations; }

int64_t HighsSearch::getStrongBranchingLpIterations() const {
  return sblpiterations + mipsolver.mipdata_->sb_lp_iterations;
}

void HighsSearch::resetLocalDomain() {
  this->lp->getLpSolver().changeColsBounds(
      0, mipsolver.numCol() - 1, mipsolver.mipdata_->domain.colLower_.data(),
      mipsolver.mipdata_->domain.colUpper_.data());
  localdom = mipsolver.mipdata_->domain;

#ifndef NDEBUG
  for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
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

HighsSearch::NodeResult HighsSearch::evaluateNode() {
  assert(!nodestack.empty());
  NodeData& currnode = nodestack.back();
  const NodeData* parent = getParentNodeData();

  const auto& domchgstack = localdom.getDomainChangeStack();

  localdom.propagate();

  if (parent != nullptr) {
    int64_t inferences = domchgstack.size() - (currnode.domgchgStackPos + 1);
    pseudocost.addInferenceObservation(
        parent->branchingdecision.column, inferences,
        parent->branchingdecision.boundtype == HighsBoundType::Lower);
  }

  NodeResult result = NodeResult::Open;

  if (localdom.infeasible()) {
    result = NodeResult::DomainInfeasible;
    localdom.clearChangedCols();
    if (parent != nullptr && parent->lp_objective != -HIGHS_CONST_INF &&
        parent->branching_point != parent->branchingdecision.boundval) {
      HighsInt col = parent->branchingdecision.column;
      bool upbranch =
          parent->branchingdecision.boundtype == HighsBoundType::Lower;
      pseudocost.addCutoffObservation(col, upbranch);
    }
  } else {
    lp->flushDomain(localdom);

#ifndef NDEBUG
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      assert(lp->getLpSolver().getLp().colLower_[i] == localdom.colLower_[i] ||
             mipsolver.variableType(i) == HighsVarType::CONTINUOUS);
      assert(lp->getLpSolver().getLp().colUpper_[i] == localdom.colUpper_[i] ||
             mipsolver.variableType(i) == HighsVarType::CONTINUOUS);
    }
#endif
    int64_t oldnumiters = lp->getNumLpIterations();
    HighsLpRelaxation::Status status = lp->resolveLp(&localdom);
    lpiterations += lp->getNumLpIterations() - oldnumiters;
    if (lp->scaledOptimal(status)) {
      lp->storeBasis();

      currnode.estimate = lp->computeBestEstimate(pseudocost);
      currnode.lp_objective = lp->getObjective();

      if (parent != nullptr && parent->lp_objective != -HIGHS_CONST_INF &&
          parent->branching_point != parent->branchingdecision.boundval) {
        HighsInt col = parent->branchingdecision.column;
        double delta =
            parent->branchingdecision.boundval - parent->branching_point;
        double objdelta =
            std::max(0.0, currnode.lp_objective - parent->lp_objective);

        pseudocost.addObservation(col, delta, objdelta);
      }

      if (lp->unscaledPrimalFeasible(status)) {
        if (lp->getFractionalIntegers().empty()) {
          result = NodeResult::BoundExceeding;
          double cutoffbnd = getCutoffBound();
          mipsolver.mipdata_->addIncumbent(
              lp->getLpSolver().getSolution().col_value, lp->getObjective(),
              inheuristic ? 'H' : 'T');
          if (mipsolver.mipdata_->upper_limit < cutoffbnd)
            lp->setObjectiveLimit(mipsolver.mipdata_->upper_limit);
          addBoundExceedingConflict();
        }
      }

      if (result == NodeResult::Open) {
        if (lp->unscaledDualFeasible(status)) {
          currnode.lower_bound =
              std::max(currnode.lp_objective, currnode.lower_bound);

          if (currnode.lower_bound > getCutoffBound()) {
            result = NodeResult::BoundExceeding;
            addBoundExceedingConflict();
          } else if (mipsolver.mipdata_->upper_limit != HIGHS_CONST_INF) {
            HighsRedcostFixing::propagateRedCost(
                mipsolver, localdom, lp->getLpSolver().getSolution().col_dual,
                lp->getObjective());
            if (localdom.infeasible()) {
              result = NodeResult::BoundExceeding;
              addBoundExceedingConflict();
              localdom.clearChangedCols();
            } else if (!localdom.getChangedCols().empty()) {
              return evaluateNode();
            }
          }
        } else if (lp->getObjective() > getCutoffBound()) {
          // the LP is not solved to dual feasibilty due to scaling/numerics
          // therefore we compute a conflict constraint as if the LP was bound
          // exceeding and propagate the local domain again. The lp relaxation
          // class will take care to consider the dual multipliers with an
          // increased zero tolerance due to the dual infeasibility when
          // computing the proof conBoundExceedingstraint.
          addBoundExceedingConflict();
          localdom.propagate();
          if (localdom.infeasible()) result = NodeResult::BoundExceeding;
        }
      }
    } else if (status == HighsLpRelaxation::Status::Infeasible) {
      result = NodeResult::LpInfeasible;
      addInfeasibleConflict();
      if (parent != nullptr && parent->lp_objective != -HIGHS_CONST_INF &&
          parent->branching_point != parent->branchingdecision.boundval) {
        HighsInt col = parent->branchingdecision.column;
        bool upbranch =
            parent->branchingdecision.boundtype == HighsBoundType::Lower;
        pseudocost.addCutoffObservation(col, upbranch);
      }
    }
  }

  if (result != NodeResult::Open) {
    mipsolver.mipdata_->debugSolution.nodePruned(localdom);
    treeweight += std::pow(0.5, getCurrentDepth() - 1);
    currnode.opensubtrees = 0;
  }

  return result;
}

HighsSearch::NodeResult HighsSearch::branch() {
  assert(localdom.getChangedCols().empty());

  NodeData& currnode = nodestack.back();
  assert(currnode.opensubtrees == 2);
  currnode.branchingdecision.column = -1;
  inbranching = true;

  HighsInt minrel = pseudocost.getMinReliable();

  NodeResult result = NodeResult::Open;
  while (currnode.opensubtrees == 2 && lp->scaledOptimal(lp->getStatus()) &&
         !lp->getFractionalIntegers().empty()) {
    int64_t sbmaxiters = 0;
    if (minrel > 0) {
      int64_t sbiters = getStrongBranchingLpIterations();
      sbmaxiters =
          100000 + ((getTotalLpIterations() - getHeuristicLpIterations() -
                     getStrongBranchingLpIterations()) >>
                    1);
      if (sbiters > sbmaxiters) {
        pseudocost.setMinReliable(0);
      } else if (sbiters > sbmaxiters / 2) {
        double reductionratio =
            (sbiters - sbmaxiters / 2) / (double)(sbmaxiters - sbmaxiters / 2);

        HighsInt minrelreduced = int(minrel - reductionratio * (minrel - 1));
        pseudocost.setMinReliable(std::min(minrel, minrelreduced));
      }
    }

    HighsInt branchcand = selectBranchingCandidate(sbmaxiters);

    if (branchcand != -1) {
      auto branching = lp->getFractionalIntegers()[branchcand];
      currnode.branchingdecision.column = branching.first;
      currnode.branching_point = branching.second;

      HighsInt col = branching.first;
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
          if (random.bit()) {
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
        case ChildSelectionRule::Disjunction: {
          int64_t numnodesup;
          int64_t numnodesdown;
          if (mipsolver.mipdata_->domain.isBinary(col)) {
            // this is faster for binary variables as the overload that gets a
            // value may be linear in the number of open nodes in the worst case
            // (if every node has tightened local bounds on this column)
            numnodesup = mipsolver.mipdata_->nodequeue.numNodesUp(col);
            numnodesdown = mipsolver.mipdata_->nodequeue.numNodesDown(col);
          } else {
            numnodesup = mipsolver.mipdata_->nodequeue.numNodesUp(
                col, currnode.branching_point);
            numnodesdown = mipsolver.mipdata_->nodequeue.numNodesDown(
                col, currnode.branching_point);
          }
          if (numnodesup > numnodesdown) {  // > -> neos*-inde sehr schnell
            currnode.branchingdecision.boundtype = HighsBoundType::Lower;
            currnode.branchingdecision.boundval =
                std::ceil(currnode.branching_point);
          } else if (numnodesdown > numnodesup) {
            currnode.branchingdecision.boundtype = HighsBoundType::Upper;
            currnode.branchingdecision.boundval =
                std::floor(currnode.branching_point);
          } else {
            if (mipsolver.colCost(col) >= 0) {
              currnode.branchingdecision.boundtype = HighsBoundType::Lower;
              currnode.branchingdecision.boundval =
                  std::ceil(currnode.branching_point);
            } else {
              currnode.branchingdecision.boundtype = HighsBoundType::Upper;
              currnode.branchingdecision.boundval =
                  std::floor(currnode.branching_point);
            }
          }
          break;
        }
      }
      result = NodeResult::Branched;
      break;
    }

    assert(!localdom.getChangedCols().empty());
    result = evaluateNode();
  }
  inbranching = false;
  pseudocost.setMinReliable(minrel);

  assert(currnode.opensubtrees == 2 || currnode.opensubtrees == 0);

  if (currnode.opensubtrees != 2) return result;

  if (currnode.branchingdecision.column == -1) {
    double bestscore = -1.0;
    // solution branching failed, so choose any integer variable to branch
    // on in case we have a different solution status could happen due to a
    // fail in the LP solution process

    for (HighsInt i : mipsolver.mipdata_->integral_cols) {
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
    lp->setIterationLimit();

    // create a fresh LP only with model rows since all integer columns are
    // fixed, the cutting planes are not required and the LP could not be solved
    // so we want to make it as easy as possible
    HighsLpRelaxation lpCopy(mipsolver);
    lpCopy.loadModel();
    lpCopy.getLpSolver().changeColsBounds(0, mipsolver.numCol() - 1,
                                          localdom.colLower_.data(),
                                          localdom.colUpper_.data());
    // temporarily use the fresh LP for the HighsSearch class
    HighsLpRelaxation* tmpLp = &lpCopy;
    std::swap(tmpLp, lp);

    // reevaluate the node with LP presolve enabled
    lp->getLpSolver().setHighsOptionValue("presolve", "on");
    result = evaluateNode();

    if (result == NodeResult::Open) {
      // LP still not solved, reevaluate with primal simplex
      lp->getLpSolver().clearSolver();
      lp->getLpSolver().setHighsOptionValue("simplex_strategy",
                                            SIMPLEX_STRATEGY_PRIMAL);
      result = evaluateNode();
      lp->getLpSolver().setHighsOptionValue("simplex_strategy",
                                            SIMPLEX_STRATEGY_DUAL);
      if (result == NodeResult::Open) {
        // LP still not solved, reevaluate with IPM instead of simplex
        lp->getLpSolver().clearSolver();
        lp->getLpSolver().setHighsOptionValue("solver", "ipm");
        result = evaluateNode();

        if (result == NodeResult::Open) {
          highsLogUser(mipsolver.options_mip_->log_options,
                       HighsLogType::WARNING,
                       "Failed to solve node with all integer columns "
                       "fixed. Declaring node infeasible.\n");
          // LP still not solved, give up and declare as infeasible
          currnode.opensubtrees = 0;
          result = NodeResult::LpInfeasible;
        }
      }
    }

    // restore old lp relaxation
    std::swap(tmpLp, lp);

    return result;
  }

  // finally open a new node with the branching decision added
  // and remember that we have one open subtree left
  HighsInt domchgPos = localdom.getDomainChangeStack().size();
  localdom.changeBound(currnode.branchingdecision);
  currnode.opensubtrees = 1;
  nodestack.emplace_back(currnode.lower_bound, currnode.estimate);
  nodestack.back().domgchgStackPos = domchgPos;

  return NodeResult::Branched;
}

bool HighsSearch::backtrack() {
  if (nodestack.empty()) return false;
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

  HighsInt domchgPos = localdom.getDomainChangeStack().size();
  localdom.changeBound(currnode.branchingdecision);
  nodestack.emplace_back(currnode.lower_bound, currnode.estimate);
  lp->flushDomain(localdom);
  nodestack.back().domgchgStackPos = domchgPos;

  return true;
}

bool HighsSearch::backtrackUntilDepth(HighsInt targetDepth) {
  if (nodestack.empty()) return false;
  assert(!nodestack.empty());
  if (getCurrentDepth() >= targetDepth) nodestack.back().opensubtrees = 0;

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

    if (getCurrentDepth() >= targetDepth) nodestack.back().opensubtrees = 0;
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

  HighsInt domchgPos = localdom.getDomainChangeStack().size();
  localdom.changeBound(currnode.branchingdecision);
  nodestack.emplace_back(currnode.lower_bound, currnode.estimate);
  lp->flushDomain(localdom);
  nodestack.back().domgchgStackPos = domchgPos;

  return true;
}

HighsSearch::NodeResult HighsSearch::dive() {
  reliableatnode.clear();

  do {
    ++nnodes;
    NodeResult result = evaluateNode();

    if (mipsolver.mipdata_->checkLimits()) return result;

    if (result != NodeResult::Open) return result;

    result = branch();
    if (result != NodeResult::Branched) return result;
  } while (true);
}

void HighsSearch::solveDepthFirst(int64_t maxbacktracks) {
  do {
    if (maxbacktracks == 0) break;

    dive();

    --maxbacktracks;

  } while (backtrack());
}
