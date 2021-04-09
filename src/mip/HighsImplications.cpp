/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsImplications.h"

#include "mip/HighsCliqueTable.h"
#include "mip/HighsMipSolverData.h"

bool HighsImplications::computeImplications(HighsInt col, bool val) {
  HighsDomain& globaldomain = mipsolver.mipdata_->domain;
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  globaldomain.propagate();
  if (globaldomain.infeasible() || !globaldomain.isBinary(col)) return true;
  const auto& domchgstack = globaldomain.getDomainChangeStack();
  HighsInt changedend = globaldomain.getChangedCols().size();

  HighsInt stackimplicstart = domchgstack.size() + 1;
  HighsInt numImplications = -stackimplicstart;
  if (val)
    globaldomain.changeBound(HighsBoundType::Lower, col, 1);
  else
    globaldomain.changeBound(HighsBoundType::Upper, col, 0);

  if (globaldomain.infeasible()) {
    globaldomain.backtrack();
    globaldomain.clearChangedCols(changedend);
    cliquetable.vertexInfeasible(globaldomain, col, val);

    return true;
  }

  globaldomain.propagate();

  if (globaldomain.infeasible()) {
    globaldomain.backtrack();
    globaldomain.clearChangedCols(changedend);

    cliquetable.vertexInfeasible(globaldomain, col, val);

    return true;
  }

  numImplications += domchgstack.size();
  mipsolver.mipdata_->pseudocost.addInferenceObservation(col, numImplications,
                                                         val);
  HighsInt stackimplicend = domchgstack.size();

  HighsInt loc = 2 * col + val;
  HighsInt implstart = implications.size();

  implications.insert(implications.end(), domchgstack.data() + stackimplicstart,
                      domchgstack.data() + stackimplicend);

  globaldomain.backtrack();
  globaldomain.clearChangedCols(changedend);

  // add the implications of binary variables to the clique table
  auto binstart =
      std::partition(implications.begin() + implstart, implications.end(),
                     [&](const HighsDomainChange& a) {
                       return !globaldomain.isBinary(a.column);
                     });

  std::sort(implications.begin() + implstart, binstart);

  HighsCliqueTable::CliqueVar clique[2];
  clique[0] = HighsCliqueTable::CliqueVar(col, val);

  for (auto i = binstart; i != implications.end(); ++i) {
    if (i->boundtype == HighsBoundType::Lower)
      clique[1] = HighsCliqueTable::CliqueVar(i->column, 0);
    else
      clique[1] = HighsCliqueTable::CliqueVar(i->column, 1);

    cliquetable.addClique(mipsolver, clique, 2);
    if (globaldomain.infeasible() || globaldomain.isFixed(col)) return true;
  }

  // store variable bounds derived from implications
  for (auto i = implications.begin() + implstart; i != binstart; ++i) {
    if (i->boundtype == HighsBoundType::Lower) {
      if (val == 1) {
        if (globaldomain.colLower_[i->column] != -HIGHS_CONST_INF)
          addVLB(i->column, col,
                 i->boundval - globaldomain.colLower_[i->column],
                 globaldomain.colLower_[i->column]);
      } else
        addVLB(i->column,
               col,  // in case the lower bound is infinite the varbound can
                     // still be tightened as soon as a finite upper bound is
                     // known because the offset is finite
               globaldomain.colLower_[i->column] - i->boundval, i->boundval);
    } else {
      if (val == 1) {
        if (globaldomain.colUpper_[i->column] != HIGHS_CONST_INF)
          addVUB(i->column, col,
                 i->boundval - globaldomain.colUpper_[i->column],
                 globaldomain.colUpper_[i->column]);
      } else
        addVUB(i->column,
               col,  // in case the upper bound is infinite the varbound can
                     // still be tightened as soon as a finite upper bound is
                     // known because the offset is finite
               globaldomain.colUpper_[i->column] - i->boundval, i->boundval);
    }
  }

  implications.erase(binstart, implications.end());

  implicationmap[loc].start = implstart;
  implicationmap[loc].num = implications.size() - implstart;

  return false;
}

bool HighsImplications::runProbing(HighsInt col, HighsInt& numboundchgs) {
  HighsDomain& globaldomain = mipsolver.mipdata_->domain;
  if (globaldomain.isBinary(col) && !implicationsCached(col, 1) &&
      !implicationsCached(col, 0) &&
      mipsolver.mipdata_->cliquetable.getSubstitution(col) == nullptr) {
    bool infeasible;

    infeasible = computeImplications(col, 1);
    if (globaldomain.infeasible()) return true;
    if (infeasible) return true;
    if (mipsolver.mipdata_->cliquetable.getSubstitution(col) != nullptr)
      return true;

    infeasible = computeImplications(col, 0);
    if (globaldomain.infeasible()) return true;
    if (infeasible) return true;
    if (mipsolver.mipdata_->cliquetable.getSubstitution(col) != nullptr)
      return true;

    // analyze implications
    const HighsDomainChange* implicsup;
    const HighsDomainChange* implicsdown;
    HighsInt nimplicsup;
    HighsInt nimplicsdown;
    nimplicsdown = getImplications(col, 0, implicsdown, infeasible);
    nimplicsup = getImplications(col, 1, implicsup, infeasible);
    HighsInt u = 0;
    HighsInt d = 0;

    while (u < nimplicsup && d < nimplicsdown) {
      if (implicsup[u].column < implicsdown[d].column)
        ++u;
      else if (implicsdown[d].column < implicsup[u].column)
        ++d;
      else {
        assert(implicsup[u].column == implicsdown[d].column);

        if (implicsup[u].boundtype == implicsdown[d].boundtype) {
          if ((implicsup[u].boundtype == HighsBoundType::Lower &&
               implicsdown[d].boundval < implicsup[u].boundval) ||
              (implicsup[u].boundtype == HighsBoundType::Upper &&
               implicsdown[d].boundval > implicsup[u].boundval))
            globaldomain.changeBound(implicsdown[d],
                                     HighsDomain::Reason::unspecified());
          else
            globaldomain.changeBound(implicsup[u],
                                     HighsDomain::Reason::unspecified());
          assert(!globaldomain.infeasible());
          ++numboundchgs;
          globaldomain.propagate();
          assert(!globaldomain.infeasible());
          ++u;
          ++d;
        } else if (!globaldomain.isFixed(implicsup[u].column) &&
                   !colsubstituted[implicsup[u].column] &&
                   globaldomain.isFixing(implicsup[u]) &&
                   globaldomain.isFixing(implicsdown[d])) {
          HighsSubstitution substitution;
          substitution.substcol = implicsup[u].column;
          substitution.staycol = col;
          substitution.offset = implicsdown[d].boundval;
          substitution.scale = implicsup[u].boundval - implicsdown[d].boundval;
          substitutions.push_back(substitution);
          colsubstituted[implicsup[u].column] = true;

        } else if ((HighsInt)implicsup[u].boundtype <
                   (HighsInt)implicsdown[d].boundtype)
          ++u;
        else
          ++d;
      }
    }

    return true;
  }

  return false;
}

void HighsImplications::addVUB(HighsInt col, HighsInt vubcol, double vubcoef,
                               double vubconstant) {
  VarBound vub{vubcoef, vubconstant};

  mipsolver.mipdata_->debugSolution.checkVub(col, vubcol, vubcoef, vubconstant);

  double minBound = vub.minValue();
  if (minBound >=
      mipsolver.mipdata_->domain.colUpper_[col] - mipsolver.mipdata_->feastol)
    return;

  auto insertresult = vubs[col].emplace(vubcol, vub);

  if (!insertresult.second) {
    VarBound& currentvub = insertresult.first->second;
    double currentMinBound = currentvub.minValue();
    if (minBound < currentMinBound - mipsolver.mipdata_->feastol) {
      currentvub.coef = vubcoef;
      currentvub.constant = vubconstant;
    }
  }
}

void HighsImplications::addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef,
                               double vlbconstant) {
  VarBound vlb{vlbcoef, vlbconstant};

  mipsolver.mipdata_->debugSolution.checkVlb(col, vlbcol, vlbcoef, vlbconstant);

  double maxBound = vlb.maxValue();
  if (vlb.maxValue() <=
      mipsolver.mipdata_->domain.colLower_[col] + mipsolver.mipdata_->feastol)
    return;

  auto insertresult = vlbs[col].emplace(vlbcol, vlb);

  if (!insertresult.second) {
    VarBound& currentvlb = insertresult.first->second;

    double currentMaxNound = currentvlb.maxValue();
    if (maxBound > currentMaxNound + mipsolver.mipdata_->feastol) {
      currentvlb.coef = vlbcoef;
      currentvlb.constant = vlbconstant;
    }
  }
}

void HighsImplications::rebuild(HighsInt ncols,
                                const std::vector<HighsInt>& orig2reducedcol,
                                const std::vector<HighsInt>& orig2reducedrow) {
  std::vector<std::map<HighsInt, VarBound>> oldvubs;
  std::vector<std::map<HighsInt, VarBound>> oldvlbs;

  oldvlbs.swap(vlbs);
  oldvubs.swap(vubs);

  colsubstituted.clear();
  colsubstituted.shrink_to_fit();
  implicationmap.clear();
  implicationmap.shrink_to_fit();

  implicationmap.resize(2 * ncols, {-1, 0});
  colsubstituted.resize(ncols);
  substitutions.clear();
  vubs.clear();
  vubs.shrink_to_fit();
  vubs.resize(ncols);
  vlbs.clear();
  vlbs.shrink_to_fit();
  vlbs.resize(ncols);
  HighsInt oldncols = oldvubs.size();

  for (HighsInt i = 0; i != oldncols; ++i) {
    HighsInt newi = orig2reducedcol[i];

    if (newi == -1) continue;

    for (const auto& oldvub : oldvubs[i]) {
      if (orig2reducedcol[oldvub.first] == -1) continue;

      if (!mipsolver.mipdata_->domain.isBinary(orig2reducedcol[oldvub.first]))
        continue;
      addVUB(newi, orig2reducedcol[oldvub.first], oldvub.second.coef,
             oldvub.second.constant);
    }

    for (const auto& oldvlb : oldvlbs[i]) {
      if (orig2reducedcol[oldvlb.first] == -1) continue;

      if (!mipsolver.mipdata_->domain.isBinary(orig2reducedcol[oldvlb.first]))
        continue;
      addVLB(newi, orig2reducedcol[oldvlb.first], oldvlb.second.coef,
             oldvlb.second.constant);
    }

    // todo also add old implications once implications can be added
    // incrementally for now we discard the old implications as they might be
    // weaker then newly computed ones and adding them would block computation
    // of new implications
  }
}

void HighsImplications::buildFrom(const HighsImplications& init) {
  return;
  HighsInt numcol = mipsolver.numCol();

  for (HighsInt i = 0; i != numcol; ++i) {
    for (const auto& vub : init.vubs[i]) {
      if (!mipsolver.mipdata_->domain.isBinary(vub.first)) continue;
      addVUB(i, vub.first, vub.second.coef, vub.second.constant);
    }

    for (const auto& vlb : init.vlbs[i]) {
      if (!mipsolver.mipdata_->domain.isBinary(vlb.first)) continue;
      addVLB(i, vlb.first, vlb.second.coef, vlb.second.constant);
    }

    // todo also add old implications once implications can be added
    // incrementally for now we discard the old implications as they might be
    // weaker then newly computed ones and adding them would block computation
    // of new implications
  }
}

void HighsImplications::separateImpliedBounds(
    const HighsLpRelaxation& lpRelaxation, const std::vector<double>& sol,
    HighsCutPool& cutpool, double feastol) {
  HighsDomain& globaldomain = mipsolver.mipdata_->domain;

  HighsInt inds[2];
  double vals[2];
  double rhs;

  HighsInt numboundchgs = 0;

  // first do probing on all candidates that have not been probed yet
  for (std::pair<HighsInt, double> fracint :
       lpRelaxation.getFractionalIntegers()) {
    HighsInt col = fracint.first;
    if (globaldomain.colLower_[col] != 0.0 ||
        globaldomain.colUpper_[col] != 1.0)
      continue;

    if (runProbing(col, numboundchgs)) {
      if (globaldomain.infeasible()) return;
    }
  }

  for (std::pair<HighsInt, double> fracint :
       lpRelaxation.getFractionalIntegers()) {
    HighsInt col = fracint.first;
    // skip non binary variables
    if (globaldomain.colLower_[col] != 0.0 ||
        globaldomain.colUpper_[col] != 1.0)
      continue;

    bool infeas;
    const HighsDomainChange* implics = nullptr;

    HighsInt nimplics = getImplications(col, 1, implics, infeas);
    if (globaldomain.infeasible()) return;

    if (infeas) {
      vals[0] = 1.0;
      inds[0] = col;
      cutpool.addCut(mipsolver, inds, vals, 1, 0.0, false, false);
      continue;
    }

    for (HighsInt i = 0; i != nimplics; ++i) {
      if (implics[i].boundtype == HighsBoundType::Upper) {
        if (implics[i].boundval + feastol >=
            globaldomain.colUpper_[implics[i].column])
          continue;

        vals[0] = 1.0;
        inds[0] = implics[i].column;
        vals[1] =
            globaldomain.colUpper_[implics[i].column] - implics[i].boundval;
        inds[1] = col;
        rhs = globaldomain.colUpper_[implics[i].column];

      } else {
        if (implics[i].boundval - feastol <=
            globaldomain.colLower_[implics[i].column])
          continue;

        vals[0] = -1.0;
        inds[0] = implics[i].column;
        vals[1] =
            globaldomain.colLower_[implics[i].column] - implics[i].boundval;
        inds[1] = col;
        rhs = -globaldomain.colLower_[implics[i].column];
      }

      double viol = sol[inds[0]] * vals[0] + sol[inds[1]] * vals[1] - rhs;

      if (viol > feastol) {
        // printf("added implied bound cut to pool\n");
        cutpool.addCut(mipsolver, inds, vals, 2, rhs,
                       mipsolver.variableType(implics[i].column) !=
                           HighsVarType::CONTINUOUS,
                       false);
      }
    }

    nimplics = getImplications(col, 0, implics, infeas);
    if (globaldomain.infeasible()) return;

    if (infeas) {
      vals[0] = -1.0;
      inds[0] = col;
      cutpool.addCut(mipsolver, inds, vals, 1, -1.0, false, false);
      continue;
    }

    for (HighsInt i = 0; i != nimplics; ++i) {
      if (implics[i].boundtype == HighsBoundType::Upper) {
        if (implics[i].boundval + feastol >=
            globaldomain.colUpper_[implics[i].column])
          continue;

        vals[0] = 1.0;
        inds[0] = implics[i].column;
        vals[1] =
            implics[i].boundval - globaldomain.colUpper_[implics[i].column];
        inds[1] = col;
        rhs = implics[i].boundval;
      } else {
        if (implics[i].boundval - feastol <=
            globaldomain.colLower_[implics[i].column])
          continue;

        vals[0] = -1.0;
        inds[0] = implics[i].column;
        vals[1] =
            globaldomain.colLower_[implics[i].column] - implics[i].boundval;
        inds[1] = col;
        rhs = -implics[i].boundval;
      }

      double viol = sol[inds[0]] * vals[0] + sol[inds[1]] * vals[1] - rhs;

      if (viol > feastol) {
        // printf("added implied bound cut to pool\n");
        cutpool.addCut(mipsolver, inds, vals, 2, rhs,
                       mipsolver.variableType(implics[i].column) !=
                           HighsVarType::CONTINUOUS,
                       false);
      }
    }
  }
}

void HighsImplications::cleanupVarbounds(HighsInt col) {
  double ub = mipsolver.mipdata_->domain.colUpper_[col];
  double lb = mipsolver.mipdata_->domain.colLower_[col];

  if (ub == lb) {
    vlbs[col].clear();
    vubs[col].clear();
    return;
  }

  auto next = vubs[col].begin();
  while (next != vubs[col].end()) {
    auto it = next++;

    mipsolver.mipdata_->debugSolution.checkVub(col, it->first, it->second.coef,
                                               it->second.constant);

    if (it->second.coef > 0) {
      double minub = it->second.constant;
      double maxub = it->second.constant + it->second.coef;
      if (minub >= ub - mipsolver.mipdata_->feastol)
        vubs[col].erase(it);  // variable bound is redundant
      else if (maxub > ub + mipsolver.mipdata_->epsilon) {
        it->second.coef =
            ub - it->second.constant;  // coefficient can be tightened
        mipsolver.mipdata_->debugSolution.checkVub(
            col, it->first, it->second.coef, it->second.constant);
      } else if (maxub < ub - mipsolver.mipdata_->epsilon) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::Upper, col, maxub,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }
    } else {
      HighsCDouble minub = HighsCDouble(it->second.constant) + it->second.coef;
      double maxub = it->second.constant;
      if (minub >= ub - mipsolver.mipdata_->feastol)
        vubs[col].erase(it);  // variable bound is redundant
      else if (maxub > ub + mipsolver.mipdata_->epsilon) {
        // variable bound can be tightened
        it->second.constant = ub;
        it->second.coef = double(minub - ub);
        mipsolver.mipdata_->debugSolution.checkVub(
            col, it->first, it->second.coef, it->second.constant);
      } else if (maxub < ub - mipsolver.mipdata_->epsilon) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::Upper, col, maxub,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }
    }
  }

  next = vlbs[col].begin();
  while (next != vlbs[col].end()) {
    auto it = next++;

    mipsolver.mipdata_->debugSolution.checkVlb(col, it->first, it->second.coef,
                                               it->second.constant);

    if (it->second.coef > 0) {
      HighsCDouble maxlb = HighsCDouble(it->second.constant) + it->second.coef;
      double minlb = it->second.constant;
      if (maxlb <= lb + mipsolver.mipdata_->feastol)
        vlbs[col].erase(it);  // variable bound is redundant
      else if (minlb < lb - mipsolver.mipdata_->epsilon) {
        // variable bound can be tightened
        it->second.constant = lb;
        it->second.coef = double(maxlb - lb);
        mipsolver.mipdata_->debugSolution.checkVlb(
            col, it->first, it->second.coef, it->second.constant);
      } else if (minlb > lb + mipsolver.mipdata_->epsilon) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::Lower, col, minlb,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }

    } else {
      double maxlb = it->second.constant;
      double minlb = it->second.constant + it->second.coef;
      if (maxlb <= lb + mipsolver.mipdata_->feastol)
        vlbs[col].erase(it);  // variable bound is redundant
      else if (minlb < lb - mipsolver.mipdata_->epsilon) {
        it->second.coef =
            lb - it->second.constant;  // variable bound can be tightened
        mipsolver.mipdata_->debugSolution.checkVlb(
            col, it->first, it->second.coef, it->second.constant);
      } else if (minlb > lb + mipsolver.mipdata_->epsilon) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::Lower, col, minlb,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }
    }
  }
}
