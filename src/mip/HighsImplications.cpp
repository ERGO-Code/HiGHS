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

bool HighsImplications::computeImplications(int col, bool val) {
  HighsDomain& globaldomain = mipsolver.mipdata_->domain;
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  globaldomain.propagate();
  if (globaldomain.infeasible()) return true;
  const auto& domchgstack = globaldomain.getDomainChangeStack();

  int changedend = globaldomain.getChangedCols().size();

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

  cliquetable.addImplications(globaldomain, col, val);

  size_t stackimplicstart = domchgstack.size();

  globaldomain.propagate();

  if (globaldomain.infeasible()) {
    globaldomain.backtrack();
    globaldomain.clearChangedCols(changedend);

    cliquetable.vertexInfeasible(globaldomain, col, val);

    return true;
  }

  size_t stackimplicend = domchgstack.size();

  int loc = 2 * col + val;
  int implstart = implications.size();

  implications.insert(implications.end(), &domchgstack[stackimplicstart],
                      &domchgstack[stackimplicend]);

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
      if (val == 1)
        addVLB(i->column, col, i->boundval - globaldomain.colLower_[i->column],
               globaldomain.colLower_[i->column]);
      else
        addVLB(i->column, col, globaldomain.colLower_[i->column] - i->boundval,
               i->boundval);
    } else {
      if (val == 1)
        addVUB(i->column, col, i->boundval - globaldomain.colUpper_[i->column],
               globaldomain.colUpper_[i->column]);
      else
        addVUB(i->column, col, globaldomain.colUpper_[i->column] - i->boundval,
               i->boundval);
    }
  }

  implications.erase(binstart, implications.end());

  implicationmap[loc].start = implstart;
  implicationmap[loc].num = implications.size() - implstart;

  return false;
}

bool HighsImplications::runProbing(int col, int& numboundchgs) {
  HighsDomain& globaldomain = mipsolver.mipdata_->domain;
  if (globaldomain.isBinary(col) && !implicationsCached(col, 1) &&
      !implicationsCached(col, 0) &&
      mipsolver.mipdata_->cliquetable.getSubstitution(col) == nullptr) {
    bool infeasible;

    infeasible = computeImplications(col, 1);
    if (globaldomain.infeasible()) return true;
    if (infeasible) return true;

    infeasible = computeImplications(col, 0);
    if (globaldomain.infeasible()) return true;
    if (infeasible) return true;

    // analyze implications
    const HighsDomainChange* implicsup;
    const HighsDomainChange* implicsdown;
    int nimplicsup;
    int nimplicsdown;
    nimplicsdown = getImplications(col, 0, implicsdown, infeasible);
    nimplicsup = getImplications(col, 1, implicsup, infeasible);
    int u = 0;
    int d = 0;

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

        } else if ((int)implicsup[u].boundtype < (int)implicsdown[d].boundtype)
          ++u;
        else
          ++d;
      }
    }
  }

  return false;
}

void HighsImplications::addVUB(int col, int vubcol, double vubcoef,
                               double vubconstant) {
  VarBound vub{vubcoef, vubconstant};

  if (vub.minValue() >=
      mipsolver.mipdata_->domain.colUpper_[col] - mipsolver.mipdata_->feastol)
    return;

  auto insertresult = vubs[col].emplace(vubcol, vub);

  if (!insertresult.second) {
    VarBound& currentvub = insertresult.first->second;
    double minbound_current = currentvub.constant - std::abs(currentvub.coef);
    double minbound_new = vubconstant - std::abs(vubcoef);
    if (minbound_new < minbound_current - mipsolver.mipdata_->feastol) {
      currentvub.coef = vubcoef;
      currentvub.constant = vubconstant;
    }
  }
}

void HighsImplications::addVLB(int col, int vlbcol, double vlbcoef,
                               double vlbconstant) {
  VarBound vlb{vlbcoef, vlbconstant};

  if (vlb.maxValue() <=
      mipsolver.mipdata_->domain.colLower_[col] + mipsolver.mipdata_->feastol)
    return;

  auto insertresult = vlbs[col].emplace(vlbcol, vlb);

  if (!insertresult.second) {
    VarBound& currentvlb = insertresult.first->second;
    double maxbound_current = currentvlb.constant + std::abs(currentvlb.coef);
    double maxbound_new = vlbconstant + std::abs(vlbcoef);
    if (maxbound_new > maxbound_current + 1e-6) {
      currentvlb.coef = vlbcoef;
      currentvlb.constant = vlbconstant;
    }
  }
}

void HighsImplications::rebuild(int ncols,
                                const std::vector<int>& orig2reducedcol,
                                const std::vector<int>& orig2reducedrow) {
  std::vector<std::map<int, VarBound>> oldvubs;
  std::vector<std::map<int, VarBound>> oldvlbs;

  oldvlbs.swap(vlbs);
  oldvubs.swap(vubs);

  colsubstituted.clear();
  colsubstituted.shrink_to_fit();
  implicationmap.clear();
  implicationmap.shrink_to_fit();

  implicationmap.resize(2 * ncols, {-1, 0});
  colsubstituted.resize(ncols);
  vubs.clear();
  vubs.shrink_to_fit();
  vubs.resize(ncols);
  vlbs.clear();
  vlbs.shrink_to_fit();
  vlbs.resize(ncols);
  int oldncols = oldvubs.size();

  for (int i = 0; i != oldncols; ++i) {
    int newi = orig2reducedcol[i];

    if (newi == -1) continue;

    for (const std::pair<int, VarBound>& oldvub : oldvubs[i]) {
      if (orig2reducedcol[oldvub.first] == -1) continue;

      addVUB(newi, orig2reducedcol[oldvub.first], oldvub.second.coef,
             oldvub.second.constant);
    }

    for (const std::pair<int, VarBound>& oldvlb : oldvlbs[i]) {
      if (orig2reducedcol[oldvlb.first] == -1) continue;

      addVLB(newi, orig2reducedcol[oldvlb.first], oldvlb.second.coef,
             oldvlb.second.constant);
    }

    // todo also add old implications once implications can be added
    // incrementally for now we discard the old implications as they might be
    // weaker then newly computed ones and adding them would block computation
    // of new implications
  }
}