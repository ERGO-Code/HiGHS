#include "mip/HighsImplications.h"

#include "mip/HighsCliqueTable.h"
#include "mip/HighsMipSolverData.h"

bool HighsImplications::computeImplications(int col, bool val) {
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

    cliquetable.addClique(globaldomain, clique, 2);
    if (globaldomain.infeasible() || globaldomain.isFixed(col)) return true;
  }

  implications.erase(binstart, implications.end());

  implicationmap[loc].start = implstart;
  implicationmap[loc].num = implications.size() - implstart;

  return false;
}

bool HighsImplications::runProbing(int col, int& numboundchgs) {
  if (globaldomain.isBinary(col) && !implicationsCached(col, 1) &&
      !implicationsCached(col, 0) &&
      cliquetable.getSubstitution(col) == nullptr) {
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