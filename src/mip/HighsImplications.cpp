#include "mip/HighsImplications.h"

#include "mip/HighsCliqueTable.h"

bool HighsImplications::computeImplications(int col, bool val) {
  globaldomain.propagate();
  if (globaldomain.infeasible()) return true;
  const auto& domchgstack = globaldomain.getDomainChangeStack();

  int changedend = globaldomain.getChangedCols().size();
  int proprowend = globaldomain.getPropagateRows().size();
  int propcutend = globaldomain.getPropagateCuts().size();

  if (val)
    globaldomain.changeBound(HighsBoundType::Lower, col, 1);
  else
    globaldomain.changeBound(HighsBoundType::Upper, col, 0);

  if (globaldomain.infeasible()) {
    globaldomain.backtrack();
    globaldomain.clearChangedCols(changedend);
    globaldomain.clearPropagateRows(proprowend);
    globaldomain.clearPropagateCuts(propcutend);

    cliquetable.vertexInfeasible(globaldomain, col, val);

    return true;
  }

  cliquetable.addImplications(globaldomain, col, val);

  size_t stackimplicstart = domchgstack.size();

  globaldomain.propagate();

  if (globaldomain.infeasible()) {
    globaldomain.backtrack();
    globaldomain.clearChangedCols(changedend);
    globaldomain.clearPropagateRows(proprowend);
    globaldomain.clearPropagateCuts(propcutend);

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
  globaldomain.clearPropagateRows(proprowend);
  globaldomain.clearPropagateCuts(propcutend);

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
    if (globaldomain.infeasible()) return true;
  }

  implications.erase(binstart, implications.end());

  implicationmap[loc].start = implstart;
  implicationmap[loc].num = implications.size() - implstart;

  return false;
}

bool HighsImplications::runProbing(int col, int& numboundchgs) {
  if (globaldomain.isBinary(col) && !implicationsCached(col, 1) &&
      !implicationsCached(col, 0)) {
    const HighsDomainChange* implicsup;
    const HighsDomainChange* implicsdown;
    int nimplicsup;
    int nimplicsdown;
    bool infeasible;
    nimplicsup = getImplications(col, 1, implicsup, infeasible);
    if (globaldomain.infeasible()) return true;

    if (infeasible) return true;

    nimplicsdown = getImplications(col, 0, implicsdown, infeasible);
    if (globaldomain.infeasible()) return true;
    if (infeasible) return true;

    // analyze implications
    int u = 0;
    int d = 0;

    while (u < nimplicsup && d < nimplicsdown) {
      if (implicsup[u] < implicsdown[d])
        ++u;
      else if (implicsdown[d] < implicsup[u])
        ++d;
      else {
        assert(implicsup[u].boundtype == implicsdown[d].boundtype);
        assert(implicsup[u].column == implicsdown[d].column);
        if ((implicsup[u].boundtype == HighsBoundType::Lower &&
             implicsdown[d].boundval < implicsup[u].boundval) ||
            (implicsup[u].boundtype == HighsBoundType::Upper &&
             implicsdown[d].boundval > implicsup[u].boundval))
          globaldomain.changeBound(implicsdown[d], -2);
        else
          globaldomain.changeBound(implicsup[u], -2);
        assert(!globaldomain.infeasible());
        ++numboundchgs;
        globaldomain.propagate();
        assert(!globaldomain.infeasible());
        ++u;
        ++d;
      }
    }
  }

  return false;
}