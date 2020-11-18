#include "mip/HighsImplications.h"

#include "mip/HighsCliqueTable.h"

bool HighsImplications::computeImplications(int col, bool val) {
  globaldomain.propagate();
  const auto& domchgstack = globaldomain.getDomainChangeStack();

  int changedend = globaldomain.getChangedCols().size();
  int proprowend = globaldomain.getPropagateRows().size();
  int propcutend = globaldomain.getPropagateCuts().size();

  if (val)
    globaldomain.changeBound(HighsBoundType::Lower, col, 1);
  else
    globaldomain.changeBound(HighsBoundType::Upper, col, 0);

  cliquetable.addImplications(globaldomain, col, val);

  size_t stackimplicstart = domchgstack.size();

  globaldomain.propagate();

  if (globaldomain.infeasible()) {
    globaldomain.backtrack();
    globaldomain.clearChangedCols(changedend);
    globaldomain.clearPropagateRows(proprowend);
    globaldomain.clearPropagateCuts(propcutend);

    if (val)
      globaldomain.changeBound(HighsBoundType::Upper, col, 0, -2);
    else
      globaldomain.changeBound(HighsBoundType::Lower, col, 1, -2);

    globaldomain.propagate();

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

  HighsCliqueTable::CliqueVar clique[2];
  clique[0] = HighsCliqueTable::CliqueVar(col, val);

  for (auto i = binstart; i != implications.end(); ++i) {
    if (i->boundtype == HighsBoundType::Lower)
      clique[1] = HighsCliqueTable::CliqueVar(i->column, 0);
    else
      clique[1] = HighsCliqueTable::CliqueVar(i->column, 1);

    cliquetable.addClique(clique, 2);
  }

  implications.erase(binstart, implications.end());

  implicationmap[loc].start = implstart;
  implicationmap[loc].num = implications.size() - implstart;

  return false;
}