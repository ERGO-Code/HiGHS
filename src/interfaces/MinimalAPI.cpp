#include "Highs.h"
#include "highs_lp_solver.h"

extern "C" void callhighs(int numcol, int numrow, int numnz, double* colcost,
                          double* collower, double* colupper, double* rowlower,
                          double* rowupper, int* astart, int* aindex,
                          double* avalue, double* col_value, double* col_dual,
                          double* row_value, double* row_dual,
                          int* col_basis_status, int* row_basis_status) {
  HighsLp lp;
  lp.numCol_ = numcol;
  lp.numRow_ = numrow;
  lp.nnz_ = numnz;

  lp.colCost_.resize(numcol);
  lp.colLower_.resize(numcol);
  lp.colUpper_.resize(numcol);

  lp.rowLower_.resize(numrow);
  lp.rowUpper_.resize(numrow);
  lp.Astart_.resize(numcol + 1);
  lp.Aindex_.resize(numnz);
  lp.Avalue_.resize(numnz);

  lp.colCost_.assign(colcost, colcost + numcol);
  lp.colLower_.assign(collower, collower + numcol);
  lp.colUpper_.assign(colupper, colupper + numcol);

  lp.rowLower_.assign(rowlower, rowlower + numrow);
  lp.rowUpper_.assign(rowupper, rowupper + numcol);
  lp.Astart_.assign(astart, astart + numcol + 1);
  lp.Aindex_.assign(aindex, aindex + numnz);
  lp.Avalue_.assign(avalue, avalue + numnz);

  HighsOptions options;
  HighsStatus status;

  Highs highs(options);

  status = highs.initializeLp(lp);
  status = highs.run();

  HighsSolution solution;
  HighsBasis basis;

  solution = highs.getSolution();
  basis = highs.getBasis();

  for (int i=0; i<numcol; i++) {
    col_value[i] = solution.col_value[i];
    col_dual[i] = solution.col_dual[i];

    col_basis_status[i] = (int)basis.col_status[i];
  }

  for (int i=0; i<numrow; i++) {
    row_value[i] = solution.row_value[i];
    row_dual[i] = solution.row_dual[i];

    row_basis_status[i] = (int)basis.row_status[i];
  }

}