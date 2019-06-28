#include "highs_c_api.h"
#include "Highs.h"

int Highs_call(int numcol, int numrow, int numnz, double *colcost,
              double *collower, double *colupper, double *rowlower,
              double *rowupper, int *astart, int *aindex, double *avalue,
              double *colvalue, double *coldual, double *rowvalue,
              double *rowdual, int *colbasisstatus, int *rowbasisstatus) {
  Highs highs;

  int status =
      Highs_loadModel(&highs, numcol, numrow, numnz, colcost, collower,
                      colupper, rowlower, rowupper, astart, aindex, avalue);
  if (status != 1) {
    return status;
  }

  status = (int)highs.run();

  if (status == 1 || status == 17) {
    HighsSolution solution;
    HighsBasis basis;

    solution = highs.getSolution();
    basis = highs.getBasis();

    for (int i = 0; i < numcol; i++) {
      colvalue[i] = solution.col_value[i];
      coldual[i] = solution.col_dual[i];

      colbasisstatus[i] = (int)basis.col_status[i];
    }

    for (int i = 0; i < numrow; i++) {
      rowvalue[i] = solution.row_value[i];
      rowdual[i] = solution.row_dual[i];

      rowbasisstatus[i] = (int)basis.row_status[i];
    }
  }

  return status;
}

void *Highs_create() { return new Highs(); }

void Highs_destroy(void *highs) { delete (Highs *)highs; }

int Highs_run(void *highs) { return (int)((Highs *)highs)->run(); }

int Highs_readFromFile(void *highs, const char *filename) {
  return (int)((Highs *)highs)->initializeFromFile(std::string(filename));
}

int Highs_writeToFile(void *highs, const char *filename) {
  return (int)((Highs *)highs)->writeToFile(std::string(filename));
}

int Highs_loadModel(void *highs, int numcol, int numrow, int numnz,
                    double *colcost, double *collower, double *colupper,
                    double *rowlower, double *rowupper, int *astart,
                    int *aindex, double *avalue) {
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

  return (int)((Highs *)highs)->initializeLp(lp);
}

int Highs_setHighsOptionValue(void *highs, const char *option,
                              const char *value) {
  return (int)((Highs *)highs)
      ->setHighsOptionValue(std::string(option), std::string(value));
}

void Highs_getSolution(void *highs, double *colvalue, double *coldual,
                       double *rowvalue, double *rowdual) {
  HighsSolution solution = ((Highs *)highs)->getSolution();

  for (int i = 0; i < solution.col_value.size(); i++) {
    colvalue[i] = solution.col_value[i];
  }

  for (int i = 0; i < solution.col_dual.size(); i++) {
    coldual[i] = solution.col_dual[i];
  }

  for (int i = 0; i < solution.row_value.size(); i++) {
    rowvalue[i] = solution.row_value[i];
  }

  for (int i = 0; i < solution.row_dual.size(); i++) {
    rowdual[i] = solution.row_dual[i];
  }
}

void Highs_getBasis(void *highs, int *colstatus, int *rowstatus) {
  HighsBasis basis = ((Highs *)highs)->getBasis();

  for (int i = 0; i < basis.col_status.size(); i++) {
    colstatus[i] = (int)basis.col_status[i];
  }

  for (int i = 0; i < basis.row_status.size(); i++) {
    rowstatus[i] = (int)basis.row_status[i];
  }
}

double Highs_getObjectiveValue(void *highs) {
  return ((Highs *)highs)->getObjectiveValue();
}

int Highs_getIterationCount(void *highs) {
  return ((Highs *)highs)->getIterationCount();
}

int Highs_addRow(void *highs, const double lower, const double upper,
                 const int num_new_nz, const int *indices,
                 const double *values) {
  return ((Highs *)highs)->addRow(lower, upper, num_new_nz, indices, values);
}

int Highs_addRows(void *highs, const int num_new_row, const double *lower,
                  const double *upper, const int num_new_nz, const int *starts,
                  const int *indices, const double *values) {
  return ((Highs *)highs)
      ->addRows(num_new_row, lower, upper, num_new_nz, starts, indices, values);
}

int Highs_addCol(void *highs, const double cost, const double lower,
                 const double upper, const int num_new_nz, const int *indices,
                 const double *values) {
  return ((Highs *)highs)
      ->addCol(cost, lower, upper, num_new_nz, indices, values);
}

int Highs_addCols(void *highs, const int num_new_col, const double *costs,
                  const double *lower, const double *upper,
                  const int num_new_nz, const int *starts, const int *indices,
                  const double *values) {
  return ((Highs *)highs)
      ->addCols(num_new_col, costs, lower, upper, num_new_nz, starts, indices,
                values);
}

int Highs_changeObjectiveSense(void *highs, const int sense) {
  return ((Highs *)highs)->changeObjectiveSense(sense);
}

int Highs_changeColCost(void *highs, const int col, const double cost) {
  return ((Highs *)highs)->changeColCost(col, cost);
}

int Highs_changeColsCostBySet(void *highs, const int num_set_entries,
                              const int *set, const double *cost) {
  return ((Highs *)highs)->changeColsCost(num_set_entries, set, cost);
}

int Highs_changeColsCostByMask(void *highs, const int *mask,
                               const double *cost) {
  return ((Highs *)highs)->changeColsCost(mask, cost);
}

int Highs_changeColBounds(void *highs, const int col, const double lower,
                          const double upper) {
  return ((Highs *)highs)->changeColBounds(col, lower, upper);
}

int Highs_changeColsBoundsByRange(void *highs, const int from_col,
                                  const int to_col, const double *lower,
                                  const double *upper) {
  return ((Highs *)highs)->changeColsBounds(from_col, to_col, lower, upper);
}

int Highs_changeColsBoundsBySet(void *highs, const int num_set_entries,
                                const int *set, const double *lower,
                                const double *upper) {
  return ((Highs *)highs)->changeColsBounds(num_set_entries, set, lower, upper);
}

int Highs_changeColsBoundsByMask(void *highs, const int *mask,
                                 const double *lower, const double *upper) {
  return ((Highs *)highs)->changeColsBounds(mask, lower, upper);
}

int Highs_changeRowBounds(void *highs, const int row, const double lower,
                          const double upper) {
  return ((Highs *)highs)->changeRowBounds(row, lower, upper);
}

int Highs_changeRowsBoundsBySet(void *highs, const int num_set_entries,
                                const int *set, const double *lower,
                                const double *upper) {
  return ((Highs *)highs)->changeRowsBounds(num_set_entries, set, lower, upper);
}

int Highs_changeRowsBoundsByMask(void *highs, const int *mask,
                                 const double *lower, const double *upper) {
  return ((Highs *)highs)->changeRowsBounds(mask, lower, upper);
}

int Highs_getColsByRange(void *highs, const int from_col, const int to_col,
                         int *num_col, double *costs, double *lower,
                         double *upper, int *num_nz, int *matrix_start,
                         int *matrix_index, double *matrix_value) {
  int numcol, numnz;
  int status = ((Highs *)highs)
                   ->getCols(from_col, to_col, numcol, costs, lower, upper,
                             numnz, matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

int Highs_getColsBySet(void *highs, const int num_set_entries, const int *set,
                       int *num_col, double *costs, double *lower,
                       double *upper, int *num_nz, int *matrix_start,
                       int *matrix_index, double *matrix_value) {
  int numcol, numnz;
  int status = ((Highs *)highs)
                   ->getCols(num_set_entries, set, numcol, costs, lower, upper,
                             numnz, matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

int Highs_getColsByMask(void *highs, const int *mask, int *num_col,
                        double *costs, double *lower, double *upper,
                        int *num_nz, int *matrix_start, int *matrix_index,
                        double *matrix_value) {
  int numcol, numnz;
  int status = ((Highs *)highs)
                   ->getCols(mask, numcol, costs, lower, upper, numnz,
                             matrix_start, matrix_index, matrix_value);
  *num_col = numcol;
  *num_nz = numnz;
  return status;
}

int Highs_getRowsByRange(void *highs, const int from_row, const int to_row,
                         int *num_row, double *lower, double *upper,
                         int *num_nz, int *matrix_start, int *matrix_index,
                         double *matrix_value) {
  int numrow, numnz;
  int status = ((Highs *)highs)
                   ->getRows(from_row, to_row, numrow, lower, upper, numnz,
                             matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

int Highs_getRowsBySet(void *highs, const int num_set_entries, const int *set,
                       int *num_row, double *lower, double *upper, int *num_nz,
                       int *matrix_start, int *matrix_index,
                       double *matrix_value) {
  int numrow, numnz;
  int status = ((Highs *)highs)
                   ->getRows(num_set_entries, set, numrow, lower, upper, numnz,
                             matrix_start, matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

int Highs_getRowsByMask(void *highs, const int *mask, int *num_row,
                        double *lower, double *upper, int *num_nz,
                        int *matrix_start, int *matrix_index,
                        double *matrix_value) {
  int numrow, numnz;
  int status = ((Highs *)highs)
                   ->getRows(mask, numrow, lower, upper, numnz, matrix_start,
                             matrix_index, matrix_value);
  *num_row = numrow;
  *num_nz = numnz;
  return status;
}

int Highs_deleteColsByRange(void *highs, const int from_col, const int to_col) {
  return ((Highs *)highs)->deleteCols(from_col, to_col);
}

int Highs_deleteColsBySet(void *highs, const int num_set_entries,
                          const int *set) {
  return ((Highs *)highs)->deleteCols(num_set_entries, set);
}

int Highs_deleteColsByMask(void *highs, int *mask) {
  return ((Highs *)highs)->deleteCols(mask);
}

int Highs_deleteRowsByRange(void *highs, const int from_row, const int to_row) {
  return ((Highs *)highs)->deleteRows(from_row, to_row);
}

int Highs_deleteRowsBySet(void *highs, const int num_set_entries,
                          const int *set) {
  return ((Highs *)highs)->deleteRows(num_set_entries, set);
}

int Highs_deleteRowsByMask(void *highs, int *mask) {
  return ((Highs *)highs)->deleteRows(mask);
}

int Highs_getNumCols(void *highs) { return ((Highs *)highs)->getLp().numCol_; }

int Highs_getNumRows(void *highs) { return ((Highs *)highs)->getLp().numRow_; }

int Highs_getNumNz(void *highs) { return ((Highs *)highs)->getLp().nnz_; }
