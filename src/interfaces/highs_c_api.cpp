#include "highs_c_api.h"
#include "Highs.h"
#include "io/Filereader.h"

void callhighs(int numcol, int numrow, int numnz, double *colcost,
               double *collower, double *colupper, double *rowlower,
               double *rowupper, int *astart, int *aindex, double *avalue,
               double *colvalue, double *coldual, double *rowvalue,
               double *rowdual, int *colbasisstatus, int *rowbasisstatus) {
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

void *Highs_create() { return new Highs(); }

void Highs_destroy(void *highs) { delete (Highs *)highs; }

int Highs_run(void *highs) { return (int)((Highs *)highs)->run(); }

int Highs_readFromFile(void *highs, const char *filename) {
  return (int)((Highs *)highs)->initializeFromFile(std::string(filename));
}

int Highs_writeToFile(void *highs, const char *filename) {
  return (int)((Highs *)highs)->writeToFile(std::string(filename));
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

int Highs_changeColCost(
    void *highs,       //!< HiGHS object reference
    const int col,     //!< The index of the column whose cost is to change
    const double cost  //!< The new cost
) {
  return ((Highs *)highs)->changeColCost(col, cost);
}

int Highs_changeColsCostBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set,     //!< Array of size num_set_entries with indices of
                        //!< columns whose costs change
    const double *cost  //!< Array of size num_set_entries with new costs
) {
  return ((Highs *)highs)->changeColsCost(num_set_entries, set, cost);
}

int Highs_changeColsCostByMask(
    void *highs,        //!< HiGHS object reference
    const int *mask,    //!< Full length array with 1 => change; 0 => not
    const double *cost  //!< Full length array of new costs
) {
  return ((Highs *)highs)->changeColsCost(mask, cost);
}

int Highs_changeColBounds(
    void *highs,         //!< HiGHS object reference
    const int col,       //!< The index of the column whose bounds are to change
    const double lower,  //!< The new lower bound
    const double upper   //!< The new upper bound
) {
  return ((Highs *)highs)->changeColBounds(col, lower, upper);
}

int Highs_changeColsBoundsByRange(
    void *highs,         //!< HiGHS object reference
    const int from_col,  //!< The index of the first column whose bounds change
    const int to_col,    //!< One more than the index of the last column whose
                         //!< bounds change
    const double
        *lower,  //!< Array of size to_col-from_col with new lower bounds
    const double
        *upper  //!< Array of size to_col-from_col with new upper bounds
) {
  return ((Highs *)highs)->changeColsBounds(from_col, to_col, lower, upper);
}

int Highs_changeColsBoundsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set,  //!< Array of size num_set_entries with indices of
                     //!< columns whose bounds change
    const double
        *lower,  //!< Array of size num_set_entries with new lower bounds
    const double
        *upper  //!< Array of size num_set_entries with new upper bounds
) {
  return ((Highs *)highs)->changeColsBounds(num_set_entries, set, lower, upper);
}

int Highs_changeColsBoundsByMask(
    void *highs,          //!< HiGHS object reference
    const int *mask,      //!< Full length array with 1 => change; 0 => not
    const double *lower,  //!< Full length array of new lower bounds
    const double *upper   //!< Full length array of new upper bounds
) {
  return ((Highs *)highs)->changeColsBounds(mask, lower, upper);
}

int Highs_changeRowBounds(
    void *highs,         //!< HiGHS object reference
    const int row,       //!< The index of the row whose bounds are to change
    const double lower,  //!< The new lower bound
    const double upper   //!< The new upper bound
) {
  return ((Highs *)highs)->changeRowBounds(row, lower, upper);
}

int Highs_changeRowsBoundsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set,  //!< Array of size num_set_entries with indices of rows
                     //!< whose bounds change
    const double
        *lower,  //!< Array of size num_set_entries with new lower bounds
    const double
        *upper  //!< Array of size num_set_entries with new upper bounds
) {
  return ((Highs *)highs)->changeRowsBounds(num_set_entries, set, lower, upper);
}

int Highs_changeRowsBoundsByMask(
    void *highs,          //!< HiGHS object reference
    const int *mask,      //!< Full length array with 1 => change; 0 => not
    const double *lower,  //!< Full length array of new lower bounds
    const double *upper   //!< Full length array of new upper bounds
) {
  return ((Highs *)highs)->changeRowsBounds(mask, lower, upper);
}

int Highs_getColsByRange(
    void *highs,          //!< HiGHS object reference
    const int from_col,   //!< The index of the first column to
                          //!< get from the model
    const int to_col,     //!< One more than the last column to get
                          //!< from the model
    int num_col,          //!< Number of columns got from the model
    double *costs,        //!< Array of size num_col with costs
    double *lower,        //!< Array of size num_col with lower bounds
    double *upper,        //!< Array of size num_col with upper bounds
    int num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!< Array of size num_col with start
                          //!< indices of the columns
    int *matrix_index,    //!< Array of size num_nz with row
                          //!< indices for the columns
    double *matrix_value  //!< Array of size num_nz with row
                          //!< values for the columns
) {
  return ((Highs *)highs)
      ->getCols(from_col, to_col, num_col, costs, lower, upper, num_nz,
                matrix_start, matrix_index, matrix_value);
}

int Highs_getColsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set,             //!< Array of size num_set_entries with indices
                                //!< of columns to get
    int num_col,                //!< Number of columns got from the model
    double *costs,              //!< Array of size num_col with costs
    double *lower,              //!< Array of size num_col with lower bounds
    double *upper,              //!< Array of size num_col with upper bounds
    int num_nz,                 //!< Number of nonzeros got from the model
    int *matrix_start,          //!< Array of size num_col with start indices
                                //!< of the columns
    int *matrix_index,          //!< Array of size num_nz with row indices
                                //!< for the columns
    double *matrix_value        //!< Array of size num_nz with row values
                                //!< for the columns
) {
  return ((Highs *)highs)
      ->getCols(num_set_entries, set, num_col, costs, lower, upper, num_nz,
                matrix_start, matrix_index, matrix_value);
}

int Highs_getColsByMask(
    void *highs,          //!< HiGHS object reference
    const int *mask,      //!< Full length array with 1 => get; 0 => not
    int num_col,          //!< Number of columns got from the model
    double *costs,        //!< Array of size num_col with costs
    double *lower,        //!< Array of size num_col with lower bounds
    double *upper,        //!< Array of size num_col with upper bounds
    int num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!<  Array of size num_col with start
                          //!<  indices of the columns
    int *matrix_index,    //!<  Array of size num_nz with row indices
                          //!<  for the columns
    double *matrix_value  //!<  Array of size num_nz with row values
                          //!<  for the columns
) {
  return ((Highs *)highs)
      ->getCols(mask, num_col, costs, lower, upper, num_nz, matrix_start,
                matrix_index, matrix_value);
}

int Highs_getRowsByRange(
    void *highs,          //!< HiGHS object reference
    const int from_row,   //!< The index of the first row to get from the model
    const int to_row,     //!< One more than the last row get from the model
    int num_row,          //!< Number of rows got from the model
    double *lower,        //!< Array of size num_row with lower bounds
    double *upper,        //!< Array of size num_row with upper bounds
    int num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!< Array of size num_row with start indices of the
                          //!< rows
    int *matrix_index,    //!< Array of size num_nz with column indices for the
                          //!< rows
    double *matrix_value  //!< Array of size num_nz with column values for the
                          //!< rows
) {
  return ((Highs *)highs)
      ->getRows(from_row, to_row, num_row, lower, upper, num_nz, matrix_start,
                matrix_index, matrix_value);
}

int Highs_getRowsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set,             //!< Array of size num_set_entries with indices
                                //!< of rows to get
    int num_row,                //!< Number of rows got from the model
    double *lower,              //!< Array of size num_row with lower bounds
    double *upper,              //!< Array of size num_row with upper bounds
    int num_nz,                 //!< Number of nonzeros got from the model
    int *matrix_start,          //!< Array of size num_row with start indices
                                //!< of the rows
    int *matrix_index,          //!< Array of size num_nz with column indices
                                //!< for the rows
    double *matrix_value        //!< Array of size num_nz with column
                                //!< values for the rows
) {
  return ((Highs *)highs)
      ->getRows(num_set_entries, set, num_row, lower, upper, num_nz,
                matrix_start, matrix_index, matrix_value);
}

int Highs_getRowsByMask(
    void *highs,          //!< HiGHS object reference
    const int *mask,      //!< Full length array with 1 => get; 0 => not
    int num_row,          //!< Number of rows got from the model
    double *lower,        //!< Array of size num_row with lower bounds
    double *upper,        //!< Array of size num_row with upper bounds
    int num_nz,           //!< Number of nonzeros got from the model
    int *matrix_start,    //!< Array of size num_row with start indices
                          //!< of the rows
    int *matrix_index,    //!< Array of size num_nz with column indices
                          //!< for the rows
    double *matrix_value  //!< Array of size num_nz with column
                          //!< values for the rows
) {
  return ((Highs *)highs)
      ->getRows(mask, num_row, lower, upper, num_nz, matrix_start, matrix_index,
                matrix_value);
}

int Highs_deleteColsByRange(
    void *highs,         //!< HiGHS object reference
    const int from_col,  //!< The index of the first column
                         //!< to delete from the model
    const int to_col     //!< One more than the last column to
                         //!< delete from the model
) {
  return ((Highs *)highs)->deleteCols(from_col, to_col);
}

int Highs_deleteColsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set  //!< Array of size num_set_entries with indices of columns
                    //!< to delete
) {
  return ((Highs *)highs)->deleteCols(num_set_entries, set);
}

int Highs_deleteColsByMask(
    void *highs,  //!< HiGHS object reference
    int *mask     //!< Full length array with 1 => delete; 0 => not
) {
  return ((Highs *)highs)->deleteCols(mask);
}

int Highs_deleteRowsByRange(
    void *highs,  //!< HiGHS object reference
    const int
        from_row,     //!< The index of the first row to delete from the model
    const int to_row  //!< One more than the last row delete from the model
) {
  return ((Highs *)highs)->deleteRows(from_row, to_row);
}

int Highs_deleteRowsBySet(
    void *highs,                //!< HiGHS object reference
    const int num_set_entries,  //!< The number of indides in the set
    const int *set  //!< Array of size num_set_entries with indices of columns
                    //!< to delete
) {
  return ((Highs *)highs)->deleteRows(num_set_entries, set);
}

int Highs_deleteRowsByMask(
    void *highs,  //!< HiGHS object reference
    int *mask     //!< Full length array with 1 => delete; 0 => not
) {
  return ((Highs *)highs)->deleteRows(mask);
}