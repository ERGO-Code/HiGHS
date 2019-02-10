/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsSimplexInterface.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HighsSimplexInterface.h"

void HighsSimplexInterface::load_from_arrays(
					     int XnumCol,
					     int Xsense,
					     const double* XcolCost,
					     const double* XcolLower,
					     const double* XcolUpper,
					     int XnumRow,
					     const double* XrowLower,
					     const double* XrowUpper,
					     int XnumNz,
					     const int* XAstart,
					     const int* XAindex,
					     const double* XAvalue
					      ) {
  HighsLp &lp = highs_model_object.lp_;
  //  printf("load_fromArrays: XnumCol = %d; XnumRow = %d; XnumNz = %d\n",
  //  XnumCol, XnumRow, XnumNz);
  assert(XnumCol > 0);
  assert(XnumRow > 0);

  lp.numCol_ = XnumCol;
  lp.numRow_ = XnumRow;
  lp.sense_ = Xsense;
  int numNz = XnumNz;
  lp.colCost_.assign(&XcolCost[0], &XcolCost[0] + lp.numCol_);
  lp.colLower_.assign(&XcolLower[0], &XcolLower[0] + lp.numCol_);
  lp.colUpper_.assign(&XcolUpper[0], &XcolUpper[0] + lp.numCol_);
  lp.rowLower_.assign(&XrowLower[0], &XrowLower[0] + lp.numRow_);
  lp.rowUpper_.assign(&XrowUpper[0], &XrowUpper[0] + lp.numRow_);
  lp.Astart_.assign(&XAstart[0], &XAstart[0] + lp.numCol_ + 1);
  lp.Aindex_.assign(&XAindex[0], &XAindex[0] + numNz);
  lp.Avalue_.assign(&XAvalue[0], &XAvalue[0] + numNz);
  // Assign and initialise the scaling factors

  // Initialise with a logical basis then allocate and populate (where
  // possible) work* arrays and allocate basis* arrays
  
  }
  
void HighsSimplexInterface::get_primal_dual_values(vector<double> &XcolValue,
						   vector<double> &XcolDual,
						   vector<double> &XrowValue,
						   vector<double> &XrowDual
						   ) {
  HighsScale &scale = highs_model_object.scale_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  // Take primal solution
  vector<double> value = simplex_info.workValue_;
  for (int iRow = 0; iRow < solver_lp.numRow_; iRow++)
    value[basis.basicIndex_[iRow]] = simplex_info.baseValue_[iRow];
  // Take dual solution
  vector<double> dual = simplex_info.workDual_;
  for (int iRow = 0; iRow < solver_lp.numRow_; iRow++) dual[basis.basicIndex_[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < solver_lp.numCol_; iCol++) {
    value[iCol] *= scale.col_[iCol];
    dual[iCol] /= (scale.col_[iCol] / scale.cost_);
  }
  for (int iRow = 0, iTot = solver_lp.numCol_; iRow < solver_lp.numRow_; iRow++, iTot++) {
    value[iTot] /= scale.row_[iRow];
    dual[iTot] *= (scale.row_[iRow] * scale.cost_);
  }

  //************** part 2: gepr and gedu
  // Now we can get the solution
  XcolValue.resize(solver_lp.numCol_);
  XcolDual.resize(solver_lp.numCol_);
  XrowValue.resize(solver_lp.numRow_);
  XrowDual.resize(solver_lp.numRow_);

  double *valuePtr = &value[0];
  for (int i = 0; i < solver_lp.numRow_; i++) XrowValue[i] = -valuePtr[i + solver_lp.numCol_];
  for (int i = 0; i < solver_lp.numCol_; i++) XcolValue[i] = valuePtr[i];
  for (int i = 0; i < solver_lp.numRow_; i++) XrowDual[i] = solver_lp.sense_ * dual[i + solver_lp.numCol_];
  for (int i = 0; i < solver_lp.numCol_; i++) XcolDual[i] = solver_lp.sense_ * dual[i];
}

double HighsSimplexInterface::get_lp_objective_value(vector<double> &XcolValue) {
  HighsLp &lp = highs_model_object.lp_;

  double lp_objective_value = 0;
  for (int i = 0; i < lp.numCol_; i++) lp_objective_value += XcolValue[i] * lp.colCost_[i];
  return lp_objective_value;
}

#ifdef HiGHSDEV
void HighsSimplexInterface::check_load_from_arrays() {
  HighsLp &lp = highs_model_object.lp_;
  // Use the arrays read from an MPS file to test the routine to
  // read a model passed by arrays. First copy the data.
  int XnumCol = lp.numCol_;
  int XnumRow = lp.numRow_;
  int XnumNz = lp.Astart_[lp.numCol_];
  int Xsense = lp.sense_;
  vector<double> XcolCost;
  vector<double> colLower;
  vector<double> XcolUpper;
  vector<double> XrowLower;
  vector<double> XrowUpper;
  vector<int> XAstart;
  vector<int> XAindex;
  vector<double> XAvalue;

  XcolCost.assign(&lp.colCost_[0], &lp.colCost_[0] + XnumCol);
  colLower.assign(&lp.colLower_[0], &lp.colLower_[0] + XnumCol);
  XcolUpper.assign(&lp.colUpper_[0], &lp.colUpper_[0] + XnumCol);
  XrowLower.assign(&lp.rowLower_[0], &lp.rowLower_[0] + XnumRow);
  XrowUpper.assign(&lp.rowUpper_[0], &lp.rowUpper_[0] + XnumRow);
  XAstart.assign(&lp.Astart_[0], &lp.Astart_[0] + XnumCol + 1);
  XAindex.assign(&lp.Aindex_[0], &lp.Aindex_[0] + XnumNz);
  XAvalue.assign(&lp.Avalue_[0], &lp.Avalue_[0] + XnumNz);

  //  clear_solver_lp(highs_model_object);
  load_from_arrays(XnumCol, Xsense, &XcolCost[0], &colLower[0],
		   &XcolUpper[0], XnumRow, &XrowLower[0], &XrowUpper[0], XnumNz,
		   &XAstart[0], &XAindex[0], &XAvalue[0]);
}

void HighsSimplexInterface::check_load_from_postsolve() {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  //  HSimplex simplex_method_;
  bool ok;

  ok = true;printf("Need to call nonbasic_flag_basic_index_ok\n"); //nonbasicFlagBasicIndex_OK(solver_lp.numCol_, solver_lp.numRow_);
  printf("check_load_from_postsolve: return from nonbasicFlagBasicIndex_OK = %d\n", ok);
  assert(ok);

  ok = true;printf("Need to call all_nonbasic_move_vs_work_arrays_ok\n"); //simplex_method_.all_nonbasic_move_vs_work_arrays_ok(highs_model_object);
  printf("check_load_from_postsolve: return from allNonbasicMoveVsWorkArrays_OK = %d\n", ok);
  assert(ok);
}
#endif

void HighsSimplexInterface::util_add_cols(
					  int ncols,
					  const double *XcolCost,
					  const double *colLower,
					  const double *XcolUpper,
					  int nnonz,
					  const int *XAstart,
					  const int *XAindex,
					  const double *XAvalue) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(ncols >= 0);
  assert(nnonz >= 0);
  // ToDo How to check that solver_lp.Astart_[solver_lp.numCol_] exists in util_addCols?
#ifdef HiGHSDEV
  printf("Called model.util_addCols(ncols=%d, nnonz = %d)\n", ncols, nnonz);
#endif

  if (ncols == 0) return;

  int nwNumCol = solver_lp.numCol_ + ncols;
  solver_lp.colCost_.resize(nwNumCol);
  solver_lp.colLower_.resize(nwNumCol);
  solver_lp.colUpper_.resize(nwNumCol);
  scale.col_.resize(nwNumCol);
  solver_lp.Astart_.resize(nwNumCol + 1);

  // Note that the new columns must have starts, even if they have no entries
  // (yet)
  for (int col = 0; col < ncols; col++) {
    solver_lp.colCost_[solver_lp.numCol_ + col] = XcolCost[col];
    solver_lp.colLower_[solver_lp.numCol_ + col] = colLower[col];
    solver_lp.colUpper_[solver_lp.numCol_ + col] = XcolUpper[col];
    scale.col_[solver_lp.numCol_ + col] = 1.0;
    //    printf("In HModel::util_addCols: column %d: setting
    //    solver_lp.Astart_[solver_lp.numCol_+col+1] = %d \n", col, solver_lp.Astart_[solver_lp.numCol_]);
    solver_lp.Astart_[solver_lp.numCol_ + col + 1] = solver_lp.Astart_[solver_lp.numCol_];
  }

  //  printf("In HModel::util_addCols: nnonz = %d; cuNnonz = %d\n", nnonz,
  //  solver_lp.Astart_[solver_lp.numCol_]); 
  if (nnonz > 0) {
    // Determine the current number of nonzeros
    int cuNnonz = solver_lp.Astart_[solver_lp.numCol_];

    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    // solver_lp.Astart_.resize(nwNumCol+1);
    solver_lp.Aindex_.resize(nwNnonz);
    solver_lp.Avalue_.resize(nwNnonz);

    // Add the new columns
    for (int col = 0; col < ncols; col++) {
      //      printf("In HModel::util_addCols: column %d: setting
      //      solver_lp.Astart_[solver_lp.numCol_+col] = %d = %d + %d\n",
      //             col, XAstart[col] + cuNnonz, XAstart[col], cuNnonz); 
      solver_lp.Astart_[solver_lp.numCol_ + col] = XAstart[col] + cuNnonz;
    }
    //    printf("In HModel::util_addCols: setting solver_lp.Astart_[solver_lp.numCol_+ncols] = %d\n",
    //    nwNnonz);
    solver_lp.Astart_[solver_lp.numCol_ + ncols] = nwNnonz;

    for (int el = 0; el < nnonz; el++) {
      int row = XAindex[el];
      assert(row >= 0);
      assert(row < solver_lp.numRow_);
      solver_lp.Aindex_[cuNnonz + el] = row;
      solver_lp.Avalue_[cuNnonz + el] = XAvalue[el];
    }
  }
  // Increase the number of columns and total number of variables in the model
  solver_lp.numCol_ += ncols;
  //  numTot += ncols;

  //  printf("In HModel::util_addCols: Model now has solver_lp.Astart_[%d] = %d
  //  nonzeros\n", solver_lp.numCol_, solver_lp.Astart_[solver_lp.numCol_]);

  // Update the basis and work vectors correponding to new nonbasic columns
  printf("Call extend_with_logical_basis(highs_model_object, solver_lp.numCol_ - ncols, solver_lp.numCol_ - 1, solver_lp.numRow_, -1);\n");
}

void HighsSimplexInterface::util_delete_cols(
					     int firstcol,
					     int lastcol
					     ) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  assert(firstcol >= 0);
  assert(lastcol < solver_lp.numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_deleteCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
  cout << flush;
#endif
  // Trivial cases are
  //
  // colStep = 0, in which case no columns are removed
  //
  // lastcol = solver_lp.numCol_-1, in which case no columns need be
  // shifted. However, this implies solver_lp.numCol_-colStep=firstcol, in which
  // case the loop is vacuous
  int colStep = lastcol - firstcol + 1;
  if (colStep) {
    for (int col = firstcol; col < solver_lp.numCol_ - colStep; col++) {
      solver_lp.colCost_[col] = solver_lp.colCost_[col + colStep];
      solver_lp.colLower_[col] = solver_lp.colLower_[col + colStep];
      solver_lp.colUpper_[col] = solver_lp.colUpper_[col + colStep];
      scale.col_[col] = scale.col_[col + colStep];
    }
  }
  // Trivial cases are
  //
  // colstep = 0, in which case no columns are removed so elStep = 0
  //
  // lastcol = solver_lp.numCol_-1, in which case no columns need be
  // shifted and the loops are vacuous
  if (colStep) {
    int elOs = solver_lp.Astart_[firstcol];
    int elStep = solver_lp.Astart_[lastcol + 1] - elOs;
    //    printf("El loop over cols %2d [%2d] to %2d [%2d]\n", lastcol+1,
    //    solver_lp.Astart_[lastcol+1], solver_lp.numCol_+1, solver_lp.Astart_[solver_lp.numCol_]-1);
    for (int el = solver_lp.Astart_[lastcol + 1]; el < solver_lp.Astart_[solver_lp.numCol_]; el++) {
      //        printf("Over-write entry %3d [%3d] by entry %3d [%3d]\n",
      //        el-elStep, solver_lp.Aindex_[el-elStep], el, solver_lp.Aindex_[el]);
      solver_lp.Aindex_[el - elStep] = solver_lp.Aindex_[el];
      solver_lp.Avalue_[el - elStep] = solver_lp.Avalue_[el];
    }
    for (int col = firstcol; col <= solver_lp.numCol_ - colStep; col++) {
      //    printf("Over-write start %3d [%3d] by entry %3d [%3d]\n", col,
      //    solver_lp.Astart_[col], col+colStep,  solver_lp.Astart_[col+colStep]-elStep);
      solver_lp.Astart_[col] = solver_lp.Astart_[col + colStep] - elStep;
    }
  }

  // Reduce the number of columns and total number of variables in the model
  solver_lp.numCol_ -= colStep;
  //  numTot -= colStep;

  // ToDo Determine consequences for basis when deleting columns
  // Invalidate matrix copies
  simplex_info.solver_lp_has_matrix_col_wise = false;
  simplex_info.solver_lp_has_matrix_row_wise = false;
}

void HighsSimplexInterface::util_delete_col_set(
						vector<int>& dstat
						) {
  printf("util_delete_col_set not implemented");
  assert(1 == 0);
}


void HighsSimplexInterface::util_extract_cols(
					      int firstcol,
					      int lastcol,
					      double* XcolLower,
					      double* XcolUpper,
					      int* nnonz,
					      int* XAstart,
					      int* XAindex,
					      double* XAvalue
					      ) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(firstcol >= 0);
  assert(lastcol < solver_lp.numCol_);
  assert(firstcol <= lastcol);
#ifdef HiGHSDEV
  printf("Called model.util_extractCols(firstcol=%d, lastcol=%d)\n", firstcol,
         lastcol);
  cout << flush;
#endif
  // Determine the number of columns to be extracted
  // int numExtractCols = lastcol-firstcol+1;
  // printf("Extracting %d columns\n", numExtractCols);cout << flush;
  int elOs = solver_lp.Astart_[firstcol];
  for (int col = firstcol; col <= lastcol; col++) {
    //    printf("Extracting column %d\n", col);cout << flush;
    XcolLower[col - firstcol] = solver_lp.colLower_[col];
    XcolUpper[col - firstcol] = solver_lp.colUpper_[col];
    XAstart[col - firstcol] = solver_lp.Astart_[col] - elOs;
  }
  for (int el = solver_lp.Astart_[firstcol]; el < solver_lp.Astart_[lastcol + 1]; el++) {
    XAindex[el - elOs] = solver_lp.Aindex_[el];
    XAvalue[el - elOs] = solver_lp.Avalue_[el];
  }
  *nnonz = solver_lp.Astart_[lastcol + 1] - elOs;
}


void HighsSimplexInterface::util_add_rows(
					  int nrows,
					  const double *XrowLower,
					  const double *XrowUpper,
					  int nnonz,
					  const int *XARstart,
					  const int *XARindex,
					  const double *XARvalue
					  ) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(nrows >= 0);
  assert(nnonz >= 0);
  assert(nnonz == 0 || solver_lp.numCol_ > 0);
#ifdef HiGHSDEV
  printf("Called model.util_addRows(nrows=%d, nnonz = %d)\n", nrows, nnonz);
#endif

  if (nrows == 0) return;

  int nwNumRow = solver_lp.numRow_ + nrows;
  solver_lp.rowLower_.resize(nwNumRow);
  solver_lp.rowUpper_.resize(nwNumRow);
  scale.row_.resize(nwNumRow);

  for (int row = 0; row < nrows; row++) {
    solver_lp.rowLower_[solver_lp.numRow_ + row] = XrowLower[row];
    solver_lp.rowUpper_[solver_lp.numRow_ + row] = XrowUpper[row];
    scale.row_[solver_lp.numRow_ + row] = 1.0;
  }
  // NB SCIP doesn't have XARstart[nrows] defined, so have to use nnonz for last
  // entry
  if (nnonz > 0) {
    int cuNnonz = solver_lp.Astart_[solver_lp.numCol_];
    vector<int> Alength;
    Alength.assign(solver_lp.numCol_, 0);
    for (int el = 0; el < nnonz; el++) {
      int col = XARindex[el];
      //      printf("El %2d: adding entry in column %2d\n", el, col); 
      assert(col >= 0);
      assert(col < solver_lp.numCol_);
      Alength[col]++;
    }
    // Determine the new number of nonzeros and resize the column-wise matrix
    // arrays
    int nwNnonz = cuNnonz + nnonz;
    solver_lp.Aindex_.resize(nwNnonz);
    solver_lp.Avalue_.resize(nwNnonz);

    // Add the new rows
    // Shift the existing columns to make space for the new entries
    int nwEl = nwNnonz;
    for (int col = solver_lp.numCol_ - 1; col >= 0; col--) {
      // printf("Column %2d has additional length %2d\n", col, Alength[col]);
      int Astart_Colp1 = nwEl;
      nwEl -= Alength[col];
      // printf("Shift: nwEl = %2d\n", nwEl);
      for (int el = solver_lp.Astart_[col + 1] - 1; el >= solver_lp.Astart_[col]; el--) {
        nwEl--;
        // printf("Shift: Over-writing solver_lp.Aindex_[%2d] with solver_lp.Aindex_[%2d]=%2d\n",
        // nwEl, el, solver_lp.Aindex_[el]);
        solver_lp.Aindex_[nwEl] = solver_lp.Aindex_[el];
        solver_lp.Avalue_[nwEl] = solver_lp.Avalue_[el];
      }
      solver_lp.Astart_[col + 1] = Astart_Colp1;
    }
    // printf("After shift: nwEl = %2d\n", nwEl);
    assert(nwEl == 0);
    // util_reportColMtx(solver_lp.numCol_, solver_lp.Astart_, solver_lp.Aindex_, solver_lp.Avalue_);

    // Insert the new entries
    for (int row = 0; row < nrows; row++) {
      int fEl = XARstart[row];
      int lEl = (row < nrows - 1 ? XARstart[row + 1] : nnonz) - 1;
      for (int el = fEl; el <= lEl; el++) {
        int col = XARindex[el];
        nwEl = solver_lp.Astart_[col + 1] - Alength[col];
        Alength[col]--;
        // printf("Insert: row = %2d; col = %2d; solver_lp.Astart_[col+1]-Alength[col] =
        // %2d; Alength[col] = %2d; nwEl = %2d\n", row, col,
        // solver_lp.Astart_[col+1]-Alength[col], Alength[col], nwEl);
        assert(nwEl >= 0);
        assert(el >= 0);
        // printf("Insert: Over-writing solver_lp.Aindex_[%2d] with solver_lp.Aindex_[%2d]=%2d\n",
        // nwEl, el, solver_lp.Aindex_[el]);
        solver_lp.Aindex_[nwEl] = solver_lp.numRow_ + row;
        solver_lp.Avalue_[nwEl] = XARvalue[el];
      }
    }
  }
  // Increase the number of rows and total number of variables in the model
  solver_lp.numRow_ += nrows;
  //  numTot += nrows;

  // Update the basis and work vectors correponding to new basic rows
  printf("Call extendWithLogicalBasis(solver_lp.numCol_, -1, solver_lp.numRow_ - nrows, solver_lp.numRow_ - 1);\n");
  
}

void HighsSimplexInterface::util_delete_rows(
					     int firstrow,
					     int lastrow
					     ) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(firstrow >= 0);
  assert(lastrow < solver_lp.numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
#endif
  // Trivial cases are
  //
  // rowStep = 0, in which case no rows are removed
  //
  // lastrow = solver_lp.numRow_-1, in which case no rows need be
  // shifted. However, this implies solver_lp.numRow_-rowStep=firstrow, in which
  // case the loop is vacuous. However, they still have to be removed
  // from the matrix unless all rows are to be removed
  int rowStep = lastrow - firstrow + 1;
  bool allRows = rowStep == solver_lp.numRow_;
#ifdef HiGHSDEV
  if (allRows) printf("In model.util_deleteRows, aa rows are being removed)\n");
#endif
  if (rowStep) {
    // Was: for (int row = firstrow; row < lastrow; row++) - surely wrong!
    for (int row = firstrow; row < solver_lp.numRow_ - rowStep; row++) {
      solver_lp.rowLower_[row] = solver_lp.rowLower_[row + rowStep];
      solver_lp.rowUpper_[row] = solver_lp.rowUpper_[row + rowStep];
      //    scale.row_[row] = scale.row_[row + rowStep];
    }
    if (!allRows) {
      int nnz = 0;
      for (int col = 0; col < solver_lp.numCol_; col++) {
        int fmEl = solver_lp.Astart_[col];
        solver_lp.Astart_[col] = nnz;
        for (int el = fmEl; el < solver_lp.Astart_[col + 1]; el++) {
          int row = solver_lp.Aindex_[el];
          if (row < firstrow || row > lastrow) {
            if (row < firstrow) {
              solver_lp.Aindex_[nnz] = row;
            } else {
              solver_lp.Aindex_[nnz] = row - rowStep;
            }
            solver_lp.Avalue_[nnz] = solver_lp.Avalue_[el];
            nnz++;
          }
        }
      }
      solver_lp.Astart_[solver_lp.numCol_] = nnz;
    }
  }

  // Reduce the number of rows and total number of variables in the model
  solver_lp.numRow_ -= rowStep;
  //  numTot -= rowStep;

  // Determine consequences for basis when deleting rows
  //  update_solver_lp_status_flags(highs_model, LpAction::DEL_ROWS);
}

void HighsSimplexInterface::util_delete_row_set(
						vector<int>& dstat
						) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  bool rp = false;
  if (rp) {
    printf("Called model.util_deleteRowSet\n");
    printf("Before\n");
  }
  //  solver_lp.reportLp();

  int newRow = 0;
  // Look through the rows removing any being deleted and shifting data
  // for the rest
  for (int row = 0; row < solver_lp.numRow_; row++) {
    if (!dstat[row]) {
      // Row is not deleted
      int var = solver_lp.numCol_ + row;
      int newVar = solver_lp.numCol_ + newRow;
      dstat[row] = newRow;
      solver_lp.rowLower_[newRow] = solver_lp.rowLower_[row];
      solver_lp.rowUpper_[newRow] = solver_lp.rowUpper_[row];
      //    scale.row_[row] = scale.row_[rowStep+row];
      basis.nonbasicFlag_[newVar] = basis.nonbasicFlag_[var];
      basis.nonbasicMove_[newVar] = basis.nonbasicMove_[var];
      simplex_info.workCost_[newVar] = simplex_info.workCost_[var];
      simplex_info.workShift_[newVar] = simplex_info.workShift_[var];
      simplex_info.workLower_[newVar] = simplex_info.workLower_[var];
      simplex_info.workUpper_[newVar] = simplex_info.workUpper_[var];
      simplex_info.workRange_[newVar] = simplex_info.workRange_[var];
      simplex_info.workValue_[newVar] = simplex_info.workValue_[var];
      if (rp)
        printf(
            "   Row %4d: dstat = %2d: Variable %2d becomes %2d; [%11g, %11g]; "
            "nonbasicFlag = %2d; nonbasicMove = %2d\n",
            row, dstat[row], var, newVar, solver_lp.rowLower_[newRow], solver_lp.rowUpper_[newRow],
            basis.nonbasicFlag_[newVar], basis.nonbasicMove_[newVar]);
      newRow++;
    } else {
      // Row is deleted
      dstat[row] = -1;
      if (rp)
        printf("   Row %4d: dstat = %2d: Variable %2d is deleted\n", row,
               dstat[row], solver_lp.numCol_ + row);
    }
  }

  if (rp) {
    printf("After\n");
    for (int row = 0; row < solver_lp.numRow_; row++)
      printf("   Row %4d: dstat = %2d\n", row, dstat[row]);
  }
  // Look through the column-wise matrix, removing entries
  // corresponding to deleted rows and shifting indices for the rest
  int nnz = 0;
  for (int col = 0; col < solver_lp.numCol_; col++) {
    int fmEl = solver_lp.Astart_[col];
    solver_lp.Astart_[col] = nnz;
    for (int el = fmEl; el < solver_lp.Astart_[col + 1]; el++) {
      int row = solver_lp.Aindex_[el];
      if (dstat[row] >= 0) {
        solver_lp.Aindex_[nnz] = dstat[row];
        solver_lp.Avalue_[nnz] = solver_lp.Avalue_[el];
        nnz++;
      }
    }
  }
  solver_lp.Astart_[solver_lp.numCol_] = nnz;

  // Reduce the number of rows and total number of variables in the model
  int dlNumRow = solver_lp.numRow_ - newRow;
#ifdef SCIP_DEV
  if (rp)
    printf("Had %d rows; removed %d rows; now %d rows\n", solver_lp.numRow_, dlNumRow,
           newRow);
#endif
  solver_lp.numRow_ -= dlNumRow;
  //  numTot -= dlNumRow;

  // Count the remaining basic variables: if there are as many as
  // there are (now) rows then the basis is OK. If there are more then some
  // columns have to be made nonbasic - but which?
  int numBasic = 0;
  bool basisOK = true;
  const int numTot = solver_lp.numCol_ + solver_lp.numRow_;
  for (int var = 0; var < numTot; var++) {
    if (!basis.nonbasicFlag_[var]) {
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      if (numBasic > newRow) {
        basisOK = false;
        break;
      }
    }
  }

  if (rp) {
    printf("Now have %d cols; %d rows and %d total\n", solver_lp.numCol_, solver_lp.numRow_, numTot);
    for (int row = 0; row < solver_lp.numRow_; row++)
      printf("Basic variable in row %2d is %2d\n", row, basis.basicIndex_[row]);
    for (int col = 0; col < solver_lp.numCol_; col++)
      printf("Col %2d has nonbasicFlag = %2d\n", col, basis.nonbasicFlag_[col]);
    for (int row = 0; row < solver_lp.numRow_; row++)
      printf("Row %2d (Variable %2d) has nonbasicFlag = %2d\n", row,
             solver_lp.numCol_ + row, basis.nonbasicFlag_[solver_lp.numCol_ + row]);
  }

  if (basisOK) {
    // All rows removed had basic slacks so basis should be OK
#ifdef SCIP_DEV
    // Check that basis is valid basis.
    basisOK = nonbasicFlagBasicIndex_OK(solver_lp.numCol_, solver_lp.numRow_);
    assert(basisOK);
    //    printf("util_deleteRowset: all rows removed are basic slacks so
    //    basisOK\n");
#endif
    // Determine consequences for basis when deleting rows to leave an OK basis
    //  update_solver_lp_status_flags(highs_model, LpAction::DEL_ROWS_BASIS_OK);
  } else {
    assert(basisOK);
#ifdef SCIP_DEV
    printf("util_deleteRowset: not all rows removed are basic slacks\n");
#endif
    // Determine consequences for basis when deleting rows to leave no basis
  //  update_solver_lp_status_flags(highs_model, LpAction::DEL_ROWS);
  }
}


void HighsSimplexInterface::util_extract_rows(
					      int firstrow,
					      int lastrow,
					      double* XrowLower,
					      double* XrowUpper,
					      int* nnonz,
					      int* XARstart,
					      int* XARindex,
					      double* XARvalue
					      ) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  HighsScale &scale = highs_model_object.scale_;
  assert(firstrow >= 0);
  assert(lastrow < solver_lp.numRow_);
  assert(firstrow <= lastrow);
#ifdef HiGHSDEV
  printf("Called model.util_extractRows(firstrow=%d, lastrow=%d)\n", firstrow,
         lastrow);
#endif
  // Determine the number of rows to be extracted
  int numExtractRows = lastrow - firstrow + 1;
  //    printf("Extracting %d rows\n", numExtractRows);
  for (int row = firstrow; row <= lastrow; row++) {
    // printf("Extracting row %d\n", row);
    XrowLower[row - firstrow] = solver_lp.rowLower_[row];
    XrowUpper[row - firstrow] = solver_lp.rowUpper_[row];
    // printf("Extracted row %d from %d with bounds [%g, %g]\n",
    //	   row-firstrow, row, XrowLower[row-firstrow],
    // XrowUpper[row-firstrow]);
  }
  // Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = solver_lp.Astart_[0]; el < solver_lp.Astart_[solver_lp.numCol_]; el++) {
    int row = solver_lp.Aindex_[el];
    if (row >= firstrow && row <= lastrow) XARlength[row - firstrow] += 1;
  }
  XARstart[0] = 0;
  // printf("Start of row %2d is %d\n", 0, XARstart[0]);
  // printf("Length of row %2d is %d\n", 0, XARlength[0]);
  for (int row = 0; row < numExtractRows - 1; row++) {
    XARstart[row + 1] = XARstart[row] + XARlength[row];
    XARlength[row] = 0;
    // printf("Start of row %2d is %d\n", row+1, XARstart[row+1]);
    // printf("Length of row %2d is %d\n", row+1, XARlength[row+1]);
  }
  XARlength[numExtractRows - 1] = 0;

  for (int col = 0; col < solver_lp.numCol_; col++) {
    for (int el = solver_lp.Astart_[col]; el < solver_lp.Astart_[col + 1]; el++) {
      int row = solver_lp.Aindex_[el];
      // printf("Is row=%d in [%d, %d]?\n", row, firstrow, lastrow);
      if (row >= firstrow && row <= lastrow) {
        int rowEl = XARstart[row - firstrow] + XARlength[row - firstrow];
        // printf("Column %2d: Extracted element %d with value %g\n", col,
        // rowEl, solver_lp.Avalue_[el]);
        XARlength[row - firstrow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = solver_lp.Avalue_[el];
      }
    }
  }
  *nnonz = XARstart[lastrow - firstrow] + XARlength[lastrow - firstrow];
  //  printf("Set nnonz = %d\n", *nnonz);
}

// Change a single coefficient in the matrix
void HighsSimplexInterface::util_change_coefficient(int row, int col, const double newval) {
  HighsLp &solver_lp = highs_model_object.solver_lp_;
  assert(row >= 0 && row < solver_lp.numRow_);
  assert(col >= 0 && col < solver_lp.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_changeCoeff(row=%d, col=%d, newval=%g)\n", row, col,
         newval);
  cout << flush;
#endif
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n",
  //  row, col, newval);cout << flush;

  //  solver_lp.reportLp();
  int cg_el = -1;
  for (int el = solver_lp.Astart_[col]; el < solver_lp.Astart_[col + 1]; el++) {
    //    printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //    solver_lp.Aindex_[el], row);
    if (solver_lp.Aindex_[el] == row) {
      cg_el = el;
      break;
    }
  }
  if (cg_el < 0) {
    //    printf("model.util_changeCoeff: Cannot find row %d in column %d\n",
    //    row, col);
    cg_el = solver_lp.Astart_[col + 1];
    int nwNnonz = solver_lp.Astart_[solver_lp.numCol_] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n",
    //    solver_lp.Astart_[solver_lp.numCol_], nwNnonz);
    solver_lp.Aindex_.resize(nwNnonz);
    solver_lp.Avalue_.resize(nwNnonz);
    for (int i = col + 1; i <= solver_lp.numCol_; i++) solver_lp.Astart_[i]++;
    for (int el = nwNnonz - 1; el > cg_el; el--) {
      solver_lp.Aindex_[el] = solver_lp.Aindex_[el - 1];
      solver_lp.Avalue_[el] = solver_lp.Avalue_[el - 1];
    }
  }
  solver_lp.Avalue_[cg_el] = newval;

  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
  //  update_solver_lp_status_flags(highs_model, LpAction::NEW_ROWS);
  //  solver_lp.reportLp();
}

// Get a single coefficient from the matrix
void HighsSimplexInterface::util_get_coefficient(HighsLp lp, int row, int col, double *val) {
  assert(row >= 0 && row < lp.numRow_);
  assert(col >= 0 && col < lp.numCol_);
#ifdef HiGHSDEV
  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);
#endif
  //  printf("Called model.util_getCoeff(row=%d, col=%d)\n", row, col);

  int get_el = -1;
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    //  printf("Column %4d: Element %4d is row %4d. Is it %4d?\n", col, el,
    //  lp.Aindex_[el], row);cout << flush;
    if (lp.Aindex_[el] == row) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    //  printf("model.util_getCoeff: Cannot find row %d in column %d\n", row,
    //  col);cout << flush;
    *val = 0;
  } else {
    //  printf("model.util_getCoeff: Found row %d in column %d as element %d:
    //  value %g\n", row, col, get_el, lp.Avalue_[get_el]);cout << flush;
    *val = lp.Avalue_[get_el];
  }
}
