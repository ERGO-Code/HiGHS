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
#include "simplex/HighsSimplexInterface.h"
#include "HConfig.h"
#include "io/HMPSIO.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
#include "simplex/HSimplex.h"
#include "util/HighsUtils.h"

HighsStatus HighsSimplexInterface::addCols(
    int XnumNewCol, const double* XcolCost, const double* XcolLower,
    const double* XcolUpper, int XnumNewNZ, const int* XAstart,
    const int* XAindex, const double* XAvalue) {
#ifdef HiGHSDEV
  printf("Called addCols(XnumNewCol=%d, XnumNewNZ = %d)\n", XnumNewCol,
         XnumNewNZ);
#endif
  HighsStatus return_status = HighsStatus::NotSet;
  if (XnumNewCol < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewCol == 0) return HighsStatus::OK;
  if (XnumNewCol > 0)
    if (isColDataNull(XcolCost, XcolLower, XcolUpper))
      return HighsStatus::Error;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(XAstart, XAindex, XAvalue)) return HighsStatus::Error;

  HighsLp& lp = highs_model_object.lp_;
  HighsOptions& options = highs_model_object.options_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsScale& scale = highs_model_object.scale_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_lp_matrix = true;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_basis = simplex_lp_status.has_basis;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
  bool apply_row_scaling = scale.is_scaled_;

  // Check that if nonzeros are to be added then the model has a positive number
  // of rows
  if (lp.numRow_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numRow_ <= 0 && XnumNewNZ > 0))
    return HighsStatus::Error;

  // Record the new number of columns
  int newNumCol = lp.numCol_ + XnumNewCol;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or
  // scaling
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
    assert(!apply_row_scaling);
  }
#endif
  HighsStatus call_status;
  call_status =
      appendLpCols(lp, XnumNewCol, XcolCost, XcolLower, XcolUpper, XnumNewNZ,
                   XAstart, XAindex, XAvalue, options, valid_lp_matrix);
  return_status = worseStatus(call_status, return_status);
  if (return_status == HighsStatus::Error) return return_status;

  if (valid_simplex_lp) {
    call_status = appendLpCols(simplex_lp, XnumNewCol, XcolCost, XcolLower,
                               XcolUpper, XnumNewNZ, XAstart, XAindex, XAvalue,
                               options, valid_simplex_matrix);
    return_status = worseStatus(call_status, return_status);
  }

  // Now consider scaling
  scale.col_.resize(newNumCol);
  for (int col = 0; col < XnumNewCol; col++) scale.col_[lp.numCol_ + col] = 1.0;

  if (apply_row_scaling) {
    // Determine scaling multipliers for this set of columns
    // Determine scale factors for this set of columns
    // Scale the simplex LP vectors for these columns
    // Scale the simplex LP matrix for these columns
  }

  // Update the basis correponding to new nonbasic columns
  if (valid_basis) append_nonbasic_cols_to_basis(lp, basis, newNumCol);
  if (valid_simplex_basis)
    append_nonbasic_cols_to_basis(simplex_lp, simplex_basis, newNumCol);

  // Deduce the consequences of adding new columns
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_COLS);

  // Increase the number of columns in the LPs
  lp.numCol_ += XnumNewCol;
  if (valid_simplex_lp) simplex_lp.numCol_ += XnumNewCol;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basisOK = highs_basis_ok();//lp, basis);
    assert(basisOK);
    report_basis(lp, basis);
  }
  if (valid_simplex_basis) {
    bool simplex_basisOK =
        nonbasic_flag_basic_index_ok(simplex_lp, simplex_basis);
    assert(simplex_basisOK);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  return return_status;
}

HighsStatus HighsSimplexInterface::deleteCols(int from_col, int to_col) {
  return deleteColsGeneral(true, from_col, to_col, false, 0, NULL, false, NULL);
}

HighsStatus HighsSimplexInterface::deleteCols(int num_set_entries,
                                              const int* col_set) {
  return deleteColsGeneral(false, 0, 0, true, num_set_entries, col_set, false,
                           NULL);
}

HighsStatus HighsSimplexInterface::deleteCols(int* col_mask) {
  return deleteColsGeneral(false, 0, 0, false, 0, NULL, true, col_mask);
}

HighsStatus HighsSimplexInterface::deleteColsGeneral(
    bool interval, int from_col, int to_col, bool set, int num_set_entries,
    const int* col_set, bool mask, int* col_mask) {
  // Uses to_col in iterator style
  HighsLp& lp = highs_model_object.lp_;

  HighsBasis& basis = highs_model_object.basis_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
  }
#endif
  bool valid_matrix = true;
  HighsStatus returnStatus;
  returnStatus =
      deleteLpCols(lp, interval, from_col, to_col, set, num_set_entries,
                   col_set, mask, col_mask, valid_matrix);
  if (returnStatus != HighsStatus::OK) return returnStatus;
  // ToDo Determine consequences for basis when deleting columns
  basis.valid_ = false;

  if (valid_simplex_lp) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    //  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
    returnStatus = deleteLpCols(simplex_lp, interval, from_col, to_col, set,
                                num_set_entries, col_set, mask, col_mask,
                                valid_simplex_matrix);
    if (returnStatus != HighsStatus::OK) return returnStatus;
    //    HighsScale& scale = highs_model_object.scale_;
    //    for (int col = from_col; col < lp.numCol_ - numDeleteCol; col++)
    //    scale.col_[col] = scale.col_[col + numDeleteCol];
    // ToDo Determine consequences for basis when deleting columns
    simplex_lp_status.has_matrix_col_wise = false;
    simplex_lp_status.has_matrix_row_wise = false;
    simplex_lp_status.has_basis = false;
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::addRows(int XnumNewRow,
                                           const double* XrowLower,
                                           const double* XrowUpper,
                                           int XnumNewNZ, const int* XARstart,
                                           const int* XARindex,
                                           const double* XARvalue) {
#ifdef HiGHSDEV
  printf("Called addRows(XnumNewRow=%d, XnumNewNZ = %d)\n", XnumNewRow,
         XnumNewNZ);
#endif
  HighsStatus return_status = HighsStatus::NotSet;
  if (XnumNewRow < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewRow == 0) return HighsStatus::OK;
  if (XnumNewRow > 0)
    if (isRowDataNull(XrowLower, XrowUpper)) return HighsStatus::Error;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(XARstart, XARindex, XARvalue))
      return HighsStatus::Error;

  HighsLp& lp = highs_model_object.lp_;
  HighsOptions& options = highs_model_object.options_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsScale& scale = highs_model_object.scale_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
  bool apply_row_scaling = scale.is_scaled_;

  // Check that if nonzeros are to be added then the model has a positive number
  // of columns
  if (lp.numCol_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numCol_ <= 0 && XnumNewNZ > 0))
    return HighsStatus::Error;

  // Record the new number of rows
  int newNumRow = lp.numRow_ + XnumNewRow;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or
  // scaling
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
    assert(!apply_row_scaling);
  }
#endif
  // Assess the bounds and matrix indices, returning on error
  bool normalise = false;
  HighsStatus call_status;
  call_status =
      assessBounds("Row", lp.numRow_, XnumNewRow, true, 0, XnumNewRow, false, 0,
                   NULL, false, NULL, (double*)XrowLower, (double*)XrowUpper,
                   options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);

  if (XnumNewNZ) {
    call_status = assessMatrix(lp.numCol_, 0, XnumNewRow, XnumNewRow, XnumNewNZ,
                               (int*)XARstart, (int*)XARindex,
                               (double*)XARvalue, options.small_matrix_value,
                               options.large_matrix_value, normalise);
    return_status = worseStatus(call_status, return_status);
    if (return_status == HighsStatus::Error) return return_status;
  }

  // Append the columns to the LP vectors and matrix
  appendRowsToLpVectors(lp, XnumNewRow, XrowLower, XrowUpper);

  // Normalise the LP row bounds
  normalise = true;
  call_status =
      assessBounds("Row", lp.numRow_, newNumRow, true, 0, newNumRow, false, 0,
                   NULL, false, NULL, &lp.rowLower_[0], &lp.rowUpper_[0],
                   options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);

  int lc_XnumNewNZ = XnumNewNZ;
  int* lc_XARstart = (int*)malloc(sizeof(int) * XnumNewRow);
  int* lc_XARindex = (int*)malloc(sizeof(int) * XnumNewNZ);
  double* lc_XARvalue = (double*)malloc(sizeof(double) * XnumNewNZ);
  if (XnumNewNZ) {
    // Copy the new row-wise matrix into a local copy that can be normalised
    std::memcpy(lc_XARstart, XARstart, sizeof(int) * XnumNewRow);
    std::memcpy(lc_XARindex, XARindex, sizeof(int) * XnumNewNZ);
    std::memcpy(lc_XARvalue, XARvalue, sizeof(double) * XnumNewNZ);
    // Normalise the new matrix columns
    normalise = true;
    call_status = assessMatrix(lp.numCol_, 0, XnumNewRow, XnumNewRow,
                               lc_XnumNewNZ, lc_XARstart, lc_XARindex,
                               lc_XARvalue, options.small_matrix_value,
                               options.large_matrix_value, normalise);
    if (lc_XnumNewNZ) {
      // Append rows to LP matrix
      appendRowsToLpMatrix(lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart,
                           lc_XARindex, lc_XARvalue);
    }
  }

  if (valid_simplex_lp) {
    appendRowsToLpVectors(simplex_lp, XnumNewRow, XrowLower, XrowUpper);
    call_status = assessBounds(
        "Row", simplex_lp.numRow_, newNumRow, true, 0, newNumRow, false, 0,
        NULL, false, NULL, &simplex_lp.colLower_[0], &simplex_lp.colUpper_[0],
        options.infinite_bound, normalise);
    return_status = worseStatus(call_status, return_status);
  }
  if (valid_simplex_matrix && lc_XnumNewNZ) {
    appendRowsToLpMatrix(simplex_lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart,
                         lc_XARindex, lc_XARvalue);
  }

  // Now consider scaling
  scale.row_.resize(newNumRow);
  for (int row = 0; row < XnumNewRow; row++) scale.row_[lp.numRow_ + row] = 1.0;

  if (apply_row_scaling) {
    // Determine scaling multipliers for this set of rows
    // Determine scale factors for this set of rows
    // Scale the simplex LP vectors for these rows
    // Scale the simplex LP matrix for these rows
  }

  // Update the basis correponding to new basic rows
  if (valid_basis) append_basic_rows_to_basis(lp, basis, newNumRow);
  //  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  //  if (simplex_lp_status.has_basis) append_basic_rows_to_basis(simplex_lp,
  //  simplex_basis, newNumRow);

  // Deduce the consequences of adding new rows
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_ROWS);

  // Increase the number of rows in the LPs
  lp.numRow_ += XnumNewRow;
  if (valid_simplex_lp) simplex_lp.numRow_ += XnumNewRow;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basisOK = highs_basis_ok();//lp, basis);
    if (!basisOK) printf("HiGHS basis not OK in addRows\n");
    assert(basisOK);
    report_basis(lp, basis);
  }
  if (simplex_lp_status.has_basis) {
    SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
    bool simplex_basisOK =
        nonbasic_flag_basic_index_ok(simplex_lp, simplex_basis);
    if (!simplex_basisOK) printf("Simplex basis not OK in addRows\n");
    assert(simplex_basisOK);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  return return_status;
}

HighsStatus HighsSimplexInterface::deleteRows(int from_row, int to_row) {
  return deleteRowsGeneral(true, from_row, to_row, false, 0, NULL, false, NULL);
}

HighsStatus HighsSimplexInterface::deleteRows(int num_set_entries,
                                              const int* row_set) {
  return deleteRowsGeneral(false, 0, 0, true, num_set_entries, row_set, false,
                           NULL);
}

HighsStatus HighsSimplexInterface::deleteRows(int* row_mask) {
  return deleteRowsGeneral(false, 0, 0, false, 0, NULL, true, row_mask);
}

HighsStatus HighsSimplexInterface::deleteRowsGeneral(
    bool interval, int from_row, int to_row, bool set, int num_set_entries,
    const int* row_set, bool mask, int* row_mask) {
  // Uses to_row in iterator style
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(from_row=%d, to_row=%d)\n", from_row,
         to_row);
#endif
  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
  }
#endif
  bool valid_matrix = true;
  HighsStatus returnStatus;
  returnStatus =
      deleteLpRows(lp, interval, from_row, to_row, set, num_set_entries,
                   row_set, mask, row_mask, valid_matrix);
  if (returnStatus != HighsStatus::OK) return returnStatus;
  // ToDo Determine consequences for basis when deleting rows
  basis.valid_ = false;

  if (valid_simplex_lp) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    //    SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
    returnStatus = deleteLpRows(simplex_lp, interval, from_row, to_row, set,
                                num_set_entries, row_set, mask, row_mask,
                                valid_simplex_matrix);
    if (returnStatus != HighsStatus::OK) return returnStatus;
    //    HighsScale& scale = highs_model_object.scale_;
    //    for (int row = from_row; row < lp.numRow_ - numDeleteRow; row++)
    //    scale.row_[row] = scale.row_[row + numDeleteRow];
    // ToDo Determine consequences for basis when deleting rows
    simplex_lp_status.has_matrix_col_wise = false;
    simplex_lp_status.has_matrix_row_wise = false;
    simplex_lp_status.has_basis = false;
  }
  // ToDo Determine consequences for basis when deleting rows
  //  updateSimplexLpStatus(simplex_lp_status, LpAction::DEL_ROWS);
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::getCols(const int from_col, const int to_col,
                                           int& num_col, double* col_cost,
                                           double* col_lower, double* col_upper,
                                           int& num_nz, int* col_matrix_start,
                                           int* col_matrix_index,
                                           double* col_matrix_value) {
  return getColsGeneral(true, from_col, to_col, false, 0, NULL, false, NULL,
                        num_col, col_cost, col_lower, col_upper, num_nz,
                        col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getCols(
    const int num_set_entries, const int* col_set, int& num_col,
    double* col_cost, double* col_lower, double* col_upper, int& num_nz,
    int* col_matrix_start, int* col_matrix_index, double* col_matrix_value) {
  return getColsGeneral(false, 0, 0, true, num_set_entries, col_set, false,
                        NULL, num_col, col_cost, col_lower, col_upper, num_nz,
                        col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getCols(const int* col_mask, int& num_col,
                                           double* col_cost, double* col_lower,
                                           double* col_upper, int& num_nz,
                                           int* col_matrix_start,
                                           int* col_matrix_index,
                                           double* col_matrix_value) {
  return getColsGeneral(false, 0, 0, false, 0, NULL, true, col_mask, num_col,
                        col_cost, col_lower, col_upper, num_nz,
                        col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getColsGeneral(
    const bool interval, const int from_col, const int to_col, const bool set,
    const int num_set_entries, const int* col_set, const bool mask,
    const int* col_mask, int& num_col, double* col_cost, double* col_lower,
    double* col_upper, int& num_nz, int* col_matrix_start,
    int* col_matrix_index, double* col_matrix_value) {
  int from_k;
  int to_k;
  HighsLp& lp = highs_model_object.lp_;
  HighsStatus return_status = assessIntervalSetMask(
      lp.numCol_, interval, from_col, to_col, set, num_set_entries, col_set,
      mask, col_mask, from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  if (from_k < 0 || to_k > lp.numCol_) return HighsStatus::Error;
  num_col = 0;
  num_nz = 0;
  if (from_k >= to_k) return HighsStatus::OK;
  int out_from_col;
  int out_to_col;
  int in_from_col;
  int in_to_col = 0;
  int current_set_entry = 0;
  int col_dim = lp.numCol_;
  for (int k = from_k; k < to_k; k++) {
    updateOutInIx(col_dim, interval, from_col, to_col, set, num_set_entries,
                  col_set, mask, col_mask, out_from_col, out_to_col,
                  in_from_col, in_to_col, current_set_entry);
    assert(out_to_col <= col_dim);
    assert(in_to_col <= col_dim);
    for (int col = out_from_col; col < out_to_col; col++) {
      col_cost[num_col] = lp.colCost_[col];
      col_lower[num_col] = lp.colLower_[col];
      col_upper[num_col] = lp.colUpper_[col];
      col_matrix_start[num_col] =
          num_nz + lp.Astart_[col] - lp.Astart_[out_from_col];
      num_col++;
    }
    for (int el = lp.Astart_[out_from_col]; el < lp.Astart_[out_to_col]; el++) {
      col_matrix_index[num_nz] = lp.Aindex_[el];
      col_matrix_value[num_nz] = lp.Avalue_[el];
      num_nz++;
    }
    if (out_to_col == col_dim || in_to_col == col_dim) break;
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::getRows(const int from_row, const int to_row,
                                           int& num_row, double* row_lower,
                                           double* row_upper, int& num_nz,
                                           int* row_matrix_start,
                                           int* row_matrix_index,
                                           double* row_matrix_value) {
  return getRowsGeneral(true, from_row, to_row, false, 0, NULL, false, NULL,
                        num_row, row_lower, row_upper, num_nz, row_matrix_start,
                        row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRows(const int num_set_entries,
                                           const int* row_set, int& num_row,
                                           double* row_lower, double* row_upper,
                                           int& num_nz, int* row_matrix_start,
                                           int* row_matrix_index,
                                           double* row_matrix_value) {
  return getRowsGeneral(false, 0, 0, true, num_set_entries, row_set, false,
                        NULL, num_row, row_lower, row_upper, num_nz,
                        row_matrix_start, row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRows(const int* row_mask, int& num_row,
                                           double* row_lower, double* row_upper,
                                           int& num_nz, int* row_matrix_start,
                                           int* row_matrix_index,
                                           double* row_matrix_value) {
  return getRowsGeneral(false, 0, 0, false, 0, NULL, true, row_mask, num_row,
                        row_lower, row_upper, num_nz, row_matrix_start,
                        row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRowsGeneral(
    const bool interval, const int from_row, const int to_row, const bool set,
    const int num_set_entries, const int* row_set, const bool mask,
    const int* row_mask, int& num_row, double* row_lower, double* row_upper,
    int& num_nz, int* row_matrix_start, int* row_matrix_index,
    double* row_matrix_value) {
  int from_k;
  int to_k;
  HighsLp& lp = highs_model_object.lp_;
  HighsStatus return_status = assessIntervalSetMask(
      lp.numRow_, interval, from_row, to_row, set, num_set_entries, row_set,
      mask, row_mask, from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  if (from_k < 0 || to_k > lp.numRow_) return HighsStatus::Error;
  num_row = 0;
  num_nz = 0;
  if (from_k >= to_k) return HighsStatus::OK;
  // "Out" means not in the set to be extrated
  // "In" means in the set to be extrated
  int out_from_row;
  int out_to_row;
  int in_from_row;
  int in_to_row = 0;
  int current_set_entry = 0;
  int row_dim = lp.numRow_;
  // Set up a row mask so that entries to be got from the column-wise
  // matrix can be identified and have their correct row index.
  int* new_index = (int*)malloc(sizeof(int) * lp.numRow_);

  if (!mask) {
    out_to_row = 0;
    current_set_entry = 0;
    for (int k = from_k; k < to_k; k++) {
      updateOutInIx(row_dim, interval, from_row, to_row, set, num_set_entries,
                    row_set, mask, row_mask, in_from_row, in_to_row,
                    out_from_row, out_to_row, current_set_entry);
      if (k == from_k) {
        // Account for any initial rows not being extracted
        for (int row = 0; row < in_from_row; row++) {
          new_index[row] = -1;
        }
      }
      for (int row = in_from_row; row < in_to_row; row++) {
        new_index[row] = num_row;
        num_row++;
      }
      for (int row = out_from_row; row < out_to_row; row++) {
        new_index[row] = -1;
      }
      if (out_to_row == row_dim) break;
    }
  } else {
    for (int row = 0; row < lp.numRow_; row++) {
      if (row_mask[row]) {
        new_index[row] = num_row;
        num_row++;
      } else {
        new_index[row] = -1;
      }
    }
  }

  // Bail out if no rows are to be extracted
  if (num_row == 0) return HighsStatus::OK;

  // Allocate an array of lengths for the row-wise matrix to be extracted
  int* row_matrix_length = (int*)malloc(sizeof(int) * num_row);

  for (int row = 0; row < lp.numRow_; row++) {
    int new_row = new_index[row];
    if (new_row >= 0) {
      assert(new_row < num_row);
      row_lower[new_row] = lp.rowLower_[row];
      row_upper[new_row] = lp.rowUpper_[row];
      row_matrix_length[new_row] = 0;
    }
  }
  // Identify the lengths of the rows in the row-wise matrix to be extracted
  for (int col = 0; col < lp.numCol_; col++) {
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      int new_row = new_index[row];
      if (new_row >= 0) row_matrix_length[new_row]++;
    }
  }

  row_matrix_start[0] = 0;
  for (int row = 0; row < num_row - 1; row++) {
    row_matrix_start[row + 1] = row_matrix_start[row] + row_matrix_length[row];
  }

  // Fill the row-wise matrix with indices and values
  for (int col = 0; col < lp.numCol_; col++) {
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      int new_row = new_index[row];
      if (new_row >= 0) {
        int row_el = row_matrix_start[new_row];
        row_matrix_index[row_el] = col;
        row_matrix_value[row_el] = lp.Avalue_[el];
        row_matrix_start[new_row]++;
      }
    }
  }
  // Restore the starts of the row-wise matrix and count the number of nonzeros
  // in it
  num_nz = 0;
  row_matrix_start[0] = 0;
  for (int row = 0; row < num_row - 1; row++) {
    row_matrix_start[row + 1] = row_matrix_start[row] + row_matrix_length[row];
    num_nz += row_matrix_length[row];
  }
  num_nz += row_matrix_length[num_row - 1];
  return HighsStatus::OK;
}

// Change a single coefficient in the matrix
HighsStatus HighsSimplexInterface::changeCoefficient(int Xrow, int Xcol,
                                                     const double XnewValue) {
#ifdef HiGHSDEV
  printf("Called util_changeCoeff(Xrow=%d, Xcol=%d, XnewValue=%g)\n", Xrow,
         Xcol, XnewValue);
#endif
  HighsLp& lp = highs_model_object.lp_;
  if (Xrow < 0 || Xrow > lp.numRow_) return HighsStatus::Error;
  if (Xcol < 0 || Xcol > lp.numCol_) return HighsStatus::Error;
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n",
  //  Xrow, Xcol, XnewValue);
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  bool valid_simplex_lp = simplex_lp_status.valid;
#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix or scaling
  if (!valid_simplex_lp) {
    assert(!simplex_lp_status.has_matrix_col_wise);
    //    assert(!apply_row_scaling);
  }
#endif
  changeLpMatrixCoefficient(lp, Xrow, Xcol, XnewValue);
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    HighsScale& scale = highs_model_object.scale_;
    double scaledXnewValue = XnewValue * scale.row_[Xrow] * scale.col_[Xcol];
    changeLpMatrixCoefficient(simplex_lp, Xrow, Xcol, scaledXnewValue);
  }
  // simplex_lp.reportLp();
  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_ROWS);
  //  simplex_lp.reportLp();
  return HighsStatus::OK;
}

void HighsSimplexInterface::shiftObjectiveValue(double Xshift) {
  printf(
      "Where is shiftObjectiveValue required - so I can interpret what's "
      "required\n");
  // Update the LP objective value with the shift
  highs_model_object.simplex_info_.dual_objective_value += Xshift;
  // Update the LP offset with the shift
  highs_model_object.lp_.offset_ += Xshift;
  if (highs_model_object.simplex_lp_status_.valid) {
    // Update the simplex LP offset with the shift
    highs_model_object.simplex_lp_.offset_ += Xshift;
  }
}

HighsStatus HighsSimplexInterface::changeObjectiveSense(int Xsense) {
  HighsLp& lp = highs_model_object.lp_;
  if ((Xsense == OBJSENSE_MINIMIZE) != (lp.sense_ == OBJSENSE_MINIMIZE)) {
    // Flip the LP objective sense
    lp.sense_ = Xsense;
  }
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    HighsLp& simplex_lp = highs_model_object.simplex_lp_;
    if ((Xsense == OBJSENSE_MINIMIZE) !=
        (simplex_lp.sense_ == OBJSENSE_MINIMIZE)) {
      // Flip the objective sense
      simplex_lp.sense_ = Xsense;
      simplex_lp_status.solution_status = SimplexSolutionStatus::UNSET;
    }
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::changeCosts(int from_col, int to_col,
                                               const double* usr_col_cost) {
  return changeCostsGeneral(true, from_col, to_col, false, 0, NULL, false, NULL,
                            usr_col_cost);
}

HighsStatus HighsSimplexInterface::changeCosts(int num_set_entries,
                                               const int* col_set,
                                               const double* usr_col_cost) {
  return changeCostsGeneral(false, 0, 0, true, num_set_entries, col_set, false,
                            NULL, usr_col_cost);
}

HighsStatus HighsSimplexInterface::changeCosts(const int* col_mask,
                                               const double* usr_col_cost) {
  return changeCostsGeneral(false, 0, 0, false, 0, NULL, true, col_mask,
                            usr_col_cost);
}

HighsStatus HighsSimplexInterface::changeCostsGeneral(
    bool interval, int from_col, int to_col, bool set, int num_set_entries,
    const int* col_set, bool mask, const int* col_mask,
    const double* usr_col_cost) {
  bool null_data = false;
  if (usr_col_cost == NULL) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "User-supplied column costs are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  HighsStatus call_status =
      changeLpCosts(highs_model_object.lp_, interval, from_col, to_col, set,
                    num_set_entries, col_set, mask, col_mask, usr_col_cost,
                    highs_model_object.options_.infinite_cost);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;
  // Deduce the consequences of new costs
  updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                        LpAction::NEW_COSTS);
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::changeColBounds(
    int from_col, int to_col, const double* usr_col_lower,
    const double* usr_col_upper) {
  return changeColBoundsGeneral(true, from_col, to_col, false, 0, NULL, false,
                                NULL, usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::changeColBounds(
    int num_set_entries, const int* col_set, const double* usr_col_lower,
    const double* usr_col_upper) {
  return changeColBoundsGeneral(false, 0, 0, true, num_set_entries, col_set,
                                false, NULL, usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::changeColBounds(
    const int* col_mask, const double* usr_col_lower,
    const double* usr_col_upper) {
  return changeColBoundsGeneral(false, 0, 0, false, 0, NULL, true, col_mask,
                                usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::changeColBoundsGeneral(
    bool interval, int from_col, int to_col, bool set, int num_set_entries,
    const int* col_set, bool mask, const int* col_mask,
    const double* usr_col_lower, const double* usr_col_upper) {
  bool null_data = false;
  if (usr_col_lower == NULL) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "User-supplied column lower bounds are NULL");
    null_data = true;
  }
  if (usr_col_upper == NULL) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "User-supplied column upper bounds are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  HighsStatus call_status = changeLpColBounds(
      highs_model_object.lp_, interval, from_col, to_col, set, num_set_entries,
      col_set, mask, col_mask, usr_col_lower, usr_col_upper,
      highs_model_object.options_.infinite_bound);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;

  if (highs_model_object.simplex_lp_status_.valid) {
    // Also change the simplex LP's column bounds
    assert(highs_model_object.lp_.numCol_ ==
           highs_model_object.simplex_lp_.numCol_);
    assert(highs_model_object.lp_.numRow_ ==
           highs_model_object.simplex_lp_.numRow_);

    call_status = changeLpColBounds(
        highs_model_object.simplex_lp_, interval, from_col, to_col, set,
        num_set_entries, col_set, mask, col_mask, usr_col_lower, usr_col_upper,
        highs_model_object.options_.infinite_bound);
    if (call_status == HighsStatus::Error) return HighsStatus::Error;
    if (highs_model_object.scale_.is_scaled_) {
      scaleLpColBounds(highs_model_object.simplex_lp_,
                       highs_model_object.scale_.col_, interval, from_col,
                       to_col, set, num_set_entries, col_set, mask, col_mask);
    }
    // Deduce the consequences of new bounds
    updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                          LpAction::NEW_BOUNDS);
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::changeRowBounds(
    int from_row, int to_row, const double* usr_row_lower,
    const double* usr_row_upper) {
  return changeRowBoundsGeneral(true, from_row, to_row, false, 0, NULL, false,
                                NULL, usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::changeRowBounds(
    int num_set_entries, const int* row_set, const double* usr_row_lower,
    const double* usr_row_upper) {
  return changeRowBoundsGeneral(false, 0, 0, true, num_set_entries, row_set,
                                false, NULL, usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::changeRowBounds(
    const int* row_mask, const double* usr_row_lower,
    const double* usr_row_upper) {
  return changeRowBoundsGeneral(false, 0, 0, false, 0, NULL, true, row_mask,
                                usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::changeRowBoundsGeneral(
    bool interval, int from_row, int to_row, bool set, int num_set_entries,
    const int* row_set, bool mask, const int* row_mask,
    const double* usr_row_lower, const double* usr_row_upper) {
  bool null_data = false;
  if (usr_row_lower == NULL) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "User-supplied row lower bounds are NULL");
    null_data = true;
  }
  if (usr_row_upper == NULL) {
    HighsLogMessage(HighsMessageType::ERROR,
                    "User-supplied row upper bounds are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  HighsStatus call_status = changeLpRowBounds(
      highs_model_object.lp_, interval, from_row, to_row, set, num_set_entries,
      row_set, mask, row_mask, usr_row_lower, usr_row_upper,
      highs_model_object.options_.infinite_bound);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;
  if (highs_model_object.simplex_lp_status_.valid) {
    // Also change the simplex LP's column bounds
    assert(highs_model_object.lp_.numCol_ ==
           highs_model_object.simplex_lp_.numCol_);
    assert(highs_model_object.lp_.numRow_ ==
           highs_model_object.simplex_lp_.numRow_);

    call_status = changeLpRowBounds(
        highs_model_object.simplex_lp_, interval, from_row, to_row, set,
        num_set_entries, row_set, mask, row_mask, usr_row_lower, usr_row_upper,
        highs_model_object.options_.infinite_bound);
    if (call_status == HighsStatus::Error) return HighsStatus::Error;
    if (highs_model_object.scale_.is_scaled_) {
      scaleLpRowBounds(highs_model_object.simplex_lp_,
                       highs_model_object.scale_.row_, interval, from_row,
                       to_row, set, num_set_entries, row_set, mask, row_mask);
    }
    // Deduce the consequences of new bounds
    updateSimplexLpStatus(highs_model_object.simplex_lp_status_,
                          LpAction::NEW_BOUNDS);
  }
  return HighsStatus::OK;
}

#ifdef HiGHSDEV
void HighsSimplexInterface::change_update_method(int updateMethod) {
  highs_model_object.factor_.change(updateMethod);
}
#endif

HighsStatus HighsSimplexInterface::lpStatusToHighsStatus(
    SimplexSolutionStatus simplex_solution_status) {
  switch (simplex_solution_status) {
    case SimplexSolutionStatus::OUT_OF_TIME:
      return HighsStatus::Timeout;
    case SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return HighsStatus::ReachedDualObjectiveUpperBound;
    case SimplexSolutionStatus::FAILED:
      return HighsStatus::SolutionError;
    case SimplexSolutionStatus::SINGULAR:
      return HighsStatus::SolutionError;
    case SimplexSolutionStatus::UNBOUNDED:
      return HighsStatus::Unbounded;
    case SimplexSolutionStatus::INFEASIBLE:
      return HighsStatus::Infeasible;
    case SimplexSolutionStatus::DUAL_FEASIBLE:
      return HighsStatus::DualFeasible;
    case SimplexSolutionStatus::PRIMAL_FEASIBLE:
      return HighsStatus::PrimalFeasible;
    case SimplexSolutionStatus::OPTIMAL:
      return HighsStatus::Optimal;
    default:
      return HighsStatus::NotImplemented;
  }
}

// Utilities to convert model basic/nonbasic status to/from SCIP-like status
int HighsSimplexInterface::convertBaseStatToHighsBasis(const int* cstat,
                                                       const int* rstat) {
  HighsBasis& basis = highs_model_object.basis_;
  HighsLp& lp = highs_model_object.lp_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;

  int numBasic = 0;
  bool error_found = false;
  basis.valid_ = false;
  for (int col = 0; col < lp.numCol_; col++) {
    if (cstat[col] == (int)HighsBasisStatus::BASIC) {
      basis.col_status[col] = HighsBasisStatus::BASIC;
      numBasic++;
      continue;
    }
    if (cstat[col] == (int)HighsBasisStatus::LOWER) {
      // Supplied basis has this column nonbasic at its lower bound: check that
      // the lower bound is finite
      error_found = highs_isInfinity(-lp.colLower_[col]);
      basis.col_status[col] = HighsBasisStatus::LOWER;
    } else if (cstat[col] == (int)HighsBasisStatus::UPPER) {
      // Supplied basis has this column nonbasic at its upper bound: check that
      // the upper bound is finite
      error_found = highs_isInfinity(lp.colUpper_[col]);
      basis.col_status[col] = HighsBasisStatus::UPPER;
    } else if (cstat[col] == (int)HighsBasisStatus::ZERO) {
      // Supplied basis has this column nonbasic at zero so free: check that
      // neither bound is finite
      error_found = !highs_isInfinity(-lp.colLower_[col]) ||
                    !highs_isInfinity(lp.colUpper_[col]);
      basis.col_status[col] = HighsBasisStatus::UPPER;
    } else {
      error_found = true;
    }
    if (error_found) {
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, cstat=%d, lower=%g, upper=%g\n",
             col, cstat[col], lp.colLower_[col], lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < lp.numRow_; row++) {
    if (rstat[row] == (int)HighsBasisStatus::BASIC) {
      basis.row_status[row] = HighsBasisStatus::BASIC;
      numBasic++;
      continue;
    }
    if (rstat[row] == (int)HighsBasisStatus::LOWER) {
      // Supplied basis has this row nonbasic at its lower bound: check that the
      // lower bound is finite
      error_found = highs_isInfinity(-lp.rowLower_[row]);
      basis.row_status[row] = HighsBasisStatus::LOWER;
    } else if (rstat[row] == (int)HighsBasisStatus::UPPER) {
      // Supplied basis has this row nonbasic at its upper bound: check that the
      // upper bound is finite
      error_found = highs_isInfinity(lp.rowUpper_[row]);
      basis.row_status[row] = HighsBasisStatus::UPPER;
    } else if (rstat[row] == (int)HighsBasisStatus::ZERO) {
      // Supplied basis has this row nonbasic at zero so free: check that
      // neither bound is finite
      error_found = !highs_isInfinity(-lp.rowLower_[row]) ||
                    !highs_isInfinity(lp.rowUpper_[row]);
      basis.row_status[row] = HighsBasisStatus::UPPER;
    } else {
      error_found = true;
    }
    if (error_found) {
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], lp.rowLower_[row], lp.rowUpper_[row]);
#endif
      return -(row + 1);
    }
  }
  assert(numBasic = lp.numRow_);
  basis.valid_ = true;
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
  return 0;
}

int HighsSimplexInterface::convertHighsBasisToBaseStat(int* cstat, int* rstat) {
  HighsBasis& basis = highs_model_object.basis_;
  HighsLp& lp = highs_model_object.lp_;
  if (cstat != NULL) {
    for (int col = 0; col < lp.numCol_; col++)
      cstat[col] = (int)basis.col_status[col];
  }
  printf("NB SCIP has row bounds [-u, -l]\n");
  if (rstat != NULL) {
    for (int row = 0; row < lp.numRow_; row++)
      rstat[row] = (int)basis.row_status[row];
  }
  return 0;
}

void HighsSimplexInterface::convertSimplexToHighsBasis() {
  HighsBasis& basis = highs_model_object.basis_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsLp& lp = highs_model_object.lp_;
  //  HighsSimplexLpStatus &simplex_lp_status =
  //  highs_model_object.simplex_lp_status_; HighsSimplexInfo &simplex_info =
  //  highs_model_object.simplex_info_;
  basis.col_status.resize(lp.numCol_);
  basis.row_status.resize(lp.numRow_);

  assert(highs_model_object.simplex_lp_status_.has_basis);
  bool permuted = highs_model_object.simplex_lp_status_.is_permuted;
  int* numColPermutation =
      &highs_model_object.simplex_info_.numColPermutation_[0];
  // numColPermutation[iCol] is the true column in column iCol
  bool error_found = false;
  basis.valid_ = false;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    int simplex_var = iCol;
    int lp_col = iCol;
    if (permuted) lp_col = numColPermutation[iCol];
    if (!simplex_basis.nonbasicFlag_[simplex_var]) {
      basis.col_status[lp_col] = HighsBasisStatus::BASIC;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_UP) {
      // Nonbasic and free to move up so should be OK to give status LOWER
#ifdef HiGHSDEV
      // Check that the lower bound isn't infinite
      error_found = highs_isInfinity(-lp.colLower_[lp_col]);
#endif
      basis.col_status[lp_col] = HighsBasisStatus::LOWER;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_DN) {
      // Nonbasic and free to move down so should be OK to give status UPPER
#ifdef HiGHSDEV
      // Check that the upper bound isn't infinite
      error_found = highs_isInfinity(lp.colUpper_[lp_col]);
#endif
      basis.col_status[lp_col] = HighsBasisStatus::UPPER;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_ZE) {
      // Column is either fixed or free, depending on the bounds
      if (lp.colLower_[lp_col] == lp.colUpper_[lp_col]) {
        // Equal bounds so should be OK to give status LOWER
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lp.colLower_[lp_col]);
#endif
        basis.col_status[lp_col] = HighsBasisStatus::LOWER;
      } else {
        // Unequal bounds so should be OK to give status ZERO
#ifdef HiGHSDEV
        // Check that neither of the bounds is finite
        error_found = !highs_isInfinity(-lp.colLower_[lp_col]) ||
                      !highs_isInfinity(lp.colUpper_[lp_col]);
#endif
        basis.col_status[lp_col] = HighsBasisStatus::ZERO;
      }
    } else {
      error_found = true;
    }
#ifdef HiGHSDEV
    if (error_found)
      printf(
          "Invalid basis status: col=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          lp_col, simplex_basis.nonbasicFlag_[simplex_var],
          simplex_basis.nonbasicMove_[simplex_var], lp.colLower_[lp_col],
          lp.colUpper_[lp_col]);
#endif
    assert(!error_found);
    if (error_found) return;
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    int simplex_var = lp.numCol_ + iRow;
    int lp_row = iRow;
    if (!simplex_basis.nonbasicFlag_[simplex_var]) {
      basis.row_status[lp_row] = HighsBasisStatus::BASIC;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_UP) {
      // Nonbasic and free to move up so should be OK to give status UPPER -
      // since simplex row bounds are flipped and negated
#ifdef HiGHSDEV
      // Check that the upper bound isn't infinite
      error_found = highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
      basis.row_status[lp_row] = HighsBasisStatus::UPPER;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_DN) {
      // Nonbasic and free to move down so should be OK to give status
      // LOWER - since simplex row bounds are flipped and negated
#ifdef HiGHSDEV
      // Check that the lower bound isn't infinite
      error_found = highs_isInfinity(-lp.rowLower_[lp_row]);
#endif
      basis.row_status[lp_row] = HighsBasisStatus::LOWER;
    } else if (simplex_basis.nonbasicMove_[simplex_var] == NONBASIC_MOVE_ZE) {
      // Row is either fixed or free, depending on the bounds
      if (lp.rowLower_[lp_row] == lp.rowUpper_[lp_row]) {
        // Equal bounds so should be OK to give status LOWER
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lp.rowLower_[lp_row]);
#endif
        basis.row_status[lp_row] = HighsBasisStatus::LOWER;
      } else {
        // Unequal bounds so should be OK to give status ZERO
#ifdef HiGHSDEV
        // Check that neither of the bounds is finite
        error_found = !highs_isInfinity(-lp.rowLower_[lp_row]) ||
                      !highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
        basis.row_status[lp_row] = HighsBasisStatus::ZERO;
      }
    } else {
      error_found = true;
    }
#ifdef HiGHSDEV
    if (error_found)
      printf(
          "Invalid basis status: row=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          lp_row, simplex_basis.nonbasicFlag_[simplex_var],
          simplex_basis.nonbasicMove_[simplex_var], lp.rowLower_[lp_row],
          lp.rowUpper_[lp_row]);
#endif
    assert(!error_found);
    if (error_found) return;
  }
  basis.valid_ = true;
}

void HighsSimplexInterface::convertHighsToSimplexBasis() {
  HighsBasis& basis = highs_model_object.basis_;
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsLp& lp = highs_model_object.lp_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  bool error_found = false;
  bool permuted = highs_model_object.simplex_lp_status_.is_permuted;
  int* numColPermutation =
      &highs_model_object.simplex_info_.numColPermutation_[0];
  // numColPermutation[iCol] is the true column in column iCol
  int num_basic = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    int simplex_var = iCol;
    int lp_col = iCol;
    if (permuted) lp_col = numColPermutation[iCol];
    if (basis.col_status[lp_col] == HighsBasisStatus::BASIC) {
      // Basic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_FALSE;
      simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      simplex_basis.basicIndex_[num_basic] = simplex_var;
      num_basic++;
    } else {
      // Nonbasic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_TRUE;
      if (basis.col_status[lp_col] == HighsBasisStatus::LOWER) {
        // HighsBasisStatus::LOWER includes fixed variables
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lp.colLower_[lp_col]);
#endif
        if (lp.colLower_[lp_col] == lp.colUpper_[lp_col]) {
          // Equal bounds so indicate that the column can't move
#ifdef HiGHSDEV
          // Check that the upper bound isn't infinite
          error_found = highs_isInfinity(lp.colUpper_[lp_col]);
#endif
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
        } else {
          // unequal bounds so indicate that the column can only move up
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_UP;
        }
      } else if (basis.col_status[lp_col] == HighsBasisStatus::UPPER) {
        // HighsBasisStatus::UPPER includes only variables at their upper bound
#ifdef HiGHSDEV
        // Check that the upper bound isn't infinite
        error_found = highs_isInfinity(lp.colUpper_[lp_col]);
#endif
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_DN;
      } else if (basis.col_status[lp_col] == HighsBasisStatus::ZERO) {
        // HighsBasisStatus::ZERO implies a free variable
#ifdef HiGHSDEV
        // Check that neither bound is finite
        error_found = !highs_isInfinity(-lp.colLower_[lp_col]) ||
                      !highs_isInfinity(lp.colUpper_[lp_col]);
#endif
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      } else {
        error_found = true;
      }
    }
#ifdef HiGHSDEV
    if (error_found)
      printf(
          "Invalid basis status: col=%d, basis.col_status=%d, lower=%g, "
          "upper=%g\n",
          lp_col, (int)basis.col_status[lp_col], lp.colLower_[lp_col],
          lp.colUpper_[lp_col]);
#endif
    assert(!error_found);
    if (error_found) return;
  }
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    int simplex_var = lp.numCol_ + iRow;
    int lp_row = iRow;
    if (basis.row_status[lp_row] == HighsBasisStatus::BASIC) {
      // Basic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_FALSE;
      simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      simplex_basis.basicIndex_[num_basic] = simplex_var;
      num_basic++;
    } else {
      // Nonbasic
      simplex_basis.nonbasicFlag_[simplex_var] = NONBASIC_FLAG_TRUE;
      if (basis.row_status[lp_row] == HighsBasisStatus::LOWER) {
        // HighsBasisStatus::LOWER includes fixed variables
#ifdef HiGHSDEV
        // Check that the lower bound isn't infinite
        error_found = highs_isInfinity(-lp.rowLower_[lp_row]);
#endif
        if (lp.rowLower_[lp_row] == lp.rowUpper_[lp_row]) {
          // Equal bounds so indicate that the row can't move
#ifdef HiGHSDEV
          // Check that the upper bound isn't infinite
          error_found = highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
        } else {
          // Unequal bounds so indicate that the row can only move
          // down - since simplex row bounds are flipped and negated
          simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_DN;
        }
      } else if (basis.row_status[lp_row] == HighsBasisStatus::UPPER) {
        // HighsBasisStatus::UPPER includes only variables at their upper bound
#ifdef HiGHSDEV
        // Check that the upper bound isn't infinite
        error_found = highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
        // Upper bounded so indicate that the row can only move
        // up - since simplex row bounds are flipped and negated
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_UP;
      } else if (basis.row_status[lp_row] == HighsBasisStatus::ZERO) {
        // HighsBasisStatus::ZERO implies a free variable
#ifdef HiGHSDEV
        // Check that neither bound is finite
        error_found = !highs_isInfinity(-lp.rowLower_[lp_row]) ||
                      !highs_isInfinity(lp.rowUpper_[lp_row]);
#endif
        simplex_basis.nonbasicMove_[simplex_var] = NONBASIC_MOVE_ZE;
      } else {
        error_found = true;
      }
    }
#ifdef HiGHSDEV
    if (error_found)
      printf(
          "Invalid basis status: row=%d, basis.row_status=%d, lower=%g, "
          "upper=%g\n",
          lp_row, (int)basis.row_status[lp_row], lp.rowLower_[lp_row],
          lp.rowUpper_[lp_row]);
#endif
    assert(!error_found);
    if (error_found) return;
  }
  assert(num_basic = lp.numRow_);
  //  populate_work_arrays(highs_model_object); // Why might this have been done
  //  here?
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_BASIS);
  simplex_lp_status.has_basis = true;
}

void HighsSimplexInterface::convertSimplexToHighsSolution() {
  HighsSolution& solution = highs_model_object.solution_;
  HighsScale& scale = highs_model_object.scale_;
  SimplexBasis& basis = highs_model_object.simplex_basis_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;

  // Take primal solution
  vector<double> value = simplex_info.workValue_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    value[basis.basicIndex_[iRow]] = simplex_info.baseValue_[iRow];
  // Take dual solution
  vector<double> dual = simplex_info.workDual_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    dual[basis.basicIndex_[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    value[iCol] *= scale.col_[iCol];
    dual[iCol] /= (scale.col_[iCol] / scale.cost_);
  }
  for (int iRow = 0, iTot = simplex_lp.numCol_; iRow < simplex_lp.numRow_;
       iRow++, iTot++) {
    value[iTot] /= scale.row_[iRow];
    dual[iTot] *= (scale.row_[iRow] * scale.cost_);
  }

  // Now we can get the solution
  solution.col_value.resize(simplex_lp.numCol_);
  solution.col_dual.resize(simplex_lp.numCol_);
  solution.row_value.resize(simplex_lp.numRow_);
  solution.row_dual.resize(simplex_lp.numRow_);

  if (highs_model_object.simplex_lp_status_.is_permuted) {
    const int* numColPermutation =
        &highs_model_object.simplex_info_.numColPermutation_[0];
    // numColPermutation[i] is the true column in column i
    for (int i = 0; i < simplex_lp.numCol_; i++) {
      int iCol = numColPermutation[i];
      solution.col_value[iCol] = value[i];
      solution.col_dual[iCol] = simplex_lp.sense_ * dual[i];
    }
  } else {
    for (int i = 0; i < simplex_lp.numCol_; i++) {
      int iCol = i;
      solution.col_value[iCol] = value[i];
      solution.col_dual[iCol] = simplex_lp.sense_ * dual[i];
    }
  }
  int row_dual_sign = 1;  //-1;
  for (int i = 0; i < simplex_lp.numRow_; i++) {
    solution.row_value[i] = -value[i + simplex_lp.numCol_];
    solution.row_dual[i] =
        row_dual_sign * simplex_lp.sense_ * dual[i + simplex_lp.numCol_];
  }
}

bool HighsSimplexInterface::analyseSingleHighsSolutionAndBasis(
    bool report, const HighsBasisStatus status, const double lower,
    const double upper, const double value, const double dual,
    int& num_non_basic_var, int& num_basic_var, double& off_bound_nonbasic,
    double& primal_infeasibility, double& dual_infeasibility) {
  double primal_feasibility_tolerance =
      highs_model_object.options_.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      highs_model_object.options_.dual_feasibility_tolerance;
  double middle = (lower + upper) * 0.5;

  bool query = false;
  bool count = !report;
  off_bound_nonbasic = 0;
  double primal_residual = max(lower - value, value - upper);
  primal_infeasibility = max(primal_residual, 0.);
  // ToDo Strange: nonbasic_flag seems to be inverted???
  if (status == HighsBasisStatus::BASIC) {
    // Basic variable: look for primal infeasibility
    if (count) num_basic_var++;
    if (primal_infeasibility > primal_feasibility_tolerance) {
      // Outside a bound
      if (value < lower) {
        query = true;
        if (report)
          printf(": Basic below lower bound by %12g", primal_residual);
      } else {
        query = true;
        if (report)
          printf(": Basic above upper bound by %12g", primal_residual);
      }
    }
    dual_infeasibility = fabs(dual);
    if (dual_infeasibility > dual_feasibility_tolerance) {
      query = true;
      if (report) printf(": Dual infeasibility of %12g", dual_infeasibility);
    }
  } else {
    // Nonbasic variable: look for primal and dual infeasibility
    if (count) num_non_basic_var++;

    if (primal_infeasibility > primal_feasibility_tolerance) {
      // Outside a bound
      dual_infeasibility = 0;
      if (value < lower) {
        query = true;
        if (report)
          printf(": Nonbasic below lower bound by %12g", primal_residual);
      } else {
        query = true;
        if (report)
          printf(": Nonbasic above upper bound by %12g", primal_residual);
      }
    } else if (primal_residual >= -primal_feasibility_tolerance) {
      // At a bound: check for dual feasibility
      if (lower < upper) {
        // Non-fixed variable
        if (value < middle) {
          // At lower
          dual_infeasibility = max(-dual, 0.);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            // Dual infeasiblility
            query = true;
            if (report)
              printf(": Dual infeasibility of %12g", dual_infeasibility);
          }
        } else {
          // At Upper
          dual_infeasibility = max(dual, 0.);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            // Dual infeasiblility
            query = true;
            if (report)
              printf(": Dual infeasibility of %12g", dual_infeasibility);
          }
        }
      } else {
        // Fixed variable
        dual_infeasibility = 0;
      }
    } else {
      // Between bounds (or free)
      if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
        // Free
        if (report) printf(": Nonbasic free");
      } else {
        query = true;
        if (report) printf(": Nonbasic off bound by %12g", -primal_residual);
        off_bound_nonbasic = -primal_residual;
      }
      dual_infeasibility = fabs(dual);
      if (dual_infeasibility > dual_feasibility_tolerance) {
        query = true;
        if (report) printf(": Dual infeasibility of %12g", dual_infeasibility);
      }
    }
  }
  query = false;
  return query;
}

SimplexSolutionStatus HighsSimplexInterface::analyseHighsSolutionAndBasis(
									  const int report_level,
									  const string message) {
  HighsLogMessage(HighsMessageType::INFO, "HiGHS basic solution: Analysis %s", message.c_str());
  HighsSolution& solution = highs_model_object.solution_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsLp& lp = highs_model_object.lp_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status =
      highs_model_object.simplex_lp_status_;
  double primal_feasibility_tolerance =
      highs_model_object.options_.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      highs_model_object.options_.dual_feasibility_tolerance;
  vector<double> primal_activities;
  vector<double> dual_activities;
  primal_activities.assign(lp.numRow_, 0);
  dual_activities.resize(lp.numCol_);
  int num_non_basic_var = 0;
  int num_basic_var = 0;
  int num_off_bound_nonbasic = 0;
  double max_off_bound_nonbasic = 0;
  double sum_off_bound_nonbasic = 0;
  bool header_written = false;
  double off_bound_nonbasic;
  double primal_infeasibility;
  double dual_infeasibility;
  int num_primal_infeasibilities = 0;
  double max_primal_infeasibility = 0;
  double sum_primal_infeasibilities = 0;
  int num_dual_infeasibilities = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibilities = 0;
  int num_nonzero_basic_duals = 0;
  int num_large_nonzero_basic_duals = 0;
  double max_nonzero_basic_dual = 0;
  double sum_nonzero_basic_duals = 0;
  double local_primal_objective_value = 0;
  double local_dual_objective_value = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double lower = lp.colLower_[iCol];
    double upper = lp.colUpper_[iCol];
    double value = solution.col_value[iCol];
    double dual = solution.col_dual[iCol];
    HighsBasisStatus status = basis.col_status[iCol];
    local_primal_objective_value += lp.colCost_[iCol] * value;
    if (status != HighsBasisStatus::BASIC)
      local_dual_objective_value += value * dual;
    bool report = false;
    bool query = analyseSingleHighsSolutionAndBasis(
        report, status, lower, upper, value, dual, num_non_basic_var,
        num_basic_var, off_bound_nonbasic, primal_infeasibility,
        dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic = max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
	num_nonzero_basic_duals++;
	if (abs_basic_dual > dual_feasibility_tolerance) num_large_nonzero_basic_duals++;
	max_nonzero_basic_dual = max(abs_basic_dual, max_nonzero_basic_dual);
	sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
	num_dual_infeasibilities++;
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
    report = report_level == 3 || (report_level == 2 && query);
    if (report) {
      if (!header_written) {
        printf(
            "\nColumns\nIndex NonBs Mv [          LB,           UB]       "
            "Primal         Dual    PrimalIfs      DualIfs\n");
        header_written = true;
      }
      printf("%5d %5d [%12g, %12g] %12g %12g", iCol, (int)status, lower, upper,
             value, dual);
      printf(" %12g %12g", primal_infeasibility, dual_infeasibility);
      analyseSingleHighsSolutionAndBasis(
          report, status, lower, upper, value, dual, num_non_basic_var,
          num_basic_var, off_bound_nonbasic, primal_infeasibility,
          dual_infeasibility);
      printf("\n");
    }
    dual_activities[iCol] = lp.colCost_[iCol];
    for (int el = lp.Astart_[iCol]; el < lp.Astart_[iCol + 1]; el++) {
      int iRow = lp.Aindex_[el];
      double Avalue = lp.Avalue_[el];
      primal_activities[iRow] += value * Avalue;
      dual_activities[iCol] += solution.row_dual[iRow] * Avalue;
    }
  }
  bool report = report_level >= 2;
  header_written = false;
  int num_primal_residual = 0;
  double max_primal_residual = 0;
  double sum_primal_residual = 0;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double primal_residual =
        fabs(primal_activities[iRow] - solution.row_value[iRow]);
    if (primal_residual > primal_feasibility_tolerance) {
      if (report) {
        if (!header_written) {
          printf(
              "\nRow primal residuals\nIndex     Activity     Solution     "
              "Residual\n");
          header_written = true;
        }
        printf("%5d %12g %12g %12g\n", iRow, primal_activities[iRow],
               solution.row_value[iRow], primal_residual);
      }
      num_primal_residual++;
    }
    max_primal_residual = max(primal_residual, max_primal_residual);
    sum_primal_residual += primal_residual;
  }
  header_written = false;
  int num_dual_residual = 0;
  double max_dual_residual = 0;
  double sum_dual_residual = 0;
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    double dual_residual =
        fabs(dual_activities[iCol] - solution.col_dual[iCol]);
    if (dual_residual > dual_feasibility_tolerance) {
      if (report) {
        if (!header_written) {
          printf(
              "\nRow dual residuals\nIndex     Activity     Solution     "
              "Residual\n");
          header_written = true;
        }
        printf("%5d %12g %12g %12g\n", iCol, dual_activities[iCol],
               solution.col_dual[iCol], dual_residual);
      }
      num_dual_residual++;
    }
    max_dual_residual = max(dual_residual, max_dual_residual);
    sum_dual_residual += dual_residual;
  }
  header_written = false;
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    double lower = lp.rowLower_[iRow];
    double upper = lp.rowUpper_[iRow];
    double value = solution.row_value[iRow];
    double dual = -solution.row_dual[iRow];
    HighsBasisStatus status = basis.row_status[iRow];
    if (status != HighsBasisStatus::BASIC)
      local_dual_objective_value += value * dual;
    bool report = false;
    bool query = analyseSingleHighsSolutionAndBasis(
        report, status, lower, upper, value, dual, num_non_basic_var,
        num_basic_var, off_bound_nonbasic, primal_infeasibility,
        dual_infeasibility);
    if (off_bound_nonbasic > 0) num_off_bound_nonbasic++;
    max_off_bound_nonbasic = max(off_bound_nonbasic, max_off_bound_nonbasic);
    sum_off_bound_nonbasic += off_bound_nonbasic;
    if (primal_infeasibility > primal_feasibility_tolerance)
      num_primal_infeasibilities++;
    max_primal_infeasibility =
        max(primal_infeasibility, max_primal_infeasibility);
    sum_primal_infeasibilities += primal_infeasibility;
    if (status == HighsBasisStatus::BASIC) {
      double abs_basic_dual = dual_infeasibility;
      if (abs_basic_dual > 0) {
	num_nonzero_basic_duals++;
	if (abs_basic_dual > dual_feasibility_tolerance) num_large_nonzero_basic_duals++;
	max_nonzero_basic_dual = max(abs_basic_dual, max_nonzero_basic_dual);
	sum_nonzero_basic_duals += abs_basic_dual;
      }
    } else {
      if (dual_infeasibility > dual_feasibility_tolerance)
	num_dual_infeasibilities++;
      max_dual_infeasibility = max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibilities += dual_infeasibility;
    }
    report = report_level == 3 || (report_level == 2 && query);
    if (report) {
      if (!header_written) {
        printf(
            "Rows\nIndex NonBs Mv [          LB,           UB]       Primal    "
            "     Dual    PrimalIfs      DualIfs\n");
        header_written = true;
      }
      printf("%5d %5d [%12g, %12g] %12g %12g", iRow, (int)status, lower, upper,
             value, dual);
      printf(" %12g %12g", primal_infeasibility, dual_infeasibility);
      analyseSingleHighsSolutionAndBasis(
          report, status, lower, upper, value, dual, num_non_basic_var,
          num_basic_var, off_bound_nonbasic, primal_infeasibility,
          dual_infeasibility);
      printf("\n");
    }
  }
  local_primal_objective_value += lp.offset_;
  local_dual_objective_value += lp.offset_;
  double primal_objective_value;
  double dual_objective_value;
  if (simplex_lp_status.has_primal_objective_value) {
    primal_objective_value = simplex_info.primal_objective_value;
  } else {
    primal_objective_value = local_primal_objective_value;
  }
    
  if (simplex_lp_status.has_dual_objective_value) {
    dual_objective_value = simplex_info.dual_objective_value;
  } else {
    dual_objective_value = local_dual_objective_value;
  }
  double primal_objective_error =
      fabs(primal_objective_value - local_primal_objective_value) /
      max(1.0, fabs(primal_objective_value));
  double dual_objective_error =
      fabs(dual_objective_value - local_primal_objective_value) /
      max(1.0, fabs(dual_objective_value));
  double relative_objective_difference =
      fabs(primal_objective_value - dual_objective_value) /
      max(fabs(primal_objective_value), max(fabs(dual_objective_value), 1.0));
  simplex_info.num_primal_infeasibilities = num_primal_infeasibilities;
  simplex_info.max_primal_infeasibility = max_primal_infeasibility;
  simplex_info.sum_primal_infeasibilities = sum_primal_infeasibilities;
  simplex_info.num_dual_infeasibilities = num_dual_infeasibilities;
  simplex_info.max_dual_infeasibility = max_dual_infeasibility;
  simplex_info.sum_dual_infeasibilities = sum_dual_infeasibilities;
  SimplexSolutionStatus solution_status;
  bool primal_feasible =
      num_primal_infeasibilities ==
      0;  // && max_primal_residual < primal_feasibility_tolerance;
  bool dual_feasible = num_dual_infeasibilities ==
                       0;  // && max_dual_residual < dual_feasibility_tolerance;
  if (primal_feasible) {
    if (dual_feasible) {
      solution_status = SimplexSolutionStatus::OPTIMAL;
    } else {
      solution_status = SimplexSolutionStatus::PRIMAL_FEASIBLE;
    }
  } else {
    if (dual_feasible) {
      solution_status = SimplexSolutionStatus::DUAL_FEASIBLE;
    } else {
      solution_status = SimplexSolutionStatus::UNSET;
    }
  }
  simplex_lp_status.solution_status = solution_status;
  if (num_nonzero_basic_duals) {
    HighsLogMessage(
		    HighsMessageType::WARNING,
		    "HiGHS basic solution: %d (%d large) nonzero basic duals; max = %g; sum = %g",
		    num_nonzero_basic_duals, num_large_nonzero_basic_duals, max_nonzero_basic_dual, sum_nonzero_basic_duals);
  }
  if (num_off_bound_nonbasic)  {
    HighsLogMessage(HighsMessageType::WARNING,
                    "Off-bound num/max/sum           %6d/%11.4g/%11.4g",
                    num_off_bound_nonbasic, max_off_bound_nonbasic, sum_off_bound_nonbasic);
  }
  if (report_level>0) {
    HighsLogMessage(HighsMessageType::INFO,
                    "Primal    num/max/sum residuals %6d/%11.4g/%11.4g: num/max/sum "
                    "infeasibilities %6d/%11.4g/%11.4g",
                    num_primal_residual, max_primal_residual,
                    sum_primal_residual, num_primal_infeasibilities,
                    max_primal_infeasibility, sum_primal_infeasibilities);
    HighsLogMessage(HighsMessageType::INFO,
                    "Dual      num/max/sum residuals %6d/%11.4g/%11.4g: num/max/sum "
                    "infeasibilities %6d/%11.4g/%11.4g",
                    num_dual_residual, max_dual_residual, sum_dual_residual,
                    num_dual_infeasibilities, max_dual_infeasibility,
                    sum_dual_infeasibilities);
    HighsLogMessage(HighsMessageType::INFO,
                    "Relative objective errors (primal = %.4g; dual = %.4g); "
                    "relative objective difference = %.4g",
                    primal_objective_error, dual_objective_error,
                    relative_objective_difference);
  } 
  HighsLogMessage(
      HighsMessageType::INFO,
      "HiGHS basic solution: Iterations = %d; Objective = %.15g; "
      "Infeasibilities Pr %d(%g); Du %d(%g); Status: %s",
      simplex_info.iteration_count, primal_objective_value,
      simplex_info.num_primal_infeasibilities,
      simplex_info.sum_primal_infeasibilities,
      simplex_info.num_dual_infeasibilities,
      simplex_info.sum_dual_infeasibilities,
      SimplexSolutionStatusToString(simplex_lp_status.solution_status).c_str());
#ifdef HiGHSDEV
  printf("grep_AnBsSol,%s,%s,%.15g,%s,%d,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g,%d,%g,%g\n",
	 lp.model_name_.c_str(), message.c_str(), primal_objective_value,
	 SimplexSolutionStatusToString(simplex_lp_status.solution_status).c_str(),
	 num_nonzero_basic_duals, num_large_nonzero_basic_duals, max_nonzero_basic_dual, sum_nonzero_basic_duals,
	 num_off_bound_nonbasic, max_off_bound_nonbasic, sum_off_bound_nonbasic,
	 num_primal_residual, max_primal_residual, sum_primal_residual,
	 num_primal_infeasibilities, max_primal_infeasibility, sum_primal_infeasibilities,
	 num_dual_residual, max_dual_residual, sum_dual_residual,
	 num_dual_infeasibilities, max_dual_infeasibility, sum_dual_infeasibilities);
#endif
  return solution_status;
}

int HighsSimplexInterface::get_basic_indices(int* bind) {
  SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = simplex_basis.basicIndex_[row];
    if (var >= simplex_lp.numCol_)
      bind[row] = -(1 + var - simplex_lp.numCol_);
    else
      bind[row] = var;
  }
  return 0;
}

#ifdef HiGHSDEV
void HighsSimplexInterface::check_load_from_postsolve() {
  HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  bool ok;

  ok = nonbasic_flag_basic_index_ok(simplex_lp,
                                    highs_model_object.simplex_basis_);
  printf(
      "check_load_from_postsolve: return from nonbasic_flag_basic_index_ok = "
      "%d\n",
      ok);
  assert(ok);

  ok = all_nonbasic_move_vs_work_arrays_ok(highs_model_object);
  printf(
      "check_load_from_postsolve: return from "
      "all_nonbasic_move_vs_work_arrays_ok = %d\n",
      ok);
  assert(ok);
}
#endif
