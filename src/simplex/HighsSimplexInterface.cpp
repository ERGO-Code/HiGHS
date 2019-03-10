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
#include "HConfig.h"
#include "lp_data/HighsLpUtils.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/HSimplex.h"
#include "io/HighsIO.h"
#include "io/HMPSIO.h"
#include "util/HighsUtils.h"

HighsStatus HighsSimplexInterface::util_add_cols(int XnumNewCol, const double *XcolCost, const double *XcolLower,  const double *XcolUpper,
						 int XnumNewNZ, const int *XAstart, const int *XAindex, const double *XAvalue) {
#ifdef HiGHSDEV
  printf("Called util_add_cols(XnumNewCol=%d, XnumNewNZ = %d)\n", XnumNewCol, XnumNewNZ);
#endif
  HighsStatus return_status = HighsStatus::NotSet;
  if (XnumNewCol < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewCol == 0) return HighsStatus::OK;
  if (XnumNewCol > 0) if (isColDataNull(XcolCost, XcolLower,  XcolUpper)) return HighsStatus::Error;
  if (XnumNewNZ > 0) if (isMatrixDataNull(XAstart, XAindex, XAvalue)) return HighsStatus::Error;

  HighsLp &lp = highs_model_object.lp_;
  HighsOptions &options = highs_model_object.options_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &simplex_basis = highs_model_object.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_lp_matrix = true;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_basis = simplex_lp_status.has_basis;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
  bool apply_row_scaling = simplex_lp_status.is_scaled;

  // Check that if nonzeros are to be added then the model has a positive number of rows
  if (lp.numRow_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numRow_ <= 0 && XnumNewNZ > 0)) return HighsStatus::Error;

  // Record the new number of columns
  int newNumCol = lp.numCol_ + XnumNewCol;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or scaling
  if (!valid_simplex_lp) {
    assert(!simplex_basis.valid_);
    assert(!valid_simplex_matrix);
    assert(!apply_row_scaling);
  }
#endif
  HighsStatus call_status;
  call_status = append_lp_cols(lp, XnumNewCol, XcolCost, XcolLower, XcolUpper,
			       XnumNewNZ, XAstart, XAindex, XAvalue,
			       options, valid_lp_matrix);
  return_status = worseStatus(call_status, return_status);
  if (return_status == HighsStatus::Error) return return_status;

  if (valid_simplex_lp) {
    call_status = append_lp_cols(simplex_lp, XnumNewCol, XcolCost, XcolLower, XcolUpper,
				 XnumNewNZ, XAstart, XAindex, XAvalue,
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
  if (valid_simplex_basis) append_nonbasic_cols_to_basis(simplex_lp, simplex_basis, newNumCol);

  // Deduce the consequences of adding new columns
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_COLS);

  // Increase the number of columns in the LPs
  lp.numCol_ += XnumNewCol;
  if (valid_simplex_lp) simplex_lp.numCol_ += XnumNewCol;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basisOK = nonbasic_flag_basic_index_ok(lp, basis);
    assert(basisOK);
    report_basis(lp, basis);
  }
  if (valid_simplex_basis) {
    bool simplex_basisOK = nonbasic_flag_basic_index_ok(simplex_lp, simplex_basis);
    assert(simplex_basisOK);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  return return_status;

}

HighsStatus HighsSimplexInterface::delete_cols(int from_col, int to_col) {
  return delete_cols_general(
			     true, from_col, to_col,
			     false, 0, NULL,
			     false, NULL
			     );
}

HighsStatus HighsSimplexInterface::delete_cols(int num_set_entries, const int* col_set) {
  return delete_cols_general(
			     false, 0, 0,
			     true, num_set_entries, col_set,
			     false, NULL
			     );
}

HighsStatus HighsSimplexInterface::delete_cols(int* col_mask) {
  return delete_cols_general(
			     false, 0, 0,
			     false, 0, NULL,
			     true, col_mask
			     );
}

HighsStatus HighsSimplexInterface::delete_cols_general(bool interval, int from_col, int to_col,
						       bool set, int num_set_entries, const int* col_set,
						       bool mask, int* col_mask) {
  // Uses to_col in iterator style
  HighsLp &lp = highs_model_object.lp_;

  HighsBasis &basis = highs_model_object.basis_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &simplex_basis = highs_model_object.simplex_basis_;

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
  returnStatus = delete_lp_cols(lp, 
				interval, from_col, to_col,
				set, num_set_entries, col_set,
				mask, col_mask,
				valid_matrix);
  if (returnStatus != HighsStatus::OK) return returnStatus;
  // ToDo Determine consequences for basis when deleting columns
  basis.valid_ = false;
  
  if (valid_simplex_lp) {
    returnStatus = delete_lp_cols(simplex_lp, 
				  interval, from_col, to_col,
				  set, num_set_entries, col_set,
				  mask, col_mask,
				  valid_simplex_matrix);
    if (returnStatus != HighsStatus::OK) return returnStatus;
    //    for (int col = from_col; col < lp.numCol_ - numDeleteCol; col++) scale.col_[col] = scale.col_[col + numDeleteCol];
    // ToDo Determine consequences for basis when deleting columns
    simplex_lp_status.has_matrix_col_wise = false;
    simplex_lp_status.has_matrix_row_wise = false;
    simplex_basis.valid_ = false;
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::util_add_rows(int XnumNewRow, const double *XrowLower, const double *XrowUpper,
						 int XnumNewNZ, const int *XARstart, const int *XARindex, const double *XARvalue) {
#ifdef HiGHSDEV
  printf("Called util_add_rows(XnumNewRow=%d, XnumNewNZ = %d)\n", XnumNewRow, XnumNewNZ);
#endif
  HighsStatus return_status = HighsStatus::NotSet;
  if (XnumNewRow < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewRow == 0) return HighsStatus::OK;
  if (XnumNewRow > 0) if (isRowDataNull(XrowLower, XrowUpper)) return HighsStatus::Error;
  if (XnumNewNZ > 0) if (isMatrixDataNull(XARstart, XARindex, XARvalue)) return HighsStatus::Error;

  HighsLp &lp = highs_model_object.lp_;
  HighsOptions &options = highs_model_object.options_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &simplex_basis = highs_model_object.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool valid_basis = basis.valid_;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_basis = simplex_lp_status.has_basis;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
  bool apply_row_scaling = simplex_lp_status.is_scaled;

  // Check that if nonzeros are to be added then the model has a positive number of columns
  if (lp.numCol_ <= 0 && XnumNewNZ > 0) return HighsStatus::Error;
  if (valid_simplex_lp && (simplex_lp.numCol_ <= 0 && XnumNewNZ > 0)) return HighsStatus::Error;

  // Record the new number of rows
  int newNumRow = lp.numRow_ + XnumNewRow;

#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no basis, matrix or scaling
  if (!valid_simplex_lp) {
    assert(!simplex_basis.valid_);
    assert(!valid_simplex_matrix);
    assert(!apply_row_scaling);
  }
#endif
  // Assess the bounds and matrix indices, returning on error
  bool normalise = false;
  HighsStatus call_status;
  call_status = assess_bounds("Row", lp.numRow_, XnumNewRow, true, 0, XnumNewRow, false, 0, NULL, false, NULL,
			      (double*)XrowLower, (double*)XrowUpper, options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);
  
  if (XnumNewNZ) {
    call_status = assessMatrix(lp.numCol_, 0, XnumNewRow, XnumNewRow,
			       XnumNewNZ, (int*)XARstart, (int*)XARindex, (double*)XARvalue,
			       options.small_matrix_value, options.large_matrix_value, normalise);
    return_status = worseStatus(call_status, return_status);
    if (return_status == HighsStatus::Error) return return_status;
  }

  // Append the columns to the LP vectors and matrix
  append_rows_to_lp_vectors(lp, XnumNewRow, XrowLower, XrowUpper);

  // Normalise the LP row bounds
  normalise = true;
  call_status = assess_bounds("Row", lp.numRow_, newNumRow, true, 0, newNumRow, false, 0, NULL, false, NULL,
			     &lp.rowLower_[0], &lp.rowUpper_[0], options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);

  int lc_XnumNewNZ = XnumNewNZ;
  int* lc_XARstart;
  int* lc_XARindex;
  double* lc_XARvalue;
  if (XnumNewNZ) {
    // Copy the new row-wise matrix into a local copy that can be normalised
    std::memcpy(lc_XARstart, XARstart, sizeof(int)*XnumNewRow);
    std::memcpy(lc_XARindex, XARindex, sizeof(int)*XnumNewNZ);
    std::memcpy(lc_XARvalue, XARvalue, sizeof(double)*XnumNewNZ);
    // Normalise the new matrix columns
    normalise = true;
    call_status = assessMatrix(lp.numCol_, 0, XnumNewRow, XnumNewRow,
			       lc_XnumNewNZ, lc_XARstart, lc_XARindex, lc_XARvalue,
			       options.small_matrix_value, options.large_matrix_value, normalise);
    if (lc_XnumNewNZ) {
      // Append rows to LP matrix
      append_rows_to_lp_matrix(lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart, lc_XARindex, lc_XARvalue);
    }
  }

  if (valid_simplex_lp) {
    append_rows_to_lp_vectors(simplex_lp, XnumNewRow, XrowLower, XrowUpper);
    call_status = assess_bounds("Row", simplex_lp.numRow_, newNumRow, true, 0, newNumRow, false, 0, NULL, false, NULL,
				&simplex_lp.colLower_[0], &simplex_lp.colUpper_[0], options.infinite_bound, normalise);
    return_status = worseStatus(call_status, return_status);
  }
  if (valid_simplex_matrix && lc_XnumNewNZ) {
    append_rows_to_lp_matrix(simplex_lp, XnumNewRow, lc_XnumNewNZ, lc_XARstart, lc_XARindex, lc_XARvalue);
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

  // Update the basis correponding to new nonbasic rows
  if (valid_basis) append_basic_rows_to_basis(lp, basis, newNumRow);
  if (valid_simplex_basis) append_basic_rows_to_basis(simplex_lp, simplex_basis, newNumRow);

  // Deduce the consequences of adding new rows
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_ROWS);

  // Increase the number of rows in the LPs
  lp.numRow_ += XnumNewRow;
  if (valid_simplex_lp) simplex_lp.numRow_ += XnumNewRow;

#ifdef HiGHSDEV
  if (valid_basis) {
    bool basisOK = nonbasic_flag_basic_index_ok(lp, basis);
    assert(basisOK);
    report_basis(lp, basis);
  }
  if (valid_simplex_basis) {
    bool simplex_basisOK = nonbasic_flag_basic_index_ok(simplex_lp, simplex_basis);
    assert(simplex_basisOK);
    report_basis(simplex_lp, simplex_basis);
  }
#endif
  return return_status;

}

HighsStatus HighsSimplexInterface::delete_rows(int from_row, int to_row) {
  return delete_rows_general(
			     true, from_row, to_row,
			     false, 0, NULL,
			     false, NULL
			     );
}

HighsStatus HighsSimplexInterface::delete_rows(int num_set_entries, const int* row_set) {
  return delete_rows_general(
			     false, 0, 0,
			     true, num_set_entries, row_set,
			     false, NULL
			     );
}

HighsStatus HighsSimplexInterface::delete_rows(int* row_mask) {
  return delete_rows_general(
			     false, 0, 0,
			     false, 0, NULL,
			     true, row_mask
			     );
}

HighsStatus HighsSimplexInterface::delete_rows_general(bool interval, int from_row, int to_row,
						       bool set, int num_set_entries, const int* row_set,
						       bool mask, int* row_mask) {
  // Uses to_row in iterator style
#ifdef HiGHSDEV
  printf("Called model.util_deleteRows(from_row=%d, to_row=%d)\n", from_row, to_row);
#endif
  HighsLp &lp = highs_model_object.lp_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsScale &scale = highs_model_object.scale_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &simplex_basis = highs_model_object.simplex_basis_;

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
  returnStatus = delete_lp_rows(lp, 
				interval, from_row, to_row,
				set, num_set_entries, row_set,
				mask, row_mask,
				valid_matrix);
  if (returnStatus != HighsStatus::OK) return returnStatus;
  // ToDo Determine consequences for basis when deleting rowumns
  basis.valid_ = false;
  
  if (valid_simplex_lp) {
    returnStatus = delete_lp_rows(simplex_lp, 
				  interval, from_row, to_row,
				  set, num_set_entries, row_set,
				  mask, row_mask,
				  valid_simplex_matrix);
    if (returnStatus != HighsStatus::OK) return returnStatus;
    //    for (int row = from_row; row < lp.numRow_ - numDeleteRow; row++) scale.row_[row] = scale.row_[row + numDeleteRow];
    // ToDo Determine consequences for basis when deleting rowumns
    simplex_lp_status.has_matrix_col_wise = false;
    simplex_lp_status.has_matrix_row_wise = false;
    simplex_basis.valid_ = false;
  }
    // ToDo Determine consequences for basis when deleting rows
  //  update_simplex_lp_status(simplex_lp_status, LpAction::DEL_ROWS);
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::getCols(const int from_col, const int to_col,
					   int &num_col, double *col_cost, double *col_lower, double *col_upper,
					   int &num_nz, int *col_matrix_start, int *col_matrix_index, double *col_matrix_value) {
  return getColsGeneral(
			true, from_col, to_col,
			false, 0, NULL,
			false, NULL,
			num_col, col_cost, col_lower, col_upper,
			num_nz, col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getCols(const int num_set_entries, const int* col_set,
					   int &num_col, double *col_cost, double *col_lower, double *col_upper,
					   int &num_nz, int *col_matrix_start, int *col_matrix_index, double *col_matrix_value) {
  return getColsGeneral(
			false, 0, 0,
			true, num_set_entries, col_set,
			false, NULL,
			num_col, col_cost, col_lower, col_upper,
			num_nz, col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getCols(const int* col_mask,
					   int &num_col, double *col_cost, double *col_lower, double *col_upper,
					   int &num_nz, int *col_matrix_start, int *col_matrix_index, double *col_matrix_value) {
  return getColsGeneral(
			false, 0, 0,
			false, 0, NULL,
			true, col_mask,
			num_col, col_cost, col_lower, col_upper,
			num_nz, col_matrix_start, col_matrix_index, col_matrix_value);
}

HighsStatus HighsSimplexInterface::getColsGeneral(const bool interval, const int from_col, const int to_col,
						  const bool set, const int num_set_entries, const int* col_set, 
						  const bool mask, const int* col_mask,
						  int &num_col, double *col_cost, double *col_lower, double *col_upper,
						  int &num_nz, int *col_matrix_start, int *col_matrix_index, double *col_matrix_value) {
  int from_k;
  int to_k;
  HighsLp &lp = highs_model_object.lp_;
  HighsStatus return_status = assess_interval_set_mask(lp.numCol_,
						     interval, from_col, to_col,
						     set, num_set_entries, col_set,
						     mask, col_mask,
						     from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  if (from_k < 0 || to_k > lp.numCol_) return HighsStatus::Error;
  if (from_k >= to_k) return HighsStatus::OK;
  int out_from_col;
  int out_to_col;
  int in_from_col;
  int in_to_col = 0;
  int current_set_entry = 0;
  int col_dim = lp.numCol_;
  num_col = 0;
  num_nz = 0;
  for (int k = from_k; k < to_k; k++) {
    update_out_in_ix(col_dim,
		     interval, from_col, to_col,
		     set, num_set_entries, col_set,
		     mask, col_mask,
		     out_from_col, out_to_col,
		     in_from_col, in_to_col,
		     current_set_entry);
    assert(out_to_col <= col_dim);
    assert(in_to_col <= col_dim);
    for (int col = out_from_col; col < out_to_col; col++) {
      col_cost[num_col] = lp.colCost_[col];
      col_lower[num_col] = lp.colLower_[col];
      col_upper[num_col] = lp.colUpper_[col];
      col_matrix_start[num_col] = num_nz + lp.Astart_[col] - lp.Astart_[out_from_col];
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
					   int &num_row, double *row_lower, double *row_upper,
					   int &num_nz, int *row_matrix_start, int *row_matrix_index, double *row_matrix_value) {
  return getRowsGeneral(
			true, from_row, to_row,
			false, 0, NULL,
			false, NULL,
			num_row, row_lower, row_upper,
			num_nz, row_matrix_start, row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRows(const int num_set_entries, const int* row_set,
					   int &num_row, double *row_lower, double *row_upper,
					   int &num_nz, int *row_matrix_start, int *row_matrix_index, double *row_matrix_value) {
  return getRowsGeneral(
			false, 0, 0,
			true, num_set_entries, row_set,
			false, NULL,
			num_row, row_lower, row_upper,
			num_nz, row_matrix_start, row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRows(const int* row_mask,
					   int &num_row, double *row_lower, double *row_upper,
					   int &num_nz, int *row_matrix_start, int *row_matrix_index, double *row_matrix_value) {
  return getRowsGeneral(
			false, 0, 0,
			false, 0, NULL,
			true, row_mask,
			num_row, row_lower, row_upper,
			num_nz, row_matrix_start, row_matrix_index, row_matrix_value);
}

HighsStatus HighsSimplexInterface::getRowsGeneral(const bool interval, const int from_row, const int to_row,
						  const bool set, const int num_set_entries, const int* row_set, 
						  const bool mask, const int* row_mask,
						  int &num_row, double *row_lower, double *row_upper,
						  int &num_nz, int *row_matrix_start, int *row_matrix_index, double *row_matrix_value) {
  /*
  HighsLp &lp = highs_model_object.lp_;
  if (XfromRow < 0 || XtoRow > lp.numRow_) return HighsStatus::Error;
  if (XfromRow >= XtoRow) return HighsStatus::OK;

  // Determine the number of rows to be extracted
  int numExtractRows = XtoRow - XfromRow;
  for (int row = XfromRow; row < XtoRow; row++) {
    XrowLower[row - XfromRow] = lp.rowLower_[row];
    XrowUpper[row - XfromRow] = lp.rowUpper_[row];
  }
  // Determine how many entries are in each row to be extracted
  vector<int> XARlength;
  XARlength.assign(numExtractRows, 0);

  for (int el = lp.Astart_[0]; el < lp.Astart_[lp.numCol_]; el++) {
    int row = lp.Aindex_[el];
    if (row >= XfromRow && row < XtoRow) XARlength[row - XfromRow] += 1;
  }
  XARstart[0] = 0;
  for (int row = 0; row < numExtractRows-1; row++) {
    XARstart[row + 1] = XARstart[row] + XARlength[row];
    XARlength[row] = 0;
  }
  XARlength[numExtractRows-1] = 0;

  for (int col = 0; col < lp.numCol_; col++) {
    for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
      int row = lp.Aindex_[el];
      if (row >= XfromRow && row < XtoRow) {
        int rowEl = XARstart[row - XfromRow] + XARlength[row - XfromRow];
        XARlength[row - XfromRow] += 1;
        XARindex[rowEl] = col;
        XARvalue[rowEl] = lp.Avalue_[el];
      }
    }
  }
  *XnumNZ = XARstart[numExtractRows-1] + XARlength[numExtractRows-1];
  //  printf("Set XnumNZ = %d\n", *XnumNZ);
*/
}

// Change a single coefficient in the matrix
HighsStatus HighsSimplexInterface::util_change_coefficient(int Xrow, int Xcol, const double XnewValue) {
#ifdef HiGHSDEV
  printf("Called util_changeCoeff(Xrow=%d, Xcol=%d, XnewValue=%g)\n", Xrow, Xcol, XnewValue);
#endif
  HighsLp &lp = highs_model_object.lp_;
  if (Xrow < 0 || Xrow > lp.numRow_) return HighsStatus::Error;
  if (Xcol < 0 || Xcol > lp.numCol_) return HighsStatus::Error;
  //  printf("\n\nCalled model.util_changeCoeff(row=%d, col=%d, newval=%g)\n\n", Xrow, Xcol, XnewValue);
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  bool valid_simplex_lp = simplex_lp_status.valid;
  bool valid_simplex_matrix = simplex_lp_status.has_matrix_col_wise;
#ifdef HiGHSDEV
  // Check that if there is no simplex LP then there is no matrix or scaling
  if (!valid_simplex_lp) {
    assert(!valid_simplex_matrix);
    //    assert(!apply_row_scaling);
  }
#endif
  change_lp_matrix_coefficient(lp, Xrow, Xcol, XnewValue);
  if (valid_simplex_lp) {
    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
    HighsScale &scale = highs_model_object.scale_;
    double scaledXnewValue = XnewValue*scale.row_[Xrow]*scale.col_[Xcol];
    change_lp_matrix_coefficient(simplex_lp, Xrow, Xcol, scaledXnewValue);
  }
  // simplex_lp.reportLp();
  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if
   update_simplex_lp_status(simplex_lp_status, LpAction::NEW_ROWS);
  //  simplex_lp.reportLp();
}

void HighsSimplexInterface::shift_objective_value(double Xshift) {
  printf("Where is shift_objective_value required - so I can interpret what's required\n");
  // Update the LP objective value with the shift
  highs_model_object.simplex_info_.dualObjectiveValue += Xshift;
  // Update the LP offset with the shift
  HighsLp &lp = highs_model_object.lp_;
  highs_model_object.lp_.offset_ += Xshift;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    // Update the simplex LP offset with the shift
    highs_model_object.simplex_lp_.offset_ += Xshift;
  }
}

HighsStatus HighsSimplexInterface::change_ObjSense(int Xsense){
  HighsLp &lp = highs_model_object.lp_;
  if ((Xsense == OBJSENSE_MINIMIZE) != (lp.sense_ == OBJSENSE_MINIMIZE)) {
    // Flip the LP objective sense
    lp.sense_ = Xsense;
  }
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  if (simplex_lp_status.valid) {
    HighsLp &simplex_lp = highs_model_object.simplex_lp_;
    HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
    if ((Xsense == OBJSENSE_MINIMIZE) != (simplex_lp.sense_ == OBJSENSE_MINIMIZE)) {
      // Flip the objective sense
      simplex_lp.sense_ = Xsense;
      simplex_lp_status.solution_status = SimplexSolutionStatus::UNSET;
    }
  }
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::change_costs(int from_col, int to_col, const double* usr_col_cost) {
  return change_costs_general(
			      true, from_col, to_col,
			      false, 0, NULL,
			      false, NULL,
			      usr_col_cost);
}

HighsStatus HighsSimplexInterface::change_costs(int num_set_entries, const int* col_set, const double* usr_col_cost) {
  return change_costs_general(
			      false, 0, 0,
			      true, num_set_entries, col_set,
			      false, NULL,
			      usr_col_cost);
}

HighsStatus HighsSimplexInterface::change_costs(const int* col_mask, const double* usr_col_cost) {
  return change_costs_general(
			      false, 0, 0,
			      false, 0, NULL,
			      true, col_mask,
			      usr_col_cost);
}

HighsStatus HighsSimplexInterface::change_costs_general(
							bool interval, int from_col, int to_col,
							bool set, int num_set_entries, const int* col_set,
							bool mask, const int* col_mask,
							const double* usr_col_cost) {
  bool null_data = false;
  if (usr_col_cost == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied column costs are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  HighsStatus call_status = change_lp_costs(highs_model_object.lp_, 
					    interval, from_col, to_col,
					    set, num_set_entries, col_set,
					    mask, col_mask,
					    usr_col_cost, highs_model_object.options_.infinite_cost);
 if (call_status == HighsStatus::Error) return HighsStatus::Error;
 // Deduce the consequences of new costs
 update_simplex_lp_status(highs_model_object.simplex_lp_status_, LpAction::NEW_COSTS);
 return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::change_col_bounds(int from_col, int to_col, const double* usr_col_lower, const double* usr_col_upper) {
  return change_col_bounds_general(
			      true, from_col, to_col,
			      false, 0, NULL,
			      false, NULL,
			      usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::change_col_bounds(int num_set_entries, const int* col_set, const double* usr_col_lower, const double* usr_col_upper) {
  return change_col_bounds_general(
			      false, 0, 0,
			      true, num_set_entries, col_set,
			      false, NULL,
			      usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::change_col_bounds(const int* col_mask, const double* usr_col_lower, const double* usr_col_upper) {
  return change_col_bounds_general(
			      false, 0, 0,
			      false, 0, NULL,
			      true, col_mask,
			      usr_col_lower, usr_col_upper);
}

HighsStatus HighsSimplexInterface::change_col_bounds_general(
							bool interval, int from_col, int to_col,
							bool set, int num_set_entries, const int* col_set,
							bool mask, const int* col_mask,
							const double* usr_col_lower, const double* usr_col_upper) {
  bool null_data = false;
  if (usr_col_lower == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied column lower bounds are NULL");
    null_data = true;
  }
  if (usr_col_upper == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied column upper bounds are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  HighsStatus call_status = change_lp_col_bounds(highs_model_object.lp_, 
						 interval, from_col, to_col,
						 set, num_set_entries, col_set,
						 mask, col_mask,
						 usr_col_lower, usr_col_upper, highs_model_object.options_.infinite_bound);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;
  // Deduce the consequences of new bounds
  update_simplex_lp_status(highs_model_object.simplex_lp_status_, LpAction::NEW_BOUNDS);
  return HighsStatus::OK;
}

HighsStatus HighsSimplexInterface::change_row_bounds(int from_row, int to_row, const double* usr_row_lower, const double* usr_row_upper) {
  return change_row_bounds_general(
			      true, from_row, to_row,
			      false, 0, NULL,
			      false, NULL,
			      usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::change_row_bounds(int num_set_entries, const int* row_set, const double* usr_row_lower, const double* usr_row_upper) {
  return change_row_bounds_general(
			      false, 0, 0,
			      true, num_set_entries, row_set,
			      false, NULL,
			      usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::change_row_bounds(const int* row_mask, const double* usr_row_lower, const double* usr_row_upper) {
  return change_row_bounds_general(
			      false, 0, 0,
			      false, 0, NULL,
			      true, row_mask,
			      usr_row_lower, usr_row_upper);
}

HighsStatus HighsSimplexInterface::change_row_bounds_general(
							bool interval, int from_row, int to_row,
							bool set, int num_set_entries, const int* row_set,
							bool mask, const int* row_mask,
							const double* usr_row_lower, const double* usr_row_upper) {
  bool null_data = false;
  if (usr_row_lower == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied row lower bounds are NULL");
    null_data = true;
  }
  if (usr_row_upper == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied row upper bounds are NULL");
    null_data = true;
  }
  if (null_data) return HighsStatus::Error;
  HighsStatus call_status = change_lp_row_bounds(highs_model_object.lp_, 
						 interval, from_row, to_row,
						 set, num_set_entries, row_set,
						 mask, row_mask,
						 usr_row_lower, usr_row_upper, highs_model_object.options_.infinite_bound);
  if (call_status == HighsStatus::Error) return HighsStatus::Error;
  // Deduce the consequences of new bounds
  update_simplex_lp_status(highs_model_object.simplex_lp_status_, LpAction::NEW_BOUNDS);
  return HighsStatus::OK;
}


#ifdef HiGHSDEV
void HighsSimplexInterface::change_update_method(int updateMethod) {
  highs_model_object.factor_.change(updateMethod);
}
#endif

void HighsSimplexInterface::report_simplex_outcome(const char *message) {
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsTimer &timer = highs_model_object.timer_;

  string ch7_status;
  if (simplex_lp_status.solution_status == SimplexSolutionStatus::OPTIMAL)
    ch7_status = "OPTIMAL";
  else
    ch7_status = "NOT-OPT";
  HighsPrintMessage(ML_ALWAYS, "%s: %7s\n", message, ch7_status.c_str());


  double dualObjectiveValue = simplex_info.dualObjectiveValue;
  double currentRunHighsTime = timer.readRunHighsClock();
#ifdef SCIP_DEV
  double prObjVal = compute_primal_objective_function_value(highs_model_object);
  double dlObjVal = abs(prObjVal - dualObjectiveValue) / max(abs(dualObjectiveValue), max(abs(prObjVal), 1.0));
  HighsPrintMessage(ML_MINIMAL, "%32s: PrObj=%20.10e; DuObj=%20.10e; DlObj=%g; Iter=%10d; %10.3f",
		    simplex_lp.model_name_.c_str(),
		    prObjVal,
		    dualObjectiveValue,
		    dlObjVal,
		    simplex_info.iteration_count,
		    currentRunHighsTime);
#endif
  HighsLogMessage(HighsMessageType::INFO, "Model name: %-s", simplex_lp.model_name_.c_str());
  HighsLogMessage(HighsMessageType::INFO, "Run status: %-s", SimplexSolutionStatusToString(simplex_lp_status.solution_status).c_str());
  HighsLogMessage(HighsMessageType::INFO, "Objective value:             %15.8g", dualObjectiveValue);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (dual phase 1): %12d", simplex_info.dual_phase1_iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (dual phase 2): %12d", simplex_info.dual_phase2_iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (primal):       %12d", simplex_info.primal_phase2_iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Iteration count (total):        %12d", simplex_info.iteration_count);
  HighsLogMessage(HighsMessageType::INFO, "Run time:                       %12.2f", currentRunHighsTime);

  // Greppable report line added
  HighsLogMessage(HighsMessageType::INFO, "grep_HiGHS,%15.8g,%d,%g,Status,%d,%16s,%d,%d,%d\n",
		    dualObjectiveValue,
		    simplex_info.iteration_count,
		    currentRunHighsTime,
		    simplex_lp_status.solution_status,
		    simplex_lp.model_name_.c_str(),
		    simplex_info.dual_phase1_iteration_count,
		    simplex_info.dual_phase2_iteration_count,
		    simplex_info.primal_phase2_iteration_count
		    );
}

double HighsSimplexInterface::get_lp_objective_value(vector<double> &XcolValue) {
  HighsLp &lp = highs_model_object.lp_;

  double lp_objective_value = 0;
  for (int i = 0; i < lp.numCol_; i++) lp_objective_value += XcolValue[i] * lp.colCost_[i];
  return lp_objective_value;
}

void HighsSimplexInterface::get_primal_dual_values(vector<double> &XcolValue,
						   vector<double> &XcolDual,
						   vector<double> &XrowValue,
						   vector<double> &XrowDual
						   ) {
  HighsScale &scale = highs_model_object.scale_;
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  // Take primal solution
  vector<double> value = simplex_info.workValue_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++)
    value[basis.basicIndex_[iRow]] = simplex_info.baseValue_[iRow];
  // Take dual solution
  vector<double> dual = simplex_info.workDual_;
  for (int iRow = 0; iRow < simplex_lp.numRow_; iRow++) dual[basis.basicIndex_[iRow]] = 0;
  // Scale back
  for (int iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    value[iCol] *= scale.col_[iCol];
    dual[iCol] /= (scale.col_[iCol] / scale.cost_);
  }
  for (int iRow = 0, iTot = simplex_lp.numCol_; iRow < simplex_lp.numRow_; iRow++, iTot++) {
    value[iTot] /= scale.row_[iRow];
    dual[iTot] *= (scale.row_[iRow] * scale.cost_);
  }

  //************** part 2: gepr and gedu
  // Now we can get the solution
  XcolValue.resize(simplex_lp.numCol_);
  XcolDual.resize(simplex_lp.numCol_);
  XrowValue.resize(simplex_lp.numRow_);
  XrowDual.resize(simplex_lp.numRow_);

  double *valuePtr = &value[0];
  for (int i = 0; i < simplex_lp.numRow_; i++) XrowValue[i] = -valuePtr[i + simplex_lp.numCol_];
  for (int i = 0; i < simplex_lp.numCol_; i++) XcolValue[i] = valuePtr[i];
  for (int i = 0; i < simplex_lp.numRow_; i++) XrowDual[i] = simplex_lp.sense_ * dual[i + simplex_lp.numCol_];
  for (int i = 0; i < simplex_lp.numCol_; i++) XcolDual[i] = simplex_lp.sense_ * dual[i];
}

void HighsSimplexInterface::get_basicIndex_nonbasicFlag(vector<int> &XbasicIndex, vector<int> &XnonbasicFlag) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  XbasicIndex.resize(simplex_lp.numRow_);
  XnonbasicFlag.resize(basis.nonbasicFlag_.size());
  int basicIndexSz = basis.basicIndex_.size();
  for (int i = 0; i < basicIndexSz; i++) XbasicIndex[i] = basis.basicIndex_[i];
  int nonbasicFlagSz = basis.nonbasicFlag_.size();
  for (int i = 0; i < nonbasicFlagSz; i++) XnonbasicFlag[i] = basis.nonbasicFlag_[i];
}

int HighsSimplexInterface::get_basic_indices(int *bind) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = basis.basicIndex_[row];
    if (var >= simplex_lp.numCol_)
      bind[row] = -(1 + var - simplex_lp.numCol_);
    else
      bind[row] = var;
  }
  return 0;
}

  // Utilities to convert model basic/nonbasic status to/from SCIP-like status
int HighsSimplexInterface::convert_baseStat_to_working(const int* cstat, const int* rstat) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  //  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;

  int numBasic = 0;
  for (int col = 0; col < simplex_lp.numCol_; col++) {
    int var = col;
    if (cstat[col] == (int) HighsBasisStatus::BASIC) {
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (cstat[col] == (int) HighsBasisStatus::LOWER) {
      // (int) HighsBasisStatus::LOWER includes fixed variables
      if (simplex_lp.colLower_[col] == simplex_lp.colUpper_[col]) {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
        continue;
      }
    } else if (cstat[col] == (int) HighsBasisStatus::UPPER) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
      continue;
    } else if (cstat[col] == (int) HighsBasisStatus::ZERO) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: col=%d, cstat=%d, lower=%g, upper=%g\n",
             col, cstat[col], simplex_lp.colLower_[col], simplex_lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  for (int row = 0; row < simplex_lp.numRow_; row++) {
    int var = simplex_lp.numCol_ + row;
    if (rstat[row] == (int) HighsBasisStatus::BASIC) {
      basis.nonbasicFlag_[var] = NONBASIC_FLAG_FALSE;
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      basis.basicIndex_[numBasic] = var;
      numBasic++;
      continue;
    }
    basis.nonbasicFlag_[var] = NONBASIC_FLAG_TRUE;
    if (rstat[row] == (int) HighsBasisStatus::LOWER) {
      // (int) HighsBasisStatus::LOWER includes fixed variables
      if (simplex_lp.rowLower_[row] == simplex_lp.rowUpper_[row]) {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
        continue;
      } else {
        basis.nonbasicMove_[var] = NONBASIC_MOVE_DN;
        continue;
      }
    } else if (rstat[row] == (int) HighsBasisStatus::UPPER) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_UP;
      continue;
    } else if (rstat[row] == (int) HighsBasisStatus::ZERO) {
      basis.nonbasicMove_[var] = NONBASIC_MOVE_ZE;
      continue;
    } else {
#ifdef HiGHSDEV
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row]);
#endif
      printf("Invalid basis status: row=%d, rstat=%d, lower=%g, upper=%g\n",
             row, rstat[row], simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row]);
      return -(row + 1);
    }
    printf(
        "convertBaseStatToWorking: row=%d, rstat=%d, lower=%g, upper=%g, "
        "nonbasicMove=%d\n",
        row, rstat[row], simplex_lp.rowLower_[row], simplex_lp.rowUpper_[row], basis.nonbasicMove_[var]);
  }
  assert(numBasic = simplex_lp.numRow_);
  populate_work_arrays(highs_model_object);
  update_simplex_lp_status(simplex_lp_status, LpAction::NEW_BASIS);
  return 0;
}

int HighsSimplexInterface::convert_Working_to_BaseStat(int* cstat, int* rstat) {
  HighsBasis &basis = highs_model_object.basis_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  if (cstat != NULL) {
    for (int col = 0; col < simplex_lp.numCol_; col++) {
      int var = col;
      if (!basis.nonbasicFlag_[var]) {
        cstat[col] = (int) HighsBasisStatus::BASIC;
        continue;
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(-simplex_lp.colLower_[col]))
#endif
        {
          cstat[col] = (int) HighsBasisStatus::LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
#ifdef HiGHSDEV
        if (!highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
        {
          cstat[col] = (int) HighsBasisStatus::UPPER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        //	printf("Var %d Move = %d [%g, %g]\n", var, basis.nonbasicMove_[var],
        // simplex_lp.colLower_[col], simplex_lp.colUpper_[col]);
        if (simplex_lp.colLower_[col] == simplex_lp.colUpper_[col]) {
#ifdef HiGHSDEV
          if (!highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
          {
            cstat[col] = (int) HighsBasisStatus::LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-simplex_lp.colLower_[col]) && highs_isInfinity(simplex_lp.colUpper_[col]))
#endif
          {
            cstat[col] = (int) HighsBasisStatus::ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: col=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          col, basis.nonbasicFlag_[var], basis.nonbasicMove_[var], simplex_lp.colLower_[col],
          simplex_lp.colUpper_[col]);
#endif
      return col + 1;
    }
  }
  if (rstat != NULL) {
    for (int row = 0; row < simplex_lp.numRow_; row++) {
      int var = simplex_lp.numCol_ + row;
      if (!basis.nonbasicFlag_[var]) {
        rstat[row] = (int) HighsBasisStatus::BASIC;
        continue;
      }
      // NB nonbasicMove for rows refers to the solver's view where the bounds
      // are switched and negated
      else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN)
      // Free to move only down from -simplex_lp.rowLower_[row]
      {
#ifdef HiGHSDEV
        if (!highs_isInfinity(-simplex_lp.rowLower_[row]))
#endif
        {
          rstat[row] = (int) HighsBasisStatus::LOWER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP)
      // Free to move only up from -simplex_lp.rowUpper_[row]
      {
#ifdef HiGHSDEV
        if (!highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
        {
          rstat[row] = (int) HighsBasisStatus::UPPER;
          continue;
        }
      } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
        if (simplex_lp.rowLower_[row] == simplex_lp.rowUpper_[row]) {
#ifdef HiGHSDEV
          if (!highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
          {
            rstat[row] = (int) HighsBasisStatus::LOWER;
            continue;
          }
        } else {
#ifdef HiGHSDEV
          if (highs_isInfinity(-simplex_lp.rowLower_[row]) && highs_isInfinity(simplex_lp.rowUpper_[row]))
#endif
          {
            rstat[row] = (int) HighsBasisStatus::ZERO;
            continue;
          }
        }
      }
#ifdef HiGHSDEV
      printf(
          "Invalid basis status: row=%d, nonbasicFlag=%d, nonbasicMove=%d, "
          "lower=%g, upper=%g\n",
          row, basis.nonbasicFlag_[var], basis.nonbasicMove_[var], simplex_lp.rowLower_[row],
          simplex_lp.rowUpper_[row]);
#endif
      return -(row + 1);
    }
  }
  return 0;
}


#ifdef HiGHSDEV
void HighsSimplexInterface::check_load_from_postsolve() {
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  bool ok;

  ok = nonbasic_flag_basic_index_ok(simplex_lp, highs_model_object.basis_);
  printf("check_load_from_postsolve: return from nonbasic_flag_basic_index_ok = %d\n", ok);
  assert(ok);

  ok = all_nonbasic_move_vs_work_arrays_ok(highs_model_object);
  printf("check_load_from_postsolve: return from all_nonbasic_move_vs_work_arrays_ok = %d\n", ok);
  assert(ok);
}
#endif

