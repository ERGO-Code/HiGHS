/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInterface.cpp
 * @brief
 */
#include "Highs.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "simplex/HSimplex.h"
#include "util/HighsMatrixUtils.h"
#include "util/HighsSort.h"

HighsStatus Highs::addColsInterface(HighsInt XnumNewCol, const double* XcolCost,
                                    const double* XcolLower,
                                    const double* XcolUpper, HighsInt XnumNewNZ,
                                    const HighsInt* XAstart,
                                    const HighsInt* XAindex,
                                    const double* XAvalue) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsOptions& options = options_;
  if (XnumNewCol < 0) return HighsStatus::kError;
  if (XnumNewNZ < 0) return HighsStatus::kError;
  if (XnumNewCol == 0) return HighsStatus::kOk;
  if (XnumNewCol > 0)
    if (isColDataNull(options.log_options, XcolCost, XcolLower, XcolUpper))
      return HighsStatus::kError;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(options.log_options, XAstart, XAindex, XAvalue))
      return HighsStatus::kError;

  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsScale& scale = lp.scale_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;
  HighsLp& simplex_lp = ekk_instance_.lp_;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;

  bool& valid_basis = basis.valid;
  bool valid_simplex_lp = false;//simplex_status.valid;
  bool& valid_simplex_basis = simplex_status.has_basis;
  bool& lp_has_scaling = lp.scale_.has_scaling;

  // Check that if nonzeros are to be added then the model has a positive number
  // of rows
  if (lp.num_row_ <= 0 && XnumNewNZ > 0) return HighsStatus::kError;
  if (valid_simplex_lp && (simplex_lp.num_row_ <= 0 && XnumNewNZ > 0))
    return HighsStatus::kError;

  // Record the new number of columns
  HighsInt newNumCol = lp.num_col_ + XnumNewCol;

  HighsIndexCollection index_collection;
  index_collection.dimension_ = XnumNewCol;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = XnumNewCol - 1;

  // Take a copy of the cost and bounds that can be normalised
  std::vector<double> local_colCost{XcolCost, XcolCost + XnumNewCol};
  std::vector<double> local_colLower{XcolLower, XcolLower + XnumNewCol};
  std::vector<double> local_colUpper{XcolUpper, XcolUpper + XnumNewCol};

  return_status =
      interpretCallStatus(assessCosts(options, lp.num_col_, index_collection,
                                      local_colCost, options.infinite_cost),
                          return_status, "assessCosts");
  if (return_status == HighsStatus::kError) return return_status;
  // Assess the column bounds
  return_status = interpretCallStatus(
      assessBounds(options, "Col", lp.num_col_, index_collection,
                   local_colLower, local_colUpper, options.infinite_bound),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;
  // Append the columns to the LP vectors and matrix
  appendColsToLpVectors(lp, XnumNewCol, local_colCost, local_colLower, local_colUpper);

  // Form a column-wise HighsSparseMatrix of the new matrix columns so
  // that is is easy to handle and, if there are nonzeros, it can be
  // normalised
  HighsSparseMatrix local_a_matrix;
  local_a_matrix.num_col_ = XnumNewCol;
  local_a_matrix.num_row_ = lp.num_row_;
  local_a_matrix.format_ = MatrixFormat::kColwise;
  if (XnumNewNZ) {
    local_a_matrix.start_ = {XAstart, XAstart + XnumNewCol};
    local_a_matrix.start_.resize(XnumNewCol + 1);
    local_a_matrix.start_[XnumNewCol] = XnumNewNZ;
    local_a_matrix.index_ = {XAindex, XAindex + XnumNewNZ};
    local_a_matrix.value_ = {XAvalue, XAvalue + XnumNewNZ};
    // Assess the matrix rows
    return_status =
      interpretCallStatus(local_a_matrix.assess(options.log_options, "LP",
						 options.small_matrix_value,
						 options.large_matrix_value),
			  return_status, "assessMatrix");
    if (return_status == HighsStatus::kError) return return_status;
  } else {
    // No nonzeros so, whether the constraint matrix is column-wise or
    // row-wise, adding the empty matrix is trivial. Complete the
    // setup of an empty column-wise HighsSparseMatrix of the new
    // matrix columns
    local_a_matrix.start_.assign(XnumNewCol + 1, 0);
  }
  // Append the columns to LP matrix
  return_status = interpretCallStatus(lp.a_matrix_.addCols(local_a_matrix), return_status, "lp.a_matrix_.addCols");
  if (return_status == HighsStatus::kError) return return_status;
  if (lp_has_scaling) {
    assert(1==0);
    // Extend the column scaling factors
    scale.col.resize(newNumCol);
    for (HighsInt iCol = 0; iCol < XnumNewCol; iCol++)
      scale.col[lp.num_col_ + iCol] = 1.0;
    scale.num_col = newNumCol;
    // Apply the existing row scaling to the new columns
    HighsSparseMatrix alt_local_a_matrix = local_a_matrix;
    applyScalingToMatrix(scale.row, XnumNewCol,
			 local_a_matrix.start_,
			 local_a_matrix.index_,
			 local_a_matrix.value_);
    alt_local_a_matrix.applyRowScale(scale);
    assert(alt_local_a_matrix==local_a_matrix);
    // Consider applying column scaling to the new columns. 
    alt_local_a_matrix.considerColScaling(options.allowed_simplex_matrix_scale_factor,
					   &scale.col[lp.num_col_]);
    // Use colScaleMatrix to take the column-wise matrix and then treat it
    // row-wise
    colScaleMatrix(options.allowed_simplex_matrix_scale_factor,
		   &scale.col[lp.num_col_], XnumNewCol,
		   local_a_matrix.start_,
		   local_a_matrix.index_,
		   local_a_matrix.value_);
    assert(alt_local_a_matrix==local_a_matrix);
  }
  return_status = interpretCallStatus(ekk_instance_.addRows(local_a_matrix), return_status, "ekk_instance_.addRows");
  if (return_status == HighsStatus::kError) return return_status;
  // Update the basis correponding to new nonbasic columns
  if (valid_basis) appendNonbasicColsToBasis(lp, basis, XnumNewCol);
  // Deduce the consequences of adding new columns
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  ekk_instance_.updateStatus(LpAction::kNewCols);

  // Increase the number of columns in the LPs
  lp.num_col_ += XnumNewCol;
  assert(lp.dimensionsOk("addCols"));
  return return_status;
}

HighsStatus Highs::addRowsInterface(HighsInt XnumNewRow,
                                    const double* XrowLower,
                                    const double* XrowUpper, HighsInt XnumNewNZ,
                                    const HighsInt* XARstart,
                                    const HighsInt* XARindex,
                                    const double* XARvalue) {
  // addRows is fundamentally different from addCols, since the new
  // matrix data are held row-wise, so we have to insert data into the
  // column-wise matrix of the LP.
  HighsStatus return_status = HighsStatus::kOk;
  HighsOptions& options = options_;
  if (XnumNewRow < 0) return HighsStatus::kError;
  if (XnumNewNZ < 0) return HighsStatus::kError;
  if (XnumNewRow == 0) return HighsStatus::kOk;
  if (XnumNewRow > 0)
    if (isRowDataNull(options.log_options, XrowLower, XrowUpper))
      return HighsStatus::kError;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(options.log_options, XARstart, XARindex, XARvalue))
      return HighsStatus::kError;

  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsScale& scale = lp.scale_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;
  HighsLp& simplex_lp = ekk_instance_.lp_;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;

  bool& valid_basis = basis.valid;
  bool valid_simplex_lp = false;//simplex_status.valid;
  bool& valid_simplex_basis = simplex_status.has_basis;
  bool& lp_has_scaling = lp.scale_.has_scaling;

  // Check that if nonzeros are to be added then the model has a positive number
  // of columns
  if (lp.num_col_ <= 0 && XnumNewNZ > 0) return HighsStatus::kError;
  if (valid_simplex_lp && (simplex_lp.num_col_ <= 0 && XnumNewNZ > 0))
    return HighsStatus::kError;

  // Record the new number of rows
  HighsInt newNumRow = lp.num_row_ + XnumNewRow;

  HighsIndexCollection index_collection;
  index_collection.dimension_ = XnumNewRow;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = XnumNewRow - 1;
  // Take a copy of the bounds that can be normalised
  std::vector<double> local_rowLower{XrowLower, XrowLower + XnumNewRow};
  std::vector<double> local_rowUpper{XrowUpper, XrowUpper + XnumNewRow};

  return_status = interpretCallStatus(
      assessBounds(options, "Row", lp.num_row_, index_collection,
                   local_rowLower, local_rowUpper, options.infinite_bound),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;

  // Append the rows to the LP vectors
  appendRowsToLpVectors(lp, XnumNewRow, local_rowLower, local_rowUpper);

  // Form a row-wise HighsSparseMatrix of the new matrix rows so that
  // is is easy to handle and, if there are nonzeros, it can be
  // normalised
  HighsSparseMatrix local_ar_matrix;
  local_ar_matrix.num_col_ = lp.num_col_;
  local_ar_matrix.num_row_ = XnumNewRow;
  local_ar_matrix.format_ = MatrixFormat::kRowwise;
  if (XnumNewNZ) {
    local_ar_matrix.start_ = {XARstart, XARstart + XnumNewRow};
    local_ar_matrix.start_.resize(XnumNewRow + 1);
    local_ar_matrix.start_[XnumNewRow] = XnumNewNZ;
    local_ar_matrix.index_ = {XARindex, XARindex + XnumNewNZ};
    local_ar_matrix.value_ = {XARvalue, XARvalue + XnumNewNZ};
    // Assess the matrix columns
    return_status =
      interpretCallStatus(local_ar_matrix.assess(options.log_options, "LP",
						 options.small_matrix_value,
						 options.large_matrix_value),
			  return_status, "assessMatrix");
    if (return_status == HighsStatus::kError) return return_status;
  } else {
    // No nonzeros so, whether the constraint matrix is row-wise or
    // column-wise, adding the empty matrix is trivial. Complete the
    // setup of an empty row-wise HighsSparseMatrix of the new matrix
    // rows
    local_ar_matrix.start_.assign(XnumNewRow + 1, 0);
  }
  // Append the rows to LP matrix
  return_status = interpretCallStatus(lp.a_matrix_.addRows(local_ar_matrix), return_status, "lp.a_matrix_.addRows");
  if (return_status == HighsStatus::kError) return return_status;
  if (lp_has_scaling) {
    assert(1==0);
    // Extend the row scaling factors
    scale.row.resize(newNumRow);
    for (HighsInt iRow = 0; iRow < XnumNewRow; iRow++)
      scale.row[lp.num_row_ + iRow] = 1.0;
    scale.num_row = newNumRow;
    // Apply the existing column scaling to the new rows
    HighsSparseMatrix alt_local_ar_matrix = local_ar_matrix;
    applyScalingToMatrix(scale.col, XnumNewRow,
			 local_ar_matrix.start_,
			 local_ar_matrix.index_,
			 local_ar_matrix.value_);
    alt_local_ar_matrix.applyColScale(scale);
    assert(alt_local_ar_matrix==local_ar_matrix);
    // Consider applying row scaling to the new rows. 
    alt_local_ar_matrix.considerRowScaling(options.allowed_simplex_matrix_scale_factor,
					   &scale.row[lp.num_row_]);
    // Use colScaleMatrix to take the row-wise matrix and then treat it
    // col-wise
    colScaleMatrix(options.allowed_simplex_matrix_scale_factor,
		   &scale.row[lp.num_row_], XnumNewRow,
		   local_ar_matrix.start_,
		   local_ar_matrix.index_,
		   local_ar_matrix.value_);
    assert(alt_local_ar_matrix==local_ar_matrix);
  }
  return_status = interpretCallStatus(ekk_instance_.addRows(local_ar_matrix), return_status, "ekk_instance_.addRows");
  if (return_status == HighsStatus::kError) return return_status;
  // Update the basis correponding to new basic rows
  if (valid_basis) appendBasicRowsToBasis(lp, basis, XnumNewRow);

  // Deduce the consequences of adding new rows
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  // ToDo clear highs_info_

  // Increase the number of rows in the LPs
  lp.num_row_ += XnumNewRow;
  assert(lp.dimensionsOk("addRows"));
  return return_status;
}

HighsStatus Highs::deleteColsInterface(HighsIndexCollection& index_collection) {
  HighsOptions& options = options_;
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;
  // Query: should simplex_status.valid be simplex_status.valid_?
  // Ensure that the LP is column-wise
  lp.ensureColWise();

  bool valid_simplex_lp = false;//simplex_status.valid;
  // Keep a copy of the original number of columns to check whether
  // any columns have been removed, and if there is mask to be updated
  HighsInt original_num_col = lp.num_col_;

  HighsStatus return_status;
  deleteLpCols(lp, index_collection);
  assert(lp.num_col_ <= original_num_col);
  if (lp.num_col_ < original_num_col) {
    // Nontrivial deletion so reset the model_status and invalidate
    // the Highs basis
    scaled_model_status_ = HighsModelStatus::kNotset;
    model_status_ =
        scaled_model_status_;
    basis.valid = false;
  }
  if (lp.scale_.has_scaling) {
    deleteScale(lp.scale_.col, index_collection);
    lp.scale_.col.resize(lp.num_col_);
  }
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = ekk_instance_.lp_;
    deleteLpCols(simplex_lp, index_collection);
    assert(simplex_lp.num_col_ <= original_num_col);
    if (simplex_lp.num_col_ < original_num_col) {
      // Nontrivial deletion so initialise the random vectors and all
      // data relating to the simplex basis
      ekk_instance_.initialiseSimplexLpRandomVectors();
      ekk_instance_.invalidateBasis();
    }
  }
  if (index_collection.is_mask_) {
    // Set the mask values to indicate the new index value of the
    // remaining columns
    HighsInt new_col = 0;
    for (HighsInt col = 0; col < original_num_col; col++) {
      if (!index_collection.mask_[col]) {
        index_collection.mask_[col] = new_col;
        new_col++;
      } else {
        index_collection.mask_[col] = -1;
      }
    }
    assert(new_col == lp.num_col_);
  }
  assert(lp.dimensionsOk("deleteCols"));
  if (valid_simplex_lp)
    assert(ekk_instance_.lp_.dimensionsOk("deleteCols - simplex"));
  return HighsStatus::kOk;
}

HighsStatus Highs::deleteRowsInterface(HighsIndexCollection& index_collection) {
  HighsOptions& options = options_;
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;
  // Query: should simplex_status.valid be simplex_status.valid_?
  // Ensure that the LP (and any simplex LP) is column-wise
  // Ensure that the LP is column-wise
  lp.ensureColWise();

  bool valid_simplex_lp = false;//simplex_status.valid;
  // Keep a copy of the original number of rows to check whether
  // any rows have been removed, and if there is mask to be updated
  HighsInt original_num_row = lp.num_row_;

  HighsStatus return_status;
  deleteLpRows(lp, index_collection);
  assert(lp.num_row_ <= original_num_row);
  if (lp.num_row_ < original_num_row) {
    // Nontrivial deletion so reset the model_status and invalidate
    // the Highs basis
    scaled_model_status_ = HighsModelStatus::kNotset;
    model_status_ = scaled_model_status_;
    basis.valid = false;
  }
  if (lp.scale_.has_scaling) {
    deleteScale(lp.scale_.row, index_collection);
    lp.scale_.row.resize(lp.num_row_);
  }
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = ekk_instance_.lp_;
    deleteLpRows(simplex_lp, index_collection);
    assert(simplex_lp.num_row_ <= original_num_row);
    if (simplex_lp.num_row_ < original_num_row) {
      // Nontrivial deletion so initialise the random vectors and all
      // data relating to the simplex basis
      ekk_instance_.initialiseSimplexLpRandomVectors();
      ekk_instance_.invalidateBasis();
    }
  }
  if (index_collection.is_mask_) {
    HighsInt new_row = 0;
    for (HighsInt row = 0; row < original_num_row; row++) {
      if (!index_collection.mask_[row]) {
        index_collection.mask_[row] = new_row;
        new_row++;
      } else {
        index_collection.mask_[row] = -1;
      }
    }
    assert(new_row == lp.num_row_);
  }
  assert(lp.dimensionsOk("deleteRows"));
  if (valid_simplex_lp)
    assert(ekk_instance_.lp_.dimensionsOk("deleteRows - simplex"));
  return HighsStatus::kOk;
}

HighsStatus Highs::getColsInterface(
    const HighsIndexCollection& index_collection, HighsInt& num_col,
    double* col_cost, double* col_lower, double* col_upper, HighsInt& num_nz,
    HighsInt* col_matrix_start, HighsInt* col_matrix_index,
    double* col_matrix_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsLp& lp = model_.lp_;
  HighsOptions& options = options_;
  // Ensure that the LP is column-wise
  lp.ensureColWise();
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k < 0 || to_k > lp.num_col_) {
    call_status = HighsStatus::kError;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  if (from_k > to_k) {
    call_status = HighsStatus::kError;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  HighsInt out_from_col;
  HighsInt out_to_col;
  HighsInt in_from_col;
  HighsInt in_to_col = -1;
  HighsInt current_set_entry = 0;
  HighsInt col_dim = lp.num_col_;
  num_col = 0;
  num_nz = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, out_from_col, out_to_col,
                                    in_from_col, in_to_col, current_set_entry);
    assert(out_to_col < col_dim);
    assert(in_to_col < col_dim);
    for (HighsInt col = out_from_col; col <= out_to_col; col++) {
      if (col_cost != NULL) col_cost[num_col] = lp.col_cost_[col];
      if (col_lower != NULL) col_lower[num_col] = lp.col_lower_[col];
      if (col_upper != NULL) col_upper[num_col] = lp.col_upper_[col];
      if (col_matrix_start != NULL)
        col_matrix_start[num_col] = num_nz + lp.a_matrix_.start_[col] -
                                    lp.a_matrix_.start_[out_from_col];
      num_col++;
    }
    for (HighsInt el = lp.a_matrix_.start_[out_from_col];
         el < lp.a_matrix_.start_[out_to_col + 1]; el++) {
      if (col_matrix_index != NULL)
        col_matrix_index[num_nz] = lp.a_matrix_.index_[el];
      if (col_matrix_value != NULL)
        col_matrix_value[num_nz] = lp.a_matrix_.value_[el];
      num_nz++;
    }
    if (out_to_col == col_dim - 1 || in_to_col == col_dim - 1) break;
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getRowsInterface(
    const HighsIndexCollection& index_collection, HighsInt& num_row,
    double* row_lower, double* row_upper, HighsInt& num_nz,
    HighsInt* row_matrix_start, HighsInt* row_matrix_index,
    double* row_matrix_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsLp& lp = model_.lp_;
  HighsOptions& options = options_;
  // Ensure that the LP is column-wise
  lp.ensureColWise();
  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  if (from_k < 0 || to_k > lp.num_row_) {
    call_status = HighsStatus::kError;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  num_row = 0;
  num_nz = 0;
  if (from_k > to_k) {
    call_status = HighsStatus::kError;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  // "Out" means not in the set to be extrated
  // "In" means in the set to be extrated
  HighsInt out_from_row;
  HighsInt out_to_row;
  HighsInt in_from_row;
  HighsInt in_to_row = -1;
  HighsInt current_set_entry = 0;
  HighsInt row_dim = lp.num_row_;
  // Ensure that the LP is column-wise
  lp.ensureColWise();
  // Set up a row mask so that entries to be got from the column-wise
  // matrix can be identified and have their correct row index.
  vector<HighsInt> new_index;
  new_index.resize(lp.num_row_);

  if (!index_collection.is_mask_) {
    out_to_row = -1;
    current_set_entry = 0;
    for (HighsInt k = from_k; k <= to_k; k++) {
      updateOutInIndex(index_collection, in_from_row, in_to_row,
                                      out_from_row, out_to_row,
                                      current_set_entry);
      if (k == from_k) {
        // Account for any initial rows not being extracted
        for (HighsInt row = 0; row < in_from_row; row++) {
          new_index[row] = -1;
        }
      }
      for (HighsInt row = in_from_row; row <= in_to_row; row++) {
        new_index[row] = num_row;
        num_row++;
      }
      for (HighsInt row = out_from_row; row <= out_to_row; row++) {
        new_index[row] = -1;
      }
      if (out_to_row >= row_dim - 1) break;
    }
  } else {
    for (HighsInt row = 0; row < lp.num_row_; row++) {
      if (index_collection.mask_[row]) {
        new_index[row] = num_row;
        num_row++;
      } else {
        new_index[row] = -1;
      }
    }
  }

  // Bail out if no rows are to be extracted
  if (num_row == 0) return HighsStatus::kOk;

  // Allocate an array of lengths for the row-wise matrix to be extracted
  vector<HighsInt> row_matrix_length;
  row_matrix_length.resize(num_row);

  for (HighsInt row = 0; row < lp.num_row_; row++) {
    HighsInt new_row = new_index[row];
    if (new_row >= 0) {
      assert(new_row < num_row);
      if (row_lower != NULL) row_lower[new_row] = lp.row_lower_[row];
      if (row_upper != NULL) row_upper[new_row] = lp.row_upper_[row];
      row_matrix_length[new_row] = 0;
    }
  }
  // Identify the lengths of the rows in the row-wise matrix to be extracted
  for (HighsInt col = 0; col < lp.num_col_; col++) {
    for (HighsInt el = lp.a_matrix_.start_[col];
         el < lp.a_matrix_.start_[col + 1]; el++) {
      HighsInt row = lp.a_matrix_.index_[el];
      HighsInt new_row = new_index[row];
      if (new_row >= 0) row_matrix_length[new_row]++;
    }
  }

  if (row_matrix_start == NULL) {
    // If the matrix start vector is null then don't get values of
    // indices, otherwise both are meaningless
    if (row_matrix_index != NULL || row_matrix_value != NULL) {
      highsLogUser(options_.log_options,
                   HighsLogType::kError,
                   "Cannot supply meaningful row matrix indices/values with "
                   "null starts\n");
      return HighsStatus::kError;
    }
  } else {
    row_matrix_start[0] = 0;
    for (HighsInt row = 0; row < num_row - 1; row++) {
      row_matrix_start[row + 1] =
          row_matrix_start[row] + row_matrix_length[row];
    }

    // Fill the row-wise matrix with indices and values
    for (HighsInt col = 0; col < lp.num_col_; col++) {
      for (HighsInt el = lp.a_matrix_.start_[col];
           el < lp.a_matrix_.start_[col + 1]; el++) {
        HighsInt row = lp.a_matrix_.index_[el];
        HighsInt new_row = new_index[row];
        if (new_row >= 0) {
          HighsInt row_el = row_matrix_start[new_row];
          if (row_matrix_index != NULL) row_matrix_index[row_el] = col;
          if (row_matrix_value != NULL)
            row_matrix_value[row_el] = lp.a_matrix_.value_[el];
          row_matrix_start[new_row]++;
        }
      }
    }
    // Restore the starts of the row-wise matrix and count the number of
    // nonzeros in it
    num_nz = 0;
    row_matrix_start[0] = 0;
    for (HighsInt row = 0; row < num_row - 1; row++) {
      row_matrix_start[row + 1] =
          row_matrix_start[row] + row_matrix_length[row];
      num_nz += row_matrix_length[row];
    }
    num_nz += row_matrix_length[num_row - 1];
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::getCoefficientInterface(const HighsInt Xrow,
                                           const HighsInt Xcol, double& value) {
  HighsLp& lp = model_.lp_;
  if (Xrow < 0 || Xrow >= lp.num_row_) return HighsStatus::kError;
  if (Xcol < 0 || Xcol >= lp.num_col_) return HighsStatus::kError;
  value = 0;
  // Ensure that the LP is column-wise
  lp.ensureColWise();
  for (HighsInt el = lp.a_matrix_.start_[Xcol];
       el < lp.a_matrix_.start_[Xcol + 1]; el++) {
    if (lp.a_matrix_.index_[el] == Xrow) {
      value = lp.a_matrix_.value_[el];
      break;
    }
  }
  return HighsStatus::kOk;
}

HighsStatus Highs::changeObjectiveSenseInterface(const ObjSense Xsense) {
  // If the sense doesn't change, just return
  if ((Xsense == ObjSense::kMinimize) ==
      (model_.lp_.sense_ == ObjSense::kMinimize))
    return HighsStatus::kOk;
  // Assume that objective sense changes
  // Set the LP objective sense
  model_.lp_.sense_ = Xsense;
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  return HighsStatus::kOk;
}

HighsStatus Highs::changeObjectiveOffsetInterface(const double Xoffset) {
  // If the offset doesn't change, just return
  if (Xoffset == model_.lp_.offset_) return HighsStatus::kOk;
  // Assume that objective offset changes
  // Update the objective value
  info_.objective_function_value += (Xoffset - model_.lp_.offset_);
  // Set the LP objective offset
  model_.lp_.offset_ = Xoffset;
  // Set any Simplex LP objective offset
  if (ekk_instance_.status_.valid)
    ekk_instance_.lp_.offset_ = Xoffset;
  return HighsStatus::kOk;
}

HighsStatus Highs::changeIntegralityInterface(
    HighsIndexCollection& index_collection,
    const HighsVarType* usr_integrality) {
  HighsOptions& options = options_;
  bool null_data = highsVarTypeUserDataNotNull(
      options.log_options, usr_integrality, "column integrality");
  if (null_data) return HighsStatus::kError;
  HighsInt num_usr_integrality = dataSize(index_collection);
  // If a non-positive number of integrality (may) need changing nothing needs
  // to be done
  if (num_usr_integrality <= 0) return HighsStatus::kOk;
  // Take a copy of the integrality that can be normalised
  std::vector<HighsVarType> local_integrality{
      usr_integrality, usr_integrality + num_usr_integrality};
  // If changing the integrality for a set of columns, ensure that the
  // set and data are in ascending order
  if (index_collection.is_set_)
    sortSetData(index_collection.set_num_entries_, index_collection.set_,
		usr_integrality,
		&local_integrality[0]);
  HighsLp& lp = model_.lp_;
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  changeLpIntegrality(lp, index_collection, local_integrality);

  // Deduce the consequences of new integrality
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  return HighsStatus::kOk;
}

HighsStatus Highs::changeCostsInterface(HighsIndexCollection& index_collection,
                                        const double* usr_col_cost) {
  HighsOptions& options = options_;
  bool null_data =
      doubleUserDataNotNull(options.log_options, usr_col_cost, "column costs");
  if (null_data) return HighsStatus::kError;
  HighsInt num_usr_col_cost = dataSize(index_collection);
  // If a non-positive number of costs (may) need changing nothing needs to be
  // done
  if (num_usr_col_cost <= 0) return HighsStatus::kOk;
  // Take a copy of the cost that can be normalised
  std::vector<double> local_colCost{usr_col_cost,
                                    usr_col_cost + num_usr_col_cost};
  // If changing the costs for a set of columns, ensure that the
  // set and data are in ascending order
  if (index_collection.is_set_)
    sortSetData(index_collection.set_num_entries_, index_collection.set_,
		usr_col_cost, NULL, NULL,
		&local_colCost[0], NULL, NULL);
  HighsLp& lp = model_.lp_;
  HighsStatus return_status = HighsStatus::kOk;
  return_status =
      interpretCallStatus(assessCosts(options, 0, index_collection,
                                      local_colCost, options.infinite_cost),
                          return_status, "assessCosts");
  if (return_status == HighsStatus::kError) return return_status;

  HighsStatus call_status;
  changeLpCosts(lp, index_collection, local_colCost);

  if (ekk_instance_.status_.valid) {
    // Also change the simplex LP's costs
    HighsLp& simplex_lp = ekk_instance_.lp_;
    assert(lp.num_col_ == simplex_lp.num_col_);
    assert(lp.num_row_ == simplex_lp.num_row_);
    changeLpCosts(simplex_lp, index_collection, local_colCost);
    if (lp.scale_.has_scaling) {
      applyScalingToLpColCost(simplex_lp,
                              lp.scale_.col, index_collection);
    }
  }
  // Deduce the consequences of new costs
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  ekk_instance_.updateStatus(LpAction::kNewCosts);
  return HighsStatus::kOk;
}

HighsStatus Highs::changeColBoundsInterface(
    HighsIndexCollection& index_collection, const double* usr_col_lower,
    const double* usr_col_upper) {
  HighsOptions& options = options_;
  bool null_data = false;
  null_data = doubleUserDataNotNull(options.log_options, usr_col_lower,
                                    "column lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(options.log_options, usr_col_upper,
                                    "column upper bounds") ||
              null_data;
  if (null_data) return HighsStatus::kError;
  HighsInt num_usr_col_bounds = dataSize(index_collection);
  // If a non-positive number of costs (may) need changing nothing needs to be
  // done
  if (num_usr_col_bounds <= 0) return HighsStatus::kOk;
  // Take a copy of the cost that can be normalised
  std::vector<double> local_colLower{usr_col_lower,
                                     usr_col_lower + num_usr_col_bounds};
  std::vector<double> local_colUpper{usr_col_upper,
                                     usr_col_upper + num_usr_col_bounds};
  // If changing the bounds for a set of columns, ensure that the
  // set and data are in ascending order
  if (index_collection.is_set_)
    sortSetData(index_collection.set_num_entries_, index_collection.set_,
                usr_col_lower, usr_col_upper, NULL,
		&local_colLower[0], &local_colUpper[0], NULL);
  HighsLp& lp = model_.lp_;
  HighsStatus return_status = HighsStatus::kOk;
  return_status = interpretCallStatus(
      assessBounds(options, "col", 0, index_collection, local_colLower,
                   local_colUpper, options.infinite_bound),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;

  HighsStatus call_status;
  changeLpColBounds(lp, index_collection,
		    local_colLower, local_colUpper);

  if (ekk_instance_.status_.valid) {
    // Also change the simplex LP's column bounds
    HighsLp& simplex_lp = ekk_instance_.lp_;
    assert(lp.num_col_ == simplex_lp.num_col_);
    assert(lp.num_row_ == simplex_lp.num_row_);
    changeLpColBounds(simplex_lp, index_collection,
		      local_colLower, local_colUpper);
    if (lp.scale_.has_scaling) {
      applyScalingToLpColBounds(simplex_lp,
                                lp.scale_.col,
                                index_collection);
    }
  }
  if (basis_.valid) {
    // Update HiGHS basis status and (any) simplex move status of
    // nonbasic variables whose bounds have changed
    return_status =
        interpretCallStatus(setNonbasicStatusInterface(index_collection, true),
                            return_status, "setNonbasicStatusInterface");
    if (return_status == HighsStatus::kError) return return_status;
  }

  // Deduce the consequences of new col bounds
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  ekk_instance_.updateStatus(LpAction::kNewBounds);
  return HighsStatus::kOk;
}

HighsStatus Highs::changeRowBoundsInterface(
    HighsIndexCollection& index_collection, const double* usr_row_lower,
    const double* usr_row_upper) {
  HighsOptions& options = options_;
  bool null_data = false;
  null_data = doubleUserDataNotNull(options.log_options, usr_row_lower,
                                    "row lower bounds") ||
              null_data;
  null_data = doubleUserDataNotNull(options.log_options, usr_row_upper,
                                    "row upper bounds") ||
              null_data;
  if (null_data) return HighsStatus::kError;
  HighsInt num_usr_row_bounds = dataSize(index_collection);
  // If a non-positive number of costs (may) need changing nothing needs to be
  // done
  if (num_usr_row_bounds <= 0) return HighsStatus::kOk;
  // Take a copy of the cost that can be normalised
  std::vector<double> local_rowLower{usr_row_lower,
                                     usr_row_lower + num_usr_row_bounds};
  std::vector<double> local_rowUpper{usr_row_upper,
                                     usr_row_upper + num_usr_row_bounds};
  // If changing the bounds for a set of rows, ensure that the
  // set and data are in ascending order
  if (index_collection.is_set_)
    sortSetData(index_collection.set_num_entries_, index_collection.set_,
                usr_row_lower, usr_row_upper, NULL, &local_rowLower[0],
                &local_rowUpper[0], NULL);
  HighsLp& lp = model_.lp_;
  HighsStatus return_status = HighsStatus::kOk;
  return_status = interpretCallStatus(
      assessBounds(options, "row", 0, index_collection, local_rowLower,
                   local_rowUpper, options.infinite_bound),
      return_status, "assessBounds");
  if (return_status == HighsStatus::kError) return return_status;

  changeLpRowBounds(lp, index_collection,
                                  local_rowLower, local_rowUpper);

  if (ekk_instance_.status_.valid) {
    // Also change the simplex LP's row bounds
    HighsLp& simplex_lp = ekk_instance_.lp_;
    assert(lp.num_col_ == simplex_lp.num_col_);
    assert(lp.num_row_ == simplex_lp.num_row_);
    changeLpRowBounds(simplex_lp, index_collection,
		      local_rowLower, local_rowUpper);
    if (lp.scale_.has_scaling) {
      applyScalingToLpRowBounds(simplex_lp,
                                lp.scale_.row,
                                index_collection);
    }
  }
  if (basis_.valid) {
    // Update HiGHS basis status and (any) simplex move status of
    // nonbasic variables whose bounds have changed
    return_status =
        interpretCallStatus(setNonbasicStatusInterface(index_collection, false),
                            return_status, "setNonbasicStatusInterface");
    if (return_status == HighsStatus::kError) return return_status;
  }
  // Deduce the consequences of new row bounds
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  ekk_instance_.updateStatus(LpAction::kNewBounds);
  return HighsStatus::kOk;
}

// Change a single coefficient in the matrix
HighsStatus Highs::changeCoefficientInterface(const HighsInt Xrow,
                                              const HighsInt Xcol,
                                              const double XnewValue) {
  HighsLp& lp = model_.lp_;
  // Ensure that the LP is column-wise
  lp.ensureColWise();
  if (Xrow < 0 || Xrow >= lp.num_row_) return HighsStatus::kError;
  if (Xcol < 0 || Xcol >= lp.num_col_) return HighsStatus::kError;
  changeLpMatrixCoefficient(lp, Xrow, Xcol, XnewValue);
  // Deduce the consequences of a changed element
  // ToDo: Can do something more intelligent if element is in nonbasic column.
  // Otherwise, treat it as if it's a new row
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  ekk_instance_.updateStatus(LpAction::kNewRows);
  return HighsStatus::kOk;
}

HighsStatus Highs::scaleColInterface(const HighsInt col,
                                     const double scaleval) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsOptions& options = options_;
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;
  HighsLp& simplex_lp = ekk_instance_.lp_;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;

  // Ensure that the LP is column-wise
  lp.ensureColWise();
  if (col < 0) return HighsStatus::kError;
  if (col >= lp.num_col_) return HighsStatus::kError;
  if (!scaleval) return HighsStatus::kError;

  return_status = interpretCallStatus(
      applyScalingToLpCol(lp, col, scaleval),
      return_status, "applyScalingToLpCol");
  if (return_status == HighsStatus::kError) return return_status;

  if (scaleval < 0 && basis.valid) {
    // Negative, so flip any nonbasic status
    if (basis.col_status[col] == HighsBasisStatus::kLower) {
      basis.col_status[col] = HighsBasisStatus::kUpper;
    } else if (basis.col_status[col] == HighsBasisStatus::kUpper) {
      basis.col_status[col] = HighsBasisStatus::kLower;
    }
  }
  if (simplex_status.valid) {
    // Apply the scaling to the simplex LP
    return_status = interpretCallStatus(
        applyScalingToLpCol(simplex_lp, col, scaleval),
        return_status, "applyScalingToLpCol");
    if (return_status == HighsStatus::kError) return return_status;
    if (scaleval < 0 && simplex_status.has_basis) {
      // Negative, so flip any nonbasic status
      if (simplex_basis.nonbasicMove_[col] == kNonbasicMoveUp) {
        simplex_basis.nonbasicMove_[col] = kNonbasicMoveDn;
      } else if (simplex_basis.nonbasicMove_[col] == kNonbasicMoveDn) {
        simplex_basis.nonbasicMove_[col] = kNonbasicMoveUp;
      }
    }
  }

  // Deduce the consequences of a scaled column
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  ekk_instance_.updateStatus(LpAction::kScaledCol);
  return HighsStatus::kOk;
}

HighsStatus Highs::scaleRowInterface(const HighsInt row,
                                     const double scaleval) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsOptions& options = options_;
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  HighsSimplexStatus& simplex_status = ekk_instance_.status_;
  HighsLp& simplex_lp = ekk_instance_.lp_;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;

  // Ensure that the LP is column-wise
  lp.ensureColWise();

  if (row < 0) return HighsStatus::kError;
  if (row >= lp.num_row_) return HighsStatus::kError;
  if (!scaleval) return HighsStatus::kError;

  return_status = interpretCallStatus(
      applyScalingToLpRow(lp, row, scaleval),
      return_status, "applyScalingToLpRow");
  if (return_status == HighsStatus::kError) return return_status;

  if (scaleval < 0 && basis.valid) {
    // Negative, so flip any nonbasic status
    if (basis.row_status[row] == HighsBasisStatus::kLower) {
      basis.row_status[row] = HighsBasisStatus::kUpper;
    } else if (basis.row_status[row] == HighsBasisStatus::kUpper) {
      basis.row_status[row] = HighsBasisStatus::kLower;
    }
  }
  if (simplex_status.valid) {
    // Apply the scaling to the simplex LP
    return_status = interpretCallStatus(
        applyScalingToLpRow(simplex_lp, row, scaleval),
        return_status, "applyScalingToLpRow");
    if (return_status == HighsStatus::kError) return return_status;
    if (scaleval < 0 && simplex_status.has_basis) {
      // Negative, so flip any nonbasic status
      const HighsInt var = simplex_lp.num_col_ + row;
      if (simplex_basis.nonbasicMove_[var] == kNonbasicMoveUp) {
        simplex_basis.nonbasicMove_[var] = kNonbasicMoveDn;
      } else if (simplex_basis.nonbasicMove_[var] == kNonbasicMoveDn) {
        simplex_basis.nonbasicMove_[var] = kNonbasicMoveUp;
      }
    }
  }

  // Deduce the consequences of a scaled row
  scaled_model_status_ = HighsModelStatus::kNotset;
  model_status_ = scaled_model_status_;
  ekk_instance_.updateStatus(LpAction::kScaledRow);
  return HighsStatus::kOk;
}

HighsStatus Highs::setNonbasicStatusInterface(
    const HighsIndexCollection& index_collection, const bool columns) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsLp& lp = model_.lp_;
  HighsBasis& basis = basis_;
  SimplexBasis& simplex_basis = ekk_instance_.basis_;
  HighsOptions& options = options_;

  assert(basis.valid);
  const bool has_simplex_basis = ekk_instance_.status_.has_basis;

  assert(ok(index_collection));
  HighsInt from_k;
  HighsInt to_k;
  limits(index_collection, from_k, to_k);
  HighsInt ix_dim;
  if (columns) {
    ix_dim = lp.num_col_;
  } else {
    ix_dim = lp.num_row_;
  }
  if (from_k < 0 || to_k > ix_dim) {
    call_status = HighsStatus::kError;
    return_status = interpretCallStatus(call_status, return_status,
                                        "setNonbasicStatusInterface");
    return return_status;
  }
  if (from_k > to_k) {
    call_status = HighsStatus::kError;
    return_status = interpretCallStatus(call_status, return_status,
                                        "setNonbasicStatusInterface");
    return return_status;
  }
  HighsInt set_from_ix;
  HighsInt set_to_ix;
  HighsInt ignore_from_ix;
  HighsInt ignore_to_ix = -1;
  HighsInt current_set_entry = 0;
  for (HighsInt k = from_k; k <= to_k; k++) {
    updateOutInIndex(index_collection, set_from_ix, set_to_ix,
                                    ignore_from_ix, ignore_to_ix,
                                    current_set_entry);
    assert(set_to_ix < ix_dim);
    assert(ignore_to_ix < ix_dim);
    if (columns) {
      for (HighsInt iCol = set_from_ix; iCol <= set_to_ix; iCol++) {
        if (basis.col_status[iCol] == HighsBasisStatus::kBasic) continue;
        // Nonbasic column
        double lower = lp.col_lower_[iCol];
        double upper = lp.col_upper_[iCol];
        if (!highs_isInfinity(-lower)) {
          basis.col_status[iCol] = HighsBasisStatus::kLower;
        } else if (!highs_isInfinity(upper)) {
          basis.col_status[iCol] = HighsBasisStatus::kUpper;
        } else {
          basis.col_status[iCol] = HighsBasisStatus::kZero;
        }
        if (has_simplex_basis) {
          // todo @ Julian this assert fails on glass4
          assert(simplex_basis.nonbasicFlag_[iCol] == kNonbasicFlagTrue);
          HighsInt move = kIllegalMoveValue;
          if (lower == upper) {
            move = kNonbasicMoveZe;
          } else if (!highs_isInfinity(-lower)) {
            // Finite lower bound so boxed or lower
            if (!highs_isInfinity(upper)) {
              // Finite upper bound so boxed
              if (fabs(lower) < fabs(upper)) {
                move = kNonbasicMoveUp;
              } else {
                move = kNonbasicMoveDn;
              }
            } else {
              // Lower (since upper bound is infinite)
              move = kNonbasicMoveUp;
            }
          } else if (!highs_isInfinity(upper)) {
            // Upper
            move = kNonbasicMoveDn;
          } else {
            // FREE
            move = kNonbasicMoveZe;
          }
          assert(move != kIllegalMoveValue);
          simplex_basis.nonbasicMove_[iCol] = move;
        }
      }
    } else {
      for (HighsInt iRow = set_from_ix; iRow <= set_to_ix; iRow++) {
        if (basis.row_status[iRow] == HighsBasisStatus::kBasic) continue;
        // Nonbasic column
        double lower = lp.row_lower_[iRow];
        double upper = lp.row_upper_[iRow];
        if (!highs_isInfinity(-lower)) {
          basis.row_status[iRow] = HighsBasisStatus::kLower;
        } else if (!highs_isInfinity(upper)) {
          basis.row_status[iRow] = HighsBasisStatus::kUpper;
        } else {
          basis.row_status[iRow] = HighsBasisStatus::kZero;
        }
        if (has_simplex_basis) {
          assert(simplex_basis.nonbasicFlag_[lp.num_col_ + iRow] ==
                 kNonbasicFlagTrue);
          HighsInt move = kIllegalMoveValue;
          if (lower == upper) {
            move = kNonbasicMoveZe;
          } else if (!highs_isInfinity(-lower)) {
            // Finite lower bound so boxed or lower
            if (!highs_isInfinity(upper)) {
              // Finite upper bound so boxed
              if (fabs(lower) < fabs(upper)) {
                move = kNonbasicMoveDn;
              } else {
                move = kNonbasicMoveUp;
              }
            } else {
              // Lower (since upper bound is infinite)
              move = kNonbasicMoveDn;
            }
          } else if (!highs_isInfinity(upper)) {
            // Upper
            move = kNonbasicMoveUp;
          } else {
            // FREE
            move = kNonbasicMoveZe;
          }
          assert(move != kIllegalMoveValue);
          simplex_basis.nonbasicMove_[lp.num_col_ + iRow] = move;
        }
      }
    }
    if (ignore_to_ix >= ix_dim - 1) break;
  }
  return return_status;
}

void Highs::clearBasisInterface() {
  ekk_instance_.updateStatus(LpAction::kNewBasis);
}

// Get the basic variables, performing INVERT if necessary
HighsStatus Highs::getBasicVariablesInterface(HighsInt* basic_variables) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsLp& ekk_lp = ekk_instance_.lp_;
  HighsInt num_row = lp.num_row_;
  HighsInt num_col = lp.num_col_;
  HighsSimplexStatus& ekk_status = ekk_instance_.status_;
  // For an LP with no rows the solution is vacuous
  if (num_row==0) return return_status;
  if (!ekk_status.has_invert) {
    // The LP has no invert to use, so have to set one up
    lp.ensureColWise();
    // Consider scaling the LP, and then move to EKK
    considerScaling(options_, lp);
    lp.moveLp(ekk_lp);
    ekk_instance_.setPointers(&options_, &timer_);
    // If the simplex LP isn't initialised, do so
    if (!ekk_status.initialised) {
      return_status = ekk_instance_.setup();
      if (return_status == HighsStatus::kError) {
	// Setup has failed - can only happen if there are excessive
	// matrix entries in an LP that should already have been
	// assessed - so move the LP back and unscale
	lp.moveLpBackAndUnapplyScaling(ekk_lp);
	return return_status;
      }
    }
    if (!ekk_status.has_basis) {
      //
      // The Ekk instance has no simplex basis, so pass the HiGHS basis
      // if it's valid, otherwise return an error for consistency with
      // old code
      //
      // Arguable that a warning should be issued and a logical basis
      // set up
      if (basis_.valid) {
	return_status = interpretCallStatus(ekk_instance_.setBasis(basis_),
					    return_status, "setBasis");
	if (return_status == HighsStatus::kError) return return_status;
      } else {
	highsLogUser(
		     options_.log_options, HighsLogType::kError,
		     "getBasicVariables called without a simplex or HiGHS basis\n");
	lp.moveLpBackAndUnapplyScaling(ekk_lp);
	return HighsStatus::kError;
      }
    }
    assert(ekk_status.has_basis);
    const bool only_from_known_basis = true;
    HighsInt return_value = ekk_instance_.initialiseSimplexLpBasisAndFactor(only_from_known_basis);
    lp.moveLpBackAndUnapplyScaling(ekk_lp);
    if (return_value) return HighsStatus::kError;
  }
  assert(ekk_status.has_invert);

  for (HighsInt row = 0; row < num_row; row++) {
    HighsInt var = ekk_instance_.basis_.basicIndex_[row];
    if (var < num_col) {
      basic_variables[row] = var;
    } else {
      basic_variables[row] = -(1 + var - num_col);
    }
  }
  return return_status;
}

// Solve (transposed) system involving the basis matrix

HighsStatus Highs::basisSolveInterface(const vector<double>& rhs,
                                       double* solution_vector,
                                       HighsInt* solution_num_nz,
                                       HighsInt* solution_indices,
                                       bool transpose) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsInt num_row = lp.num_row_;
  HighsInt num_col = lp.num_col_;
  // For an LP with no rows the solution is vacuous
  if (num_row==0) return return_status;
  assert(ekk_instance_.status_.has_invert);
  assert(!lp.is_moved_);
  // Set up solve vector with suitably scaled RHS
  HVector solve_vector;
  solve_vector.setup(num_row);
  solve_vector.clear();
  HighsScale& scale = lp.scale_;
  HighsInt rhs_num_nz = 0;
  if (transpose) {
    for (HighsInt row = 0; row < num_row; row++) {
      if (rhs[row]) {
        solve_vector.index[rhs_num_nz++] = row;
        double rhs_value = rhs[row];
	if (scale.has_scaling) {
	  HighsInt col = ekk_instance_.basis_.basicIndex_[row];
	  if (col < num_col) {
	    rhs_value *= scale.col[col];
	  } else {
	    double scale_value = scale.row[col - num_col];
	    rhs_value /= scale_value;
	  }
	}
        solve_vector.array[row] = rhs_value;
      }
    }
  } else {
    for (HighsInt row = 0; row < num_row; row++) {
      if (rhs[row]) {
        solve_vector.index[rhs_num_nz++] = row;
        solve_vector.array[row] = rhs[row];
	if (scale.has_scaling)
	  solve_vector.array[row] *= scale.row[row];
      }
    }
  }
  solve_vector.count = rhs_num_nz;
  //
  // Note that solve_vector.count is just used to determine whether
  // hyper-sparse solves should be used. The indices of the nonzeros
  // in the solution are always accumulated. There's no switch (such
  // as setting solve_vector.count = num_row+1) to not do this.
  //
  // Get expected_density from analysis during simplex solve.
  const double expected_density = 1;
  if (transpose) {
    ekk_instance_.simplex_nla_.btran(solve_vector, expected_density);
  } else {
    ekk_instance_.simplex_nla_.ftran(solve_vector, expected_density);
  }
  // Extract the solution
  if (solution_indices == NULL) {
    // Nonzeros in the solution not required
    if (solve_vector.count > num_row) {
      // Solution nonzeros not known
      for (HighsInt row = 0; row < num_row; row++) {
        solution_vector[row] = solve_vector.array[row];
      }
    } else {
      // Solution nonzeros are known
      for (HighsInt row = 0; row < num_row; row++) solution_vector[row] = 0;
      for (HighsInt ix = 0; ix < solve_vector.count; ix++) {
        HighsInt row = solve_vector.index[ix];
        solution_vector[row] = solve_vector.array[row];
      }
    }
  } else {
    // Nonzeros in the solution are required
    if (solve_vector.count > num_row) {
      // Solution nonzeros not known
      solution_num_nz = 0;
      for (HighsInt row = 0; row < num_row; row++) {
        solution_vector[row] = 0;
        if (solve_vector.array[row]) {
          solution_vector[row] = solve_vector.array[row];
          solution_indices[*solution_num_nz++] = row;
        }
      }
    } else {
      // Solution nonzeros are known
      for (HighsInt row = 0; row < num_row; row++) solution_vector[row] = 0;
      for (HighsInt ix = 0; ix < solve_vector.count; ix++) {
        HighsInt row = solve_vector.index[ix];
        solution_vector[row] = solve_vector.array[row];
        solution_indices[ix] = row;
      }
      *solution_num_nz = solve_vector.count;
    }
  }
  // Scale the solution
  if (scale.has_scaling) {
    if (transpose) {
      if (solve_vector.count > num_row) {
	// Solution nonzeros not known
	for (HighsInt row = 0; row < num_row; row++) {
	  double scale_value = scale.row[row];
	  solution_vector[row] *= scale_value;
	}
      } else {
	for (HighsInt ix = 0; ix < solve_vector.count; ix++) {
	  HighsInt row = solve_vector.index[ix];
	  double scale_value = scale.row[row];
	  solution_vector[row] *= scale_value;
	}
      }
    } else {
      if (solve_vector.count > num_row) {
	// Solution nonzeros not known
	for (HighsInt row = 0; row < num_row; row++) {
	  HighsInt col = ekk_instance_.basis_.basicIndex_[row];
	  if (col < num_col) {
	    solution_vector[row] *= scale.col[col];
	  } else {
	    double scale_value = scale.row[col - num_col];
	    solution_vector[row] /= scale_value;
	  }
	}
      } else {
	for (HighsInt ix = 0; ix < solve_vector.count; ix++) {
	  HighsInt row = solve_vector.index[ix];
	  HighsInt col = ekk_instance_.basis_.basicIndex_[row];
	  if (col < num_col) {
	    solution_vector[row] *= scale.col[col];
	  } else {
	    double scale_value = scale.row[col - num_col];
	    solution_vector[row] /= scale_value;
	  }
	}
      }
    }
  }
  return HighsStatus::kOk;
}

void Highs::zeroIterationCounts() {
  info_.simplex_iteration_count = 0;
  info_.ipm_iteration_count = 0;
  info_.crossover_iteration_count = 0;
  info_.qp_iteration_count = 0;
}

HighsStatus Highs::getDualRayInterface(bool& has_dual_ray,
                                       double* dual_ray_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsInt num_row = lp.num_row_;
  // For an LP with no rows the dual ray is vacuous
  if (num_row==0) return return_status;
  assert(ekk_instance_.status_.has_invert);
  assert(!lp.is_moved_);
  has_dual_ray = ekk_instance_.status_.has_dual_ray;
  if (has_dual_ray && dual_ray_value != NULL) {
    vector<double> rhs;
    HighsInt iRow = ekk_instance_.info_.dual_ray_row_;
    rhs.assign(num_row, 0);
    rhs[iRow] = ekk_instance_.info_.dual_ray_sign_;
    HighsInt* dual_ray_num_nz = 0;
    basisSolveInterface(rhs, dual_ray_value, dual_ray_num_nz, NULL, true);
  }
  return return_status;
}

HighsStatus Highs::getPrimalRayInterface(bool& has_primal_ray,
                                         double* primal_ray_value) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsLp& lp = model_.lp_;
  HighsInt num_row = lp.num_row_;
  HighsInt num_col = lp.num_col_;
  // For an LP with no rows the primal ray is vacuous
  if (num_row==0) return return_status;
  assert(ekk_instance_.status_.has_invert);
  assert(!lp.is_moved_);
  has_primal_ray = ekk_instance_.status_.has_primal_ray;
  if (has_primal_ray && primal_ray_value != NULL) {
    HighsInt col = ekk_instance_.info_.primal_ray_col_;
    assert(ekk_instance_.basis_.nonbasicFlag_[col] == kNonbasicFlagTrue);
    // Get this pivotal column
    vector<double> rhs;
    vector<double> column;
    column.assign(num_row, 0);
    rhs.assign(num_row, 0);
    lp.ensureColWise();
    HighsInt primal_ray_sign = ekk_instance_.info_.primal_ray_sign_;
    if (col < num_col) {
      for (HighsInt iEl = lp.a_matrix_.start_[col];
           iEl < lp.a_matrix_.start_[col + 1]; iEl++)
        rhs[lp.a_matrix_.index_[iEl]] =
            primal_ray_sign * lp.a_matrix_.value_[iEl];
    } else {
      rhs[col - num_col] = primal_ray_sign;
    }
    HighsInt* column_num_nz = 0;
    basisSolveInterface(rhs, &column[0], column_num_nz, NULL, false);
    // Now zero primal_ray_value and scatter the column according to
    // the basic variables.
    for (HighsInt iCol = 0; iCol < num_col; iCol++) primal_ray_value[iCol] = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      HighsInt iCol = ekk_instance_.basis_.basicIndex_[iRow];
      if (iCol < num_col) primal_ray_value[iCol] = column[iRow];
    }
    if (col < num_col) primal_ray_value[col] = -primal_ray_sign;
  }
  return return_status;
}

bool Highs::qFormatOk(const HighsInt num_nz, const HighsInt format) {
  if (!num_nz) return true;
  const bool ok_format = format == (HighsInt)MatrixFormat::kColwise;
  assert(ok_format);
  if (!ok_format)
    highsLogUser(
        options_.log_options, HighsLogType::kError,
        "Non-empty Hessian matrix has illegal format = %" HIGHSINT_FORMAT "\n",
        format);
  return ok_format;
}

void Highs::clearZeroHessian() {
  HighsHessian& hessian = model_.hessian_;
  if (hessian.dim_) {
    // Clear any zero Hessian
    if (hessian.numNz() == 0) {
      highsLogUser(options_.log_options, HighsLogType::kInfo,
                   "Hessian has dimension %" HIGHSINT_FORMAT
                   " but no nonzeros, so is ignored\n",
                   hessian.dim_);
      hessian.clear();
    }
  }
}

HighsStatus Highs::checkOptimality(const std::string solver_type, HighsStatus return_status) {
  // Check for infeasibility measures incompatible with optimality
  assert(return_status != HighsStatus::kError);
  // Cannot expect to have no dual_infeasibilities since the QP solver
  // (and, of course, the MIP solver) give no dual information
  if (info_.num_primal_infeasibilities == 0 &&
      info_.num_dual_infeasibilities <=0) return HighsStatus::kOk;
  HighsLogType log_type = HighsLogType::kWarning;
  return_status = HighsStatus::kWarning;
  if (info_.max_primal_infeasibility > sqrt(options_.primal_feasibility_tolerance) ||
      info_.max_dual_infeasibility > sqrt(options_.dual_feasibility_tolerance)) {
    // Check for gross errors
    log_type = HighsLogType::kError;
    return_status = HighsStatus::kError;
  }
  highsLogUser(options_.log_options, log_type,
                     "%s solver claims optimality, but with num/sum/max "
                     "primal(%" HIGHSINT_FORMAT
                     "/%g/%g) and dual(%" HIGHSINT_FORMAT
                     "/%g/%g) infeasibilities\n", solver_type.c_str(),
                     info_.num_primal_infeasibilities,
                     info_.sum_primal_infeasibilities,
                     info_.max_primal_infeasibility,
                     info_.num_dual_infeasibilities,
                     info_.sum_dual_infeasibilities,
                     info_.max_dual_infeasibility);
  return return_status;
}

HighsStatus Highs::invertRequirementError(std::string method_name) {
  assert(!ekk_instance_.status_.has_invert);
  highsLogUser(options_.log_options, HighsLogType::kError,
	       "No invertible representation for %s\n", method_name.c_str());
  return HighsStatus::kError;
}
