/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsInterface.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Highs.h"

#include "HConfig.h"
//#include "io/HMPSIO.h"
//#include "io/HighsIO.h"
#include "lp_data/HighsLpUtils.h"
//#include "lp_data/HighsModelUtils.h"
#include "simplex/HSimplex.h"
//#include "simplex/HSimplexDebug.h"
//#include "util/HighsSort.h"
//#include "util/HighsUtils.h"

HighsStatus Highs::addColsInterface(int XnumNewCol, const double* XcolCost, const double* XcolLower,
				    const double* XcolUpper, int XnumNewNZ, const int* XAstart,
				    const int* XAindex, const double* XAvalue) {
  HighsModelObject& highs_model_object = hmos_[0];
  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsStatus return_status = HighsStatus::OK;
  HighsOptions& options = highs_model_object.options_;
  if (XnumNewCol < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewCol == 0) return HighsStatus::OK;
  if (XnumNewCol > 0)
    if (isColDataNull(options, XcolCost, XcolLower, XcolUpper))
      return HighsStatus::Error;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(options, XAstart, XAindex, XAvalue))
      return HighsStatus::Error;

  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsScale& scale = highs_model_object.scale_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;
  HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;

  bool& valid_basis = basis.valid_;
  bool& valid_simplex_lp = simplex_lp_status.valid;
  bool& valid_simplex_basis = simplex_lp_status.has_basis;
  bool& scaled_simplex_lp = scale.is_scaled_;

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
    assert(!scaled_simplex_lp);
  }
#endif

  HighsIndexCollection index_collection;
  index_collection.dimension_ = XnumNewCol;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = XnumNewCol - 1;

  // Take a copy of the cost and bounds that can be normalised
  std::vector<double> local_colCost{XcolCost, XcolCost + XnumNewCol};
  std::vector<double> local_colLower{XcolLower, XcolLower + XnumNewCol};
  std::vector<double> local_colUpper{XcolUpper, XcolUpper + XnumNewCol};

  // There are sure to be new columns since XnumNewCol <= 0 is handled above
  // Assess the column costs
  assert(XnumNewCol > 0);
  return_status =
      interpretCallStatus(assessCosts(options, lp.numCol_, index_collection,
                                      local_colCost, options.infinite_cost),
                          return_status, "assessCosts");
  if (return_status == HighsStatus::Error) return return_status;
  // Assess the column bounds
  return_status = interpretCallStatus(
      assessBounds(options, "Col", lp.numCol_, index_collection, local_colLower,
                   local_colUpper, options.infinite_bound),
      return_status, "assessBounds");
  if (return_status == HighsStatus::Error) return return_status;
  // Append the columns to the LP vectors and matrix
  return_status =
      interpretCallStatus(appendColsToLpVectors(lp, XnumNewCol, local_colCost,
                                                local_colLower, local_colUpper),
                          return_status, "appendColsToLpVectors");
  if (return_status == HighsStatus::Error) return return_status;

  if (valid_simplex_lp) {
    // Append the columns to the Simplex LP vectors and matrix
    return_status = interpretCallStatus(
        appendColsToLpVectors(simplex_lp, XnumNewCol, local_colCost,
                              local_colLower, local_colUpper),
        return_status, "appendColsToLpVectors");
    if (return_status == HighsStatus::Error) return return_status;
  }

  // Now consider scaling. First resize the scaling factors and
  // initialise the new components
  scale.col_.resize(newNumCol);
  for (int col = 0; col < XnumNewCol; col++)
    scale.col_[lp.numCol_ + col] = 1.0;

  // Now consider any new matrix columns
  if (XnumNewNZ) {
    // There are nonzeros, so take a copy of the matrix that can be
    // normalised
    int local_num_new_nz = XnumNewNZ;
    std::vector<int> local_Astart{XAstart, XAstart + XnumNewCol};
    std::vector<int> local_Aindex{XAindex, XAindex + XnumNewNZ};
    std::vector<double> local_Avalue{XAvalue, XAvalue + XnumNewNZ};
    local_Astart.resize(XnumNewCol + 1);
    local_Astart[XnumNewCol] = XnumNewNZ;
    // Assess the matrix columns
    return_status = interpretCallStatus(
        assessMatrix(options, lp.numRow_, XnumNewCol, local_Astart,
                     local_Aindex, local_Avalue, options.small_matrix_value,
                     options.large_matrix_value),
        return_status, "assessMatrix");
    if (return_status == HighsStatus::Error) return return_status;
    local_num_new_nz = local_Astart[XnumNewCol];
    // Append the columns to the LP matrix
    return_status = interpretCallStatus(
        appendColsToLpMatrix(lp, XnumNewCol, local_num_new_nz, &local_Astart[0],
                             &local_Aindex[0], &local_Avalue[0]),
        return_status, "appendColsToLpMatrix");
    if (return_status == HighsStatus::Error) return return_status;
    if (valid_simplex_lp) {
      if (scaled_simplex_lp) {
        // Apply the row scaling to the new columns
        applyRowScalingToMatrix(scale.row_, XnumNewCol, local_Astart,
                                local_Aindex, local_Avalue);
        // Determine and apply the column scaling for the new columns
        colScaleMatrix(options.allowed_simplex_matrix_scale_factor,
                       &scale.col_[lp.numCol_], XnumNewCol,
                       local_Astart, local_Aindex, local_Avalue);
      }
      // Append the columns to the Simplex LP matrix
      return_status = interpretCallStatus(
          appendColsToLpMatrix(simplex_lp, XnumNewCol, local_num_new_nz,
                               &local_Astart[0], &local_Aindex[0],
                               &local_Avalue[0]),
          return_status, "appendColsToLpMatrix");
      if (return_status == HighsStatus::Error) return return_status;
      if (scaled_simplex_lp) {
        // Apply the column scaling to the costs and bounds
        HighsIndexCollection scaling_index_collection;
        scaling_index_collection.dimension_ = newNumCol;
        scaling_index_collection.is_interval_ = true;
        scaling_index_collection.from_ = lp.numCol_;
        scaling_index_collection.to_ = newNumCol - 1;
        return_status = interpretCallStatus(
            applyScalingToLpColCost(options, simplex_lp, scale.col_,
                                    scaling_index_collection),
            return_status, "applyScalingToLpColCost");
        if (return_status == HighsStatus::Error) return return_status;
        return_status = interpretCallStatus(
            applyScalingToLpColBounds(options, simplex_lp, scale.col_,
                                      scaling_index_collection),
            return_status, "applyScalingToLpColBounds");
        if (return_status == HighsStatus::Error) return return_status;
      }
    }
  } else {
    // There are no nonzeros, so XAstart/XAindex/XAvalue may be null. Have to
    // set up starts for empty columns
    assert(XnumNewCol > 0);
    appendColsToLpMatrix(lp, XnumNewCol, 0, NULL, NULL, NULL);
    if (valid_simplex_lp) {
      appendColsToLpMatrix(simplex_lp, XnumNewCol, 0, NULL, NULL, NULL);
      // Should be extendSimplexLpRandomVectors here
    }
  }
  // Update the basis correponding to new nonbasic columns
  if (valid_basis) appendNonbasicColsToBasis(lp, basis, XnumNewCol);
  if (valid_simplex_basis)
    appendNonbasicColsToBasis(simplex_lp, simplex_basis, XnumNewCol);

  // Deduce the consequences of adding new columns
  highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
  highs_model_object.unscaled_model_status_ =
      highs_model_object.scaled_model_status_;
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_COLS);

  // Increase the number of columns in the LPs
  lp.numCol_ += XnumNewCol;
  if (valid_simplex_lp) {
    simplex_lp.numCol_ += XnumNewCol;
    ekk_instance.initialiseSimplexLpRandomVectors();
  }

  return return_status;
}

HighsStatus Highs::addRowsInterface(int XnumNewRow,
				    const double* XrowLower,
				    const double* XrowUpper,
				    int XnumNewNZ, const int* XARstart,
				    const int* XARindex,
				    const double* XARvalue) {
  // addRows is fundamentally different from addCols, since the new
  // matrix data are held row-wise, so we have to insert data into the
  // column-wise matrix of the LP.
  HighsModelObject& highs_model_object = hmos_[0];
  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsStatus return_status = HighsStatus::OK;
  HighsOptions& options = highs_model_object.options_;
  if (XnumNewRow < 0) return HighsStatus::Error;
  if (XnumNewNZ < 0) return HighsStatus::Error;
  if (XnumNewRow == 0) return HighsStatus::OK;
  if (XnumNewRow > 0)
    if (isRowDataNull(options, XrowLower, XrowUpper)) return HighsStatus::Error;
  if (XnumNewNZ > 0)
    if (isMatrixDataNull(options, XARstart, XARindex, XARvalue))
      return HighsStatus::Error;

  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsScale& scale = highs_model_object.scale_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;
  HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  SimplexBasis& simplex_basis = ekk_instance.simplex_basis_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool& valid_basis = basis.valid_;
  bool& valid_simplex_lp = simplex_lp_status.valid;
  bool& valid_simplex_basis = simplex_lp_status.has_basis;
  bool& scaled_simplex_lp = scale.is_scaled_;

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
    assert(!scaled_simplex_lp);
  }
#endif

  HighsIndexCollection index_collection;
  index_collection.dimension_ = XnumNewRow;
  index_collection.is_interval_ = true;
  index_collection.from_ = 0;
  index_collection.to_ = XnumNewRow - 1;
  // Take a copy of the bounds that can be normalised
  std::vector<double> local_rowLower{XrowLower, XrowLower + XnumNewRow};
  std::vector<double> local_rowUpper{XrowUpper, XrowUpper + XnumNewRow};

  return_status = interpretCallStatus(
      assessBounds(options, "Row", lp.numRow_, index_collection, local_rowLower,
                   local_rowUpper, options.infinite_bound),
      return_status, "assessBounds");
  if (return_status == HighsStatus::Error) return return_status;

  // Append the rows to the LP vectors
  return_status = interpretCallStatus(
      appendRowsToLpVectors(lp, XnumNewRow, local_rowLower, local_rowUpper),
      return_status, "appendRowsToLpVectors");
  if (return_status == HighsStatus::Error) return return_status;

  if (valid_simplex_lp) {
    // Append the rows to the Simplex LP vectors
    return_status = interpretCallStatus(
        appendRowsToLpVectors(simplex_lp, XnumNewRow, local_rowLower,
                              local_rowUpper),
        return_status, "appendRowsToLpVectors");
    if (return_status == HighsStatus::Error) return return_status;
  }

  // Now consider scaling. First resize the scaling factors and
  // initialise the new components
  scale.row_.resize(newNumRow);
  for (int row = 0; row < XnumNewRow; row++)
    scale.row_[lp.numRow_ + row] = 1.0;

  // Now consider any new matrix rows
  if (XnumNewNZ) {
    // There are nonzeros, so take a copy of the matrix that can be
    // normalised
    int local_num_new_nz = XnumNewNZ;
    std::vector<int> local_ARstart{XARstart, XARstart + XnumNewRow};
    std::vector<int> local_ARindex{XARindex, XARindex + XnumNewNZ};
    std::vector<double> local_ARvalue{XARvalue, XARvalue + XnumNewNZ};
    local_ARstart.resize(XnumNewRow + 1);
    local_ARstart[XnumNewRow] = XnumNewNZ;
    // Assess the matrix columns
    return_status = interpretCallStatus(
        assessMatrix(options, lp.numCol_, XnumNewRow, local_ARstart,
                     local_ARindex, local_ARvalue, options.small_matrix_value,
                     options.large_matrix_value),
        return_status, "assessMatrix");
    if (return_status == HighsStatus::Error) return return_status;
    local_num_new_nz = local_ARstart[XnumNewRow];
    // Append the rows to LP matrix
    return_status = interpretCallStatus(
        appendRowsToLpMatrix(lp, XnumNewRow, local_num_new_nz,
                             &local_ARstart[0], &local_ARindex[0],
                             &local_ARvalue[0]),
        return_status, "appendRowsToLpMatrix");
    if (return_status == HighsStatus::Error) return return_status;
    if (valid_simplex_lp) {
      if (scaled_simplex_lp) {
        // Apply the column scaling to the new rows
        applyRowScalingToMatrix(scale.col_, XnumNewRow, local_ARstart,
                                local_ARindex, local_ARvalue);
        // Determine and apply the row scaling for the new rows. Using
        // colScaleMatrix to take the row-wise matrix and then treat
        // it col-wise
        colScaleMatrix(options.allowed_simplex_matrix_scale_factor,
                       &scale.row_[lp.numRow_], XnumNewRow,
                       local_ARstart, local_ARindex, local_ARvalue);
      }
      // Append the rows to the Simplex LP matrix
      return_status = interpretCallStatus(
          appendRowsToLpMatrix(simplex_lp, XnumNewRow, local_num_new_nz,
                               &local_ARstart[0], &local_ARindex[0],
                               &local_ARvalue[0]),
          return_status, "appendRowsToLpMatrix");
      if (return_status == HighsStatus::Error) return return_status;
      // Should be extendSimplexLpRandomVectors
      if (scaled_simplex_lp) {
        // Apply the row scaling to the bounds
        HighsIndexCollection scaling_index_collection;
        scaling_index_collection.dimension_ = newNumRow;
        scaling_index_collection.is_interval_ = true;
        scaling_index_collection.from_ = lp.numRow_;
        scaling_index_collection.to_ = newNumRow - 1;
        return_status = interpretCallStatus(
            applyScalingToLpRowBounds(options, simplex_lp, scale.row_,
                                      scaling_index_collection),
            return_status, "applyScalingToLpRowBounds");
        if (return_status == HighsStatus::Error) return return_status;
      }
    }
  }
  // Update the basis correponding to new basic rows
  if (valid_basis) appendBasicRowsToBasis(lp, basis, XnumNewRow);
  if (valid_simplex_basis)
    appendBasicRowsToBasis(simplex_lp, simplex_basis, XnumNewRow);

  // Deduce the consequences of adding new rows
  highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
  highs_model_object.unscaled_model_status_ =
      highs_model_object.scaled_model_status_;
  updateSimplexLpStatus(simplex_lp_status, LpAction::NEW_ROWS);

  // Increase the number of rows in the LPs
  lp.numRow_ += XnumNewRow;
  if (valid_simplex_lp) {
    simplex_lp.numRow_ += XnumNewRow;
    ekk_instance.initialiseSimplexLpRandomVectors();
  }

  return return_status;
}

HighsStatus Highs::deleteColsInterface(HighsIndexCollection& index_collection) {
  HighsModelObject& highs_model_object = hmos_[0];
  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsOptions& options = highs_model_object.options_;
  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;
  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool& valid_simplex_lp = simplex_lp_status.valid;
  // Keep a copy of the original number of columns to check whether
  // any columns have been removed, and if there is mask to be updated
  int original_num_col = lp.numCol_;

  HighsStatus return_status;
  return_status = deleteLpCols(options, lp, index_collection);
  if (return_status != HighsStatus::OK) return return_status;
  assert(lp.numCol_ <= original_num_col);
  if (lp.numCol_ < original_num_col) {
    // Nontrivial deletion so reset the model_status and invalidate
    // the Highs basis
    highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
    highs_model_object.unscaled_model_status_ =
        highs_model_object.scaled_model_status_;
    basis.valid_ = false;
  }
  return_status = interpretCallStatus(
      deleteScale(options, highs_model_object.scale_.col_, index_collection),
      return_status, "deleteScale");
  if (return_status == HighsStatus::Error) return return_status;
  highs_model_object.scale_.col_.resize(lp.numCol_);
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = ekk_instance.simplex_lp_;
    return_status = deleteLpCols(options, simplex_lp, index_collection);
    if (return_status != HighsStatus::OK) return return_status;
    assert(simplex_lp.numCol_ <= original_num_col);
    if (simplex_lp.numCol_ < original_num_col) {
      // Nontrivial deletion so initialise the random vectors and all
      // data relating to the simplex basis
      ekk_instance.initialiseSimplexLpRandomVectors();
      invalidateSimplexLpBasis(simplex_lp_status);
    }
  }
  if (index_collection.is_mask_) {
    // Set the mask values to indicate the new index value of the
    // remaining columns
    int new_col = 0;
    for (int col = 0; col < original_num_col; col++) {
      if (!index_collection.mask_[col]) {
        index_collection.mask_[col] = new_col;
        new_col++;
      } else {
        index_collection.mask_[col] = -1;
      }
    }
    assert(new_col == lp.numCol_);
  }
  return HighsStatus::OK;
}

HighsStatus Highs::deleteRowsInterface(HighsIndexCollection& index_collection) {
  HighsModelObject& highs_model_object = hmos_[0];
  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsOptions& options = highs_model_object.options_;
  HighsLp& lp = highs_model_object.lp_;
  HighsBasis& basis = highs_model_object.basis_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;

  // Query: should simplex_lp_status.valid be simplex_lp_status.valid_?
  bool& valid_simplex_lp = simplex_lp_status.valid;
  // Keep a copy of the original number of rows to check whether
  // any rows have been removed, and if there is mask to be updated
  int original_num_row = lp.numRow_;

  HighsStatus return_status;
  return_status = deleteLpRows(options, lp, index_collection);
  if (return_status != HighsStatus::OK) return return_status;
  assert(lp.numRow_ <= original_num_row);
  if (lp.numRow_ < original_num_row) {
    // Nontrivial deletion so reset the model_status and invalidate
    // the Highs basis
    highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
    highs_model_object.unscaled_model_status_ =
        highs_model_object.scaled_model_status_;
    basis.valid_ = false;
  }
  return_status = interpretCallStatus(
      deleteScale(options, highs_model_object.scale_.row_, index_collection),
      return_status, "deleteScale");
  if (return_status == HighsStatus::Error) return return_status;

  highs_model_object.scale_.row_.resize(lp.numRow_);
  if (valid_simplex_lp) {
    HighsLp& simplex_lp = ekk_instance.simplex_lp_;
    return_status = deleteLpRows(options, simplex_lp, index_collection);
    if (return_status != HighsStatus::OK) return return_status;
    assert(simplex_lp.numRow_ <= original_num_row);
    if (simplex_lp.numRow_ < original_num_row) {
      // Nontrivial deletion so initialise the random vectors and all
      // data relating to the simplex basis
      ekk_instance.initialiseSimplexLpRandomVectors();
      invalidateSimplexLpBasis(simplex_lp_status);
    }
  }
  if (index_collection.is_mask_) {
    int new_row = 0;
    for (int row = 0; row < original_num_row; row++) {
      if (!index_collection.mask_[row]) {
        index_collection.mask_[row] = new_row;
        new_row++;
      } else {
        index_collection.mask_[row] = -1;
      }
    }
    assert(new_row == lp.numRow_);
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getColsInterface(
    const HighsIndexCollection& index_collection, int& num_col,
    double* col_cost, double* col_lower, double* col_upper, int& num_nz,
    int* col_matrix_start, int* col_matrix_index, double* col_matrix_value) {
  HighsModelObject& highs_model_object = hmos_[0];
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp& lp = highs_model_object.lp_;
  HighsOptions& options = highs_model_object.options_;
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k < 0 || to_k > lp.numCol_) {
    call_status = HighsStatus::Error;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  if (from_k > to_k) {
    call_status = HighsStatus::Error;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  int out_from_col;
  int out_to_col;
  int in_from_col;
  int in_to_col = -1;
  int current_set_entry = 0;
  int col_dim = lp.numCol_;

  num_col = 0;
  num_nz = 0;
  for (int k = from_k; k <= to_k; k++) {
    updateIndexCollectionOutInIndex(index_collection, out_from_col, out_to_col,
                                    in_from_col, in_to_col, current_set_entry);
    assert(out_to_col < col_dim);
    assert(in_to_col < col_dim);
    for (int col = out_from_col; col <= out_to_col; col++) {
      if (col_cost != NULL) col_cost[num_col] = lp.colCost_[col];
      if (col_lower != NULL) col_lower[num_col] = lp.colLower_[col];
      if (col_upper != NULL) col_upper[num_col] = lp.colUpper_[col];
      if (col_matrix_start != NULL)
        col_matrix_start[num_col] =
            num_nz + lp.Astart_[col] - lp.Astart_[out_from_col];
      num_col++;
    }
    if (col_matrix_index != NULL || col_matrix_value != NULL) {
      for (int el = lp.Astart_[out_from_col]; el < lp.Astart_[out_to_col + 1];
           el++) {
        if (col_matrix_index != NULL) col_matrix_index[num_nz] = lp.Aindex_[el];
        if (col_matrix_value != NULL) col_matrix_value[num_nz] = lp.Avalue_[el];
        num_nz++;
      }
    }
    if (out_to_col == col_dim - 1 || in_to_col == col_dim - 1) break;
  }
  return HighsStatus::OK;
}

HighsStatus Highs::getRowsInterface(
    const HighsIndexCollection& index_collection, int& num_row,
    double* row_lower, double* row_upper, int& num_nz, int* row_matrix_start,
    int* row_matrix_index, double* row_matrix_value) {
  HighsModelObject& highs_model_object = hmos_[0];
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsLp& lp = highs_model_object.lp_;
  HighsOptions& options = highs_model_object.options_;
  if (!assessIndexCollection(options, index_collection))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "assessIndexCollection");
  int from_k;
  int to_k;
  if (!limitsForIndexCollection(options, index_collection, from_k, to_k))
    return interpretCallStatus(HighsStatus::Error, return_status,
                               "limitsForIndexCollection");
  if (from_k < 0 || to_k > lp.numRow_) {
    call_status = HighsStatus::Error;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  num_row = 0;
  num_nz = 0;
  if (from_k > to_k) {
    call_status = HighsStatus::Error;
    return_status = interpretCallStatus(call_status, return_status, "getCols");
    return return_status;
  }
  // "Out" means not in the set to be extrated
  // "In" means in the set to be extrated
  int out_from_row;
  int out_to_row;
  int in_from_row;
  int in_to_row = -1;
  int current_set_entry = 0;
  int row_dim = lp.numRow_;

  // Set up a row mask so that entries to be got from the column-wise
  // matrix can be identified and have their correct row index.
  vector<int> new_index;
  new_index.resize(lp.numRow_);

  if (!index_collection.is_mask_) {
    out_to_row = -1;
    current_set_entry = 0;
    for (int k = from_k; k <= to_k; k++) {
      updateIndexCollectionOutInIndex(index_collection, in_from_row, in_to_row,
                                      out_from_row, out_to_row,
                                      current_set_entry);
      if (k == from_k) {
        // Account for any initial rows not being extracted
        for (int row = 0; row < in_from_row; row++) {
          new_index[row] = -1;
        }
      }
      for (int row = in_from_row; row <= in_to_row; row++) {
        new_index[row] = num_row;
        num_row++;
      }
      for (int row = out_from_row; row <= out_to_row; row++) {
        new_index[row] = -1;
      }
      if (out_to_row >= row_dim - 1) break;
    }
  } else {
    for (int row = 0; row < lp.numRow_; row++) {
      if (index_collection.mask_[row]) {
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
  vector<int> row_matrix_length;
  row_matrix_length.resize(num_row);

  for (int row = 0; row < lp.numRow_; row++) {
    int new_row = new_index[row];
    if (new_row >= 0) {
      assert(new_row < num_row);
      if (row_lower != NULL) row_lower[new_row] = lp.rowLower_[row];
      if (row_upper != NULL) row_upper[new_row] = lp.rowUpper_[row];
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

  if (row_matrix_start == NULL) {
    // If the matrix start vector is null then don't get values of
    // indices, otherwise both are meaningless
    if (row_matrix_index != NULL || row_matrix_value != NULL) {
      HighsLogMessage(highs_model_object.options_.logfile,
                      HighsMessageType::ERROR,
                      "Cannot supply meaningful row matrix indices/values with "
                      "null starts");
      return HighsStatus::Error;
    }
  } else {
    row_matrix_start[0] = 0;
    for (int row = 0; row < num_row - 1; row++) {
      row_matrix_start[row + 1] =
          row_matrix_start[row] + row_matrix_length[row];
    }

    // Fill the row-wise matrix with indices and values
    for (int col = 0; col < lp.numCol_; col++) {
      for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
        int row = lp.Aindex_[el];
        int new_row = new_index[row];
        if (new_row >= 0) {
          int row_el = row_matrix_start[new_row];
          if (row_matrix_index != NULL) row_matrix_index[row_el] = col;
          if (row_matrix_value != NULL)
            row_matrix_value[row_el] = lp.Avalue_[el];
          row_matrix_start[new_row]++;
        }
      }
    }
    // Restore the starts of the row-wise matrix and count the number of
    // nonzeros in it
    num_nz = 0;
    row_matrix_start[0] = 0;
    for (int row = 0; row < num_row - 1; row++) {
      row_matrix_start[row + 1] =
          row_matrix_start[row] + row_matrix_length[row];
      num_nz += row_matrix_length[row];
    }
    num_nz += row_matrix_length[num_row - 1];
  }
  return HighsStatus::OK;
}

