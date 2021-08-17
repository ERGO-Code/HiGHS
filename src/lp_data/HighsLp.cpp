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
/**@file lp_data/HighsLp.cpp
 * @brief
 */
#include "lp_data/HighsLp.h"

#include <cassert>

#include "util/HighsMatrixUtils.h"

bool HighsLp::isMip() const {
  HighsInt integrality_size = this->integrality_.size();
  if (integrality_size) {
    assert(integrality_size == this->num_col_);
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
      if (this->integrality_[iCol] != HighsVarType::kContinuous) return true;
  }
  return false;
}

bool HighsLp::operator==(const HighsLp& lp) {
  bool equal = equalButForNames(lp);
  equal = this->row_names_ == lp.row_names_ && equal;
  equal = this->col_names_ == lp.col_names_ && equal;
  return equal;
}

bool HighsLp::equalButForNames(const HighsLp& lp) const {
  bool equal = true;
  equal = this->num_col_ == lp.num_col_ && equal;
  equal = this->num_row_ == lp.num_row_ && equal;
  equal = this->sense_ == lp.sense_ && equal;
  equal = this->offset_ == lp.offset_ && equal;
  equal = this->model_name_ == lp.model_name_ && equal;
  equal = this->col_cost_ == lp.col_cost_ && equal;
  equal = this->col_upper_ == lp.col_upper_ && equal;
  equal = this->col_lower_ == lp.col_lower_ && equal;
  equal = this->row_upper_ == lp.row_upper_ && equal;
  equal = this->row_lower_ == lp.row_lower_ && equal;

  equal = this->a_matrix_ == lp.a_matrix_;

  equal = this->scale_.strategy == lp.scale_.strategy && equal;
  equal = this->scale_.has_scaling == lp.scale_.has_scaling && equal;
  equal = this->scale_.num_col == lp.scale_.num_col && equal;
  equal = this->scale_.num_row == lp.scale_.num_row && equal;
  equal = this->scale_.cost == lp.scale_.cost && equal;
  equal = this->scale_.col == lp.scale_.col && equal;
  equal = this->scale_.row == lp.scale_.row && equal;
  return equal;
}

double HighsLp::objectiveValue(const std::vector<double>& solution) const {
  assert((int)solution.size() >= this->num_col_);
  double objective_function_value = this->offset_;
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
    objective_function_value += this->col_cost_[iCol] * solution[iCol];
  return objective_function_value;
}

void HighsLp::setMatrixDimensions() {
  this->a_matrix_.num_col_ = this->num_col_;
  this->a_matrix_.num_row_ = this->num_row_;
}

void HighsLp::resetScale() {
  // Should allow user-supplied scale to be retained
  //  const bool dimensions_ok =
  //    this->scale_.num_col_ == this->num_col_ &&
  //    this->scale_.num_row_ == this->num_row_;
  this->clearScale();
}

bool HighsLp::dimensionsOk(std::string message) const {
  bool ok = true;
  const HighsInt num_col = this->num_col_;
  const HighsInt num_row = this->num_row_;
  if (!(num_col >= 0))
    printf("HighsLp::dimensionsOK (%s) fails on num_col >= 0\n",
           message.c_str());
  ok = num_col >= 0 && ok;
  if (!(num_row >= 0))
    printf("HighsLp::dimensionsOK (%s) fails on num_row >= 0\n",
           message.c_str());
  ok = num_row >= 0 && ok;
  if (!ok) return ok;

  HighsInt col_cost_size = this->col_cost_.size();
  HighsInt col_lower_size = this->col_lower_.size();
  HighsInt col_upper_size = this->col_upper_.size();
  HighsInt matrix_start_size = this->a_matrix_.start_.size();
  bool legal_col_cost_size = col_cost_size >= num_col;
  bool legal_col_lower_size = col_lower_size >= num_col;
  bool legal_col_upper_size = col_lower_size >= num_col;
  if (!legal_col_cost_size)
    printf("HighsLp::dimensionsOK (%s) fails on legal_col_cost_size\n",
           message.c_str());
  ok = legal_col_cost_size && ok;
  if (!legal_col_lower_size)
    printf("HighsLp::dimensionsOK (%s) fails on legal_col_lower_size\n",
           message.c_str());
  ok = legal_col_lower_size && ok;
  if (!legal_col_upper_size)
    printf("HighsLp::dimensionsOK (%s) fails on legal_col_upper_size\n",
           message.c_str());
  ok = legal_col_upper_size && ok;

  bool legal_format = this->a_matrix_.format_ == MatrixFormat::kColwise ||
                      this->a_matrix_.format_ == MatrixFormat::kRowwise;
  if (!legal_format)
    printf("HighsLp::dimensionsOK (%s) fails on legal_format\n",
           message.c_str());
  ok = legal_format && ok;
  HighsInt num_vec;
  if (this->a_matrix_.isColwise()) {
    num_vec = num_col;
  } else {
    num_vec = num_row;
  }
  const bool partitioned = false;
  vector<HighsInt> a_matrix_p_end;
  bool legal_matrix_dimensions =
      assessMatrixDimensions(num_vec, partitioned, this->a_matrix_.start_,
                             a_matrix_p_end, this->a_matrix_.index_,
                             this->a_matrix_.value_) == HighsStatus::kOk;
  if (!legal_matrix_dimensions)
    printf("HighsLp::dimensionsOK (%s) fails on legal_matrix_dimensions\n",
           message.c_str());
  ok = legal_matrix_dimensions && ok;

  HighsInt row_lower_size = this->row_lower_.size();
  HighsInt row_upper_size = this->row_upper_.size();
  bool legal_row_lower_size = row_lower_size >= num_row;
  bool legal_row_upper_size = row_lower_size >= num_row;
  if (!legal_row_lower_size)
    printf("HighsLp::dimensionsOK (%s) fails on legal_row_lower_size\n",
           message.c_str());
  ok = legal_row_lower_size && ok;
  if (!legal_row_upper_size)
    printf("HighsLp::dimensionsOK (%s) fails on legal_row_upper_size\n",
           message.c_str());
  ok = legal_row_upper_size && ok;

  bool legal_a_matrix_num_col = this->a_matrix_.num_col_ == num_col;
  bool legal_a_matrix_num_row = this->a_matrix_.num_row_ == num_row;
  if (!legal_a_matrix_num_col)
    printf("HighsLp::dimensionsOK (%s) fails on legal_a_matrix_num_col\n",
           message.c_str());
  ok = legal_a_matrix_num_col && ok;
  if (!legal_a_matrix_num_row)
    printf("HighsLp::dimensionsOK (%s) fails on legal_a_matrix_num_row\n",
           message.c_str());
  ok = legal_a_matrix_num_row && ok;

  HighsInt scale_strategy = (HighsInt)this->scale_.strategy;
  bool legal_scale_strategy = scale_strategy >= 0;
  if (!legal_scale_strategy)
    printf("HighsLp::dimensionsOK (%s) fails on legal_scale_strategy\n",
           message.c_str());
  ok = legal_scale_strategy && ok;
  HighsInt scale_row_size = (HighsInt)this->scale_.row.size();
  HighsInt scale_col_size = (HighsInt)this->scale_.col.size();
  bool legal_scale_num_col = false;
  bool legal_scale_num_row = false;
  bool legal_scale_row_size = false;
  bool legal_scale_col_size = false;
  if (this->scale_.has_scaling) {
    legal_scale_num_col = this->scale_.num_col == num_col;
    legal_scale_num_row = this->scale_.num_row == num_row;
    legal_scale_row_size = scale_row_size >= num_row;
    legal_scale_col_size = scale_col_size >= num_col;
  } else {
    legal_scale_num_col = this->scale_.num_col == 0;
    legal_scale_num_row = this->scale_.num_row == 0;
    legal_scale_row_size = scale_row_size == 0;
    legal_scale_col_size = scale_col_size == 0;
  }
  if (!legal_scale_num_col)
    printf("HighsLp::dimensionsOK (%s) fails on legal_scale_num_col\n",
           message.c_str());
  ok = legal_scale_num_col && ok;
  if (!legal_scale_num_row)
    printf("HighsLp::dimensionsOK (%s) fails on legal_scale_num_row\n",
           message.c_str());
  ok = legal_scale_num_row && ok;
  if (!legal_scale_row_size)
    printf("HighsLp::dimensionsOK (%s) fails on legal_scale_row_size\n",
           message.c_str());
  ok = legal_scale_row_size && ok;
  if (!legal_scale_col_size)
    printf("HighsLp::dimensionsOK (%s) fails on legal_scale_col_size\n",
           message.c_str());
  ok = legal_scale_col_size && ok;
  if (!ok) {
    printf("HighsLp::dimensionsOK (%s) fails \n", message.c_str());
  }

  return ok;
}

void HighsLp::setFormat(const MatrixFormat format) {
  this->a_matrix_.setFormat(format);
}

void HighsLp::clear() {
  this->num_col_ = 0;
  this->num_row_ = 0;

  this->col_cost_.clear();
  this->col_lower_.clear();
  this->col_upper_.clear();
  this->row_lower_.clear();
  this->row_upper_.clear();

  this->a_matrix_.clear();

  this->sense_ = ObjSense::kMinimize;
  this->offset_ = 0;

  this->model_name_ = "";

  this->col_names_.clear();
  this->row_names_.clear();

  this->integrality_.clear();

  this->clearScale();
  this->is_scaled_ = false;
  this->is_moved_ = false;
}

void HighsLp::clearScale() {
  this->scale_.strategy = kSimplexScaleStrategyOff;
  this->scale_.has_scaling = false;
  this->scale_.num_col = 0;
  this->scale_.num_row = 0;
  this->scale_.cost = 0;
  this->scale_.col.clear();
  this->scale_.row.clear();
}

void HighsLp::clearScaling() {
  this->unapplyScale();
  this->clearScale();
}

void HighsLp::applyScale() {
  // Ensure that any scaling is applied
  const HighsScale& scale = this->scale_;
  if (this->is_scaled_) {
    // Already scaled - so check that there is scaling and return
    assert(scale.has_scaling);
    return;
  }
  // No scaling currently applied
  this->is_scaled_ = false;
  if (scale.has_scaling) {
    // Apply the scaling to the bounds, costs and matrix, and record
    // that it has been applied
    for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
      this->col_lower_[iCol] /= scale.col[iCol];
      this->col_upper_[iCol] /= scale.col[iCol];
      this->col_cost_[iCol] *= scale.col[iCol];
    }
    for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
      this->row_lower_[iRow] *= scale.row[iRow];
      this->row_upper_[iRow] *= scale.row[iRow];
    }
    this->a_matrix_.applyScale(scale);
    this->is_scaled_ = true;
  }
}

void HighsLp::unapplyScale() {
  // Ensure that any scaling is not applied
  const HighsScale& scale = this->scale_;
  if (!this->is_scaled_) {
    // Not scaled so return
    return;
  }
  // Already scaled - so check that there is scaling
  assert(scale.has_scaling);
  // Unapply the scaling to the bounds, costs and matrix, and record
  // that it has been unapplied
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++) {
    this->col_lower_[iCol] *= scale.col[iCol];
    this->col_upper_[iCol] *= scale.col[iCol];
    this->col_cost_[iCol] /= scale.col[iCol];
  }
  for (HighsInt iRow = 0; iRow < this->num_row_; iRow++) {
    this->row_lower_[iRow] /= scale.row[iRow];
    this->row_upper_[iRow] /= scale.row[iRow];
  }
  this->a_matrix_.unapplyScale(scale);
  this->is_scaled_ = false;
}

void HighsLp::moveBackLpAndUnapplyScaling(HighsLp lp) {
  assert(this->is_moved_ == true);
  *this = std::move(lp);
  this->unapplyScale();
  assert(this->is_moved_ == false);
}
