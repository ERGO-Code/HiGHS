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
  equal = this->colCost_ == lp.colCost_ && equal;
  equal = this->colUpper_ == lp.colUpper_ && equal;
  equal = this->colLower_ == lp.colLower_ && equal;
  equal = this->rowUpper_ == lp.rowUpper_ && equal;
  equal = this->rowLower_ == lp.rowLower_ && equal;
  equal = this->a_start_ == lp.a_start_ && equal;
  equal = this->a_index_ == lp.a_index_ && equal;
  equal = this->a_value_ == lp.a_value_ && equal;
  equal = this->format_ == lp.format_ && equal;
  return equal;
}

double HighsLp::objectiveValue(const std::vector<double>& solution) const {
  assert((int)solution.size() >= this->num_col_);
  double objective_function_value = this->offset_;
  for (HighsInt iCol = 0; iCol < this->num_col_; iCol++)
    objective_function_value += this->colCost_[iCol] * solution[iCol];
  return objective_function_value;
}

void HighsLp::clear() {
  this->num_col_ = 0;
  this->num_row_ = 0;

  this->a_start_.clear();
  this->a_index_.clear();
  this->a_value_.clear();
  this->colCost_.clear();
  this->colLower_.clear();
  this->colUpper_.clear();
  this->rowLower_.clear();
  this->rowUpper_.clear();

  this->sense_ = ObjSense::kMinimize;
  this->offset_ = 0;
  this->format_ = MatrixFormat::kNone;

  this->model_name_ = "";

  this->col_names_.clear();
  this->row_names_.clear();

  this->integrality_.clear();
}
