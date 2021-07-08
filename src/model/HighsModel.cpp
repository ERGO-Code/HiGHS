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
/**@file lp_data/HighsModel.cpp
 * @brief
 */
#include "model/HighsModel.h"

#include <cassert>

bool HighsModel::isQp() const {
  HighsInt dim = this->hessian_.dim_;
  HighsInt num_col = this->lp_.numCol_;
  if (dim) {
    // If there's a Hessian then it's a QP, but its dimension should
    // be the same as the number of LP columns
    assert(dim == num_col);
    return true;
  }
  return false;
}

void HighsModel::clear() {
  this->lp_.clear();
  this->hessian_.clear();
}

double HighsModel::objectiveValue(const std::vector<double>& solution) const {
  assert((int)solution.size() >= this->lp_.numCol_);
  double objective_function_value = this->lp_.offset_;
  for (HighsInt iCol = 0; iCol < this->lp_.numCol_; iCol++)
    objective_function_value += this->lp_.colCost_[iCol] * solution[iCol];
  if (this->hessian_.dim_ == 0) return objective_function_value;
  for (HighsInt iCol = 0; iCol < this->hessian_.dim_; iCol++) {
    HighsInt iEl = this->hessian_.q_start_[iCol];
    assert(this->hessian_.q_index_[iEl] == iCol);
    objective_function_value +=
        0.5 * solution[iCol] * this->hessian_.q_value_[iEl] * solution[iCol];
    for (HighsInt iEl = this->hessian_.q_start_[iCol] + 1;
         iEl < this->hessian_.q_start_[iCol + 1]; iEl++)
      objective_function_value += solution[iCol] *
                                  this->hessian_.q_value_[iEl] *
                                  solution[this->hessian_.q_index_[iEl]];
  }
  return objective_function_value;
}
