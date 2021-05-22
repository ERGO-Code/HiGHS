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

bool HighsModel::isQp() {
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

