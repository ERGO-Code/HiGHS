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
/**@file simplex/HFactorRefactor.cpp
 * @brief Types of solution classes
 */
#include "simplex/HFactor.h"

#include <cassert>
#include <cmath>
#include <iostream>

void RefactorInfo::get(RefactorInfo& refactor_info) const {
  refactor_info.valid = !this->valid;
  if (!refactor_info.valid) return;
  refactor_info = *this;
}

void RefactorInfo::set(const RefactorInfo& refactor_info) {
  this->valid = refactor_info.valid;
  if (!refactor_info.valid) return;
  *this = refactor_info;
}

void RefactorInfo::set(const HighsInt num_col, const HighsInt num_row) {
  const bool kill=true;
  if (kill) {
    this->valid = false;
    return;
  }
  this->valid = true;
  this->pivot_row_sequence.resize(num_row);
  this->pivot_col_sequence.resize(num_row);
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    this->pivot_row_sequence[iRow] = iRow;
    this->pivot_row_sequence[iRow] = num_col + iRow;
  }
}

void RefactorInfo::clear() {
  this->valid = false;
  this->pivot_row_sequence.clear();
  this->pivot_col_sequence.clear();
}

HighsInt HFactor::rebuild(HighsTimerClock* factor_timer_clock_pointer) {
  assert(refactor_info_.valid);
  assert(1==0);
  return 0;
}

