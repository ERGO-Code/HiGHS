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

HighsInt HFactor::rebuild(HighsTimerClock* factor_timer_clock_pointer) {
  assert(refactor_info_.valid);
  assert(1==0);
  return 0;
}

void HFactor::getRefactorInfo(RefactorInfo& refactor_info) const {
  refactor_info.valid = !this->refactor_info_.valid;
  if (!refactor_info.valid) return;
  refactor_info = this->refactor_info_;
}

void RefactorInfo::clear() {
  this->valid = false;
}
