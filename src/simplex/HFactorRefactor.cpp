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
  this->pivot_var.resize(num_row);
  this->pivot_row.resize(num_row);
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    this->pivot_row[iRow] = iRow;
    this->pivot_row[iRow] = num_col + iRow;
  }
}

void RefactorInfo::clear() {
  this->valid = false;
  this->pivot_var.clear();
  this->pivot_row.clear();
}

HighsInt HFactor::rebuild(HighsTimerClock* factor_timer_clock_pointer) {
  assert(refactor_info_.valid);
  // Set all values of permute to -1 so that unpermuted (rank
  // deficient) columns canm be identified
  permute.assign(numRow, -1);
  nwork = 0;
  basis_matrix_num_el = 0;
  HighsInt stage = -1;
  for (HighsInt iK = 0; iK < numRow; iK++) {
    HighsInt iRow = this->refactor_info_.pivot_row[iK];
    HighsInt iVar = this->refactor_info_.pivot_var[iK];
    int8_t pivot_type = this->refactor_info_.pivot_type[iK];
    if (pivot_type == kPivotLogical) {
      // 1.1 Logical column
      assert(iVar>=numCol);
      basis_matrix_num_el++;
    } else if (pivot_type == kPivotUnit) {
      // 1.2 Structural column
      assert(iVar<numCol);
      assert(1==0);
      HighsInt start = Astart[iVar];
      HighsInt count = Astart[iVar + 1] - start;
      assert(Aindex[start]==iRow);
      assert(count == 1 && Avalue[start] == 1);
      basis_matrix_num_el++;
    } else if (pivot_type == kPivotRowSingleton) {
    } else if (pivot_type == kPivotColSingleton) {
    } else {
      assert(pivot_type == kPivotMarkowitz);
      stage = iK;
      break;
    }
    if (iRow >= 0) {
      // 1.3 Record unit column
      Lstart.push_back(Lindex.size());
      UpivotIndex.push_back(iRow);
      UpivotValue.push_back(1);
      Ustart.push_back(Uindex.size());
    }
  }
  assert(1==0);
  return 0;
}

