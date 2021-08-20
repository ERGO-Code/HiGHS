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
  this->pivot_type.resize(num_row);
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    this->pivot_var[iRow] = num_col + iRow;
    this->pivot_row[iRow] = iRow;
    this->pivot_type[iRow] = kPivotLogical;
  }
}

void RefactorInfo::clear() {
  this->valid = false;
  this->pivot_var.clear();
  this->pivot_row.clear();
  this->pivot_type.clear();
}

HighsInt HFactor::rebuild(HighsTimerClock* factor_timer_clock_pointer) {
  assert(refactor_info_.valid);
  /**
   * 0. Clear L and U factor
   */
  LuClear();

  // Set all values of permute to -1 so that unpermuted (rank
  // deficient) columns canm be identified
  //  permute.assign(numRow, -1);
  nwork = 0;
  basis_matrix_num_el = 0;
  HighsInt stage = numRow;
  vector<bool> has_pivot;
  has_pivot.assign(numRow, false);
  printf("\nRefactor\n");
  for (HighsInt iK = 0; iK < numRow; iK++) {
    HighsInt iRow = this->refactor_info_.pivot_row[iK];
    HighsInt iVar = this->refactor_info_.pivot_var[iK];
    int8_t pivot_type = this->refactor_info_.pivot_type[iK];
    assert(!has_pivot[iRow]);
    if (pivot_type == kPivotLogical) {
      //
      // 1.1 Logical column
      printf("Stage %d: Logical\n", (int)iK);
      assert(iVar>=numCol);
      basis_matrix_num_el++;
    } else if (pivot_type == kPivotUnit) {
      //
      // 1.2 (Structural) unit column 
      printf("Stage %d: Unit\n", (int)iK);
      assert(iVar<numCol);
      assert(1==0);
      HighsInt start = Astart[iVar];
      HighsInt count = Astart[iVar + 1] - start;
      assert(Aindex[start]==iRow);
      assert(count == 1 && Avalue[start] == 1);
      basis_matrix_num_el++;
    } else if (pivot_type == kPivotRowSingleton || pivot_type == kPivotColSingleton) {
      //
      // Row or column singleton
      assert(iVar<numCol);
      const HighsInt start = Astart[iVar];
      const HighsInt end = Astart[iVar + 1];
      // Find where the pivot is
      HighsInt pivot_k = -1;
      for (HighsInt k = start; k < end; k++) {
	if (Aindex[k] == iRow) {
	  pivot_k = k;
          break;
	}
      }
      assert(pivot_k>=0);
      if (pivot_type == kPivotRowSingleton) {
	//
	// 2.2 Deal with row singleton
	const double pivotX = 1 / Avalue[pivot_k];
	printf("Stage %d: Row singleton (%4d, %g)\n", (int)iK, (int)pivot_k, pivotX);
	for (HighsInt section = 0; section < 2; section++) {
	  HighsInt p0 = section == 0 ? start : pivot_k + 1;
	  HighsInt p1 = section == 0 ? pivot_k : end;
	  for (HighsInt k = p0; k < p1; k++) {
	    HighsInt local_iRow = Aindex[k];
	    if (!has_pivot[local_iRow]) {
	      printf("Row singleton: L En (%4d, %11.4g)\n", (int)local_iRow, Avalue[k] * pivotX);
	      Lindex.push_back(local_iRow);
	      Lvalue.push_back(Avalue[k] * pivotX);
	    } else {
	      printf("Row singleton: U En (%4d, %11.4g)\n", (int)local_iRow, Avalue[k]);
	      Uindex.push_back(local_iRow);
	      Uvalue.push_back(Avalue[k]);
	    }
	  }
	}
	Lstart.push_back(Lindex.size());
	printf("Row singleton: U Pv (%4d, %11.4g)\n", (int)iRow, Avalue[pivot_k]);
	UpivotIndex.push_back(iRow);
	UpivotValue.push_back(Avalue[pivot_k]);
	Ustart.push_back(Uindex.size());
      } else {
	//
        // 2.3 Deal with column singleton
	printf("Stage %d: Col singleton\n", (int)iK);
        for (HighsInt k = start; k < pivot_k; k++) {
	  printf("Col singleton: U En (%4d, %11.4g)\n", (int)Aindex[k], Avalue[k]);
          Uindex.push_back(Aindex[k]);
          Uvalue.push_back(Avalue[k]);
        }
        for (HighsInt k = pivot_k + 1; k < end; k++) {
	  printf("Col singleton: U En (%4d, %11.4g)\n", (int)Aindex[k], Avalue[k]);
          Uindex.push_back(Aindex[k]);
          Uvalue.push_back(Avalue[k]);
        }
        Lstart.push_back(Lindex.size());
	printf("Col singleton: U Pv (%4d, %11.4g)\n", (int)iRow, Avalue[pivot_k]);
        UpivotIndex.push_back(iRow);
        UpivotValue.push_back(Avalue[pivot_k]);
        Ustart.push_back(Uindex.size());
      }
    } else {
      assert(pivot_type == kPivotMarkowitz);
      stage = iK;
      break;
    }
    has_pivot[iRow] = true;
    if (pivot_type == kPivotLogical || pivot_type == kPivotUnit) {
      // 1.3 Record unit column
      //      permute[iCol] = iRow;
      Lstart.push_back(Lindex.size());
      UpivotIndex.push_back(iRow);
      UpivotValue.push_back(1);
      Ustart.push_back(Uindex.size());
    }
  }
  if (stage<numRow) {
    // Handle the remaining Markowitz pivots
    HVector column;
    column.setup(numRow);
    for (HighsInt iK = stage; iK < numRow; iK++) {
      HighsInt iRow = this->refactor_info_.pivot_row[iK];
      HighsInt iVar = this->refactor_info_.pivot_var[iK];
      int8_t pivot_type = this->refactor_info_.pivot_type[iK];
      assert(!has_pivot[iRow]);
      assert(pivot_type == kPivotMarkowitz);
      column.clear();
      const HighsInt start = Astart[iVar];
      const HighsInt end = Astart[iVar + 1];
      for (HighsInt iEl = start; iEl < end; iEl++) {
	HighsInt local_iRow = Aindex[iEl];
	if (has_pivot[local_iRow]) continue;
	column.index[column.count++] = local_iRow;
	column.array[local_iRow] = Avalue[iEl];
      }
      const double expected_density = 1.0;
      ftranL(column, expected_density, factor_timer_clock_pointer);
      assert(1==0);
    }
  }
  buildFinish();
  return 0;
}

