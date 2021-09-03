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
#include <cassert>
#include <cmath>
#include <iostream>

#include "simplex/HFactor.h"

// std::max and std::min used in HFactor.h for local in-line
// functions, so HFactor.h has #include <algorithm>
using std::fabs;

void RefactorInfo::clear() {
  this->use = false;
  this->build_synthetic_tick = 0.0;
  this->pivot_var.clear();
  this->pivot_row.clear();
  this->pivot_type.clear();
}

HighsInt HFactor::rebuild(HighsTimerClock* factor_timer_clock_pointer) {
  const bool report_lu = false;
  // Check that the refactorzation information should be used
  assert(refactor_info_.use);
  /**
   * 0. Clear L and U factor
   */
  luClear();

  nwork = 0;
  basis_matrix_num_el = 0;
  HighsInt stage = numRow;
  HighsInt rank_deficiency = 0;
  vector<bool> has_pivot;
  has_pivot.assign(numRow, false);
  const bool report_unit = false;
  const bool report_singletons = false;
  const bool report_markowitz = false;
  const bool report_anything =
      report_unit || report_singletons || report_markowitz;
  if (report_anything) printf("\nRefactor\n");
  // Take build_synthetic_tick from the refactor info so that this
  // refactorization doesn't look unrealistically cheap.
  this->build_synthetic_tick = this->refactor_info_.build_synthetic_tick;
  for (HighsInt iK = 0; iK < numRow; iK++) {
    HighsInt iRow = this->refactor_info_.pivot_row[iK];
    HighsInt iVar = this->refactor_info_.pivot_var[iK];
    int8_t pivot_type = this->refactor_info_.pivot_type[iK];
    assert(!has_pivot[iRow]);

    if (pivot_type == kPivotLogical || pivot_type == kPivotUnit) {
      if (pivot_type == kPivotLogical) {
        //
        // 1.1 Logical column
        if (report_unit) printf("Stage %d: Logical\n", (int)iK);
        assert(iVar >= numCol);
        basis_matrix_num_el++;
      } else if (pivot_type == kPivotUnit) {
        //
        // 1.2 (Structural) unit column
        if (report_unit) printf("Stage %d: Unit\n", (int)iK);
        assert(iVar < numCol);
        HighsInt start = Astart[iVar];
        HighsInt count = Astart[iVar + 1] - start;
        assert(Aindex[start] == iRow);
        assert(count == 1 && Avalue[start] == 1);
        basis_matrix_num_el++;
      }
      // 1.3 Record unit column
      Lstart.push_back(Lindex.size());
      UpivotIndex.push_back(iRow);
      UpivotValue.push_back(1);
      Ustart.push_back(Uindex.size());
    } else if (pivot_type == kPivotRowSingleton ||
               pivot_type == kPivotColSingleton) {
      //
      // Row or column singleton
      assert(iVar < numCol);
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
      assert(pivot_k >= 0);
      // Check that the pivot isn't too small. Shouldn't happen since
      // this is refactorization
      double abs_pivot = std::fabs(Avalue[pivot_k]);
      assert(abs_pivot >= pivot_tolerance);
      if (abs_pivot < pivot_tolerance) {
        rank_deficiency = nwork + 1;
        return rank_deficiency;
      }
      if (pivot_type == kPivotRowSingleton) {
        //
        // 2.2 Deal with row singleton
        const double pivotX = 1 / Avalue[pivot_k];
        if (report_singletons)
          printf("Stage %d: Row singleton (%4d, %g)\n", (int)iK, (int)pivot_k,
                 pivotX);
        for (HighsInt section = 0; section < 2; section++) {
          HighsInt p0 = section == 0 ? start : pivot_k + 1;
          HighsInt p1 = section == 0 ? pivot_k : end;
          for (HighsInt k = p0; k < p1; k++) {
            HighsInt local_iRow = Aindex[k];
            if (!has_pivot[local_iRow]) {
              if (report_singletons)
                printf("Row singleton: L En (%4d, %11.4g)\n", (int)local_iRow,
                       Avalue[k] * pivotX);
              Lindex.push_back(local_iRow);
              Lvalue.push_back(Avalue[k] * pivotX);
            } else {
              if (report_singletons)
                printf("Row singleton: U En (%4d, %11.4g)\n", (int)local_iRow,
                       Avalue[k]);
              Uindex.push_back(local_iRow);
              Uvalue.push_back(Avalue[k]);
            }
          }
        }
        Lstart.push_back(Lindex.size());
        if (report_singletons)
          printf("Row singleton: U Pv (%4d, %11.4g)\n", (int)iRow,
                 Avalue[pivot_k]);
        UpivotIndex.push_back(iRow);
        UpivotValue.push_back(Avalue[pivot_k]);
        Ustart.push_back(Uindex.size());
      } else {
        //
        // 2.3 Deal with column singleton
        if (report_singletons) printf("Stage %d: Col singleton\n", (int)iK);
        for (HighsInt k = start; k < pivot_k; k++) {
          if (report_singletons)
            printf("Col singleton: U En (%4d, %11.4g)\n", (int)Aindex[k],
                   Avalue[k]);
          Uindex.push_back(Aindex[k]);
          Uvalue.push_back(Avalue[k]);
        }
        for (HighsInt k = pivot_k + 1; k < end; k++) {
          if (report_singletons)
            printf("Col singleton: U En (%4d, %11.4g)\n", (int)Aindex[k],
                   Avalue[k]);
          Uindex.push_back(Aindex[k]);
          Uvalue.push_back(Avalue[k]);
        }
        Lstart.push_back(Lindex.size());
        if (report_singletons)
          printf("Col singleton: U Pv (%4d, %11.4g)\n", (int)iRow,
                 Avalue[pivot_k]);
        UpivotIndex.push_back(iRow);
        UpivotValue.push_back(Avalue[pivot_k]);
        Ustart.push_back(Uindex.size());
      }
    } else {
      assert(pivot_type == kPivotMarkowitz);
      stage = iK;
      break;
    }
    baseIndex[iRow] = iVar;
    has_pivot[iRow] = true;
  }
  if (report_lu) {
    printf("\nAfter units and singletons\n");
    reportLu(kReportLuBoth, false);
  }
  if (stage < numRow) {
    // Handle the remaining Markowitz pivots
    //
    // First of all complete the L factor with identity columns so
    // that FtranL counts the RHS entries in rows that don't yet have
    // picots by running to completion. In the hyper-sparse code,
    // these will HOPEFULLY be skipped
    //
    // There are already Lstart entries for the first stage rows, but
    // LpivotIndex is not assigned, as UpivotIndex gets copied into it
    Lstart.resize(numRow + 1);
    for (HighsInt iK = stage; iK < numRow; iK++) Lstart[iK + 1] = Lstart[iK];
    LpivotIndex.resize(numRow);
    for (HighsInt iK = 0; iK < numRow; iK++)
      LpivotIndex[iK] = this->refactor_info_.pivot_row[iK];
    // To do hyper-sparse FtranL operations, have to set up LpivotLookup.
    LpivotLookup.resize(numRow);
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      if (iRow < stage) {
        if (LpivotLookup[LpivotIndex[iRow]] != iRow) {
          //	  printf("Strange: Thought that LpivotLookup[LpivotIndex[iRow]]
          //== iRow\n");
        }
      }
      LpivotLookup[LpivotIndex[iRow]] = iRow;
    }
    // Need to know whether to consider matrix entries for FtranL
    // operation. Initially these correspond to all the rows without
    // pivots
    vector<bool> not_in_bump = has_pivot;
    // Monitor density of FtranL result to possibly switch from exploiting
    // hyper-sparsity
    double expected_density = 0.0;
    // Initialise a HVector in which the L and U entries of the
    // pivotal column will be formed
    HVector column;
    column.setup(numRow);
    for (HighsInt iK = stage; iK < numRow; iK++) {
      HighsInt iRow = this->refactor_info_.pivot_row[iK];
      HighsInt iVar = this->refactor_info_.pivot_var[iK];
      int8_t pivot_type = this->refactor_info_.pivot_type[iK];
      assert(!has_pivot[iRow]);
      assert(pivot_type == kPivotMarkowitz);
      // Set up the column for the FtranL. It contains the matrix
      // entries in rows without pivots, and the remaining entries
      // start forming the U column
      column.clear();
      HighsInt start = Astart[iVar];
      HighsInt end = Astart[iVar + 1];
      for (HighsInt iEl = start; iEl < end; iEl++) {
        HighsInt local_iRow = Aindex[iEl];
        if (not_in_bump[local_iRow]) {
          Uindex.push_back(local_iRow);
          Uvalue.push_back(Avalue[iEl]);
        } else {
          column.index[column.count++] = local_iRow;
          column.array[local_iRow] = Avalue[iEl];
        }
      }
      // Perform FtranL, but don't time it!
      ftranL(column, expected_density);
      // Update the running average density
      double local_density = (1.0 * column.count) / numRow;
      expected_density = kRunningAverageMultiplier * local_density +
                         (1 - kRunningAverageMultiplier) * expected_density;
      // Strip out small values
      column.tight();
      // Now form the column of L
      //
      // Find the pivot
      HighsInt pivot_k = -1;
      start = 0;
      end = column.count;
      for (HighsInt k = start; k < end; k++) {
        if (column.index[k] == iRow) {
          pivot_k = k;
          break;
        }
      }
      assert(pivot_k >= 0);
      // Check that the pivot isn't too small. Shouldn't happen since
      // this is refactorization
      double abs_pivot = std::fabs(column.array[iRow]);
      assert(abs_pivot >= pivot_tolerance);
      if (abs_pivot < pivot_tolerance) {
        rank_deficiency = numRow - iK;
        return rank_deficiency;
      }
      const double pivotX = 1 / column.array[iRow];
      for (HighsInt section = 0; section < 2; section++) {
        HighsInt p0 = section == 0 ? start : pivot_k + 1;
        HighsInt p1 = section == 0 ? pivot_k : end;
        for (HighsInt k = p0; k < p1; k++) {
          HighsInt local_iRow = column.index[k];
          if (!has_pivot[local_iRow]) {
            Lindex.push_back(local_iRow);
            Lvalue.push_back(column.array[local_iRow] * pivotX);
          } else {
            Uindex.push_back(local_iRow);
            Uvalue.push_back(column.array[local_iRow]);
          }
        }
      }
      Lstart[iK + 1] = Lindex.size();
      UpivotIndex.push_back(iRow);
      UpivotValue.push_back(column.array[iRow]);
      Ustart.push_back(Uindex.size());
      baseIndex[iRow] = iVar;
      has_pivot[iRow] = true;
      if (report_lu) {
        printf("\nAfter Markowitz %d\n", (int)(iK - stage));
        reportLu(kReportLuBoth, false);
      }
    }
  }
  if (report_lu) {
    printf("\nRefactored INVERT\n");
    reportLu(kReportLuBoth, false);
  }
  buildFinish();
  return 0;
}
