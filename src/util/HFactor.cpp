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
/**@file util/HFactor.cpp
 * @brief Types of solution classes
 */
#include "util/HFactor.h"

#include <cassert>
#include <iostream>

#include "lp_data/HConst.h"
#include "pdqsort/pdqsort.h"
#include "util/FactorTimer.h"
#include "util/HFactorDebug.h"
#include "util/HVector.h"
#include "util/HVectorBase.h"
#include "util/HighsTimer.h"

// std::vector, std::max and std::min used in HFactor.h for local
// in-line functions, so HFactor.h has #include <algorithm>
using std::fabs;

using std::copy;
using std::fill_n;
using std::make_pair;
using std::min;
using std::pair;

void solveMatrixT(const HighsInt Xstart, const HighsInt Xend,
                  const HighsInt Ystart, const HighsInt Yend,
                  const HighsInt* Tindex, const double* Tvalue,
                  const double Tpivot, HighsInt* RHScount, HighsInt* RHSindex,
                  double* RHSarray) {
  // Collect by X
  double pivotX = 0;
  for (HighsInt k = Xstart; k < Xend; k++)
    pivotX += Tvalue[k] * RHSarray[Tindex[k]];

  // Scatter by Y
  if (fabs(pivotX) > kHighsTiny) {
    HighsInt workCount = *RHScount;

    pivotX /= Tpivot;
    for (HighsInt k = Ystart; k < Yend; k++) {
      const HighsInt index = Tindex[k];
      const double value0 = RHSarray[index];
      const double value1 = value0 - pivotX * Tvalue[k];
      if (value0 == 0) RHSindex[workCount++] = index;
      RHSarray[index] = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
    }

    *RHScount = workCount;
  }
}

void solveHyper(const HighsInt Hsize, const HighsInt* Hlookup,
                const HighsInt* HpivotIndex, const double* HpivotValue,
                const HighsInt* Hstart, const HighsInt* Hend,
                const HighsInt* Hindex, const double* Hvalue, HVector* rhs) {
  HighsInt RHScount = rhs->count;
  HighsInt* RHSindex = &rhs->index[0];
  double* RHSarray = &rhs->array[0];

  // Take count

  // Build list
  char* listMark = &rhs->cwork[0];
  HighsInt* listIndex = &rhs->iwork[0];
  HighsInt* listStack = &rhs->iwork[Hsize];
  HighsInt listCount = 0;

  HighsInt countPivot = 0;
  HighsInt countEntry = 0;

  for (HighsInt i = 0; i < RHScount; i++) {
    // Skip touched index
    HighsInt iTrans = Hlookup[RHSindex[i]];  // XXX: this contains a bug iTran
    if (listMark[iTrans])                    // XXX bug here
      continue;

    HighsInt Hi = iTrans;      // H matrix pivot index
    HighsInt Hk = Hstart[Hi];  // H matrix non zero position
    HighsInt nStack = -1;      // Usage of the stack (-1 not used)

    listMark[Hi] = 1;  // Mark this as touched

    for (;;) {
      if (Hk < Hend[Hi]) {
        HighsInt Hi_sub = Hlookup[Hindex[Hk++]];
        if (listMark[Hi_sub] == 0) {  // Go to a child
          listMark[Hi_sub] = 1;       // Mark as touched
          listStack[++nStack] = Hi;   // Store current into stack
          listStack[++nStack] = Hk;
          Hi = Hi_sub;  // Replace current with child
          Hk = Hstart[Hi];
          if (Hi >= Hsize) {
            countPivot++;
            countEntry += Hend[Hi] - Hstart[Hi];
          }
        }
      } else {
        listIndex[listCount++] = Hi;
        if (nStack == -1)  // Quit on empty stack
          break;
        Hk = listStack[nStack--];  // Back to last in stack
        Hi = listStack[nStack--];
      }
    }
  }

  rhs->synthetic_tick += countPivot * 20 + countEntry * 10;

  // Solve with list
  if (HpivotValue == 0) {
    RHScount = 0;
    for (HighsInt iList = listCount - 1; iList >= 0; iList--) {
      HighsInt i = listIndex[iList];
      listMark[i] = 0;
      HighsInt pivotRow = HpivotIndex[i];
      double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > kHighsTiny) {
        RHSindex[RHScount++] = pivotRow;
        const HighsInt start = Hstart[i];
        const HighsInt end = Hend[i];
        for (HighsInt k = start; k < end; k++)
          RHSarray[Hindex[k]] -= pivotX * Hvalue[k];
      } else
        RHSarray[pivotRow] = 0;
    }
    rhs->count = RHScount;
  } else {
    RHScount = 0;
    for (HighsInt iList = listCount - 1; iList >= 0; iList--) {
      HighsInt i = listIndex[iList];
      listMark[i] = 0;
      HighsInt pivotRow = HpivotIndex[i];
      double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > kHighsTiny) {
        pivotX /= HpivotValue[i];
        RHSarray[pivotRow] = pivotX;
        RHSindex[RHScount++] = pivotRow;
        const HighsInt start = Hstart[i];
        const HighsInt end = Hend[i];
        for (HighsInt k = start; k < end; k++)
          RHSarray[Hindex[k]] -= pivotX * Hvalue[k];
      } else
        RHSarray[pivotRow] = 0;
    }
    rhs->count = RHScount;
  }
}

void HFactor::setupGeneral(const HighsSparseMatrix* a_matrix,
                           HighsInt num_basic_variables, HighsInt* baseIndex,
                           const double pivot_threshold,
                           const double pivot_tolerance,
                           const HighsInt highs_debug_level,
                           const HighsLogOptions* log_options) {
  this->setupGeneral(a_matrix->num_col_, a_matrix->num_row_,
                     &a_matrix->start_[0], &a_matrix->index_[0],
                     &a_matrix->value_[0], num_basic_variables, baseIndex,
                     pivot_threshold, pivot_tolerance, highs_debug_level,
                     log_options, true, kUpdateMethodFt);
}

void HFactor::setup(
    const HighsInt num_col_, const HighsInt num_row_, const HighsInt* a_start_,
    const HighsInt* a_index_, const double* a_value_, HighsInt* baseIndex_,
    const double pivot_threshold_, const double pivot_tolerance_,
    const HighsInt highs_debug_level_, const HighsLogOptions* log_options_,
    const bool use_original_HFactor_logic_, const HighsInt updateMethod_) {
  setupGeneral(num_col_, num_row_, a_start_, a_index_, a_value_, num_row_, baseIndex_,
               pivot_threshold_, pivot_tolerance_, highs_debug_level_,
               log_options_, use_original_HFactor_logic_, updateMethod_);
}

void HFactor::setupGeneral(const HighsInt num_col_, const HighsInt num_row_,
                           const HighsInt* a_start_, const HighsInt* a_index_,
                           const double* a_value_, HighsInt num_basic_variables,
                           HighsInt* baseIndex_, const double pivot_threshold_,
                           const double pivot_tolerance_,
                           const HighsInt highs_debug_level_,
                           const HighsLogOptions* log_options_,
                           const bool use_original_HFactor_logic_,
                           const HighsInt updateMethod_) {
  // Copy Problem size and (pointer to) coefficient matrix
  num_row = num_row_;
  num_col = num_col_;
  this->a_matrix_valid = true;
  a_start = a_start_;
  a_index = a_index_;
  a_value = a_value_;
  baseIndex = baseIndex_;
  pivot_threshold =
      max(kMinPivotThreshold, min(pivot_threshold_, kMaxPivotThreshold));
  pivot_tolerance =
      max(kMinPivotTolerance, min(pivot_tolerance_, kMaxPivotTolerance));
  highs_debug_level = highs_debug_level_;
  if (!log_options_) {
    output_flag = false;
    log_file_stream = NULL;
    log_to_console = true;
    log_dev_level = 0;
  } else {
    output_flag = *(log_options_->output_flag);
    log_file_stream = log_options_->log_file_stream;
    log_to_console = *(log_options_->log_to_console);
    log_dev_level = *(log_options_->log_dev_level);
    log_options = *log_options_;
  }
  use_original_HFactor_logic = use_original_HFactor_logic_;
  updateMethod = updateMethod_;

  // Allocate for working buffer
  iwork.reserve(num_row * 2);
  dwork.assign(num_row, 0);

  // Find Basis matrix limit size
  BlimitX = 0;
  iwork.assign(num_row + 1, 0);
  for (HighsInt i = 0; i < num_col; i++) iwork[a_start[i + 1] - a_start[i]]++;
  for (HighsInt i = num_row, counted = 0; i >= 0 && counted < num_row; i--)
    BlimitX += i * iwork[i], counted += iwork[i];
  BlimitX += num_row;

  // Allocate space for basis matrix, L, U factor and Update buffer
  Bvar.resize(num_row);
  Bstart.resize(num_row + 1, 0);
  Bindex.resize(BlimitX);
  Bvalue.resize(BlimitX);

  // Allocate space for pivot records
  permute.resize(num_row);

  // Allocate space for Markowitz matrices
  MCvar.resize(num_row);
  MCstart.resize(num_row);
  MCcountA.resize(num_row);
  MCcountN.resize(num_row);
  MCspace.resize(num_row);
  MCminpivot.resize(num_row);
  MCindex.resize(BlimitX * kMCExtraEntriesMultiplier);
  MCvalue.resize(BlimitX * kMCExtraEntriesMultiplier);

  MRstart.resize(num_row);
  MRcount.resize(num_row);
  MRspace.resize(num_row);
  MRcountb4.resize(num_row);
  MRindex.resize(BlimitX * kMRExtraEntriesMultiplier);

  mwz_column_mark.assign(num_row, 0);
  mwz_column_index.resize(num_row);
  mwz_column_array.assign(num_row, 0);

  // Allocate space for count-link-list
  col_link_first.assign(num_row + 1, -1);
  col_link_next.resize(num_row);
  col_link_last.resize(num_row);

  row_link_first.assign(num_row + 1, -1);
  row_link_next.resize(num_row);
  row_link_last.resize(num_row);

  // Allocate space for L factor
  l_pivot_lookup.resize(num_row);
  l_pivot_index.reserve(num_row);
  l_start.reserve(num_row + 1);
  l_index.reserve(BlimitX * kLFactorExtraEntriesMultiplier);
  l_value.reserve(BlimitX * kLFactorExtraEntriesMultiplier);

  lr_start.reserve(num_row + 1);
  lr_index.reserve(BlimitX * kLFactorExtraEntriesMultiplier);
  lr_value.reserve(BlimitX * kLFactorExtraEntriesMultiplier);

  // Allocate space for U factor
  u_pivot_lookup.resize(num_row);
  u_pivot_index.reserve(num_row + kUFactorExtraVectors);
  u_pivot_value.reserve(num_row + kUFactorExtraVectors);

  u_start.reserve(num_row + kUFactorExtraVectors + 1);
  u_last_p.reserve(num_row + kUFactorExtraVectors);
  u_index.reserve(BlimitX * kUFactorExtraEntriesMultiplier);
  u_value.reserve(BlimitX * kUFactorExtraEntriesMultiplier);

  ur_start.reserve(num_row + kUFactorExtraVectors + 1);
  ur_lastp.reserve(num_row + kUFactorExtraVectors);
  ur_space.reserve(num_row + kUFactorExtraVectors);
  ur_index.reserve(BlimitX * kUFactorExtraEntriesMultiplier);
  ur_value.reserve(BlimitX * kUFactorExtraEntriesMultiplier);

  // Allocate spaces for Update buffer
  PFpivotValue.reserve(kPFFPivotEntries);
  PFpivotIndex.reserve(kPFFPivotEntries);
  PFstart.reserve(kPFVectors + 1);
  PFindex.reserve(BlimitX * kPFEntriesMultiplier);
  PFvalue.reserve(BlimitX * kPFEntriesMultiplier);
}

void HFactor::setupMatrix(const HighsInt* a_start_, const HighsInt* a_index_,
                          const double* a_value_) {
  a_start = a_start_;
  a_index = a_index_;
  a_value = a_value_;
  this->a_matrix_valid = true;
}

void HFactor::setupMatrix(const HighsSparseMatrix* a_matrix) {
  setupMatrix(&a_matrix->start_[0], &a_matrix->index_[0], &a_matrix->value_[0]);
}

HighsInt HFactor::build(HighsTimerClock* factor_timer_clock_pointer) {
  const bool report_lu = false;
  // Ensure that the A matrix is valid for factorization
  assert(this->a_matrix_valid);
  FactorTimer factor_timer;
  // Possibly use the refactorization information!
  if (refactor_info_.use) {
    factor_timer.start(FactorReinvert, factor_timer_clock_pointer);
    rank_deficiency = rebuild(factor_timer_clock_pointer);
    factor_timer.stop(FactorReinvert, factor_timer_clock_pointer);
    if (!rank_deficiency) return 0;
  }
  // Refactoring from just the list of basic variables. Initialise the
  // refactorization information.
  refactor_info_.clear();
  // Start the timer
  factor_timer.start(FactorInvert, factor_timer_clock_pointer);
  build_synthetic_tick = 0;
  factor_timer.start(FactorInvertSimple, factor_timer_clock_pointer);
  // Build the L, U factor
  buildSimple();
  if (report_lu) {
    printf("\nAfter units and singletons\n");
    reportLu(kReportLuBoth, false);
  }
  factor_timer.stop(FactorInvertSimple, factor_timer_clock_pointer);
  factor_timer.start(FactorInvertKernel, factor_timer_clock_pointer);
  rank_deficiency = buildKernel();
  factor_timer.stop(FactorInvertKernel, factor_timer_clock_pointer);
  if (rank_deficiency) {
    factor_timer.start(FactorInvertDeficient, factor_timer_clock_pointer);
    highsLogDev(log_options, HighsLogType::kWarning,
                "Rank deficiency of %" HIGHSINT_FORMAT
                " identified in basis matrix\n",
                rank_deficiency);
    // Singular matrix B: reorder the basic variables so that the
    // singular columns are in the position corresponding to the
    // logical which replaces them
    buildHandleRankDeficiency();
    buildMarkSingC();
    factor_timer.stop(FactorInvertDeficient, factor_timer_clock_pointer);
  }
  // Complete INVERT
  if (report_lu) {
    printf("\nFactored INVERT\n");
    reportLu(kReportLuBoth, false);
  }
  factor_timer.start(FactorInvertFinish, factor_timer_clock_pointer);
  buildFinish();
  factor_timer.stop(FactorInvertFinish, factor_timer_clock_pointer);
  //
  // Indicate that the refactorization information is known unless the basis was
  // rank deficient
  if (rank_deficiency) {
    this->refactor_info_.clear();
  } else {
    // Check that the refactorization information is not (yet) flagged
    // to be used in a future call
    assert(!this->refactor_info_.use);
    // Record build_synthetic_tick to use as the value of
    // build_synthetic_tick if refactorization is performed. Not only
    // is there no build_synthetic_tick measure for refactorization,
    // if there were it would give an unrealistic underestimate of the
    // cost of factorization from scratch
    this->refactor_info_.build_synthetic_tick = this->build_synthetic_tick;
  }

  // Record the number of entries in the INVERT
  invert_num_el = l_start[num_row] + u_last_p[num_row - 1] + num_row;

  kernel_dim -= rank_deficiency;
  debugLogRankDeficiency(highs_debug_level, log_options, rank_deficiency,
                         basis_matrix_num_el, invert_num_el, kernel_dim,
                         kernel_num_el, nwork);
  factor_timer.stop(FactorInvert, factor_timer_clock_pointer);
  return rank_deficiency;
}

void HFactor::ftranCall(HVector& vector, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorFtran, factor_timer_clock_pointer);
  ftranL(vector, expected_density, factor_timer_clock_pointer);
  ftranU(vector, expected_density, factor_timer_clock_pointer);
  factor_timer.stop(FactorFtran, factor_timer_clock_pointer);
}

void HFactor::btranCall(HVector& vector, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorBtran, factor_timer_clock_pointer);
  btranU(vector, expected_density, factor_timer_clock_pointer);
  btranL(vector, expected_density, factor_timer_clock_pointer);
  factor_timer.stop(FactorBtran, factor_timer_clock_pointer);
}

void HFactor::update(HVector* aq, HVector* ep, HighsInt* iRow, HighsInt* hint) {
  // Updating implies a change of basis. Since the refactorizaion info
  // no longer corresponds to the current basis, it must be
  // invalidated
  this->refactor_info_.clear();
  // Special case
  if (aq->next) {
    updateCFT(aq, ep, iRow);
    return;
  }

  if (updateMethod == kUpdateMethodFt) updateFT(aq, ep, *iRow);
  if (updateMethod == kUpdateMethodPf) updatePF(aq, *iRow, hint);
  if (updateMethod == kUpdateMethodMpf) updateMPF(aq, ep, *iRow, hint);
  if (updateMethod == kUpdateMethodApf) updateAPF(aq, ep, *iRow);
}

bool HFactor::setPivotThreshold(const double new_pivot_threshold) {
  if (new_pivot_threshold < kMinPivotThreshold) return false;
  if (new_pivot_threshold > kMaxPivotThreshold) return false;
  pivot_threshold = new_pivot_threshold;
  return true;
}

void HFactor::luClear() {
  l_start.clear();
  l_start.push_back(0);
  l_index.clear();
  l_value.clear();

  u_pivot_index.clear();
  u_pivot_value.clear();
  u_start.clear();
  u_start.push_back(0);
  u_index.clear();
  u_value.clear();
}

void HFactor::buildSimple() {
  /**
   * 0. Clear L and U factor
   */
  luClear();

  // Set all values of permute to -1 so that unpermuted (rank
  // deficient) columns canm be identified
  permute.assign(num_row, -1);

  /**
   * 1. Prepare basis matrix and deal with unit columns
   */

  const bool report_unit = false;
  const bool report_singletons = false;
  const bool report_markowitz = false;
  const bool report_anything =
      report_unit || report_singletons || report_markowitz;
  HighsInt BcountX = 0;
  fill_n(&MRcountb4[0], num_row, 0);
  nwork = 0;
  if (report_anything) printf("\nFactor\n");
  // Compile a vector iwork of the nwork non-unit structural columns
  // in baseindex: they will be formed into the B matrix as the kernel
  for (HighsInt iCol = 0; iCol < num_row; iCol++) {
    HighsInt iMat = baseIndex[iCol];
    HighsInt iRow = -1;
    int8_t pivot_type = kPivotIllegal;
    if (iMat >= num_col) {
      if (report_unit) printf("Stage %d: Logical\n", (int)(l_start.size() - 1));
      // 1.1 Logical column
      pivot_type = kPivotLogical;
      // Check for double pivot
      HighsInt lc_iRow = iMat - num_col;
      if (MRcountb4[lc_iRow] >= 0) {
        iRow = lc_iRow;
      } else {
        highsLogDev(log_options, HighsLogType::kError,
                    "INVERT Error: Found a logical column with pivot "
                    "already in row %" HIGHSINT_FORMAT "\n",
                    lc_iRow);
        // Treat this as a column to be handled in the kernel, so that
        // the rank deficiency is detected as late as possible.
        MRcountb4[lc_iRow]++;
        Bindex[BcountX] = lc_iRow;
        Bvalue[BcountX++] = 1.0;
        iwork[nwork++] = iCol;
      }
    } else {
      // 1.2 Structural column
      HighsInt start = a_start[iMat];
      HighsInt count = a_start[iMat + 1] - start;
      HighsInt lc_iRow = a_index[start];
      // Check for unit column with double pivot
      bool unit_col = count == 1 && a_value[start] == 1;
      if (unit_col && MRcountb4[lc_iRow] >= 0) {
        if (report_unit) printf("Stage %d: Unit\n", (int)(l_start.size() - 1));
        // Don't exploit this special case in case the matrix is
        // re-factorized after scaling has been applied, making this
        // column non-unit.
        pivot_type = kPivotColSingleton;  //;kPivotUnit;//
        iRow = lc_iRow;
      } else {
        if (unit_col)
          highsLogDev(log_options, HighsLogType::kError,
                      "INVERT Error: Found a second unit column with pivot in "
                      "row %" HIGHSINT_FORMAT "\n",
                      lc_iRow);
        for (HighsInt k = start; k < start + count; k++) {
          MRcountb4[a_index[k]]++;
          if (BcountX >= BlimitX) {
            printf("BcountX >= BlimitX\n");
            fflush(stdout);
          }

          assert(BcountX < BlimitX);
          Bindex[BcountX] = a_index[k];
          Bvalue[BcountX++] = a_value[k];
        }
        iwork[nwork++] = iCol;
      }
    }

    if (iRow >= 0) {
      // 1.3 Record unit column
      permute[iCol] = iRow;
      l_start.push_back(l_index.size());
      u_pivot_index.push_back(iRow);
      u_pivot_value.push_back(1);
      u_start.push_back(u_index.size());
      MRcountb4[iRow] = -num_row;
      assert(pivot_type != kPivotIllegal);
      this->refactor_info_.pivot_row.push_back(iRow);
      this->refactor_info_.pivot_var.push_back(iMat);
      this->refactor_info_.pivot_type.push_back(pivot_type);
    }
    Bstart[iCol + 1] = BcountX;
    Bvar[iCol] = iMat;
  }
  // Record the number of elements in the basis matrix
  basis_matrix_num_el = num_row - nwork + BcountX;

  // count1 = 0;
  // Comments: for pds-20, dfl001: 60 / 80
  // Comments: when system is large: enlarge
  // Comments: when system is small: decrease
  build_synthetic_tick += BcountX * 60 + (num_row - nwork) * 80;

  /**
   * 2. Search for and deal with singletons
   */
  double t2_search = 0;
  double t2_storeL = l_index.size();
  double t2_storeU = u_index.size();
  double t2_storep = nwork;
  while (nwork > 0) {
    HighsInt nworkLast = nwork;
    nwork = 0;
    for (HighsInt i = 0; i < nworkLast; i++) {
      const HighsInt iCol = iwork[i];
      const HighsInt start = Bstart[iCol];
      const HighsInt end = Bstart[iCol + 1];
      HighsInt pivot_k = -1;
      HighsInt found_row_singleton = 0;
      HighsInt count = 0;

      // 2.1 Search for singleton
      t2_search += end - start;
      for (HighsInt k = start; k < end; k++) {
        const HighsInt iRow = Bindex[k];
        if (MRcountb4[iRow] == 1) {
          pivot_k = k;
          found_row_singleton = 1;
          break;
        }
        if (MRcountb4[iRow] > 1) {
          pivot_k = k;
          count++;
        }
      }

      if (found_row_singleton) {
        // 2.2 Deal with row singleton
        const double pivotX = 1 / Bvalue[pivot_k];
        if (report_singletons)
          printf("Stage %d: Row singleton (%4d, %g)\n",
                 (int)(l_start.size() - 1), (int)pivot_k, pivotX);
        for (HighsInt section = 0; section < 2; section++) {
          HighsInt p0 = section == 0 ? start : pivot_k + 1;
          HighsInt p1 = section == 0 ? pivot_k : end;
          for (HighsInt k = p0; k < p1; k++) {
            HighsInt iRow = Bindex[k];
            if (MRcountb4[iRow] > 0) {
              if (report_singletons)
                printf("Row singleton: L En (%4d, %11.4g)\n", (int)iRow,
                       Bvalue[k] * pivotX);
              l_index.push_back(iRow);
              l_value.push_back(Bvalue[k] * pivotX);
            } else {
              if (report_singletons)
                printf("Row singleton: U En (%4d, %11.4g)\n", (int)iRow,
                       Bvalue[k]);
              u_index.push_back(iRow);
              u_value.push_back(Bvalue[k]);
            }
            MRcountb4[iRow]--;
          }
        }
        HighsInt iRow = Bindex[pivot_k];
        MRcountb4[iRow] = 0;
        permute[iCol] = iRow;
        l_start.push_back(l_index.size());

        if (report_singletons)
          printf("Row singleton: U Pv (%4d, %11.4g)\n", (int)iRow,
                 Bvalue[pivot_k]);
        u_pivot_index.push_back(iRow);
        u_pivot_value.push_back(Bvalue[pivot_k]);
        u_start.push_back(u_index.size());
        assert(Bvar[iCol] == baseIndex[iCol]);

        this->refactor_info_.pivot_row.push_back(iRow);
        this->refactor_info_.pivot_var.push_back(baseIndex[iCol]);
        this->refactor_info_.pivot_type.push_back(kPivotRowSingleton);
      } else if (count == 1) {
        if (report_singletons)
          printf("Stage %d: Col singleton \n", (int)(l_start.size() - 1));
        // 2.3 Deal with column singleton
        for (HighsInt k = start; k < pivot_k; k++) {
          if (report_singletons)
            printf("Col singleton: U En (%4d, %11.4g)\n", (int)Bindex[k],
                   Bvalue[k]);
          u_index.push_back(Bindex[k]);
          u_value.push_back(Bvalue[k]);
        }
        for (HighsInt k = pivot_k + 1; k < end; k++) {
          if (report_singletons)
            printf("Col singleton: U En (%4d, %11.4g)\n", (int)Bindex[k],
                   Bvalue[k]);
          u_index.push_back(Bindex[k]);
          u_value.push_back(Bvalue[k]);
        }

        HighsInt iRow = Bindex[pivot_k];
        MRcountb4[iRow] = 0;
        permute[iCol] = iRow;
        l_start.push_back(l_index.size());

        if (report_singletons)
          printf("Col singleton: U Pv (%4d, %11.4g)\n", (int)iRow,
                 Bvalue[pivot_k]);
        u_pivot_index.push_back(iRow);
        u_pivot_value.push_back(Bvalue[pivot_k]);
        u_start.push_back(u_index.size());
        assert(Bvar[iCol] == baseIndex[iCol]);
        this->refactor_info_.pivot_row.push_back(iRow);
        this->refactor_info_.pivot_var.push_back(baseIndex[iCol]);
        this->refactor_info_.pivot_type.push_back(kPivotColSingleton);
      } else {
        iwork[nwork++] = iCol;
      }
    }

    // No singleton found in the last pass
    if (nworkLast == nwork) break;
  }
  if (report_anything) reportLu(kReportLuBoth, false);
  t2_storeL = l_index.size() - t2_storeL;
  t2_storeU = u_index.size() - t2_storeU;
  t2_storep = t2_storep - nwork;

  build_synthetic_tick +=
      t2_search * 20 + (t2_storep + t2_storeL + t2_storeU) * 80;

  /**
   * 3. Prepare the kernel parts
   */
  // 3.1 Prepare row links, row matrix spaces
  row_link_first.assign(num_row + 1, -1);
  MRcount.assign(num_row, 0);
  HighsInt MRcountX = 0;
  // Determine the number of entries in the kernel
  kernel_num_el = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt count = MRcountb4[iRow];
    if (count > 0) {
      MRstart[iRow] = MRcountX;
      MRspace[iRow] = count * 2;
      MRcountX += count * 2;
      rlinkAdd(iRow, count);
      kernel_num_el += count + 1;
    }
  }
  MRindex.resize(MRcountX);

  // 3.2 Prepare column links, kernel matrix
  col_link_first.assign(num_row + 1, -1);
  MCindex.clear();
  MCvalue.clear();
  MCcountA.assign(num_row, 0);
  MCcountN.assign(num_row, 0);
  HighsInt MCcountX = 0;
  for (HighsInt i = 0; i < nwork; i++) {
    HighsInt iCol = iwork[i];
    MCvar[iCol] = Bvar[iCol];
    MCstart[iCol] = MCcountX;
    MCspace[iCol] = (Bstart[iCol + 1] - Bstart[iCol]) * 2;
    MCcountX += MCspace[iCol];
    MCindex.resize(MCcountX);
    MCvalue.resize(MCcountX);
    for (HighsInt k = Bstart[iCol]; k < Bstart[iCol + 1]; k++) {
      const HighsInt iRow = Bindex[k];
      const double value = Bvalue[k];
      if (MRcountb4[iRow] > 0) {
        colInsert(iCol, iRow, value);
        rowInsert(iCol, iRow);
      } else {
        colStoreN(iCol, iRow, value);
      }
    }
    colFixMax(iCol);
    clinkAdd(iCol, MCcountA[iCol]);
  }
  build_synthetic_tick += (num_row + nwork + MCcountX) * 40 + MRcountX * 20;
  // Record the kernel dimension
  kernel_dim = nwork;
  assert((HighsInt)this->refactor_info_.pivot_row.size() == num_row - nwork);
}

HighsInt HFactor::buildKernel() {
  // Deal with the kernel part by 'n-work' pivoting

  double fake_search = 0;
  double fake_fill = 0;
  double fake_eliminate = 0;

  while (nwork-- > 0) {
    /**
     * 1. Search for the pivot
     */
    HighsInt jColPivot = -1;
    HighsInt iRowPivot = -1;
    //    int8_t pivot_type = kPivotIllegal;
    // 1.1. Setup search merits
    HighsInt searchLimit = min(nwork, HighsInt{8});
    HighsInt searchCount = 0;
    double meritLimit = 1.0 * num_row * num_row;
    double meritPivot = meritLimit;

    // 1.2. Search for local singletons
    bool foundPivot = false;
    if (!foundPivot && col_link_first[1] != -1) {
      jColPivot = col_link_first[1];
      iRowPivot = MCindex[MCstart[jColPivot]];
      foundPivot = true;
    }
    if (!foundPivot && row_link_first[1] != -1) {
      iRowPivot = row_link_first[1];
      jColPivot = MRindex[MRstart[iRowPivot]];
      foundPivot = true;
    }
    const bool singleton_pivot = foundPivot;
    // 1.3. Major search loop
    double candidate_pivot_value = 0;
    for (HighsInt count = 2; !foundPivot && count <= num_row; count++) {
      // 1.3.1 Search for columns
      for (HighsInt j = col_link_first[count]; j != -1; j = col_link_next[j]) {
        double minpivot = MCminpivot[j];
        HighsInt start = MCstart[j];
        HighsInt end = start + MCcountA[j];
        for (HighsInt k = start; k < end; k++) {
          if (fabs(MCvalue[k]) >= minpivot) {
            HighsInt i = MCindex[k];
            HighsInt rowCount = MRcount[i];
            double meritLocal = 1.0 * (count - 1) * (rowCount - 1);
            if (meritPivot > meritLocal) {
              candidate_pivot_value = fabs(MCvalue[k]);
              meritPivot = meritLocal;
              jColPivot = j;
              iRowPivot = i;
              foundPivot = foundPivot || (rowCount < count);
            }
          }
        }

        if (searchCount++ >= searchLimit && meritPivot < meritLimit)
          foundPivot = true;
        if (foundPivot) break;

        fake_search += count;
      }

      // 1.3.2 Search for rows
      for (HighsInt i = row_link_first[count]; i != -1; i = row_link_next[i]) {
        HighsInt start = MRstart[i];
        HighsInt end = start + MRcount[i];
        for (HighsInt k = start; k < end; k++) {
          HighsInt j = MRindex[k];
          HighsInt columnCount = MCcountA[j];
          double meritLocal = 1.0 * (count - 1) * (columnCount - 1);
          if (meritLocal < meritPivot) {
            HighsInt ifind = MCstart[j];
            while (MCindex[ifind] != i) ifind++;
            if (fabs(MCvalue[ifind]) >= MCminpivot[j]) {
              candidate_pivot_value = fabs(MCvalue[ifind]);
              meritPivot = meritLocal;
              jColPivot = j;
              iRowPivot = i;
              foundPivot = foundPivot || (columnCount <= count);
            }
          }
        }
        if (searchCount++ >= searchLimit && meritPivot < meritLimit)
          foundPivot = true;
        if (foundPivot) break;
      }

      fake_search += count;
    }

    // 1.4. If we found nothing: tell singular
    if (!foundPivot) {
      rank_deficiency = nwork + 1;
      return rank_deficiency;
    }

    /**
     * 2. Elimination other elements by the pivot
     */
    // 2.1. Delete the pivot
    double pivotX = colDelete(jColPivot, iRowPivot);
    if (!singleton_pivot) assert(candidate_pivot_value == fabs(pivotX));
    if (fabs(pivotX) < pivot_tolerance) {
      highsLogDev(log_options, HighsLogType::kWarning,
                  "Small |pivot| = %g when nwork = %" HIGHSINT_FORMAT "\n",
                  fabs(pivotX), nwork);
      rank_deficiency = nwork + 1;
      assert((HighsInt)this->refactor_info_.pivot_row.size() +
                 rank_deficiency ==
             num_row);
      return rank_deficiency;
    }
    rowDelete(jColPivot, iRowPivot);
    clinkDel(jColPivot);
    rlinkDel(iRowPivot);
    permute[jColPivot] = iRowPivot;
    //    printf("Mwz pivot %3d; MCvar[%3d] = %d; baseIndex[%3d] = %d\n",
    //    (int)(num_row-nwork),
    //	   (int)jColPivot, (int)MCvar[jColPivot],
    //	   (int)jColPivot, (int)baseIndex[jColPivot]);
    assert(MCvar[jColPivot] == baseIndex[jColPivot]);

    this->refactor_info_.pivot_row.push_back(iRowPivot);
    this->refactor_info_.pivot_var.push_back(baseIndex[jColPivot]);
    this->refactor_info_.pivot_type.push_back(kPivotMarkowitz);

    // 2.2. Store active pivot column to L
    HighsInt start_A = MCstart[jColPivot];
    HighsInt end_A = start_A + MCcountA[jColPivot];
    HighsInt mwz_column_count = 0;
    for (HighsInt k = start_A; k < end_A; k++) {
      const HighsInt iRow = MCindex[k];
      const double value = MCvalue[k] / pivotX;
      mwz_column_index[mwz_column_count++] = iRow;
      mwz_column_array[iRow] = value;
      mwz_column_mark[iRow] = 1;
      l_index.push_back(iRow);
      l_value.push_back(value);
      MRcountb4[iRow] = MRcount[iRow];
      rowDelete(jColPivot, (int)iRow);
    }
    l_start.push_back(l_index.size());
    fake_fill += 2 * MCcountA[jColPivot];

    // 2.3. Store non active pivot column to U
    HighsInt end_N = start_A + MCspace[jColPivot];
    HighsInt start_N = end_N - MCcountN[jColPivot];
    for (HighsInt i = start_N; i < end_N; i++) {
      u_index.push_back(MCindex[i]);
      u_value.push_back(MCvalue[i]);
    }
    u_pivot_index.push_back(iRowPivot);
    u_pivot_value.push_back(pivotX);
    u_start.push_back(u_index.size());
    fake_fill += end_N - start_N;

    // 2.4. Loop over pivot row to eliminate other column
    const HighsInt row_start = MRstart[iRowPivot];
    const HighsInt row_end = row_start + MRcount[iRowPivot];
    for (HighsInt row_k = row_start; row_k < row_end; row_k++) {
      // 2.4.1. My pointer
      HighsInt iCol = MRindex[row_k];
      const HighsInt my_count = MCcountA[iCol];
      const HighsInt my_start = MCstart[iCol];
      const HighsInt my_end = my_start + my_count - 1;
      double my_pivot = colDelete(iCol, iRowPivot);
      colStoreN(iCol, iRowPivot, my_pivot);

      // 2.4.2. Elimination on the overlapping part
      HighsInt nFillin = mwz_column_count;
      HighsInt nCancel = 0;
      for (HighsInt my_k = my_start; my_k < my_end; my_k++) {
        HighsInt iRow = MCindex[my_k];
        double value = MCvalue[my_k];
        if (mwz_column_mark[iRow]) {
          mwz_column_mark[iRow] = 0;
          nFillin--;
          value -= my_pivot * mwz_column_array[iRow];
          if (fabs(value) < kHighsTiny) {
            value = 0;
            nCancel++;
          }
          MCvalue[my_k] = value;
        }
      }
      fake_eliminate += mwz_column_count;
      fake_eliminate += nFillin * 2;

      // 2.4.3. Remove cancellation gaps
      if (nCancel > 0) {
        HighsInt new_end = my_start;
        for (HighsInt my_k = my_start; my_k < my_end; my_k++) {
          if (MCvalue[my_k] != 0) {
            MCindex[new_end] = MCindex[my_k];
            MCvalue[new_end++] = MCvalue[my_k];
          } else {
            rowDelete(iCol, MCindex[my_k]);
          }
        }
        MCcountA[iCol] = new_end - my_start;
      }

      // 2.4.4. Insert fill-in
      if (nFillin > 0) {
        // 2.4.4.1 Check column size
        if (MCcountA[iCol] + MCcountN[iCol] + nFillin > MCspace[iCol]) {
          // p1&2=active, p3&4=non active, p5=new p1, p7=new p3
          HighsInt p1 = MCstart[iCol];
          HighsInt p2 = p1 + MCcountA[iCol];
          HighsInt p3 = p1 + MCspace[iCol] - MCcountN[iCol];
          HighsInt p4 = p1 + MCspace[iCol];
          MCspace[iCol] += max(MCspace[iCol], nFillin);
          HighsInt p5 = MCstart[iCol] = MCindex.size();
          HighsInt p7 = p5 + MCspace[iCol] - MCcountN[iCol];
          MCindex.resize(p5 + MCspace[iCol]);
          MCvalue.resize(p5 + MCspace[iCol]);
          copy(&MCindex[p1], &MCindex[p2], &MCindex[p5]);
          copy(&MCvalue[p1], &MCvalue[p2], &MCvalue[p5]);
          copy(&MCindex[p3], &MCindex[p4], &MCindex[p7]);
          copy(&MCvalue[p3], &MCvalue[p4], &MCvalue[p7]);
        }

        // 2.4.4.2 Fill into column copy
        for (HighsInt i = 0; i < mwz_column_count; i++) {
          HighsInt iRow = mwz_column_index[i];
          if (mwz_column_mark[iRow])
            colInsert(iCol, iRow, -my_pivot * mwz_column_array[iRow]);
        }

        // 2.4.4.3 Fill into the row copy
        for (HighsInt i = 0; i < mwz_column_count; i++) {
          HighsInt iRow = mwz_column_index[i];
          if (mwz_column_mark[iRow]) {
            // Expand row space
            if (MRcount[iRow] == MRspace[iRow]) {
              HighsInt p1 = MRstart[iRow];
              HighsInt p2 = p1 + MRcount[iRow];
              HighsInt p3 = MRstart[iRow] = MRindex.size();
              MRspace[iRow] *= 2;
              MRindex.resize(p3 + MRspace[iRow]);
              copy(&MRindex[p1], &MRindex[p2], &MRindex[p3]);
            }
            rowInsert(iCol, iRow);
          }
        }
      }

      // 2.4.5. Reset pivot column mark
      for (HighsInt i = 0; i < mwz_column_count; i++)
        mwz_column_mark[mwz_column_index[i]] = 1;

      // 2.4.6. Fix max value and link list
      colFixMax(iCol);
      if (my_count != MCcountA[iCol]) {
        clinkDel(iCol);
        clinkAdd(iCol, MCcountA[iCol]);
      }
    }

    // 2.5. Clear pivot column buffer
    for (HighsInt i = 0; i < mwz_column_count; i++)
      mwz_column_mark[mwz_column_index[i]] = 0;

    // 2.6. Correct row links for the remain active part
    for (HighsInt i = start_A; i < end_A; i++) {
      HighsInt iRow = MCindex[i];
      if (MRcountb4[iRow] != MRcount[iRow]) {
        rlinkDel(iRow);
        rlinkAdd(iRow, MRcount[iRow]);
      }
    }
  }
  build_synthetic_tick +=
      fake_search * 20 + fake_fill * 160 + fake_eliminate * 80;
  rank_deficiency = 0;
  return rank_deficiency;
}

void HFactor::buildHandleRankDeficiency() {
  debugReportRankDeficiency(0, highs_debug_level, log_options, num_row, permute,
                            iwork, baseIndex, rank_deficiency,
                            row_with_no_pivot, col_with_no_pivot);
  // iwork can now be used as workspace: use it to accumulate the new
  // baseIndex. iwork is set to -1 and baseIndex is permuted into it.
  // Indices of iwork corresponding to missing indices in permute
  // remain -1. Hence the -1's become markers for the logicals which
  // will replace singular columns. Once baseIndex[i] is read, it can
  // be used to pack up the entries in baseIndex which are not
  // permuted anywhere - and so will be singular columns.
  row_with_no_pivot.resize(rank_deficiency);
  col_with_no_pivot.resize(rank_deficiency);
  HighsInt lc_rank_deficiency = 0;
  for (HighsInt i = 0; i < num_row; i++) iwork[i] = -1;
  for (HighsInt i = 0; i < num_row; i++) {
    HighsInt perm_i = permute[i];
    if (perm_i >= 0) {
      iwork[perm_i] = baseIndex[i];
    } else {
      col_with_no_pivot[lc_rank_deficiency++] = i;
    }
  }
  assert(lc_rank_deficiency == rank_deficiency);
  lc_rank_deficiency = 0;
  for (HighsInt i = 0; i < num_row; i++) {
    if (iwork[i] < 0) {
      // Record the rows with no pivots in row_with_no_pivot and indicate them
      // within iwork by storing the negation of one more than their
      // rank deficiency counter [since we can't have -0].
      row_with_no_pivot[lc_rank_deficiency] = i;
      iwork[i] = -(lc_rank_deficiency + 1);
      lc_rank_deficiency++;
    }
  }
  assert(lc_rank_deficiency == rank_deficiency);
  debugReportRankDeficiency(1, highs_debug_level, log_options, num_row, permute,
                            iwork, baseIndex, rank_deficiency,
                            row_with_no_pivot, col_with_no_pivot);
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    HighsInt iRow = row_with_no_pivot[k];
    HighsInt iCol = col_with_no_pivot[k];
    assert(permute[iCol] == -1);
    permute[iCol] = iRow;
    l_start.push_back(l_index.size());
    u_pivot_index.push_back(iRow);
    u_pivot_value.push_back(1);
    u_start.push_back(u_index.size());
  }
  debugReportRankDeficiency(2, highs_debug_level, log_options, num_row, permute,
                            iwork, baseIndex, rank_deficiency,
                            row_with_no_pivot, col_with_no_pivot);
  debugReportRankDeficientASM(
      highs_debug_level, log_options, num_row, MCstart, MCcountA, MCindex,
      MCvalue, iwork, rank_deficiency, col_with_no_pivot, row_with_no_pivot);
}

void HFactor::buildMarkSingC() {
  // Singular matrix B: reorder the basic variables so that the
  // singular columns are in the position corresponding to the
  // logical which replaces them
  debugReportMarkSingC(0, highs_debug_level, log_options, num_row, iwork,
                       baseIndex);

  var_with_no_pivot.resize(rank_deficiency);
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    HighsInt ASMrow = row_with_no_pivot[k];
    HighsInt ASMcol = col_with_no_pivot[k];
    assert(-iwork[ASMrow] - 1 >= 0 && -iwork[ASMrow] - 1 < rank_deficiency);
    // Store negation of 1+ASMcol so that removing column 0 can be
    // identified!
    iwork[ASMrow] = -(ASMcol + 1);
    var_with_no_pivot[k] = baseIndex[ASMcol];
    baseIndex[ASMcol] = num_col + ASMrow;
  }
  debugReportMarkSingC(1, highs_debug_level, log_options, num_row, iwork,
                       baseIndex);
}

void HFactor::buildFinish() {
  //  debugPivotValueAnalysis(highs_debug_level, log_options, num_row,
  //  u_pivot_value);
  // The look up table
  for (HighsInt i = 0; i < num_row; i++) u_pivot_lookup[u_pivot_index[i]] = i;
  l_pivot_index = u_pivot_index;
  l_pivot_lookup = u_pivot_lookup;

  // LR space
  HighsInt LcountX = l_index.size();
  lr_index.resize(LcountX);
  lr_value.resize(LcountX);

  // LR pointer
  iwork.assign(num_row, 0);
  for (HighsInt k = 0; k < LcountX; k++) iwork[l_pivot_lookup[l_index[k]]]++;

  lr_start.assign(num_row + 1, 0);
  for (HighsInt i = 1; i <= num_row; i++)
    lr_start[i] = lr_start[i - 1] + iwork[i - 1];

  // LR elements
  iwork.assign(&lr_start[0], &lr_start[num_row]);
  for (HighsInt i = 0; i < num_row; i++) {
    const HighsInt index = l_pivot_index[i];
    for (HighsInt k = l_start[i]; k < l_start[i + 1]; k++) {
      HighsInt iRow = l_pivot_lookup[l_index[k]];
      HighsInt iPut = iwork[iRow]++;
      lr_index[iPut] = index;
      lr_value[iPut] = l_value[k];
    }
  }

  // U pointer
  u_start.push_back(0);
  u_last_p.assign(&u_start[1], &u_start[num_row + 1]);
  u_start.resize(num_row);

  // UR space
  HighsInt UcountX = u_index.size();
  HighsInt URstuffX = updateMethod == kUpdateMethodFt ? 5 : 0;
  HighsInt URcountX = UcountX + URstuffX * num_row;
  ur_index.resize(URcountX);
  ur_value.resize(URcountX);

  // UR pointer
  //
  // NB ur_lastp just being used as temporary storage here
  ur_start.assign(num_row + 1, 0);
  ur_lastp.assign(num_row, 0);
  ur_space.assign(num_row, URstuffX);
  for (HighsInt k = 0; k < UcountX; k++) ur_lastp[u_pivot_lookup[u_index[k]]]++;
  for (HighsInt i = 1; i <= num_row; i++)
    ur_start[i] = ur_start[i - 1] + ur_lastp[i - 1] + URstuffX;
  ur_start.resize(num_row);

  // UR element
  //
  // NB ur_lastp initialised here!
  ur_lastp = ur_start;
  for (HighsInt i = 0; i < num_row; i++) {
    const HighsInt index = u_pivot_index[i];
    for (HighsInt k = u_start[i]; k < u_last_p[i]; k++) {
      HighsInt iRow = u_pivot_lookup[u_index[k]];
      HighsInt iPut = ur_lastp[iRow]++;
      ur_index[iPut] = index;
      ur_value[iPut] = u_value[k];
    }
  }

  // Re-factor merit
  UmeritX = num_row + (LcountX + UcountX) * 1.5;
  UtotalX = UcountX;
  if (updateMethod == kUpdateMethodPf) UmeritX = num_row + UcountX * 4;
  if (updateMethod == kUpdateMethodMpf) UmeritX = num_row + UcountX * 3;

  // Clear update buffer
  PFpivotValue.clear();
  PFpivotIndex.clear();
  PFstart.clear();
  PFstart.push_back(0);
  PFindex.clear();
  PFvalue.clear();

  if (!this->refactor_info_.use) {
    // Finally, if not calling buildFinish after refactorizing,
    // permute the base index
    iwork.assign(baseIndex, baseIndex + num_row);
    for (HighsInt i = 0; i < num_row; i++) baseIndex[permute[i]] = iwork[i];
    // Add cost of buildFinish to build_synthetic_tick
    build_synthetic_tick += num_row * 80 + (LcountX + UcountX) * 60;
  }
}

void HFactor::ftranL(HVector& rhs, const double expected_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorFtranLower, factor_timer_clock_pointer);
  if (updateMethod == kUpdateMethodApf) {
    factor_timer.start(FactorFtranLowerAPF, factor_timer_clock_pointer);
    rhs.tight();
    rhs.pack();
    ftranAPF(rhs);
    factor_timer.stop(FactorFtranLowerAPF, factor_timer_clock_pointer);
    rhs.tight();
  }

  double current_density = 1.0 * rhs.count / num_row;
  if (current_density > kHyperCancel || expected_density > kHyperFtranL) {
    factor_timer.start(FactorFtranLowerSps, factor_timer_clock_pointer);
    // Alias to RHS
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to factor L
    const HighsInt* l_start = &this->l_start[0];
    const HighsInt* l_index = this->l_index.size() > 0 ? &this->l_index[0] : NULL;
    const double* l_value = this->l_value.size() > 0 ? &this->l_value[0] : NULL;

    // Transform
    for (HighsInt i = 0; i < num_row; i++) {
      HighsInt pivotRow = l_pivot_index[i];
      const double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > kHighsTiny) {
        RHSindex[RHScount++] = pivotRow;
        const HighsInt start = l_start[i];
        const HighsInt end = l_start[i + 1];
        for (HighsInt k = start; k < end; k++)
          RHSarray[l_index[k]] -= pivotX * l_value[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    factor_timer.stop(FactorFtranLowerSps, factor_timer_clock_pointer);
  } else {
    factor_timer.start(FactorFtranLowerHyper, factor_timer_clock_pointer);
    const HighsInt* l_index = this->l_index.size() > 0 ? &this->l_index[0] : NULL;
    const double* l_value = this->l_value.size() > 0 ? &this->l_value[0] : NULL;
    solveHyper(num_row, &l_pivot_lookup[0], &l_pivot_index[0], 0, &l_start[0],
               &l_start[1], &l_index[0], &l_value[0], &rhs);
    factor_timer.stop(FactorFtranLowerHyper, factor_timer_clock_pointer);
  }
  factor_timer.stop(FactorFtranLower, factor_timer_clock_pointer);
}

void HFactor::btranL(HVector& rhs, const double expected_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorBtranLower, factor_timer_clock_pointer);
  double current_density = 1.0 * rhs.count / num_row;
  if (current_density > kHyperCancel || expected_density > kHyperBtranL) {
    // Alias to RHS
    factor_timer.start(FactorBtranLowerSps, factor_timer_clock_pointer);
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to factor L
    const HighsInt* lr_start = &this->lr_start[0];
    const HighsInt* lr_index =
        this->lr_index.size() > 0 ? &this->lr_index[0] : NULL;
    const double* lr_value = this->lr_value.size() > 0 ? &this->lr_value[0] : NULL;

    // Transform
    for (HighsInt i = num_row - 1; i >= 0; i--) {
      HighsInt pivotRow = l_pivot_index[i];
      const double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > kHighsTiny) {
        RHSindex[RHScount++] = pivotRow;
        RHSarray[pivotRow] = pivotX;
        const HighsInt start = lr_start[i];
        const HighsInt end = lr_start[i + 1];
        for (HighsInt k = start; k < end; k++)
          RHSarray[lr_index[k]] -= pivotX * lr_value[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    factor_timer.stop(FactorBtranLowerSps, factor_timer_clock_pointer);
  } else {
    factor_timer.start(FactorBtranLowerHyper, factor_timer_clock_pointer);
    const HighsInt* lr_index =
        this->lr_index.size() > 0 ? &this->lr_index[0] : NULL;
    const double* lr_value = this->lr_value.size() > 0 ? &this->lr_value[0] : NULL;
    solveHyper(num_row, &l_pivot_lookup[0], &l_pivot_index[0], 0, &lr_start[0],
               &lr_start[1], &lr_index[0], &lr_value[0], &rhs);
    factor_timer.stop(FactorBtranLowerHyper, factor_timer_clock_pointer);
  }

  if (updateMethod == kUpdateMethodApf) {
    factor_timer.start(FactorBtranLowerAPF, factor_timer_clock_pointer);
    btranAPF(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorBtranLowerAPF, factor_timer_clock_pointer);
  }
  factor_timer.stop(FactorBtranLower, factor_timer_clock_pointer);
}

void HFactor::ftranU(HVector& rhs, const double expected_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorFtranUpper, factor_timer_clock_pointer);
  // The update part
  if (updateMethod == kUpdateMethodFt) {
    factor_timer.start(FactorFtranUpperFT, factor_timer_clock_pointer);
    //    const double current_density = 1.0 * rhs.count / num_row;
    ftranFT(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorFtranUpperFT, factor_timer_clock_pointer);
  }
  if (updateMethod == kUpdateMethodMpf) {
    factor_timer.start(FactorFtranUpperMPF, factor_timer_clock_pointer);
    ftranMPF(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorFtranUpperMPF, factor_timer_clock_pointer);
  }

  // The regular part
  const double current_density = 1.0 * rhs.count / num_row;
  if (current_density > kHyperCancel || expected_density > kHyperFtranU) {
    const bool report_ftran_upper_sparse =
        false;  // current_density < kHyperCancel;
    HighsInt use_clock;
    if (current_density < 0.1)
      use_clock = FactorFtranUpperSps2;
    else if (current_density < 0.5)
      use_clock = FactorFtranUpperSps1;
    else
      use_clock = FactorFtranUpperSps0;
    factor_timer.start(use_clock, factor_timer_clock_pointer);
    // Alias to non constant
    double RHS_synthetic_tick = 0;
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to the factor
    const HighsInt* u_start = &this->u_start[0];
    const HighsInt* Uend = &this->u_last_p[0];
    const HighsInt* u_index = this->u_index.size() > 0 ? &this->u_index[0] : NULL;
    const double* u_value = this->u_value.size() > 0 ? &this->u_value[0] : NULL;

    // Transform
    HighsInt UpivotCount = u_pivot_index.size();
    for (HighsInt iLogic = UpivotCount - 1; iLogic >= 0; iLogic--) {
      // Skip void
      if (u_pivot_index[iLogic] == -1) continue;

      // Normal part
      const HighsInt pivotRow = u_pivot_index[iLogic];
      double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > kHighsTiny) {
        pivotX /= u_pivot_value[iLogic];
        RHSindex[RHScount++] = pivotRow;
        RHSarray[pivotRow] = pivotX;
        const HighsInt start = u_start[iLogic];
        const HighsInt end = Uend[iLogic];
        if (iLogic >= num_row) {
          RHS_synthetic_tick += (end - start);
        }
        for (HighsInt k = start; k < end; k++)
          RHSarray[u_index[k]] -= pivotX * u_value[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    rhs.synthetic_tick += RHS_synthetic_tick * 15 + (UpivotCount - num_row) * 10;
    factor_timer.stop(use_clock, factor_timer_clock_pointer);
    if (report_ftran_upper_sparse) {
      const double final_density = 1.0 * rhs.count / num_row;
      printf(
          "FactorFtranUpperSps: expected_density = %10.4g; current_density = "
          "%10.4g; final_density = %10.4g\n",
          expected_density, current_density, final_density);
    }
  } else {
    HighsInt use_clock = -1;
    if (current_density < 5e-6)
      use_clock = FactorFtranUpperHyper5;
    else if (current_density < 1e-5)
      use_clock = FactorFtranUpperHyper4;
    else if (current_density < 1e-4)
      use_clock = FactorFtranUpperHyper3;
    else if (current_density < 1e-3)
      use_clock = FactorFtranUpperHyper2;
    else if (current_density < 1e-2)
      use_clock = FactorFtranUpperHyper1;
    else
      use_clock = FactorFtranUpperHyper0;
    factor_timer.start(use_clock, factor_timer_clock_pointer);
    const HighsInt* u_index = this->u_index.size() > 0 ? &this->u_index[0] : NULL;
    const double* u_value = this->u_value.size() > 0 ? &this->u_value[0] : NULL;
    solveHyper(num_row, &u_pivot_lookup[0], &u_pivot_index[0], &u_pivot_value[0],
               &u_start[0], &u_last_p[0], &u_index[0], &u_value[0], &rhs);
    factor_timer.stop(use_clock, factor_timer_clock_pointer);
  }
  if (updateMethod == kUpdateMethodPf) {
    factor_timer.start(FactorFtranUpperPF, factor_timer_clock_pointer);
    ftranPF(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorFtranUpperPF, factor_timer_clock_pointer);
  }
  factor_timer.stop(FactorFtranUpper, factor_timer_clock_pointer);
}

void HFactor::btranU(HVector& rhs, const double expected_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorBtranUpper, factor_timer_clock_pointer);
  if (updateMethod == kUpdateMethodPf) {
    factor_timer.start(FactorBtranUpperPF, factor_timer_clock_pointer);
    btranPF(rhs);
    factor_timer.stop(FactorBtranUpperPF, factor_timer_clock_pointer);
  }

  // The regular part
  double current_density = 1.0 * rhs.count / num_row;
  if (current_density > kHyperCancel || expected_density > kHyperBtranU) {
    factor_timer.start(FactorBtranUpperSps, factor_timer_clock_pointer);
    // Alias to non constant
    double RHS_synthetic_tick = 0;
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to the factor
    const HighsInt* ur_start = &this->ur_start[0];
    const HighsInt* URend = &this->ur_lastp[0];
    const HighsInt* ur_index = &this->ur_index[0];
    const double* ur_value = &this->ur_value[0];

    // Transform
    HighsInt UpivotCount = u_pivot_index.size();
    for (HighsInt iLogic = 0; iLogic < UpivotCount; iLogic++) {
      // Skip void
      if (u_pivot_index[iLogic] == -1) continue;

      // Normal part
      const HighsInt pivotRow = u_pivot_index[iLogic];
      double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > kHighsTiny) {
        pivotX /= u_pivot_value[iLogic];
        RHSindex[RHScount++] = pivotRow;
        RHSarray[pivotRow] = pivotX;
        const HighsInt start = ur_start[iLogic];
        const HighsInt end = URend[iLogic];
        if (iLogic >= num_row) {
          RHS_synthetic_tick += (end - start);
        }
        for (HighsInt k = start; k < end; k++)
          RHSarray[ur_index[k]] -= pivotX * ur_value[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    rhs.synthetic_tick += RHS_synthetic_tick * 15 + (UpivotCount - num_row) * 10;
    factor_timer.stop(FactorBtranUpperSps, factor_timer_clock_pointer);
  } else {
    factor_timer.start(FactorBtranUpperHyper, factor_timer_clock_pointer);
    solveHyper(num_row, &u_pivot_lookup[0], &u_pivot_index[0], &u_pivot_value[0],
               &ur_start[0], &ur_lastp[0], &ur_index[0], &ur_value[0], &rhs);
    factor_timer.stop(FactorBtranUpperHyper, factor_timer_clock_pointer);
  }

  // The update part
  if (updateMethod == kUpdateMethodFt) {
    factor_timer.start(FactorBtranUpperFT, factor_timer_clock_pointer);
    rhs.tight();
    rhs.pack();
    //    const double current_density = 1.0 * rhs.count / num_row;
    btranFT(rhs);
    rhs.tight();
    factor_timer.stop(FactorBtranUpperFT, factor_timer_clock_pointer);
  }
  if (updateMethod == kUpdateMethodMpf) {
    factor_timer.start(FactorBtranUpperMPF, factor_timer_clock_pointer);
    rhs.tight();
    rhs.pack();
    btranMPF(rhs);
    rhs.tight();
    factor_timer.stop(FactorBtranUpperMPF, factor_timer_clock_pointer);
  }
  factor_timer.stop(FactorBtranUpper, factor_timer_clock_pointer);
}

void HFactor::ftranFT(HVector& vector) const {
  // Alias to PF buffer
  const HighsInt PFpivotCount = PFpivotIndex.size();
  HighsInt* PFpivotIndex = NULL;
  if (this->PFpivotIndex.size() > 0)
    PFpivotIndex = (HighsInt*)&this->PFpivotIndex[0];

  const HighsInt* PFstart = this->PFstart.size() > 0 ? &this->PFstart[0] : NULL;
  const HighsInt* PFindex = this->PFindex.size() > 0 ? &this->PFindex[0] : NULL;
  const double* PFvalue = this->PFvalue.size() > 0 ? &this->PFvalue[0] : NULL;

  // Alias to non constant
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Forwardly apply row ETA
  for (HighsInt i = 0; i < PFpivotCount; i++) {
    HighsInt iRow = PFpivotIndex[i];
    double value0 = RHSarray[iRow];
    double value1 = value0;
    const HighsInt start = PFstart[i];
    const HighsInt end = PFstart[i + 1];
    for (HighsInt k = start; k < end; k++)
      value1 -= RHSarray[PFindex[k]] * PFvalue[k];
    // This would skip the situation where they are both zeros
    if (value0 || value1) {
      if (value0 == 0) RHSindex[RHScount++] = iRow;
      RHSarray[iRow] = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
    }
  }

  // Save count back
  vector.count = RHScount;
  vector.synthetic_tick += PFpivotCount * 20 + PFstart[PFpivotCount] * 5;
  if (PFstart[PFpivotCount] / (PFpivotCount + 1) < 5) {
    vector.synthetic_tick += PFstart[PFpivotCount] * 5;
  }
}

void HFactor::btranFT(HVector& vector) const {
  // Alias to PF buffer
  const HighsInt PFpivotCount = PFpivotIndex.size();
  const HighsInt* PFpivotIndex =
      this->PFpivotIndex.size() > 0 ? &this->PFpivotIndex[0] : NULL;
  const HighsInt* PFstart = this->PFstart.size() > 0 ? &this->PFstart[0] : NULL;
  const HighsInt* PFindex = this->PFindex.size() > 0 ? &this->PFindex[0] : NULL;
  const double* PFvalue = this->PFvalue.size() > 0 ? &this->PFvalue[0] : NULL;

  // Alias to non constant
  double RHS_synthetic_tick = 0;
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Backwardly apply row ETA
  for (HighsInt i = PFpivotCount - 1; i >= 0; i--) {
    HighsInt pivotRow = PFpivotIndex[i];
    double pivotX = RHSarray[pivotRow];
    if (pivotX) {
      const HighsInt start = PFstart[i];
      const HighsInt end = PFstart[i + 1];
      RHS_synthetic_tick += (end - start);
      for (HighsInt k = start; k < end; k++) {
        HighsInt iRow = PFindex[k];
        double value0 = RHSarray[iRow];
        double value1 = value0 - pivotX * PFvalue[k];
        if (value0 == 0) RHSindex[RHScount++] = iRow;
        RHSarray[iRow] = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
      }
    }
  }

  vector.synthetic_tick += RHS_synthetic_tick * 15 + PFpivotCount * 10;

  // Save count back
  vector.count = RHScount;
}

void HFactor::ftranPF(HVector& vector) const {
  // Alias to PF buffer
  const HighsInt PFpivotCount = PFpivotIndex.size();
  const HighsInt* PFpivotIndex = &this->PFpivotIndex[0];
  const double* PFpivotValue = &this->PFpivotValue[0];
  const HighsInt* PFstart = &this->PFstart[0];
  const HighsInt* PFindex = &this->PFindex[0];
  const double* PFvalue = &this->PFvalue[0];

  // Alias to non constant
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Forwardly
  for (HighsInt i = 0; i < PFpivotCount; i++) {
    HighsInt pivotRow = PFpivotIndex[i];
    double pivotX = RHSarray[pivotRow];
    if (fabs(pivotX) > kHighsTiny) {
      pivotX /= PFpivotValue[i];
      RHSarray[pivotRow] = pivotX;
      for (HighsInt k = PFstart[i]; k < PFstart[i + 1]; k++) {
        const HighsInt index = PFindex[k];
        const double value0 = RHSarray[index];
        const double value1 = value0 - pivotX * PFvalue[k];
        if (value0 == 0) RHSindex[RHScount++] = index;
        RHSarray[index] = (fabs(value1) < kHighsTiny) ? kHighsZero : value1;
      }
    }
  }

  // Save count
  vector.count = RHScount;
}

void HFactor::btranPF(HVector& vector) const {
  // Alias to PF buffer
  const HighsInt PFpivotCount = PFpivotIndex.size();
  const HighsInt* PFpivotIndex = &this->PFpivotIndex[0];
  const double* PFpivotValue = &this->PFpivotValue[0];
  const HighsInt* PFstart = &this->PFstart[0];
  const HighsInt* PFindex = &this->PFindex[0];
  const double* PFvalue = &this->PFvalue[0];

  // Alias to non constant
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Backwardly
  for (HighsInt i = PFpivotCount - 1; i >= 0; i--) {
    HighsInt pivotRow = PFpivotIndex[i];
    double pivotX = RHSarray[pivotRow];
    for (HighsInt k = PFstart[i]; k < PFstart[i + 1]; k++)
      pivotX -= PFvalue[k] * RHSarray[PFindex[k]];
    pivotX /= PFpivotValue[i];

    if (RHSarray[pivotRow] == 0) RHSindex[RHScount++] = pivotRow;
    RHSarray[pivotRow] = (fabs(pivotX) < kHighsTiny) ? 1e-100 : pivotX;
  }

  // Save count
  vector.count = RHScount;
}

void HFactor::ftranMPF(HVector& vector) const {
  // Alias to non constant
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Forwardly
  HighsInt PFpivotCount = PFpivotValue.size();
  for (HighsInt i = 0; i < PFpivotCount; i++) {
    solveMatrixT(PFstart[i * 2 + 1], PFstart[i * 2 + 2], PFstart[i * 2],
                 PFstart[i * 2 + 1], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                 &RHScount, RHSindex, RHSarray);
  }

  // Remove cancellation
  vector.count = RHScount;
}

void HFactor::btranMPF(HVector& vector) const {
  // Alias to non constant
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Backwardly
  for (HighsInt i = PFpivotValue.size() - 1; i >= 0; i--) {
    solveMatrixT(PFstart[i * 2], PFstart[i * 2 + 1], PFstart[i * 2 + 1],
                 PFstart[i * 2 + 2], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                 &RHScount, RHSindex, RHSarray);
  }

  // Remove cancellation
  vector.count = RHScount;
}

void HFactor::ftranAPF(HVector& vector) const {
  // Alias to non constant
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Backwardly
  HighsInt PFpivotCount = PFpivotValue.size();
  for (HighsInt i = PFpivotCount - 1; i >= 0; i--) {
    solveMatrixT(PFstart[i * 2 + 1], PFstart[i * 2 + 2], PFstart[i * 2],
                 PFstart[i * 2 + 1], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                 &RHScount, RHSindex, RHSarray);
  }

  // Remove cancellation
  vector.count = RHScount;
}

void HFactor::btranAPF(HVector& vector) const {
  // Alias to non constant
  HighsInt RHScount = vector.count;
  HighsInt* RHSindex = &vector.index[0];
  double* RHSarray = &vector.array[0];

  // Forwardly
  HighsInt PFpivotCount = PFpivotValue.size();
  for (HighsInt i = 0; i < PFpivotCount; i++) {
    solveMatrixT(PFstart[i * 2], PFstart[i * 2 + 1], PFstart[i * 2 + 1],
                 PFstart[i * 2 + 2], &PFindex[0], &PFvalue[0], PFpivotValue[i],
                 &RHScount, RHSindex, RHSarray);
  }
  vector.count = RHScount;
}

void HFactor::updateCFT(HVector* aq, HVector* ep, HighsInt* iRow
                        //, HighsInt* hint
) {
  /*
   * In the major update loop, the prefix
   *
   * c(p) = current working pivot
   * p(p) = previous pivot  (0 =< pp < cp)
   */

  HighsInt numUpdate = 0;
  for (HVector* vec = aq; vec != 0; vec = vec->next) numUpdate++;

  HVector** aqWork = new HVector*[numUpdate];
  HVector** epWork = new HVector*[numUpdate];

  for (HighsInt i = 0; i < numUpdate; i++) {
    aqWork[i] = aq;
    epWork[i] = ep;
    aq = aq->next;
    ep = ep->next;
  }

  // Pivot related buffers
  HighsInt PFnp0 = PFpivotIndex.size();
  HighsInt* pLogic = new HighsInt[numUpdate];
  double* pValue = new double[numUpdate];
  double* pAlpha = new double[numUpdate];
  for (HighsInt cp = 0; cp < numUpdate; cp++) {
    HighsInt cRow = iRow[cp];
    HighsInt iLogic = u_pivot_lookup[cRow];
    pLogic[cp] = iLogic;
    pValue[cp] = u_pivot_value[iLogic];
    pAlpha[cp] = aqWork[cp]->array[cRow];
  }

  // Temporary U pointers
  HighsInt* Tstart = new HighsInt[numUpdate + 1];
  double* Tpivot = new double[numUpdate];
  Tstart[0] = u_index.size();

  // Logically sorted previous row_ep
  vector<pair<HighsInt, int> > sorted_pp;

  // Major update loop
  for (HighsInt cp = 0; cp < numUpdate; cp++) {
    // 1. Expand partial FTRAN result to buffer
    iwork.clear();
    for (HighsInt i = 0; i < aqWork[cp]->packCount; i++) {
      HighsInt index = aqWork[cp]->packIndex[i];
      double value = aqWork[cp]->packValue[i];
      iwork.push_back(index);
      dwork[index] = value;
    }

    // 2. Update partial FTRAN result by recent FT matrix
    for (HighsInt pp = 0; pp < cp; pp++) {
      HighsInt pRow = iRow[pp];
      double value = dwork[pRow];
      HighsInt PFpp = pp + PFnp0;
      for (HighsInt i = PFstart[PFpp]; i < PFstart[PFpp + 1]; i++)
        value -= dwork[PFindex[i]] * PFvalue[i];
      iwork.push_back(pRow);  // OK to duplicate
      dwork[pRow] = value;
    }

    // 3. Store the partial FTRAN result to matirx U
    double ppaq = dwork[iRow[cp]];  // pivot of the partial aq
    dwork[iRow[cp]] = 0;
    HighsInt UcountX = Tstart[cp];
    HighsInt u_startX = UcountX;
    for (unsigned i = 0; i < iwork.size(); i++) {
      HighsInt index = iwork[i];
      double value = dwork[index];
      dwork[index] = 0;  // This effectively removes all duplication
      if (fabs(value) > kHighsTiny) {
        u_index.push_back(index);
        u_value.push_back(value);
      }
    }
    UcountX = u_index.size();
    Tstart[cp + 1] = UcountX;
    Tpivot[cp] = pValue[cp] * pAlpha[cp];

    // 4. Expand partial BTRAN result to buffer
    iwork.clear();
    for (HighsInt i = 0; i < epWork[cp]->packCount; i++) {
      HighsInt index = epWork[cp]->packIndex[i];
      double value = epWork[cp]->packValue[i];
      iwork.push_back(index);
      dwork[index] = value;
    }

    // 5. Delete logical later rows (in logical order)
    for (HighsInt isort = 0; isort < cp; isort++) {
      HighsInt pp = sorted_pp[isort].second;
      HighsInt pRow = iRow[pp];
      double multiplier = -pValue[pp] * dwork[pRow];
      if (fabs(dwork[pRow]) > kHighsTiny) {
        for (HighsInt i = 0; i < epWork[pp]->packCount; i++) {
          HighsInt index = epWork[pp]->packIndex[i];
          double value = epWork[pp]->packValue[i];
          iwork.push_back(index);
          dwork[index] += value * multiplier;
        }
      }
      dwork[pRow] = 0;  // Force to be 0
    }

    // 6. Update partial BTRAN result by recent U columns
    for (HighsInt pp = 0; pp < cp; pp++) {
      HighsInt kpivot = iRow[pp];
      double value = dwork[kpivot];
      for (HighsInt k = Tstart[pp]; k < Tstart[pp + 1]; k++)
        value -= dwork[u_index[k]] * u_value[k];
      value /= Tpivot[pp];
      iwork.push_back(kpivot);
      dwork[kpivot] = value;  // Again OK to duplicate
    }

    // 6.x compute current alpha
    double thex = 0;
    for (HighsInt k = u_startX; k < UcountX; k++) {
      HighsInt index = u_index[k];
      double value = u_value[k];
      thex += dwork[index] * value;
    }
    Tpivot[cp] = ppaq + thex * pValue[cp];

    // 7. Store BTRAN result to FT elimination, update logic helper
    dwork[iRow[cp]] = 0;
    double pivotX = -pValue[cp];
    for (unsigned i = 0; i < iwork.size(); i++) {
      HighsInt index = iwork[i];
      double value = dwork[index];
      dwork[index] = 0;
      if (fabs(value) > kHighsTiny) {
        PFindex.push_back(index);
        PFvalue.push_back(value * pivotX);
      }
    }
    PFpivotIndex.push_back(iRow[cp]);
    UtotalX += PFindex.size() - PFstart.back();
    PFstart.push_back(PFindex.size());

    // 8. Update the sorted ep
    sorted_pp.push_back(make_pair(pLogic[cp], cp));
    pdqsort(sorted_pp.begin(), sorted_pp.end());
  }

  // Now modify the U matrix
  for (HighsInt cp = 0; cp < numUpdate; cp++) {
    // 1. Delete pivotal row from U
    HighsInt cIndex = iRow[cp];
    HighsInt cLogic = pLogic[cp];
    UtotalX -= ur_lastp[cLogic] - ur_start[cLogic];
    for (HighsInt k = ur_start[cLogic]; k < ur_lastp[cLogic]; k++) {
      // Find the pivotal position
      HighsInt iLogic = u_pivot_lookup[ur_index[k]];
      HighsInt iFind = u_start[iLogic];
      HighsInt iLast = --u_last_p[iLogic];
      for (; iFind <= iLast; iFind++)
        if (u_index[iFind] == cIndex) break;
      // Put last to find, and delete last
      u_index[iFind] = u_index[iLast];
      u_value[iFind] = u_value[iLast];
    }

    // 2. Delete pivotal column from UR
    UtotalX -= u_last_p[cLogic] - u_start[cLogic];
    for (HighsInt k = u_start[cLogic]; k < u_last_p[cLogic]; k++) {
      // Find the pivotal position
      HighsInt iLogic = u_pivot_lookup[u_index[k]];
      HighsInt iFind = ur_start[iLogic];
      HighsInt iLast = --ur_lastp[iLogic];
      for (; iFind <= iLast; iFind++)
        if (ur_index[iFind] == cIndex) break;
      // Put last to find, and delete last
      ur_space[iLogic]++;
      ur_index[iFind] = ur_index[iLast];
      ur_value[iFind] = ur_value[iLast];
    }

    // 3. Insert the (stored) partial FTRAN to the row matrix
    HighsInt u_startX = Tstart[cp];
    HighsInt UendX = Tstart[cp + 1];
    UtotalX += UendX - u_startX;
    // Store column as UR elements
    for (HighsInt k = u_startX; k < UendX; k++) {
      // Which ETA file
      HighsInt iLogic = u_pivot_lookup[u_index[k]];

      // Move row to the end if necessary
      if (ur_space[iLogic] == 0) {
        // Make pointers
        HighsInt row_start = ur_start[iLogic];
        HighsInt row_count = ur_lastp[iLogic] - row_start;
        HighsInt new_start = ur_index.size();
        HighsInt new_space = row_count * 1.1 + 5;

        // Check matrix UR
        ur_index.resize(new_start + new_space);
        ur_value.resize(new_start + new_space);

        // Move elements
        HighsInt iFrom = row_start;
        HighsInt iEnd = row_start + row_count;
        HighsInt iTo = new_start;
        copy(&ur_index[iFrom], &ur_index[iEnd], &ur_index[iTo]);
        copy(&ur_value[iFrom], &ur_value[iEnd], &ur_value[iTo]);

        // Save new pointers
        ur_start[iLogic] = new_start;
        ur_lastp[iLogic] = new_start + row_count;
        ur_space[iLogic] = new_space - row_count;
      }

      // Put into the next available space
      ur_space[iLogic]--;
      HighsInt iPut = ur_lastp[iLogic]++;
      ur_index[iPut] = cIndex;
      ur_value[iPut] = u_value[k];
    }

    // 4. Save pointers
    u_start.push_back(u_startX);
    u_last_p.push_back(UendX);

    ur_start.push_back(ur_start[cLogic]);
    ur_lastp.push_back(ur_start[cLogic]);
    ur_space.push_back(ur_space[cLogic] + ur_lastp[cLogic] - ur_start[cLogic]);

    u_pivot_lookup[cIndex] = u_pivot_index.size();
    u_pivot_index[cLogic] = -1;
    u_pivot_index.push_back(cIndex);
    u_pivot_value.push_back(Tpivot[cp]);
  }

  //    // See if we want refactor
  //    if (UtotalX > UmeritX && PFpivotIndex.size() > 100)
  //        *hint = 1;
  delete[] aqWork;
  delete[] epWork;
  delete[] pLogic;
  delete[] pValue;
  delete[] pAlpha;
  delete[] Tstart;
  delete[] Tpivot;
}

void HFactor::updateFT(HVector* aq, HVector* ep, HighsInt iRow
                       //, HighsInt* hint
) {
  // Store pivot
  HighsInt pLogic = u_pivot_lookup[iRow];
  double pivot = u_pivot_value[pLogic];
  double alpha = aq->array[iRow];
  u_pivot_index[pLogic] = -1;

  // Delete pivotal row from U
  for (HighsInt k = ur_start[pLogic]; k < ur_lastp[pLogic]; k++) {
    // Find the pivotal position
    HighsInt iLogic = u_pivot_lookup[ur_index[k]];
    HighsInt iFind = u_start[iLogic];
    HighsInt iLast = --u_last_p[iLogic];
    for (; iFind <= iLast; iFind++)
      if (u_index[iFind] == iRow) break;
    // Put last to find, and delete last
    u_index[iFind] = u_index[iLast];
    u_value[iFind] = u_value[iLast];
  }

  // Delete pivotal column from UR
  for (HighsInt k = u_start[pLogic]; k < u_last_p[pLogic]; k++) {
    // Find the pivotal position
    HighsInt iLogic = u_pivot_lookup[u_index[k]];
    HighsInt iFind = ur_start[iLogic];
    HighsInt iLast = --ur_lastp[iLogic];
    for (; iFind <= iLast; iFind++)
      if (ur_index[iFind] == iRow) break;
    // Put last to find, and delete last
    ur_space[iLogic]++;
    ur_index[iFind] = ur_index[iLast];
    ur_value[iFind] = ur_value[iLast];
  }

  // Store column to U
  u_start.push_back(u_index.size());
  for (HighsInt i = 0; i < aq->packCount; i++)
    if (aq->packIndex[i] != iRow) {
      u_index.push_back(aq->packIndex[i]);
      u_value.push_back(aq->packValue[i]);
    }
  u_last_p.push_back(u_index.size());
  HighsInt u_startX = u_start.back();
  HighsInt UendX = u_last_p.back();
  UtotalX += UendX - u_startX + 1;

  // Store column as UR elements
  for (HighsInt k = u_startX; k < UendX; k++) {
    // Which ETA file
    HighsInt iLogic = u_pivot_lookup[u_index[k]];

    // Move row to the end if necessary
    if (ur_space[iLogic] == 0) {
      // Make pointers
      HighsInt row_start = ur_start[iLogic];
      HighsInt row_count = ur_lastp[iLogic] - row_start;
      HighsInt new_start = ur_index.size();
      HighsInt new_space = row_count * 1.1 + 5;

      // Check matrix UR
      ur_index.resize(new_start + new_space);
      ur_value.resize(new_start + new_space);

      // Move elements
      HighsInt iFrom = row_start;
      HighsInt iEnd = row_start + row_count;
      HighsInt iTo = new_start;
      copy(&ur_index[iFrom], &ur_index[iEnd], &ur_index[iTo]);
      copy(&ur_value[iFrom], &ur_value[iEnd], &ur_value[iTo]);

      // Save new pointers
      ur_start[iLogic] = new_start;
      ur_lastp[iLogic] = new_start + row_count;
      ur_space[iLogic] = new_space - row_count;
    }

    // Put into the next available space
    ur_space[iLogic]--;
    HighsInt iPut = ur_lastp[iLogic]++;
    ur_index[iPut] = iRow;
    ur_value[iPut] = u_value[k];
  }

  // Store UR pointers
  ur_start.push_back(ur_start[pLogic]);
  ur_lastp.push_back(ur_start[pLogic]);
  ur_space.push_back(ur_space[pLogic] + ur_lastp[pLogic] - ur_start[pLogic]);

  // Update pivot count
  u_pivot_lookup[iRow] = u_pivot_index.size();
  u_pivot_index.push_back(iRow);
  u_pivot_value.push_back(pivot * alpha);

  // Store row_ep as R matrix
  for (HighsInt i = 0; i < ep->packCount; i++) {
    if (ep->packIndex[i] != iRow) {
      PFindex.push_back(ep->packIndex[i]);
      PFvalue.push_back(-ep->packValue[i] * pivot);
    }
  }
  UtotalX += PFindex.size() - PFstart.back();

  // Store R matrix pivot
  PFpivotIndex.push_back(iRow);
  PFstart.push_back(PFindex.size());

  // Update total countX
  UtotalX -= u_last_p[pLogic] - u_start[pLogic];
  UtotalX -= ur_lastp[pLogic] - ur_start[pLogic];

  //    // See if we want refactor
  //    if (UtotalX > UmeritX && PFpivotIndex.size() > 100)
  //        *hint = 1;
}

void HFactor::updatePF(HVector* aq, HighsInt iRow, HighsInt* hint) {
  // Check space
  const HighsInt columnCount = aq->packCount;
  const HighsInt* variable_index = &aq->packIndex[0];
  const double* columnArray = &aq->packValue[0];

  // Copy the pivotal column
  for (HighsInt i = 0; i < columnCount; i++) {
    HighsInt index = variable_index[i];
    double value = columnArray[i];
    if (index != iRow) {
      PFindex.push_back(index);
      PFvalue.push_back(value);
    }
  }

  // Save pivot
  PFpivotIndex.push_back(iRow);
  PFpivotValue.push_back(aq->array[iRow]);
  PFstart.push_back(PFindex.size());

  // Check refactor
  UtotalX += aq->packCount;
  if (UtotalX > UmeritX) *hint = 1;
}

void HFactor::updateMPF(HVector* aq, HVector* ep, HighsInt iRow,
                        HighsInt* hint) {
  // Store elements
  for (HighsInt i = 0; i < aq->packCount; i++) {
    PFindex.push_back(aq->packIndex[i]);
    PFvalue.push_back(aq->packValue[i]);
  }
  HighsInt pLogic = u_pivot_lookup[iRow];
  HighsInt u_startX = u_start[pLogic];
  HighsInt UendX = u_start[pLogic + 1];
  for (HighsInt k = u_startX; k < UendX; k++) {
    PFindex.push_back(u_index[k]);
    PFvalue.push_back(-u_value[k]);
  }
  PFindex.push_back(iRow);
  PFvalue.push_back(-u_pivot_value[pLogic]);
  PFstart.push_back(PFindex.size());

  for (HighsInt i = 0; i < ep->packCount; i++) {
    PFindex.push_back(ep->packIndex[i]);
    PFvalue.push_back(ep->packValue[i]);
  }
  PFstart.push_back(PFindex.size());

  // Store pivot
  PFpivotValue.push_back(aq->array[iRow]);

  // Refactor or not
  UtotalX += aq->packCount + ep->packCount;
  if (UtotalX > UmeritX) *hint = 1;
}

void HFactor::updateAPF(HVector* aq, HVector* ep, HighsInt iRow
                        //, HighsInt* hint
) {
  // Store elements
  for (HighsInt i = 0; i < aq->packCount; i++) {
    PFindex.push_back(aq->packIndex[i]);
    PFvalue.push_back(aq->packValue[i]);
  }

  HighsInt variable_out = baseIndex[iRow];
  if (variable_out >= num_col) {
    PFindex.push_back(variable_out - num_col);
    PFvalue.push_back(-1);
  } else {
    for (HighsInt k = a_start[variable_out]; k < a_start[variable_out + 1]; k++) {
      PFindex.push_back(a_index[k]);
      PFvalue.push_back(-a_value[k]);
    }
  }
  PFstart.push_back(PFindex.size());

  for (HighsInt i = 0; i < ep->packCount; i++) {
    PFindex.push_back(ep->packIndex[i]);
    PFvalue.push_back(ep->packValue[i]);
  }
  PFstart.push_back(PFindex.size());

  // Store pivot
  PFpivotValue.push_back(aq->array[iRow]);
}
