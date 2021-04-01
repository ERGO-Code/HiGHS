/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HFactor.cpp
 * @brief Types of solution classes
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HFactor.h"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "lp_data/HConst.h"
//#include "io/HighsIO.h"
#include "simplex/FactorTimer.h"
#include "simplex/HFactorDebug.h"
#include "simplex/HVector.h"
#include "util/HighsTimer.h"

using std::copy;
using std::fill_n;
using std::make_pair;
using std::pair;
using std::vector;

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
  if (fabs(pivotX) > HIGHS_CONST_TINY) {
    HighsInt workCount = *RHScount;

    pivotX /= Tpivot;
    for (HighsInt k = Ystart; k < Yend; k++) {
      const HighsInt index = Tindex[k];
      const double value0 = RHSarray[index];
      const double value1 = value0 - pivotX * Tvalue[k];
      if (value0 == 0) RHSindex[workCount++] = index;
      RHSarray[index] =
          (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
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

  rhs->syntheticTick += countPivot * 20 + countEntry * 10;

  // Solve with list
  if (HpivotValue == 0) {
    RHScount = 0;
    for (HighsInt iList = listCount - 1; iList >= 0; iList--) {
      HighsInt i = listIndex[iList];
      listMark[i] = 0;
      HighsInt pivotRow = HpivotIndex[i];
      double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > HIGHS_CONST_TINY) {
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
      if (fabs(pivotX) > HIGHS_CONST_TINY) {
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

void HFactor::setup(HighsInt numCol_, HighsInt numRow_, const HighsInt* Astart_,
                    const HighsInt* Aindex_, const double* Avalue_,
                    HighsInt* baseIndex_, double pivot_threshold_,
                    double pivot_tolerance_, HighsInt highs_debug_level_,
                    bool output_flag_, FILE* log_file_stream_,
                    bool log_to_console_, int log_dev_level_,
                    const bool use_original_HFactor_logic_,
                    const HighsInt updateMethod_) {
  // Copy Problem size and (pointer to) coefficient matrix
  numRow = numRow_;
  numCol = numCol_;
  Astart = Astart_;
  Aindex = Aindex_;
  Avalue = Avalue_;
  baseIndex = baseIndex_;
  pivot_threshold =
      max(min_pivot_threshold, min(pivot_threshold_, max_pivot_threshold));
  pivot_tolerance =
      max(min_pivot_tolerance, min(pivot_tolerance_, max_pivot_tolerance));
  highs_debug_level = highs_debug_level_;
  log_options.log_file_stream = log_file_stream_;
  log_data = decltype(log_data)(new std::tuple<bool, bool, HighsInt>(
      output_flag_, log_to_console_, log_dev_level_));
  log_options.output_flag = &std::get<0>(*log_data);
  log_options.log_to_console = &std::get<1>(*log_data);
  log_options.log_dev_level = &std::get<2>(*log_data);
  use_original_HFactor_logic = use_original_HFactor_logic_;
  updateMethod = updateMethod_;

  // Allocate for working buffer
  iwork.reserve(numRow * 2);
  dwork.assign(numRow, 0);

  // Find Basis matrix limit size
  HighsInt BlimitX = 0;
  iwork.assign(numRow + 1, 0);
  for (HighsInt i = 0; i < numCol; i++) iwork[Astart[i + 1] - Astart[i]]++;
  for (HighsInt i = numRow, counted = 0; i >= 0 && counted < numRow; i--)
    BlimitX += i * iwork[i], counted += iwork[i];
  BlimitX += numRow;

  // Allocate space for basis matrix, L, U factor and Update buffer
  Bstart.resize(numRow + 1, 0);
  Bindex.resize(BlimitX);
  Bvalue.resize(BlimitX);

  // Allocate space for pivot records
  permute.resize(numRow);

  // Allocate space for Markowitz matrices
  MCstart.resize(numRow);
  MCcountA.resize(numRow);
  MCcountN.resize(numRow);
  MCspace.resize(numRow);
  MCminpivot.resize(numRow);
  MCindex.resize(BlimitX * 2);
  MCvalue.resize(BlimitX * 2);

  MRstart.resize(numRow);
  MRcount.resize(numRow);
  MRspace.resize(numRow);
  MRcountb4.resize(numRow);
  MRindex.resize(BlimitX * 2);

  mwz_column_mark.assign(numRow, 0);
  mwz_column_index.resize(numRow);
  mwz_column_array.assign(numRow, 0);

  // Allocate space for count-link-list
  clinkFirst.assign(numRow + 1, -1);
  clinkNext.resize(numRow);
  clinkLast.resize(numRow);

  rlinkFirst.assign(numRow + 1, -1);
  rlinkNext.resize(numRow);
  rlinkLast.resize(numRow);

  // Allocate space for L factor
  LpivotLookup.resize(numRow);
  LpivotIndex.reserve(numRow);
  Lstart.reserve(numRow + 1);
  Lindex.reserve(BlimitX * 3);
  Lvalue.reserve(BlimitX * 3);

  LRstart.reserve(numRow + 1);
  LRindex.reserve(BlimitX * 3);
  LRvalue.reserve(BlimitX * 3);

  // Allocate space for U factor
  UpivotLookup.resize(numRow);
  UpivotIndex.reserve(numRow + 1000);
  UpivotValue.reserve(numRow + 1000);

  Ustart.reserve(numRow + 1000 + 1);
  Ulastp.reserve(numRow + 1000);
  Uindex.reserve(BlimitX * 3);
  Uvalue.reserve(BlimitX * 3);

  URstart.reserve(numRow + 1000 + 1);
  URlastp.reserve(numRow + 1000);
  URspace.reserve(numRow + 1000);
  URindex.reserve(BlimitX * 3);
  URvalue.reserve(BlimitX * 3);

  // Allocate spaces for Update buffer
  PFpivotValue.reserve(1000);
  PFpivotIndex.reserve(1000);
  PFstart.reserve(2000 + 1);
  PFindex.reserve(BlimitX * 4);
  PFvalue.reserve(BlimitX * 4);
}

HighsInt HFactor::build(HighsTimerClock* factor_timer_clock_pointer) {
  FactorTimer factor_timer;
  factor_timer.start(FactorInvert, factor_timer_clock_pointer);
  build_syntheticTick = 0;
  factor_timer.start(FactorInvertSimple, factor_timer_clock_pointer);
  // Build the L, U factor
  buildSimple();
  factor_timer.stop(FactorInvertSimple, factor_timer_clock_pointer);
  factor_timer.start(FactorInvertKernel, factor_timer_clock_pointer);
  rank_deficiency = buildKernel();
  factor_timer.stop(FactorInvertKernel, factor_timer_clock_pointer);
  if (rank_deficiency) {
    factor_timer.start(FactorInvertDeficient, factor_timer_clock_pointer);
    highsLogUser(log_options, HighsLogType::WARNING,
                 "Rank deficiency of %d identified in basis matrix\n",
                 rank_deficiency);
    // Singular matrix B: reorder the basic variables so that the
    // singular columns are in the position corresponding to the
    // logical which replaces them
    buildHandleRankDeficiency();
    buildMarkSingC();
    factor_timer.stop(FactorInvertDeficient, factor_timer_clock_pointer);
  }
  // Complete INVERT
  factor_timer.start(FactorInvertFinish, factor_timer_clock_pointer);
  buildFinish();
  factor_timer.stop(FactorInvertFinish, factor_timer_clock_pointer);
  // Record the number of entries in the INVERT
  invert_num_el = Lstart[numRow] + Ulastp[numRow - 1] + numRow;

  kernel_dim -= rank_deficiency;
  debugLogRankDeficiency(highs_debug_level, log_options, rank_deficiency,
                         basis_matrix_num_el, invert_num_el, kernel_dim,
                         kernel_num_el, nwork);
  factor_timer.stop(FactorInvert, factor_timer_clock_pointer);
  return rank_deficiency;
}

void HFactor::ftran(HVector& vector, double historical_density,
                    HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorFtran, factor_timer_clock_pointer);
  ftranL(vector, historical_density, factor_timer_clock_pointer);
  ftranU(vector, historical_density, factor_timer_clock_pointer);
  factor_timer.stop(FactorFtran, factor_timer_clock_pointer);
}

void HFactor::btran(HVector& vector, double historical_density,
                    HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorBtran, factor_timer_clock_pointer);
  btranU(vector, historical_density, factor_timer_clock_pointer);
  btranL(vector, historical_density, factor_timer_clock_pointer);
  factor_timer.stop(FactorBtran, factor_timer_clock_pointer);
}

void HFactor::update(HVector* aq, HVector* ep, HighsInt* iRow, HighsInt* hint) {
  // Special case
  if (aq->next) {
    updateCFT(aq, ep, iRow);
    return;
  }

  if (updateMethod == UPDATE_METHOD_FT) updateFT(aq, ep, *iRow);
  if (updateMethod == UPDATE_METHOD_PF) updatePF(aq, *iRow, hint);
  if (updateMethod == UPDATE_METHOD_MPF) updateMPF(aq, ep, *iRow, hint);
  if (updateMethod == UPDATE_METHOD_APF) updateAPF(aq, ep, *iRow);
}

bool HFactor::setPivotThreshold(const double new_pivot_threshold) {
  if (new_pivot_threshold < min_pivot_threshold) return false;
  if (new_pivot_threshold > max_pivot_threshold) return false;
  pivot_threshold = new_pivot_threshold;
  return true;
}

void HFactor::buildSimple() {
  /**
   * 0. Clear L and U factor
   */
  Lstart.clear();
  Lstart.push_back(0);
  Lindex.clear();
  Lvalue.clear();

  UpivotIndex.clear();
  UpivotValue.clear();
  Ustart.clear();
  Ustart.push_back(0);
  Uindex.clear();
  Uvalue.clear();

  // Set all values of permute to -1 so that unpermuted (rank
  // deficient) columns canm be identified
  permute.assign(numRow, -1);

  /**
   * 1. Prepare basis matrix and deal with unit columns
   */

  HighsInt BcountX = 0;
  fill_n(&MRcountb4[0], numRow, 0);
  nwork = 0;
  for (HighsInt iCol = 0; iCol < numRow; iCol++) {
    HighsInt iMat = baseIndex[iCol];
    HighsInt iRow = -1;
    if (iMat >= numCol) {
      // 1.1 Logical column
      // Check for double pivot
      HighsInt lc_iRow = iMat - numCol;
      if (MRcountb4[lc_iRow] >= 0) {
        iRow = lc_iRow;
      } else {
        highsLogUser(log_options, HighsLogType::ERROR,
                     "INVERT Error: Found a logical column with pivot "
                     "already in row %d\n",
                     lc_iRow);
        MRcountb4[lc_iRow]++;
        Bindex[BcountX] = lc_iRow;
        Bvalue[BcountX++] = 1.0;
        iwork[nwork++] = iCol;
      }
    } else {
      // 1.2 Structural column
      HighsInt start = Astart[iMat];
      HighsInt count = Astart[iMat + 1] - start;
      HighsInt lc_iRow = Aindex[start];
      // Check for unit column with double pivot
      bool unit_col = count == 1 && Avalue[start] == 1;
      if (unit_col && MRcountb4[lc_iRow] >= 0) {
        iRow = lc_iRow;
      } else {
        if (unit_col)
          highsLogUser(
              log_options, HighsLogType::ERROR,
              "INVERT Error: Found a second unit column with pivot in row %d\n",
              lc_iRow);
        for (HighsInt k = start; k < start + count; k++) {
          MRcountb4[Aindex[k]]++;
          Bindex[BcountX] = Aindex[k];
          Bvalue[BcountX++] = Avalue[k];
        }
        iwork[nwork++] = iCol;
      }
    }

    if (iRow >= 0) {
      // 1.3 Record unit column
      // Uindex.size());
      permute[iCol] = iRow;
      Lstart.push_back(Lindex.size());
      UpivotIndex.push_back(iRow);
      UpivotValue.push_back(1);
      Ustart.push_back(Uindex.size());
      MRcountb4[iRow] = -numRow;
    }
    Bstart[iCol + 1] = BcountX;
  }
  // Record the number of elements in the basis matrix
  basis_matrix_num_el = numRow - nwork + BcountX;

  // count1 = 0;
  // Comments: for pds-20, dfl001: 60 / 80
  // Comments: when system is large: enlarge
  // Comments: when system is small: decrease
  build_syntheticTick += BcountX * 60 + (numRow - nwork) * 80;

  /**
   * 2. Search for and deal with singletons
   */
  double t2_search = 0;
  double t2_storeL = Lindex.size();
  double t2_storeU = Uindex.size();
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
        for (HighsInt section = 0; section < 2; section++) {
          HighsInt p0 = section == 0 ? start : pivot_k + 1;
          HighsInt p1 = section == 0 ? pivot_k : end;
          for (HighsInt k = p0; k < p1; k++) {
            HighsInt iRow = Bindex[k];
            if (MRcountb4[iRow] > 0) {
              Lindex.push_back(iRow);
              Lvalue.push_back(Bvalue[k] * pivotX);
            } else {
              Uindex.push_back(iRow);
              Uvalue.push_back(Bvalue[k]);
            }
            MRcountb4[iRow]--;
          }
        }
        HighsInt iRow = Bindex[pivot_k];
        MRcountb4[iRow] = 0;
        permute[iCol] = iRow;
        Lstart.push_back(Lindex.size());

        UpivotIndex.push_back(iRow);
        UpivotValue.push_back(Bvalue[pivot_k]);
        Ustart.push_back(Uindex.size());
      } else if (count == 1) {
        // 2.3 Deal with column singleton
        for (HighsInt k = start; k < pivot_k; k++) {
          Uindex.push_back(Bindex[k]);
          Uvalue.push_back(Bvalue[k]);
        }
        for (HighsInt k = pivot_k + 1; k < end; k++) {
          Uindex.push_back(Bindex[k]);
          Uvalue.push_back(Bvalue[k]);
        }

        HighsInt iRow = Bindex[pivot_k];
        MRcountb4[iRow] = 0;
        permute[iCol] = iRow;
        Lstart.push_back(Lindex.size());

        UpivotIndex.push_back(iRow);
        UpivotValue.push_back(Bvalue[pivot_k]);
        Ustart.push_back(Uindex.size());
      } else {
        iwork[nwork++] = iCol;
      }
    }

    // No singleton found in the last pass
    if (nworkLast == nwork) break;
  }
  t2_storeL = Lindex.size() - t2_storeL;
  t2_storeU = Uindex.size() - t2_storeU;
  t2_storep = t2_storep - nwork;

  build_syntheticTick +=
      t2_search * 20 + (t2_storep + t2_storeL + t2_storeU) * 80;

  /**
   * 3. Prepare the kernel parts
   */
  // 3.1 Prepare row links, row matrix spaces
  rlinkFirst.assign(numRow + 1, -1);
  MRcount.assign(numRow, 0);
  HighsInt MRcountX = 0;
  // Determine the number of entries in the kernel
  kernel_num_el = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
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
  clinkFirst.assign(numRow + 1, -1);
  MCindex.clear();
  MCvalue.clear();
  MCcountA.assign(numRow, 0);
  MCcountN.assign(numRow, 0);
  HighsInt MCcountX = 0;
  for (HighsInt i = 0; i < nwork; i++) {
    HighsInt iCol = iwork[i];
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
  build_syntheticTick += (numRow + nwork + MCcountX) * 40 + MRcountX * 20;
  // Record the kernel dimension
  kernel_dim = nwork;
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
    /*
    bool rp_r_k = false;
    if (rp_r_k) {
      printf("Row counts:");
      bool f_k = true;
      for (HighsInt k = 0; k < numRow; k++) {
        if (rlinkFirst[k] >= 0) {
          if (f_k) {
            printf(" (%2d:", k);
            f_k = false;
          } else {
            printf("; (%2d:", k);
          }
          for (HighsInt i = rlinkFirst[k]; i != -1; i = rlinkNext[i]) {
            printf(" %2d", i);
          }
          printf(")");
        }
      }
      printf("\n");
    }
    bool rp_permute = false;
    if (rp_permute) {
      printf("Permute:\n");
      for (HighsInt i = 0; i < numRow; i++) {
        printf(" %2d", i);
      }
      printf("\n");
      for (HighsInt i = 0; i < numRow; i++) {
        printf(" %2d", permute[i]);
      }
      printf("\n");
    }
    */
    HighsInt jColPivot = -1;
    HighsInt iRowPivot = -1;
    // 1.1. Setup search merits
    HighsInt searchLimit = min(nwork, HighsInt{8});
    HighsInt searchCount = 0;
    double meritLimit = 1.0 * numRow * numRow;
    double meritPivot = meritLimit;

    // 1.2. Search for local singletons
    bool foundPivot = false;
    if (!foundPivot && clinkFirst[1] != -1) {
      jColPivot = clinkFirst[1];
      iRowPivot = MCindex[MCstart[jColPivot]];
      foundPivot = true;
    }
    if (!foundPivot && rlinkFirst[1] != -1) {
      iRowPivot = rlinkFirst[1];
      jColPivot = MRindex[MRstart[iRowPivot]];
      foundPivot = true;
    }
    const bool singleton_pivot = foundPivot;
    // 1.3. Major search loop
    double candidate_pivot_value = 0;
    for (HighsInt count = 2; !foundPivot && count <= numRow; count++) {
      // 1.3.1 Search for columns
      for (HighsInt j = clinkFirst[count]; j != -1; j = clinkNext[j]) {
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
      for (HighsInt i = rlinkFirst[count]; i != -1; i = rlinkNext[i]) {
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
      highsLogUser(log_options, HighsLogType::WARNING,
                   "Small |pivot| = %g when nwork = %d\n", fabs(pivotX), nwork);
      rank_deficiency = nwork + 1;
      return rank_deficiency;
    }
    rowDelete(jColPivot, iRowPivot);
    clinkDel(jColPivot);
    rlinkDel(iRowPivot);
    permute[jColPivot] = iRowPivot;

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
      Lindex.push_back(iRow);
      Lvalue.push_back(value);
      MRcountb4[iRow] = MRcount[iRow];
      rowDelete(jColPivot, iRow);
    }
    Lstart.push_back(Lindex.size());
    fake_fill += 2 * MCcountA[jColPivot];

    // 2.3. Store non active pivot column to U
    HighsInt end_N = start_A + MCspace[jColPivot];
    HighsInt start_N = end_N - MCcountN[jColPivot];
    for (HighsInt i = start_N; i < end_N; i++) {
      Uindex.push_back(MCindex[i]);
      Uvalue.push_back(MCvalue[i]);
    }
    UpivotIndex.push_back(iRowPivot);
    UpivotValue.push_back(pivotX);
    Ustart.push_back(Uindex.size());
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
          if (fabs(value) < HIGHS_CONST_TINY) {
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
  build_syntheticTick +=
      fake_search * 20 + fake_fill * 160 + fake_eliminate * 80;
  rank_deficiency = 0;
  return rank_deficiency;
}

void HFactor::buildHandleRankDeficiency() {
  debugReportRankDeficiency(0, highs_debug_level, log_options, numRow, permute,
                            iwork, baseIndex, rank_deficiency, noPvR, noPvC);
  // iwork can now be used as workspace: use it to accumulate the new
  // baseIndex. iwork is set to -1 and baseIndex is permuted into it.
  // Indices of iwork corresponding to missing indices in permute
  // remain -1. Hence the -1's become markers for the logicals which
  // will replace singular columns. Once baseIndex[i] is read, it can
  // be used to pack up the entries in baseIndex which are not
  // permuted anywhere - and so will be singular columns.
  noPvR.resize(rank_deficiency);
  noPvC.resize(rank_deficiency);
  HighsInt lc_rank_deficiency = 0;
  for (HighsInt i = 0; i < numRow; i++) iwork[i] = -1;
  for (HighsInt i = 0; i < numRow; i++) {
    HighsInt perm_i = permute[i];
    if (perm_i >= 0) {
      iwork[perm_i] = baseIndex[i];
    } else {
      noPvC[lc_rank_deficiency++] = i;
    }
  }
  assert(lc_rank_deficiency == rank_deficiency);
  lc_rank_deficiency = 0;
  for (HighsInt i = 0; i < numRow; i++) {
    if (iwork[i] < 0) {
      // Record the rows with no pivots in noPvR and indicate them
      // within iwork by storing the negation of one more than their
      // rank deficiency counter [since we can't have -0].
      noPvR[lc_rank_deficiency] = i;
      iwork[i] = -(lc_rank_deficiency + 1);
      lc_rank_deficiency++;
    }
  }
  assert(lc_rank_deficiency == rank_deficiency);
  debugReportRankDeficiency(1, highs_debug_level, log_options, numRow, permute,
                            iwork, baseIndex, rank_deficiency, noPvR, noPvC);
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    HighsInt iRow = noPvR[k];
    HighsInt iCol = noPvC[k];
    assert(permute[iCol] == -1);
    permute[iCol] = iRow;
    Lstart.push_back(Lindex.size());
    UpivotIndex.push_back(iRow);
    UpivotValue.push_back(1);
    Ustart.push_back(Uindex.size());
  }
  debugReportRankDeficiency(2, highs_debug_level, log_options, numRow, permute,
                            iwork, baseIndex, rank_deficiency, noPvR, noPvC);
  debugReportRankDeficientASM(highs_debug_level, log_options, numRow, MCstart,
                              MCcountA, MCindex, MCvalue, iwork,
                              rank_deficiency, noPvC, noPvR);
}

void HFactor::buildMarkSingC() {
  // Singular matrix B: reorder the basic variables so that the
  // singular columns are in the position corresponding to the
  // logical which replaces them
  debugReportMarkSingC(0, highs_debug_level, log_options, numRow, iwork,
                       baseIndex);

  for (HighsInt k = 0; k < rank_deficiency; k++) {
    HighsInt ASMrow = noPvR[k];
    HighsInt ASMcol = noPvC[k];
    assert(-iwork[ASMrow] - 1 >= 0 && -iwork[ASMrow] - 1 < rank_deficiency);
    // Store negation of 1+ASMcol so that removing column 0 can be
    // identified!
    iwork[ASMrow] = -(ASMcol + 1);
    noPvC[k] = baseIndex[ASMcol];
    baseIndex[ASMcol] = numCol + ASMrow;
  }
  debugReportMarkSingC(1, highs_debug_level, log_options, numRow, iwork,
                       baseIndex);
}

void HFactor::buildFinish() {
  //  debugPivotValueAnalysis(highs_debug_level, log_options, numRow,
  //  UpivotValue);
  // The look up table
  for (HighsInt i = 0; i < numRow; i++) UpivotLookup[UpivotIndex[i]] = i;
  LpivotIndex = UpivotIndex;
  LpivotLookup = UpivotLookup;

  // LR space
  HighsInt LcountX = Lindex.size();
  LRindex.resize(LcountX);
  LRvalue.resize(LcountX);

  // LR pointer
  iwork.assign(numRow, 0);
  for (HighsInt k = 0; k < LcountX; k++) iwork[LpivotLookup[Lindex[k]]]++;

  LRstart.assign(numRow + 1, 0);
  for (HighsInt i = 1; i <= numRow; i++)
    LRstart[i] = LRstart[i - 1] + iwork[i - 1];

  // LR elements
  iwork.assign(&LRstart[0], &LRstart[numRow]);
  for (HighsInt i = 0; i < numRow; i++) {
    const HighsInt index = LpivotIndex[i];
    for (HighsInt k = Lstart[i]; k < Lstart[i + 1]; k++) {
      HighsInt iRow = LpivotLookup[Lindex[k]];
      HighsInt iPut = iwork[iRow]++;
      LRindex[iPut] = index;
      LRvalue[iPut] = Lvalue[k];
    }
  }

  // U pointer
  Ustart.push_back(0);
  Ulastp.assign(&Ustart[1], &Ustart[numRow + 1]);
  Ustart.resize(numRow);

  // UR space
  HighsInt UcountX = Uindex.size();
  HighsInt URstuffX = updateMethod == UPDATE_METHOD_FT ? 5 : 0;
  HighsInt URcountX = UcountX + URstuffX * numRow;
  URindex.resize(URcountX);
  URvalue.resize(URcountX);

  // UR pointer
  URstart.assign(numRow + 1, 0);
  URlastp.assign(numRow, 0);
  URspace.assign(numRow, URstuffX);
  for (HighsInt k = 0; k < UcountX; k++) URlastp[UpivotLookup[Uindex[k]]]++;
  for (HighsInt i = 1; i <= numRow; i++)
    URstart[i] = URstart[i - 1] + URlastp[i - 1] + URstuffX;
  URstart.resize(numRow);

  // UR element
  URlastp = URstart;
  for (HighsInt i = 0; i < numRow; i++) {
    const HighsInt index = UpivotIndex[i];
    for (HighsInt k = Ustart[i]; k < Ulastp[i]; k++) {
      HighsInt iRow = UpivotLookup[Uindex[k]];
      HighsInt iPut = URlastp[iRow]++;
      URindex[iPut] = index;
      URvalue[iPut] = Uvalue[k];
    }
  }

  // Re-factor merit
  UmeritX = numRow + (LcountX + UcountX) * 1.5;
  UtotalX = UcountX;
  if (updateMethod == UPDATE_METHOD_PF) UmeritX = numRow + UcountX * 4;
  if (updateMethod == UPDATE_METHOD_MPF) UmeritX = numRow + UcountX * 3;

  // Clear update buffer
  PFpivotValue.clear();
  PFpivotIndex.clear();
  PFstart.clear();
  PFstart.push_back(0);
  PFindex.clear();
  PFvalue.clear();

  // Finally, permute the base index
  iwork.assign(baseIndex, baseIndex + numRow);
  for (HighsInt i = 0; i < numRow; i++) baseIndex[permute[i]] = iwork[i];

  build_syntheticTick += numRow * 80 + (LcountX + UcountX) * 60;
}

void HFactor::ftranL(HVector& rhs, double historical_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorFtranLower, factor_timer_clock_pointer);
  if (updateMethod == UPDATE_METHOD_APF) {
    factor_timer.start(FactorFtranLowerAPF, factor_timer_clock_pointer);
    rhs.tight();
    rhs.pack();
    ftranAPF(rhs);
    factor_timer.stop(FactorFtranLowerAPF, factor_timer_clock_pointer);
    rhs.tight();
  }

  double current_density = 1.0 * rhs.count / numRow;
  if (current_density > hyperCANCEL || historical_density > hyperFTRANL) {
    factor_timer.start(FactorFtranLowerSps, factor_timer_clock_pointer);
    // Alias to RHS
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to factor L
    const HighsInt* Lstart = &this->Lstart[0];
    const HighsInt* Lindex = this->Lindex.size() > 0 ? &this->Lindex[0] : NULL;
    const double* Lvalue = this->Lvalue.size() > 0 ? &this->Lvalue[0] : NULL;

    // Transform
    for (HighsInt i = 0; i < numRow; i++) {
      HighsInt pivotRow = LpivotIndex[i];
      const double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > HIGHS_CONST_TINY) {
        RHSindex[RHScount++] = pivotRow;
        const HighsInt start = Lstart[i];
        const HighsInt end = Lstart[i + 1];
        for (HighsInt k = start; k < end; k++)
          RHSarray[Lindex[k]] -= pivotX * Lvalue[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    factor_timer.stop(FactorFtranLowerSps, factor_timer_clock_pointer);
  } else {
    factor_timer.start(FactorFtranLowerHyper, factor_timer_clock_pointer);
    const HighsInt* Lindex = this->Lindex.size() > 0 ? &this->Lindex[0] : NULL;
    const double* Lvalue = this->Lvalue.size() > 0 ? &this->Lvalue[0] : NULL;
    solveHyper(numRow, &LpivotLookup[0], &LpivotIndex[0], 0, &Lstart[0],
               &Lstart[1], &Lindex[0], &Lvalue[0], &rhs);
    factor_timer.stop(FactorFtranLowerHyper, factor_timer_clock_pointer);
  }
  factor_timer.stop(FactorFtranLower, factor_timer_clock_pointer);
}

void HFactor::btranL(HVector& rhs, double historical_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorBtranLower, factor_timer_clock_pointer);
  double current_density = 1.0 * rhs.count / numRow;
  if (current_density > hyperCANCEL || historical_density > hyperBTRANL) {
    // Alias to RHS
    factor_timer.start(FactorBtranLowerSps, factor_timer_clock_pointer);
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to factor L
    const HighsInt* LRstart = &this->LRstart[0];
    const HighsInt* LRindex =
        this->LRindex.size() > 0 ? &this->LRindex[0] : NULL;
    const double* LRvalue = this->LRvalue.size() > 0 ? &this->LRvalue[0] : NULL;

    // Transform
    for (HighsInt i = numRow - 1; i >= 0; i--) {
      HighsInt pivotRow = LpivotIndex[i];
      const double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > HIGHS_CONST_TINY) {
        RHSindex[RHScount++] = pivotRow;
        RHSarray[pivotRow] = pivotX;
        const HighsInt start = LRstart[i];
        const HighsInt end = LRstart[i + 1];
        for (HighsInt k = start; k < end; k++)
          RHSarray[LRindex[k]] -= pivotX * LRvalue[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    factor_timer.stop(FactorBtranLowerSps, factor_timer_clock_pointer);
  } else {
    factor_timer.start(FactorBtranLowerHyper, factor_timer_clock_pointer);
    const HighsInt* LRindex =
        this->LRindex.size() > 0 ? &this->LRindex[0] : NULL;
    const double* LRvalue = this->LRvalue.size() > 0 ? &this->LRvalue[0] : NULL;
    solveHyper(numRow, &LpivotLookup[0], &LpivotIndex[0], 0, &LRstart[0],
               &LRstart[1], &LRindex[0], &LRvalue[0], &rhs);
    factor_timer.stop(FactorBtranLowerHyper, factor_timer_clock_pointer);
  }

  if (updateMethod == UPDATE_METHOD_APF) {
    factor_timer.start(FactorBtranLowerAPF, factor_timer_clock_pointer);
    btranAPF(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorBtranLowerAPF, factor_timer_clock_pointer);
  }
  factor_timer.stop(FactorBtranLower, factor_timer_clock_pointer);
}

void HFactor::ftranU(HVector& rhs, double historical_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorFtranUpper, factor_timer_clock_pointer);
  // The update part
  if (updateMethod == UPDATE_METHOD_FT) {
    factor_timer.start(FactorFtranUpperFT, factor_timer_clock_pointer);
    //    const double current_density = 1.0 * rhs.count / numRow;
    ftranFT(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorFtranUpperFT, factor_timer_clock_pointer);
  }
  if (updateMethod == UPDATE_METHOD_MPF) {
    factor_timer.start(FactorFtranUpperMPF, factor_timer_clock_pointer);
    ftranMPF(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorFtranUpperMPF, factor_timer_clock_pointer);
  }

  // The regular part
  const double current_density = 1.0 * rhs.count / numRow;
  if (current_density > hyperCANCEL || historical_density > hyperFTRANU) {
    const bool report_ftran_upper_sparse =
        false;  // current_density < hyperCANCEL;
    HighsInt use_clock;
    if (current_density < 0.1)
      use_clock = FactorFtranUpperSps2;
    else if (current_density < 0.5)
      use_clock = FactorFtranUpperSps1;
    else
      use_clock = FactorFtranUpperSps0;
    factor_timer.start(use_clock, factor_timer_clock_pointer);
    // Alias to non constant
    double RHS_syntheticTick = 0;
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to the factor
    const HighsInt* Ustart = &this->Ustart[0];
    const HighsInt* Uend = &this->Ulastp[0];
    const HighsInt* Uindex = this->Uindex.size() > 0 ? &this->Uindex[0] : NULL;
    const double* Uvalue = this->Uvalue.size() > 0 ? &this->Uvalue[0] : NULL;

    // Transform
    HighsInt UpivotCount = UpivotIndex.size();
    for (HighsInt iLogic = UpivotCount - 1; iLogic >= 0; iLogic--) {
      // Skip void
      if (UpivotIndex[iLogic] == -1) continue;

      // Normal part
      const HighsInt pivotRow = UpivotIndex[iLogic];
      double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > HIGHS_CONST_TINY) {
        pivotX /= UpivotValue[iLogic];
        RHSindex[RHScount++] = pivotRow;
        RHSarray[pivotRow] = pivotX;
        const HighsInt start = Ustart[iLogic];
        const HighsInt end = Uend[iLogic];
        if (iLogic >= numRow) {
          RHS_syntheticTick += (end - start);
        }
        for (HighsInt k = start; k < end; k++)
          RHSarray[Uindex[k]] -= pivotX * Uvalue[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    rhs.syntheticTick += RHS_syntheticTick * 15 + (UpivotCount - numRow) * 10;
    factor_timer.stop(use_clock, factor_timer_clock_pointer);
    if (report_ftran_upper_sparse) {
      const double final_density = 1.0 * rhs.count / numRow;
      printf(
          "FactorFtranUpperSps: historical_density = %10.4g; current_density = "
          "%10.4g; final_density = %10.4g\n",
          historical_density, current_density, final_density);
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
    const HighsInt* Uindex = this->Uindex.size() > 0 ? &this->Uindex[0] : NULL;
    const double* Uvalue = this->Uvalue.size() > 0 ? &this->Uvalue[0] : NULL;
    solveHyper(numRow, &UpivotLookup[0], &UpivotIndex[0], &UpivotValue[0],
               &Ustart[0], &Ulastp[0], &Uindex[0], &Uvalue[0], &rhs);
    factor_timer.stop(use_clock, factor_timer_clock_pointer);
  }
  if (updateMethod == UPDATE_METHOD_PF) {
    factor_timer.start(FactorFtranUpperPF, factor_timer_clock_pointer);
    ftranPF(rhs);
    rhs.tight();
    rhs.pack();
    factor_timer.stop(FactorFtranUpperPF, factor_timer_clock_pointer);
  }
  factor_timer.stop(FactorFtranUpper, factor_timer_clock_pointer);
}

void HFactor::btranU(HVector& rhs, double historical_density,
                     HighsTimerClock* factor_timer_clock_pointer) const {
  FactorTimer factor_timer;
  factor_timer.start(FactorBtranUpper, factor_timer_clock_pointer);
  if (updateMethod == UPDATE_METHOD_PF) {
    factor_timer.start(FactorBtranUpperPF, factor_timer_clock_pointer);
    btranPF(rhs);
    factor_timer.stop(FactorBtranUpperPF, factor_timer_clock_pointer);
  }

  // The regular part
  double current_density = 1.0 * rhs.count / numRow;
  if (current_density > hyperCANCEL || historical_density > hyperBTRANU) {
    factor_timer.start(FactorBtranUpperSps, factor_timer_clock_pointer);
    // Alias to non constant
    double RHS_syntheticTick = 0;
    HighsInt RHScount = 0;
    HighsInt* RHSindex = &rhs.index[0];
    double* RHSarray = &rhs.array[0];

    // Alias to the factor
    const HighsInt* URstart = &this->URstart[0];
    const HighsInt* URend = &this->URlastp[0];
    const HighsInt* URindex = &this->URindex[0];
    const double* URvalue = &this->URvalue[0];

    // Transform
    HighsInt UpivotCount = UpivotIndex.size();
    for (HighsInt iLogic = 0; iLogic < UpivotCount; iLogic++) {
      // Skip void
      if (UpivotIndex[iLogic] == -1) continue;

      // Normal part
      const HighsInt pivotRow = UpivotIndex[iLogic];
      double pivotX = RHSarray[pivotRow];
      if (fabs(pivotX) > HIGHS_CONST_TINY) {
        pivotX /= UpivotValue[iLogic];
        RHSindex[RHScount++] = pivotRow;
        RHSarray[pivotRow] = pivotX;
        const HighsInt start = URstart[iLogic];
        const HighsInt end = URend[iLogic];
        if (iLogic >= numRow) {
          RHS_syntheticTick += (end - start);
        }
        for (HighsInt k = start; k < end; k++)
          RHSarray[URindex[k]] -= pivotX * URvalue[k];
      } else
        RHSarray[pivotRow] = 0;
    }

    // Save the count
    rhs.count = RHScount;
    rhs.syntheticTick += RHS_syntheticTick * 15 + (UpivotCount - numRow) * 10;
    factor_timer.stop(FactorBtranUpperSps, factor_timer_clock_pointer);
  } else {
    factor_timer.start(FactorBtranUpperHyper, factor_timer_clock_pointer);
    solveHyper(numRow, &UpivotLookup[0], &UpivotIndex[0], &UpivotValue[0],
               &URstart[0], &URlastp[0], &URindex[0], &URvalue[0], &rhs);
    factor_timer.stop(FactorBtranUpperHyper, factor_timer_clock_pointer);
  }

  // The update part
  if (updateMethod == UPDATE_METHOD_FT) {
    factor_timer.start(FactorBtranUpperFT, factor_timer_clock_pointer);
    rhs.tight();
    rhs.pack();
    //    const double current_density = 1.0 * rhs.count / numRow;
    btranFT(rhs);
    rhs.tight();
    factor_timer.stop(FactorBtranUpperFT, factor_timer_clock_pointer);
  }
  if (updateMethod == UPDATE_METHOD_MPF) {
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
      RHSarray[iRow] =
          (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
    }
  }

  // Save count back
  vector.count = RHScount;
  vector.syntheticTick += PFpivotCount * 20 + PFstart[PFpivotCount] * 5;
  if (PFstart[PFpivotCount] / (PFpivotCount + 1) < 5) {
    vector.syntheticTick += PFstart[PFpivotCount] * 5;
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
  double RHS_syntheticTick = 0;
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
      RHS_syntheticTick += (end - start);
      for (HighsInt k = start; k < end; k++) {
        HighsInt iRow = PFindex[k];
        double value0 = RHSarray[iRow];
        double value1 = value0 - pivotX * PFvalue[k];
        if (value0 == 0) RHSindex[RHScount++] = iRow;
        RHSarray[iRow] =
            (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
      }
    }
  }

  vector.syntheticTick += RHS_syntheticTick * 15 + PFpivotCount * 10;

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
    if (fabs(pivotX) > HIGHS_CONST_TINY) {
      pivotX /= PFpivotValue[i];
      RHSarray[pivotRow] = pivotX;
      for (HighsInt k = PFstart[i]; k < PFstart[i + 1]; k++) {
        const HighsInt index = PFindex[k];
        const double value0 = RHSarray[index];
        const double value1 = value0 - pivotX * PFvalue[k];
        if (value0 == 0) RHSindex[RHScount++] = index;
        RHSarray[index] =
            (fabs(value1) < HIGHS_CONST_TINY) ? HIGHS_CONST_ZERO : value1;
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
    RHSarray[pivotRow] = (fabs(pivotX) < HIGHS_CONST_TINY) ? 1e-100 : pivotX;
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
    HighsInt iLogic = UpivotLookup[cRow];
    pLogic[cp] = iLogic;
    pValue[cp] = UpivotValue[iLogic];
    pAlpha[cp] = aqWork[cp]->array[cRow];
  }

  // Temporary U pointers
  HighsInt* Tstart = new HighsInt[numUpdate + 1];
  double* Tpivot = new double[numUpdate];
  Tstart[0] = Uindex.size();

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
    HighsInt UstartX = UcountX;
    for (unsigned i = 0; i < iwork.size(); i++) {
      HighsInt index = iwork[i];
      double value = dwork[index];
      dwork[index] = 0;  // This effectively removes all duplication
      if (fabs(value) > HIGHS_CONST_TINY) {
        Uindex.push_back(index);
        Uvalue.push_back(value);
      }
    }
    UcountX = Uindex.size();
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
      if (fabs(dwork[pRow]) > HIGHS_CONST_TINY) {
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
        value -= dwork[Uindex[k]] * Uvalue[k];
      value /= Tpivot[pp];
      iwork.push_back(kpivot);
      dwork[kpivot] = value;  // Again OK to duplicate
    }

    // 6.x compute current alpha
    double thex = 0;
    for (HighsInt k = UstartX; k < UcountX; k++) {
      HighsInt index = Uindex[k];
      double value = Uvalue[k];
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
      if (fabs(value) > HIGHS_CONST_TINY) {
        PFindex.push_back(index);
        PFvalue.push_back(value * pivotX);
      }
    }
    PFpivotIndex.push_back(iRow[cp]);
    UtotalX += PFindex.size() - PFstart.back();
    PFstart.push_back(PFindex.size());

    // 8. Update the sorted ep
    sorted_pp.push_back(make_pair(pLogic[cp], cp));
    sort(sorted_pp.begin(), sorted_pp.end());
  }

  // Now modify the U matrix
  for (HighsInt cp = 0; cp < numUpdate; cp++) {
    // 1. Delete pivotal row from U
    HighsInt cIndex = iRow[cp];
    HighsInt cLogic = pLogic[cp];
    UtotalX -= URlastp[cLogic] - URstart[cLogic];
    for (HighsInt k = URstart[cLogic]; k < URlastp[cLogic]; k++) {
      // Find the pivotal position
      HighsInt iLogic = UpivotLookup[URindex[k]];
      HighsInt iFind = Ustart[iLogic];
      HighsInt iLast = --Ulastp[iLogic];
      for (; iFind <= iLast; iFind++)
        if (Uindex[iFind] == cIndex) break;
      // Put last to find, and delete last
      Uindex[iFind] = Uindex[iLast];
      Uvalue[iFind] = Uvalue[iLast];
    }

    // 2. Delete pivotal column from UR
    UtotalX -= Ulastp[cLogic] - Ustart[cLogic];
    for (HighsInt k = Ustart[cLogic]; k < Ulastp[cLogic]; k++) {
      // Find the pivotal position
      HighsInt iLogic = UpivotLookup[Uindex[k]];
      HighsInt iFind = URstart[iLogic];
      HighsInt iLast = --URlastp[iLogic];
      for (; iFind <= iLast; iFind++)
        if (URindex[iFind] == cIndex) break;
      // Put last to find, and delete last
      URspace[iLogic]++;
      URindex[iFind] = URindex[iLast];
      URvalue[iFind] = URvalue[iLast];
    }

    // 3. Insert the (stored) partial FTRAN to the row matrix
    HighsInt UstartX = Tstart[cp];
    HighsInt UendX = Tstart[cp + 1];
    UtotalX += UendX - UstartX;
    // Store column as UR elements
    for (HighsInt k = UstartX; k < UendX; k++) {
      // Which ETA file
      HighsInt iLogic = UpivotLookup[Uindex[k]];

      // Move row to the end if necessary
      if (URspace[iLogic] == 0) {
        // Make pointers
        HighsInt row_start = URstart[iLogic];
        HighsInt row_count = URlastp[iLogic] - row_start;
        HighsInt new_start = URindex.size();
        HighsInt new_space = row_count * 1.1 + 5;

        // Check matrix UR
        URindex.resize(new_start + new_space);
        URvalue.resize(new_start + new_space);

        // Move elements
        HighsInt iFrom = row_start;
        HighsInt iEnd = row_start + row_count;
        HighsInt iTo = new_start;
        copy(&URindex[iFrom], &URindex[iEnd], &URindex[iTo]);
        copy(&URvalue[iFrom], &URvalue[iEnd], &URvalue[iTo]);

        // Save new pointers
        URstart[iLogic] = new_start;
        URlastp[iLogic] = new_start + row_count;
        URspace[iLogic] = new_space - row_count;
      }

      // Put into the next available space
      URspace[iLogic]--;
      HighsInt iPut = URlastp[iLogic]++;
      URindex[iPut] = cIndex;
      URvalue[iPut] = Uvalue[k];
    }

    // 4. Save pointers
    Ustart.push_back(UstartX);
    Ulastp.push_back(UendX);

    URstart.push_back(URstart[cLogic]);
    URlastp.push_back(URstart[cLogic]);
    URspace.push_back(URspace[cLogic] + URlastp[cLogic] - URstart[cLogic]);

    UpivotLookup[cIndex] = UpivotIndex.size();
    UpivotIndex[cLogic] = -1;
    UpivotIndex.push_back(cIndex);
    UpivotValue.push_back(Tpivot[cp]);
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
  HighsInt pLogic = UpivotLookup[iRow];
  double pivot = UpivotValue[pLogic];
  double alpha = aq->array[iRow];
  UpivotIndex[pLogic] = -1;

  // Delete pivotal row from U
  for (HighsInt k = URstart[pLogic]; k < URlastp[pLogic]; k++) {
    // Find the pivotal position
    HighsInt iLogic = UpivotLookup[URindex[k]];
    HighsInt iFind = Ustart[iLogic];
    HighsInt iLast = --Ulastp[iLogic];
    for (; iFind <= iLast; iFind++)
      if (Uindex[iFind] == iRow) break;
    // Put last to find, and delete last
    Uindex[iFind] = Uindex[iLast];
    Uvalue[iFind] = Uvalue[iLast];
  }

  // Delete pivotal column from UR
  for (HighsInt k = Ustart[pLogic]; k < Ulastp[pLogic]; k++) {
    // Find the pivotal position
    HighsInt iLogic = UpivotLookup[Uindex[k]];
    HighsInt iFind = URstart[iLogic];
    HighsInt iLast = --URlastp[iLogic];
    for (; iFind <= iLast; iFind++)
      if (URindex[iFind] == iRow) break;
    // Put last to find, and delete last
    URspace[iLogic]++;
    URindex[iFind] = URindex[iLast];
    URvalue[iFind] = URvalue[iLast];
  }

  // Store column to U
  Ustart.push_back(Uindex.size());
  for (HighsInt i = 0; i < aq->packCount; i++)
    if (aq->packIndex[i] != iRow) {
      Uindex.push_back(aq->packIndex[i]);
      Uvalue.push_back(aq->packValue[i]);
    }
  Ulastp.push_back(Uindex.size());
  HighsInt UstartX = Ustart.back();
  HighsInt UendX = Ulastp.back();
  UtotalX += UendX - UstartX + 1;

  // Store column as UR elements
  for (HighsInt k = UstartX; k < UendX; k++) {
    // Which ETA file
    HighsInt iLogic = UpivotLookup[Uindex[k]];

    // Move row to the end if necessary
    if (URspace[iLogic] == 0) {
      // Make pointers
      HighsInt row_start = URstart[iLogic];
      HighsInt row_count = URlastp[iLogic] - row_start;
      HighsInt new_start = URindex.size();
      HighsInt new_space = row_count * 1.1 + 5;

      // Check matrix UR
      URindex.resize(new_start + new_space);
      URvalue.resize(new_start + new_space);

      // Move elements
      HighsInt iFrom = row_start;
      HighsInt iEnd = row_start + row_count;
      HighsInt iTo = new_start;
      copy(&URindex[iFrom], &URindex[iEnd], &URindex[iTo]);
      copy(&URvalue[iFrom], &URvalue[iEnd], &URvalue[iTo]);

      // Save new pointers
      URstart[iLogic] = new_start;
      URlastp[iLogic] = new_start + row_count;
      URspace[iLogic] = new_space - row_count;
    }

    // Put into the next available space
    URspace[iLogic]--;
    HighsInt iPut = URlastp[iLogic]++;
    URindex[iPut] = iRow;
    URvalue[iPut] = Uvalue[k];
  }

  // Store UR pointers
  URstart.push_back(URstart[pLogic]);
  URlastp.push_back(URstart[pLogic]);
  URspace.push_back(URspace[pLogic] + URlastp[pLogic] - URstart[pLogic]);

  // Update pivot count
  UpivotLookup[iRow] = UpivotIndex.size();
  UpivotIndex.push_back(iRow);
  UpivotValue.push_back(pivot * alpha);

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
  UtotalX -= Ulastp[pLogic] - Ustart[pLogic];
  UtotalX -= URlastp[pLogic] - URstart[pLogic];

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
  HighsInt pLogic = UpivotLookup[iRow];
  HighsInt UstartX = Ustart[pLogic];
  HighsInt UendX = Ustart[pLogic + 1];
  for (HighsInt k = UstartX; k < UendX; k++) {
    PFindex.push_back(Uindex[k]);
    PFvalue.push_back(-Uvalue[k]);
  }
  PFindex.push_back(iRow);
  PFvalue.push_back(-UpivotValue[pLogic]);
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
  if (variable_out >= numCol) {
    PFindex.push_back(variable_out - numCol);
    PFvalue.push_back(-1);
  } else {
    for (HighsInt k = Astart[variable_out]; k < Astart[variable_out + 1]; k++) {
      PFindex.push_back(Aindex[k]);
      PFvalue.push_back(-Avalue[k]);
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
