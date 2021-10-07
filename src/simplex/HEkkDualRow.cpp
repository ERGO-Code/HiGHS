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
/**@file simplex/HEkkDualRow.cpp
 * @brief
 */
#include "simplex/HEkkDualRow.h"

#include <cassert>
#include <iostream>

#include "simplex/HSimplexDebug.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsSort.h"

using std::make_pair;
using std::pair;
using std::set;

void HEkkDualRow::setupSlice(HighsInt size) {
  workSize = size;
  workMove = &ekk_instance_.basis_.nonbasicMove_[0];
  workDual = &ekk_instance_.info_.workDual_[0];
  workRange = &ekk_instance_.info_.workRange_[0];
  work_devex_index = &ekk_instance_.info_.devex_index_[0];

  // Allocate spaces
  packCount = 0;
  packIndex.resize(workSize);
  packValue.resize(workSize);

  workCount = 0;
  workData.resize(workSize);
  analysis = &ekk_instance_.analysis_;
}

void HEkkDualRow::setup() {
  // Setup common vectors
  const HighsInt numTot =
      ekk_instance_.lp_.num_col_ + ekk_instance_.lp_.num_row_;
  setupSlice(numTot);
  workNumTotPermutation = &ekk_instance_.info_.numTotPermutation_[0];

  // deleteFreelist() is being called in Phase 1 and Phase 2 since
  // it's in updatePivots(), but create_Freelist() is only called in
  // Phase 2. Hence freeList is not initialised when freeList.empty()
  // is used in deleteFreelist(), clear freeList now.
  freeList.clear();
}

void HEkkDualRow::clear() {
  packCount = 0;
  workCount = 0;
}

void HEkkDualRow::chooseMakepack(const HVector* row, const HighsInt offset) {
  /**
   * Pack the indices and values for the row
   *
   * Offset of numCol is used when packing row_ep
   */
  const HighsInt rowCount = row->count;
  const HighsInt* rowIndex = &row->index[0];
  const double* rowArray = &row->array[0];

  for (HighsInt i = 0; i < rowCount; i++) {
    const HighsInt index = rowIndex[i];
    const double value = rowArray[index];
    packIndex[packCount] = index + offset;
    packValue[packCount++] = value;
  }
}

void HEkkDualRow::choosePossible() {
  /**
   * Determine the possible variables - candidates for CHUZC
   * TODO: Check with Qi what this is doing
   */
  const double Ta = ekk_instance_.info_.update_count < 10
                        ? 1e-9
                        : ekk_instance_.info_.update_count < 20 ? 3e-8 : 1e-6;
  const double Td = ekk_instance_.options_->dual_feasibility_tolerance;
  const HighsInt move_out = workDelta < 0 ? -1 : 1;
  workTheta = kHighsInf;
  workCount = 0;
  for (HighsInt i = 0; i < packCount; i++) {
    const HighsInt iCol = packIndex[i];
    const HighsInt move = workMove[iCol];
    const double alpha = packValue[i] * move_out * move;
    if (alpha > Ta) {
      workData[workCount++] = make_pair(iCol, alpha);
      const double relax = workDual[iCol] * move + Td;
      if (workTheta * alpha > relax) workTheta = relax / alpha;
    }
  }
}

void HEkkDualRow::chooseJoinpack(const HEkkDualRow* otherRow) {
  /**
   * Join pack of possible candidates in this row with possible
   * candidates in otherRow
   */
  const HighsInt otherCount = otherRow->workCount;
  const pair<HighsInt, double>* otherData = &otherRow->workData[0];
  copy(otherData, otherData + otherCount, &workData[workCount]);
  workCount = workCount + otherCount;
  workTheta = min(workTheta, otherRow->workTheta);
}

HighsInt HEkkDualRow::chooseFinal() {
  /**
   * Chooses the entering variable via BFRT and EXPAND
   *
   * It will
   * (1) reduce the candidates as a small collection
   * (2) choose by BFRT by going over break points
   * (3) choose final by alpha
   * (4) determine final flip variables
   */

  // 1. Reduce by large step BFRT
  analysis->simplexTimerStart(Chuzc2Clock);
  HighsInt fullCount = workCount;
  workCount = 0;
  double totalChange = 0;
  const double totalDelta = fabs(workDelta);
  double selectTheta = 10 * workTheta + 1e-7;
  for (;;) {
    for (HighsInt i = workCount; i < fullCount; i++) {
      HighsInt iCol = workData[i].first;
      double alpha = workData[i].second;
      double tight = workMove[iCol] * workDual[iCol];
      if (alpha * selectTheta >= tight) {
        swap(workData[workCount++], workData[i]);
        totalChange += workRange[iCol] * alpha;
      }
    }
    selectTheta *= 10;
    if (totalChange >= totalDelta || workCount == fullCount) break;
  }
  analysis->simplexTimerStop(Chuzc2Clock);
  HighsInt i109 = findInWorkData(109, fullCount, workData);
  // 2. Choose by small step BFRT

  bool use_quad_sort = false;
  bool use_heap_sort = false;
  // Use the quadratic cost sort for smaller values of workCount,
  // otherwise use the heap-based sort
  use_quad_sort = workCount < 100;
  const HighsInt check_iter = 26;
  if (ekk_instance_.iteration_count_ == -check_iter) {
    printf("Overruling use_quad_sort\n");
    use_quad_sort = true;
  }
  
  use_heap_sort = !use_quad_sort;
  assert(use_heap_sort || use_quad_sort);

  if (use_heap_sort) {
    // Take a copy of workData and workCount for the independent
    // heap-based code
    original_workData = workData;
    alt_workCount = workCount;
  }
  analysis->simplexTimerStart(Chuzc3Clock);
  bool choose_ok;
  if (use_quad_sort) {
    // Use the O(n^2) quadratic sort for the candidates
    analysis->simplexTimerStart(Chuzc3a0Clock);
    choose_ok = chooseFinalWorkGroupQuad();
    analysis->simplexTimerStop(Chuzc3a0Clock);
  }
  if (use_heap_sort) {
    // Use the O(n log n) heap sort for the candidates
    analysis->simplexTimerStart(Chuzc3a1Clock);
    choose_ok = chooseFinalWorkGroupHeap();
    analysis->simplexTimerStop(Chuzc3a1Clock);
  }
  if (!choose_ok) {
    analysis->simplexTimerStop(Chuzc3Clock);
    return -1;
  }
  // Make sure that there is at least one group according to sorting procedure
  if (use_quad_sort) assert((HighsInt)workGroup.size() > 1);
  if (use_heap_sort) assert((HighsInt)alt_workGroup.size() > 1);

  // 3. Choose large alpha
  analysis->simplexTimerStart(Chuzc3bClock);
  HighsInt breakIndex;
  HighsInt breakGroup;
  HighsInt alt_breakIndex;
  HighsInt alt_breakGroup;
  if (use_quad_sort)
    chooseFinalLargeAlpha(breakIndex, breakGroup, workCount, workData,
                          workGroup);
  if (use_heap_sort)
    chooseFinalLargeAlpha(alt_breakIndex, alt_breakGroup, alt_workCount,
                          sorted_workData, alt_workGroup);
  analysis->simplexTimerStop(Chuzc3bClock);

  if (!use_quad_sort) {
    // If the quadratic sort is not being used, revert to the heap
    // sort results
    breakIndex = alt_breakIndex;
    breakGroup = alt_breakGroup;
  }
  analysis->simplexTimerStart(Chuzc3cClock);

  HighsInt move_out = workDelta < 0 ? -1 : 1;
  assert(breakIndex >= 0);
  if (use_quad_sort) {
    workPivot = workData[breakIndex].first;
    workAlpha = workData[breakIndex].second * move_out * workMove[workPivot];
  } else {
    workPivot = sorted_workData[breakIndex].first;
    workAlpha =
        sorted_workData[breakIndex].second * move_out * workMove[workPivot];
  }
  if (workDual[workPivot] * workMove[workPivot] > 0) {
    workTheta = workDual[workPivot] / workAlpha;
  } else {
    workTheta = 0;
  }
  printf("breakIndex = %3d: workPivot = %4d; workAlpha = %11.4g; Dual = %11.4g; workTheta = %11.4g\n",
	 (int)breakIndex, (int)workPivot, workAlpha, workDual[workPivot], workTheta);
  analysis->simplexTimerStop(Chuzc3cClock);

  analysis->simplexTimerStart(Chuzc3dClock);

  // 4. Determine BFRT flip index: flip all
  fullCount = breakIndex;// Not used
  workCount = 0;
  const bool report = true;//ekk_instance_.iteration_count_ == check_iter;
  if (use_quad_sort) {
    for (HighsInt i = 0; i < workGroup[breakGroup]; i++) {
      const HighsInt iCol = workData[i].first;
      const HighsInt move = workMove[iCol];
      workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
    }
    if (report) printf("Final(quad): workGroup[%d] = %d with workCount = %d\n",
		       (int)breakGroup, (int)workGroup[breakGroup], (int)workCount);
  } else {
    for (HighsInt i = 0; i < alt_workGroup[breakGroup]; i++) {
      const HighsInt iCol = sorted_workData[i].first;
      const HighsInt move = workMove[iCol];
      workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
    }
    if (report) printf("Final(heap): workGroup[%d] = %d with workCount = %d\n",
		       (int)breakGroup, (int)alt_workGroup[breakGroup], (int)workCount);
    // Look at all entries of final group to see what dual
    // infeasibilities might be created
    assert(breakGroup+1 < (int)alt_workGroup.size());
    const HighsInt to_i = alt_workGroup[breakGroup+1];
    assert(to_i <= (int)sorted_workData.size());
    HighsInt num_infeasibility = 0;
    const double Td = ekk_instance_.options_->dual_feasibility_tolerance;
    for (HighsInt i = alt_workGroup[breakGroup]; i < to_i; i++) {
      const HighsInt iCol = sorted_workData[i].first;
      const double value = sorted_workData[i].second;
      const HighsInt move = workMove[iCol];
      const double dual = workDual[iCol];
      const double new_dual = dual - move_out * move * workTheta * value;
      const double new_dual_infeasibility = move*new_dual;
      const bool infeasible = new_dual_infeasibility < -Td;
      if (infeasible) num_infeasibility++;
      printf("%3d: iCol = %4d; dual = %11.4g; value = %11.4g; move = %2d; new_dual = %11.4g; infeasibility = %11.4g: %d\n",
	     (int)i, (int)iCol, dual, value, (int)move, new_dual, new_dual_infeasibility, infeasible);
      if (infeasible) workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
    }
    if (num_infeasibility) printf("Would create %d dual infeasibilities\n", (int)num_infeasibility);
  }
  if (workTheta == 0) workCount = 0;
  analysis->simplexTimerStop(Chuzc3dClock);

  analysis->simplexTimerStart(Chuzc3eClock);
  // Sort workData by .first (iCol) so that columns of A are accessed in order?
  sort(workData.begin(), workData.begin() + workCount);
  analysis->simplexTimerStop(Chuzc3eClock);
  analysis->simplexTimerStop(Chuzc3Clock);
  
  HighsInt num_infeasibility = debugChooseColumnInfeasibilities();
  if (num_infeasibility) printf("CHUZC creates %d dual infeasibilities\n", (int)num_infeasibility);
  
  if (ekk_instance_.iteration_count_ == check_iter) {
    printf("breakIndex = %d; workCount = %d; workTheta = %g; workAlpha = %g\n",
	   (int)breakIndex, (int)workCount,
	   workTheta, workAlpha);
  }
  return 0;
}

bool HEkkDualRow::chooseFinalWorkGroupQuad() {
  const double Td = ekk_instance_.options_->dual_feasibility_tolerance;
  HighsInt fullCount = workCount;
  workCount = 0;
  double totalChange = kInitialTotalChange;
  double selectTheta = workTheta;
  const double totalDelta = fabs(workDelta);
  workGroup.clear();
  workGroup.push_back(0);
  HighsInt prev_workCount = workCount;
  double prev_remainTheta = kInitialRemainTheta;
  double prev_selectTheta = selectTheta;
  HighsInt debug_num_loop = 0;

  while (selectTheta < kMaxSelectTheta) {
    double remainTheta = kInitialRemainTheta;
    debug_num_loop++;
    HighsInt debug_loop_ln = 0;
    for (HighsInt i = workCount; i < fullCount; i++) {
      HighsInt iCol = workData[i].first;
      double value = workData[i].second;
      double dual = workMove[iCol] * workDual[iCol];
      // Tight satisfy
      if (dual <= selectTheta * value) {
        swap(workData[workCount++], workData[i]);
        totalChange += value * (workRange[iCol]);
      } else if (dual + Td < remainTheta * value) {
        remainTheta = (dual + Td) / value;
      }
      debug_loop_ln++;
    }
    workGroup.push_back(workCount);

    // Update selectTheta with the value of remainTheta;
    selectTheta = remainTheta;
    // Check for no change in this loop - to prevent infinite loop
    if ((workCount == prev_workCount) && (prev_selectTheta == selectTheta) &&
        (prev_remainTheta == remainTheta)) {
      HighsInt num_var =
          ekk_instance_.lp_.num_col_ + ekk_instance_.lp_.num_row_;
      debugDualChuzcFailQuad0(*ekk_instance_.options_, workCount, workData,
                              num_var, workDual, selectTheta, remainTheta,
                              true);
      return false;
    }
    // Record the initial values of workCount, remainTheta and selectTheta for
    // the next pass through the loop - to check for infinite loop condition
    prev_workCount = workCount;
    prev_remainTheta = remainTheta;
    prev_selectTheta = selectTheta;
    if (totalChange >= totalDelta || workCount == fullCount) break;
  }
  const HighsInt check_iter = 26;
  if (ekk_instance_.iteration_count_ == check_iter) {
    printf("Quad search: workCount = %d\n", (int)workCount);
    HighsInt from_ix = 0;
    for (HighsInt iGroup=0; iGroup<workGroup.size(); iGroup++) {
      printf("Group %d\n", (int)iGroup);
      for(HighsInt ix=from_ix; ix<workGroup[iGroup]; ix++) {
	HighsInt iCol = workData[ix].first;
	double value = workData[ix].second;
	printf("%4d: iCol = %4d: value = %11.4g; workDual = %11.4g; ratio = %17.10g\n",
	       (int)ix, (int)iCol, value, workDual[iCol], fabs(workDual[iCol]/value));
      }
      from_ix = workGroup[iGroup];
    }
  }
  // Check that at least one group has been identified
  if ((HighsInt)workGroup.size() <= 1) {
    HighsInt num_var = ekk_instance_.lp_.num_col_ + ekk_instance_.lp_.num_row_;
    debugDualChuzcFailQuad1(*ekk_instance_.options_, workCount, workData,
                            num_var, workDual, selectTheta, true);
    return false;
  }
  return true;
}

bool HEkkDualRow::chooseFinalWorkGroupHeap() {
  const double Td = ekk_instance_.options_->dual_feasibility_tolerance;
  HighsInt fullCount = alt_workCount;
  double totalChange = kInitialTotalChange;
  double selectTheta = workTheta;
  const double totalDelta = fabs(workDelta);
  HighsInt heap_num_en = 0;
  std::vector<HighsInt> heap_i;
  std::vector<double> heap_v;
  heap_i.resize(fullCount + 1);
  heap_v.resize(fullCount + 1);
  HighsInt check_iCol = 109;
  HighsInt check_i = findInWorkData(check_iCol, fullCount, original_workData);
  for (HighsInt i = 0; i < fullCount; i++) {
    HighsInt iCol = original_workData[i].first;
    double value = original_workData[i].second;
    double dual = workMove[iCol] * workDual[iCol];
    double ratio = dual / value;
    if (ratio < kMaxSelectTheta) {
      heap_num_en++;
      heap_i[heap_num_en] = i;
      heap_v[heap_num_en] = ratio;
    }
  }
  maxheapsort(&heap_v[0], &heap_i[0], heap_num_en);

  alt_workCount = 0;
  alt_workGroup.clear();
  alt_workGroup.push_back(alt_workCount);
  if (heap_num_en <= 0) {
    HighsInt num_var = ekk_instance_.lp_.num_col_ + ekk_instance_.lp_.num_row_;
    // No entries in heap = > failure
    debugDualChuzcFailHeap(*ekk_instance_.options_, alt_workCount,
                           original_workData, num_var, workDual, selectTheta,
                           true);
    return false;
  }
  HighsInt this_group_first_entry = alt_workCount;
  sorted_workData.resize(heap_num_en);
  const HighsInt check_iter = 26;
  const bool report = true;//ekk_instance_.iteration_count_ == check_iter;
  for (HighsInt en = 1; en <= heap_num_en; en++) {
    HighsInt i = heap_i[en];
    HighsInt iCol = original_workData[i].first;
    double value = original_workData[i].second;
    double dual = workMove[iCol] * workDual[iCol];
    if (iCol == check_iCol) {
      if (report) printf("Found %d as entry %d\n", (int)check_iCol, (int)i);
      if (report) printf("dual = %g; selectTheta = %g; value = %g; selectTheta * value = %g; dual/value = %g\n",
	     dual,selectTheta,value,selectTheta * value, dual/value);
    }
    if (dual > selectTheta * value) {
      // Breakpoint is in the next group, so record the pointer to its
      // first entry
      alt_workGroup.push_back(alt_workCount);
      this_group_first_entry = alt_workCount;
      HighsInt alt_workGroup_size = alt_workGroup.size();
      selectTheta = (dual + Td) / value;
      if (report) printf("Next group (%2d) begins at %2d: selectTheta = %g; totalChange = %g; totalDelta = %g\n",
	     (int)(alt_workGroup_size-1), (int)alt_workCount, selectTheta, totalChange, totalDelta);
      // End loop if all permitted groups have been identified
      if (totalChange >= totalDelta) break;
    }
    // Store the breakpoint
    sorted_workData[alt_workCount].first = iCol;
    sorted_workData[alt_workCount].second = value;
    totalChange += value * (workRange[iCol]);
    alt_workCount++;
    if (report) printf("alt_workCount = %d: iCol = %3d; totalChange = %g \n",
						   (int)alt_workCount, (int)iCol, totalChange);
  }
  if (report) printf("alt_workCount = %d, this_group_first_entry = %d\n", (int)alt_workCount, (int)this_group_first_entry);
  if (alt_workCount > this_group_first_entry)
    alt_workGroup.push_back(alt_workCount);
  if (report) {
    printf("Quad search: workCount = %d\n", (int)workCount);
    HighsInt from_ix = 0;
    for (HighsInt iGroup=0; iGroup<alt_workGroup.size(); iGroup++) {
      printf("Group %d\n", (int)iGroup);
      for(HighsInt ix=from_ix; ix<alt_workGroup[iGroup]; ix++) {
	HighsInt iCol = sorted_workData[ix].first;
	double value = sorted_workData[ix].second;
	printf("%4d: iCol = %4d: value = %11.4g; workDual = %11.4g; ratio = %17.10g\n",
	       (int)ix, (int)iCol, value, workDual[iCol], fabs(workDual[iCol]/value));
      }
      from_ix = alt_workGroup[iGroup];
    }
  }
  return true;
}

void HEkkDualRow::chooseFinalLargeAlpha(
    HighsInt& breakIndex, HighsInt& breakGroup, HighsInt pass_workCount,
    const std::vector<std::pair<HighsInt, double>>& pass_workData,
    const std::vector<HighsInt>& pass_workGroup) {
  double finalCompare = 0;
  for (HighsInt i = 0; i < pass_workCount; i++)
    finalCompare = max(finalCompare, pass_workData[i].second);
  finalCompare = min(0.1 * finalCompare, 1.0);
  HighsInt countGroup = pass_workGroup.size() - 1;
  breakGroup = -1;
  breakIndex = -1;
  for (HighsInt iGroup = countGroup - 1; iGroup >= 0; iGroup--) {
    double dMaxFinal = 0;
    HighsInt iMaxFinal = -1;
    for (HighsInt i = pass_workGroup[iGroup]; i < pass_workGroup[iGroup + 1];
         i++) {
      if (dMaxFinal < pass_workData[i].second) {
        dMaxFinal = pass_workData[i].second;
        iMaxFinal = i;
      } else if (dMaxFinal == pass_workData[i].second) {
        HighsInt jCol = pass_workData[iMaxFinal].first;
        HighsInt iCol = pass_workData[i].first;
        if (workNumTotPermutation[iCol] < workNumTotPermutation[jCol]) {
          iMaxFinal = i;
        }
      }
    }

    if (pass_workData[iMaxFinal].second > finalCompare) {
      breakIndex = iMaxFinal;
      breakGroup = iGroup;
      break;
    }
  }
}

void HEkkDualRow::updateFlip(HVector* bfrtColumn) {
  double* workDual = &ekk_instance_.info_.workDual_[0];
  double dual_objective_value_change = 0;
  bfrtColumn->clear();
  const HighsInt check_iter = 26;
  const bool report = false;//ekk_instance_.iteration_count_ == check_iter;
  if (report) printf("HEkkDualRow::updateFlip workCount = %d\n", (int)workCount);
  for (HighsInt i = 0; i < workCount; i++) {
    const HighsInt iCol = workData[i].first;
    const double change = workData[i].second;
    double local_dual_objective_change = change * workDual[iCol];
    local_dual_objective_change *= ekk_instance_.cost_scale_;
    dual_objective_value_change += local_dual_objective_change;
    if (report) printf("Flip bound %5d with dual %11.4g (%2d)\n", (int)iCol, workDual[iCol],
		       ekk_instance_.basis_.nonbasicMove_[iCol]);
    ekk_instance_.flipBound(iCol);
    ekk_instance_.lp_.a_matrix_.collectAj(*bfrtColumn, iCol, change);
  }
  ekk_instance_.info_.updated_dual_objective_value +=
      dual_objective_value_change;
}

void HEkkDualRow::updateDual(double theta) {
  analysis->simplexTimerStart(UpdateDualClock);
  double* workDual = &ekk_instance_.info_.workDual_[0];
  double dual_objective_value_change = 0;
  for (HighsInt i = 0; i < packCount; i++) {
    workDual[packIndex[i]] -= theta * packValue[i];
    // Identify the change to the dual objective
    HighsInt iCol = packIndex[i];
    const double delta_dual = theta * packValue[i];
    const double local_value = ekk_instance_.info_.workValue_[iCol];
    double local_dual_objective_change =
        ekk_instance_.basis_.nonbasicFlag_[iCol] * (-local_value * delta_dual);
    local_dual_objective_change *= ekk_instance_.cost_scale_;
    dual_objective_value_change += local_dual_objective_change;
  }
  const HighsInt check_col = 2529;
  printf("Dual for workDual[%d] = %g; move[2529] = %d\n",
	 (int)check_col, workDual[check_col], (int)ekk_instance_.basis_.nonbasicMove_[check_col]);
  ekk_instance_.info_.updated_dual_objective_value +=
      dual_objective_value_change;
  analysis->simplexTimerStop(UpdateDualClock);
}

void HEkkDualRow::createFreelist() {
  freeList.clear();
  for (HighsInt i = 0;
       i < ekk_instance_.lp_.num_col_ + ekk_instance_.lp_.num_row_; i++) {
    if (ekk_instance_.basis_.nonbasicFlag_[i] &&
        highs_isInfinity(-ekk_instance_.info_.workLower_[i]) &&
        highs_isInfinity(ekk_instance_.info_.workUpper_[i]))
      freeList.insert(i);
  }
  //  debugFreeListNumEntries(ekk_instance_, freeList);
}

void HEkkDualRow::createFreemove(HVector* row_ep) {
  // TODO: Check with Qi what this is doing and why it's expensive
  if (!freeList.empty()) {
    double Ta = ekk_instance_.info_.update_count < 10
                    ? 1e-9
                    : ekk_instance_.info_.update_count < 20 ? 3e-8 : 1e-6;
    HighsInt move_out = workDelta < 0 ? -1 : 1;
    set<HighsInt>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) {
      HighsInt iCol = *sit;
      assert(iCol < ekk_instance_.lp_.num_col_);
      double alpha = ekk_instance_.lp_.a_matrix_.computeDot(*row_ep, iCol);
      if (fabs(alpha) > Ta) {
        if (alpha * move_out > 0)
          ekk_instance_.basis_.nonbasicMove_[iCol] = 1;
        else
          ekk_instance_.basis_.nonbasicMove_[iCol] = -1;
      }
    }
  }
}
void HEkkDualRow::deleteFreemove() {
  if (!freeList.empty()) {
    set<HighsInt>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) {
      HighsInt iCol = *sit;
      assert(iCol < ekk_instance_.lp_.num_col_);
      ekk_instance_.basis_.nonbasicMove_[iCol] = 0;
    }
  }
}

void HEkkDualRow::deleteFreelist(HighsInt iColumn) {
  if (!freeList.empty()) {
    if (freeList.count(iColumn)) freeList.erase(iColumn);
  }
}

void HEkkDualRow::computeDevexWeight(const HighsInt slice) {
  const bool rp_computed_edge_weight = false;
  computed_edge_weight = 0;
  for (HighsInt el_n = 0; el_n < packCount; el_n++) {
    HighsInt vr_n = packIndex[el_n];
    if (!ekk_instance_.basis_.nonbasicFlag_[vr_n]) {
      //      printf("Basic variable %" HIGHSINT_FORMAT " in packIndex is
      //      skipped\n", vr_n);
      continue;
    }
    double pv = work_devex_index[vr_n] * packValue[el_n];
    if (pv) {
      computed_edge_weight += pv * pv;
    }
  }
  if (rp_computed_edge_weight) {
    if (slice >= 0)
      printf("HEkkDualRow::computeDevexWeight: Slice %1" HIGHSINT_FORMAT
             "; computed_edge_weight = "
             "%11.4g\n",
             slice, computed_edge_weight);
  }
}

HighsInt HEkkDualRow::findInWorkData(const HighsInt iCol,
				     const HighsInt count,
				     const std::vector<std::pair<HighsInt, double>>& workData_) {
  for(HighsInt Ix=0; Ix<count; Ix++) 
    if (workData_[Ix].first == iCol) return Ix;
  return -1;
}

HighsInt HEkkDualRow::debugChooseColumnInfeasibilities() const {
  const HighsInt move_out = workDelta < 0 ? -1 : 1;
  std::vector<double> unpack_value;
  HighsLp& lp = ekk_instance_.lp_;
  unpack_value.resize(lp.num_col_ + lp.num_row_);
  for (HighsInt ix = 0; ix < packCount; ix++) 
    unpack_value[packIndex[ix]] = packValue[ix];
  HighsInt num_infeasibility = 0;
  const double Td = ekk_instance_.options_->dual_feasibility_tolerance;
  printf("HEkkDualRow::debugChooseColumnInfeasibilities: workCount = %d; move_out = %d; workTheta = %g\n",
	 (int)workCount, (int)move_out, workTheta);
  for (HighsInt i = 0; i < workCount; i++) {
    const HighsInt iCol = workData[i].first;
    const double delta = workData[i].second;
    const double value = unpack_value[iCol];
    const HighsInt move = workMove[iCol];
    const double dual = workDual[iCol];
    const double delta_dual = fabs(workTheta * value);
    const double new_dual = dual - workTheta * value;
    const double infeasibility_after_flip = -move*new_dual;
    const bool infeasible = infeasibility_after_flip < -Td;
    if (infeasible) {
      printf("%3d: iCol = %4d; dual = %11.4g; value = %11.4g; move = %2d; delta = %11.4g; new_dual = %11.4g; infeasibility = %11.4g: %d\n",
	     (int)i, (int)iCol, dual, value, (int)move, delta_dual, new_dual, infeasibility_after_flip, infeasible);

      num_infeasibility++;
    }
  }
  assert(!num_infeasibility);
  return num_infeasibility;
}

