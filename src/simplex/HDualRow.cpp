/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDualRow.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HDualRow.h"
#include "HighsModelObject.h"

#include <cassert>
#include <iostream>

#include "HConst.h"
#include "HModel.h"
#include "SimplexTimer.h"
#include "HVector.h"

using std::make_pair;
using std::pair;
using std::set;

void HDualRow::setupSlice(HighsModelObject& highs_model_object, int size) {
  // Copy pointer
  workHMO = highs_model_object;
  HModel *model;
  model = &workHMO->hmodel_[0];
  workModel = model;
  workSize = size;
  workMove = &workHMO->basis_.nonbasicMove_[0];
  workDual = &workHMO->simplex_.workDual_[0];
  workRange = &workHMO->simplex_.workRange_[0];

  // Allocate spaces
  packCount = 0;
  packIndex.resize(workSize);
  packValue.resize(workSize);

  workCount = 0;
  workData.resize(workSize);
}

void HDualRow::setup() {
  // Setup common vectors
  HModel *model;
  model = &highs_model_object->hmodel_[0];
  const int numTot = model->lp_scaled_->numCol_ + model->lp_scaled_->numRow_;
  setupSlice(highs_model_object, numTot);
  workColPermutation = &model->colPermutation[0];

  
 // delete_Freelist() is being called in Phase 1 and Phase 2 since
 // it's in updatePivots(), but create_Freelist() is only called in
 // Phase 2. Hence freeList and freeListSize are not initialised when
 // freeList.empty() is used to identify that freeListSize should be
 // tested for zero. Suddenly freeListSize is 1212631365 rather than
 // zero when uninitialised, triggering a warning. So, let's set
 // clear freeList and set freeListSize = 0.
  freeList.clear();
  freeListSize = 0;

}

void HDualRow::clear() {
  packCount = 0;
  workCount = 0;
}

void HDualRow::choose_makepack(const HVector *row, const int offset) {
  /**
   * Pack the indices and values for the row
   *
   * Offset of numCol is used when packing row_ep
   */
  const int rowCount = row->count;
  const int *rowIndex = &row->index[0];
  const double *rowArray = &row->array[0];
  const double *rowPackValue = &row->packValue[0];
  const int rowPWd = row->pWd;

  if (rowPWd == row->dfSparseDaStr) {
    for (int i = 0; i < rowCount; i++) {
      const int index = rowIndex[i];
      const double value = rowArray[index];
      packIndex[packCount] = index + offset;
      packValue[packCount++] = value;
    }
  } else if (rowPWd >= row->p0SparseDaStr) {
    for (int i = 0; i < rowCount; i++) {
      const int index = rowIndex[i];
      const double value = rowPackValue[i];
      packIndex[packCount] = index + offset;
      packValue[packCount++] = value;
    }
  } else {
    printf("HDualRow::choose_makepack: Cannot handle rowPWd = %d\n", rowPWd);
  }
}

void HDualRow::choose_possible() {
  /**
   * Determine the possible variables - candidates for CHUZC
   * TODO: Check with Qi what this is doing
   */
  const double Ta = workModel->countUpdate < 10
                        ? 1e-9
                        : workModel->countUpdate < 20 ? 3e-8 : 1e-6;
  const double Td = workModel->dblOption[DBLOPT_DUAL_TOL];
  const int sourceOut = workDelta < 0 ? -1 : 1;
  workTheta = HIGHS_CONST_INF;
  workCount = 0;
  for (int i = 0; i < packCount; i++) {
    const int iCol = packIndex[i];
    const int move = workMove[iCol];
    const double alpha = packValue[i] * sourceOut * move;
    if (alpha > Ta) {
      workData[workCount++] = make_pair(iCol, alpha);
      const double relax = workDual[iCol] * move + Td;
      if (workTheta * alpha > relax) workTheta = relax / alpha;
    }
  }
}

void HDualRow::choose_joinpack(const HDualRow *otherRow) {
  /**
   * Join pack of possible candidates in this row with possible
   * candidates in otherRow
   */
  const int otherCount = otherRow->workCount;
  const pair<int, double> *otherData = &otherRow->workData[0];
  copy(otherData, otherData + otherCount, &workData[workCount]);
  workCount = workCount + otherCount;
  workTheta = min(workTheta, otherRow->workTheta);
}

bool HDualRow::choose_final() {
  HighsTimer &timer = workHMO->timer_;
  HighsSimplexInfo &simplex = workHMO->simplex_;
  /**
   * Chooses the entering variable via BFRT and EXPAND
   *
   * It will
   * (1) reduce the candidates as a small collection
   * (2) choose by BFRT by going over break points
   * (3) choose final by alpha
   * (4) determine final flip variables
   */

#ifdef HiGHSDEV
  bool rp_Choose_final = false;
  //   rp_Choose_final = true;
#endif
  // 1. Reduce by large step BFRT
  timer.start(simplex.clock_[Chuzc2Clock]);
  int fullCount = workCount;
  workCount = 0;
  double totalChange = 0;
  double totalDelta = fabs(workDelta);
  double selectTheta = 10 * workTheta + 1e-7;
  for (;;) {
    for (int i = workCount; i < fullCount; i++) {
      int iCol = workData[i].first;
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
  timer.stop(simplex.clock_[Chuzc2Clock]);

#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 1\n");
#endif
  // 2. Choose by small step BFRT
  timer.start(simplex.clock_[Chuzc3Clock]);
  const double Td = workModel->dblOption[DBLOPT_DUAL_TOL];
  fullCount = workCount;
  workCount = 0;
  totalChange = 1e-12;
  selectTheta = workTheta;
  workGroup.clear();
  workGroup.push_back(0);
  const double iz_remainTheta = 1e100;
  int prev_workCount = workCount;
  double prev_remainTheta = iz_remainTheta;
  double prev_selectTheta = selectTheta;
  while (selectTheta < 1e18) {
    double remainTheta = iz_remainTheta;
#ifdef HiGHSDEV
    if (rp_Choose_final)
      printf(
          "Performing choose_final 2; selectTheta = %11.4g; workCount=%d; "
          "fullCount=%d\n",
          selectTheta, workCount, fullCount);
#endif
    for (int i = workCount; i < fullCount; i++) {
      int iCol = workData[i].first;
      double value = workData[i].second;
      double dual = workMove[iCol] * workDual[iCol];
#ifdef HiGHSDEV
      if (rp_Choose_final)
        printf("iCol=%4d; v=%11.4g; d=%11.4g |", iCol, value, dual);
#endif
        // Tight satisfy
#ifdef HiGHSDEV
      if (rp_Choose_final)
        printf(" %11.4g = dual ?<=? sTh * v = %11.4g; workCount=%2d", dual,
               selectTheta * value, workCount);
#endif
      if (dual <= selectTheta * value) {
        swap(workData[workCount++], workData[i]);
        totalChange += value * (workRange[iCol]);
      } else if (dual + Td < remainTheta * value) {
        remainTheta = (dual + Td) / value;
      }
#ifdef HiGHSDEV
      if (rp_Choose_final)
        printf(": totCg=%11.4g; rmTh=%11.4g\n", totalChange, remainTheta);
#endif
    }
    workGroup.push_back(workCount);
    // Update selectTheta with the value of remainTheta;
    selectTheta = remainTheta;
    // Check for no change in this loop - to prevent infinite loop
    if ((workCount == prev_workCount) && (prev_selectTheta == selectTheta) &&
        (prev_remainTheta == remainTheta)) {
#ifdef HiGHSDEV
      printf("In choose_final: No change in loop 2 so return error\n");
      double workDataNorm = 0;
      double dualNorm = 0;
      for (int i = 0; i < workCount; i++) {
        int iCol = workData[i].first;
        double value = workData[i].second;
        workDataNorm += value * value;
        value = workDual[iCol];
        dualNorm += value * value;
      }
      workDataNorm += sqrt(workDataNorm);
      dualNorm += sqrt(dualNorm);
      printf("   workCount = %d; selectTheta=%g; remainTheta=%g\n", workCount,
             selectTheta, remainTheta);
      printf("workDataNorm = %g; dualNorm = %g\n", workDataNorm, dualNorm);
#endif
      return true;
    }
    // Record the initial values of workCount, remainTheta and selectTheta for
    // the next pass through the loop
    prev_workCount = workCount;
    prev_remainTheta = remainTheta;
    prev_selectTheta = selectTheta;
    if (totalChange >= totalDelta || workCount == fullCount) break;
  }

#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 2\n");
#endif
  // 3. Choose large alpha
  double finalCompare = 0;
  for (int i = 0; i < workCount; i++)
    finalCompare = max(finalCompare, workData[i].second);
  finalCompare = min(0.1 * finalCompare, 1.0);
  int countGroup = workGroup.size() - 1;
  int breakGroup = -1;
  int breakIndex = -1;
  for (int iGroup = countGroup - 1; iGroup >= 0; iGroup--) {
    double dMaxFinal = 0;
    int iMaxFinal = -1;
    for (int i = workGroup[iGroup]; i < workGroup[iGroup + 1]; i++) {
      if (dMaxFinal < workData[i].second) {
        dMaxFinal = workData[i].second;
        iMaxFinal = i;
      } else if (dMaxFinal == workData[i].second) {
        int jCol = workData[iMaxFinal].first;
        int iCol = workData[i].first;
        if (workColPermutation[iCol] < workColPermutation[jCol]) {
          iMaxFinal = i;
        }
      }
    }

    if (workData[iMaxFinal].second > finalCompare) {
      breakIndex = iMaxFinal;
      breakGroup = iGroup;
      break;
    }
  }

#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 3\n");
#endif
  int sourceOut = workDelta < 0 ? -1 : 1;
  workPivot = workData[breakIndex].first;
  workAlpha = workData[breakIndex].second * sourceOut * workMove[workPivot];
  if (workDual[workPivot] * workMove[workPivot] > 0) {
    workTheta = workDual[workPivot] / workAlpha;
  } else {
    workTheta = 0;
  }

  // 4. Determine BFRT flip index: flip all
  fullCount = breakIndex;
  workCount = 0;
  for (int i = 0; i < workGroup[breakGroup]; i++) {
    const int iCol = workData[i].first;
    const int move = workMove[iCol];
    workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
  }
  if (workTheta == 0) workCount = 0;
  sort(workData.begin(), workData.begin() + workCount);
  timer.stop(simplex.clock_[Chuzc3Clock]);
#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 4\n");
#endif
  return false;
}

void HDualRow::update_flip(HVector *bfrtColumn) {
  //  workModel->checkDualObjectiveValue("Before update_flip");
  double *workDual = &workHMO->simplex_.workDual_[0];//
  //  double *workLower = &workHMO->simplex_.workLower_[0];
  //  double *workUpper = &workHMO->simplex_.workUpper_[0];
  //  double *workValue = &workHMO->simplex_.workValue_[0];
  double dualObjectiveValueChange = 0;
  bfrtColumn->clear();
  for (int i = 0; i < workCount; i++) {
    const int iCol = workData[i].first;
    const double change = workData[i].second;

    double lcDualObjectiveValueChange = change*workDual[iCol];
    //    printf("%6d: [%11.4g, %11.4g, %11.4g], (%11.4g) DlObj = %11.4g dualObjectiveValueChange = %11.4g\n",
    //	   iCol, workLower[iCol], workValue[iCol], workUpper[iCol], change, lcDualObjectiveValueChange, dualObjectiveValueChange);
    dualObjectiveValueChange += lcDualObjectiveValueChange;
    workModel->flipBound(iCol);
    workModel->matrix_->collect_aj(*bfrtColumn, iCol, change);
  }
  workModel->updatedDualObjectiveValue += dualObjectiveValueChange;
  //  workModel->checkDualObjectiveValue("After  update_flip");
}

void HDualRow::update_dual(double theta, int columnOut) {
  HighsTimer &timer = workHMO->timer_;
  HighsSimplexInfo &simplex = workHMO->simplex_;
  //  workModel->checkDualObjectiveValue("Before update_dual");
  timer.start(simplex.clock_[UpdateDualClock]);
  double *workDual = &workHMO->simplex_.workDual_[0];
  //  int columnOut_i = -1;
  for (int i = 0; i < packCount; i++) {
    workDual[packIndex[i]] -= theta * packValue[i];
    // Identify the change to the dual objective
    int iCol = packIndex[i];
    //    if (iCol == columnOut) columnOut_i = i;
    double dlDual = theta * packValue[i];
    double iColWorkValue = workHMO->simplex_.workValue_[iCol];
    double dlDuObj = workHMO->basis_.nonbasicFlag_[iCol] * (-iColWorkValue * dlDual);
    dlDuObj *= workHMO->scale_.cost_;
    workModel->updatedDualObjectiveValue += dlDuObj;
  }
  /*
  // Apply correction because the updated dual objective value may
  // contain a rogue contribution corresponding to the leaving column
  double duObjCorrection = 0;
  if (columnOut_i >= 0) {
    //    double nonUnitPackValue = abs(packValue[columnOut_i] - 1.0);
    //    if (nonUnitPackValue > 1e-8) printf("STRANGE: packValue[columnOut_i] = %11.4g has nonUnitPackValue = %11.4g\n", packValue[columnOut_i], nonUnitPackValue);
    int iCol = columnOut;
    double dlDual = theta * packValue[columnOut_i];
    double iColWorkValue = workModel->workValue[iCol];
    double dlDuObj = -iColWorkValue * dlDual;
    duObjCorrection = dlDuObj;
    //    if (abs(duObjCorrection)>1e-8) printf("Dual objective correction (columnOut_i>=0) = %11.4g\n", duObjCorrection);
  }
  
  double duObjEr = workModel->checkDualObjectiveValue("After  update_dual");
  double rlvErDen = max(1.0, abs(workModel->dualObjectiveValue));
  double rlvDuObjEr = abs(duObjEr)/rlvErDen;
  if (rlvDuObjEr > 1e-8) {
    if (columnOut_i < 0) printf("ERROR: columnOut_i = %d in update_dual\n", columnOut_i);
    double nwDuObjEr = duObjEr+duObjCorrection;
    printf("theta = %11.4g, duObjEr = %11.4g; duObjCorrection = %11.4g; nwDuObjEr = %11.4g\n",
	   theta, duObjEr, duObjCorrection, nwDuObjEr);
    if (abs(nwDuObjEr)/rlvErDen > 1e-12)
      printf("ERROR:  duObjEr+duObjCorrection = %11.4g + %11.4g = %11.4g\n", duObjEr, duObjCorrection, nwDuObjEr);
    for (int i = 0; i < packCount; i++) {
      int iCol = packIndex[i];
      double dlDual = theta * packValue[i];
      double iColWorkValue = workModel->workValue[iCol];
      double dlDuObj = -iColWorkValue * dlDual;
      dlDuObj *= workModel->scale.cost_;
      if (!workModel->nonbasicFlag[iCol])
	printf("Column %5d: packValue = %11.4g Fg = %2d; dlDual = %11.4g; iColWorkValue = %11.4g; dlDuObj = %11.4g: DuObj = %11.4g\n", 
	       iCol, packValue[i], workModel->nonbasicFlag[iCol], dlDual, iColWorkValue, dlDuObj, workModel->dualObjectiveValue);
    }
  }
  */
  timer.stop(simplex.clock_[UpdateDualClock]);
}

void HDualRow::create_Freelist() {
  freeList.clear();
  const int *nonbasicFlag = &workHMO->basis_.nonbasicFlag_[0];
  int ckFreeListSize = 0;
  const int numTot = workModel->lp_scaled_->numCol_ + workModel->lp_scaled_->numRow_;
  for (int i = 0; i < numTot; i++) {
    if (nonbasicFlag[i] && workRange[i] > 1.5 * HIGHS_CONST_INF) {
      freeList.insert(i);
      ckFreeListSize++;
    }
  }
  //  int freeListSa = *freeList.begin();
  //  int freeListE = *freeList.end();
  freeListSize = *freeList.end();
  if (freeListSize != ckFreeListSize) {
    printf("!! STRANGE: freeListSize != ckFreeListSize\n");
  }
  // const int numTot = workModel->lp_scaled_->numCol_ + workModel->lp_scaled_->numRow_;
  //  printf("Create Freelist %d:%d has size %d (%3d%%)\n", freeListSa,
  //  freeListE, freeListSize, 100*freeListSize/numTot);
}

void HDualRow::create_Freemove(HVector *row_ep) {
  // TODO: Check with Qi what this is doing and why it's expensive
  if (!freeList.empty()) {
    double Ta = workModel->countUpdate < 10
                    ? 1e-9
                    : workModel->countUpdate < 20 ? 3e-8 : 1e-6;
    int sourceOut = workDelta < 0 ? -1 : 1;
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) {
      int iCol = *sit;
      assert(iCol < workModel->lp_scaled_->numCol_);
      double alpha = workModel->matrix_->compute_dot(*row_ep, iCol);
      if (fabs(alpha) > Ta) {
        if (alpha * sourceOut > 0)
          workHMO->basis_.nonbasicMove_[iCol] = 1;
        else
          workHMO->basis_.nonbasicMove_[iCol] = -1;
      }
    }
  }
}
void HDualRow::delete_Freemove() {
  if (!freeList.empty()) {
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) {
      int iCol = *sit;
      assert(iCol < workModel->lp_scaled_->numCol_);
      workHMO->basis_.nonbasicMove_[iCol] = 0;
    }
  }
}

void HDualRow::delete_Freelist(int iColumn) {
  if (!freeList.empty()) {
    if (freeList.count(iColumn)) freeList.erase(iColumn);
    //  int freeListSa = *freeList.begin();
    //  int freeListE = *freeList.end();
    int ckFreeListSize = 0;
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) ckFreeListSize++;
    freeListSize = *freeList.end();
    if (freeListSize != ckFreeListSize) {
      printf("!! STRANGE: freeListSize != ckFreeListSize\n");
    }
    // const int numTot = workModel->lp_scaled_->numCol_ + workModel->lp_scaled_->numRow_;
    //  printf("Update Freelist %d:%d has size %d (%3d%%)\n", freeListSa,
    //  freeListE, freeListSize, 100*freeListSize/numTot); if
    //  (freeList.empty()) {
    //    printf("Empty  Freelist\n");
    //  } else {
    //    printf("\n");
    //  }
  } else {
    if (freeListSize > 0)
      printf("!! STRANGE: Empty Freelist has size %d\n", freeListSize);
  }
}
