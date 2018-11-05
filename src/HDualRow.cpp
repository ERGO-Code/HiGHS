#include "HDualRow.h"
#include "HModel.h"
#include "HConst.h"

#include <cassert>
#include <iostream>
using namespace std;

void HDualRow::setupSlice(HModel *model, int size)
{
  // Copy pointer
  workModel = model;
  workSize = size;
  workMove = model->getNonbasicMove();
  workDual = model->getWorkDual();
  workRange = model->getWorkRange();

  // Allocate spaces
  packCount = 0;
  packIndex.resize(workSize);
  packValue.resize(workSize);

  workCount = 0;
  workData.resize(workSize);
}

void HDualRow::setup(HModel *model)
{
  // Setup common vectors
  setupSlice(model, model->getNumTot());
  workRand = model->getWorkIntBreak();
}

void HDualRow::clear()
{
  packCount = 0;
  workCount = 0;
}

void HDualRow::choose_makepack(const HVector *row, const int offset)
{
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

void HDualRow::choose_possible()
{
  /**
   * Determine the possible variables - candidates for CHUZC
   * TODO: Check with Qi what this is doing
   */
  const double Ta = workModel->countUpdate < 10 ? 1e-9 : workModel->countUpdate < 20 ? 3e-8 : 1e-6;
  const double Td = workModel->dblOption[DBLOPT_DUAL_TOL];
  const int sourceOut = workDelta < 0 ? -1 : 1;
  workTheta = HSOL_CONST_INF;
  workCount = 0;
  for (int i = 0; i < packCount; i++)
  {
    const int iCol = packIndex[i];
    const int move = workMove[iCol];
    const double alpha = packValue[i] * sourceOut * move;
    if (alpha > Ta)
    {
      workData[workCount++] = make_pair(iCol, alpha);
      const double relax = workDual[iCol] * move + Td;
      if (workTheta * alpha > relax)
        workTheta = relax / alpha;
    }
  }
}

void HDualRow::choose_joinpack(const HDualRow *otherRow)
{
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

bool HDualRow::choose_final()
{
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
  workModel->timer.recordStart(HTICK_CHUZC2);
  int fullCount = workCount;
  workCount = 0;
  double totalChange = 0;
  double totalDelta = fabs(workDelta);
  double selectTheta = 10 * workTheta + 1e-7;
  for (;;)
  {
    for (int i = workCount; i < fullCount; i++)
    {
      int iCol = workData[i].first;
      double alpha = workData[i].second;
      double tight = workMove[iCol] * workDual[iCol];
      if (alpha * selectTheta >= tight)
      {
        swap(workData[workCount++], workData[i]);
        totalChange += workRange[iCol] * alpha;
      }
    }
    selectTheta *= 10;
    if (totalChange >= totalDelta || workCount == fullCount)
      break;
  }
  workModel->timer.recordFinish(HTICK_CHUZC2);

#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 1\n");
#endif
  // 2. Choose by small step BFRT
  workModel->timer.recordStart(HTICK_CHUZC3);
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
  while (selectTheta < 1e18)
  {
    double remainTheta = iz_remainTheta;
#ifdef HiGHSDEV
    if (rp_Choose_final) printf("Performing choose_final 2; selectTheta = %11.4g; workCount=%d; fullCount=%d\n", selectTheta, workCount, fullCount);
#endif
    for (int i = workCount; i < fullCount; i++)
    {
      int iCol = workData[i].first;
      double value = workData[i].second;
      double dual = workMove[iCol] * workDual[iCol];
#ifdef HiGHSDEV
      if (rp_Choose_final) printf("iCol=%4d; v=%11.4g; d=%11.4g |", iCol, value, dual);
#endif
      // Tight satisfy
#ifdef HiGHSDEV
      if (rp_Choose_final) printf(" %11.4g = dual ?<=? sTh * v = %11.4g; workCount=%2d", dual, selectTheta * value, workCount);
#endif
      if (dual <= selectTheta * value)
      {
        swap(workData[workCount++], workData[i]);
        totalChange += value * (workRange[iCol]);
      }
      else if (dual + Td < remainTheta * value)
      {
        remainTheta = (dual + Td) / value;
      }
#ifdef HiGHSDEV
      if (rp_Choose_final) printf(": totCg=%11.4g; rmTh=%11.4g\n", totalChange, remainTheta);
#endif
    }
    workGroup.push_back(workCount);
    // Update selectTheta with the value of remainTheta;
    selectTheta = remainTheta;
    //Check for no change in this loop - to prevent infinite loop
    if ((workCount == prev_workCount) && (prev_selectTheta == selectTheta) && (prev_remainTheta == remainTheta)) {
#ifdef HiGHSDEV
      printf("In choose_final: No change in loop 2 so return error\n");
      double workDataNorm = 0;
      double dualNorm = 0;
      for (int i = 0; i < workCount; i++) {
	int iCol = workData[i].first;
	double value = workData[i].second;
	workDataNorm += value*value;
	value = workDual[iCol];
	dualNorm += value*value;
      }
      workDataNorm += sqrt(workDataNorm);
      dualNorm += sqrt(dualNorm);
      printf("   workCount = %d; selectTheta=%g; remainTheta=%g\n", workCount, selectTheta, remainTheta);
      printf(   "workDataNorm = %g; dualNorm = %g\n", workDataNorm, dualNorm);      
#endif
      return true;
    }
    // Record the initial values of workCount, remainTheta and selectTheta for the next pass through the loop
    prev_workCount = workCount;
    prev_remainTheta = remainTheta;
    prev_selectTheta = selectTheta;
    if (totalChange >= totalDelta || workCount == fullCount)
      break;
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
  for (int iGroup = countGroup - 1; iGroup >= 0; iGroup--)
  {
    double dMaxFinal = 0;
    int iMaxFinal = -1;
    for (int i = workGroup[iGroup]; i < workGroup[iGroup + 1]; i++)
    {
      if (dMaxFinal < workData[i].second)
      {
        dMaxFinal = workData[i].second;
        iMaxFinal = i;
      }
      else if (dMaxFinal == workData[i].second)
      {
        int jCol = workData[iMaxFinal].first;
        int iCol = workData[i].first;
        if (workRand[iCol] < workRand[jCol])
        {
          iMaxFinal = i;
        }
      }
    }

    if (workData[iMaxFinal].second > finalCompare)
    {
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
  if (workDual[workPivot] * workMove[workPivot] > 0)
    {
      workTheta = workDual[workPivot] / workAlpha;
    }
  else
    {
      workTheta = 0;
    }

  // 4. Determine BFRT flip index: flip all
  fullCount = breakIndex;
  workCount = 0;
  for (int i = 0; i < workGroup[breakGroup]; i++)
  {
    const int iCol = workData[i].first;
    const int move = workMove[iCol];
    workData[workCount++] = make_pair(iCol, move * workRange[iCol]);
  }
  if (workTheta == 0)
    workCount = 0;
  sort(workData.begin(), workData.begin() + workCount);
  workModel->timer.recordFinish(HTICK_CHUZC3);
#ifdef HiGHSDEV
  if (rp_Choose_final) printf("Completed  choose_final 4\n");
#endif
  return false;
}

void HDualRow::update_flip(HVector *bfrtColumn)
{
  bfrtColumn->clear();
  for (int i = 0; i < workCount; i++)
  {
    const int iCol = workData[i].first;
    const double change = workData[i].second;
    workModel->flipBound(iCol);
    workModel->getMatrix()->collect_aj(*bfrtColumn, iCol, change);
  }
}

void HDualRow::update_dual(double theta)
{
  workModel->timer.recordStart(HTICK_UPDATE_DUAL);
  double *workDual = workModel->getWorkDual();
  for (int i = 0; i < packCount; i++) {
    workDual[packIndex[i]] -= theta * packValue[i];
    // Identify the change to the dual objective
    //JAJH10/10
    /*
    int iCol = packIndex[i];
    double dlDual = theta * packValue[i];
    double iColWorkValue = workModel->workValue[iCol];
    double dlDuObj = -iColWorkValue * dlDual;
    //    dlDuObj *= costScale;
    workModel->objective += dlDuObj;
    printf("Column %2d: Fg = %2d; dlDual = %11.4g; iColWorkValue = %11.4g; dlDuObj = %11.4g: DuObj = %11.4g\n",
    	   iCol, workModel->nonbasicFlag[i], dlDual, iColWorkValue, dlDuObj, workModel->objective);
    */
  }
  workModel->timer.recordFinish(HTICK_UPDATE_DUAL);
}

void HDualRow::create_Freelist()
{
  freeList.clear();
  const int *nonbasicFlag = workModel->getNonbasicFlag();
  int ckFreeListSize = 0;
  for (int i = 0; i < workModel->getNumTot(); i++)
  {
    if (nonbasicFlag[i] && workRange[i] > 1.5 * HSOL_CONST_INF) {
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
  //  printf("Create Freelist %d:%d has size %d (%3d%%)\n", freeListSa, freeListE, freeListSize, 100*freeListSize/workModel->getNumTot());
}

void HDualRow::create_Freemove(HVector *row_ep)
{
  // TODO: Check with Qi what this is doing and why it's expensive
  if (!freeList.empty())
  {
    double Ta = workModel->countUpdate < 10 ? 1e-9 : workModel->countUpdate < 20 ? 3e-8 : 1e-6;
    int sourceOut = workDelta < 0 ? -1 : 1;
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++)
    {
      int iCol = *sit;
      assert(iCol < workModel->getNumCol());
      double alpha = workModel->getMatrix()->compute_dot(*row_ep, iCol);
      if (fabs(alpha) > Ta)
      {
        if (alpha * sourceOut > 0)
          workModel->getNonbasicMove()[iCol] = 1;
        else
          workModel->getNonbasicMove()[iCol] = -1;
      }
    }
  }
}
void HDualRow::delete_Freemove()
{
  if (!freeList.empty())
  {
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++)
    {
      int iCol = *sit;
      assert(iCol < workModel->getNumCol());
      workModel->getNonbasicMove()[iCol] = 0;
    }
  }
}

void HDualRow::delete_Freelist(int iColumn)
{
  if (!freeList.empty())
  {
    if (freeList.count(iColumn))
      freeList.erase(iColumn);
    //  int freeListSa = *freeList.begin();
    //  int freeListE = *freeList.end();
    int ckFreeListSize = 0;
    set<int>::iterator sit;
    for (sit = freeList.begin(); sit != freeList.end(); sit++) ckFreeListSize++;
    freeListSize = *freeList.end();
    if (freeListSize != ckFreeListSize) {
      printf("!! STRANGE: freeListSize != ckFreeListSize\n");
    }
    //  printf("Update Freelist %d:%d has size %d (%3d%%)\n", freeListSa, freeListE, freeListSize, 100*freeListSize/workModel->getNumTot());
    //  if (freeList.empty()) {
    //    printf("Empty  Freelist\n");
    //  } else {
    //    printf("\n");
    //  }
  } else {
    if (freeListSize > 0)
      printf("!! STRANGE: Empty Freelist has size %d\n", freeListSize);
  }
}
