#include "HDualRHS.h"
#include "HConst.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
using namespace std;

void HDualRHS::setup(HModel *model)
{
  workModel = model;
  workMark.resize(workModel->getNumRow());
  workIndex.resize(workModel->getNumRow());
  workArray.resize(workModel->getNumRow());
  workEdWt.assign(workModel->getNumRow(), 1);
  workEdWtFull.resize(workModel->getNumTot());
  partNum = 0;
  partSwitch = 0;
}

void HDualRHS::setup_partition(const char *filename)
{
  ifstream reader(filename);
  reader >> partNum >> partNumRow >> partNumCol >> partNumCut;
  if (partNumRow != workModel->getNumRow())
  {
    partNum = 0;
    workModel->util_reportMessage("wrong-partition-file");
    reader.close();
    return;
  }
  workPartition.assign(workModel->getNumRow(), 0);
  for (int i = 0; i < partNumRow; i++)
    reader >> workPartition[i];
  reader.close();
  partSwitch = 1;
}

void HDualRHS::choose_normal(int *chIndex)
{
  // Moved the following to the top to avoid starting the clock for a trivial call.
  if (workCount == 0)
    {
      *chIndex = -1;
      return;
    }

  // Since choose_normal calls itself, only start the clock if it's not currently running 
  bool keepTimerRunning = workModel->timer.itemStart[HTICK_CHUZR1] < 0;
  if (!keepTimerRunning) workModel->timer.recordStart(HTICK_CHUZR1);

  int random = workModel->random.intRandom();
  if (workCount < 0)
  {
    // DENSE mode
    const int numRow = -workCount;
    int randomStart = random % numRow;
    double bestMerit = 0;
    int bestIndex = -1;
    for (int section = 0; section < 2; section++)
    {
      const int start = (section == 0) ? randomStart : 0;
      const int end = (section == 0) ? numRow : randomStart;
      for (int iRow = start; iRow < end; iRow++)
      {
        if (workArray[iRow] > HSOL_CONST_ZERO)
        {
          const double myInfeas = workArray[iRow];
          const double myWeight = workEdWt[iRow];
          if (bestMerit * myWeight < myInfeas)
          {
            bestMerit = myInfeas / myWeight;
            bestIndex = iRow;
          }
        }
      }
    }
    *chIndex = bestIndex;
  }
  else
  {
    // SPARSE mode
    // Moved the following to the top to avoid starting the clock for a trivial call.
    //    if (workCount == 0)
    //    {
    //      *chIndex = -1;
    //      return;
    //    }

    int randomStart = random % workCount;
    double bestMerit = 0;
    int bestIndex = -1;
    for (int section = 0; section < 2; section++)
    {
      const int start = (section == 0) ? randomStart : 0;
      const int end = (section == 0) ? workCount : randomStart;
      for (int i = start; i < end; i++)
      {
        int iRow = workIndex[i];
        if (workArray[iRow] > HSOL_CONST_ZERO)
        {
          const double myInfeas = workArray[iRow];
          const double myWeight = workEdWt[iRow];
          if (bestMerit * myWeight < myInfeas)
          {
            bestMerit = myInfeas / myWeight;
            bestIndex = iRow;
          }
        }
      }
    }

    int createListAgain = 0;
    if (bestIndex == -1)
    {
      createListAgain = workCutoff > 0;
    }
    else if (bestMerit <= workCutoff * 0.99)
    {
      createListAgain = 1;
    }
    if (createListAgain)
    {
      create_infeasList(0);
      choose_normal(&bestIndex);
    }
    *chIndex = bestIndex;
  }
  // Since choose_normal calls itself, only stop the clock if it's not currently running 
  if (!keepTimerRunning) workModel->timer.recordFinish(HTICK_CHUZR1);
}

void HDualRHS::choose_multi_global(int *chIndex, int *chCount, int chLimit)
{
  workModel->timer.recordStart(HTICK_CHUZR1);

  for (int i = 0; i < chLimit; i++)
    chIndex[i] = -1;

  const unsigned int chooseCHECK = chLimit * 2;
  vector<pair<double, int>> setP;
  setP.reserve(chooseCHECK);

  if (workCount < 0)
  {
    // DENSE mode
    const int numRow = -workCount;
    double cutoffMerit = 0;
    for (int iRow = 0; iRow < numRow; iRow++)
    {
      if (workArray[iRow] > HSOL_CONST_ZERO)
      {
        const double myInfeas = workArray[iRow];
        const double myWeight = workEdWt[iRow];
        if (cutoffMerit * myWeight < myInfeas)
        {
          // Save
          setP.push_back(make_pair(-myInfeas / myWeight, iRow));
          // Shrink
          if (setP.size() >= chooseCHECK)
          {
            sort(setP.begin(), setP.end());
            setP.resize(chLimit);
            cutoffMerit = -setP.back().first;
          }
        }
      }
    }
  }
  else
  {
    // SPARSE Mode
    double cutoffMerit = 0;
    for (int i = 0; i < workCount; i++)
    {
      int iRow = workIndex[i];
      if (workArray[iRow] > HSOL_CONST_ZERO)
      {
        const double myInfeas = workArray[iRow];
        const double myWeight = workEdWt[iRow];
        if (cutoffMerit * myWeight < myInfeas)
        {
          // Save
          setP.push_back(make_pair(-myInfeas / myWeight, iRow));
          // Shrink
          if (setP.size() >= chooseCHECK)
          {
            sort(setP.begin(), setP.end());
            setP.resize(chLimit);
            cutoffMerit = -setP.back().first;
          }
        }
      }
    }
  }

  // Store the setP
  sort(setP.begin(), setP.end());
  if ((int)(setP.size()) > chLimit)
    setP.resize(chLimit);
  *chCount = setP.size();
  for (unsigned i = 0; i < setP.size(); i++)
    chIndex[i] = setP[i].second;
  workModel->timer.recordFinish(HTICK_CHUZR1);
}

void HDualRHS::choose_multi_HGauto(int *chIndex, int *chCount, int chLimit)
{
  // Automatically decide to use partition or not
  if (partSwitch)
    choose_multi_HGpart(chIndex, chCount, chLimit);
  else
    choose_multi_global(chIndex, chCount, chLimit);
}

void HDualRHS::choose_multi_HGpart(int *chIndex, int *chCount, int chLimit)
{
  workModel->timer.recordStart(HTICK_CHUZR1);

  // Force to use partition method, unless doesn't exist
  if (partNum != chLimit)
  {
    choose_multi_global(chIndex, chCount, chLimit);
    partSwitch = 0;
    workModel->timer.recordFinish(HTICK_CHUZR1);
    return;
  }

  // Initialise
  for (int i = 0; i < chLimit; i++)
    chIndex[i] = -1;
  *chCount = 0;

  int random = workModel->random.intRandom();
  if (workCount < 0)
  {
    // DENSE mode
    const int numRow = -workCount;
    int randomStart = random % numRow;
    vector<double> bestMerit(chLimit, 0);
    vector<int> bestIndex(chLimit, -1);
    for (int section = 0; section < 2; section++)
    {
      const int start = (section == 0) ? randomStart : 0;
      const int end = (section == 0) ? numRow : randomStart;
      for (int iRow = start; iRow < end; iRow++)
      {
        if (workArray[iRow] > HSOL_CONST_ZERO)
        {
          int iPart = workPartition[iRow];
          const double myInfeas = workArray[iRow];
          const double myWeight = workEdWt[iRow];
          if (bestMerit[iPart] * myWeight < myInfeas)
          {
            bestMerit[iPart] = myInfeas / myWeight;
            bestIndex[iPart] = iRow;
          }
        }
      }
    }
    int count = 0;
    for (int i = 0; i < chLimit; i++)
    {
      if (bestIndex[i] != -1)
      {
        chIndex[count++] = bestIndex[i];
      }
    }
    *chCount = count;
  }
  else
  {
    // SPARSE mode
    if (workCount == 0)
    {
      workModel->timer.recordFinish(HTICK_CHUZR1);
      return;
    }

    int randomStart = random % workCount;
    vector<double> bestMerit(chLimit, 0);
    vector<int> bestIndex(chLimit, -1);
    for (int section = 0; section < 2; section++)
    {
      const int start = (section == 0) ? randomStart : 0;
      const int end = (section == 0) ? workCount : randomStart;
      for (int i = start; i < end; i++)
      {
        int iRow = workIndex[i];
        if (workArray[iRow] > HSOL_CONST_ZERO)
        {
          int iPart = workPartition[iRow];
          const double myInfeas = workArray[iRow];
          const double myWeight = workEdWt[iRow];
          if (bestMerit[iPart] * myWeight < myInfeas)
          {
            bestMerit[iPart] = myInfeas / myWeight;
            bestIndex[iPart] = iRow;
          }
        }
      }
    }
    int count = 0;
    for (int i = 0; i < chLimit; i++)
    {
      if (bestIndex[i] != -1)
      {
        chIndex[count++] = bestIndex[i];
      }
    }
    *chCount = count;
  }

  workModel->timer.recordFinish(HTICK_CHUZR1);
}

void HDualRHS::update_primal(HVector *column, double theta)
{
  workModel->timer.recordStart(HTICK_UPDATE_PRIMAL);

  const int numRow = workModel->getNumRow();
  const int columnCount = column->count;
  const int *columnIndex = &column->index[0];
  const double *columnArray = &column->array[0];

  const double *baseLower = workModel->getBaseLower();
  const double *baseUpper = workModel->getBaseUpper();
  const double Tp = workModel->dblOption[DBLOPT_PRIMAL_TOL];

  double *baseValue = workModel->getBaseValue();

  bool updatePrimal_inDense = columnCount < 0 || columnCount > 0.4 * numRow;

  if (updatePrimal_inDense)
  {
    for (int iRow = 0; iRow < numRow; iRow++)
    {
      baseValue[iRow] -= theta * columnArray[iRow];
      const double value = baseValue[iRow];
      const double less = baseLower[iRow] - value;
      const double more = value - baseUpper[iRow];
      double infeas = less > Tp ? less : (more > Tp ? more : 0);
      workArray[iRow] = infeas * infeas;
    }
  }
  else
  {
    for (int i = 0; i < columnCount; i++)
    {
      int iRow = columnIndex[i];
      baseValue[iRow] -= theta * columnArray[iRow];
      const double value = baseValue[iRow];
      const double less = baseLower[iRow] - value;
      const double more = value - baseUpper[iRow];
      double infeas = less > Tp ? less : (more > Tp ? more : 0);
      workArray[iRow] = infeas * infeas;
    }
  }

  workModel->timer.recordFinish(HTICK_UPDATE_PRIMAL);
}

//This is the DSE update weight, but keep its name unchanged, just to be safe.
//Could do with changing the name "devex", though!
void HDualRHS::update_weight(HVector *column, double devex, double Kai,
                             double *dseArray)
{
  workModel->timer.recordStart(HTICK_UPDATE_WEIGHT);

  const int numRow = workModel->getNumRow();
  const int columnCount = column->count;
  const int *columnIndex = &column->index[0];
  const double *columnArray = &column->array[0];

  bool updateWeight_inDense = columnCount < 0 || columnCount > 0.4 * numRow;
  if (updateWeight_inDense)
  {
    for (int iRow = 0; iRow < numRow; iRow++)
    {
      const double val = columnArray[iRow];
      workEdWt[iRow] += val * (devex * val + Kai * dseArray[iRow]);
      if (workEdWt[iRow] < 1e-4)
        workEdWt[iRow] = 1e-4;
    }
  }
  else
  {
    for (int i = 0; i < columnCount; i++)
    {
      const int iRow = columnIndex[i];
      const double val = columnArray[iRow];
      workEdWt[iRow] += val * (devex * val + Kai * dseArray[iRow]);
      if (workEdWt[iRow] < 1e-4)
        workEdWt[iRow] = 1e-4;
    }
  }
  workModel->timer.recordFinish(HTICK_UPDATE_WEIGHT);
}
//This is the Devex update weight
void HDualRHS::update_weight_Dvx(HVector *column, double dvx_wt_o_rowOut)
{
  const bool rp_dvx = false;
  workModel->timer.recordStart(HTICK_UPDATE_WEIGHT);

  const int numRow = workModel->getNumRow();
  const int columnCount = column->count;
  const int *columnIndex = &column->index[0];
  const double *columnArray = &column->array[0];

  bool updateWeight_inDense = columnCount < 0 || columnCount > 0.4 * numRow;
  if (updateWeight_inDense)
  {
    for (int iRow = 0; iRow < numRow; iRow++)
    {
      double aa_iRow = columnArray[iRow];
      double nw_wt = max(workEdWt[iRow], dvx_wt_o_rowOut * aa_iRow * aa_iRow);
      if (abs(nw_wt - workEdWt[iRow]) > 1e-1)
        if (rp_dvx)
          printf("Row %2d: Edge weight changes from %g to %g\n", iRow, workEdWt[iRow], nw_wt);
      workEdWt[iRow] = nw_wt;
    }
  }
  else
  {
    for (int i = 0; i < columnCount; i++)
    {
      int iRow = columnIndex[i];
      double aa_iRow = columnArray[iRow];
      double nw_wt = max(workEdWt[iRow], dvx_wt_o_rowOut * aa_iRow * aa_iRow);
      if (abs(nw_wt - workEdWt[iRow]) > 1e-1)
        if (rp_dvx)
          printf("Row %2d: Edge weight changes from %g to %g\n", iRow, workEdWt[iRow], nw_wt);
      workEdWt[iRow] = nw_wt;
    }
  }
  workModel->timer.recordFinish(HTICK_UPDATE_WEIGHT);
}

void HDualRHS::update_pivots(int iRow, double value)
{
  const double *baseLower = workModel->getBaseLower();
  const double *baseUpper = workModel->getBaseUpper();
  const double Tp = workModel->dblOption[DBLOPT_PRIMAL_TOL];
  double *baseValue = workModel->getBaseValue();
  baseValue[iRow] = value;
  double pivotInfeas = 0;
  if (baseValue[iRow] < baseLower[iRow] - Tp)
    pivotInfeas = baseValue[iRow] - baseLower[iRow];
  if (baseValue[iRow] > baseUpper[iRow] + Tp)
    pivotInfeas = baseValue[iRow] - baseUpper[iRow];
  workArray[iRow] = pivotInfeas * pivotInfeas;
}

void HDualRHS::update_infeasList(HVector *column)
{
  const int columnCount = column->count;
  const int *columnIndex = &column->index[0];

  // DENSE mode: disabled
  if (workCount < 0)
    return;

  workModel->timer.recordStart(HTICK_UPDATE_PRIMAL);

  if (workCutoff <= 0)
  {
    // The regular sparse way
    for (int i = 0; i < columnCount; i++)
    {
      int iRow = columnIndex[i];
      if (workMark[iRow] == 0)
      {
        if (workArray[iRow])
        {
          workIndex[workCount++] = iRow;
          workMark[iRow] = 1;
        }
      }
    }
  }
  else
  {
    // The hyper sparse way
    for (int i = 0; i < columnCount; i++)
    {
      int iRow = columnIndex[i];
      if (workMark[iRow] == 0)
      {
        if (workArray[iRow] > workEdWt[iRow] * workCutoff)
        {
          workIndex[workCount++] = iRow;
          workMark[iRow] = 1;
        }
      }
    }
  }

  workModel->timer.recordFinish(HTICK_UPDATE_PRIMAL);
}

void HDualRHS::create_infeasArray()
{
  int numRow = workModel->getNumRow();
  const double *baseValue = workModel->getBaseValue();
  const double *baseLower = workModel->getBaseLower();
  const double *baseUpper = workModel->getBaseUpper();
  const double Tp = workModel->dblOption[DBLOPT_PRIMAL_TOL];
  for (int i = 0; i < numRow; i++)
  {
    const double value = baseValue[i];
    const double less = baseLower[i] - value;
    const double more = value - baseUpper[i];
    double infeas = less > Tp ? less : (more > Tp ? more : 0);
    workArray[i] = infeas * infeas;
  }
}

void HDualRHS::create_infeasList(double columnDensity)
{
  int numRow = workModel->getNumRow();
  double *dwork = &workEdWtFull[0];

  // 1. Build the full list
  fill_n(&workMark[0], numRow, 0);
  workCount = 0;
  workCutoff = 0;
  for (int iRow = 0; iRow < numRow; iRow++)
  {
    if (workArray[iRow])
    {
      workMark[iRow] = 1;
      workIndex[workCount++] = iRow;
    }
  }

  // 2. See if it worth to try to go sparse
  //    (Many candidates, really sparse RHS)
  if (workCount > max(numRow * 0.01, 500.0) && columnDensity < 0.05)
  {
    int icutoff = max(workCount * 0.001, 500.0);
    double maxMerit = 0;
    for (int iRow = 0, iPut = 0; iRow < numRow; iRow++)
      if (workMark[iRow])
      {
        double myMerit = workArray[iRow] / workEdWt[iRow];
        if (maxMerit < myMerit)
          maxMerit = myMerit;
        dwork[iPut++] = -myMerit;
      }
    nth_element(dwork, dwork + icutoff, dwork + workCount);
    double cutMerit = -dwork[icutoff];
    workCutoff = min(maxMerit * 0.99999, cutMerit * 1.00001);

    // Create again
    fill_n(&workMark[0], numRow, 0);
    workCount = 0;
    for (int iRow = 0; iRow < numRow; iRow++)
    {
      if (workArray[iRow] >= workEdWt[iRow] * workCutoff)
      {
        workIndex[workCount++] = iRow;
        workMark[iRow] = 1;
      }
    }

    // Reduce by drop smaller
    if (workCount > icutoff * 1.5)
    {
      // Firstly take up "icutoff" number of elements
      int fullCount = workCount;
      workCount = icutoff;
      for (int i = icutoff; i < fullCount; i++)
      {
        int iRow = workIndex[i];
        if (workArray[iRow] > workEdWt[iRow] * cutMerit)
        {
          workIndex[workCount++] = iRow;
        }
        else
        {
          workMark[iRow] = 0;
        }
      }
    }

    //        cout
    //                << "======================================================> WORK COUNT = "
    //                << workCount << "\t icutoff = " << icutoff << "\t maxMerit = "
    //                << maxMerit << "\t cutMerit = " << cutMerit << endl;
  }

  // 3. If there is still too much candidates: disable them
  if (workCount > 0.2 * numRow)
  {
    workCount = -numRow;
    workCutoff = 0;
  }
}
