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
/**@file simplex/HEkkDualRHS.cpp
 * @brief
 */
#include "HEkkDualRHS.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>

#include "lp_data/HConst.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HVector.h"
#include "simplex/SimplexTimer.h"

using std::fill_n;
using std::make_pair;
using std::nth_element;
using std::pair;

void HEkkDualRHS::setup() {
  const HighsInt numRow = ekk_instance_.simplex_lp_.numRow_;
  const HighsInt numTot =
      ekk_instance_.simplex_lp_.numCol_ + ekk_instance_.simplex_lp_.numRow_;
  workMark.resize(numRow);
  workIndex.resize(numRow);
  work_infeasibility.resize(numRow);
  workEdWt.assign(numRow, 1);
  workEdWtFull.resize(numTot);
  partNum = 0;
  partSwitch = 0;
  analysis = &ekk_instance_.analysis_;
}

void HEkkDualRHS::chooseNormal(HighsInt* chIndex) {
  // Moved the following to the top to avoid starting the clock for a trivial
  // call. NB Must still call HighsInt to maintain sequence of random numbers
  // for code reproducibility!! Never mind if we're not timing the random number
  // call!!
  // HighsInt random = ekk_instance_.random_.integer();
  if (workCount == 0) {
    *chIndex = -1;
    return;
  }

  // Since chooseNormal calls itself, only start the clock if it's not
  // currently running
  bool keep_timer_running = analysis->simplexTimerRunning(ChuzrDualClock);
  //      timer.clock_start[simplex_info.clock_[ChuzrDualClock]] < 0;
  if (!keep_timer_running) {
    analysis->simplexTimerStart(ChuzrDualClock);
  }

  if (workCount < 0) {
    // DENSE mode
    const HighsInt numRow = -workCount;
    HighsInt randomStart = ekk_instance_.random_.integer(numRow);
    double bestMerit = 0;
    HighsInt bestIndex = -1;
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? numRow : randomStart;
      for (HighsInt iRow = start; iRow < end; iRow++) {
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = workEdWt[iRow];
          //	  printf("Dense: Row %4" HIGHSINT_FORMAT " weight = %g\n", iRow,
          // myWeight);
          if (bestMerit * myWeight < myInfeas) {
            bestMerit = myInfeas / myWeight;
            bestIndex = iRow;
          }
        }
      }
    }
    *chIndex = bestIndex;
  } else {
    // SPARSE mode
    // Moved the following to the top to avoid starting the clock for a trivial
    // call.
    //    if (workCount == 0)
    //    {
    //      *chIndex = -1;
    //      return;
    //    }

    HighsInt randomStart = ekk_instance_.random_.integer(workCount);
    double bestMerit = 0;
    HighsInt bestIndex = -1;
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? workCount : randomStart;
      for (HighsInt i = start; i < end; i++) {
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = workEdWt[iRow];
          /*
          const double myMerit = myInfeas / myWeight;
          printf("CHUZR: iRow = %6" HIGHSINT_FORMAT "; Infeas = %11.4g; Weight =
          %11.4g; Merit = %11.4g\n", iRow, myInfeas, myWeight, myMerit);
          */
          if (bestMerit * myWeight < myInfeas) {
            bestMerit = myInfeas / myWeight;
            bestIndex = iRow;
          }
        }
      }
    }

    HighsInt createListAgain = 0;
    if (bestIndex == -1) {
      createListAgain = workCutoff > 0;
    } else if (bestMerit <= workCutoff * 0.99) {
      createListAgain = 1;
    }
    if (createListAgain) {
      createInfeasList(0);
      chooseNormal(&bestIndex);
    }
    *chIndex = bestIndex;
  }
  // Since chooseNormal calls itself, only stop the clock if it's not currently
  // running
  if (!keep_timer_running) analysis->simplexTimerStop(ChuzrDualClock);
}

void HEkkDualRHS::chooseMultiGlobal(HighsInt* chIndex, HighsInt* chCount,
                                    HighsInt chLimit) {
  analysis->simplexTimerStart(ChuzrDualClock);

  for (HighsInt i = 0; i < chLimit; i++) chIndex[i] = -1;

  const HighsUInt chooseCHECK = chLimit * 2;
  vector<pair<double, int>> setP;
  setP.reserve(chooseCHECK);

  if (workCount < 0) {
    // DENSE mode
    const HighsInt numRow = -workCount;
    HighsInt randomStart = ekk_instance_.random_.integer(numRow);
    double cutoffMerit = 0;
    // Now
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? numRow : randomStart;
      for (HighsInt iRow = start; iRow < end; iRow++) {
        // Was
        //    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
        // Continue
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = workEdWt[iRow];
          if (cutoffMerit * myWeight < myInfeas) {
            // Save
            setP.push_back(make_pair(-myInfeas / myWeight, iRow));
            // Shrink
            if (setP.size() >= chooseCHECK) {
              sort(setP.begin(), setP.end());
              setP.resize(chLimit);
              cutoffMerit = -setP.back().first;
            }
          }
        }
      }
    }
  } else {
    // SPARSE Mode
    HighsInt randomStart;
    if (workCount) {
      randomStart = ekk_instance_.random_.integer(workCount);
    } else {
      // workCount = 0
      randomStart = 0;
    }
    double cutoffMerit = 0;
    // Now
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? workCount : randomStart;
      for (HighsInt i = start; i < end; i++) {
        // Was
        //    for (HighsInt i = 0; i < workCount; i++) {
        // Continue
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > kHighsZero) {
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = workEdWt[iRow];
          /*
          const double myMerit = myInfeas / myWeight;
          printf("CHUZR: iRow = %6" HIGHSINT_FORMAT "; Infeas = %11.4g; Weight =
          %11.4g; Merit = %11.4g\n", iRow, myInfeas, myWeight, myMerit);
          */
          if (cutoffMerit * myWeight < myInfeas) {
            // Save
            setP.push_back(make_pair(-myInfeas / myWeight, iRow));
            // Shrink
            if (setP.size() >= chooseCHECK) {
              sort(setP.begin(), setP.end());
              setP.resize(chLimit);
              cutoffMerit = -setP.back().first;
            }
          }
        }
      }
    }
  }

  // Store the setP
  sort(setP.begin(), setP.end());
  if ((HighsInt)(setP.size()) > chLimit) setP.resize(chLimit);
  *chCount = setP.size();
  for (unsigned i = 0; i < setP.size(); i++) chIndex[i] = setP[i].second;
  analysis->simplexTimerStop(ChuzrDualClock);
}

void HEkkDualRHS::chooseMultiHyperGraphAuto(HighsInt* chIndex,
                                            HighsInt* chCount,
                                            HighsInt chLimit) {
  // Automatically decide to use partition or not
  if (partSwitch)
    chooseMultiHyperGraphPart(chIndex, chCount, chLimit);
  else
    chooseMultiGlobal(chIndex, chCount, chLimit);
}

void HEkkDualRHS::chooseMultiHyperGraphPart(HighsInt* chIndex,
                                            HighsInt* chCount,
                                            HighsInt chLimit) {
  analysis->simplexTimerStart(ChuzrDualClock);

  // Force to use partition method, unless doesn't exist
  if (partNum != chLimit) {
    chooseMultiGlobal(chIndex, chCount, chLimit);
    partSwitch = 0;
    analysis->simplexTimerStop(ChuzrDualClock);
    return;
  }

  // Initialise
  for (HighsInt i = 0; i < chLimit; i++) chIndex[i] = -1;
  *chCount = 0;

  if (workCount < 0) {
    // DENSE mode
    const HighsInt numRow = -workCount;
    HighsInt randomStart = ekk_instance_.random_.integer(numRow);
    vector<double> bestMerit(chLimit, 0);
    vector<HighsInt> bestIndex(chLimit, -1);
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? numRow : randomStart;
      for (HighsInt iRow = start; iRow < end; iRow++) {
        if (work_infeasibility[iRow] > kHighsZero) {
          HighsInt iPart = workPartition[iRow];
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = workEdWt[iRow];
          if (bestMerit[iPart] * myWeight < myInfeas) {
            bestMerit[iPart] = myInfeas / myWeight;
            bestIndex[iPart] = iRow;
          }
        }
      }
    }
    HighsInt count = 0;
    for (HighsInt i = 0; i < chLimit; i++) {
      if (bestIndex[i] != -1) {
        chIndex[count++] = bestIndex[i];
      }
    }
    *chCount = count;
  } else {
    // SPARSE mode
    if (workCount == 0) {
      analysis->simplexTimerStop(ChuzrDualClock);
      return;
    }

    HighsInt randomStart = ekk_instance_.random_.integer(workCount);
    vector<double> bestMerit(chLimit, 0);
    vector<HighsInt> bestIndex(chLimit, -1);
    for (HighsInt section = 0; section < 2; section++) {
      const HighsInt start = (section == 0) ? randomStart : 0;
      const HighsInt end = (section == 0) ? workCount : randomStart;
      for (HighsInt i = start; i < end; i++) {
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > kHighsZero) {
          HighsInt iPart = workPartition[iRow];
          const double myInfeas = work_infeasibility[iRow];
          const double myWeight = workEdWt[iRow];
          if (bestMerit[iPart] * myWeight < myInfeas) {
            bestMerit[iPart] = myInfeas / myWeight;
            bestIndex[iPart] = iRow;
          }
        }
      }
    }
    HighsInt count = 0;
    for (HighsInt i = 0; i < chLimit; i++) {
      if (bestIndex[i] != -1) {
        chIndex[count++] = bestIndex[i];
      }
    }
    *chCount = count;
  }

  analysis->simplexTimerStop(ChuzrDualClock);
}

void HEkkDualRHS::updatePrimal(HVector* column, double theta) {
  analysis->simplexTimerStart(UpdatePrimalClock);

  const HighsInt numRow = ekk_instance_.simplex_lp_.numRow_;
  const HighsInt columnCount = column->count;
  const HighsInt* variable_index = &column->index[0];
  const double* columnArray = &column->array[0];

  const double* baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  const double* baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  const double Tp = ekk_instance_.options_.primal_feasibility_tolerance;
  double* baseValue = &ekk_instance_.simplex_info_.baseValue_[0];

  bool updatePrimal_inDense = columnCount < 0 || columnCount > 0.4 * numRow;

  if (updatePrimal_inDense) {
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      baseValue[iRow] -= theta * columnArray[iRow];
      const double value = baseValue[iRow];
      const double less = baseLower[iRow] - value;
      const double more = value - baseUpper[iRow];
      double infeas = less > Tp ? less : (more > Tp ? more : 0);
      //    work_infeasibility[iRow] = infeas * infeas;
      if (ekk_instance_.simplex_info_.store_squared_primal_infeasibility)
        work_infeasibility[iRow] = infeas * infeas;
      else
        work_infeasibility[iRow] = fabs(infeas);
    }
  } else {
    for (HighsInt i = 0; i < columnCount; i++) {
      HighsInt iRow = variable_index[i];
      baseValue[iRow] -= theta * columnArray[iRow];
      const double value = baseValue[iRow];
      const double less = baseLower[iRow] - value;
      const double more = value - baseUpper[iRow];
      double infeas = less > Tp ? less : (more > Tp ? more : 0);
      if (ekk_instance_.simplex_info_.store_squared_primal_infeasibility)
        work_infeasibility[iRow] = infeas * infeas;
      else
        work_infeasibility[iRow] = fabs(infeas);
    }
  }

  analysis->simplexTimerStop(UpdatePrimalClock);
}

// Update the DSE weights
void HEkkDualRHS::updateWeightDualSteepestEdge(
    HVector* column, const double new_pivotal_edge_weight, double Kai,
    double* dseArray) {
  analysis->simplexTimerStart(DseUpdateWeightClock);

  const HighsInt numRow = ekk_instance_.simplex_lp_.numRow_;
  const HighsInt columnCount = column->count;
  const HighsInt* variable_index = &column->index[0];
  const double* columnArray = &column->array[0];

  bool updateWeight_inDense = columnCount < 0 || columnCount > 0.4 * numRow;
  if (updateWeight_inDense) {
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      const double aa_iRow = columnArray[iRow];
      workEdWt[iRow] +=
          aa_iRow * (new_pivotal_edge_weight * aa_iRow + Kai * dseArray[iRow]);
      if (workEdWt[iRow] < min_dual_steepest_edge_weight)
        workEdWt[iRow] = min_dual_steepest_edge_weight;
    }
  } else {
    for (HighsInt i = 0; i < columnCount; i++) {
      const HighsInt iRow = variable_index[i];
      const double aa_iRow = columnArray[iRow];
      workEdWt[iRow] +=
          aa_iRow * (new_pivotal_edge_weight * aa_iRow + Kai * dseArray[iRow]);
      if (workEdWt[iRow] < min_dual_steepest_edge_weight)
        workEdWt[iRow] = min_dual_steepest_edge_weight;
    }
  }
  analysis->simplexTimerStop(DseUpdateWeightClock);
}
// Update the Devex weights
void HEkkDualRHS::updateWeightDevex(HVector* column,
                                    const double new_pivotal_edge_weight) {
  analysis->simplexTimerStart(DevexUpdateWeightClock);

  const HighsInt numRow = ekk_instance_.simplex_lp_.numRow_;
  const HighsInt columnCount = column->count;
  const HighsInt* variable_index = &column->index[0];
  const double* columnArray = &column->array[0];

  bool updateWeight_inDense = columnCount < 0 || columnCount > 0.4 * numRow;
  if (updateWeight_inDense) {
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      double aa_iRow = columnArray[iRow];
      workEdWt[iRow] =
          max(workEdWt[iRow], new_pivotal_edge_weight * aa_iRow * aa_iRow);
    }
  } else {
    for (HighsInt i = 0; i < columnCount; i++) {
      HighsInt iRow = variable_index[i];
      double aa_iRow = columnArray[iRow];
      workEdWt[iRow] =
          max(workEdWt[iRow], new_pivotal_edge_weight * aa_iRow * aa_iRow);
    }
  }
  analysis->simplexTimerStop(DevexUpdateWeightClock);
}

void HEkkDualRHS::updatePivots(HighsInt iRow, double value) {
  // Update the primal value for the row (iRow) where the basis change
  // has occurred, and set the corresponding squared primal
  // infeasibility value in work_infeasibility
  //
  const double* baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  const double* baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  const double Tp = ekk_instance_.options_.primal_feasibility_tolerance;
  double* baseValue = &ekk_instance_.simplex_info_.baseValue_[0];
  baseValue[iRow] = value;
  double pivotInfeas = 0;
  if (baseValue[iRow] < baseLower[iRow] - Tp)
    pivotInfeas = baseValue[iRow] - baseLower[iRow];
  if (baseValue[iRow] > baseUpper[iRow] + Tp)
    pivotInfeas = baseValue[iRow] - baseUpper[iRow];
  // work_infeasibility[iRow] = pivotInfeas * pivotInfeas;
  if (ekk_instance_.simplex_info_.store_squared_primal_infeasibility)
    work_infeasibility[iRow] = pivotInfeas * pivotInfeas;
  else
    work_infeasibility[iRow] = fabs(pivotInfeas);
}

void HEkkDualRHS::updateInfeasList(HVector* column) {
  const HighsInt columnCount = column->count;
  const HighsInt* variable_index = &column->index[0];

  // DENSE mode: disabled
  if (workCount < 0) return;

  analysis->simplexTimerStart(UpdatePrimalClock);

  if (workCutoff <= 0) {
    // The regular sparse way
    for (HighsInt i = 0; i < columnCount; i++) {
      HighsInt iRow = variable_index[i];
      if (workMark[iRow] == 0) {
        if (work_infeasibility[iRow]) {
          workIndex[workCount++] = iRow;
          workMark[iRow] = 1;
        }
      }
    }
  } else {
    // The hyper sparse way
    for (HighsInt i = 0; i < columnCount; i++) {
      HighsInt iRow = variable_index[i];
      if (workMark[iRow] == 0) {
        if (work_infeasibility[iRow] > workEdWt[iRow] * workCutoff) {
          workIndex[workCount++] = iRow;
          workMark[iRow] = 1;
        }
      }
    }
  }

  analysis->simplexTimerStop(UpdatePrimalClock);
}

void HEkkDualRHS::createArrayOfPrimalInfeasibilities() {
  HighsInt numRow = ekk_instance_.simplex_lp_.numRow_;
  const double* baseValue = &ekk_instance_.simplex_info_.baseValue_[0];
  const double* baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  const double* baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  const double Tp = ekk_instance_.options_.primal_feasibility_tolerance;
  for (HighsInt i = 0; i < numRow; i++) {
    const double value = baseValue[i];
    const double less = baseLower[i] - value;
    const double more = value - baseUpper[i];
    double infeas = less > Tp ? less : (more > Tp ? more : 0);
    //    work_infeasibility[i] = infeas * infeas;
    if (ekk_instance_.simplex_info_.store_squared_primal_infeasibility)
      work_infeasibility[i] = infeas * infeas;
    else
      work_infeasibility[i] = fabs(infeas);
  }
}

void HEkkDualRHS::createInfeasList(double columnDensity) {
  HighsInt numRow = ekk_instance_.simplex_lp_.numRow_;
  double* dwork = &workEdWtFull[0];

  // 1. Build the full list
  fill_n(&workMark[0], numRow, 0);
  workCount = 0;
  workCutoff = 0;
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    if (work_infeasibility[iRow]) {
      workMark[iRow] = 1;
      workIndex[workCount++] = iRow;
    }
  }

  // 2. See if it worth to try to go sparse
  //    (Many candidates, really sparse RHS)
  if (workCount > max(numRow * 0.01, 500.0) && columnDensity < 0.05) {
    HighsInt icutoff = max(workCount * 0.001, 500.0);
    double maxMerit = 0;
    for (HighsInt iRow = 0, iPut = 0; iRow < numRow; iRow++)
      if (workMark[iRow]) {
        double myMerit = work_infeasibility[iRow] / workEdWt[iRow];
        if (maxMerit < myMerit) maxMerit = myMerit;
        dwork[iPut++] = -myMerit;
      }
    nth_element(dwork, dwork + icutoff, dwork + workCount);
    double cutMerit = -dwork[icutoff];
    workCutoff = min(maxMerit * 0.99999, cutMerit * 1.00001);

    // Create again
    fill_n(&workMark[0], numRow, 0);
    workCount = 0;
    for (HighsInt iRow = 0; iRow < numRow; iRow++) {
      if (work_infeasibility[iRow] >= workEdWt[iRow] * workCutoff) {
        workIndex[workCount++] = iRow;
        workMark[iRow] = 1;
      }
    }

    // Reduce by drop smaller
    if (workCount > icutoff * 1.5) {
      // Firstly take up "icutoff" number of elements
      HighsInt fullCount = workCount;
      workCount = icutoff;
      for (HighsInt i = icutoff; i < fullCount; i++) {
        HighsInt iRow = workIndex[i];
        if (work_infeasibility[iRow] > workEdWt[iRow] * cutMerit) {
          workIndex[workCount++] = iRow;
        } else {
          workMark[iRow] = 0;
        }
      }
    }
  }

  // 3. If there are still too many candidates: disable them
  if (workCount > 0.2 * numRow) {
    workCount = -numRow;
    workCutoff = 0;
  }
}
