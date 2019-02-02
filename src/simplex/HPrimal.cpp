/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HPrimal.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HPrimal.h"
#include "HModel.h"
#include "HConst.h"
#include "HighsIO.h"
#include "HSimplex.h"
#include "SimplexTimer.h"

#include <cassert>
#include <cstdio>
#include <iostream>

using std::runtime_error;

void HPrimal::solvePhase2() {
  model = &workHMO.hmodel_[0]; // Pointer to model within workHMO: defined in HDual.h
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;

  numCol = model->solver_lp_->numCol_;
  numRow = model->solver_lp_->numRow_;
  numTot = model->solver_lp_->numCol_ + model->solver_lp_->numRow_;

#ifdef HiGHSDEV
  printf("************************************\n");
  printf("Performing primal simplex iterations\n");
  printf("************************************\n");
#endif
  // Setup update limits
  limitUpdate = min(100 + numRow / 100, 1000);
  countUpdate = 0;

  // Setup local vectors
  column.setup(numRow);
  row_ep.setup(numRow);
  row_ap.setup(numCol);
  columnDensity = 0;
  row_epDensity = 0;

  // Setup other buffers

  HighsPrintMessage(ML_DETAILED, "primal-start\n");

  HighsTimer &timer = *(model->timer_);

  // Main solving structure
  for (;;) {
    timer.start(simplex_info.clock_[IteratePrimalRebuildClock]);
    primalRebuild();
    timer.stop(simplex_info.clock_[IteratePrimalRebuildClock]);

    for (;;) {
      primalChooseColumn();
      if (columnIn == -1) {
        invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
        break;
      }
      primalChooseRow();
      if (rowOut == -1) {
        invertHint = INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED;
        break;
      }
      primalUpdate();
      if (invertHint) {
        break;
      }
      double current_dual_objective_value = simplex_info.updatedDualObjectiveValue;
      // printf("HPrimal::solve_phase2: Iter = %d; Objective = %g\n",
      // model->numberIteration, current_dual_objective_value);
      if (current_dual_objective_value > simplex_info.dual_objective_value_upper_bound) {
#ifdef SCIP_DEV
        printf("HPrimal::solve_phase2: %12g = Objective > ObjectiveUB\n",
	       current_dual_objective_value, simplex_info.dual_objective_value_upper_bound);
#endif
        simplex_info.solution_status = SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND;
        break;
      }
    }

    double currentRunHighsTime = timer.readRunHighsClock();
    if (currentRunHighsTime > simplex_info.highs_run_time_limit) {
      simplex_info.solution_status = SimplexSolutionStatus::OUT_OF_TIME;
      break;
    }
    if (simplex_info.solution_status == SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND) break;
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (countUpdate == 0) break;
    if (model->mlFg_haveFreshRebuild) break;
  }

  if (simplex_info.solution_status == SimplexSolutionStatus::OUT_OF_TIME ||
      simplex_info.solution_status == SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND)
    return;

  if (columnIn == -1) {
    HighsPrintMessage(ML_DETAILED, "primal-optimal\n");
    HighsPrintMessage(ML_DETAILED, "problem-optimal\n");
    simplex_info.solution_status = SimplexSolutionStatus::OPTIMAL;
  } else {
    HighsPrintMessage(ML_MINIMAL, "primal-unbounded\n");
    simplex_info.solution_status = SimplexSolutionStatus::UNBOUNDED;
  }
}

void HPrimal::primalRebuild() {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsTimer &timer = workHMO.timer_;
  // Move this to Simplex class once it's created
  //  simplex_method.record_pivots(-1, -1, 0);  // Indicate REINVERT

  // Rebuild model->factor - only if we got updates
  int sv_invertHint = invertHint;
  invertHint = INVERT_HINT_NO;
  // Possibly Rebuild model->factor
  bool reInvert = model->countUpdate > 0;
  if (!model->InvertIfRowOutNeg) {
    // Don't reinvert if columnIn is negative [equivalently, if sv_invertHint ==
    // INVERT_HINT_POSSIBLY_OPTIMAL]
    if (sv_invertHint == INVERT_HINT_POSSIBLY_OPTIMAL) {
      assert(columnIn == -1);
      reInvert = false;
    }
  }
  if (reInvert) {
    int rankDeficiency = model->computeFactor();
    if (rankDeficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    countUpdate = 0;
  }
  model->computeDual();
  model->computePrimal();
  simplex_method_.computeDualObjectiveValue(workHMO);
  model->util_reportNumberIterationObjectiveValue(sv_invertHint);

#ifdef HiGHSDEV
  if (simplex_info.analyseRebuildTime) {
    int iClock = simplex_info.clock_[IteratePrimalRebuildClock];
    int totalRebuilds = timer.clockNumCall[iClock];
    double totalRebuildTime = timer.read(iClock);
    printf(
        "Primal     rebuild %d (%1d) on iteration %9d: Total rebuild time %g\n",
        totalRebuilds, sv_invertHint, model->numberIteration, totalRebuildTime);
  }
#endif
  // Data are fresh from rebuild
  model->mlFg_haveFreshRebuild = 1;
}

void HPrimal::primalChooseColumn() {
  columnIn = -1;
  double bestInfeas = 0;
  const int *jFlag = &workHMO.basis_.nonbasicFlag_[0];
  const int *jMove = &workHMO.basis_.nonbasicMove_[0];
  double *workDual = &workHMO.simplex_info_.workDual_[0];
  const double *workLower = &workHMO.simplex_info_.workLower_[0];
  const double *workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double dualTolerance = workHMO.simplex_info_.dual_feasibility_tolerance;

  const int numTot = model->solver_lp_->numCol_ + model->solver_lp_->numRow_;
  for (int iCol = 0; iCol < numTot; iCol++) {
    if (jFlag[iCol] && fabs(workDual[iCol]) > dualTolerance) {
      // Always take free
      // TODO: if we found free,
      // Then deal with it in dual phase 1
      if (workLower[iCol] == -HIGHS_CONST_INF &&
          workUpper[iCol] == HIGHS_CONST_INF) {
        columnIn = iCol;
        break;
      }
      // Then look at dual infeasible
      if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
        if (bestInfeas < fabs(workDual[iCol])) {
          bestInfeas = fabs(workDual[iCol]);
          columnIn = iCol;
        }
      }
    }
  }
}

void HPrimal::primalChooseRow() {
  const double *baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double *baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double *baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance = workHMO.simplex_info_.primal_feasibility_tolerance;

  // Compute pivot column
  column.clear();
  column.packFlag = true;
  model->matrix_->collect_aj(column, columnIn, 1);
  model->factor_->ftran(column, columnDensity);
  columnDensity = 0.95 * columnDensity + 0.05 * column.count / numRow;

  // Initialize
  rowOut = -1;

  // Choose column pass 1
  double alphaTol = countUpdate < 10 ? 1e-9 : countUpdate < 20 ? 1e-8 : 1e-7;
  const int *jMove = &workHMO.basis_.nonbasicMove_[0];
  int moveIn = jMove[columnIn];
  if (moveIn == 0) {
    // If there's still free in the N
    // We would report not-solved
    // Need to handle free
  }
  double relaxTheta = 1e100;
  for (int i = 0; i < column.count; i++) {
    int index = column.index[i];
    double alpha = column.array[index] * moveIn;
    if (alpha > alphaTol) {
      double relaxSpace = baseValue[index] - baseLower[index] + primalTolerance;
      if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    } else if (alpha < -alphaTol) {
      double relaxSpace = baseValue[index] - baseUpper[index] - primalTolerance;
      if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    }
  }

  // Choose column pass 2
  double bestAlpha = 0;
  for (int i = 0; i < column.count; i++) {
    int index = column.index[i];
    double alpha = column.array[index] * moveIn;
    if (alpha > alphaTol) {
      double tightSpace = baseValue[index] - baseLower[index];
      if (tightSpace < relaxTheta * alpha) {
        if (bestAlpha < alpha) {
          bestAlpha = alpha;
          rowOut = index;
        }
      }
    } else if (alpha < -alphaTol) {
      double tightSpace = baseValue[index] - baseUpper[index];
      if (tightSpace > relaxTheta * alpha) {
        if (bestAlpha < -alpha) {
          bestAlpha = -alpha;
          rowOut = index;
        }
      }
    }
  }
}

void HPrimal::primalUpdate() {
  int *jMove = &workHMO.basis_.nonbasicMove_[0];
  double *workDual = &workHMO.simplex_info_.workDual_[0];
  const double *workLower = &workHMO.simplex_info_.workLower_[0];
  const double *workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double *baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double *baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double *workValue = &workHMO.simplex_info_.workValue_[0];
  double *baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance = workHMO.simplex_info_.primal_feasibility_tolerance;

  // Compute thetaPrimal
  int moveIn = jMove[columnIn];
  int columnOut = workHMO.basis_.basicIndex_[rowOut];
  double alpha = column.array[rowOut];
  double thetaPrimal = 0;
  if (alpha * moveIn > 0) {
    // Lower bound
    thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alpha;
  } else {
    // Upper bound
    thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alpha;
  }

  // 1. Make sure it is inside bounds or just flip bound
  double lowerIn = workLower[columnIn];
  double upperIn = workUpper[columnIn];
  double valueIn = workValue[columnIn] + thetaPrimal;
  bool flipped = false;
  if (jMove[columnIn] == 1) {
    if (valueIn > upperIn + primalTolerance) {
      // Flip to upper
      workValue[columnIn] = upperIn;
      thetaPrimal = upperIn - lowerIn;
      flipped = true;
      jMove[columnIn] = -1;
    }
  } else if (jMove[columnIn] == -1) {
    if (valueIn < lowerIn - primalTolerance) {
      // Flip to lower
      workValue[columnIn] = lowerIn;
      thetaPrimal = lowerIn - upperIn;
      flipped = true;
      jMove[columnIn] = 1;
    }
  }

  for (int i = 0; i < column.count; i++) {
    int index = column.index[i];
    baseValue[index] -= thetaPrimal * column.array[index];
  }

  // If flipped, then no need touch the pivots
  if (flipped) {
    return;
  }

  // Pivot in
  int sourceOut = alpha * moveIn > 0 ? -1 : 1;
  model->updatePivots(columnIn, rowOut, sourceOut);

  baseValue[rowOut] = valueIn;

  // Check for any possible infeasible
  for (int iRow = 0; iRow < numRow; iRow++) {
    if (baseValue[iRow] < baseLower[iRow] - primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    } else if (baseValue[iRow] > baseUpper[iRow] + primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    }
  }

  // 2. Now we can update the dual
  row_ep.clear();
  row_ap.clear();
  row_ep.count = 1;
  row_ep.index[0] = rowOut;
  row_ep.array[rowOut] = 1;
  row_ep.packFlag = true;
  model->factor_->btran(row_ep, row_epDensity);
  model->matrix_->price_by_row(row_ap, row_ep);
  row_epDensity = 0.95 * row_epDensity + 0.05 * row_ep.count / numRow;

  double thetaDual = workDual[columnIn] / alpha;
  for (int i = 0; i < row_ap.count; i++) {
    int iCol = row_ap.index[i];
    workDual[iCol] -= thetaDual * row_ap.array[iCol];
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iGet = row_ep.index[i];
    int iCol = iGet + numCol;
    workDual[iCol] -= thetaDual * row_ep.array[iGet];
  }

  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  // Update model->factor basis
  model->updateFactor(&column, &row_ep, &rowOut, &invertHint);
  model->updateMatrix(columnIn, columnOut);
  if (++countUpdate >= limitUpdate)
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;  // Was true;

  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  model->numberIteration++;
}
