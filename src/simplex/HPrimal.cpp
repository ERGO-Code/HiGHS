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

#include <cassert>
#include <cstdio>
#include <iostream>

using std::runtime_error;

void HPrimal::solvePhase2(HighsModelObject *ptr_highs_model_object) {
  // Copy size
  highs_model_object = ptr_highs_model_object; // Pointer to highs_model_object: defined in HPrimal.h
  model = &highs_model_object->hmodel_[0];
  //  model->basis_ = &highs_model_object->basis_;
  numCol = model->lp_scaled_->numCol_;
  numRow = model->lp_scaled_->numRow_;
  numTot = model->lp_scaled_->numCol_ + model->lp_scaled_->numRow_;

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

#ifdef HiGHSDEV
  // Initialise rebuild count and time
  totalRebuildTime = 0;
  totalRebuilds = 0;
#endif

  // Setup other buffers

  model->util_reportMessage("primal-start");

  double lc_totalTime = model->totalTime + model->timer.getTime();

  // Main solving structure
  for (;;) {
    primalRebuild();

    for (;;) {
      primalChooseColumn();
      if (columnIn == -1) {
        invertHint = invertHint_possiblyOptimal;
        break;
      }
      primalChooseRow();
      if (rowOut == -1) {
        invertHint = invertHint_possiblyPrimalUnbounded;
        break;
      }
      primalUpdate();
      if (invertHint) {
        break;
      }
      // printf("HPrimal::solve_phase2: Iter = %d; Objective = %g\n",
      // model->numberIteration, model->dualObjectiveValue);
      if (model->dualObjectiveValue > model->dblOption[DBLOPT_OBJ_UB]) {
#ifdef SCIP_DEV
        printf("HPrimal::solve_phase2: Objective = %g > %g = dblOption[DBLOPT_OBJ_UB]\n",
	       model->dualObjectiveValue, model->dblOption[DBLOPT_OBJ_UB]);
#endif
        model->problemStatus = LP_Status_ObjUB;
        break;
      }
    }

    lc_totalTime = model->totalTime + model->timer.getTime();
    //	printf("Primal Ph2: lc_totalTime = %5.2f\n", lc_totalTime);
    if (lc_totalTime > TimeLimitValue) {
      model->problemStatus = LP_Status_OutOfTime;
      break;
    }
    if (model->problemStatus == LP_Status_ObjUB) break;
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (countUpdate == 0) break;
    if (model->mlFg_haveFreshRebuild) break;
  }

  if (model->problemStatus == LP_Status_OutOfTime ||
      model->problemStatus == LP_Status_ObjUB)
    return;

  if (columnIn == -1) {
    model->util_reportMessage("primal-optimal");
    model->util_reportMessage("problem-optimal");
    model->setProblemStatus(LP_Status_Optimal);
  } else {
    model->util_reportMessage("primal-unbounded");
    model->setProblemStatus(LP_Status_Unbounded);
  }
}

void HPrimal::primalRebuild() {
  model->recordPivots(-1, -1, 0);  // Indicate REINVERT
#ifdef HiGHSDEV
  double tt0 = 0;
  if (anRebuildTime) tt0 = model->timer.getTime();
#endif
  // Rebuild model->factor - only if we got updates
  int sv_invertHint = invertHint;
  invertHint = invertHint_no;  // Was 0
  // Possibly Rebuild model->factor
  bool reInvert = model->countUpdate > 0;
  if (!model->InvertIfRowOutNeg) {
    // Don't reinvert if rowOut is negative [equivalently, if sv_invertHint ==
    // invertHint_possiblyOptimal]
    if (sv_invertHint == invertHint_possiblyOptimal) {
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
  model->computeDualObjectiveValue();
  model->util_reportNumberIterationObjectiveValue(sv_invertHint);

#ifdef HiGHSDEV
  if (anRebuildTime) {
    double rebuildTime = model->timer.getTime() - tt0;
    totalRebuilds++;
    totalRebuildTime += rebuildTime;
    printf(
        "Primal     rebuild %d (%1d) on iteration %9d: Rebuild time = %g; "
        "Total rebuild time %g\n",
        totalRebuilds, sv_invertHint, model->numberIteration, rebuildTime,
        totalRebuildTime);
  }
#endif
  // Data are fresh from rebuild
  model->mlFg_haveFreshRebuild = 1;
}

void HPrimal::primalChooseColumn() {
  columnIn = -1;
  double bestInfeas = 0;
  const int *jFlag = &highs_model_object->basis_.nonbasicFlag_[0];
  const int *jMove = &highs_model_object->basis_.nonbasicMove_[0];
  double *workDual = &highs_model_object->simplex_.workDual_[0];
  const double *workLower = &highs_model_object->simplex_.workLower_[0];
  const double *workUpper = &highs_model_object->simplex_.workUpper_[0];
  const double dualTolerance = model->dblOption[DBLOPT_DUAL_TOL];

  const int numTot = model->lp_scaled_->numCol_ + model->lp_scaled_->numRow_;
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
  const double *baseLower = &highs_model_object->simplex_.baseLower_[0];
  const double *baseUpper = &highs_model_object->simplex_.baseUpper_[0];
  double *baseValue = &highs_model_object->simplex_.baseValue_[0];
  const double primalTolerance = model->dblOption[DBLOPT_PRIMAL_TOL];

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
  const int *jMove = &highs_model_object->basis_.nonbasicMove_[0];
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
  int *jMove = &highs_model_object->basis_.nonbasicMove_[0];
  double *workDual = &highs_model_object->simplex_.workDual_[0];
  const double *workLower = &highs_model_object->simplex_.workLower_[0];
  const double *workUpper = &highs_model_object->simplex_.workUpper_[0];
  const double *baseLower = &highs_model_object->simplex_.baseLower_[0];
  const double *baseUpper = &highs_model_object->simplex_.baseUpper_[0];
  double *workValue = &highs_model_object->simplex_.workValue_[0];
  double *baseValue = &highs_model_object->simplex_.baseValue_[0];
  const double primalTolerance = model->dblOption[DBLOPT_PRIMAL_TOL];

  // Compute thetaPrimal
  int moveIn = jMove[columnIn];
  int columnOut = highs_model_object->basis_.basicIndex_[rowOut];
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
      invertHint = invertHint_primalInfeasibleInPrimalSimplex;  // Was 1
    } else if (baseValue[iRow] > baseUpper[iRow] + primalTolerance) {
      invertHint = invertHint_primalInfeasibleInPrimalSimplex;  // Was 1
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
    invertHint = invertHint_updateLimitReached;  // Was true;

  model->recordPivots(columnIn, columnOut, alpha);
}
