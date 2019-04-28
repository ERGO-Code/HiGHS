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
#include "simplex/HPrimal.h"
#include "lp_data/HConst.h"
#include "io/HighsIO.h"
#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsUtils.h"
#include "util/HighsRandom.h"

#include <cassert>
#include <cstdio>
#include <iostream>

using std::runtime_error;

void HPrimal::solvePhase2() {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer &timer = workHMO.timer_;

  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

#ifdef HiGHSDEV
  printf("************************************\n");
  printf("Performing primal simplex iterations\n");
  printf("************************************\n");
#endif
  // Setup update limits
  simplex_info.update_limit = min(100 + solver_num_row / 100, 1000); // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  column.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);
  columnDensity = 0;
  row_epDensity = 0;

  num_tabu_col = 0;
  tabu_col_p.resize(solver_num_tot);
  tabu_col.assign(solver_num_tot, 1);
  
  no_free_columns = true;
  for (int iCol = 0; iCol < solver_num_tot; iCol++) {
    if (highs_isInfinity(-workHMO.simplex_info_.workLower_[iCol])) {
      if (highs_isInfinity(workHMO.simplex_info_.workUpper_[iCol])) {
	// Free column
	no_free_columns = false;
	break;
      }
    }
  }
  if (no_free_columns) {
    printf("Model has no free columns\n");
  } else {
    printf("Model has free columns\n");
  }

  // Setup other buffers

  HighsPrintMessage(ML_DETAILED, "primal-start\n");

  // Main solving structure
  for (;;) {
    timer.start(simplex_info.clock_[IteratePrimalRebuildClock]);
    primalRebuild();
    timer.stop(simplex_info.clock_[IteratePrimalRebuildClock]);

    for (;;) {
      for (;;) {
	if (num_tabu_col == solver_num_tot) {
	  printf("All columns tabu!\n");
	  invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
	  break;
	}
	primalChooseColumn();
	if (columnIn == -1) {
	  printf("Possibly optimal: num_tabu_col = %d\n", num_tabu_col);
	  invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
	  break;
	}
	primalChooseRow();
	if (alpha) break;
	tabu_col_p[num_tabu_col] = columnIn;
	tabu_col[columnIn] = 0;
	num_tabu_col++;
      }
      if (invertHint) break;
      if (num_tabu_col) {
	for (int p = 0; p < num_tabu_col; p++) {
	  tabu_col[tabu_col_p[p]] = 1;
	}
	//	printf("Removed %2d tabu columns\n", num_tabu_col);
	num_tabu_col = 0;
      }
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
      // simplex_info.iteration_count, current_dual_objective_value);
      if (current_dual_objective_value > simplex_info.dual_objective_value_upper_bound) {
#ifdef SCIP_DEV
        printf("HPrimal::solve_phase2: %12g = Objective > ObjectiveUB\n",
	       current_dual_objective_value, simplex_info.dual_objective_value_upper_bound);
#endif
        simplex_lp_status.solution_status = SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND;
        break;
      }
    }

    double currentRunHighsTime = timer.readRunHighsClock();
    if (currentRunHighsTime > simplex_info.highs_run_time_limit) {
      simplex_lp_status.solution_status = SimplexSolutionStatus::OUT_OF_TIME;
      break;
    }
    if (simplex_lp_status.solution_status == SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND) break;
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) break;
  }

  if (simplex_lp_status.solution_status == SimplexSolutionStatus::OUT_OF_TIME ||
      simplex_lp_status.solution_status == SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND)
    return;

  if (columnIn == -1) {
    HighsPrintMessage(ML_DETAILED, "primal-optimal\n");
    HighsPrintMessage(ML_DETAILED, "problem-optimal\n");
    simplex_lp_status.solution_status = SimplexSolutionStatus::OPTIMAL;
#ifdef HiGHSDEV
    util_analyse_lp_solution(workHMO);
#endif    
  } else {
    HighsPrintMessage(ML_MINIMAL, "primal-unbounded\n");
    simplex_lp_status.solution_status = SimplexSolutionStatus::UNBOUNDED;
  }
  compute_dual_objective_value(workHMO);
}

void HPrimal::primalRebuild() {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer &timer = workHMO.timer_;
  // Move this to Simplex class once it's created
  //  simplex_method.record_pivots(-1, -1, 0);  // Indicate REINVERT

  // Rebuild workHMO.factor_ - only if we got updates
  int sv_invertHint = invertHint;
  invertHint = INVERT_HINT_NO;
  // Possibly Rebuild workHMO.factor_
  bool reInvert = simplex_info.update_count > 0;
  if (!invert_if_row_out_negative) {
    // Don't reinvert if columnIn is negative [equivalently, if sv_invertHint ==
    // INVERT_HINT_POSSIBLY_OPTIMAL]
    if (sv_invertHint == INVERT_HINT_POSSIBLY_OPTIMAL) {
      assert(columnIn == -1);
      reInvert = false;
    }
  }
  if (reInvert) {
    timer.start(simplex_info.clock_[InvertClock]);
    int rankDeficiency = compute_factor(workHMO);
    timer.stop(simplex_info.clock_[InvertClock]);
    if (rankDeficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    simplex_info.update_count = 0;
  }
  timer.start(simplex_info.clock_[ComputeDualClock]);
  compute_dual(workHMO);
  timer.stop(simplex_info.clock_[ComputeDualClock]);

  timer.start(simplex_info.clock_[ComputePrimalClock]);
  compute_primal(workHMO);
  timer.stop(simplex_info.clock_[ComputePrimalClock]);

  timer.start(simplex_info.clock_[CollectPrIfsClock]);
  int numPrimalInfeas = computePrimalInfeasible(workHMO);
  int numBinaryColumnValues = computeNumBinaryColumnValues(workHMO);
  timer.stop(simplex_info.clock_[CollectPrIfsClock]);

  if (num_tabu_col) {
    for (int p = 0; p < num_tabu_col; p++) {
      tabu_col[tabu_col_p[p]] = 1;
    }
    num_tabu_col = 0;
  }

  // Primal objective section
  bool checkPrimalObjectiveValue = simplex_lp_status.has_primal_objective_value;
  timer.start(simplex_info.clock_[ComputeProbjClock]);
  compute_primal_objective_value(workHMO);
  timer.stop(simplex_info.clock_[ComputeProbjClock]);
  report_iteration_count_primal_objective_value(workHMO, numBinaryColumnValues);//sv_invertHint);

  double primalObjectiveValue = simplex_info.primalObjectiveValue;
  if (checkPrimalObjectiveValue) {
    double absPrimalObjectiveError = fabs(simplex_info.updatedPrimalObjectiveValue - primalObjectiveValue);
    double rlvPrimalObjectiveError = absPrimalObjectiveError/max(1.0, fabs(primalObjectiveValue));
#ifdef HiGHSDEV
    // TODO Investigate these Primal objective value errors
    if (rlvPrimalObjectiveError >= 1e-8) {
      HighsLogMessage(HighsMessageType::WARNING, "Primal objective value error abs(rel) = %12g (%12g)",
			absPrimalObjectiveError, rlvPrimalObjectiveError);
    }
#endif
  }
  simplex_info.updatedPrimalObjectiveValue = primalObjectiveValue;

#ifdef HiGHSDEV
  if (simplex_info.analyseRebuildTime) {
    int iClock = simplex_info.clock_[IteratePrimalRebuildClock];
    int totalRebuilds = timer.clock_num_call[iClock];
    double totalRebuildTime = timer.read(iClock);
    printf(
        "Primal     rebuild %d (%1d) on iteration %9d: Total rebuild time %g\n",
        totalRebuilds, sv_invertHint, simplex_info.iteration_count, totalRebuildTime);
  }
#endif
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HPrimal::primalChooseColumn() {
  HighsRandom &random = workHMO.random_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsTimer &timer = workHMO.timer_;
  const int *jFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  const int *jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double *workDual = &workHMO.simplex_info_.workDual_[0];
  const double *workLower = &workHMO.simplex_info_.workLower_[0];
  const double *workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double dualTolerance = workHMO.simplex_info_.dual_feasibility_tolerance;

  timer.start(simplex_info.clock_[ChuzcClock]);
  columnIn = -1;
  double bestInfeas = 0;
  if (no_free_columns) {
    const int numSection = 1;
    int startSection = random.integer() % numSection;
    int deltaCol = (solver_num_tot+numSection-1)/numSection;
    int fromCol = startSection*deltaCol;
    int toCol = min(fromCol+deltaCol, solver_num_tot);
    int numPass = 1;
    //    printf("\nstartSection = %1d; deltaCol = %d\n", startSection, deltaCol);
    for (;;) {
      //      printf("CHUZC: %1d [%6d, %6d] %6d\n", numPass, fromCol, toCol, solver_num_tot);
      for (int iCol = fromCol; iCol < toCol; iCol++) {
	// Then look at dual infeasible
	if (tabu_col[iCol] * jMove[iCol] * workDual[iCol] < -dualTolerance) {
	  if (bestInfeas < fabs(workDual[iCol])) {
	    bestInfeas = fabs(workDual[iCol]);
	    columnIn = iCol;
	  }
	}
      }
      if (columnIn >= 0 || numPass == numSection) {
	//	printf("Break from CHUZC after %d passes\n", numPass);
	break;
      }
      if (toCol == solver_num_tot) {
	fromCol = 0;
	toCol = deltaCol;
      } else {
	fromCol = toCol;
	toCol = min(fromCol+deltaCol, solver_num_tot);
      }
      numPass++;
    }
  } else {
    for (int iCol = 0; iCol < solver_num_tot; iCol++) {
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
  timer.stop(simplex_info.clock_[ChuzcClock]);
}

void HPrimal::primalChooseRow() {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsTimer &timer = workHMO.timer_;
  const double *baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double *baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double *baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance = workHMO.simplex_info_.primal_feasibility_tolerance;
  const bool require_unit_pivot = false;
  const bool allow_extra_pass = false;
    
  // Compute pivot column
  timer.start(simplex_info.clock_[FtranClock]);
  column.clear();
  column.packFlag = true;
  workHMO.matrix_.collect_aj(column, columnIn, 1);
  workHMO.factor_.ftran(column, columnDensity);
  timer.stop(simplex_info.clock_[FtranClock]);
  columnDensity = 0.95 * columnDensity + 0.05 * column.count / solver_num_row;

  timer.start(simplex_info.clock_[Chuzr1Clock]);
  // Initialize
  rowOut = -1;
  int simplexIteration = simplex_info.iteration_count;

  bool report = false;
  // Choose row pass 1
  double alphaTol = workHMO.simplex_info_.update_count < 10 ? 1e-9 : workHMO.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  const int *jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  int moveIn = jMove[columnIn];
  if (moveIn == 0) {
    // If there's still free in the N
    // We would report not-solved
    // Need to handle free
  }
  int numPass = 1;
  for (;;) {
    report = numPass == 2;
    double relaxTheta = 1e100;
    double relaxSpace;
    for (int i = 0; i < column.count; i++) {
      int index = column.index[i];
      //    double
      alpha = column.array[index] * moveIn;
      bool report_alpha = false;
      if (alpha > alphaTol) {
	report_alpha = abs(alpha-1.0) < primalTolerance;
	relaxSpace = baseValue[index] - baseLower[index] + primalTolerance;
	if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
      } else if (alpha < -alphaTol) {
	report_alpha = abs(alpha+1.0) < primalTolerance;
	relaxSpace = baseValue[index] - baseUpper[index] - primalTolerance;
	if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
      }
      double ratio = relaxSpace / alpha;
      if (report && report_alpha && ratio < 10*relaxTheta)
	printf("Row %6d: entry %12g: relaxSpace = %12g; ratio = %12g; minRatio = %12g\n", index, alpha, relaxSpace, ratio, relaxTheta);
    }
    if (numPass == 1) timer.stop(simplex_info.clock_[Chuzr1Clock]);
    
    if (numPass == 1) timer.start(simplex_info.clock_[Chuzr2Clock]);
    // Choose row pass 2
    double bestAlpha = 0;
    int numCandidate = 0;
    int unitPivotRowOut = -1;
    for (int i = 0; i < column.count; i++) {
      int index = column.index[i];
      //    double
      alpha = column.array[index] * moveIn;
      if (alpha > alphaTol) {
	// Positive pivotal column entry
	double tightSpace = baseValue[index] - baseLower[index];
	if (tightSpace <= relaxTheta * alpha) {// NB JAJH Changed < to <= 
	  numCandidate++;
	  if (numPass == 2) {
	    printf("CHUZR: %3d has tightSpace = %12g < %12g = %12g * %12g (relaxTheta * alpha)\n",
		   numCandidate, tightSpace, relaxTheta * alpha, relaxTheta, alpha);
	  }
	  if (alpha == 1.0) unitPivotRowOut = index;
	  if (bestAlpha < alpha) {
	    bestAlpha = alpha;
	    rowOut = index;
	  }
	}
      } else if (alpha < -alphaTol) {
	// Negative pivotal column entry
	double tightSpace = baseValue[index] - baseUpper[index];
	if (tightSpace >= relaxTheta * alpha) {// NB JAJH Changed > to >= 
	  numCandidate++;
	  if (numPass == 2) {
	    printf("CHUZR: %3d has tightSpace = %12g > %12g = %12g * %12g (relaxTheta * alpha)\n",
		   numCandidate, tightSpace, relaxTheta * alpha, relaxTheta, alpha);
	  }
	  if (alpha == -1.0) unitPivotRowOut = index;
	  if (bestAlpha < -alpha) {
	    bestAlpha = -alpha;
	    rowOut = index;
	  }
	}
      }
    }
    if (unitPivotRowOut >= 0) {
      rowOut = unitPivotRowOut;
      alpha = column.array[rowOut];
    } else {
      if (require_unit_pivot) {
	rowOut = -1;
	alpha = 0;
      } else {
	alpha = column.array[rowOut];
      }
    }
    if (numPass == 1) timer.stop(simplex_info.clock_[Chuzr2Clock]);
    if (rowOut >= 0) break;
    if (numPass > 1) break;
    if (!allow_extra_pass) break;
    numPass++;
    printf("Iteration %6d: Cannot find unit pivot: relaxTheta = %12g\n", simplexIteration, relaxTheta);
  }
}

void HPrimal::primalUpdate() {
  HighsTimer &timer = workHMO.timer_;
  int *jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double *workDual = &workHMO.simplex_info_.workDual_[0];
  const double *workLower = &workHMO.simplex_info_.workLower_[0];
  const double *workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double *baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double *baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double *workValue = &workHMO.simplex_info_.workValue_[0];
  double *baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance = workHMO.simplex_info_.primal_feasibility_tolerance;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;

  // Compute thetaPrimal
  int moveIn = jMove[columnIn];
  //  int
  columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  //  double
  alpha = column.array[rowOut];
  //  double
  thetaPrimal = 0;
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

  timer.start(simplex_info.clock_[UpdatePrimalClock]);
  for (int i = 0; i < column.count; i++) {
    int index = column.index[i];
    baseValue[index] -= thetaPrimal * column.array[index];
  }
  timer.stop(simplex_info.clock_[UpdatePrimalClock]);

  simplex_info.updatedPrimalObjectiveValue += workDual[columnIn]*thetaPrimal;

  // If flipped, then no need touch the pivots
  if (flipped) {
    return;
  }

  // Pivot in
  int sourceOut = alpha * moveIn > 0 ? -1 : 1;
  timer.start(simplex_info.clock_[IteratePivotsClock]);
  update_pivots(workHMO, columnIn, rowOut, sourceOut);
  timer.stop(simplex_info.clock_[IteratePivotsClock]);

  baseValue[rowOut] = valueIn;

  timer.start(simplex_info.clock_[CollectPrIfsClock]);
  // Check for any possible infeasible
  for (int iRow = 0; iRow < solver_num_row; iRow++) {
    if (baseValue[iRow] < baseLower[iRow] - primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    } else if (baseValue[iRow] > baseUpper[iRow] + primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    }
  }
  timer.stop(simplex_info.clock_[CollectPrIfsClock]);

  // 2. Now we can update the dual
  
  timer.start(simplex_info.clock_[BtranClock]);
  row_ep.clear();
  row_ap.clear();
  row_ep.count = 1;
  row_ep.index[0] = rowOut;
  row_ep.array[rowOut] = 1;
  row_ep.packFlag = true;
#ifdef HiGHSDEV
  //  if (simplex_info.analyseSimplexIterations) iterateOpRecBf(AnIterOpTy_Btran, row_ep, row_epDensity);
#endif
  workHMO.factor_.btran(row_ep, row_epDensity);
#ifdef HiGHSDEV
  //  if (simplex_info.analyseSimplexIterations) iterateOpRecAf(AnIterOpTy_Btran, row_ep);
#endif
  timer.stop(simplex_info.clock_[BtranClock]);

  timer.start(simplex_info.clock_[PriceClock]);
  workHMO.matrix_.price_by_row(row_ap, row_ep);
  timer.stop(simplex_info.clock_[PriceClock]);
  row_epDensity = 0.95 * row_epDensity + 0.05 * row_ep.count / solver_num_row;

  timer.start(simplex_info.clock_[UpdateDualClock]);
  //  double
  thetaDual = workDual[columnIn] / alpha;
  for (int i = 0; i < row_ap.count; i++) {
    int iCol = row_ap.index[i];
    workDual[iCol] -= thetaDual * row_ap.array[iCol];
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iGet = row_ep.index[i];
    int iCol = iGet + solver_num_col;
    workDual[iCol] -= thetaDual * row_ep.array[iGet];
  }
  timer.stop(simplex_info.clock_[UpdateDualClock]);

  /*
  // updateVerify for primal
  double aCol = fabs(alpha);
  double alphaRow;
  if (columnIn < workHMO.simplex_lp_.numCol_) {
    alphaRow = row_ap.array[columnIn];
  } else {
    alphaRow = row_ep.array[rowOut];
  }
  double aRow = fabs(alphaRow);
  double aDiff = fabs(aCol - aRow);
  numericalTrouble = aDiff / min(aCol, aRow);
  if (numericalTrouble > 1e-7) 
    printf("Numerical check: alphaCol = %12g, alphaRow = a%12g, aDiff = a%12g: measure = %12g\n", alpha, alphaRow, aDiff, numericalTrouble);
  // Reinvert if the relative difference is large enough, and updates have been performed
  //  if (numericalTrouble > 1e-7 && workHMO.simplex_info_.update_count > 0) invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
  */
  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  // Update workHMO.factor_ basis
  update_factor(workHMO, &column, &row_ep, &rowOut, &invertHint);
  update_matrix(workHMO, columnIn, columnOut);
  if (simplex_info.update_count >= simplex_info.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }
  // Report on the iteration
  iterateRp();

  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  simplex_info.iteration_count++;

}

void HPrimal::iterateRp() {
  int numIter = workHMO.simplex_info_.iteration_count;
  bool header = numIter % 10 == 1;
  //  header = true;  // JAJH10/10
  if (header) iterateRpFull(header);
  iterateRpFull(false);
}

void HPrimal::iterateRpFull(bool header) {
  if (header) {
    iterateRpIterPh(ML_DETAILED, true);
    iterateRpPrObj(ML_DETAILED, true);
#ifdef HiGHSDEV
    iterateRpIterDa(ML_DETAILED, true);
    //    iterateRpDsty(ML_DETAILED, true);
    //    HighsPrintMessage(ML_DETAILED, " FreeLsZ");
#endif
    HighsPrintMessage(ML_DETAILED, "\n");
  } else {
    iterateRpIterPh(ML_DETAILED, false);
    iterateRpPrObj(ML_DETAILED, false);
#ifdef HiGHSDEV
    iterateRpIterDa(ML_DETAILED, false);
    //    iterateRpDsty(ML_DETAILED, false);
    //    HighsPrintMessage(ML_DETAILED, " %7d", dualRow.freeListSize);
#endif
    if (workHMO.options_.messageLevel & ML_DETAILED) {
      int numBinaryColumnValues = computeNumBinaryColumnValues(workHMO);
      HighsPrintMessage(ML_DETAILED, " %7d", numBinaryColumnValues);
    }
    HighsPrintMessage(ML_DETAILED, "\n");
  }
}

void HPrimal::iterateRpIterPh(int iterate_log_level, bool header) {
  int solvePhase=2;
  if (header) {
    HighsPrintMessage(iterate_log_level, " Iteration Ph");
  } else {
    int numIter = workHMO.simplex_info_.iteration_count;
    HighsPrintMessage(iterate_log_level, " %9d %2d", numIter, solvePhase);
  }
}

void HPrimal::iterateRpPrObj(int iterate_log_level, bool header) {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  if (header) {
    HighsPrintMessage(iterate_log_level, "  PrimalObjective    ");
  } else {
    HighsPrintMessage(iterate_log_level, " %20.10e", simplex_info.updatedPrimalObjectiveValue);
  }
}

void HPrimal::iterateRpIterDa(int iterate_log_level, bool header) {
  if (header) {
    HighsPrintMessage(iterate_log_level, " Inv       NumCk     EnC     LvR     LvC        ThDu        ThPr          Aa");
  } else {
    HighsPrintMessage(iterate_log_level, " %3d %11.4g %7d %7d %7d %11.4g %11.4g %11.4g", 
		      invertHint, numericalTrouble, columnIn, rowOut, columnOut, 
		      thetaDual, thetaPrimal, alpha);
  }
}

/*
void HPrimal::iterateRpDsty(int iterate_log_level, bool header) {
  bool rp_dual_steepest_edge = dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE;
  if (header) {
    HighsPrintMessage(iterate_log_level, "  Col R_Ep R_Ap");
    if (rp_dual_steepest_edge) {
      HighsPrintMessage(iterate_log_level, "  DSE");
    } else {
      HighsPrintMessage(iterate_log_level, "     ");
    }
  } else {
    int l10ColDse = intLog10(columnDensity);
    int l10REpDse = intLog10(row_epDensity);
    int l10RapDse = intLog10(row_apDensity);
    HighsPrintMessage(iterate_log_level, " %4d %4d %4d", l10ColDse, l10REpDse, l10RapDse);
    if (rp_dual_steepest_edge) {
      int l10DseDse = intLog10(rowdseDensity);
      HighsPrintMessage(iterate_log_level, " %4d", l10DseDse);
    } else {
      HighsPrintMessage(iterate_log_level, "     ");
    }
  }
}
int HPrimal::intLog10(double v) {
  int intLog10V = -99;
  if (v > 0) intLog10V = log(v) / log(10.0);
  return intLog10V;
}

*/
void HPrimal::iterateRpInvert(int i_v) {
#ifdef HiGHSDEV
  HighsPrintMessage(ML_MINIMAL, "Iter %10d:", workHMO.simplex_info_.iteration_count);
  //  iterateRpDsty(ML_MINIMAL, true);
  //  iterateRpDsty(ML_MINIMAL, false);
  iterateRpPrObj(ML_MINIMAL, false);
  HighsPrintMessage(ML_MINIMAL, " %2d\n", i_v);
#else
  report_iteration_count_primal_objective_value(workHMO, i_v);
#endif
}

