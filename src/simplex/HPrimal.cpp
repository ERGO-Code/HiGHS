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

#include <cassert>
#include <cstdio>
#include <iostream>

using std::runtime_error;

void HPrimal::solvePhase2() {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;

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

  // Setup other buffers

  HighsPrintMessage(ML_DETAILED, "primal-start\n");

  HighsTimer &timer = workHMO.timer_;

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
  } else {
    HighsPrintMessage(ML_MINIMAL, "primal-unbounded\n");
    simplex_lp_status.solution_status = SimplexSolutionStatus::UNBOUNDED;
  }
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
    int rankDeficiency = compute_factor(workHMO);
    if (rankDeficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    simplex_info.update_count = 0;
  }
  compute_dual(workHMO);
  compute_primal(workHMO);
  int numPrimalInfeas = computePrimalInfeasible(workHMO);

  compute_dual_objective_value(workHMO);
  report_iteration_count_dual_objective_value(workHMO, sv_invertHint);

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
  columnIn = -1;
  double bestInfeas = 0;
  const int *jFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  const int *jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double *workDual = &workHMO.simplex_info_.workDual_[0];
  const double *workLower = &workHMO.simplex_info_.workLower_[0];
  const double *workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double dualTolerance = workHMO.simplex_info_.dual_feasibility_tolerance;

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

void HPrimal::primalChooseRow() {
  const double *baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double *baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double *baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance = workHMO.simplex_info_.primal_feasibility_tolerance;

  // Compute pivot column
  column.clear();
  column.packFlag = true;
  workHMO.matrix_.collect_aj(column, columnIn, 1);
  workHMO.factor_.ftran(column, columnDensity);
  columnDensity = 0.95 * columnDensity + 0.05 * column.count / solver_num_row;

  // Initialize
  rowOut = -1;

  // Choose column pass 1
  double alphaTol = workHMO.simplex_info_.update_count < 10 ? 1e-9 : workHMO.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  const int *jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
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
  update_pivots(workHMO, columnIn, rowOut, sourceOut);

  baseValue[rowOut] = valueIn;

  // Check for any possible infeasible
  for (int iRow = 0; iRow < solver_num_row; iRow++) {
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
  workHMO.factor_.btran(row_ep, row_epDensity);
  workHMO.matrix_.price_by_row(row_ap, row_ep);
  row_epDensity = 0.95 * row_epDensity + 0.05 * row_ep.count / solver_num_row;

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
    iterateRpDuObj(ML_DETAILED, true);
#ifdef HiGHSDEV
    iterateRpIterDa(ML_DETAILED, true);
    //    iterateRpDsty(ML_DETAILED, true);
    //    HighsPrintMessage(ML_DETAILED, " FreeLsZ");
#endif
    HighsPrintMessage(ML_DETAILED, "\n");
  } else {
    iterateRpIterPh(ML_DETAILED, false);
    iterateRpDuObj(ML_DETAILED, false);
#ifdef HiGHSDEV
    iterateRpIterDa(ML_DETAILED, false);
    //    iterateRpDsty(ML_DETAILED, false);
    //    HighsPrintMessage(ML_DETAILED, " %7d", dualRow.freeListSize);
#endif
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

void HPrimal::iterateRpDuObj(int iterate_log_level, bool header) {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  if (header) {
    HighsPrintMessage(iterate_log_level, "    DualObjective    ");
  } else {
    HighsPrintMessage(iterate_log_level, " %20.10e", simplex_info.updatedDualObjectiveValue);
  }
}

void HPrimal::iterateRpIterDa(int iterate_log_level, bool header) {
  if (header) {
    HighsPrintMessage(iterate_log_level, " Inv       NumCk     EnC     LvR     LvC        DlPr        ThDu        ThPr          Aa");
  } else {
    HighsPrintMessage(iterate_log_level, " %3d %11.4g %7d %7d %7d %11.4g %11.4g %11.4g %11.4g", 
		      invertHint, numericalTrouble, columnIn, rowOut, columnOut, deltaPrimal,
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
  iterateRpDuObj(ML_MINIMAL, false);
  HighsPrintMessage(ML_MINIMAL, " %2d\n", i_v);
#else
  report_iteration_count_dual_objective_value(workHMO, i_v);
#endif
}

