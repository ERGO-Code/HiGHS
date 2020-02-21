/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HPrimal.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HPrimal.h"

#include <cassert>
#include <cstdio>
#include <iostream>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsRandom.h"
#include "util/HighsUtils.h"

using std::runtime_error;

HighsStatus HPrimal::solve() {
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  workHMO.scaled_model_status_ = HighsModelStatus::NOTSET;
  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = workHMO.simplex_lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(workHMO.options_.logfile, HighsMessageType::ERROR,
                    "HPrimal::solve called for LP with non-positive (%d) "
                    "number of constraints",
                    workHMO.simplex_lp_.numRow_);
    return HighsStatus::Error;
  }

  HighsTimer& timer = workHMO.timer_;
  invertHint = INVERT_HINT_NO;

  // Setup aspects of the model data which are needed for solve() but better
  // left until now for efficiency reasons.
  // ToDo primal simplex version
  // setup_for_solve(workHMO);

  // Set SolveBailout to be true if control is to be returned immediately to
  // calling function
  // ToDo Move to Simplex
  //  SolveBailout = false;

  // Initialise working environment
  // Does LOTS, including initialisation of edge weights. Should only
  // be called if model dimension changes
  // ToDo primal simplex version
  // init();
  // initParallel();

  // ToDo primal simplex version
  // initialise_cost(workHMO, 1); //  model->initCost(1);
  assert(simplex_lp_status.has_fresh_invert);
  if (!simplex_lp_status.has_fresh_invert) {
    printf(
        "ERROR: Should enter with fresh INVERT - unless no_invert_on_optimal "
        "is set\n");
  }
  // Consider initialising edge weights - create Primal variants
  //
#ifdef HiGHSDEV
  //  printf("simplex_lp_status.has_dual_steepest_edge_weights 2 = %d;
  //  dual_edge_weight_mode = %d; DualEdgeWeightMode::STEEPEST_EDGE = %d\n",
  //	 simplex_lp_status.has_dual_steepest_edge_weights,
  // dual_edge_weight_mode, DualEdgeWeightMode::STEEPEST_EDGE);cout<<flush;
  //  printf("Edge weights known? %d\n",
  //  !simplex_lp_status.has_dual_steepest_edge_weights);cout<<flush;
#endif
  /*
  if (!simplex_lp_status.has_dual_steepest_edge_weights) {
    // Edge weights are not known
    // Set up edge weights according to dual_edge_weight_mode and
  initialise_dual_steepest_edge_weights
    // Using dual Devex edge weights
    // Zero the number of Devex frameworks used and set up the first one
    devex_index.assign(solver_num_tot, 0);
    initialiseDevexFramework();
    // Indicate that edge weights are known
    simplex_lp_status.has_dual_steepest_edge_weights = true;
  }
  */

  // ToDo Determine primal simplex phase from initial primal values
  //
  /*
  computePrimal(workHMO);
  compute_primal_infeasible_in_??(workHMO, &dualInfeasCount);
  solvePhase = ??InfeasCount > 0 ? 1 : 2;
  */
  solvePhase = 0;  // Frig to skip while (solvePhase) {*}

  // Check that the model is OK to solve:
  //
  // Level 0 just checks the flags
  //
  // Level 1 also checks that the basis is OK and that the necessary
  // data in work* is populated.
  //
  // Level 2 (will) checks things like the nonbasic duals and basic
  // primal values
  //
  // Level 3 (will) checks expensive things like the INVERT and
  // steepeest edge weights
  //
  // ToDo Write primal simplex equivalent
  /*
  bool ok = ok_to_solve(workHMO, 1, solvePhase);
  if (!ok) {printf("NOT OK TO SOLVE???\n"); cout << flush;}
  assert(ok);
  */
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "Before HPrimal major solving
  //  loop");
#endif
  // The major solving loop

  while (solvePhase) {
    /*
    int it0 = scaled_solution_params.simplex_iteration_count;
    switch (solvePhase) {
      case 1:
        timer.start(simplex_info.clock_[SimplexPrimalPhase1Clock]);
        solvePhase1();
        timer.stop(simplex_info.clock_[SimplexPrimalPhase1Clock]);
        simplex_info.primal_phase1_iteration_count +=
    (scaled_solution_params.simplex_iteration_count - it0); break; case 2:
        timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
        solvePhase2();
        timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);
        simplex_info.primal_phase2_iteration_count +=
    (scaled_solution_params.simplex_iteration_count - it0); break; case 4:
    break; default: solvePhase = 0; break;
    }
    // Jump for primal
    if (solvePhase == 4) break;
    // Possibly bail out
    if (SolveBailout) break;
    */
  }
  solvePhase = 2;
  if (workHMO.scaled_model_status_ != HighsModelStatus::REACHED_TIME_LIMIT) {
    if (solvePhase == 2) {
      int it0 = scaled_solution_params.simplex_iteration_count;

      timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
      solvePhase2();
      timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);

      simplex_info.primal_phase2_iteration_count +=
          (scaled_solution_params.simplex_iteration_count - it0);
    }
  }
  /*
  // ToDo Adapt ok_to_solve to be used by primal
  bool ok = ok_to_solve(workHMO, 1, solvePhase);// model->OKtoSolve(1,
  solvePhase); if (!ok) {printf("NOT OK After Solve???\n"); cout << flush;}
  assert(ok);
  */
  return HighsStatus::OK;
}

void HPrimal::solvePhase2() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;

  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 2;
  // Set up local copies of model dimensions
  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  analysis = &workHMO.simplex_analysis_;

  // Setup update limits
  simplex_info.update_limit =
      min(100 + solver_num_row / 100,
          1000);  // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);

#ifdef HiGHSDEV
  printf(
      "HPrimal::solvePhase2 - WARNING: Not setting analysis->col_aq_density = "
      "0\n");
  printf(
      "HPrimal::solvePhase2 - WARNING: Not setting analysis->row_ep_density = "
      "0\n");
#endif
  //  analysis->col_aq_density = 0;
  //  analysis->row_ep_density = 0;

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
#ifdef HiGHSDEV
  if (no_free_columns) {
    printf("Model has no free columns\n");
  } else {
    printf("Model has free columns\n");
  }
#endif

  // Setup other buffers

  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                    ML_DETAILED, "primal-phase2-start\n");
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
    }

    double currentRunHighsTime = timer.readRunHighsClock();
    if (currentRunHighsTime > workHMO.options_.time_limit) {
      workHMO.scaled_model_status_ = HighsModelStatus::REACHED_TIME_LIMIT;
      break;
    }
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) {
#ifdef HiGHSDEV
      if (num_flip_since_rebuild)
        printf("Consider doing a primal rebuild if flips have occurred\n");
#endif
      //      if (num_flip_since_rebuild == 0)
      break;
    }
  }

  if (workHMO.scaled_model_status_ == HighsModelStatus::REACHED_TIME_LIMIT) {
    return;
  }

  if (columnIn == -1) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_DETAILED, "primal-optimal\n");
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_DETAILED, "problem-optimal\n");
    workHMO.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level,
                      ML_MINIMAL, "primal-unbounded\n");
    workHMO.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  computeDualObjectiveValue(workHMO);
}

void HPrimal::primalRebuild() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
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
  computeDual(workHMO);
  timer.stop(simplex_info.clock_[ComputeDualClock]);

  timer.start(simplex_info.clock_[ComputePrimalClock]);
  computePrimal(workHMO);
  timer.stop(simplex_info.clock_[ComputePrimalClock]);

  // Primal objective section
  timer.start(simplex_info.clock_[ComputePrObjClock]);
  computePrimalObjectiveValue(workHMO);
  timer.stop(simplex_info.clock_[ComputePrObjClock]);

  double primal_objective_value = simplex_info.primal_objective_value;
#ifdef HiGHSDEV
  // Check the objective value maintained by updating against the
  // value when computed exactly - so long as there is a value to
  // check against
  if (simplex_lp_status.has_primal_objective_value) {
    double absPrimalObjectiveError = fabs(
        simplex_info.updated_primal_objective_value - primal_objective_value);
    double rlvPrimalObjectiveError =
        absPrimalObjectiveError / max(1.0, fabs(primal_objective_value));
    // TODO Investigate these Primal objective value errors
    if (rlvPrimalObjectiveError >= 1e-8) {
      HighsLogMessage(workHMO.options_.logfile, HighsMessageType::WARNING,
                      "Primal objective value error |rel| = %12g (%12g)",
                      absPrimalObjectiveError, rlvPrimalObjectiveError);
    }
  }
#endif
  simplex_info.updated_primal_objective_value = primal_objective_value;

  timer.start(simplex_info.clock_[ComputePrIfsClock]);
  computePrimalInfeasible(workHMO);
  timer.stop(simplex_info.clock_[ComputePrIfsClock]);

  timer.start(simplex_info.clock_[ComputeDuIfsClock]);
  computeDualInfeasible(workHMO);
  timer.stop(simplex_info.clock_[ComputeDuIfsClock]);

  timer.start(simplex_info.clock_[ReportRebuildClock]);
  reportRebuild(sv_invertHint);
  timer.stop(simplex_info.clock_[ReportRebuildClock]);

#ifdef HiGHSDEV
  if (simplex_info.analyse_rebuild_time) {
    int iClock = simplex_info.clock_[IteratePrimalRebuildClock];
    int totalRebuilds = timer.clock_num_call[iClock];
    double totalRebuildTime = timer.read(iClock);
    printf(
        "Primal     rebuild %d (%1d) on iteration %9d: Total rebuild time %g\n",
        totalRebuilds, sv_invertHint,
        workHMO.scaled_solution_params_.simplex_iteration_count,
        totalRebuildTime);
  }
#endif
  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HPrimal::primalChooseColumn() {
  HighsRandom& random = workHMO.random_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  const int* jFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  const int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double dualTolerance =
      workHMO.scaled_solution_params_.dual_feasibility_tolerance;

  timer.start(simplex_info.clock_[ChuzcPrimalClock]);
  columnIn = -1;
  double bestInfeas = 0;
  if (no_free_columns) {
    const int numSection = 1;
    int startSection = random.integer() % numSection;
    int deltaCol = (solver_num_tot + numSection - 1) / numSection;
    int fromCol = startSection * deltaCol;
    int toCol = min(fromCol + deltaCol, solver_num_tot);
    int numPass = 1;
    //    printf("\nstartSection = %1d; deltaCol = %d\n", startSection,
    //    deltaCol);
    for (;;) {
      //      printf("CHUZC: %1d [%6d, %6d] %6d\n", numPass, fromCol, toCol,
      //      solver_num_tot);
      for (int iCol = fromCol; iCol < toCol; iCol++) {
        // Then look at dual infeasible
        if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
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
        toCol = min(fromCol + deltaCol, solver_num_tot);
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
  timer.stop(simplex_info.clock_[ChuzcPrimalClock]);
}

void HPrimal::primalChooseRow() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance =
      workHMO.scaled_solution_params_.primal_feasibility_tolerance;

  // Compute pivot column
  timer.start(simplex_info.clock_[FtranClock]);
  col_aq.clear();
  col_aq.packFlag = true;
  workHMO.matrix_.collect_aj(col_aq, columnIn, 1);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq,
                                    analysis->col_aq_density);
#endif
  workHMO.factor_.ftran(col_aq, analysis->col_aq_density);
  timer.stop(simplex_info.clock_[FtranClock]);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
#endif

  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density,
                                         analysis->col_aq_density);

  const bool check_dual = false;
  if (check_dual) {
    const double* workCost = &workHMO.simplex_info_.workCost_[0];
    const double* workDual = &workHMO.simplex_info_.workDual_[0];
    const int* basicIndex = &workHMO.simplex_basis_.basicIndex_[0];
    double check_dual_value = workCost[columnIn];
    for (int i = 0; i < col_aq.count; i++) {
      int row = col_aq.index[i];
      int col = basicIndex[row];
      double value = col_aq.array[row];
      double cost = workCost[col];
      check_dual_value -= value * cost;
      //    printf("Entry %2d: [%2d, %12g] Cost = %12g; check_dual_value =
      //    %12g\n", i, row, value, cost, check_dual_value);
    }
    thetaDual = workDual[columnIn];
    double dual_error =
        fabs(check_dual_value - thetaDual) / max(1.0, fabs(thetaDual));
    if (dual_error > 1e-8)
      printf("Checking dual: updated = %12g; direct = %12g; error = %12g\n",
             thetaDual, check_dual_value, dual_error);
  }

  timer.start(simplex_info.clock_[Chuzr1Clock]);
  // Initialize
  rowOut = -1;

  // Choose row pass 1
  double alphaTol = workHMO.simplex_info_.update_count < 10
                        ? 1e-9
                        : workHMO.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  const int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  int moveIn = jMove[columnIn];
  if (moveIn == 0) {
    // If there's still free in the N
    // We would report not-solved
    // Need to handle free
  }
  double relaxTheta = 1e100;
  double relaxSpace;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    alpha = col_aq.array[index] * moveIn;
    if (alpha > alphaTol) {
      relaxSpace = baseValue[index] - baseLower[index] + primalTolerance;
      if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    } else if (alpha < -alphaTol) {
      relaxSpace = baseValue[index] - baseUpper[index] - primalTolerance;
      if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    }
  }
  timer.stop(simplex_info.clock_[Chuzr1Clock]);

  timer.start(simplex_info.clock_[Chuzr2Clock]);
  double bestAlpha = 0;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    alpha = col_aq.array[index] * moveIn;
    if (alpha > alphaTol) {
      // Positive pivotal column entry
      double tightSpace = baseValue[index] - baseLower[index];
      if (tightSpace < relaxTheta * alpha) {
        if (bestAlpha < alpha) {
          bestAlpha = alpha;
          rowOut = index;
        }
      }
    } else if (alpha < -alphaTol) {
      // Negative pivotal column entry
      double tightSpace = baseValue[index] - baseUpper[index];
      if (tightSpace > relaxTheta * alpha) {
        if (bestAlpha < -alpha) {
          bestAlpha = -alpha;
          rowOut = index;
        }
      }
    }
  }
  timer.stop(simplex_info.clock_[Chuzr2Clock]);
}

void HPrimal::primalUpdate() {
  HighsTimer& timer = workHMO.timer_;
  int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* workValue = &workHMO.simplex_info_.workValue_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance =
      workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;

  // Compute thetaPrimal
  int moveIn = jMove[columnIn];
  //  int
  columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  //  double
  alpha = col_aq.array[rowOut];
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
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    baseValue[index] -= thetaPrimal * col_aq.array[index];
  }
  timer.stop(simplex_info.clock_[UpdatePrimalClock]);

  simplex_info.updated_primal_objective_value +=
      workDual[columnIn] * thetaPrimal;

  computePrimalInfeasible(workHMO);

  // If flipped, then no need touch the pivots
  if (flipped) {
    rowOut = -1;
    numericalTrouble = 0;
    thetaDual = workDual[columnIn];
    iterationAnalysis();
    num_flip_since_rebuild++;
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
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep,
                                    analysis->row_ep_density);
#endif
  workHMO.factor_.btran(row_ep, analysis->row_ep_density);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep);
#endif
  timer.stop(simplex_info.clock_[BtranClock]);
  //
  // PRICE
  //
  computeTableauRowFromPiP(workHMO, row_ep, row_ap);
  /*
  #ifdef HiGHSDEV
    if (simplex_info.analyse_iterations) {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
  analysis->row_ap_density); analysis->num_row_price++;
    }
  #endif
    timer.start(simplex_info.clock_[PriceClock]);
    workHMO.matrix_.priceByRowSparseResult(row_ap, row_ep);
    timer.stop(simplex_info.clock_[PriceClock]);
  #ifdef HiGHSDEV
    if (simplex_info.analyse_iterations)
      analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep);
  #endif

    const double local_row_ep_density = (double)row_ep.count / solver_num_row;
    analysis->updateOperationResultDensity(local_row_ep_density,
  analysis->row_ep_density);
  */
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

  // After dual update in primal simplex the dual objective value is not known
  workHMO.simplex_lp_status_.has_dual_objective_value = false;

  // updateVerify for primal
  numericalTrouble = 0;
  /*
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
    printf("Numerical check: alphaCol = %12g, alphaRow = a%12g, aDiff = a%12g:
  measure = %12g\n", alpha, alphaRow, aDiff, numericalTrouble);
  // Reinvert if the relative difference is large enough, and updates have been
  performed
  //  if (numericalTrouble > 1e-7 && workHMO.simplex_info_.update_count > 0)
  invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
  */
  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  // Update workHMO.factor_ basis
  update_factor(workHMO, &col_aq, &row_ep, &rowOut, &invertHint);
  update_matrix(workHMO, columnIn, columnOut);
  if (simplex_info.update_count >= simplex_info.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }
  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  workHMO.scaled_solution_params_.simplex_iteration_count++;

  // Report on the iteration
  iterationAnalysis();
}

void HPrimal::iterationAnalysisData() {
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  analysis->simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
  analysis->edge_weight_mode = DualEdgeWeightMode::DANTZIG;
  analysis->solve_phase = solvePhase;
  analysis->simplex_iteration_count =
      scaled_solution_params.simplex_iteration_count;
  analysis->devex_iteration_count = 0;
  analysis->pivotal_row_index = rowOut;
  analysis->leaving_variable = columnOut;
  analysis->entering_variable = columnIn;
  analysis->invert_hint = invertHint;
  analysis->freelist_size = 0;
  analysis->reduced_rhs_value = 0;
  analysis->reduced_cost_value = 0;
  analysis->edge_weight = 0;
  analysis->primal_delta = 0;
  analysis->primal_step = thetaPrimal;
  analysis->dual_step = thetaDual;
  analysis->pivot_value_from_column = alpha;
  analysis->pivot_value_from_row = alpha;  // Row;
  analysis->numerical_trouble = numericalTrouble;
  analysis->objective_value = simplex_info.updated_primal_objective_value;
  analysis->num_primal_infeasibilities =
      scaled_solution_params.num_primal_infeasibilities;
  analysis->num_dual_infeasibilities =
      scaled_solution_params.num_dual_infeasibilities;
  analysis->sum_primal_infeasibilities =
      scaled_solution_params.sum_primal_infeasibilities;
  analysis->sum_dual_infeasibilities =
      scaled_solution_params.sum_dual_infeasibilities;
#ifdef HiGHSDEV
  analysis->basis_condition = simplex_info.invert_condition;
#endif
}

void HPrimal::iterationAnalysis() {
  // Possibly report on the iteration
  iterationAnalysisData();
  analysis->iterationReport();

#ifdef HiGHSDEV
  analysis->iterationRecord();
#endif
}

void HPrimal::reportRebuild(const int rebuild_invert_hint) {
  iterationAnalysisData();
  analysis->invert_hint = rebuild_invert_hint;
  analysis->invertReport();
}
