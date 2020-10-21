/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkPrimal.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HEkkPrimal.h"

#include "simplex/HEkkDebug.h"

//#include <cassert>
//#include <cstdio>
//#include <iostream>

//#include "simplex/HSimplexDebug.h"
#include "simplex/SimplexTimer.h"
//#include "util/HighsUtils.h"

using std::runtime_error;

HighsStatus HEkkPrimal::solve() {
  HighsOptions& options = ekk_instance_.options_;
  HighsLp& simplex_lp = ekk_instance_.simplex_lp_;
  HighsLogMessage(
      options.logfile, HighsMessageType::INFO,
      "HEkkPrimal::solve called for LP with %d columns, %d rows and %d entries",
      simplex_lp.numCol_, simplex_lp.numRow_,
      simplex_lp.Astart_[simplex_lp.numCol_]);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = ekk_instance_.simplex_lp_.numRow_ > 0;
  if (!positive_num_row) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HEkkPrimal::solve called for LP with non-positive (%d) "
                    "number of constraints",
                    ekk_instance_.simplex_lp_.numRow_);
    assert(positive_num_row);
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }
  if (ekk_instance_.bailoutOnTimeIterations())
    return ekk_instance_.returnFromSolve(HighsStatus::Warning);

  // Set up bound perturbation as cost perturbation in HDual
  if (!use_bound_perturbation)
    HighsLogMessage(options.logfile, HighsMessageType::INFO,
                    "HEkkPrimal::solve not using bound perturbation");

  if (!simplex_lp_status.has_invert) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HEkkPrimal::solve called without INVERT");
    assert(simplex_lp_status.has_fresh_invert);
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }

  if (use_bound_perturbation) {
    ekk_instance_.computePrimal();
    ekk_instance_.computeSimplexPrimalInfeasible();
  }
  solvePhase =
      ekk_instance_.simplex_info_.num_primal_infeasibilities > 0 ? 1 : 2;

  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         use_bound_perturbation) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);

  //  if (ekkDebugSimplex(ekk_instance_, algorithm, solvePhase) ==
  //  HighsDebugStatus::LOGICAL_ERROR) return
  //  ekk_instance_.returnFromSolve(HighsStatus::Error);

  // The major solving loop
  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  // ToDo Move to simplex and adapt so it's OK for primal and dual
  //  iterationAnalysisInitialise();

  //  while (solvePhase) {
  //    int it0 = ekk_instance_.iteration_count_;
  //    switch (solvePhase) {
  //      case 1:
  //        analysis->simplexTimerStart(SimplexPrimalPhase1Clock);
  //        solvePhase1();
  //        analysis->simplexTimerStop(SimplexPrimalPhase1Clock);
  //        simplex_info.primal_phase1_iteration_count +=
  //        (ekk_instance_.iteration_count_ - it0);
  //	break;
  //    case 2:
  //      analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
  //      solvePhase2();
  //      analysis->simplexTimerStop(SimplexPrimalPhase2Clock);
  //      simplex_info.primal_phase2_iteration_count +=
  //      (ekk_instance_.iteration_count_ - it0); break;
  //    case 4:
  //      break;
  //    default:
  //      solvePhase = 0;
  //      break;
  //    }
  //    // Jump for primal
  //    if (solvePhase == 4) break;
  //    // Possibly bail out
  //    if (SolveBailout) break;
  //  }

  solvePhase = 2;
  assert(ekk_instance_.solve_bailout_ == false);
  //  HighsSimplexAnalysis& analysis = ekk_instance_.simplex_analysis_;
  if (solvePhase == 2) {
    int it0 = ekk_instance_.iteration_count_;

    analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
    solvePhase2();
    analysis->simplexTimerStop(SimplexPrimalPhase2Clock);

    simplex_info.primal_phase2_iteration_count +=
        (ekk_instance_.iteration_count_ - it0);
    if (ekk_instance_.bailoutReturn()) return HighsStatus::Warning;
  }
  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         use_bound_perturbation) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  return HighsStatus::OK;
}

void HEkkPrimal::solvePhase2() {
  //  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 2;
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailoutReturn()) return;
  // Setup update limits
  //  simplex_info.update_limit =
  //      min(100 + num_row / 100,
  //          1000);  // TODO: Consider allowing the dual limit to be used
  //  simplex_info.update_count = 0;

  HighsPrintMessage(ekk_instance_.options_.output,
                    ekk_instance_.options_.message_level, ML_DETAILED,
                    "primal-phase2-start\n");
  // Main solving structure
  for (;;) {
    analysis->simplexTimerStart(IteratePrimalRebuildClock);
    primalRebuild();
    analysis->simplexTimerStop(IteratePrimalRebuildClock);

    if (isPrimalPhase1) {
      for (;;) {
        // Primal phase 1 choose column
        phase1ChooseColumn();
        if (columnIn == -1) {
          invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
          break;
        }

        // Primal phase 1 choose row
        phase1ChooseRow();
        if (rowOut == -1) {
          HighsLogMessage(ekk_instance_.options_.logfile,
                          HighsMessageType::ERROR,
                          "Primal phase 1 choose row failed");
          exit(0);
        }

        // Primal phase 1 update
        phase1Update();
        if (invertHint) {
          break;
        }
        if (ekk_instance_.bailoutReturn()) return;
      }
      // Go to the next rebuild
      if (invertHint) {
        // Stop when the invert is new
        if (simplex_lp_status.has_fresh_rebuild) {
          break;
        }
        continue;
      }
    }
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
      if (ekk_instance_.bailoutReturn()) return;
      if (invertHint) {
        break;
      }
    }
    // If the data are fresh from rebuild() and no flips have occurred, break
    // out of the outer loop to see what's ocurred
    if (simplex_lp_status.has_fresh_rebuild && num_flip_since_rebuild == 0)
      break;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);

  if (columnIn == -1) {
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_DETAILED,
                      "primal-optimal\n");
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_DETAILED,
                      "problem-optimal\n");
    ekk_instance_.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_MINIMAL,
                      "primal-unbounded\n");
    ekk_instance_.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  ekk_instance_.computeDualObjectiveValue();
}

void HEkkPrimal::initialise() {
  analysis = &ekk_instance_.analysis_;

  num_col = ekk_instance_.simplex_lp_.numCol_;
  num_row = ekk_instance_.simplex_lp_.numRow_;
  num_tot = num_col + num_row;

  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance =
      ekk_instance_.scaled_solution_params_.primal_feasibility_tolerance;
  dual_feasibility_tolerance =
      ekk_instance_.scaled_solution_params_.dual_feasibility_tolerance;

  ekk_instance_.simplex_lp_status_.has_primal_objective_value = false;
  ekk_instance_.simplex_lp_status_.has_dual_objective_value = false;
  ekk_instance_.scaled_model_status_ = HighsModelStatus::NOTSET;
  ekk_instance_.solve_bailout_ = false;

  // Setup local vectors
  col_aq.setup(num_row);
  row_ep.setup(num_row);
  row_ap.setup(num_col);

  ph1SorterR.reserve(num_row);
  ph1SorterT.reserve(num_row);

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < num_tot; iCol++) {
    if (ekk_instance_.simplex_info_.workLower_[iCol] == -HIGHS_CONST_INF &&
        ekk_instance_.simplex_info_.workUpper_[iCol] == HIGHS_CONST_INF) {
      // Free column
      no_free_columns = false;
      break;
    }
  }
  if (!no_free_columns)
    HighsLogMessage(ekk_instance_.options_.logfile, HighsMessageType::INFO,
                    "HEkkPrimal:: LP has free columns");
}

void HEkkPrimal::primalRebuild() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;

  // Record whether the update objective value should be tested. If
  // the objective value is known, then the updated objective value
  // should be correct - once the correction due to recomputing the
  // dual values has been applied.
  //
  // Note that computePrimalObjectiveValue sets
  // has_primal_objective_value
  //
  // Have to do this before INVERT, as this permutes the indices of
  // basic variables, and baseValue only corresponds to the new
  // ordering once computePrimal has been called
  const bool check_updated_objective_value =
      simplex_lp_status.has_primal_objective_value;
  double previous_primal_objective_value;
  if (check_updated_objective_value) {
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase,
    //    "Before INVERT");
    previous_primal_objective_value =
        simplex_info.updated_primal_objective_value;
  } else {
    // Reset the knowledge of previous objective values
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, -1, "");
  }

  // Rebuild ekk_instance_.factor_ - only if we got updates
  int sv_invertHint = invertHint;
  invertHint = INVERT_HINT_NO;
  // Possibly Rebuild factor
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
    analysis->simplexTimerStart(InvertClock);
    int rank_deficiency = ekk_instance_.computeFactor();
    analysis->simplexTimerStop(InvertClock);
    if (rank_deficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    simplex_info.update_count = 0;
  }
  ekk_instance_.computeDual();
  ekk_instance_.computePrimal();
  ekk_instance_.computePrimalObjectiveValue();
  if (check_updated_objective_value) {
    // Apply the objective value correction due to computing primal
    // values from scratch.
    const double primal_objective_value_correction =
        simplex_info.primal_objective_value - previous_primal_objective_value;
    simplex_info.updated_primal_objective_value +=
        primal_objective_value_correction;
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm);
  }
  // Now that there's a new dual_objective_value, reset the updated
  // value
  simplex_info.updated_primal_objective_value =
      simplex_info.primal_objective_value;

  ekk_instance_.computeSimplexInfeasible();
  // Determine whether simplex_info.num_primal_infeasibilities and
  // simplex_info.num_dual_infeasibilities can be used
  //  ekk_instance_.copySimplexInfeasible();

  // Whether to switch to primal phase 1
  isPrimalPhase1 = 0;
  if (simplex_info.num_primal_infeasibilities > 0) {
    isPrimalPhase1 = 1;
    phase1ComputeDual();
  }

  reportRebuild(sv_invertHint);
  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HEkkPrimal::primalChooseColumn() {
  HighsRandom& random = ekk_instance_.random_;
  const int* jFlag = &ekk_instance_.simplex_basis_.nonbasicFlag_[0];
  const int* jMove = &ekk_instance_.simplex_basis_.nonbasicMove_[0];
  double* workDual = &ekk_instance_.simplex_info_.workDual_[0];
  const double* workLower = &ekk_instance_.simplex_info_.workLower_[0];
  const double* workUpper = &ekk_instance_.simplex_info_.workUpper_[0];
  const double dualTolerance =
      ekk_instance_.scaled_solution_params_.dual_feasibility_tolerance;

  analysis->simplexTimerStart(ChuzcPrimalClock);
  columnIn = -1;
  double bestInfeas = 0;
  if (no_free_columns) {
    const int numSection = 1;
    int startSection = random.integer() % numSection;
    int deltaCol = (num_tot + numSection - 1) / numSection;
    int fromCol = startSection * deltaCol;
    int toCol = min(fromCol + deltaCol, num_tot);
    int numPass = 1;
    //    printf("\nstartSection = %1d; deltaCol = %d\n", startSection,
    //    deltaCol);
    for (;;) {
      //      printf("CHUZC: %1d [%6d, %6d] %6d\n", numPass, fromCol, toCol,
      //      num_tot);
      for (int iCol = fromCol; iCol < toCol; iCol++) {
        // Then look at dual infeasible
        if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
          if (bestInfeas * devex_weight[iCol] < fabs(workDual[iCol])) {
            bestInfeas = fabs(workDual[iCol]) / devex_weight[iCol];
            columnIn = iCol;
          }
        }
      }
      if (columnIn >= 0 || numPass == numSection) {
        //	printf("Break from CHUZC after %d passes\n", numPass);
        break;
      }
      if (toCol == num_tot) {
        fromCol = 0;
        toCol = deltaCol;
      } else {
        fromCol = toCol;
        toCol = min(fromCol + deltaCol, num_tot);
      }
      numPass++;
    }
  } else {
    for (int iCol = 0; iCol < num_tot; iCol++) {
      if (jFlag[iCol] && fabs(workDual[iCol]) > dualTolerance) {
        // Always take free
        // TODO: if we found free,
        // Then deal with it in dual phase 1
        if (workLower[iCol] <= -HIGHS_CONST_INF &&
            workUpper[iCol] >= HIGHS_CONST_INF) {
          columnIn = iCol;
          break;
        }
        // Then look at dual infeasible
        if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
          if (bestInfeas * devex_weight[iCol] < fabs(workDual[iCol])) {
            bestInfeas = fabs(workDual[iCol]) / devex_weight[iCol];
            columnIn = iCol;
          }
        }
      }
    }
  }
  analysis->simplexTimerStop(ChuzcPrimalClock);
  columnIn = -1;
}

void HEkkPrimal::primalChooseRow() {
  const double* baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  const double* baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  double* baseValue = &ekk_instance_.simplex_info_.baseValue_[0];
  const double primalTolerance =
      ekk_instance_.scaled_solution_params_.primal_feasibility_tolerance;

  // Compute pivot column
  analysis->simplexTimerStart(FtranClock);
  col_aq.clear();
  col_aq.packFlag = true;
  ekk_instance_.matrix_.collect_aj(col_aq, columnIn, 1);
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq,
                                    analysis->col_aq_density);
#endif
  ekk_instance_.factor_.ftran(col_aq, analysis->col_aq_density,
                              analysis->pointer_serial_factor_clocks);
  analysis->simplexTimerStop(FtranClock);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
#endif

  const double local_col_aq_density = (double)col_aq.count / num_row;
  analysis->updateOperationResultDensity(local_col_aq_density,
                                         analysis->col_aq_density);

  const bool check_dual = false;
  if (check_dual) {
    const double* workCost = &ekk_instance_.simplex_info_.workCost_[0];
    const double* workDual = &ekk_instance_.simplex_info_.workDual_[0];
    const int* basicIndex = &ekk_instance_.simplex_basis_.basicIndex_[0];
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

  analysis->simplexTimerStart(Chuzr1Clock);
  // Initialize
  rowOut = -1;

  // Choose row pass 1
  double alphaTol =
      ekk_instance_.simplex_info_.update_count < 10
          ? 1e-9
          : ekk_instance_.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  const int* jMove = &ekk_instance_.simplex_basis_.nonbasicMove_[0];
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
  analysis->simplexTimerStop(Chuzr1Clock);

  analysis->simplexTimerStart(Chuzr2Clock);
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
  analysis->simplexTimerStop(Chuzr2Clock);
}

void HEkkPrimal::primalUpdate() {
  int* jMove = &ekk_instance_.simplex_basis_.nonbasicMove_[0];
  double* workDual = &ekk_instance_.simplex_info_.workDual_[0];
  const double* workLower = &ekk_instance_.simplex_info_.workLower_[0];
  const double* workUpper = &ekk_instance_.simplex_info_.workUpper_[0];
  const double* baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  const double* baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  double* workValue = &ekk_instance_.simplex_info_.workValue_[0];
  double* baseValue = &ekk_instance_.simplex_info_.baseValue_[0];
  const double primalTolerance =
      ekk_instance_.scaled_solution_params_.primal_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;

  // Compute thetaPrimal
  int moveIn = jMove[columnIn];
  //  int
  columnOut = ekk_instance_.simplex_basis_.basicIndex_[rowOut];
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

  analysis->simplexTimerStart(UpdatePrimalClock);
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    baseValue[index] -= thetaPrimal * col_aq.array[index];
  }
  analysis->simplexTimerStop(UpdatePrimalClock);

  simplex_info.updated_primal_objective_value +=
      workDual[columnIn] * thetaPrimal;

  // Why is the detailed primal infeasibility information needed?
  ekk_instance_.computeSimplexPrimalInfeasible();
  // copySimplexPrimalInfeasible();

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
  ekk_instance_.updatePivots(columnIn, rowOut, sourceOut);

  baseValue[rowOut] = valueIn;

  analysis->simplexTimerStart(CollectPrIfsClock);
  // Check for any possible infeasible
  for (int iRow = 0; iRow < num_row; iRow++) {
    if (baseValue[iRow] < baseLower[iRow] - primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    } else if (baseValue[iRow] > baseUpper[iRow] + primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    }
  }
  analysis->simplexTimerStop(CollectPrIfsClock);

  // 2. Now we can update the dual

  analysis->simplexTimerStart(BtranClock);
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
  ekk_instance_.factor_.btran(row_ep, analysis->row_ep_density,
                              analysis->pointer_serial_factor_clocks);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep);
#endif
  analysis->simplexTimerStop(BtranClock);
  //
  // PRICE
  //
  ekk_instance_.computeTableauRowFromPiP(row_ep, row_ap);
  /*
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) {
  analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep,
analysis->row_ap_density); analysis->num_row_price++;
  }
#endif
  analysis->simplexTimerStart(PriceClock);
  ekk_instance_.matrix_.priceByRowSparseResult(row_ap, row_ep);
  analysis->simplexTimerStop(PriceClock);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
  analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ep);
#endif

  const double local_row_ep_density = (double)row_ep.count / num_row;
  analysis->updateOperationResultDensity(local_row_ep_density,
analysis->row_ep_density);
  */
  analysis->simplexTimerStart(UpdateDualClock);
  //  double
  thetaDual = workDual[columnIn] / alpha;
  for (int i = 0; i < row_ap.count; i++) {
    int iCol = row_ap.index[i];
    workDual[iCol] -= thetaDual * row_ap.array[iCol];
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iGet = row_ep.index[i];
    int iCol = iGet + num_col;
    workDual[iCol] -= thetaDual * row_ep.array[iGet];
  }
  analysis->simplexTimerStop(UpdateDualClock);

  /* Update the devex weight */
  devexUpdate();

  // After dual update in primal simplex the dual objective value is not known
  ekk_instance_.simplex_lp_status_.has_dual_objective_value = false;

  // updateVerify for primal
  numericalTrouble = 0;
  /*
  double aCol = fabs(alpha);
  double alphaRow;
  if (columnIn < ekk_instance_.simplex_lp_.numCol_) {
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
  //  if (numericalTrouble > 1e-7 && ekk_instance_.simplex_info_.update_count >
  0) invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
  */
  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  // Update ekk_instance_.factor_ basis
  ekk_instance_.updateFactor(&col_aq, &row_ep, &rowOut, &invertHint);
  ekk_instance_.updateMatrix(columnIn, columnOut);
  if (simplex_info.update_count >= simplex_info.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }
  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  ekk_instance_.iteration_count_++;

  /* Reset the devex when there are too many errors */
  if (num_bad_devex_weight > 3) {
    devexReset();
  }

  // Report on the iteration
  iterationAnalysis();
}

void HEkkPrimal::phase1ComputeDual() {
  /* Alias to problem size, tolerance and work arrays */
  const int nRow = num_row;
  const int nCol = num_col;
  const double dFeasTol = dual_feasibility_tolerance;
  const double* baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  const double* baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  const double* baseValue = &ekk_instance_.simplex_info_.baseValue_[0];

  analysis->simplexTimerStart(BtranClock);
  /* Setup artificial cost and compute pi with BTran */
  HVector buffer;
  buffer.setup(nRow);
  buffer.clear();
  for (int iRow = 0; iRow < nRow; iRow++) {
    buffer.index[iRow] = iRow;
    if (baseValue[iRow] < baseLower[iRow] - dFeasTol) {
      buffer.array[iRow] = -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
      buffer.array[iRow] = 1.0;
    } else {
      buffer.array[iRow] = 0.0;
    }
  }
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, buffer,
                                    analysis->row_ep_density);
#endif
  ekk_instance_.factor_.btran(buffer, 1,
                              analysis->pointer_serial_factor_clocks);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, buffer);
#endif
  analysis->simplexTimerStop(BtranClock);

  analysis->simplexTimerStart(PriceClock);
  /* The compute the phase 1 reduced cost for all variables by SpMV */
  HVector bufferLong;
  bufferLong.setup(nCol);
  bufferLong.clear();
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) {
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP, buffer,
                                    0.0);
    analysis->num_col_price++;
  }
#endif
  ekk_instance_.matrix_.priceByColumn(bufferLong, buffer);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP, row_ap);
#endif
  analysis->simplexTimerStop(PriceClock);

  const int* nbFlag = &ekk_instance_.simplex_basis_.nonbasicFlag_[0];
  double* workDual = &ekk_instance_.simplex_info_.workDual_[0];
  for (int iSeq = 0; iSeq < nCol + nRow; iSeq++) {
    workDual[iSeq] = 0.0;
  }
  for (int iSeq = 0; iSeq < nCol; iSeq++) {
    if (nbFlag[iSeq]) workDual[iSeq] = -bufferLong.array[iSeq];
  }
  for (int iRow = 0, iSeq = nCol; iRow < nRow; iRow++, iSeq++) {
    if (nbFlag[iSeq]) workDual[iSeq] = -buffer.array[iRow];
  }

  /* Recompute number of dual infeasible variables with the phase 1 cost */
  ekk_instance_.computeSimplexDualInfeasible();
  // Determine whether simplex_info.num_dual_infeasibilities can be used
  //  copySimplexDualInfeasible();
}

void HEkkPrimal::phase1ChooseColumn() {
  const int nSeq = num_col + num_row;
  const int* nbMove = &ekk_instance_.simplex_basis_.nonbasicMove_[0];
  const double* workDual = &ekk_instance_.simplex_info_.workDual_[0];
  const double dDualTol =
      ekk_instance_.scaled_solution_params_.dual_feasibility_tolerance;
  analysis->simplexTimerStart(ChuzcPrimalClock);
  double dBestScore = 0;
  columnIn = -1;
  for (int iSeq = 0; iSeq < nSeq; iSeq++) {
    double dMyDual = nbMove[iSeq] * workDual[iSeq];
    double dMyScore = dMyDual / devex_weight[iSeq];
    if (dMyDual < -dDualTol && dMyScore < dBestScore) {
      dBestScore = dMyScore;
      columnIn = iSeq;
    }
  }
  analysis->simplexTimerStop(ChuzcPrimalClock);
  columnIn = -1;
}

void HEkkPrimal::phase1ChooseRow() {
  /* Alias to work arrays */
  const double dFeasTol =
      ekk_instance_.scaled_solution_params_.primal_feasibility_tolerance;
  const double* baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  const double* baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  const double* baseValue = &ekk_instance_.simplex_info_.baseValue_[0];

  /* Compute the transformed pivot column and update its density */
  analysis->simplexTimerStart(FtranClock);
  col_aq.clear();
  col_aq.packFlag = true;
  ekk_instance_.matrix_.collect_aj(col_aq, columnIn, 1);
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq,
                                    analysis->col_aq_density);
#endif
  ekk_instance_.factor_.ftran(col_aq, analysis->col_aq_density,
                              analysis->pointer_serial_factor_clocks);
  analysis->simplexTimerStop(FtranClock);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
#endif

  const double local_col_aq_density = (double)col_aq.count / num_row;
  analysis->updateOperationResultDensity(local_col_aq_density,
                                         analysis->col_aq_density);

  /* Compute the reduced cost for the pivot column and compare it with the kept
   * value */
  double dCompDual = 0.0;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    if (baseValue[iRow] < baseLower[iRow] - dFeasTol) {
      dCompDual -= col_aq.array[iRow] * -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
      dCompDual -= col_aq.array[iRow] * +1.0;
    }
  }
  if (fabs(ekk_instance_.simplex_info_.workDual_[columnIn] - dCompDual) >
      (fabs(dCompDual) + 1.0) * 1e-9) {
    printf("==> Phase 1 reduced cost. Updated %g, Computed %g\n",
           ekk_instance_.simplex_info_.workDual_[columnIn], dCompDual);
  }

  analysis->simplexTimerStart(Chuzr1Clock);
  /* Collect phase 1 theta lists */
  int nRow = num_row;
  const int iMoveIn = ekk_instance_.simplex_basis_.nonbasicMove_[columnIn];
  const double dPivotTol =
      ekk_instance_.simplex_info_.update_count < 10
          ? 1e-9
          : ekk_instance_.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  ph1SorterR.clear();
  ph1SorterT.clear();
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    double dAlpha = col_aq.array[iRow] * iMoveIn;

    /* When the basic variable x[i] decrease */
    if (dAlpha > +dPivotTol) {
      /* Whether it can become feasible by going below its upper bound */
      if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
        double dFeasTheta =
            (baseValue[iRow] - baseUpper[iRow] - dFeasTol) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow));
      }
      /* Whether it can become infeasible (again) by going below its lower bound
       */
      if (baseValue[iRow] > baseLower[iRow] - dFeasTol &&
          baseLower[iRow] > -HIGHS_CONST_INF) {
        double dRelaxTheta =
            (baseValue[iRow] - baseLower[iRow] + dFeasTol) / dAlpha;
        double dTightTheta = (baseValue[iRow] - baseLower[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow - nRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow - nRow));
      }
    }

    /* When the basic variable x[i] increase */
    if (dAlpha < -dPivotTol) {
      /* Whether it can become feasible by going above its lower bound */
      if (baseValue[iRow] < baseLower[iRow] - dFeasTol) {
        double dFeasTheta =
            (baseValue[iRow] - baseLower[iRow] + dFeasTol) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow - nRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow - nRow));
      }

      /* Whether it can become infeasible (again) by going above its upper bound
       */
      if (baseValue[iRow] < baseUpper[iRow] + dFeasTol &&
          baseUpper[iRow] < +HIGHS_CONST_INF) {
        double dRelaxTheta =
            (baseValue[iRow] - baseUpper[iRow] - dFeasTol) / dAlpha;
        double dTightTheta = (baseValue[iRow] - baseUpper[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow));
      }
    }
  }

  analysis->simplexTimerStop(Chuzr1Clock);
  /* When there is no candidates at all, we can leave it here */
  if (ph1SorterR.empty()) {
    rowOut = -1;
    columnOut = -1;
    return;
  }

  /*
   * Now sort the relaxed theta to find the final break point.
   * TODO: Consider partial sort.
   *       Or heapify [O(n)] and then pop k points [kO(log(n))].
   */
  analysis->simplexTimerStart(Chuzr2Clock);
  std::sort(ph1SorterR.begin(), ph1SorterR.end());
  double dMaxTheta = ph1SorterR.at(0).first;
  double dGradient = fabs(ekk_instance_.simplex_info_.workDual_[columnIn]);
  for (unsigned int i = 0; i < ph1SorterR.size(); i++) {
    double dMyTheta = ph1SorterR.at(i).first;
    int index = ph1SorterR.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    dGradient -= fabs(col_aq.array[iRow]);
    /* Stop when the gradient start to decrease */
    if (dGradient <= 0) {
      break;
    }
    dMaxTheta = dMyTheta;
  }

  /* Find out the biggest possible alpha for pivot */
  std::sort(ph1SorterT.begin(), ph1SorterT.end());
  double dMaxAlpha = 0.0;
  unsigned int iLast = ph1SorterT.size();
  for (unsigned int i = 0; i < ph1SorterT.size(); i++) {
    double dMyTheta = ph1SorterT.at(i).first;
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    /* Stop when the theta is too large */
    if (dMyTheta > dMaxTheta) {
      iLast = i;
      break;
    }
    /* Update the maximal possible alpha */
    if (dMaxAlpha < dAbsAlpha) {
      dMaxAlpha = dAbsAlpha;
    }
  }

  /* Finally choose a pivot with good enough alpha, backwardly */
  rowOut = -1;
  columnOut = -1;
  phase1OutBnd = 0;
  for (int i = iLast - 1; i >= 0; i--) {
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    if (dAbsAlpha > dMaxAlpha * 0.1) {
      rowOut = iRow;
      phase1OutBnd = index > 0 ? 1 : -1;
      break;
    }
  }
  if (rowOut != -1) {
    columnOut = ekk_instance_.simplex_basis_.basicIndex_[rowOut];
  }
  analysis->simplexTimerStop(Chuzr2Clock);
}

void HEkkPrimal::phase1Update() {}

void HEkkPrimal::devexReset() {
  devex_weight.assign(num_tot, 1.0);
  devex_index.assign(num_tot, 0);
  for (int iVar = 0; iVar < num_tot; iVar++) {
    const int nonbasicFlag = ekk_instance_.simplex_basis_.nonbasicFlag_[iVar];
    devex_index[iVar] = nonbasicFlag * nonbasicFlag;
  }
  num_devex_iterations = 0;
  num_bad_devex_weight = 0;
}

void HEkkPrimal::devexUpdate() {}

void HEkkPrimal::iterationAnalysisData() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  analysis->simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
  analysis->edge_weight_mode = DualEdgeWeightMode::DEVEX;
  analysis->solve_phase = solvePhase;
  analysis->simplex_iteration_count = ekk_instance_.iteration_count_;
  analysis->devex_iteration_count = num_devex_iterations;
  analysis->pivotal_row_index = rowOut;
  analysis->leaving_variable = columnOut;
  analysis->entering_variable = columnIn;
  analysis->invert_hint = invertHint;
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
      simplex_info.num_primal_infeasibilities;
  analysis->num_dual_infeasibilities = simplex_info.num_dual_infeasibilities;
  analysis->sum_primal_infeasibilities =
      simplex_info.sum_primal_infeasibilities;
  analysis->sum_dual_infeasibilities = simplex_info.sum_dual_infeasibilities;
#ifdef HiGHSDEV
  analysis->basis_condition = simplex_info.invert_condition;
#endif
  if ((analysis->edge_weight_mode == DualEdgeWeightMode::DEVEX) &&
      (num_devex_iterations == 0))
    analysis->num_devex_framework++;
}

void HEkkPrimal::iterationAnalysis() {
  iterationAnalysisData();
  analysis->iterationReport();
#ifdef HiGHSDEV
  analysis->iterationRecord();
#endif
}

void HEkkPrimal::reportRebuild(const int rebuild_invert_hint) {
  analysis->simplexTimerStart(ReportRebuildClock);
  iterationAnalysisData();
  analysis->invert_hint = rebuild_invert_hint;
  analysis->invertReport();
  analysis->simplexTimerStop(ReportRebuildClock);
}
