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

  // Get the nonabsic free column set
  getNonbasicFreeColumnSet();

  if (use_bound_perturbation) {
    ekk_instance_.computePrimal();
    ekk_instance_.computeSimplexPrimalInfeasible();
  }
  int num_primal_infeasibilities =
      ekk_instance_.simplex_info_.num_primal_infeasibilities;
  solvePhase = num_primal_infeasibilities > 0 ? SOLVE_PHASE_1 : SOLVE_PHASE_2;

  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         use_bound_perturbation) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);

  // The major solving loop
  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  // ToDo Move to simplex and adapt so it's OK for primal and dual
  //  iterationAnalysisInitialise();

  while (solvePhase) {
    int it0 = ekk_instance_.iteration_count_;
    // When starting a new phase the (updated) primal objective function
    // value isn't known. Indicate this so that when the value
    // computed from scratch in rebuild() isn't checked against the the
    // updated value
    simplex_lp_status.has_primal_objective_value = false;
    if (solvePhase == SOLVE_PHASE_UNKNOWN) {
      // Reset the phase 2 bounds so that true number of dual
      // infeasibilities can be determined
      ekk_instance_.initialiseBound();
      // Determine the number of primal infeasibilities, and hence the solve
      // phase
      ekk_instance_.computeSimplexPrimalInfeasible();
      num_primal_infeasibilities =
          ekk_instance_.simplex_info_.num_primal_infeasibilities;
      solvePhase =
          num_primal_infeasibilities > 0 ? SOLVE_PHASE_1 : SOLVE_PHASE_2;
      /*
      if (simplex_info.backtracking_) {
        // Backtracking, so set the bounds and primal values
        ekk_instance_.initialiseBound(solvePhase);
        ekk_instance_.initialiseValueAndNonbasicMove();
        // Can now forget that we might have been backtracking
        simplex_info.backtracking_ = false;
      }
      */
    }
    assert(solvePhase == SOLVE_PHASE_1 || solvePhase == SOLVE_PHASE_2);
    if (solvePhase == SOLVE_PHASE_1) {
      // Phase 1
      analysis->simplexTimerStart(SimplexPrimalPhase1Clock);
      solvePhase1();
      analysis->simplexTimerStop(SimplexPrimalPhase1Clock);
      simplex_info.primal_phase1_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else if (solvePhase == SOLVE_PHASE_2) {
      // Phase 2
      analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
      solvePhase2();
      analysis->simplexTimerStop(SimplexPrimalPhase2Clock);
      simplex_info.primal_phase2_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else {
      // Should only be SOLVE_PHASE_1 or SOLVE_PHASE_2
      ekk_instance_.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
      return ekk_instance_.returnFromSolve(HighsStatus::Error);
    }
    // Return if bailing out from solve
    if (ekk_instance_.solve_bailout_)
      return ekk_instance_.returnFromSolve(HighsStatus::Warning);
    // Can have all possible cases of solvePhase
    assert(solvePhase >= SOLVE_PHASE_MIN && solvePhase <= SOLVE_PHASE_MAX);
    // Look for scenarios when the major solving loop ends
    if (solvePhase == SOLVE_PHASE_ERROR) {
      // Solver error so return HighsStatus::Error
      ekk_instance_.scaled_model_status_ = HighsModelStatus::SOLVE_ERROR;
      return ekk_instance_.returnFromSolve(HighsStatus::Error);
    }
    if (solvePhase == SOLVE_PHASE_EXIT) {
      // LP identified as not having an optimal solution
      assert(ekk_instance_.scaled_model_status_ ==
                 HighsModelStatus::PRIMAL_DUAL_INFEASIBLE ||
             ekk_instance_.scaled_model_status_ ==
                 HighsModelStatus::PRIMAL_INFEASIBLE ||
             ekk_instance_.scaled_model_status_ ==
                 HighsModelStatus::PRIMAL_UNBOUNDED);
      break;
    }
    if (solvePhase == SOLVE_PHASE_1 && ekk_instance_.scaled_model_status_ ==
                                           HighsModelStatus::DUAL_INFEASIBLE) {
      // Dual infeasibilities after phase 2 for a problem known to be dual
      // infeasible.
      break;
    }
    if (solvePhase == SOLVE_PHASE_CLEANUP) {
      // Primal infeasibilities after phase 2 for a problem not known
      // to be primal infeasible. Dual feasible with primal
      // infeasibilities so use dual simplex to clean up
      break;
    }
    // If solvePhase == SOLVE_PHASE_OPTIMAL == 0 then major solving
    // loop ends naturally since solvePhase is false
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  // Should only have these cases
  assert(solvePhase == SOLVE_PHASE_EXIT || solvePhase == SOLVE_PHASE_UNKNOWN ||
         solvePhase == SOLVE_PHASE_OPTIMAL || solvePhase == SOLVE_PHASE_1 ||
         solvePhase == SOLVE_PHASE_CLEANUP);
  if (solvePhase == SOLVE_PHASE_OPTIMAL)
    ekk_instance_.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         use_bound_perturbation) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  return ekk_instance_.returnFromSolve(HighsStatus::OK);
}

void HEkkPrimal::solvePhase1() {
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailoutReturn()) return;
  HighsPrintMessage(ekk_instance_.options_.output,
                    ekk_instance_.options_.message_level, ML_DETAILED,
                    "primal-phase1-start\n");
  // Main solving structure
  for (;;) {
    analysis->simplexTimerStart(IteratePrimalRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IteratePrimalRebuildClock);
    if (solvePhase == SOLVE_PHASE_ERROR) return;
    if (!isPrimalPhase1) {
      // No primal infeasibilities found in rebuild() so break and
      // return to phase 2
      solvePhase = SOLVE_PHASE_2;
      break;
    }
    if (ekk_instance_.bailoutOnTimeIterations()) return;

    for (;;) {
      if (debugPrimalSimplex("Before phase 1 iteration") ==
          HighsDebugStatus::LOGICAL_ERROR) {
        solvePhase = SOLVE_PHASE_ERROR;
        return;
      }
      // Primal phase 1 choose column
      chooseColumn();
      if (columnIn == -1) {
        invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
        break;
      }

      // Primal phase 1 choose row
      phase1ChooseRow();
      if (rowOut == -1) {
        HighsLogMessage(ekk_instance_.options_.logfile, HighsMessageType::ERROR,
                        "Primal phase 1 choose row failed");
        solvePhase = SOLVE_PHASE_ERROR;
        return;
      }

      // Primal phase 1 update
      phase1Update();
      if (invertHint) {
        break;
      }
      if (ekk_instance_.bailoutOnTimeIterations()) return;
    }
    // Go to the next rebuild
    if (invertHint) {
      // Stop when the invert is new
      if (simplex_lp_status.has_fresh_rebuild) {
        break;
      }
      continue;
    }
    // If the data are fresh from rebuild() and no flips have occurred, break
    // out of the outer loop to see what's ocurred
    if (simplex_lp_status.has_fresh_rebuild && num_flip_since_rebuild == 0)
      break;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  if (debugPrimalSimplex("End of solvePhase1") ==
      HighsDebugStatus::LOGICAL_ERROR) {
    solvePhase = SOLVE_PHASE_ERROR;
    return;
  }
  if (ekk_instance_.simplex_info_.num_primal_infeasibilities == 0) {
    // Primal feasible so switch to phase 2
    solvePhase = SOLVE_PHASE_2;
  }
}

void HEkkPrimal::solvePhase2() {
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailoutReturn()) return;
  HighsPrintMessage(ekk_instance_.options_.output,
                    ekk_instance_.options_.message_level, ML_DETAILED,
                    "primal-phase2-start\n");
  // Main solving structure
  for (;;) {
    analysis->simplexTimerStart(IteratePrimalRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IteratePrimalRebuildClock);
    if (solvePhase == SOLVE_PHASE_ERROR) return;
    if (ekk_instance_.bailoutOnTimeIterations()) return;

    if (isPrimalPhase1) {
      // Primal infeasibilities found in rebuild() Should be
      // shifted but, for now, break and return to phase 1
      solvePhase = SOLVE_PHASE_1;
      break;
    }

    for (;;) {
      if (debugPrimalSimplex("Before phase 2 iteration") ==
          HighsDebugStatus::LOGICAL_ERROR) {
        solvePhase = SOLVE_PHASE_ERROR;
        return;
      }
      chooseColumn();
      if (columnIn == -1) {
        invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
        break;
      }
      chooseRow();
      phase2Update();
      if (ekk_instance_.bailoutOnTimeIterations()) return;
      if (invertHint) {
        break;
      }
    }
    // If the data are fresh from rebuild() and no flips have occurred, break
    // out of the outer loop to see what's ocurred
    if (simplex_lp_status.has_fresh_rebuild && num_flip_since_rebuild == 0)
      break;
  }
  if (debugPrimalSimplex("End of solvePhase2") ==
      HighsDebugStatus::LOGICAL_ERROR) {
    solvePhase = SOLVE_PHASE_ERROR;
    return;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);

  if (isPrimalPhase1) {
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_DETAILED,
                      "primal-return-phase1\n");
  } else {
    if (columnIn == -1) {
      HighsPrintMessage(ekk_instance_.options_.output,
                        ekk_instance_.options_.message_level, ML_DETAILED,
                        "primal-optimal\n");
      solvePhase = SOLVE_PHASE_OPTIMAL;
    } else {
      HighsPrintMessage(ekk_instance_.options_.output,
                        ekk_instance_.options_.message_level, ML_MINIMAL,
                        "primal-unbounded\n");
      solvePhase = SOLVE_PHASE_EXIT;
      ekk_instance_.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
    }
    ekk_instance_.computeDualObjectiveValue();
  }
}

void HEkkPrimal::initialise() {
  analysis = &ekk_instance_.analysis_;

  num_col = ekk_instance_.simplex_lp_.numCol_;
  num_row = ekk_instance_.simplex_lp_.numRow_;
  num_tot = num_col + num_row;

  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance =
      ekk_instance_.options_.primal_feasibility_tolerance;
  dual_feasibility_tolerance =
      ekk_instance_.options_.dual_feasibility_tolerance;

  invertHint = INVERT_HINT_NO;

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

  num_free_col = 0;
  for (int iCol = 0; iCol < num_tot; iCol++) {
    if (ekk_instance_.simplex_info_.workLower_[iCol] == -HIGHS_CONST_INF &&
        ekk_instance_.simplex_info_.workUpper_[iCol] == HIGHS_CONST_INF) {
      // Free column
      num_free_col++;
    }
  }
  if (num_free_col) {
    HighsLogMessage(ekk_instance_.options_.logfile, HighsMessageType::INFO,
                    "HEkkPrimal:: LP has %d free columns", num_free_col);
    bool debug =
        ekk_instance_.options_.highs_debug_level > HIGHS_DEBUG_LEVEL_CHEAP;
    nonbasic_free_col_set.setup(num_free_col, num_tot,
                                ekk_instance_.options_.output, debug);
  }
}

void HEkkPrimal::rebuild() {
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
  ekk_instance_.computePrimal();
  ekk_instance_.computeSimplexPrimalInfeasible();
  isPrimalPhase1 = 0;
  if (simplex_info.num_primal_infeasibilities > 0) {
    // Primal infeasibilities so in phase 1
    isPrimalPhase1 = 1;
    if (solvePhase == SOLVE_PHASE_2)
      HighsLogMessage(
          ekk_instance_.options_.logfile, HighsMessageType::WARNING,
          "HEkkPrimal::rebuild switching back to phase 1 from phase 2");
    phase1ComputeDual();
  } else {
    // No primal infeasibilities so in phase 2. Reset costs if was
    // previously in phase 1
    if (solvePhase == SOLVE_PHASE_1) ekk_instance_.initialiseCost();
    ekk_instance_.computeDual();
  }
  ekk_instance_.computeSimplexDualInfeasible();

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

  reportRebuild(sv_invertHint);

  ekk_instance_.build_syntheticTick_ =
      ekk_instance_.factor_.build_syntheticTick;
  ekk_instance_.total_syntheticTick_ = 0;

  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HEkkPrimal::chooseColumn() {
  const vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  analysis->simplexTimerStart(ChuzcPrimalClock);
  double dBestScore = 0;
  columnIn = -1;

  // Choose any attractive nonbasic free column
  const int& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (num_nonbasic_free_col) {
    const vector<int>& nonbasic_free_col_set_entry =
        nonbasic_free_col_set.entry();
    for (int ix = 0; ix < num_nonbasic_free_col; ix++) {
      int iCol = nonbasic_free_col_set_entry[ix];
      if (fabs(workDual[iCol]) > dual_feasibility_tolerance) {
        columnIn = iCol;
        analysis->simplexTimerStop(ChuzcPrimalClock);
        return;
      }
    }
  }
  // Now look at other columns
  for (int iCol = 0; iCol < num_tot; iCol++) {
    double dMyDual = nonbasicMove[iCol] * workDual[iCol];
    double dMyScore = dMyDual / devex_weight[iCol];
    if (dMyDual < -dual_feasibility_tolerance && dMyScore < dBestScore) {
      dBestScore = dMyScore;
      columnIn = iCol;
    }
  }
  analysis->simplexTimerStop(ChuzcPrimalClock);
}

void HEkkPrimal::chooseRow() {
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<double>& baseValue = simplex_info.baseValue_;
  const vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  //
  // FTRAN
  //
  // Compute pivot column
  ekk_instance_.pivotColumnFtran(columnIn, col_aq);
  // Compute the reduced cost for the pivot column and compare it with
  // the kept value
  thetaDual = simplex_info.workDual_[columnIn];
  analysis->dualValueSignOk(ekk_instance_.options_, thetaDual, columnIn, col_aq,
                            simplex_info.workCost_,
                            ekk_instance_.simplex_basis_.basicIndex_);

  analysis->simplexTimerStart(Chuzr1Clock);
  // Initialize
  rowOut = -1;

  // Choose row pass 1
  double alphaTol = simplex_info.update_count < 10
                        ? 1e-9
                        : simplex_info.update_count < 20 ? 1e-8 : 1e-7;
  const int moveIn = thetaDual > 0 ? -1 : 1;
  if (nonbasicMove[columnIn]) assert(nonbasicMove[columnIn] == moveIn);

  double relaxTheta = 1e100;
  double relaxSpace;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    double alpha = col_aq.array[index] * moveIn;
    if (alpha > alphaTol) {
      relaxSpace =
          baseValue[index] - baseLower[index] + primal_feasibility_tolerance;
      if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    } else if (alpha < -alphaTol) {
      relaxSpace =
          baseValue[index] - baseUpper[index] - primal_feasibility_tolerance;
      if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    }
  }
  analysis->simplexTimerStop(Chuzr1Clock);

  analysis->simplexTimerStart(Chuzr2Clock);
  double bestAlpha = 0;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    double alpha = col_aq.array[index] * moveIn;
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

void HEkkPrimal::phase2Update() {
  vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  const vector<double>& workLower = ekk_instance_.simplex_info_.workLower_;
  const vector<double>& workUpper = ekk_instance_.simplex_info_.workUpper_;
  const vector<double>& baseLower = ekk_instance_.simplex_info_.baseLower_;
  const vector<double>& baseUpper = ekk_instance_.simplex_info_.baseUpper_;
  vector<double>& workValue = ekk_instance_.simplex_info_.workValue_;
  vector<double>& baseValue = ekk_instance_.simplex_info_.baseValue_;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;

  // Compute thetaPrimal
  int moveIn = thetaDual > 0 ? -1 : 1;
  if (nonbasicMove[columnIn]) assert(nonbasicMove[columnIn] == moveIn);

  if (rowOut < 0) {
    // No binding ratio in CHUZR, so flip or unbounded
    thetaPrimal = moveIn * HIGHS_CONST_INF;
  } else {
    columnOut = ekk_instance_.simplex_basis_.basicIndex_[rowOut];
    alphaCol = col_aq.array[rowOut];
    thetaPrimal = 0;
    if (alphaCol * moveIn > 0) {
      // Lower bound
      thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alphaCol;
    } else {
      // Upper bound
      thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alphaCol;
    }
  }

  // 1. Make sure it is inside bounds or just flip bound
  double lowerIn = workLower[columnIn];
  double upperIn = workUpper[columnIn];
  double valueIn = workValue[columnIn] + thetaPrimal;
  bool flipped = false;
  if (moveIn > 0) {
    if (valueIn > upperIn + primal_feasibility_tolerance) {
      // Flip to upper
      workValue[columnIn] = upperIn;
      thetaPrimal = upperIn - lowerIn;
      flipped = true;
      nonbasicMove[columnIn] = NONBASIC_MOVE_DN;
    }
  } else {
    if (valueIn < lowerIn - primal_feasibility_tolerance) {
      // Flip to lower
      workValue[columnIn] = lowerIn;
      thetaPrimal = lowerIn - upperIn;
      flipped = true;
      nonbasicMove[columnIn] = NONBASIC_MOVE_UP;
    }
  }

  if (rowOut < 0 && !flipped) {
    invertHint = INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED;
    return;
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

  // If flipped, then no need touch the pivots
  if (flipped) {
    rowOut = -1;
    alphaCol = 0;
    columnOut = columnIn;
    numericalTrouble = 0;
    iterationAnalysis();
    num_flip_since_rebuild++;
    // Update the synthetic clock
    ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
    return;
  }

  bool remove_nonbasic_free_column = nonbasicMove[columnIn] == 0;
  if (remove_nonbasic_free_column) {
    bool removed_nonbasic_free_column = nonbasic_free_col_set.remove(columnIn);
    if (!removed_nonbasic_free_column) {
      assert(removed_nonbasic_free_column);
      HighsLogMessage(
          ekk_instance_.options_.logfile, HighsMessageType::ERROR,
          "HEkkPrimal::phase2update failed to remove Nonbasic free column %d",
          columnIn);
    }
  }

  baseValue[rowOut] = valueIn;

  analysis->simplexTimerStart(CollectPrIfsClock);
  // Check for any possible infeasible
  for (int iRow = 0; iRow < num_row; iRow++) {
    if (baseValue[iRow] < baseLower[iRow] - primal_feasibility_tolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    } else if (baseValue[iRow] >
               baseUpper[iRow] + primal_feasibility_tolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    }
  }
  analysis->simplexTimerStop(CollectPrIfsClock);

  // 2. Now we can update the dual
  //
  // BTRAN
  //
  // Compute unit BTran for tableau row and FT update
  ekk_instance_.unitBtran(rowOut, row_ep);
  //
  // PRICE
  //
  ekk_instance_.computeTableauRowFromPiP(row_ep, row_ap);
  // Now update the dual values
  analysis->simplexTimerStart(UpdateDualClock);
  thetaDual = workDual[columnIn] / alphaCol;
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

  // Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  updateVerify();

  // Update the devex weight
  devexUpdate();

  // After dual update in primal simplex the dual objective value is not known
  ekk_instance_.simplex_lp_status_.has_dual_objective_value = false;

  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  // Update ekk_instance_.factor_ basis
  // Pivot in
  int sourceOut = alphaCol * moveIn > 0 ? -1 : 1;
  ekk_instance_.updatePivots(columnIn, rowOut, sourceOut);
  ekk_instance_.updateFactor(&col_aq, &row_ep, &rowOut, &invertHint);
  ekk_instance_.updateMatrix(columnIn, columnOut);
  if (simplex_info.update_count >= simplex_info.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }
  ekk_instance_.iteration_count_++;

  // Reset the devex when there are too many errors
  if (num_bad_devex_weight > 3) devexReset();

  // Report on the iteration
  iterationAnalysis();

  // Update the synthetic clock
  ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
  ekk_instance_.total_syntheticTick_ += row_ep.syntheticTick;
}

void HEkkPrimal::phase1ComputeDual() {
  const vector<double>& baseLower = ekk_instance_.simplex_info_.baseLower_;
  const vector<double>& baseUpper = ekk_instance_.simplex_info_.baseUpper_;
  const vector<double>& baseValue = ekk_instance_.simplex_info_.baseValue_;

  // Accumulate costs for checking
  vector<double>& workCost = ekk_instance_.simplex_info_.workCost_;
  workCost.assign(num_tot, 0);

  HVector buffer;
  buffer.setup(num_row);
  buffer.clear();
  for (int iRow = 0; iRow < num_row; iRow++) {
    buffer.index[iRow] = iRow;
    double cost = 0;
    if (baseValue[iRow] < baseLower[iRow] - dual_feasibility_tolerance) {
      cost = -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dual_feasibility_tolerance) {
      cost = 1.0;
    }
    buffer.array[iRow] = cost;
    workCost[ekk_instance_.simplex_basis_.basicIndex_[iRow]] = cost;
  }
  // 
  // Full BTRAN
  // 
  ekk_instance_.fullBtran(buffer);

  analysis->simplexTimerStart(PriceClock);
  HVector bufferLong;
  bufferLong.setup(num_col);
  bufferLong.clear();
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  if (simplex_info.analyse_iterations) {
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_FULL, buffer,
                                    0.0);
    analysis->num_col_price++;
  }
#endif
  ekk_instance_.matrix_.priceByColumn(bufferLong, buffer);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_FULL, row_ap);
#endif
  analysis->simplexTimerStop(PriceClock);

  const vector<int>& nonbasicFlag = ekk_instance_.simplex_basis_.nonbasicFlag_;
  vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  for (int iSeq = 0; iSeq < num_tot; iSeq++) {
    workDual[iSeq] = 0.0;
  }
  for (int iSeq = 0; iSeq < num_col; iSeq++) {
    if (nonbasicFlag[iSeq]) workDual[iSeq] = -bufferLong.array[iSeq];
  }
  for (int iRow = 0, iSeq = num_col; iRow < num_row; iRow++, iSeq++) {
    if (nonbasicFlag[iSeq]) workDual[iSeq] = -buffer.array[iRow];
  }

  ekk_instance_.computeSimplexDualInfeasible();
}

void HEkkPrimal::phase1ChooseRow() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<double>& baseValue = simplex_info.baseValue_;
  const vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  //
  // FTRAN
  //
  // Compute pivot column
  ekk_instance_.pivotColumnFtran(columnIn, col_aq);
  // Compute the reduced cost for the pivot column and compare it with
  // the kept value
  thetaDual = simplex_info.workDual_[columnIn];
  analysis->dualValueSignOk(ekk_instance_.options_, thetaDual, columnIn, col_aq,
                            simplex_info.workCost_,
                            ekk_instance_.simplex_basis_.basicIndex_);

  analysis->simplexTimerStart(Chuzr1Clock);
  // Collect phase 1 theta lists
  //
  // Determine the move direction - can't use nonbasicMove_[columnIn]
  // due to free columns
  const int moveIn = thetaDual > 0 ? -1 : 1;
  if (nonbasicMove[columnIn]) assert(nonbasicMove[columnIn] == moveIn);

  const double dPivotTol = simplex_info.update_count < 10
                               ? 1e-9
                               : simplex_info.update_count < 20 ? 1e-8 : 1e-7;
  ph1SorterR.clear();
  ph1SorterT.clear();
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    double dAlpha = col_aq.array[iRow] * moveIn;

    // When the basic variable x[i] decrease
    if (dAlpha > +dPivotTol) {
      // Whether it can become feasible by going below its upper bound
      if (baseValue[iRow] > baseUpper[iRow] + primal_feasibility_tolerance) {
        double dFeasTheta =
            (baseValue[iRow] - baseUpper[iRow] - primal_feasibility_tolerance) /
            dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow));
      }
      // Whether it can become infeasible (again) by going below its
      // lower bound
      if (baseValue[iRow] > baseLower[iRow] - primal_feasibility_tolerance &&
          baseLower[iRow] > -HIGHS_CONST_INF) {
        double dRelaxTheta =
            (baseValue[iRow] - baseLower[iRow] + primal_feasibility_tolerance) /
            dAlpha;
        double dTightTheta = (baseValue[iRow] - baseLower[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow - num_row));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow - num_row));
      }
    }

    // When the basic variable x[i] increase
    if (dAlpha < -dPivotTol) {
      // Whether it can become feasible by going above its lower bound
      if (baseValue[iRow] < baseLower[iRow] - primal_feasibility_tolerance) {
        double dFeasTheta =
            (baseValue[iRow] - baseLower[iRow] + primal_feasibility_tolerance) /
            dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow - num_row));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow - num_row));
      }
      // Whether it can become infeasible (again) by going above its
      // upper bound
      if (baseValue[iRow] < baseUpper[iRow] + primal_feasibility_tolerance &&
          baseUpper[iRow] < +HIGHS_CONST_INF) {
        double dRelaxTheta =
            (baseValue[iRow] - baseUpper[iRow] - primal_feasibility_tolerance) /
            dAlpha;
        double dTightTheta = (baseValue[iRow] - baseUpper[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow));
      }
    }
  }

  analysis->simplexTimerStop(Chuzr1Clock);
  // When there are no candidates at all, we can leave it here
  if (ph1SorterR.empty()) {
    rowOut = -1;
    columnOut = -1;
    return;
  }

  // Now sort the relaxed theta to find the final break point. TODO:
  // Consider partial sort. Or heapify [O(n)] and then pop k points
  // [kO(log(n))].

  analysis->simplexTimerStart(Chuzr2Clock);
  std::sort(ph1SorterR.begin(), ph1SorterR.end());
  double dMaxTheta = ph1SorterR.at(0).first;
  double dGradient = fabs(thetaDual);
  for (unsigned int i = 0; i < ph1SorterR.size(); i++) {
    double dMyTheta = ph1SorterR.at(i).first;
    int index = ph1SorterR.at(i).second;
    int iRow = index >= 0 ? index : index + num_row;
    dGradient -= fabs(col_aq.array[iRow]);
    // Stop when the gradient start to decrease
    if (dGradient <= 0) {
      break;
    }
    dMaxTheta = dMyTheta;
  }

  // Find out the biggest possible alpha for pivot
  std::sort(ph1SorterT.begin(), ph1SorterT.end());
  double dMaxAlpha = 0.0;
  unsigned int iLast = ph1SorterT.size();
  for (unsigned int i = 0; i < ph1SorterT.size(); i++) {
    double dMyTheta = ph1SorterT.at(i).first;
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + num_row;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    // Stop when the theta is too large
    if (dMyTheta > dMaxTheta) {
      iLast = i;
      break;
    }
    // Update the maximal possible alpha
    if (dMaxAlpha < dAbsAlpha) {
      dMaxAlpha = dAbsAlpha;
    }
  }

  // Finally choose a pivot with good enough alpha, working backwards
  rowOut = -1;
  columnOut = -1;
  phase1OutBnd = 0;
  for (int i = iLast - 1; i >= 0; i--) {
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + num_row;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    if (dAbsAlpha > dMaxAlpha * 0.1) {
      rowOut = iRow;
      phase1OutBnd = index >= 0 ? 1 : -1;
      break;
    }
  }
  if (rowOut != -1) {
    columnOut = ekk_instance_.simplex_basis_.basicIndex_[rowOut];
  }
  analysis->simplexTimerStop(Chuzr2Clock);
}

void HEkkPrimal::phase1Update() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& workLower = simplex_info.workLower_;
  const vector<double>& workUpper = simplex_info.workUpper_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  vector<double>& workDual = simplex_info.workDual_;
  vector<double>& workValue = simplex_info.workValue_;
  vector<double>& baseValue = simplex_info.baseValue_;
  vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;

  // Identify the direction of movement
  const int moveIn = thetaDual > 0 ? -1 : 1;
  if (nonbasicMove[columnIn]) assert(nonbasicMove[columnIn] == moveIn);

  // Compute the primal theta and see if we should have done a bound
  // flip instead
  alphaCol = col_aq.array[rowOut];
  thetaPrimal = 0.0;
  if (phase1OutBnd == 1) {
    thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alphaCol;
  } else {
    thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alphaCol;
  }
  assert(thetaPrimal > -HIGHS_CONST_INF && thetaPrimal < HIGHS_CONST_INF);
  double lowerIn = workLower[columnIn];
  double upperIn = workUpper[columnIn];
  double valueIn = workValue[columnIn] + thetaPrimal;
  int ifFlip = 0;
  if (moveIn == +1 && valueIn > upperIn + primal_feasibility_tolerance) {
    workValue[columnIn] = upperIn;
    thetaPrimal = upperIn - lowerIn;
    ifFlip = 1;
    nonbasicMove[columnIn] = NONBASIC_MOVE_DN;
  }
  if (moveIn == -1 && valueIn < lowerIn - primal_feasibility_tolerance) {
    workValue[columnIn] = lowerIn;
    thetaPrimal = lowerIn - upperIn;
    ifFlip = 1;
    nonbasicMove[columnIn] = NONBASIC_MOVE_UP;
  }

  // Update for the flip case
  if (ifFlip) {
    rowOut = -1;
    alphaCol = 0;
    columnOut = columnIn;
    numericalTrouble = 0;
    iterationAnalysis();
    num_flip_since_rebuild++;
    // Recompute things on flip
    if (invertHint == 0) {
      ekk_instance_.computePrimal();
      ekk_instance_.computeSimplexPrimalInfeasible();
      if (simplex_info.num_primal_infeasibilities > 0) {
        isPrimalPhase1 = 1;
        phase1ComputeDual();
      } else {
        invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
      }
    }
    // Update the synthetic clock
    ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
    return;
  }

  //
  // BTRAN
  //
  // Compute unit BTran for tableau row and FT update
  ekk_instance_.unitBtran(rowOut, row_ep);
  //
  // PRICE
  //
  ekk_instance_.computeTableauRowFromPiP(row_ep, row_ap);

  // Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  updateVerify();

  // Update the devex weight
  devexUpdate();

  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  bool remove_nonbasic_free_column = nonbasicMove[columnIn] == 0;
  if (remove_nonbasic_free_column) {
    bool removed_nonbasic_free_column = nonbasic_free_col_set.remove(columnIn);
    if (!removed_nonbasic_free_column) {
      assert(removed_nonbasic_free_column);
      HighsLogMessage(
          ekk_instance_.options_.logfile, HighsMessageType::ERROR,
          "HEkkPrimal::phase1update failed to remove Nonbasic free column %d",
          columnIn);
    }
  }

  // Update other things
  ekk_instance_.updatePivots(columnIn, rowOut, phase1OutBnd);
  ekk_instance_.updateFactor(&col_aq, &row_ep, &rowOut, &invertHint);
  ekk_instance_.updateMatrix(columnIn, columnOut);
  if (simplex_info.update_count >= simplex_info.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }
  ekk_instance_.iteration_count_++;

  // Reset the devex framework when necessary
  if (num_bad_devex_weight > 3) devexReset();

  // Report on the iteration
  iterationAnalysis();

  // Update the synthetic clock
  ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
  ekk_instance_.total_syntheticTick_ += row_ep.syntheticTick;

  // Recompute dual and primal
  if (invertHint == 0) {
    ekk_instance_.computePrimal();
    ekk_instance_.computeSimplexPrimalInfeasible();
    if (simplex_info.num_primal_infeasibilities > 0) {
      isPrimalPhase1 = 1;
      phase1ComputeDual();
    } else {
      // Crude way to force rebuild
      invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
    }
  }
}

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

void HEkkPrimal::devexUpdate() {
  // Compute the pivot weight from the reference set
  analysis->simplexTimerStart(DevexUpdateWeightClock);
  double dPivotWeight = 0.0;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    int iSeq = ekk_instance_.simplex_basis_.basicIndex_[iRow];
    double dAlpha = devex_index[iSeq] * col_aq.array[iRow];
    dPivotWeight += dAlpha * dAlpha;
  }
  dPivotWeight += devex_index[columnIn] * 1.0;
  dPivotWeight = sqrt(dPivotWeight);

  // Check if the saved weight is too large
  if (devex_weight[columnIn] > 3.0 * dPivotWeight) num_bad_devex_weight++;

  // Update the devex weight for all
  double dPivot = col_aq.array[rowOut];
  dPivotWeight /= fabs(dPivot);

  for (int i = 0; i < row_ap.count; i++) {
    int iSeq = row_ap.index[i];
    double alpha = row_ap.array[iSeq];
    double devex = dPivotWeight * fabs(alpha);
    devex += devex_index[iSeq] * 1.0;
    if (devex_weight[iSeq] < devex) {
      devex_weight[iSeq] = devex;
    }
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iPtr = row_ep.index[i];
    int iSeq = row_ep.index[i] + num_col;
    double alpha = row_ep.array[iPtr];
    double devex = dPivotWeight * fabs(alpha);
    devex += devex_index[iSeq] * 1.0;
    if (devex_weight[iSeq] < devex) {
      devex_weight[iSeq] = devex;
    }
  }

  // Update devex weight for the pivots
  devex_weight[columnOut] = max(1.0, dPivotWeight);
  devex_weight[columnIn] = 1.0;
  num_devex_iterations++;
  analysis->simplexTimerStop(DevexUpdateWeightClock);
}

void HEkkPrimal::updateVerify() {
  // updateVerify for primal
  numericalTrouble = 0;
  double abs_alpha_from_col = fabs(alphaCol);
  bool column_in = columnIn < num_col;
  std::string alphaRow_source;
  if (column_in) {
    alphaRow = row_ap.array[columnIn];
    alphaRow_source = "Col";
  } else {
    alphaRow = row_ep.array[columnIn - num_col];
    alphaRow_source = "Row";
  }
  double abs_alpha_from_row = fabs(alphaRow);
  double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  double min_abs_alpha = min(abs_alpha_from_col, abs_alpha_from_row);
  numericalTrouble = abs_alpha_diff / min_abs_alpha;
  if (numericalTrouble > 1e-7)
    printf(
        "Numerical check: Iter %4d: alphaCol = %12g, (From %3s alphaRow = "
        "%12g), aDiff = %12g: measure = %12g\n",
        ekk_instance_.iteration_count_, alphaCol, alphaRow_source.c_str(),
        alphaRow, abs_alpha_diff, numericalTrouble);
  assert(numericalTrouble < 1e-3);
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  //
  //  if (numericalTrouble > 1e-7 && ekk_instance_.simplex_info_.update_count >
  //  0) invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
}

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
  analysis->pivot_value_from_column = alphaCol;
  analysis->pivot_value_from_row = alphaRow;
  analysis->numerical_trouble = 0;  // numericalTrouble;
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

void HEkkPrimal::getNonbasicFreeColumnSet() {
  if (!num_free_col) return;
  assert(num_free_col > 0);
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const SimplexBasis& simplex_basis = ekk_instance_.simplex_basis_;
  nonbasic_free_col_set.clear();
  for (int iVar = 0; iVar < num_tot; iVar++) {
    bool nonbasic_free =
        simplex_basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_TRUE &&
        simplex_info.workLower_[iVar] <= -HIGHS_CONST_INF &&
        simplex_info.workUpper_[iVar] >= HIGHS_CONST_INF;
    if (nonbasic_free) nonbasic_free_col_set.add(iVar);
  }
}

HighsDebugStatus HEkkPrimal::debugPrimalSimplex(const std::string message) {
  HighsDebugStatus return_status =
      ekkDebugSimplex(message, ekk_instance_, algorithm, solvePhase);
  if (return_status == HighsDebugStatus::LOGICAL_ERROR) return return_status;
  return_status = ekkDebugNonbasicFreeColumnSet(ekk_instance_, num_free_col,
                                                nonbasic_free_col_set);
  if (return_status == HighsDebugStatus::LOGICAL_ERROR) return return_status;
  return HighsDebugStatus::OK;
}
