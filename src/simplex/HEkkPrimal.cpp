/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*      <                                                                 */
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
#include "simplex/SimplexTimer.h"
#include "util/HighsSort.h"

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

  localReportIter(true);
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
      solvePhase1();
      simplex_info.primal_phase1_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else if (solvePhase == SOLVE_PHASE_2) {
      // Phase 2
      solvePhase2();
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

  rebuild_reason = REBUILD_REASON_NO;

  ekk_instance_.simplex_lp_status_.has_primal_objective_value = false;
  ekk_instance_.simplex_lp_status_.has_dual_objective_value = false;
  ekk_instance_.scaled_model_status_ = HighsModelStatus::NOTSET;
  ekk_instance_.solve_bailout_ = false;

  // Setup local vectors
  col_aq.setup(num_row);
  row_ep.setup(num_row);
  row_ap.setup(num_col);
  col_basic_feasibility_change.setup(num_row);
  row_basic_feasibility_change.setup(num_col);

  ph1SorterR.reserve(num_row);
  ph1SorterT.reserve(num_row);

  resetDevex();

  num_free_col = 0;
  for (int iCol = 0; iCol < num_tot; iCol++) {
    if (ekk_instance_.simplex_info_.workLower_[iCol] == -HIGHS_CONST_INF &&
        ekk_instance_.simplex_info_.workUpper_[iCol] == HIGHS_CONST_INF) {
      // Free column
      num_free_col++;
    }
  }
  // Set up the HSet instances, possibly using the internal error reporting and
  // debug option
  const bool debug =
      ekk_instance_.options_.highs_debug_level > HIGHS_DEBUG_LEVEL_CHEAP;
  FILE* output = ekk_instance_.options_.output;
  if (num_free_col) {
    HighsLogMessage(ekk_instance_.options_.logfile, HighsMessageType::INFO,
                    "HEkkPrimal:: LP has %d free columns", num_free_col);
    nonbasic_free_col_set.setup(num_free_col, num_tot, output, debug);
  }
  // Set up the hyper-sparse CHUZC data
  hyper_chuzc_candidate.resize(1 + max_num_hyper_chuzc_candidates);
  hyper_chuzc_measure.resize(1 + max_num_hyper_chuzc_candidates);
  hyper_chuzc_candidate_set.setup(max_num_hyper_chuzc_candidates, num_tot,
                                  output, debug);
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
    rebuild();
    if (solvePhase == SOLVE_PHASE_ERROR) return;
    if (!is_primal_phase1_) {
      // No primal infeasibilities found in rebuild() so break and
      // return to phase 2
      solvePhase = SOLVE_PHASE_2;
      break;
    }
    if (ekk_instance_.bailoutOnTimeIterations()) return;

    for (;;) {
      iterate();
      if (ekk_instance_.bailoutOnTimeIterations()) return;
      if (solvePhase == SOLVE_PHASE_ERROR) return;
      if (rebuild_reason) break;
    }
    // Go to the next rebuild
    if (rebuild_reason) {
      // Stop when the invert is new
      if (simplex_lp_status.has_fresh_rebuild) break;
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
    rebuild();
    if (ekk_instance_.bailoutOnTimeIterations()) return;
    if (solvePhase == SOLVE_PHASE_ERROR) return;

    if (is_primal_phase1_) {
      // Primal infeasibilities found in rebuild() Should be
      // shifted but, for now, break and return to phase 1
      solvePhase = SOLVE_PHASE_1;
      break;
    }

    for (;;) {
      iterate();
      if (ekk_instance_.bailoutOnTimeIterations()) return;
      if (solvePhase == SOLVE_PHASE_ERROR) return;
      if (rebuild_reason) break;
    }
    // If the data are fresh from rebuild() and no flips have occurred, break
    // out of the outer loop to see what's ocurred
    if (simplex_lp_status.has_fresh_rebuild && num_flip_since_rebuild == 0)
      break;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  if (debugPrimalSimplex("End of solvePhase2") ==
      HighsDebugStatus::LOGICAL_ERROR) {
    solvePhase = SOLVE_PHASE_ERROR;
    return;
  }
  if (is_primal_phase1_) {
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_DETAILED,
                      "primal-return-phase1\n");
  } else {
    if (variable_in == -1) {
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
  int sv_rebuild_reason = rebuild_reason;
  rebuild_reason = REBUILD_REASON_NO;
  // Possibly Rebuild factor
  bool reInvert = simplex_info.update_count > 0;
  if (!invert_if_row_out_negative) {
    // Don't reinvert if variable_in is negative [equivalently, if sv_rebuild_reason ==
    // REBUILD_REASON_POSSIBLY_OPTIMAL]
    if (sv_rebuild_reason == REBUILD_REASON_POSSIBLY_OPTIMAL) {
      assert(variable_in == -1);
      reInvert = false;
    }
  }
  if (reInvert) {
    int rank_deficiency = ekk_instance_.computeFactor();
    if (rank_deficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    simplex_info.update_count = 0;
  }
  ekk_instance_.computePrimal();
  getBasicPrimalInfeasibility();
  is_primal_phase1_ = false;
  if (simplex_info.num_primal_infeasibilities > 0) {
    // Primal infeasibilities so in phase 1
    is_primal_phase1_ = true;
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
    ekk_instance_.computeSimplexDualInfeasible();
  }
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

  reportRebuild(sv_rebuild_reason);

  // Record the synthetic clock for INVERT, and zero it for UPDATE
  ekk_instance_.build_syntheticTick_ =
      ekk_instance_.factor_.build_syntheticTick;
  ekk_instance_.total_syntheticTick_ = 0;

  // Determine whether to use hyper-sparse CHUZC
  if (solvePhase == SOLVE_PHASE_1) {
    use_hyper_chuzc = true;
  } else {
    use_hyper_chuzc = true;
  }
  hyperChooseColumnClear();

  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HEkkPrimal::iterate() {
  //
  bool check = ekk_instance_.iteration_count_ >= check_iter;
  if (check) {
    printf("Iter %d\n", ekk_instance_.iteration_count_);
    ekk_instance_.options_.highs_debug_level = HIGHS_DEBUG_LEVEL_EXPENSIVE;
  }
  if (debugPrimalSimplex("Before iteration") ==
      HighsDebugStatus::LOGICAL_ERROR) {
    solvePhase = SOLVE_PHASE_ERROR;
    return;
  }

  // Perform CHUZC
  if (solvePhase == SOLVE_PHASE_1) {
    chooseColumn();
    if (variable_in == -1) {
      rebuild_reason = REBUILD_REASON_CHOOSE_COLUMN_FAIL;
      return;
    }
  } else {
    chuzc();
    if (variable_in == -1) {
      rebuild_reason = REBUILD_REASON_POSSIBLY_OPTIMAL;
      return;
    }
  }

  // Perform FTRAN - and dual value cross-check
  assessVariableIn();
  if (solvePhase == SOLVE_PHASE_ERROR) return;

  // Perform CHUZR
  if (solvePhase == SOLVE_PHASE_1) {
    phase1ChooseRow();
    if (row_out == -1) {
      HighsLogMessage(ekk_instance_.options_.logfile, HighsMessageType::ERROR,
		      "Primal phase 1 choose row failed");
      solvePhase = SOLVE_PHASE_ERROR;
      return;
    }
  } else {
    chooseRow();
  }
  if (use_considerBoundSwap) considerBoundSwap();
  // If pivoting, perform unit BTRAN and PRICE to get pivotal row - and do numerical check
  if (row_out >= 0) {
    assessPivot();
    if (rebuild_reason) {
      assert(rebuild_reason==REBUILD_REASON_POSSIBLY_SINGULAR_BASIS);
      return;
    }
  }
  // Any pivoting is numerically acceptable, so perform update
  update();
}

void HEkkPrimal::assessVariableIn() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  vector<double>& workDual = simplex_info.workDual_;
  vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  theta_dual = workDual[variable_in];
  // Determine the move direction - can't use nonbasicMove_[variable_in]
  // due to free columns
  move_in = theta_dual > 0 ? -1 : 1;
  // Unless the variable is free, nonbasicMove[variable_in] should be the same as move_in
  if (nonbasicMove[variable_in]) assert(nonbasicMove[variable_in] == move_in);
  //
  // FTRAN
  //
  // Compute pivot column
  ekk_instance_.pivotColumnFtran(variable_in, col_aq);
  // Compute the dual for the pivot column and compare it with the
  // updated value
  double computed_theta_dual;
  bool theta_dual_sign_ok = 
    analysis->dualValueSignOk(ekk_instance_.options_, theta_dual, variable_in, col_aq,
			      simplex_info.workCost_,
			      ekk_instance_.simplex_basis_.basicIndex_,
			      computed_theta_dual);
  // Really should do something if theta_dual_sign_ok is false
  if (!theta_dual_sign_ok) {
    HighsLogMessage(ekk_instance_.options_.logfile, HighsMessageType::ERROR,
		    "Primal simplex: Computed / update dual of entering variable are %g / %g: sign error",
		    computed_theta_dual, theta_dual);
    solvePhase = SOLVE_PHASE_ERROR;
    return;
  }
  // Feed in the computed dual value
  //  simplex_info.workDual_[variable_in] = computed_theta_dual;
  //  theta_dual = simplex_info.workDual_[variable_in];
  
}

void HEkkPrimal::considerBoundSwap() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& workDual = simplex_info.workDual_;
  const vector<double>& workLower = simplex_info.workLower_;
  const vector<double>& workUpper = simplex_info.workUpper_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  vector<double>& workValue = simplex_info.workValue_;
  vector<double>& baseValue = simplex_info.baseValue_;
  vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  // Compute the primal theta and see if we should have done a bound
  // flip instead
  if (row_out < 0) {
    assert(solvePhase == SOLVE_PHASE_2);
    // No binding ratio in CHUZR, so flip or unbounded
    theta_primal = move_in * HIGHS_CONST_INF;
    move_out = 0;
  } else {
    // Determine the step to the leaving bound
    //
    alpha_col = col_aq.array[row_out];
    // In Phase 1, move_out depends on whether the leaving variable is
    // becoming feasible - moves up to lower (down to upper) - or
    // remaining feasible - moves down to lower (up to upper) - so
    // can't be set so easily as in phase 2
    if (solvePhase == SOLVE_PHASE_2) 
      move_out = alpha_col * move_in > 0 ? -1 : 1;
    theta_primal = 0;
    if (move_out == 1) {
      theta_primal = (baseValue[row_out] - baseUpper[row_out]) / alpha_col;
    } else {
      theta_primal = (baseValue[row_out] - baseLower[row_out]) / alpha_col;
    }
    assert(theta_primal > -HIGHS_CONST_INF && theta_primal < HIGHS_CONST_INF);
  }

  // Look to see if there is a bound flip
  bool flipped = false;
  double lower_in = workLower[variable_in];
  double upper_in = workUpper[variable_in];
  value_in = workValue[variable_in] + theta_primal;
  if (move_in > 0) {
    if (value_in > upper_in + primal_feasibility_tolerance) {
      // Flip to upper
      value_in = upper_in;
      if (!flip_update_later) {
	workValue[variable_in] = value_in;
	nonbasicMove[variable_in] = NONBASIC_MOVE_DN;
      }
      theta_primal = upper_in - lower_in;
      flipped = true;
    }
  } else {
    if (value_in < lower_in - primal_feasibility_tolerance) {
      // Flip to lower
      value_in = lower_in;
      if (!flip_update_later) {
	workValue[variable_in] = value_in;
	nonbasicMove[variable_in] = NONBASIC_MOVE_UP;
      }
      theta_primal = lower_in - upper_in;
      flipped = true;
    }
  }
  if (flipped) {
    if (!flip_update_later) assert(nonbasicMove[variable_in] = -move_in);
    row_out = -1;
    variable_out = variable_in;
    alpha_col = 0;
    numericalTrouble = 0;
    if (flip_update_later) {
      workValue[variable_in] = value_in;
      assert(nonbasicMove[variable_in] = move_in);
      nonbasicMove[variable_in] = -move_in;
    }
  }

  const bool pivot_or_flipped = row_out >= 0 || flipped;
  if (solvePhase == SOLVE_PHASE_2) {
    // Check for possible unboundedness
    if (!pivot_or_flipped) {
      rebuild_reason = REBUILD_REASON_POSSIBLY_PRIMAL_UNBOUNDED;
      return;
    }
  }
  // Check for possible error
  assert(pivot_or_flipped);
  assert(flipped == (row_out == -1));

}

void HEkkPrimal::assessPivot() {
  assert(row_out >= 0);
  // Record the pivot entry
  alpha_col = col_aq.array[row_out];
  variable_out = ekk_instance_.simplex_basis_.basicIndex_[row_out];
  
  // Compute the tableau row
  //
  // BTRAN
  //
  // Compute unit BTran for tableau row and FT update
  ekk_instance_.unitBtran(row_out, row_ep);
  //
  // PRICE
  //
  ekk_instance_.tableauRowPrice(row_ep, row_ap);
  
  // Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  updateVerify();
  
}

void HEkkPrimal::update() {
  // Perform update operations that are independent of phase
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& workDual = simplex_info.workDual_;
  const vector<double>& workLower = simplex_info.workLower_;
  const vector<double>& workUpper = simplex_info.workUpper_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  vector<double>& workValue = simplex_info.workValue_;
  vector<double>& baseValue = simplex_info.baseValue_;
  vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;

  if (!use_considerBoundSwap) {
  // Compute the primal theta and see if we should have done a bound
  // flip instead
  if (row_out < 0) {
    assert(solvePhase == SOLVE_PHASE_2);
    // No binding ratio in CHUZR, so flip or unbounded
    theta_primal = move_in * HIGHS_CONST_INF;
    move_out = 0;
  } else {
    // Determine the step to the leaving bound
    //
    // Check that the pivot entry is known
    assert(alpha_col == col_aq.array[row_out]);
    // In Phase 1, move_out depends on whether the leaving variable is
    // becoming feasible - moves up to lower (down to upper) - or
    // remaining feasible - moves down to lower (up to upper) - so
    // can't be set so easily as in phase 2
    if (solvePhase == SOLVE_PHASE_2) 
      move_out = alpha_col * move_in > 0 ? -1 : 1;
    theta_primal = 0;
    if (move_out == 1) {
      theta_primal = (baseValue[row_out] - baseUpper[row_out]) / alpha_col;
    } else {
      theta_primal = (baseValue[row_out] - baseLower[row_out]) / alpha_col;
    }
    assert(theta_primal > -HIGHS_CONST_INF && theta_primal < HIGHS_CONST_INF);
  }

  // Look to see if there is a bound flip
  bool flipped = false;
  double lower_in = workLower[variable_in];
  double upper_in = workUpper[variable_in];
  value_in = workValue[variable_in] + theta_primal;
  if (move_in > 0) {
    if (value_in > upper_in + primal_feasibility_tolerance) {
      // Flip to upper
      value_in = upper_in;
      if (!flip_update_later) {
	workValue[variable_in] = value_in;
	nonbasicMove[variable_in] = NONBASIC_MOVE_DN;
      }
      theta_primal = upper_in - lower_in;
      flipped = true;
    }
  } else {
    if (value_in < lower_in - primal_feasibility_tolerance) {
      // Flip to lower
      value_in = lower_in;
      if (!flip_update_later) {
	workValue[variable_in] = value_in;
	nonbasicMove[variable_in] = NONBASIC_MOVE_UP;
      }
      theta_primal = lower_in - upper_in;
      flipped = true;
    }
  }
  if (flipped) {
    if (!flip_update_later) assert(nonbasicMove[variable_in] = -move_in);
    row_out = -1;
    variable_out = variable_in;
    alpha_col = 0;
    numericalTrouble = 0;
    if (flip_update_later) {
      workValue[variable_in] = value_in;
      assert(nonbasicMove[variable_in] = move_in);
      nonbasicMove[variable_in] = -move_in;
    }
  }

  const bool pivot_or_flipped = row_out >= 0 || flipped;
  if (solvePhase == SOLVE_PHASE_2) {
    // Check for possible unboundedness
    if (!pivot_or_flipped) {
      rebuild_reason = REBUILD_REASON_POSSIBLY_PRIMAL_UNBOUNDED;
      return;
    }
  }
  // Check for possible error
  assert(pivot_or_flipped);
  assert(flipped == (row_out == -1));
  }

  bool flipped = row_out < 0;
  // Start hyper-sparse CHUZC, that takes place through phase1Update()
  hyperChooseColumnStart();

  if (solvePhase == SOLVE_PHASE_1) {
    
    // Update primal values
    phase1UpdatePrimal();
    
    // Update the duals with respect to feasibility changes
    basicFeasibilityChangeUpdateDual();
    
    // For hyper-sparse CHUZC, analyse the duals that have just changed
    hyperChooseColumnBasicFeasibilityChange();
    
    // Update for the flip case
    if (flipped) {
      simplex_info.primal_bound_swap++;
      ekk_instance_.invalidateDualInfeasibilityRecord();
      iterationAnalysis();
      localReportIter();
      num_flip_since_rebuild++;
      // Recompute things on flip
      if (rebuild_reason == 0) {
	if (ekk_instance_.simplex_info_.num_primal_infeasibilities > 0) {
	  is_primal_phase1_ = true;
	} else {
	  // Crude way to force rebuild
	  rebuild_reason = REBUILD_REASON_UPDATE_LIMIT_REACHED;
	}
      }
      // Update the synthetic clock for UPDATE
      ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
      return;
    }
    
    // Now set the value of the entering variable
    assert(row_out >= 0);
    baseValue[row_out] = value_in;
    // Consider whether the entering value is feasible and, if not, take
    // action
    considerInfeasibleValueIn();
    
  } else {
    // Phase 2
    
    //
    // Update primal values, and identify any infeasibilities
    //
    phase2UpdatePrimal();
    
    // Why is the detailed primal infeasibility information needed?
    //  ekk_instance_.computeSimplexPrimalInfeasible();
    
    // If flipped, then no need touch the pivots
    if (flipped) {
      simplex_info.primal_bound_swap++;
      ekk_instance_.invalidateDualInfeasibilityRecord();
      iterationAnalysis();
      localReportIter();
      num_flip_since_rebuild++;
      // Update the synthetic clock for UPDATE
      ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
      return;
    }
    // Remaining update operations are independent of phase
  }

  // Update the dual values
  theta_dual = workDual[variable_in];
  updateDual();
  
  // Update the devex weight
  updateDevex();
  
  // If entering column was nonbasic free, remove it from the set
  removeNonbasicFreeColumn();
  
  // For hyper-sparse CHUZC, analyse the duals and weights that have
  // just changed
  hyperChooseColumnDualChange();
  
  // Perform pivoting
  ekk_instance_.updatePivots(variable_in, row_out, move_out);
  ekk_instance_.updateFactor(&col_aq, &row_ep, &row_out, &rebuild_reason);
  ekk_instance_.updateMatrix(variable_in, variable_out);
  if (simplex_info.update_count >= simplex_info.update_limit)
    rebuild_reason = REBUILD_REASON_UPDATE_LIMIT_REACHED;
  
  // Update the iteration count
  ekk_instance_.iteration_count_++;
  
  // Reset the devex when there are too many errors
  if (num_bad_devex_weight > allowed_num_bad_devex_weight) resetDevex();
  
  // Report on the iteration
  iterationAnalysis();
  localReportIter();
  
  // Update the synthetic clock for UPDATE
  ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
  ekk_instance_.total_syntheticTick_ += row_ep.syntheticTick;
  
  // Perform hyper-sparse CHUZC
  hyperChooseColumn();

  if (solvePhase == SOLVE_PHASE_1) {
    // Possibly rebuild if not already doing it and there are primal
    // infeasibilities
    if (!rebuild_reason) {
      if (ekk_instance_.simplex_info_.num_primal_infeasibilities > 0) {
	is_primal_phase1_ = true;
      } else {
	// Crude way to force rebuild
	rebuild_reason = REBUILD_REASON_UPDATE_LIMIT_REACHED;
      }
    }
  }
}

void HEkkPrimal::chuzc() {
  if (done_next_chuzc) assert(use_hyper_chuzc);
  if (!done_next_chuzc) chooseColumn(true);
  if (use_hyper_chuzc) {
    const bool check_hyper_chuzc = true;
    if (check_hyper_chuzc) {
      int hyper_sparse_variable_in = variable_in;
      chooseColumn(false);
      HighsSimplexInfo& info = ekk_instance_.simplex_info_;
      double hyper_sparse_measure = 0;
      if (hyper_sparse_variable_in >= 0)
        hyper_sparse_measure = fabs(info.workDual_[hyper_sparse_variable_in]) /
	  devex_weight[hyper_sparse_variable_in];
      double measure = 0;
      if (variable_in >= 0)
        measure = fabs(info.workDual_[variable_in]) / devex_weight[variable_in];
      if (hyper_sparse_measure != measure) {
        printf(
	       "Iteration %d: Hyper-sparse CHUZC measure %g != %g = Full "
	       "CHUZC measure (%d, %d)\n",
	       ekk_instance_.iteration_count_, hyper_sparse_measure, measure,
	       hyper_sparse_variable_in, variable_in);
        assert(hyper_sparse_measure == measure);
      }
      variable_in = hyper_sparse_variable_in;
    }
  }
}

void HEkkPrimal::chooseColumn(const bool hyper_sparse) {
  assert(!hyper_sparse || !done_next_chuzc);
  const vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  double best_measure = 0;
  variable_in = -1;
  
  const bool local_use_hyper_chuzc = hyper_sparse;
  // Consider nonbasic free columns first
  const int& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (local_use_hyper_chuzc) {
    if (!initialise_hyper_chuzc) hyperChooseColumn();
    if (initialise_hyper_chuzc) {
      analysis->simplexTimerStart(ChuzcHyperInitialiselClock);
      num_hyper_chuzc_candidates = 0;
      if (num_nonbasic_free_col) {
        const vector<int>& nonbasic_free_col_set_entry =
	  nonbasic_free_col_set.entry();
        for (int ix = 0; ix < num_nonbasic_free_col; ix++) {
          int iCol = nonbasic_free_col_set_entry[ix];
          double dual_infeasibility = fabs(workDual[iCol]);
          if (dual_infeasibility > dual_feasibility_tolerance) {
            double measure = dual_infeasibility / devex_weight[iCol];
            addToDecreasingHeap(
				num_hyper_chuzc_candidates, max_num_hyper_chuzc_candidates,
				hyper_chuzc_measure, hyper_chuzc_candidate, measure, iCol);
          }
        }
      }
      // Now look at other columns
      for (int iCol = 0; iCol < num_tot; iCol++) {
        double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
        if (dual_infeasibility > dual_feasibility_tolerance) {
          double measure = dual_infeasibility / devex_weight[iCol];
          addToDecreasingHeap(
			      num_hyper_chuzc_candidates, max_num_hyper_chuzc_candidates,
			      hyper_chuzc_measure, hyper_chuzc_candidate, measure, iCol);
        }
      }
      // Sort the heap
      sortDecreasingHeap(num_hyper_chuzc_candidates, hyper_chuzc_measure,
                         hyper_chuzc_candidate);
      initialise_hyper_chuzc = false;
      analysis->simplexTimerStop(ChuzcHyperInitialiselClock);
      // Choose the first entry - if there is one
      if (num_hyper_chuzc_candidates) {
        variable_in = hyper_chuzc_candidate[1];
        best_measure = hyper_chuzc_measure[1];
        max_hyper_chuzc_non_candidate_measure =
	  hyper_chuzc_measure[num_hyper_chuzc_candidates];
        if (report_hyper_chuzc)
          printf(
		 "Full CHUZC: Max         measure is %9.4g for column %4d, and "
		 "max non-candiate measure of  %9.4g\n",
		 best_measure, variable_in, max_hyper_chuzc_non_candidate_measure);
      }
    }
  } else {
    analysis->simplexTimerStart(ChuzcPrimalClock);
    // Choose any attractive nonbasic free column
    if (num_nonbasic_free_col) {
      const vector<int>& nonbasic_free_col_set_entry =
	nonbasic_free_col_set.entry();
      for (int ix = 0; ix < num_nonbasic_free_col; ix++) {
        int iCol = nonbasic_free_col_set_entry[ix];
        double dual_infeasibility = fabs(workDual[iCol]);
        if (dual_infeasibility > dual_feasibility_tolerance &&
            dual_infeasibility > best_measure * devex_weight[iCol]) {
          variable_in = iCol;
          best_measure = dual_infeasibility / devex_weight[iCol];
        }
      }
    }
    // Now look at other columns
    for (int iCol = 0; iCol < num_tot; iCol++) {
      double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
      if (dual_infeasibility > dual_feasibility_tolerance &&
          dual_infeasibility > best_measure * devex_weight[iCol]) {
        variable_in = iCol;
        best_measure = dual_infeasibility / devex_weight[iCol];
      }
    }
    analysis->simplexTimerStop(ChuzcPrimalClock);
  }
  //  printf("ChooseColumn: Iteration %d, choose column %d with measure %g\n",
  //	 ekk_instance_.iteration_count_, variable_in, best_measure);
}

void HEkkPrimal::hyperChooseColumn() {
  if (!use_hyper_chuzc) return;
  if (initialise_hyper_chuzc) return;
  analysis->simplexTimerStart(ChuzcHyperClock);
  const vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  const vector<int>& nonbasicFlag = ekk_instance_.simplex_basis_.nonbasicFlag_;
  const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  if (report_hyper_chuzc)
    printf("H-S  CHUZC: Max changed measure is %9.4g for column %4d",
           max_changed_measure_value, max_changed_measure_column);
  double best_measure = max_changed_measure_value;
  variable_in = max_changed_measure_column;
  const bool consider_nonbasic_free_column = nonbasic_free_col_set.count();
  if (num_hyper_chuzc_candidates) {
    for (int iEntry = 1; iEntry <= num_hyper_chuzc_candidates; iEntry++) {
      int iCol = hyper_chuzc_candidate[iEntry];
      if (nonbasicFlag[iCol] == NONBASIC_FLAG_FALSE) {
        assert(!nonbasicMove[iCol]);
        continue;
      }
      // Assess any dual infeasibility
      double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
      if (consider_nonbasic_free_column) {
        if (nonbasic_free_col_set.in(iCol))
          dual_infeasibility = fabs(workDual[iCol]);
      }
      if (dual_infeasibility > dual_feasibility_tolerance) {
        if (dual_infeasibility > best_measure * devex_weight[iCol]) {
          best_measure = dual_infeasibility / devex_weight[iCol];
          variable_in = iCol;
        }
      }
    }
  }
  if (variable_in != max_changed_measure_column) {
    if (report_hyper_chuzc)
      printf(", and after HS CHUZC set it is now %9.4g for column %4d",
             best_measure, variable_in);
    max_hyper_chuzc_non_candidate_measure =
      max(max_changed_measure_value, max_hyper_chuzc_non_candidate_measure);
  }
  if (best_measure >= max_hyper_chuzc_non_candidate_measure) {
    // Candidate is at least as good as any unknown column, so accept it
    done_next_chuzc = true;
    if (report_hyper_chuzc)
      printf(", and no       has  measure >  %9.4g\n",
             max_hyper_chuzc_non_candidate_measure);
  } else {
    // Candidate isn't as good as best unknown column, so do a full CHUZC
    // Shouldn't claim to have done the next CHUZC
    assert(!done_next_chuzc);
    done_next_chuzc = false;
    initialise_hyper_chuzc = true;
    if (report_hyper_chuzc)
      printf(", but some may have measure >= %9.4g\n",
             max_hyper_chuzc_non_candidate_measure);
  }
  analysis->simplexTimerStop(ChuzcHyperClock);
}

void HEkkPrimal::hyperChooseColumnStart() {
  max_changed_measure_value = 0;
  max_changed_measure_column = -1;
  done_next_chuzc = false;
}

void HEkkPrimal::hyperChooseColumnClear() {
  initialise_hyper_chuzc = use_hyper_chuzc;
  max_hyper_chuzc_non_candidate_measure = -1;
  done_next_chuzc = false;
}

void HEkkPrimal::hyperChooseColumnChangedInfeasibility(
						       const double infeasibility, const int iCol) {
  if (infeasibility > max_changed_measure_value * devex_weight[iCol]) {
    max_hyper_chuzc_non_candidate_measure =
      max(max_changed_measure_value, max_hyper_chuzc_non_candidate_measure);
    max_changed_measure_value = infeasibility / devex_weight[iCol];
    max_changed_measure_column = iCol;
  } else if (infeasibility >
             max_hyper_chuzc_non_candidate_measure * devex_weight[iCol]) {
    max_hyper_chuzc_non_candidate_measure = infeasibility / devex_weight[iCol];
  }
}

void HEkkPrimal::hyperChooseColumnBasicFeasibilityChange() {
  if (!use_hyper_chuzc) return;
  analysis->simplexTimerStart(ChuzcHyperBasicFeasibilityChangeClock);
  const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  const vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  int to_entry;
  const bool use_row_indices = ekk_instance_.sparseLoopStyle(
							     row_basic_feasibility_change.count, num_col, to_entry);
  for (int iEntry = 0; iEntry < to_entry; iEntry++) {
    int iCol;
    if (use_row_indices) {
      iCol = row_basic_feasibility_change.index[iEntry];
    } else {
      iCol = iEntry;
    }
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  const bool use_col_indices = ekk_instance_.sparseLoopStyle(
							     col_basic_feasibility_change.count, num_row, to_entry);
  for (int iEntry = 0; iEntry < to_entry; iEntry++) {
    int iRow;
    if (use_col_indices) {
      iRow = col_basic_feasibility_change.index[iEntry];
    } else {
      iRow = iEntry;
    }
    int iCol = num_col + iRow;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Any nonbasic free columns will be handled explicitly in
  // hyperChooseColumnDualChange, so only look at them here if not
  // flipping
  const int& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (row_out < 0 && num_nonbasic_free_col) {
    const vector<int>& nonbasic_free_col_set_entry =
      nonbasic_free_col_set.entry();
    for (int iEntry = 0; iEntry < num_nonbasic_free_col; iEntry++) {
      int iCol = nonbasic_free_col_set_entry[iEntry];
      double dual_infeasibility = fabs(workDual[iCol]);
      if (dual_infeasibility > dual_feasibility_tolerance)
        hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
    }
  }
  analysis->simplexTimerStop(ChuzcHyperBasicFeasibilityChangeClock);
}

void HEkkPrimal::hyperChooseColumnDualChange() {
  if (!use_hyper_chuzc) return;
  analysis->simplexTimerStart(ChuzcHyperDualClock);
  const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  const vector<int>& nonbasicMove = ekk_instance_.simplex_basis_.nonbasicMove_;
  int to_entry;
  // Look at changes in the columns and assess any dual infeasibility
  const bool use_row_indices =
    ekk_instance_.sparseLoopStyle(row_ap.count, num_col, to_entry);
  for (int iEntry = 0; iEntry < to_entry; iEntry++) {
    int iCol;
    if (use_row_indices) {
      iCol = row_ap.index[iEntry];
    } else {
      iCol = iEntry;
    }
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (iCol == check_column && ekk_instance_.iteration_count_ >= check_iter) {
      double measure = dual_infeasibility / devex_weight[iCol];
      if (report_hyper_chuzc) {
        printf("Changing column %d: measure = %g \n", check_column, measure);
      }
    }
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Look at changes in the rows and assess any dual infeasibility
  const bool use_col_indices =
    ekk_instance_.sparseLoopStyle(row_ep.count, num_row, to_entry);
  for (int iEntry = 0; iEntry < to_entry; iEntry++) {
    int iRow;
    if (use_col_indices) {
      iRow = row_ep.index[iEntry];
    } else {
      iRow = iEntry;
    }
    int iCol = iRow + num_col;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (iCol == check_column && ekk_instance_.iteration_count_ >= check_iter) {
      double measure = dual_infeasibility / devex_weight[iCol];
      if (report_hyper_chuzc) {
        printf("Changing column %d: measure = %g \n", check_column, measure);
      }
    }
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Look for measure changes in any nonbasic free columns and assess
  // any dual infeasibility
  const int& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (num_nonbasic_free_col) {
    const vector<int>& nonbasic_free_col_set_entry =
      nonbasic_free_col_set.entry();
    for (int iEntry = 0; iEntry < num_nonbasic_free_col; iEntry++) {
      int iCol = nonbasic_free_col_set_entry[iEntry];
      double dual_infeasibility = fabs(workDual[iCol]);
      if (dual_infeasibility > dual_feasibility_tolerance)
        hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
    }
  }
  // Assess any dual infeasibility for the leaving column - should be dual
  // feasible!
  int iCol = variable_out;
  double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
  if (dual_infeasibility > dual_feasibility_tolerance) {
    printf("Dual infeasibility %g for leaving column!\n", dual_infeasibility);
    assert(dual_infeasibility <= dual_feasibility_tolerance);
    hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  analysis->simplexTimerStop(ChuzcHyperDualClock);
}

void HEkkPrimal::chooseRow() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<double>& baseValue = simplex_info.baseValue_;
  analysis->simplexTimerStart(Chuzr1Clock);
  // Initialize
  row_out = -1;
  
  // Choose row pass 1
  double alphaTol = simplex_info.update_count < 10
						? 1e-9
						: simplex_info.update_count < 20 ? 1e-8 : 1e-7;
  
  double relaxTheta = 1e100;
  double relaxSpace;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    double alpha = col_aq.array[iRow] * move_in;
    if (alpha > alphaTol) {
      relaxSpace =
	baseValue[iRow] - baseLower[iRow] + primal_feasibility_tolerance;
      if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    } else if (alpha < -alphaTol) {
      relaxSpace =
	baseValue[iRow] - baseUpper[iRow] - primal_feasibility_tolerance;
      if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    }
  }
  analysis->simplexTimerStop(Chuzr1Clock);
  
  analysis->simplexTimerStart(Chuzr2Clock);
  double bestAlpha = 0;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    double alpha = col_aq.array[iRow] * move_in;
    if (alpha > alphaTol) {
      // Positive pivotal column entry
      double tightSpace = baseValue[iRow] - baseLower[iRow];
      if (tightSpace < relaxTheta * alpha) {
        if (bestAlpha < alpha) {
          bestAlpha = alpha;
          row_out = iRow;
        }
      }
    } else if (alpha < -alphaTol) {
      // Negative pivotal column entry
      double tightSpace = baseValue[iRow] - baseUpper[iRow];
      if (tightSpace > relaxTheta * alpha) {
        if (bestAlpha < -alpha) {
          bestAlpha = -alpha;
          row_out = iRow;
        }
      }
    }
  }
  analysis->simplexTimerStop(Chuzr2Clock);
}

void HEkkPrimal::updateDual() {
  analysis->simplexTimerStart(UpdateDualClock);
  assert(alpha_col);
  assert(row_out >= 0);
  vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  //  const vector<int>& nonbasicMove =
  //  ekk_instance_.simplex_basis_.nonbasicMove_;
  // Update the duals
  theta_dual = workDual[variable_in] / alpha_col;
  for (int iEl = 0; iEl < row_ap.count; iEl++) {
    int iCol = row_ap.index[iEl];
    workDual[iCol] -= theta_dual * row_ap.array[iCol];
  }
  for (int iEl = 0; iEl < row_ep.count; iEl++) {
    int iRow = row_ep.index[iEl];
    int iCol = iRow + num_col;
    workDual[iCol] -= theta_dual * row_ep.array[iRow];
  }
  // Dual for the pivot
  workDual[variable_in] = 0;
  workDual[variable_out] = -theta_dual;
  
  ekk_instance_.invalidateDualInfeasibilityRecord();
  // After dual update in primal simplex the dual objective value is not known
  ekk_instance_.simplex_lp_status_.has_dual_objective_value = false;
  analysis->simplexTimerStop(UpdateDualClock);
}

void HEkkPrimal::phase1ComputeDual() {
  const vector<double>& baseLower = ekk_instance_.simplex_info_.baseLower_;
  const vector<double>& baseUpper = ekk_instance_.simplex_info_.baseUpper_;
  const vector<double>& baseValue = ekk_instance_.simplex_info_.baseValue_;
  const vector<int>& nonbasicFlag = ekk_instance_.simplex_basis_.nonbasicFlag_;
  vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  
  vector<double>& workCost = ekk_instance_.simplex_info_.workCost_;
  // Accumulate costs for checking
  workCost.assign(num_tot, 0);
  
  HVector buffer;
  buffer.setup(num_row);
  buffer.clear();
  buffer.count = 0;
  for (int iRow = 0; iRow < num_row; iRow++) {
    double cost = 0;
    if (baseValue[iRow] < baseLower[iRow] - dual_feasibility_tolerance) {
      cost = -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dual_feasibility_tolerance) {
      cost = 1.0;
    }
    buffer.array[iRow] = cost;
    if (cost) buffer.index[buffer.count++] = iRow;
    workCost[ekk_instance_.simplex_basis_.basicIndex_[iRow]] = cost;
  }
  //
  // Full BTRAN
  //
  ekk_instance_.fullBtran(buffer);
  //
  // Full PRICE
  //
  HVector bufferLong;
  bufferLong.setup(num_col);
  ekk_instance_.fullPrice(buffer, bufferLong);
  
  for (int iCol = 0; iCol < num_col; iCol++)
    workDual[iCol] = -nonbasicFlag[iCol] * bufferLong.array[iCol];
  for (int iRow = 0, iCol = num_col; iRow < num_row; iRow++, iCol++)
    workDual[iCol] = -nonbasicFlag[iCol] * buffer.array[iRow];
  
  ekk_instance_.computeSimplexDualInfeasible();
}

void HEkkPrimal::phase1ChooseRow() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<double>& baseValue = simplex_info.baseValue_;
  analysis->simplexTimerStart(Chuzr1Clock);
  // Collect phase 1 theta lists
  //
  
  const double dPivotTol = simplex_info.update_count < 10
						       ? 1e-9
						       : simplex_info.update_count < 20 ? 1e-8 : 1e-7;
  ph1SorterR.clear();
  ph1SorterT.clear();
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    double dAlpha = col_aq.array[iRow] * move_in;
    
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
    row_out = -1;
    variable_out = -1;
    return;
  }

  // Now sort the relaxed theta to find the final break point. TODO:
  // Consider partial sort. Or heapify [O(n)] and then pop k points
  // [kO(log(n))].

  analysis->simplexTimerStart(Chuzr2Clock);
  std::sort(ph1SorterR.begin(), ph1SorterR.end());
  double dMaxTheta = ph1SorterR.at(0).first;
  double dGradient = fabs(theta_dual);
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
  row_out = -1;
  variable_out = -1;
  move_out = 0;
  for (int i = iLast - 1; i >= 0; i--) {
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + num_row;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    if (dAbsAlpha > dMaxAlpha * 0.1) {
      row_out = iRow;
      move_out = index >= 0 ? 1 : -1;
      break;
    }
  }
  analysis->simplexTimerStop(Chuzr2Clock);
}

void HEkkPrimal::phase1UpdatePrimal() {
  analysis->simplexTimerStart(UpdatePrimalClock);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<int>& basicIndex = ekk_instance_.simplex_basis_.basicIndex_;
  vector<double>& workDual = simplex_info.workDual_;
  vector<double>& workCost = simplex_info.workCost_;
  vector<double>& baseValue = simplex_info.baseValue_;
  col_basic_feasibility_change.clear();
  //
  // Update basic primal values, identifying all the feasibility
  // changes giving a value to col_basic_feasibility_change so that the duals
  // can be updated.
  //  if (ekk_instance_.sparseLoopStyle(col_aq.count, num_row, to_entry)) {
  for (int iEl = 0; iEl < col_aq.count; iEl++) {
    int iRow = col_aq.index[iEl];
    baseValue[iRow] -= theta_primal * col_aq.array[iRow];
    int iCol = basicIndex[iRow];
    double was_cost = workCost[iCol];
    double cost = 0;
    if (baseValue[iRow] < baseLower[iRow] - primal_feasibility_tolerance) {
      cost = -1.0;
    } else if (baseValue[iRow] >
               baseUpper[iRow] + primal_feasibility_tolerance) {
      cost = 1.0;
    }
    workCost[iCol] = cost;
    if (was_cost) {
      if (!cost) simplex_info.num_primal_infeasibilities--;
    } else {
      if (cost) simplex_info.num_primal_infeasibilities++;
    }
    double delta_cost = cost - was_cost;
    if (delta_cost) {
      col_basic_feasibility_change.array[iRow] = delta_cost;
      col_basic_feasibility_change.index[col_basic_feasibility_change.count++] =
          iRow;
      if (iCol >= num_col) workDual[iCol] += delta_cost;
    }
  }
  // Don't set baseValue[row_out] yet so that dual update due to
  // feasibility changes is done correctly
  ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
  analysis->simplexTimerStop(UpdatePrimalClock);
}

void HEkkPrimal::basicFeasibilityChangeUpdateDual() {
  analysis->simplexTimerStart(UpdateDualBasicFeasibilityChangeClock);
  vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  // For basic logicals, the change in the basic cost will be a
  // component in col_basic_feasibility_change. This will lead to it being
  // subtracted from workDual in the loop below over the
  // nonzeros in col_basic_feasibility_change, so add it in now. For basic
  // structurals, there will be no corresponding component in
  // row_basic_feasibility_change, since only the nonbasic components are
  // computed (avoided using row pricing, and basic components
  // zeroed after column pricing). Hence there will be no
  // subtraction in the loop below over the nonzeros in
  // row_basic_feasibility_change. Hence, only add in the basic cost change
  // for logicals.
  //
  // Assumes that row_basic_feasibility_change has been set up in
  // phase1UpdatePrimal()

  basicFeasibilityChangeBtran();
  basicFeasibilityChangePrice();
  int to_entry;
  const bool use_row_indices = ekk_instance_.sparseLoopStyle(
      row_basic_feasibility_change.count, num_col, to_entry);
  for (int iEntry = 0; iEntry < to_entry; iEntry++) {
    int iCol;
    if (use_row_indices) {
      iCol = row_basic_feasibility_change.index[iEntry];
    } else {
      iCol = iEntry;
    }
    workDual[iCol] -= row_basic_feasibility_change.array[iCol];
  }
  const bool use_col_indices = ekk_instance_.sparseLoopStyle(
      col_basic_feasibility_change.count, num_row, to_entry);
  for (int iEntry = 0; iEntry < to_entry; iEntry++) {
    int iRow;
    if (use_col_indices) {
      iRow = col_basic_feasibility_change.index[iEntry];
    } else {
      iRow = iEntry;
    }
    int iCol = num_col + iRow;
    workDual[iCol] -= col_basic_feasibility_change.array[iRow];
  }
  ekk_instance_.invalidateDualInfeasibilityRecord();
  analysis->simplexTimerStop(UpdateDualBasicFeasibilityChangeClock);
}

void HEkkPrimal::basicFeasibilityChangeBtran() {
  // Performs BTRAN on col_basic_feasibility_change. Make sure that
  // col_basic_feasibility_change.count is large (>simplex_lp_.numRow_ to be
  // sure) rather than 0 if the indices of the RHS (and true value of
  // col_basic_feasibility_change.count) isn't known.
  analysis->simplexTimerStart(BtranBasicFeasibilityChangeClock);
  const int solver_num_row = ekk_instance_.simplex_lp_.numRow_;
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  if (simplex_info.analyse_iterations)
    analysis->operationRecordBefore(
        ANALYSIS_OPERATION_TYPE_BTRAN_BASIC_FEASIBILITY_CHANGE,
        col_basic_feasibility_change,
        analysis->col_basic_feasibility_change_density);
#endif
  ekk_instance_.factor_.btran(col_basic_feasibility_change,
                              analysis->col_basic_feasibility_change_density,
                              analysis->pointer_serial_factor_clocks);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(
        ANALYSIS_OPERATION_TYPE_BTRAN_BASIC_FEASIBILITY_CHANGE,
        col_basic_feasibility_change);
#endif
  const double local_col_basic_feasibility_change_density =
      (double)col_basic_feasibility_change.count / solver_num_row;
  analysis->updateOperationResultDensity(
      local_col_basic_feasibility_change_density,
      analysis->col_basic_feasibility_change_density);
  analysis->simplexTimerStop(BtranBasicFeasibilityChangeClock);
}

void HEkkPrimal::basicFeasibilityChangePrice() {
  analysis->simplexTimerStart(PriceBasicFeasibilityChangeClock);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const int solver_num_row = ekk_instance_.simplex_lp_.numRow_;
  const int solver_num_col = ekk_instance_.simplex_lp_.numCol_;
  const double local_density =
      1.0 * col_basic_feasibility_change.count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  ekk_instance_.choosePriceTechnique(simplex_info.price_strategy, local_density,
                                     use_col_price, use_row_price_w_switch);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) {
    if (use_col_price) {
      const double historical_density_for_non_hypersparse_operation = 1;
      analysis->operationRecordBefore(
          ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE,
          col_basic_feasibility_change,
          historical_density_for_non_hypersparse_operation);
      analysis->num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis->operationRecordBefore(
          ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE,
          col_basic_feasibility_change,
          analysis->col_basic_feasibility_change_density);
      analysis->num_row_price_with_switch++;
    } else {
      analysis->operationRecordBefore(
          ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE,
          col_basic_feasibility_change,
          analysis->col_basic_feasibility_change_density);
      analysis->num_row_price++;
    }
  }
#endif
  row_basic_feasibility_change.clear();
  if (use_col_price) {
    // Perform column-wise PRICE
    ekk_instance_.matrix_.priceByColumn(row_basic_feasibility_change,
                                        col_basic_feasibility_change);
  } else if (use_row_price_w_switch) {
    // Perform hyper-sparse row-wise PRICE, but switch if the density of
    // row_basic_feasibility_change becomes extreme
    const double switch_density = ekk_instance_.matrix_.hyperPRICE;
    ekk_instance_.matrix_.priceByRowSparseResultWithSwitch(
        row_basic_feasibility_change, col_basic_feasibility_change,
        analysis->row_basic_feasibility_change_density, 0, switch_density);
  } else {
    // Perform hyper-sparse row-wise PRICE
    ekk_instance_.matrix_.priceByRowSparseResult(row_basic_feasibility_change,
                                                 col_basic_feasibility_change);
  }
  if (use_col_price) {
    // Column-wise PRICE computes components corresponding to basic
    // variables, so zero these by exploiting the fact that, for basic
    // variables, nonbasicFlag[*]=0
    const int* nonbasicFlag = &ekk_instance_.simplex_basis_.nonbasicFlag_[0];
    for (int iCol = 0; iCol < solver_num_col; iCol++)
      row_basic_feasibility_change.array[iCol] *= nonbasicFlag[iCol];
  }
  // Update the record of average row_basic_feasibility_change density
  const double local_row_basic_feasibility_change_density =
      (double)row_basic_feasibility_change.count / solver_num_col;
  analysis->updateOperationResultDensity(
      local_row_basic_feasibility_change_density,
      analysis->row_basic_feasibility_change_density);
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations)
    analysis->operationRecordAfter(
        ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE,
        row_basic_feasibility_change);
#endif
  analysis->simplexTimerStop(PriceBasicFeasibilityChangeClock);
}

void HEkkPrimal::phase2UpdatePrimal() {
  analysis->simplexTimerStart(UpdatePrimalClock);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<double>& workLower = simplex_info.workLower_;
  const vector<double>& workUpper = simplex_info.workUpper_;
  const vector<double>& workDual = simplex_info.workDual_;
  vector<double>& baseValue = simplex_info.baseValue_;

  bool primal_infeasible = false;
  //  if (ekk_instance_.sparseLoopStyle(col_aq.count, num_row, to_entry)) {
  for (int iEl = 0; iEl < col_aq.count; iEl++) {
    int iRow = col_aq.index[iEl];
    baseValue[iRow] -= theta_primal * col_aq.array[iRow];
    double lower = baseLower[iRow];
    double upper = baseUpper[iRow];
    double value = baseValue[iRow];
    double primal_infeasibility = 0;
    if (value < lower - primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > primal_feasibility_tolerance) {
      simplex_info.num_primal_infeasibilities++;
      primal_infeasible = true;
    }
  }
  if (row_out >= 0) {
    // Check the feasibility of the entering variable - using work
    // vector values since its base vector values haven't been updated
    baseValue[row_out] = value_in;
    double lower = workLower[variable_in];
    double upper = workUpper[variable_in];
    double value = value_in;
    double primal_infeasibility = 0;
    if (value < lower - primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > primal_feasibility_tolerance) {
      printf(
          "Entering varible has primal infeasibility of %g for [%g, %g, %g]\n",
          primal_infeasibility, workLower[variable_in], value_in,
          workUpper[variable_in]);
      simplex_info.num_primal_infeasibilities++;
      primal_infeasible = true;
    }
  }

  if (primal_infeasible)
    rebuild_reason = REBUILD_REASON_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;

  simplex_info.updated_primal_objective_value +=
      workDual[variable_in] * theta_primal;

  ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
  analysis->simplexTimerStop(UpdatePrimalClock);
}

void HEkkPrimal::considerInfeasibleValueIn() {
  assert(row_out >= 0);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& workLower = simplex_info.workLower_;
  const vector<double>& workUpper = simplex_info.workUpper_;
  vector<double>& workCost = simplex_info.workCost_;
  vector<double>& workDual = simplex_info.workDual_;
  double cost = 0;
  if (value_in < workLower[variable_in] - primal_feasibility_tolerance) {
    cost = -1.0;
  } else if (value_in > workUpper[variable_in] + primal_feasibility_tolerance) {
    cost = 1.0;
  }
  workCost[variable_in] = cost;
  if (cost) simplex_info.num_primal_infeasibilities++;
  workDual[variable_in] += cost;
}

void HEkkPrimal::resetDevex() {
  devex_weight.assign(num_tot, 1.0);
  devex_index.assign(num_tot, 0);
  for (int iCol = 0; iCol < num_tot; iCol++) {
    const int nonbasicFlag = ekk_instance_.simplex_basis_.nonbasicFlag_[iCol];
    devex_index[iCol] = nonbasicFlag * nonbasicFlag;
  }
  num_devex_iterations = 0;
  num_bad_devex_weight = 0;
  if (report_hyper_chuzc) printf("resetDevex\n");
  hyperChooseColumnClear();
}

void HEkkPrimal::updateDevex() {
  analysis->simplexTimerStart(DevexUpdateWeightClock);
  // Compute the pivot weight from the reference set
  double dPivotWeight = 0.0;
  int to_entry;
  const bool use_col_indices =
      ekk_instance_.sparseLoopStyle(col_aq.count, num_row, to_entry);
  for (int iEntry = 0; iEntry < to_entry; iEntry++) {
    int iRow;
    if (use_col_indices) {
      iRow = col_aq.index[iEntry];
    } else {
      iRow = iEntry;
    }
    int iCol = ekk_instance_.simplex_basis_.basicIndex_[iRow];
    double dAlpha = devex_index[iCol] * col_aq.array[iRow];
    dPivotWeight += dAlpha * dAlpha;
  }
  dPivotWeight += devex_index[variable_in] * 1.0;
  dPivotWeight = sqrt(dPivotWeight);

  // Check if the saved weight is too large
  if (devex_weight[variable_in] > bad_devex_weight_factor * dPivotWeight)
    num_bad_devex_weight++;

  // Update the devex weight for all
  double dPivot = col_aq.array[row_out];
  dPivotWeight /= fabs(dPivot);

  for (int iEl = 0; iEl < row_ap.count; iEl++) {
    int iCol = row_ap.index[iEl];
    double alpha = row_ap.array[iCol];
    double devex = dPivotWeight * fabs(alpha);
    devex += devex_index[iCol] * 1.0;
    if (devex_weight[iCol] < devex) {
      devex_weight[iCol] = devex;
    }
  }
  for (int iEl = 0; iEl < row_ep.count; iEl++) {
    int iRow = row_ep.index[iEl];
    int iCol = iRow + num_col;
    double alpha = row_ep.array[iRow];
    double devex = dPivotWeight * fabs(alpha);
    devex += devex_index[iCol] * 1.0;
    if (devex_weight[iCol] < devex) {
      devex_weight[iCol] = devex;
    }
  }

  // Update devex weight for the pivots
  devex_weight[variable_out] = max(1.0, dPivotWeight);
  devex_weight[variable_in] = 1.0;
  num_devex_iterations++;
  analysis->simplexTimerStop(DevexUpdateWeightClock);
}

void HEkkPrimal::updateVerify() {
  // updateVerify for primal
  const double numerical_trouble_tolerance = 1e-7;
  numericalTrouble = 0;
  double abs_alpha_from_col = fabs(alpha_col);
  std::string alpha_row_source;
  if (variable_in < num_col) {
    alpha_row = row_ap.array[variable_in];
    alpha_row_source = "Col";
  } else {
    alpha_row = row_ep.array[variable_in - num_col];
    alpha_row_source = "Row";
  }
  double abs_alpha_from_row = fabs(alpha_row);
  double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  double min_abs_alpha = min(abs_alpha_from_col, abs_alpha_from_row);
  numericalTrouble = abs_alpha_diff / min_abs_alpha;
  if (numericalTrouble > numerical_trouble_tolerance)
    printf(
        "Numerical check: Iter %4d: alpha_col = %12g, (From %3s alpha_row = "
        "%12g), aDiff = %12g: measure = %12g\n",
        ekk_instance_.iteration_count_, alpha_col, alpha_row_source.c_str(),
        alpha_row, abs_alpha_diff, numericalTrouble);
  assert(numericalTrouble < 1e-3);
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  //
  if (numericalTrouble > 1e-7 && ekk_instance_.simplex_info_.update_count > 0)
    rebuild_reason = REBUILD_REASON_POSSIBLY_SINGULAR_BASIS;
}

void HEkkPrimal::iterationAnalysisData() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  analysis->simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
  analysis->edge_weight_mode = DualEdgeWeightMode::DEVEX;
  analysis->solve_phase = solvePhase;
  analysis->simplex_iteration_count = ekk_instance_.iteration_count_;
  analysis->devex_iteration_count = num_devex_iterations;
  analysis->pivotal_row_index = row_out;
  analysis->leaving_variable = variable_out;
  analysis->entering_variable = variable_in;
  analysis->rebuild_reason = rebuild_reason;
  analysis->reduced_rhs_value = 0;
  analysis->reduced_cost_value = 0;
  analysis->edge_weight = 0;
  analysis->primal_delta = 0;
  analysis->primal_step = theta_primal;
  analysis->dual_step = theta_dual;
  analysis->pivot_value_from_column = alpha_col;
  analysis->pivot_value_from_row = alpha_row;
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

void HEkkPrimal::localReportIterHeader() {
  printf(" Iter ColIn Row_Out ColOut\n");
}

void HEkkPrimal::localReportIter(const bool header) {
  if (!report_hyper_chuzc) return;
  static int last_header_iteration_count;
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  int iteration_count = ekk_instance_.iteration_count_;
  if (header) {
    localReportIterHeader();
    last_header_iteration_count = iteration_count;
  } else {
    if (ekk_instance_.iteration_count_ > last_header_iteration_count + 10) {
      localReportIterHeader();
      last_header_iteration_count = iteration_count;
    }
    if (row_out >= 0) {
      printf("%5d %5d  %5d  %5d", iteration_count, variable_in, row_out, variable_out);
    } else {
      printf("%5d %5d Bound flip   ", iteration_count, variable_in);
    }
    if (check_column >= 0 && iteration_count >= check_iter) {
      int flag = ekk_instance_.simplex_basis_.nonbasicFlag_[check_column];
      int move = ekk_instance_.simplex_basis_.nonbasicMove_[check_column];
      double lower = simplex_info.workLower_[check_column];
      double upper = simplex_info.workUpper_[check_column];
      double value;
      if (flag == NONBASIC_FLAG_TRUE) {
        value = simplex_info.workValue_[check_column];
      } else {
        int iRow;
        for (iRow = 0; iRow < num_row; iRow++) {
          if (ekk_instance_.simplex_basis_.basicIndex_[iRow] == check_column)
            break;
        }
        assert(iRow < num_row);
        value = simplex_info.baseValue_[iRow];
      }
      printf(": Var %2d (%1d, %2d) [%9.4g, %9.4g, %9.4g]", check_column, flag,
             move, lower, value, upper);
      if (flag == NONBASIC_FLAG_TRUE) {
        double dual = simplex_info.workDual_[check_column];
        double weight = devex_weight[check_column];
        double infeasibility = -move * dual;
        if (lower == -HIGHS_CONST_INF && upper == HIGHS_CONST_INF)
          infeasibility = fabs(dual);
        if (infeasibility < dual_feasibility_tolerance) infeasibility = 0;
        double measure = infeasibility / weight;
        printf(" Du = %9.4g; Wt = %9.4g; Ms = %9.4g", dual, weight, measure);
      }
    }
    printf("\n");
  }
}

void HEkkPrimal::reportRebuild(const int reason_for_rebuild) {
  analysis->simplexTimerStart(ReportRebuildClock);
  iterationAnalysisData();
  analysis->rebuild_reason = reason_for_rebuild;
  analysis->invertReport();
  analysis->simplexTimerStop(ReportRebuildClock);
}

void HEkkPrimal::getNonbasicFreeColumnSet() {
  if (!num_free_col) return;
  assert(num_free_col > 0);
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const SimplexBasis& simplex_basis = ekk_instance_.simplex_basis_;
  nonbasic_free_col_set.clear();
  for (int iCol = 0; iCol < num_tot; iCol++) {
    bool nonbasic_free =
        simplex_basis.nonbasicFlag_[iCol] == NONBASIC_FLAG_TRUE &&
        simplex_info.workLower_[iCol] <= -HIGHS_CONST_INF &&
        simplex_info.workUpper_[iCol] >= HIGHS_CONST_INF;
    if (nonbasic_free) nonbasic_free_col_set.add(iCol);
  }
  nonbasic_free_col_set.print();
}

void HEkkPrimal::removeNonbasicFreeColumn() {
  bool remove_nonbasic_free_column =
      ekk_instance_.simplex_basis_.nonbasicMove_[variable_in] == 0;
  if (remove_nonbasic_free_column) {
    bool removed_nonbasic_free_column = nonbasic_free_col_set.remove(variable_in);
    if (!removed_nonbasic_free_column) {
      HighsLogMessage(
          ekk_instance_.options_.logfile, HighsMessageType::ERROR,
          "HEkkPrimal::phase1update failed to remove nonbasic free column %d",
          variable_in);
      assert(removed_nonbasic_free_column);
    }
  }
}
void HEkkPrimal::getBasicPrimalInfeasibility() {
  // Gets the num/max/sum of basic primal infeasibliities,
  analysis->simplexTimerStart(ComputePrIfsClock);
  const double primal_feasibility_tolerance =
      ekk_instance_.options_.primal_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<double>& baseValue = simplex_info.baseValue_;
  int& num_primal_infeasibilities = simplex_info.num_primal_infeasibilities;
  double& max_primal_infeasibility = simplex_info.max_primal_infeasibility;
  double& sum_primal_infeasibilities = simplex_info.sum_primal_infeasibilities;
  const int updated_num_primal_infeasibilities = num_primal_infeasibilities;
  num_primal_infeasibilities = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibilities = 0;

  for (int iRow = 0; iRow < num_row; iRow++) {
    double value = baseValue[iRow];
    double lower = baseLower[iRow];
    double upper = baseUpper[iRow];
    double primal_infeasibility = 0;
    if (value < lower - primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > primal_feasibility_tolerance)
        num_primal_infeasibilities++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibilities += primal_infeasibility;
    }
  }
  if (updated_num_primal_infeasibilities >= 0) {
    // The number of primal infeasibliities should be correct
    bool num_primal_infeasibilities_ok =
        num_primal_infeasibilities == updated_num_primal_infeasibilities;
    if (!num_primal_infeasibilities_ok) {
      printf(
          "In iteration %d: num_primal_infeasibilities = %d != %d = "
          "updated_num_primal_infeasibilities\n",
          ekk_instance_.iteration_count_, num_primal_infeasibilities,
          updated_num_primal_infeasibilities);
      assert(num_primal_infeasibilities_ok);
    }
  }
  analysis->simplexTimerStop(ComputePrIfsClock);
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
