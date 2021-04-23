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
/**@file simplex/HEkkPrimal.cpp
 * @brief
 */
#include "simplex/HEkkPrimal.h"

#include "simplex/HEkkDebug.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsSort.h"

HighsStatus HEkkPrimal::solve() {
  HighsOptions& options = ekk_instance_.options_;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // Assumes that the LP has a positive number of rows
  bool positive_num_row = ekk_instance_.simplex_lp_.numRow_ > 0;
  if (!positive_num_row) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "HEkkPrimal::solve called for LP with non-positive (%" HIGHSINT_FORMAT
        ") "
        "number of constraints\n",
        ekk_instance_.simplex_lp_.numRow_);
    assert(positive_num_row);
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  }
  if (ekk_instance_.bailoutOnTimeIterations())
    return ekk_instance_.returnFromSolve(HighsStatus::kWarning);

  if (!simplex_lp_status.has_invert) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HEkkPrimal::solve called without INVERT\n");
    assert(simplex_lp_status.has_fresh_invert);
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  }

  if (debugPrimalSimplex("Initialise", true) == HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);

  // Get the nonabsic free column set
  getNonbasicFreeColumnSet();

  // Determine whether the solution is near-optimal. Value 1 is
  // unimportant, as the sum of dual infeasiblilities for near-optimal
  // solutions is typically many orders of magnitude smaller than 1,
  // and the sum of dual infeasiblilities will be very much larger for
  // non-trivial LPs that are primal feasible for a logical or crash
  // basis.
  const bool near_optimal = simplex_info.num_primal_infeasibility == 0 &&
                            simplex_info.sum_dual_infeasibility < 1;
  if (near_optimal)
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "Primal feasible and num / max / sum dual infeasibilities are "
                "%" HIGHSINT_FORMAT
                " / %g "
                "/ %g, so near-optimal\n",
                simplex_info.num_dual_infeasibility,
                simplex_info.max_dual_infeasibility,
                simplex_info.sum_dual_infeasibility);

  // Perturb bounds according to whether the solution is near-optimnal
  const bool perturb_bounds = !near_optimal;
  if (!perturb_bounds)
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "Near-optimal, so don't use bound perturbation\n");
  if (perturb_bounds &&
      simplex_info.primal_simplex_bound_perturbation_multiplier) {
    ekk_instance_.initialiseBound(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown,
                                  perturb_bounds);
    ekk_instance_.initialiseNonbasicValueAndMove();
    ekk_instance_.computePrimal();
    ekk_instance_.computeSimplexPrimalInfeasible();
  }
  HighsInt num_primal_infeasibility =
      ekk_instance_.simplex_info_.num_primal_infeasibility;
  solvePhase = num_primal_infeasibility > 0 ? kSolvePhase1 : kSolvePhase2;

  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         ekk_instance_.scaled_model_status_) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);

  // The major solving loop
  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  // ToDo Move to simplex and adapt so it's OK for primal and dual
  //  iterationAnalysisInitialise();

  localReportIter(true);
  correctPrimal(true);
  while (solvePhase) {
    HighsInt it0 = ekk_instance_.iteration_count_;
    // When starting a new phase the (updated) primal objective function
    // value isn't known. Indicate this so that when the value
    // computed from scratch in rebuild() isn't checked against the the
    // updated value
    simplex_lp_status.has_primal_objective_value = false;
    if (solvePhase == kSolvePhaseUnknown) {
      // Determine the number of primal infeasibilities, and hence the solve
      // phase
      ekk_instance_.computeSimplexPrimalInfeasible();
      num_primal_infeasibility =
          ekk_instance_.simplex_info_.num_primal_infeasibility;
      solvePhase = num_primal_infeasibility > 0 ? kSolvePhase1 : kSolvePhase2;
      if (simplex_info.backtracking_) {
        // Backtracking
        ekk_instance_.initialiseCost(SimplexAlgorithm::kPrimal, solvePhase);
        ekk_instance_.initialiseNonbasicValueAndMove();
        // Can now forget that we might have been backtracking
        simplex_info.backtracking_ = false;
      }
    }
    assert(solvePhase == kSolvePhase1 || solvePhase == kSolvePhase2);
    if (solvePhase == kSolvePhase1) {
      //
      // Phase 1
      //
      // solvePhase = kSolvePhase1 if the iteration or time limit has
      // been reached
      //
      // solvePhase = kSolvePhase2 if there are no primal infeasibilities
      //
      // solvePhase = kSolvePhaseUnknown if backtracking
      //
      // solvePhase = kSolvePhaseExit if primal infeasiblilty is
      // detected, in which case scaled_model_status_ =
      // HighsModelStatus::kInfeasible is set
      //
      // solvePhase = kSolvePhaseError is set if an error occurs
      solvePhase1();
      assert(solvePhase == kSolvePhase1 || solvePhase == kSolvePhase2 ||
             solvePhase == kSolvePhaseUnknown ||
             solvePhase == kSolvePhaseExit || solvePhase == kSolvePhaseError);
      simplex_info.primal_phase1_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else if (solvePhase == kSolvePhase2) {
      //
      // Phase 2
      //
      // solvePhase = kSolvePhaseOptimal if there are no dual
      // infeasibilities
      //
      // solvePhase = kSolvePhase1 if there are primal
      // infeasibilities
      //
      // solvePhase = kSolvePhase2 if the iteration or time limit has
      // been reached
      //
      // solvePhase = kSolvePhaseCleanup if there are primal
      // infeasiblilities to clean up after removing bound shifts
      //
      // solvePhase = kSolvePhaseUnknown if backtracking
      //
      // solvePhase = kSolvePhaseExit if primal unboundedness is
      // detected, in which case scaled_model_status_ =
      // HighsModelStatus::kUnbounded is set
      //
      // solvePhase = kSolvePhaseError is set if an error occurs
      solvePhase2();
      assert(solvePhase == kSolvePhaseOptimal || solvePhase == kSolvePhase1 ||
             solvePhase == kSolvePhase2 || solvePhase == kSolvePhaseCleanup ||
             solvePhase == kSolvePhaseUnknown ||
             solvePhase == kSolvePhaseExit || solvePhase == kSolvePhaseError);
      assert(solvePhase != kSolvePhaseExit ||
             ekk_instance_.scaled_model_status_ ==
                 HighsModelStatus::kUnbounded);
      simplex_info.primal_phase2_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else {
      // Should only be kSolvePhase1 or kSolvePhase2
      ekk_instance_.scaled_model_status_ = HighsModelStatus::kSolveError;
      return ekk_instance_.returnFromSolve(HighsStatus::kError);
    }
    // Return if bailing out from solve
    if (ekk_instance_.solve_bailout_)
      return ekk_instance_.returnFromSolve(HighsStatus::kWarning);
    // Can have all possible cases of solvePhase
    assert(solvePhase >= kSolvePhaseMin && solvePhase <= kSolvePhaseMax);
    // Look for scenarios when the major solving loop ends
    if (solvePhase == kSolvePhaseError) {
      // Solver error so return HighsStatus::kError
      ekk_instance_.scaled_model_status_ = HighsModelStatus::kSolveError;
      return ekk_instance_.returnFromSolve(HighsStatus::kError);
    }
    if (solvePhase == kSolvePhaseExit) {
      // LP identified as not having an optimal solution
      assert(
          ekk_instance_.scaled_model_status_ == HighsModelStatus::kInfeasible ||
          ekk_instance_.scaled_model_status_ == HighsModelStatus::kUnbounded);
      break;
    }
    if (solvePhase == kSolvePhaseCleanup) {
      // Primal infeasibilities after phase 2 for a problem not known
      // to be primal infeasible. Dual feasible with primal
      // infeasibilities so use dual simplex to clean up
      break;
    }
    // If solvePhase == kSolvePhaseOptimal == 0 then major solving
    // loop ends naturally since solvePhase is false
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  // Should only have these cases
  assert(solvePhase == kSolvePhaseExit || solvePhase == kSolvePhaseUnknown ||
         solvePhase == kSolvePhaseOptimal || solvePhase == kSolvePhase1 ||
         solvePhase == kSolvePhaseCleanup);
  if (solvePhase == kSolvePhaseOptimal)
    ekk_instance_.scaled_model_status_ = HighsModelStatus::kOptimal;
  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         ekk_instance_.scaled_model_status_) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::kError);
  return ekk_instance_.returnFromSolve(HighsStatus::kOk);
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

  rebuild_reason = kRebuildReasonNo;

  ekk_instance_.simplex_lp_status_.has_primal_objective_value = false;
  ekk_instance_.simplex_lp_status_.has_dual_objective_value = false;
  ekk_instance_.scaled_model_status_ = HighsModelStatus::kNotset;
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
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (ekk_instance_.simplex_info_.workLower_[iCol] == -kHighsInf &&
        ekk_instance_.simplex_info_.workUpper_[iCol] == kHighsInf) {
      // Free column
      num_free_col++;
    }
  }
  // Set up the HSet instances, possibly using the internal error reporting and
  // debug option
  const bool debug =
      ekk_instance_.options_.highs_debug_level > kHighsDebugLevelCheap;
  if (num_free_col) {
    highsLogUser(ekk_instance_.options_.log_options, HighsLogType::kInfo,
                 "HEkkPrimal:: LP has %" HIGHSINT_FORMAT " free columns\n",
                 num_free_col);
    nonbasic_free_col_set.setup(num_free_col, num_tot,
                                ekk_instance_.options_.output_flag,
                                ekk_instance_.options_.log_file_stream, debug);
  }
  // Set up the hyper-sparse CHUZC data
  hyper_chuzc_candidate.resize(1 + max_num_hyper_chuzc_candidates);
  hyper_chuzc_measure.resize(1 + max_num_hyper_chuzc_candidates);
  hyper_chuzc_candidate_set.setup(max_num_hyper_chuzc_candidates, num_tot,
                                  ekk_instance_.options_.output_flag,
                                  ekk_instance_.options_.log_file_stream,
                                  debug);
}

void HEkkPrimal::solvePhase1() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailoutReturn()) return;
  highsLogDev(ekk_instance_.options_.log_options, HighsLogType::kDetailed,
              "primal-phase1-start\n");
  // If there's no backtracking basis, save the initial basis in case of
  // backtracking
  if (!simplex_info.valid_backtracking_basis_)
    ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  for (;;) {
    //
    // Rebuild
    //
    // solvePhase = kSolvePhaseError is set if the basis matrix is singular
    rebuild();
    if (solvePhase == kSolvePhaseError) return;
    if (solvePhase == kSolvePhaseUnknown) return;
    if (ekk_instance_.bailoutOnTimeIterations()) return;
    assert(solvePhase == kSolvePhase1 || solvePhase == kSolvePhase2);
    //
    // solvePhase = kSolvePhase2 is set if no primal infeasibilities
    // are found in rebuild(), in which case return for phase 2
    if (solvePhase == kSolvePhase2) break;

    for (;;) {
      iterate();
      if (ekk_instance_.bailoutOnTimeIterations()) return;
      if (solvePhase == kSolvePhaseError) return;
      assert(solvePhase == kSolvePhase1);
      if (rebuild_reason) break;
    }
    // If the data are fresh from rebuild() and no flips have
    // occurred, break out of the outer loop to see what's ocurred
    if (simplex_lp_status.has_fresh_rebuild && num_flip_since_rebuild == 0)
      break;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  // Will only have accurate simplex info if moving to phase 2 - but
  // should check primal feasiblilty and residual information if LP
  // is primal infeasible
  if (debugPrimalSimplex("End of solvePhase1") ==
      HighsDebugStatus::kLogicalError) {
    solvePhase = kSolvePhaseError;
    return;
  }
  // Possible to have switched to phase 2 without any iterations - if
  // LP with inconsistent bounds is being solved. Not that this should
  // be allowed to happen! ToDo add inconsistent bounds check in run()
  if (solvePhase == kSolvePhase1) {
    // Determine whether primal infeasiblility has been identified
    if (variable_in < 0) {
      // Optimal in phase 1, so should have primal infeasiblilities
      assert(simplex_info.num_primal_infeasibility > 0);
      ekk_instance_.scaled_model_status_ = HighsModelStatus::kInfeasible;
      solvePhase = kSolvePhaseExit;
    }
  }
  if (solvePhase == kSolvePhase2) {
    // Moving to phase 2 so comment if bound perturbation is not permitted
    //
    // It may have been prevented to avoid cleanup-perturbation loops
    if (!simplex_info.allow_bound_perturbation)
      highsLogDev(ekk_instance_.options_.log_options, HighsLogType::kWarning,
                  "Moving to phase 2, but not allowing bound perturbation\n");
  }
}

void HEkkPrimal::solvePhase2() {
  HighsOptions& options = ekk_instance_.options_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  HighsModelStatus& scaled_model_status = ekk_instance_.scaled_model_status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailoutReturn()) return;
  highsLogDev(options.log_options, HighsLogType::kDetailed,
              "primal-phase2-start\n");
  phase2UpdatePrimal(true);

  // If there's no backtracking basis Save the initial basis in case of
  // backtracking
  if (!ekk_instance_.simplex_info_.valid_backtracking_basis_)
    ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  for (;;) {
    //
    // Rebuild
    //
    // solvePhase = kSolvePhaseError is set if the basis matrix is singular
    rebuild();
    if (solvePhase == kSolvePhaseError) return;
    if (solvePhase == kSolvePhaseUnknown) return;
    if (ekk_instance_.bailoutOnTimeIterations()) return;
    assert(solvePhase == kSolvePhase1 || solvePhase == kSolvePhase2);
    //
    // solvePhase = kSolvePhase1 is set if primal infeasibilities
    // are found in rebuild(), in which case return for phase 1
    if (solvePhase == kSolvePhase1) break;

    for (;;) {
      iterate();
      if (ekk_instance_.bailoutOnTimeIterations()) return;
      if (solvePhase == kSolvePhaseError) return;
      assert(solvePhase == kSolvePhase2);
      if (rebuild_reason) break;
    }
    // If the data are fresh from rebuild() and no flips have
    // occurred, break out of the outer loop to see what's ocurred
    if (simplex_lp_status.has_fresh_rebuild && num_flip_since_rebuild == 0)
      break;
  }
  // If bailing out, should have returned already
  assert(!ekk_instance_.solve_bailout_);
  if (debugPrimalSimplex("End of solvePhase2") ==
      HighsDebugStatus::kLogicalError) {
    solvePhase = kSolvePhaseError;
    return;
  }
  if (solvePhase == kSolvePhase1) {
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "primal-return-phase1\n");
  } else if (variable_in == -1) {
    // There is no candidate in CHUZC, even after rebuild so probably optimal
    highsLogDev(options.log_options, HighsLogType::kDetailed,
                "primal-phase-2-optimal\n");
    // Remove any bound perturbations and see if basis is still primal feasible
    cleanup();
    if (ekk_instance_.simplex_info_.num_primal_infeasibility > 0) {
      // There are primal infeasiblities, so consider performing dual
      // simplex iterations to get primal feasibility
      solvePhase = kSolvePhaseCleanup;
    } else {
      // There are no dual infeasiblities so optimal!
      solvePhase = kSolvePhaseOptimal;
      highsLogDev(options.log_options, HighsLogType::kDetailed,
                  "problem-optimal\n");
      scaled_model_status = HighsModelStatus::kOptimal;
      ekk_instance_.computeDualObjectiveValue();  // Why?
    }
  } else {
    assert(row_out < 0);

    // There is no candidate in CHUZR, so probably primal unbounded
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "primal-phase-2-unbounded\n");
    if (ekk_instance_.simplex_info_.bounds_perturbed) {
      // If the bounds have been perturbed, clean up and return
      cleanup();
      // If there are primal infeasiblities, go back to phase 1
      if (ekk_instance_.simplex_info_.num_primal_infeasibility > 0)
        solvePhase = kSolvePhase1;
    } else {
      // The bounds have not been perturbed, so primal unbounded
      solvePhase = kSolvePhaseExit;
      // Primal unbounded, so save primal ray
      savePrimalRay();
      // Model status should be unset
      assert(scaled_model_status == HighsModelStatus::kNotset);
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "problem-primal-unbounded\n");
      scaled_model_status = HighsModelStatus::kUnbounded;
    }
  }
}

void HEkkPrimal::cleanup() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  if (!simplex_info.bounds_perturbed) return;
  highsLogDev(ekk_instance_.options_.log_options, HighsLogType::kDetailed,
              "primal-cleanup-shift\n");
  // Remove perturbation and don't permit further perturbation
  ekk_instance_.initialiseBound(SimplexAlgorithm::kPrimal, solvePhase, false);
  ekk_instance_.initialiseNonbasicValueAndMove();
  simplex_info.allow_bound_perturbation = false;
  // Possibly take a copy of the original duals before recomputing them
  /*
  vector<double> original_baseValue;
  if (ekk_instance_.options_.highs_debug_level > kHighsDebugLevelCheap)
    original_baseValue = simplex_info.baseValue_;
  */
  // Compute the primal values
  ekk_instance_.computePrimal();
  // Possibly analyse the change in duals
  /*  debugCleanup(ekk_instance_, original_baseValue); */
  // Compute the primal infeasibilities
  ekk_instance_.computeSimplexPrimalInfeasible();

  // Compute the primal objective value
  ekk_instance_.computePrimalObjectiveValue();
  // Now that there's a new primal_objective_value, reset the updated
  // value
  simplex_info.updated_primal_objective_value =
      simplex_info.primal_objective_value;

  //  if (!simplex_info.run_quiet) {
  // Report the dual infeasiblities
  ekk_instance_.computeSimplexDualInfeasible();
  // In phase 1, report the simplex LP dual infeasiblities
  // In phase 2, report the simplex dual infeasiblities (known)
  //    if (solvePhase == kSolvePhase1)
  //    computeSimplexLpDualInfeasible(ekk_instance_);
  reportRebuild();
  //  }
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
  HighsInt reason_for_rebuild = rebuild_reason;
  rebuild_reason = kRebuildReasonNo;
  // Possibly Rebuild factor
  bool reInvert = simplex_info.update_count > 0;
  if (reInvert) {
    // Get a nonsingular inverse if possible. One of three things
    // happens: Current basis is nonsingular; Current basis is
    // singular and last nonsingular basis is refactorized as
    // nonsingular - or found singular. Latter is code failure.
    if (!ekk_instance_.getNonsingularInverse(solvePhase)) {
      solvePhase = kSolvePhaseError;
      return;
    }
  }
  if (!ekk_instance_.simplex_lp_status_.has_matrix) {
    // Don't have the matrix either row-wise or col-wise, so
    // reinitialise it
    assert(simplex_info.backtracking_);
    HighsLp& simplex_lp = ekk_instance_.simplex_lp_;
    analysis->simplexTimerStart(matrixSetupClock);
    ekk_instance_.matrix_.setup(simplex_lp.numCol_, simplex_lp.numRow_,
                                &simplex_lp.Astart_[0], &simplex_lp.Aindex_[0],
                                &simplex_lp.Avalue_[0],
                                &ekk_instance_.simplex_basis_.nonbasicFlag_[0]);
    simplex_lp_status.has_matrix = true;
    analysis->simplexTimerStop(matrixSetupClock);
  }

  if (simplex_info.backtracking_) {
    // If backtracking, may change phase, so drop out
    solvePhase = kSolvePhaseUnknown;
    return;
  }

  ekk_instance_.computePrimal();
  if (solvePhase == kSolvePhase2) {
    bool correct_primal_ok = correctPrimal();
    assert(correct_primal_ok);
  }
  getBasicPrimalInfeasibility();
  if (simplex_info.num_primal_infeasibility > 0) {
    // Primal infeasibilities so should be in phase 1
    if (solvePhase == kSolvePhase2) {
      highsLogUser(
          ekk_instance_.options_.log_options, HighsLogType::kWarning,
          "HEkkPrimal::rebuild switching back to phase 1 from phase 2\n");
      solvePhase = kSolvePhase1;
    }
    phase1ComputeDual();
  } else {
    // No primal infeasibilities so in phase 2. Reset costs if was
    // previously in phase 1
    if (solvePhase == kSolvePhase1) {
      ekk_instance_.initialiseCost(SimplexAlgorithm::kPrimal, solvePhase);
      solvePhase = kSolvePhase2;
    }
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

  reportRebuild(reason_for_rebuild);

  // Record the synthetic clock for INVERT, and zero it for UPDATE
  ekk_instance_.build_syntheticTick_ =
      ekk_instance_.factor_.build_syntheticTick;
  ekk_instance_.total_syntheticTick_ = 0;

  // Determine whether to use hyper-sparse CHUZC
  if (solvePhase == kSolvePhase1) {
    use_hyper_chuzc = false;
  } else {
    use_hyper_chuzc = true;
  }
  hyperChooseColumnClear();

  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
  assert(solvePhase == kSolvePhase1 || solvePhase == kSolvePhase2);
}

void HEkkPrimal::iterate() {
  bool check = ekk_instance_.iteration_count_ >= check_iter;
  if (check) {
    printf("Iter %" HIGHSINT_FORMAT "\n", ekk_instance_.iteration_count_);
    ekk_instance_.options_.highs_debug_level = kHighsDebugLevelExpensive;
  }
  if (debugPrimalSimplex("Before iteration") ==
      HighsDebugStatus::kLogicalError) {
    solvePhase = kSolvePhaseError;
    return;
  }

  // Perform CHUZC
  //
  chuzc();
  if (variable_in == -1) {
    rebuild_reason = kRebuildReasonPossiblyOptimal;
    return;
  }

  // Perform FTRAN - and dual value cross-check to decide whether to use the
  // variable
  //
  // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
  // numerical trouble is detected
  if (!useVariableIn()) {
    if (rebuild_reason)
      assert(rebuild_reason == kRebuildReasonPossiblySingularBasis);
    return;
  }
  assert(!rebuild_reason);

  // Perform CHUZR
  if (solvePhase == kSolvePhase1) {
    phase1ChooseRow();
    if (row_out < 0) {
      highsLogUser(ekk_instance_.options_.log_options, HighsLogType::kError,
                   "Primal phase 1 choose row failed\n");
      solvePhase = kSolvePhaseError;
      return;
    }
  } else {
    chooseRow();
  }
  assert(!rebuild_reason);

  // Consider whether to perform a bound swap - either because it's
  // shorter than the pivoting step or, in the case of Phase 1,
  // because it's cheaper than pivoting - which may be questionable
  //
  // rebuild_reason = kRebuildReasonPossiblyPrimalUnbounded is set
  // in phase 2 if there's no pivot or bound swap. In phase 1 there is
  // always a pivot at this stage since row_out < 0 is trapped (above)
  // as an error.
  considerBoundSwap();
  if (rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded) return;
  assert(!rebuild_reason);

  if (row_out >= 0) {
    // Perform unit BTRAN and PRICE to get pivotal row - and do a
    // numerical check.
    //
    // rebuild_reason = kRebuildReasonPossiblySingularBasis is set
    // if numerical trouble is detected
    assessPivot();
    if (rebuild_reason) {
      assert(rebuild_reason == kRebuildReasonPossiblySingularBasis);
      return;
    }
  }
  // Any pivoting is numerically acceptable, so perform update.
  //
  // rebuild_reason =
  // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set if a
  // primal infeasiblility is found in phase 2
  //
  // rebuild_reason = kRebuildReasonUpdateLimitReached is set in
  // phase 1 if the number of primal infeasiblilities is reduced to
  // zero, or in either phase if the update count reaches the limit!
  //
  // rebuild_reason = kRebuildReasonSyntheticClockSaysInvert is
  // set in updateFactor() if it is considered to be more efficient to
  // reinvert.
  update();
  // Crude way to force rebuild if there are no infeasibilities in phase 1
  if (!ekk_instance_.simplex_info_.num_primal_infeasibility &&
      solvePhase == kSolvePhase1)
    rebuild_reason = kRebuildReasonUpdateLimitReached;

  assert(rebuild_reason == kRebuildReasonNo ||
         rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex ||
         rebuild_reason == kRebuildReasonSyntheticClockSaysInvert ||
         rebuild_reason == kRebuildReasonUpdateLimitReached);
  assert(solvePhase == kSolvePhase1 || solvePhase == kSolvePhase2);
}

void HEkkPrimal::chuzc() {
  if (done_next_chuzc) assert(use_hyper_chuzc);
  if (use_hyper_chuzc) {
    // Perform hyper-sparse CHUZC and then check result using full CHUZC
    if (!done_next_chuzc) chooseColumn(true);
    const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
    const bool check_hyper_chuzc = true;
    if (check_hyper_chuzc) {
      HighsInt hyper_sparse_variable_in = variable_in;
      chooseColumn(false);
      double hyper_sparse_measure = 0;
      if (hyper_sparse_variable_in >= 0)
        hyper_sparse_measure = fabs(workDual[hyper_sparse_variable_in]) /
                               devex_weight[hyper_sparse_variable_in];
      double measure = 0;
      if (variable_in >= 0)
        measure = fabs(workDual[variable_in]) / devex_weight[variable_in];
      double abs_measure_error = fabs(hyper_sparse_measure - measure);
      bool measure_error = abs_measure_error > 1e-12;
      if (measure_error) {
        printf("Iteration %" HIGHSINT_FORMAT
               ": Hyper-sparse CHUZC measure %g != %g = Full "
               "CHUZC measure (%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
               "): error %g\n",
               ekk_instance_.iteration_count_, hyper_sparse_measure, measure,
               hyper_sparse_variable_in, variable_in, abs_measure_error);
        assert(!measure_error);
      }
      variable_in = hyper_sparse_variable_in;
    }
  } else {
    chooseColumn(false);
  }
}

void HEkkPrimal::chooseColumn(const bool hyper_sparse) {
  assert(!hyper_sparse || !done_next_chuzc);
  const vector<int8_t>& nonbasicMove =
      ekk_instance_.simplex_basis_.nonbasicMove_;
  const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  double best_measure = 0;
  variable_in = -1;

  const bool local_use_hyper_chuzc = hyper_sparse;
  // Consider nonbasic free columns first
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (local_use_hyper_chuzc) {
    if (!initialise_hyper_chuzc) hyperChooseColumn();
    if (initialise_hyper_chuzc) {
      analysis->simplexTimerStart(ChuzcHyperInitialiselClock);
      num_hyper_chuzc_candidates = 0;
      if (num_nonbasic_free_col) {
        const vector<HighsInt>& nonbasic_free_col_set_entry =
            nonbasic_free_col_set.entry();
        for (HighsInt ix = 0; ix < num_nonbasic_free_col; ix++) {
          HighsInt iCol = nonbasic_free_col_set_entry[ix];
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
      for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
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
              "Full CHUZC: Max         measure is %9.4g for column "
              "%4" HIGHSINT_FORMAT
              ", and "
              "max non-candiate measure of  %9.4g\n",
              best_measure, variable_in, max_hyper_chuzc_non_candidate_measure);
      }
    }
  } else {
    analysis->simplexTimerStart(ChuzcPrimalClock);
    // Choose any attractive nonbasic free column
    if (num_nonbasic_free_col) {
      const vector<HighsInt>& nonbasic_free_col_set_entry =
          nonbasic_free_col_set.entry();
      for (HighsInt ix = 0; ix < num_nonbasic_free_col; ix++) {
        HighsInt iCol = nonbasic_free_col_set_entry[ix];
        double dual_infeasibility = fabs(workDual[iCol]);
        if (dual_infeasibility > dual_feasibility_tolerance &&
            dual_infeasibility > best_measure * devex_weight[iCol]) {
          variable_in = iCol;
          best_measure = dual_infeasibility / devex_weight[iCol];
        }
      }
    }
    // Now look at other columns
    for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
      double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
      if (dual_infeasibility > dual_feasibility_tolerance &&
          dual_infeasibility > best_measure * devex_weight[iCol]) {
        variable_in = iCol;
        best_measure = dual_infeasibility / devex_weight[iCol];
      }
    }
    analysis->simplexTimerStop(ChuzcPrimalClock);
  }
  //  printf("ChooseColumn: Iteration %" HIGHSINT_FORMAT ", choose column %"
  //  HIGHSINT_FORMAT " with measure %g\n",
  //	 ekk_instance_.iteration_count_, variable_in, best_measure);
}

bool HEkkPrimal::useVariableIn() {
  // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
  // numerical trouble is detected
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  vector<double>& workDual = simplex_info.workDual_;
  const vector<int8_t>& nonbasicMove =
      ekk_instance_.simplex_basis_.nonbasicMove_;
  const double updated_theta_dual = workDual[variable_in];
  // Determine the move direction - can't use nonbasicMove_[variable_in]
  // due to free columns
  move_in = updated_theta_dual > 0 ? -1 : 1;
  // Unless the variable is free, nonbasicMove[variable_in] should be the same
  // as move_in
  if (nonbasicMove[variable_in]) assert(nonbasicMove[variable_in] == move_in);
  //
  // FTRAN
  //
  // Compute pivot column
  ekk_instance_.pivotColumnFtran(variable_in, col_aq);
  // Compute the dual for the pivot column and compare it with the
  // updated value
  double computed_theta_dual =
      ekk_instance_.computeDualForTableauColumn(variable_in, col_aq);
  ekkDebugUpdatedDual(ekk_instance_.options_, updated_theta_dual,
                      computed_theta_dual);

  // Feed in the computed dual value
  simplex_info.workDual_[variable_in] = computed_theta_dual;
  // Reassign theta_dual to be the computed value
  theta_dual = simplex_info.workDual_[variable_in];
  // Determine whether theta_dual is too small or has changed sign
  const bool theta_dual_small = fabs(theta_dual) <= dual_feasibility_tolerance;
  const bool theta_dual_sign_error =
      updated_theta_dual * computed_theta_dual <= 0;

  if (theta_dual_small || theta_dual_sign_error) {
    // The computed dual is small or has a sign error, so don't use it
    std::string theta_dual_size = "";
    if (theta_dual_small) theta_dual_size = "; too small";
    std::string theta_dual_sign = "";
    if (theta_dual_sign_error) theta_dual_sign = "; sign error";
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::kInfo,
                "Chosen entering variable %" HIGHSINT_FORMAT
                " (Iter = %" HIGHSINT_FORMAT "; Update = %" HIGHSINT_FORMAT
                ") has computed "
                "(updated) dual of %10.4g (%10.4g) so don't use it%s%s\n",
                variable_in, ekk_instance_.iteration_count_,
                simplex_info.update_count, computed_theta_dual,
                updated_theta_dual, theta_dual_size.c_str(),
                theta_dual_sign.c_str());
    // If a significant computed dual has sign error, consider reinverting
    if (!theta_dual_small && simplex_info.update_count > 0)
      rebuild_reason = kRebuildReasonPossiblySingularBasis;
    hyperChooseColumnClear();
    return false;
  }
  return true;
}

void HEkkPrimal::phase1ChooseRow() {
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
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
  for (HighsInt i = 0; i < col_aq.count; i++) {
    HighsInt iRow = col_aq.index[i];
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
          baseLower[iRow] > -kHighsInf) {
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
          baseUpper[iRow] < +kHighsInf) {
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
  for (HighsUInt i = 0; i < ph1SorterR.size(); i++) {
    double dMyTheta = ph1SorterR.at(i).first;
    HighsInt index = ph1SorterR.at(i).second;
    HighsInt iRow = index >= 0 ? index : index + num_row;
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
  HighsUInt iLast = ph1SorterT.size();
  for (HighsUInt i = 0; i < ph1SorterT.size(); i++) {
    double dMyTheta = ph1SorterT.at(i).first;
    HighsInt index = ph1SorterT.at(i).second;
    HighsInt iRow = index >= 0 ? index : index + num_row;
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
  for (HighsInt i = iLast - 1; i >= 0; i--) {
    HighsInt index = ph1SorterT.at(i).second;
    HighsInt iRow = index >= 0 ? index : index + num_row;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    if (dAbsAlpha > dMaxAlpha * 0.1) {
      row_out = iRow;
      move_out = index >= 0 ? 1 : -1;
      break;
    }
  }
  analysis->simplexTimerStop(Chuzr2Clock);
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
  for (HighsInt i = 0; i < col_aq.count; i++) {
    HighsInt iRow = col_aq.index[i];
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
  for (HighsInt i = 0; i < col_aq.count; i++) {
    HighsInt iRow = col_aq.index[i];
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

void HEkkPrimal::considerBoundSwap() {
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<double>& workLower = simplex_info.workLower_;
  const vector<double>& workUpper = simplex_info.workUpper_;
  const vector<double>& baseLower = simplex_info.baseLower_;
  const vector<double>& baseUpper = simplex_info.baseUpper_;
  const vector<double>& workValue = simplex_info.workValue_;
  const vector<double>& baseValue = simplex_info.baseValue_;

  // Compute the primal theta and see if we should have done a bound
  // flip instead
  if (row_out < 0) {
    assert(solvePhase == kSolvePhase2);
    // No binding ratio in CHUZR, so flip or unbounded
    theta_primal = move_in * kHighsInf;
    move_out = 0;
  } else {
    // Determine the step to the leaving bound
    //
    alpha_col = col_aq.array[row_out];
    // In Phase 1, move_out depends on whether the leaving variable is
    // becoming feasible - moves up to lower (down to upper) - or
    // remaining feasible - moves down to lower (up to upper) - so
    // can't be set so easily as in phase 2
    if (solvePhase == kSolvePhase2) move_out = alpha_col * move_in > 0 ? -1 : 1;
    theta_primal = 0;
    if (move_out == 1) {
      theta_primal = (baseValue[row_out] - baseUpper[row_out]) / alpha_col;
    } else {
      theta_primal = (baseValue[row_out] - baseLower[row_out]) / alpha_col;
    }
    assert(theta_primal > -kHighsInf && theta_primal < kHighsInf);
  }

  // Look to see if there is a bound flip
  bool flipped = false;
  double lower_in = workLower[variable_in];
  double upper_in = workUpper[variable_in];
  value_in = workValue[variable_in] + theta_primal;
  if (move_in > 0) {
    if (value_in > upper_in + primal_feasibility_tolerance) {
      flipped = true;
      row_out = -1;
      value_in = upper_in;
      theta_primal = upper_in - lower_in;
    }
  } else {
    if (value_in < lower_in - primal_feasibility_tolerance) {
      flipped = true;
      row_out = -1;
      value_in = lower_in;
      theta_primal = lower_in - upper_in;
    }
  }
  const bool pivot_or_flipped = row_out >= 0 || flipped;
  if (solvePhase == kSolvePhase2) {
    // Check for possible unboundedness
    if (!pivot_or_flipped) {
      rebuild_reason = kRebuildReasonPossiblyPrimalUnbounded;
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
  //
  // rebuild_reason = kRebuildReasonPossiblySingularBasis is set if
  // numerical trouble is detected
  updateVerify();
}

void HEkkPrimal::update() {
  // Perform update operations that are independent of phase
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  assert(!rebuild_reason);
  bool flipped = row_out < 0;
  if (flipped) {
    variable_out = variable_in;
    alpha_col = 0;
    numericalTrouble = 0;
    simplex_info.workValue_[variable_in] = value_in;
    assert(ekk_instance_.simplex_basis_.nonbasicMove_[variable_in] = move_in);
    ekk_instance_.simplex_basis_.nonbasicMove_[variable_in] = -move_in;
  } else {
    // Adjust perturbation if leaving equation
    adjustPerturbedEquationOut();
  }

  // Start hyper-sparse CHUZC, that takes place through phase1Update()
  hyperChooseColumnStart();

  if (solvePhase == kSolvePhase1) {
    // Update primal values
    phase1UpdatePrimal();

    // Update the duals with respect to feasibility changes
    basicFeasibilityChangeUpdateDual();

    // For hyper-sparse CHUZC, analyse the duals that have just changed
    hyperChooseColumnBasicFeasibilityChange();

  } else {
    // Update primal values, and identify any infeasibilities
    //
    // rebuild_reason =
    // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set if a
    // primal infeasiblility is found
    phase2UpdatePrimal();
  }

  assert(rebuild_reason == kRebuildReasonNo ||
         rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex);

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

  assert(row_out >= 0);
  // Now set the value of the entering variable
  simplex_info.baseValue_[row_out] = value_in;
  // Consider whether the entering value is feasible and, if not, take
  // action
  //
  // rebuild_reason =
  // kRebuildReasonPrimalInfeasibleInPrimalSimplex is set in
  // phase 2 if a primal infeasiblility is found
  considerInfeasibleValueIn();

  // Update the dual values
  theta_dual = simplex_info.workDual_[variable_in];
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
    rebuild_reason = kRebuildReasonUpdateLimitReached;

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
}

void HEkkPrimal::hyperChooseColumn() {
  if (!use_hyper_chuzc) return;
  if (initialise_hyper_chuzc) return;
  analysis->simplexTimerStart(ChuzcHyperClock);
  const vector<int8_t>& nonbasicMove =
      ekk_instance_.simplex_basis_.nonbasicMove_;
  const vector<int8_t>& nonbasicFlag =
      ekk_instance_.simplex_basis_.nonbasicFlag_;
  const vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  if (report_hyper_chuzc)
    printf(
        "H-S  CHUZC: Max changed measure is %9.4g for column %4" HIGHSINT_FORMAT
        "",
        max_changed_measure_value, max_changed_measure_column);
  double best_measure = max_changed_measure_value;
  variable_in = max_changed_measure_column;
  const bool consider_nonbasic_free_column = nonbasic_free_col_set.count();
  if (num_hyper_chuzc_candidates) {
    for (HighsInt iEntry = 1; iEntry <= num_hyper_chuzc_candidates; iEntry++) {
      HighsInt iCol = hyper_chuzc_candidate[iEntry];
      if (nonbasicFlag[iCol] == kNonbasicFlagFalse) {
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
      printf(
          ", and after HS CHUZC set it is now %9.4g for column "
          "%4" HIGHSINT_FORMAT "",
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
    const double infeasibility, const HighsInt iCol) {
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
  const vector<int8_t>& nonbasicMove =
      ekk_instance_.simplex_basis_.nonbasicMove_;
  HighsInt to_entry;
  const bool use_row_indices = ekk_instance_.sparseLoopStyle(
      row_basic_feasibility_change.count, num_col, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iCol;
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
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iRow;
    if (use_col_indices) {
      iRow = col_basic_feasibility_change.index[iEntry];
    } else {
      iRow = iEntry;
    }
    HighsInt iCol = num_col + iRow;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Any nonbasic free columns will be handled explicitly in
  // hyperChooseColumnDualChange, so only look at them here if not
  // flipping
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (row_out < 0 && num_nonbasic_free_col) {
    const vector<HighsInt>& nonbasic_free_col_set_entry =
        nonbasic_free_col_set.entry();
    for (HighsInt iEntry = 0; iEntry < num_nonbasic_free_col; iEntry++) {
      HighsInt iCol = nonbasic_free_col_set_entry[iEntry];
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
  const vector<int8_t>& nonbasicMove =
      ekk_instance_.simplex_basis_.nonbasicMove_;
  HighsInt to_entry;
  // Look at changes in the columns and assess any dual infeasibility
  const bool use_row_indices =
      ekk_instance_.sparseLoopStyle(row_ap.count, num_col, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iCol;
    if (use_row_indices) {
      iCol = row_ap.index[iEntry];
    } else {
      iCol = iEntry;
    }
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (iCol == check_column && ekk_instance_.iteration_count_ >= check_iter) {
      double measure = dual_infeasibility / devex_weight[iCol];
      if (report_hyper_chuzc) {
        printf("Changing column %" HIGHSINT_FORMAT ": measure = %g \n",
               check_column, measure);
      }
    }
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Look at changes in the rows and assess any dual infeasibility
  const bool use_col_indices =
      ekk_instance_.sparseLoopStyle(row_ep.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iRow;
    if (use_col_indices) {
      iRow = row_ep.index[iEntry];
    } else {
      iRow = iEntry;
    }
    HighsInt iCol = iRow + num_col;
    double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
    if (iCol == check_column && ekk_instance_.iteration_count_ >= check_iter) {
      double measure = dual_infeasibility / devex_weight[iCol];
      if (report_hyper_chuzc) {
        printf("Changing column %" HIGHSINT_FORMAT ": measure = %g \n",
               check_column, measure);
      }
    }
    if (dual_infeasibility > dual_feasibility_tolerance)
      hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  // Look for measure changes in any nonbasic free columns and assess
  // any dual infeasibility
  const HighsInt& num_nonbasic_free_col = nonbasic_free_col_set.count();
  if (num_nonbasic_free_col) {
    const vector<HighsInt>& nonbasic_free_col_set_entry =
        nonbasic_free_col_set.entry();
    for (HighsInt iEntry = 0; iEntry < num_nonbasic_free_col; iEntry++) {
      HighsInt iCol = nonbasic_free_col_set_entry[iEntry];
      double dual_infeasibility = fabs(workDual[iCol]);
      if (dual_infeasibility > dual_feasibility_tolerance)
        hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
    }
  }
  // Assess any dual infeasibility for the leaving column - should be dual
  // feasible!
  HighsInt iCol = variable_out;
  double dual_infeasibility = -nonbasicMove[iCol] * workDual[iCol];
  if (dual_infeasibility > dual_feasibility_tolerance) {
    printf("Dual infeasibility %g for leaving column!\n", dual_infeasibility);
    assert(dual_infeasibility <= dual_feasibility_tolerance);
    hyperChooseColumnChangedInfeasibility(dual_infeasibility, iCol);
  }
  analysis->simplexTimerStop(ChuzcHyperDualClock);
}

void HEkkPrimal::updateDual() {
  analysis->simplexTimerStart(UpdateDualClock);
  assert(alpha_col);
  assert(row_out >= 0);
  vector<double>& workDual = ekk_instance_.simplex_info_.workDual_;
  //  const vector<HighsInt>& nonbasicMove =
  //  ekk_instance_.simplex_basis_.nonbasicMove_;
  // Update the duals
  theta_dual = workDual[variable_in] / alpha_col;
  for (HighsInt iEl = 0; iEl < row_ap.count; iEl++) {
    HighsInt iCol = row_ap.index[iEl];
    workDual[iCol] -= theta_dual * row_ap.array[iCol];
  }
  for (HighsInt iEl = 0; iEl < row_ep.count; iEl++) {
    HighsInt iRow = row_ep.index[iEl];
    HighsInt iCol = iRow + num_col;
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
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const vector<int8_t>& nonbasicFlag =
      ekk_instance_.simplex_basis_.nonbasicFlag_;

  HVector buffer;
  buffer.setup(num_row);
  buffer.clear();
  buffer.count = 0;
  // Accumulate costs for checking
  simplex_info.workCost_.assign(num_tot, 0);
  // Zero the dual values
  simplex_info.workDual_.assign(num_tot, 0);
  // Determine the base value for cost perturbation
  const double base =
      simplex_info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    const double value = simplex_info.baseValue_[iRow];
    const double lower = simplex_info.baseLower_[iRow];
    const double upper = simplex_info.baseUpper_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (!bound_violated) continue;
    double cost = bound_violated;
    if (base) cost *= 1 + base * simplex_info.numTotRandomValue_[iRow];
    buffer.array[iRow] = cost;
    buffer.index[buffer.count++] = iRow;
  }
  if (buffer.count <= 0) {
    // Strange, should be a non-trivial RHS
    assert(buffer.count > 0);
    return;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    simplex_info.workCost_[ekk_instance_.simplex_basis_.basicIndex_[iRow]] =
        buffer.array[iRow];
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

  for (HighsInt iCol = 0; iCol < num_col; iCol++)
    simplex_info.workDual_[iCol] = -nonbasicFlag[iCol] * bufferLong.array[iCol];
  for (HighsInt iRow = 0, iCol = num_col; iRow < num_row; iRow++, iCol++)
    simplex_info.workDual_[iCol] = -nonbasicFlag[iCol] * buffer.array[iRow];
}

void HEkkPrimal::phase1UpdatePrimal() {
  analysis->simplexTimerStart(UpdatePrimalClock);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  col_basic_feasibility_change.clear();
  //
  // Update basic primal values, identifying all the feasibility
  // changes giving a value to col_basic_feasibility_change so that the duals
  // can be updated.
  //
  // Determine the base value for cost perturbation
  const double base =
      simplex_info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  //  if (ekk_instance_.sparseLoopStyle(col_aq.count, num_row, to_entry)) {
  for (HighsInt iEl = 0; iEl < col_aq.count; iEl++) {
    HighsInt iRow = col_aq.index[iEl];
    simplex_info.baseValue_[iRow] -= theta_primal * col_aq.array[iRow];
    HighsInt iCol = ekk_instance_.simplex_basis_.basicIndex_[iRow];
    double was_cost = simplex_info.workCost_[iCol];
    const double value = simplex_info.baseValue_[iRow];
    const double lower = simplex_info.baseLower_[iRow];
    const double upper = simplex_info.baseUpper_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1.0;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1.0;
    }
    double cost = bound_violated;
    if (base) cost *= 1 + base * simplex_info.numTotRandomValue_[iRow];
    simplex_info.workCost_[iCol] = cost;
    if (was_cost) {
      if (!cost) simplex_info.num_primal_infeasibility--;
    } else {
      if (cost) simplex_info.num_primal_infeasibility++;
    }
    double delta_cost = cost - was_cost;
    if (delta_cost) {
      col_basic_feasibility_change.array[iRow] = delta_cost;
      col_basic_feasibility_change.index[col_basic_feasibility_change.count++] =
          iRow;
      if (iCol >= num_col) simplex_info.workDual_[iCol] += delta_cost;
    }
  }
  // Don't set baseValue[row_out] yet so that dual update due to
  // feasibility changes is done correctly
  ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
  analysis->simplexTimerStop(UpdatePrimalClock);
}

void HEkkPrimal::considerInfeasibleValueIn() {
  assert(row_out >= 0);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  // Determine the base value for cost perturbation
  const double base =
      simplex_info.primal_simplex_phase1_cost_perturbation_multiplier * 5e-7;
  const double lower = simplex_info.workLower_[variable_in];
  const double upper = simplex_info.workUpper_[variable_in];
  HighsInt bound_violated = 0;
  if (value_in < lower - primal_feasibility_tolerance) {
    bound_violated = -1;
  } else if (value_in > upper + primal_feasibility_tolerance) {
    bound_violated = 1;
  }
  if (!bound_violated) return;
  // The primal value of the entering variable is not feasible
  if (solvePhase == kSolvePhase1) {
    simplex_info.num_primal_infeasibility++;
    double cost = bound_violated;
    if (base) cost *= 1 + base * simplex_info.numTotRandomValue_[row_out];
    simplex_info.workCost_[variable_in] = cost;
    simplex_info.workDual_[variable_in] += cost;
  } else if (primal_correction_strategy ==
             kSimplexPrimalCorrectionStrategyNone) {
    // @primal_infeasibility calculation
    double primal_infeasibility;
    if (bound_violated < 0) {
      primal_infeasibility = lower - value_in;
    } else {
      primal_infeasibility = value_in - upper;
    }
    simplex_info.num_primal_infeasibility++;
    printf(
        "Entering variable has primal infeasibility of %g for [%g, %g, %g]\n",
        primal_infeasibility, lower, value_in, upper);
    rebuild_reason = kRebuildReasonPrimalInfeasibleInPrimalSimplex;
  } else {
    double bound_shift;
    if (bound_violated > 0) {
      // Perturb the upper bound to accommodate the infeasiblilty
      shiftBound(false, variable_in, value_in,
                 simplex_info.numTotRandomValue_[variable_in],
                 simplex_info.workUpper_[variable_in], bound_shift, true);
      simplex_info.workUpperShift_[variable_in] += bound_shift;
    } else {
      // Perturb the lower bound to accommodate the infeasiblilty
      shiftBound(true, variable_in, value_in,
                 simplex_info.numTotRandomValue_[variable_in],
                 simplex_info.workLower_[variable_in], bound_shift, true);
      simplex_info.workLowerShift_[variable_in] += bound_shift;
    }
    simplex_info.bounds_perturbed = true;
  }
  ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
}

void HEkkPrimal::phase2UpdatePrimal(const bool initialise) {
  static double max_max_local_primal_infeasibility;
  static double max_max_ignored_violation;
  if (initialise) {
    max_max_local_primal_infeasibility = 0;
    max_max_ignored_violation = 0;
    return;
  }
  analysis->simplexTimerStart(UpdatePrimalClock);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  bool primal_infeasible = false;
  double max_local_primal_infeasibility = 0;
  double max_ignored_violation = 0;
  // If shifts are only identified in rebuild() the bounds can be
  // ignored. If they aren't ignored, then violations lead to either
  // identification of infeasiblilities (and return to Phase 1) or
  // shifting of bounds to accommodate them.
  const bool ignore_bounds =
      primal_correction_strategy == kSimplexPrimalCorrectionStrategyInRebuild;
  HighsInt to_entry;
  const bool use_col_indices =
      ekk_instance_.sparseLoopStyle(col_aq.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iRow;
    if (use_col_indices) {
      iRow = col_aq.index[iEntry];
    } else {
      iRow = iEntry;
    }
    simplex_info.baseValue_[iRow] -= theta_primal * col_aq.array[iRow];
    //    if (ignore_bounds) continue;
    // Determine whether a bound is violated and take action
    double lower = simplex_info.baseLower_[iRow];
    double upper = simplex_info.baseUpper_[iRow];
    double value = simplex_info.baseValue_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (!bound_violated) continue;
    // A bound is violated
    if (primal_correction_strategy == kSimplexPrimalCorrectionStrategyNone) {
      // @primal_infeasibility calculation
      double primal_infeasibility;
      if (bound_violated < 0) {
        primal_infeasibility = lower - value;
      } else {
        primal_infeasibility = value - upper;
      }
      max_local_primal_infeasibility =
          max(primal_infeasibility, max_local_primal_infeasibility);
      if (primal_infeasibility > primal_feasibility_tolerance) {
        simplex_info.num_primal_infeasibility++;
        primal_infeasible = true;
      }
    } else if (ignore_bounds) {
      double ignored_violation;
      if (bound_violated < 0) {
        ignored_violation = lower - value;
      } else {
        ignored_violation = value - upper;
      }
      max_ignored_violation = max(ignored_violation, max_ignored_violation);
    } else {
      HighsInt iCol = ekk_instance_.simplex_basis_.basicIndex_[iRow];
      double bound_shift;
      if (bound_violated > 0) {
        // Perturb the upper bound to accommodate the infeasiblilty
        shiftBound(false, iCol, simplex_info.baseValue_[iRow],
                   simplex_info.numTotRandomValue_[iCol],
                   simplex_info.workUpper_[iCol], bound_shift, true);
        simplex_info.baseUpper_[iRow] = simplex_info.workUpper_[iCol];
        simplex_info.workUpperShift_[iCol] += bound_shift;
      } else {
        // Perturb the lower bound to accommodate the infeasiblilty
        shiftBound(true, iCol, simplex_info.baseValue_[iRow],
                   simplex_info.numTotRandomValue_[iCol],
                   simplex_info.workLower_[iCol], bound_shift, true);
        simplex_info.baseLower_[iRow] = simplex_info.workLower_[iCol];
        simplex_info.workLowerShift_[iCol] += bound_shift;
      }
      assert(bound_shift > 0);
    }
  }
  if (primal_infeasible) {
    rebuild_reason = kRebuildReasonPrimalInfeasibleInPrimalSimplex;
    if (max_local_primal_infeasibility >
        max_max_local_primal_infeasibility * 2) {
      max_max_local_primal_infeasibility = max_local_primal_infeasibility;
      printf("phase2UpdatePrimal: max_local_primal_infeasibility = %g\n",
             max_local_primal_infeasibility);
    }
    ekk_instance_.invalidatePrimalMaxSumInfeasibilityRecord();
  }
  if (max_ignored_violation > max_max_ignored_violation * 2) {
    max_max_ignored_violation = max_ignored_violation;
    printf("phase2UpdatePrimal: max_ignored_violation = %g\n",
           max_ignored_violation);
  }
  simplex_info.updated_primal_objective_value +=
      simplex_info.workDual_[variable_in] * theta_primal;

  analysis->simplexTimerStop(UpdatePrimalClock);
}

bool HEkkPrimal::correctPrimal(const bool initialise) {
  if (primal_correction_strategy == kSimplexPrimalCorrectionStrategyNone)
    return true;
  static double max_max_primal_correction;
  if (initialise) {
    max_max_primal_correction = 0;
    return true;
  }
  assert(solvePhase == kSolvePhase2);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsInt num_primal_correction = 0;
  double max_primal_correction = 0;
  double sum_primal_correction = 0;
  HighsInt num_primal_correction_skipped = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double lower = simplex_info.baseLower_[iRow];
    double upper = simplex_info.baseUpper_[iRow];
    double value = simplex_info.baseValue_[iRow];
    HighsInt bound_violated = 0;
    if (value < lower - primal_feasibility_tolerance) {
      bound_violated = -1;
    } else if (value > upper + primal_feasibility_tolerance) {
      bound_violated = 1;
    }
    if (bound_violated) {
      if (simplex_info.allow_bound_perturbation) {
	HighsInt iCol = ekk_instance_.simplex_basis_.basicIndex_[iRow];
	double bound_shift;
	if (bound_violated > 0) {
          // Perturb the upper bound to accommodate the infeasiblilty
          shiftBound(false, iCol, simplex_info.baseValue_[iRow],
                     simplex_info.numTotRandomValue_[iCol],
                     simplex_info.workUpper_[iCol], bound_shift, true);
          simplex_info.baseUpper_[iRow] = simplex_info.workUpper_[iCol];
          simplex_info.workUpperShift_[iCol] += bound_shift;
        } else {
          // Perturb the lower bound to accommodate the infeasiblilty
          shiftBound(true, iCol, simplex_info.baseValue_[iRow],
                     simplex_info.numTotRandomValue_[iCol],
                     simplex_info.workLower_[iCol], bound_shift, true);
          simplex_info.baseLower_[iRow] = simplex_info.workLower_[iCol];
          simplex_info.workLowerShift_[iCol] += bound_shift;
        }
        assert(bound_shift > 0);
        num_primal_correction++;
        max_primal_correction = max(bound_shift, max_primal_correction);
        sum_primal_correction += bound_shift;
        simplex_info.bounds_perturbed = true;
      } else {
        // Bound perturbation is not permitted
        num_primal_correction_skipped++;
      }
    }
  }
  if (num_primal_correction_skipped) {
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::kError,
                "correctPrimal: Missed %d bound shifts\n",
                num_primal_correction_skipped);
    assert(!num_primal_correction_skipped);
    return false;
  }
  if (max_primal_correction > 2 * max_max_primal_correction) {
    printf(
        "phase2CorrectPrimal: num / max / sum primal corrections = "
        "%" HIGHSINT_FORMAT
        " / %g / "
        "%g\n",
        num_primal_correction, max_primal_correction, sum_primal_correction);
    max_max_primal_correction = max_primal_correction;
  }
  return true;
}

void HEkkPrimal::basicFeasibilityChangeUpdateDual() {
  analysis->simplexTimerStart(UpdateDualBasicFeasibilityChangeClock);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
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
  HighsInt to_entry;
  const bool use_row_indices = ekk_instance_.sparseLoopStyle(
      row_basic_feasibility_change.count, num_col, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iCol;
    if (use_row_indices) {
      iCol = row_basic_feasibility_change.index[iEntry];
    } else {
      iCol = iEntry;
    }
    simplex_info.workDual_[iCol] -= row_basic_feasibility_change.array[iCol];
  }
  const bool use_col_indices = ekk_instance_.sparseLoopStyle(
      col_basic_feasibility_change.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iRow;
    if (use_col_indices) {
      iRow = col_basic_feasibility_change.index[iEntry];
    } else {
      iRow = iEntry;
    }
    HighsInt iCol = num_col + iRow;
    simplex_info.workDual_[iCol] -= col_basic_feasibility_change.array[iRow];
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
  const HighsInt solver_num_row = ekk_instance_.simplex_lp_.numRow_;
  if (analysis->analyse_simplex_data)
    analysis->operationRecordBefore(
        ANALYSIS_OPERATION_TYPE_BTRAN_BASIC_FEASIBILITY_CHANGE,
        col_basic_feasibility_change,
        analysis->col_basic_feasibility_change_density);
  ekk_instance_.factor_.btran(col_basic_feasibility_change,
                              analysis->col_basic_feasibility_change_density,
                              analysis->pointer_serial_factor_clocks);
  if (analysis->analyse_simplex_data)
    analysis->operationRecordAfter(
        ANALYSIS_OPERATION_TYPE_BTRAN_BASIC_FEASIBILITY_CHANGE,
        col_basic_feasibility_change);
  const double local_col_basic_feasibility_change_density =
      (double)col_basic_feasibility_change.count / solver_num_row;
  analysis->updateOperationResultDensity(
      local_col_basic_feasibility_change_density,
      analysis->col_basic_feasibility_change_density);
  analysis->simplexTimerStop(BtranBasicFeasibilityChangeClock);
}

void HEkkPrimal::basicFeasibilityChangePrice() {
  analysis->simplexTimerStart(PriceBasicFeasibilityChangeClock);
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  const double local_density =
      1.0 * col_basic_feasibility_change.count / num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  ekk_instance_.choosePriceTechnique(simplex_info.price_strategy, local_density,
                                     use_col_price, use_row_price_w_switch);
  if (analysis->analyse_simplex_data) {
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
    const int8_t* nonbasicFlag = &ekk_instance_.simplex_basis_.nonbasicFlag_[0];
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      row_basic_feasibility_change.array[iCol] *= nonbasicFlag[iCol];
  }
  // Update the record of average row_basic_feasibility_change density
  const double local_row_basic_feasibility_change_density =
      (double)row_basic_feasibility_change.count / num_col;
  analysis->updateOperationResultDensity(
      local_row_basic_feasibility_change_density,
      analysis->row_basic_feasibility_change_density);
  if (analysis->analyse_simplex_data)
    analysis->operationRecordAfter(
        ANALYSIS_OPERATION_TYPE_PRICE_BASIC_FEASIBILITY_CHANGE,
        row_basic_feasibility_change);
  analysis->simplexTimerStop(PriceBasicFeasibilityChangeClock);
}

void HEkkPrimal::resetDevex() {
  devex_weight.assign(num_tot, 1.0);
  devex_index.assign(num_tot, 0);
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    const HighsInt nonbasicFlag =
        ekk_instance_.simplex_basis_.nonbasicFlag_[iCol];
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
  HighsInt to_entry;
  const bool use_col_indices =
      ekk_instance_.sparseLoopStyle(col_aq.count, num_row, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iRow;
    if (use_col_indices) {
      iRow = col_aq.index[iEntry];
    } else {
      iRow = iEntry;
    }
    HighsInt iCol = ekk_instance_.simplex_basis_.basicIndex_[iRow];
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

  for (HighsInt iEl = 0; iEl < row_ap.count; iEl++) {
    HighsInt iCol = row_ap.index[iEl];
    double alpha = row_ap.array[iCol];
    double devex = dPivotWeight * fabs(alpha);
    devex += devex_index[iCol] * 1.0;
    if (devex_weight[iCol] < devex) {
      devex_weight[iCol] = devex;
    }
  }
  for (HighsInt iEl = 0; iEl < row_ep.count; iEl++) {
    HighsInt iRow = row_ep.index[iEl];
    HighsInt iCol = iRow + num_col;
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
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
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
    printf("Numerical check: Iter %4" HIGHSINT_FORMAT
           ": alpha_col = %12g, (From %3s alpha_row = "
           "%12g), aDiff = %12g: measure = %12g\n",
           ekk_instance_.iteration_count_, alpha_col, alpha_row_source.c_str(),
           alpha_row, abs_alpha_diff, numericalTrouble);
  assert(numericalTrouble < 1e-3);
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  //
  if (numericalTrouble > 1e-7 && simplex_info.update_count > 0)
    rebuild_reason = kRebuildReasonPossiblySingularBasis;
}

void HEkkPrimal::iterationAnalysisData() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  analysis->simplex_strategy = kSimplexStrategyPrimal;
  analysis->edge_weight_mode = DualEdgeWeightMode::kDevex;
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
  analysis->numerical_trouble = numericalTrouble;
  analysis->objective_value = simplex_info.updated_primal_objective_value;
  analysis->num_primal_infeasibility = simplex_info.num_primal_infeasibility;
  analysis->num_dual_infeasibility = simplex_info.num_dual_infeasibility;
  analysis->sum_primal_infeasibility = simplex_info.sum_primal_infeasibility;
  analysis->sum_dual_infeasibility = simplex_info.sum_dual_infeasibility;
  if ((analysis->edge_weight_mode == DualEdgeWeightMode::kDevex) &&
      (num_devex_iterations == 0))
    analysis->num_devex_framework++;
}

void HEkkPrimal::iterationAnalysis() {
  iterationAnalysisData();
  analysis->iterationReport();
  if (analysis->analyse_simplex_data) analysis->iterationRecord();
}

void HEkkPrimal::localReportIterHeader() {
  printf(" Iter ColIn Row_Out ColOut\n");
}

void HEkkPrimal::localReportIter(const bool header) {
  if (!report_hyper_chuzc) return;
  static HighsInt last_header_iteration_count;
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsInt iteration_count = ekk_instance_.iteration_count_;
  if (header) {
    localReportIterHeader();
    last_header_iteration_count = iteration_count;
  } else {
    if (ekk_instance_.iteration_count_ > last_header_iteration_count + 10) {
      localReportIterHeader();
      last_header_iteration_count = iteration_count;
    }
    if (row_out >= 0) {
      printf("%5" HIGHSINT_FORMAT " %5" HIGHSINT_FORMAT "  %5" HIGHSINT_FORMAT
             "  %5" HIGHSINT_FORMAT "",
             iteration_count, variable_in, row_out, variable_out);
    } else {
      printf("%5" HIGHSINT_FORMAT " %5" HIGHSINT_FORMAT " Bound flip   ",
             iteration_count, variable_in);
    }
    if (check_column >= 0 && iteration_count >= check_iter) {
      HighsInt flag = ekk_instance_.simplex_basis_.nonbasicFlag_[check_column];
      HighsInt move = ekk_instance_.simplex_basis_.nonbasicMove_[check_column];
      double lower = simplex_info.workLower_[check_column];
      double upper = simplex_info.workUpper_[check_column];
      double value;
      if (flag == kNonbasicFlagTrue) {
        value = simplex_info.workValue_[check_column];
      } else {
        HighsInt iRow;
        for (iRow = 0; iRow < num_row; iRow++) {
          if (ekk_instance_.simplex_basis_.basicIndex_[iRow] == check_column)
            break;
        }
        assert(iRow < num_row);
        value = simplex_info.baseValue_[iRow];
      }
      printf(": Var %2" HIGHSINT_FORMAT " (%1" HIGHSINT_FORMAT
             ", %2" HIGHSINT_FORMAT ") [%9.4g, %9.4g, %9.4g]",
             check_column, flag, move, lower, value, upper);
      if (flag == kNonbasicFlagTrue) {
        double dual = simplex_info.workDual_[check_column];
        double weight = devex_weight[check_column];
        double infeasibility = -move * dual;
        if (lower == -kHighsInf && upper == kHighsInf)
          infeasibility = fabs(dual);
        if (infeasibility < dual_feasibility_tolerance) infeasibility = 0;
        double measure = infeasibility / weight;
        printf(" Du = %9.4g; Wt = %9.4g; Ms = %9.4g", dual, weight, measure);
      }
    }
    printf("\n");
  }
}

void HEkkPrimal::reportRebuild(const HighsInt reason_for_rebuild) {
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
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    bool nonbasic_free =
        simplex_basis.nonbasicFlag_[iCol] == kNonbasicFlagTrue &&
        simplex_info.workLower_[iCol] <= -kHighsInf &&
        simplex_info.workUpper_[iCol] >= kHighsInf;
    if (nonbasic_free) nonbasic_free_col_set.add(iCol);
  }
  //  nonbasic_free_col_set.print();
}

void HEkkPrimal::removeNonbasicFreeColumn() {
  bool remove_nonbasic_free_column =
      ekk_instance_.simplex_basis_.nonbasicMove_[variable_in] == 0;
  if (remove_nonbasic_free_column) {
    bool removed_nonbasic_free_column =
        nonbasic_free_col_set.remove(variable_in);
    if (!removed_nonbasic_free_column) {
      highsLogUser(ekk_instance_.options_.log_options, HighsLogType::kError,
                   "HEkkPrimal::phase1update failed to remove nonbasic free "
                   "column %" HIGHSINT_FORMAT "\n",
                   variable_in);
      assert(removed_nonbasic_free_column);
    }
  }
}

void HEkkPrimal::adjustPerturbedEquationOut() {
  if (!ekk_instance_.simplex_info_.bounds_perturbed) return;
  const HighsLp& simplex_lp = ekk_instance_.simplex_lp_;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  double lp_lower;
  double lp_upper;
  if (variable_out < num_col) {
    lp_lower = simplex_lp.colLower_[variable_out];
    lp_upper = simplex_lp.colUpper_[variable_out];
  } else {
    lp_lower = -simplex_lp.rowUpper_[variable_out - num_col];
    lp_upper = -simplex_lp.rowLower_[variable_out - num_col];
  }
  if (lp_lower < lp_upper) return;
  // Leaving variable is fixed
  //  double save_theta_primal = theta_primal;
  double true_fixed_value = lp_lower;
  // Modify theta_primal so that variable leaves at true fixed value
  theta_primal =
      (simplex_info.baseValue_[row_out] - true_fixed_value) / alpha_col;
  /*
    printf("For equation %4" HIGHSINT_FORMAT " to be nonbasic at RHS %10.4g
    requires theta_primal to change by %10.4g from %10.4g to %10.4g\n",
    variable_out, true_fixed_value, theta_primal-save_theta_primal,
    save_theta_primal, theta_primal);
  */
  simplex_info.workLower_[variable_out] = true_fixed_value;
  simplex_info.workUpper_[variable_out] = true_fixed_value;
  simplex_info.workRange_[variable_out] = 0;
  value_in = simplex_info.workValue_[variable_in] + theta_primal;
}

void HEkkPrimal::getBasicPrimalInfeasibility() {
  // Gets the num/max/sum of basic primal infeasibliities,
  analysis->simplexTimerStart(ComputePrIfsClock);
  const double primal_feasibility_tolerance =
      ekk_instance_.options_.primal_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsInt& num_primal_infeasibility = simplex_info.num_primal_infeasibility;
  double& max_primal_infeasibility = simplex_info.max_primal_infeasibility;
  double& sum_primal_infeasibility = simplex_info.sum_primal_infeasibility;
  const HighsInt updated_num_primal_infeasibility = num_primal_infeasibility;
  num_primal_infeasibility = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibility = 0;

  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double value = simplex_info.baseValue_[iRow];
    double lower = simplex_info.baseLower_[iRow];
    double upper = simplex_info.baseUpper_[iRow];
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (value < lower - primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > primal_feasibility_tolerance)
        num_primal_infeasibility++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibility += primal_infeasibility;
    }
  }
  if (updated_num_primal_infeasibility >= 0) {
    // The number of primal infeasibliities should be correct
    bool num_primal_infeasibility_ok =
        num_primal_infeasibility == updated_num_primal_infeasibility;
    if (!num_primal_infeasibility_ok) {
      printf("In iteration %" HIGHSINT_FORMAT
             ": num_primal_infeasibility = %" HIGHSINT_FORMAT
             " != %" HIGHSINT_FORMAT
             " = "
             "updated_num_primal_infeasibility\n",
             ekk_instance_.iteration_count_, num_primal_infeasibility,
             updated_num_primal_infeasibility);
      assert(num_primal_infeasibility_ok);
    }
  }
  analysis->simplexTimerStop(ComputePrIfsClock);
}

void HEkkPrimal::shiftBound(const bool lower, const HighsInt iVar,
                            const double value, const double random_value,
                            double& bound, double& shift, const bool report) {
  double feasibility = (1 + random_value) * primal_feasibility_tolerance;
  double old_bound = bound;
  std::string type;
  double infeasibility;
  double new_infeasibility;
  if (lower) {
    // Bound to shift is lower
    type = "lower";
    assert(value < bound - primal_feasibility_tolerance);
    infeasibility = bound - value;
    assert(infeasibility > 0);
    // Determine the amount by which value will be feasible - so that
    // it's not degenerate
    shift = infeasibility + feasibility;
    bound -= shift;
    new_infeasibility = bound - value;
    assert(new_infeasibility < 0);
  } else {
    // Bound to shift is upper
    type = "upper";
    assert(value > bound + primal_feasibility_tolerance);
    infeasibility = value - bound;
    assert(infeasibility > 0);
    // Determine the amount by which value will be feasible - so that
    // it's not degenerate
    shift = infeasibility + feasibility;
    bound += shift;
    new_infeasibility = value - bound;
    assert(new_infeasibility < 0);
  }
  double error = fabs(-new_infeasibility - feasibility);
  if (report)
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::kVerbose,
                "Value(%4" HIGHSINT_FORMAT
                ") = %10.4g exceeds %s = %10.4g by %9.4g, so shift bound by "
                "%9.4g to %10.4g: infeasibility %10.4g with error %g\n",
                iVar, value, type.c_str(), old_bound, infeasibility, shift,
                bound, new_infeasibility, error);
}

void HEkkPrimal::savePrimalRay() {
  ekk_instance_.simplex_lp_status_.has_primal_ray = true;
  ekk_instance_.simplex_info_.primal_ray_col_ = variable_in;
  ekk_instance_.simplex_info_.primal_ray_sign_ = -move_in;
}

HighsDebugStatus HEkkPrimal::debugPrimalSimplex(const std::string message,
                                                const bool initialise) {
  HighsDebugStatus return_status = ekkDebugSimplex(
      message, ekk_instance_, algorithm, solvePhase, initialise);
  if (return_status == HighsDebugStatus::kLogicalError) return return_status;
  if (initialise) return return_status;
  return_status = ekkDebugNonbasicFreeColumnSet(ekk_instance_, num_free_col,
                                                nonbasic_free_col_set);
  if (return_status == HighsDebugStatus::kLogicalError) return return_status;
  return HighsDebugStatus::kOk;
}
