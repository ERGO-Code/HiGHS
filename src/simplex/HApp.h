/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef SIMPLEX_HAPP_H_
#define SIMPLEX_HAPP_H_

// todo: clear includes.
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
#include "lp_data/HighsStatus.h"
#include "simplex/HDual.h"
#include "simplex/HEkk.h"
#include "simplex/HPrimal.h"
#include "simplex/HSimplex.h"
#include "simplex/HSimplexDebug.h"
#include "simplex/HSimplexReport.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/SimplexConst.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsUtils.h"

#ifdef OPENMP
#include "omp.h"
#endif

// Single function to solve the (scaled) LP according to
// options. Assumes that the LP has a positive number of rows, since
// unconstrained LPs should be solved in solveLpSimplex
//
// Also sets the solution parameters for the unscaled LP
HighsStatus runHmoSimplexSolver(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  FILE* logfile = highs_model_object.options_.logfile;

  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLp
  bool positive_num_row = highs_model_object.lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "runSimplexSolver called for LP with non-positive (%d) "
                    "number of constraints",
                    highs_model_object.lp_.numRow_);
    return HighsStatus::Error;
  }
#ifdef HiGHSDEV
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  analysis.simplexTimerStart(SimplexTotalClock);
#endif

  // Indicate that dual and primal rays are not known
  highs_model_object.simplex_lp_status_.has_dual_ray = false;
  highs_model_object.simplex_lp_status_.has_primal_ray = false;

  // Transition to the best possible simplex basis and solution
  call_status = transition(highs_model_object);
  return_status = interpretCallStatus(call_status, return_status, "transition");
  if (return_status == HighsStatus::Error) return return_status;

  if (highs_model_object.options_.simplex_class_ekk) {
    HEkk& ekk = highs_model_object.ekk_instance_;
    // Initialise the phase iteration count data
    reportSimplexPhaseIterations(logfile, ekk.iteration_count_,
                                 ekk.simplex_info_, true);
    ekk.passLp(highs_model_object.simplex_lp_);
    call_status = ekk.solve();
    return_status =
        interpretCallStatus(call_status, return_status, "HEkk::solve");
    if (return_status == HighsStatus::Error) return return_status;
    reportSimplexPhaseIterations(
        logfile, ekk.iteration_count_,
        ekk.simplex_info_);  // Check num_basic_logicals
    highs_model_object.scaled_model_status_ = ekk.scaled_model_status_;
    highs_model_object.scaled_solution_params_ = ekk.getSolutionParams();
    if (highs_model_object.scaled_model_status_ ==
            HighsModelStatus::REACHED_TIME_LIMIT ||
        highs_model_object.scaled_model_status_ ==
            HighsModelStatus::REACHED_ITERATION_LIMIT ||
        highs_model_object.scaled_model_status_ ==
            HighsModelStatus::PRIMAL_INFEASIBLE ||
        highs_model_object.scaled_model_status_ ==
            HighsModelStatus::PRIMAL_UNBOUNDED) {
      return return_status;
    }
    simplex_info = ekk.simplex_info_;
    highs_model_object.simplex_basis_ = ekk.simplex_basis_;
    highs_model_object.iteration_counts_.simplex = ekk.iteration_count_;
    int num_row = highs_model_object.simplex_lp_.numRow_;
    int num_col = highs_model_object.simplex_lp_.numCol_;
    int* Astart = &highs_model_object.simplex_lp_.Astart_[0];
    int* Aindex = &highs_model_object.simplex_lp_.Aindex_[0];
    double* Avalue = &highs_model_object.simplex_lp_.Avalue_[0];
    int* nonbasicFlag = &highs_model_object.simplex_basis_.nonbasicFlag_[0];
    int* basicIndex = &highs_model_object.simplex_basis_.basicIndex_[0];
    HighsOptions& options = highs_model_object.options_;
    highs_model_object.matrix_.setup(num_col, num_row, Astart, Aindex, Avalue,
                                     nonbasicFlag);
    highs_model_object.factor_.setup(
        num_col, num_row, Astart, Aindex, Avalue, basicIndex,
        options.highs_debug_level, options.logfile, options.output,
        options.message_level, simplex_info.factor_pivot_threshold,
        options.factor_pivot_tolerance);
    highs_model_object.simplex_lp_status_.has_invert = false;
    int rank_deficiency = computeFactor(highs_model_object);
    if (rank_deficiency) {
      HighsLogMessage(
          highs_model_object.options_.logfile, HighsMessageType::ERROR,
          "Rank deficiency of %d when reinverting in HMO after EkkPrimal",
          rank_deficiency);
      assert(!rank_deficiency);
      return HighsStatus::Error;
    }
    // Use this since data in Ekk now computed without max function
    computeSimplexInfeasible(highs_model_object);
    copySimplexInfeasible(highs_model_object);
    if (simplex_info.num_primal_infeasibilities) {
      // Have primal infeasibilities to clean up
      printf(
          "Leaving EkkPrimal with %d primal and %d dual infeasibilities: "
          "Status %s\n",
          simplex_info.num_primal_infeasibilities,
          simplex_info.num_dual_infeasibilities,
          utilHighsModelStatusToString(highs_model_object.scaled_model_status_)
              .c_str());
    }
    // Should have no dual infeasibilities
    assert(!simplex_info.num_dual_infeasibilities);
    // Use dual simplex (phase 2) with Devex pricing
    highs_model_object.simplex_info_.simplex_strategy = SIMPLEX_STRATEGY_DUAL;
    highs_model_object.simplex_info_.dual_simplex_cost_perturbation_multiplier =
        0;
    highs_model_object.simplex_info_.dual_edge_weight_strategy =
        SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX;
    HighsSimplexAnalysis& simplex_analysis =
        highs_model_object.simplex_analysis_;
    simplex_analysis.setup(highs_model_object.lp_, highs_model_object.options_,
                           highs_model_object.iteration_counts_.simplex);
    HDual dual_solver(highs_model_object);
    dual_solver.options();
    HighsLogMessage(logfile, HighsMessageType::INFO,
                    "Using dual simplex solver - serial");
    call_status = dual_solver.solve();
    return_status =
        interpretCallStatus(call_status, return_status, "HDual::solve");
    if (return_status == HighsStatus::Error) return return_status;
    computeSimplexInfeasible(highs_model_object);
    copySimplexInfeasible(highs_model_object);
    return return_status;
  }

#ifdef HiGHSDEV
  // reportSimplexLpStatus(simplex_lp_status, "After transition");
#endif
  HighsSolutionParams& scaled_solution_params =
      highs_model_object.scaled_solution_params_;
  // Determine whether the scaled solution is optimal
  if (scaled_solution_params.num_primal_infeasibilities == 0 &&
      scaled_solution_params.num_dual_infeasibilities == 0) {
    // No scaled primal or dual infeasiblities => Optimal
    highs_model_object.scaled_model_status_ = HighsModelStatus::OPTIMAL;
    scaled_solution_params.primal_status =
        PrimalDualStatus::STATUS_FEASIBLE_POINT;
    scaled_solution_params.dual_status =
        PrimalDualStatus::STATUS_FEASIBLE_POINT;
  } else {
    // Not optimal
    //
    assert(simplex_info.num_primal_infeasibilities ==
           scaled_solution_params.num_primal_infeasibilities);
    assert(simplex_info.num_dual_infeasibilities ==
           scaled_solution_params.num_dual_infeasibilities);
    highs_model_object.ekk_instance_.chooseSimplexStrategyThreads(
        highs_model_object.options_, simplex_info);
    int simplex_strategy = simplex_info.simplex_strategy;
    reportSimplexPhaseIterations(logfile,
                                 highs_model_object.iteration_counts_.simplex,
                                 highs_model_object.simplex_info_, true);
    if (simplex_strategy == SIMPLEX_STRATEGY_PRIMAL) {
      // Use primal simplex solver
      HighsLogMessage(logfile, HighsMessageType::ERROR,
                      "Primal simplex solver unavailable");
      return HighsStatus::Error;
    } else {
      // Use dual simplex solver
      HDual dual_solver(highs_model_object);
      dual_solver.options();
      // Solve, depending on the particular strategy
      if (simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS) {
        // Parallel - SIP
        HighsLogMessage(logfile, HighsMessageType::INFO,
                        "Using parallel simplex solver - SIP with %d threads",
                        simplex_info.num_threads);
        // writePivots("tasks");
        call_status = dual_solver.solve();
        return_status =
            interpretCallStatus(call_status, return_status, "HDual::solve");
        if (return_status == HighsStatus::Error) return return_status;
      } else if (simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
        // Parallel - PAMI
        HighsLogMessage(logfile, HighsMessageType::INFO,
                        "Using parallel simplex solver - PAMI with %d threads",
                        simplex_info.num_threads);
        call_status = dual_solver.solve();
        return_status =
            interpretCallStatus(call_status, return_status, "HDual::solve");
        if (return_status == HighsStatus::Error) return return_status;
      } else {
        // Serial
        HighsLogMessage(logfile, HighsMessageType::INFO,
                        "Using dual simplex solver - serial");
        call_status = dual_solver.solve();
        return_status =
            interpretCallStatus(call_status, return_status, "HDual::solve");
        if (return_status == HighsStatus::Error) return return_status;
      }

      int& num_scaled_primal_infeasibilities =
          simplex_info.num_primal_infeasibilities;
      if (highs_model_object.scaled_model_status_ ==
              HighsModelStatus::OPTIMAL &&
          num_scaled_primal_infeasibilities) {
        // If Phase 2 primal simplex solver creates primal
        // infeasibilities it doesn't check and may claim
        // optimality. Try again with serial dual solver
        HighsLogMessage(
            logfile, HighsMessageType::WARNING,
            "Phase 2 primal simplex clean-up infeasibilities: Pr %d(Max %9.4g, "
            "Sum %9.4g) so re-solving",
            num_scaled_primal_infeasibilities,
            highs_model_object.scaled_solution_params_.max_primal_infeasibility,
            highs_model_object.scaled_solution_params_
                .sum_primal_infeasibilities);
        call_status = dual_solver.solve();
        return_status =
            interpretCallStatus(call_status, return_status, "HDual::solve");
        if (return_status == HighsStatus::Error) return return_status;
        if (highs_model_object.scaled_model_status_ ==
                HighsModelStatus::OPTIMAL &&
            num_scaled_primal_infeasibilities) {
          // Still optimal with primal infeasibilities
          highs_model_object.scaled_model_status_ = HighsModelStatus::NOTSET;
        }
      }
    }

    computeSimplexInfeasible(highs_model_object);
    copySimplexInfeasible(highs_model_object);

    scaled_solution_params.objective_function_value =
        simplex_info.primal_objective_value;

    if (highs_model_object.scaled_model_status_ == HighsModelStatus::OPTIMAL) {
      highs_model_object.scaled_solution_params_.primal_status =
          PrimalDualStatus::STATUS_FEASIBLE_POINT;
      highs_model_object.scaled_solution_params_.dual_status =
          PrimalDualStatus::STATUS_FEASIBLE_POINT;
    }
    debugBasisCondition(highs_model_object, "Final");

    // Official finish of solver
    reportSimplexPhaseIterations(logfile,
                                 highs_model_object.iteration_counts_.simplex,
                                 highs_model_object.simplex_info_);
  }

#ifdef HiGHSDEV
  analysis.simplexTimerStop(SimplexTotalClock);
#endif
  // Reaches here whether optimal or not
  if (debugSimplexBasicSolution("After runSimplexSolver", highs_model_object) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return HighsStatus::Error;

  return_status =
      highsStatusFromHighsModelStatus(highs_model_object.scaled_model_status_);
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "After running the simplex
  //  solver");
#endif
  return return_status;
}

HighsStatus tryToSolveUnscaledLp(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  for (int pass = 0; pass < 2; pass++) {
    double new_primal_feasibility_tolerance;
    double new_dual_feasibility_tolerance;
#ifdef HiGHSDEV
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "tryToSolveUnscaledLp pass %1d:", pass);
#endif
    // Deduce the unscaled solution parameters, and new fasibility tolerances if
    // not primal and/or dual feasible
    call_status = getNewInfeasibilityTolerancesFromSimplexBasicSolution(
        highs_model_object, highs_model_object.unscaled_solution_params_,
        new_primal_feasibility_tolerance, new_dual_feasibility_tolerance);
    return_status = interpretCallStatus(
        call_status, return_status,
        "getNewInfeasibilityTolerancesFromSimplexBasicSolution");
    if (return_status == HighsStatus::Error) return return_status;
    int num_unscaled_primal_infeasibilities =
        highs_model_object.unscaled_solution_params_.num_primal_infeasibilities;
    int num_unscaled_dual_infeasibilities =
        highs_model_object.unscaled_solution_params_.num_dual_infeasibilities;
    // Set the model and solution status according to the unscaled solution
    // parameters
    if (num_unscaled_primal_infeasibilities == 0 &&
        num_unscaled_dual_infeasibilities == 0) {
      highs_model_object.unscaled_model_status_ = HighsModelStatus::OPTIMAL;
      highs_model_object.unscaled_solution_params_.primal_status =
          PrimalDualStatus::STATUS_FEASIBLE_POINT;
      highs_model_object.unscaled_solution_params_.dual_status =
          PrimalDualStatus::STATUS_FEASIBLE_POINT;
      return HighsStatus::OK;
    }

    // Not optimal
    assert(num_unscaled_primal_infeasibilities > 0 ||
           num_unscaled_dual_infeasibilities > 0);

    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Have %d primal and %d dual unscaled infeasibilities",
                    num_unscaled_primal_infeasibilities,
                    num_unscaled_dual_infeasibilities);
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Possibly re-solve with feasibility tolerances of %g "
                    "primal and %g dual",
                    new_primal_feasibility_tolerance,
                    new_dual_feasibility_tolerance);
    const bool refinement = false;
    if (refinement) {
      HighsLogMessage(highs_model_object.options_.logfile,
                      HighsMessageType::INFO,
                      "Re-solving with refined tolerances");
      highs_model_object.scaled_solution_params_.primal_feasibility_tolerance =
          new_primal_feasibility_tolerance;
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance =
          new_dual_feasibility_tolerance;

      HighsOptions save_options = highs_model_object.options_;
      HighsOptions& options = highs_model_object.options_;
      options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
      call_status = runHmoSimplexSolver(highs_model_object);
      options = save_options;
      return_status =
          interpretCallStatus(call_status, return_status, "runSimplexSolver");
      if (return_status == HighsStatus::Error) return return_status;
      // Assess success according to the scaled model status, unless
      // something worse has happened earlier
      call_status = highsStatusFromHighsModelStatus(
          highs_model_object.scaled_model_status_);
      return_status = interpretCallStatus(call_status, return_status);
      if (return_status == HighsStatus::Error) return return_status;
    } else {
      HighsLogMessage(highs_model_object.options_.logfile,
                      HighsMessageType::INFO,
                      "Not re-solving with refined tolerances");
      return return_status;
    }
  }
  return return_status;
}

// Single method to solve an LP with the simplex method. Solves the
// scaled LP then analyses the unscaled solution. If it doesn't satisfy
// the required tolerances, tolerances for the scaled LP are
// identified which, if used, might yield an unscaled solution that
// satisfies the required tolerances.
//
// This method and tryToSolveUnscaledLp may make mutiple calls to
// runSimplexSolver
//
// It sets the HiGHS basis within highs_model_object and, if optimal,
// the HiGHS solution, too
HighsStatus solveLpHmoSimplex(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Reset unscaled and scaled model status and solution params
  resetTwoModelStatusAndSolutionParams(highs_model_object);
  // Set the value of simplex_info_.run_quiet to suppress computation
  // that is just for reporting
  //  setRunQuiet(highs_model_object);
  //  printf("Forcing simplex_info_.run_quiet true for testing\n");
  //  highs_model_object.simplex_info_.run_quiet = true;

  // Assumes that the LP has a positive number of rows
  bool positive_num_row = highs_model_object.lp_.numRow_ > 0;
  if (!positive_num_row) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::ERROR,
                    "solveLpSimplex called for LP with non-positive (%d) "
                    "number of constraints",
                    highs_model_object.lp_.numRow_);
    assert(positive_num_row);
    return HighsStatus::Error;
  }
  HighsSimplexAnalysis& simplex_analysis = highs_model_object.simplex_analysis_;
  simplex_analysis.setup(highs_model_object.lp_, highs_model_object.options_,
                         highs_model_object.iteration_counts_.simplex);
  //  SimplexTimer simplex_timer;
  //  simplex_timer.initialiseSimplexClocks(highs_model_object);
  // (Try to) solve the scaled LP
  call_status = runHmoSimplexSolver(highs_model_object);
  return_status =
      interpretCallStatus(call_status, return_status, "runSimplexSolver");
  if (return_status == HighsStatus::Error) return return_status;

  double cost_scale = highs_model_object.scale_.cost_;
#ifdef HiGHSDEV
  if (cost_scale != 1) printf("solveLpSimplex: Can't handle cost scaling\n");
#endif
  assert(cost_scale == 1);
  if (cost_scale != 1) return HighsStatus::Error;

  if (highs_model_object.scaled_model_status_ == HighsModelStatus::OPTIMAL) {
    // (Scaled) LP solved to optimality
    if (highs_model_object.scale_.is_scaled_) {
      // LP solved was scaled, so see whether the scaled problem has
      // been solved
      //
      // Analyse the unscaled solution and, if it doesn't satisfy the
      // required tolerances, tolerances for the scaled LP are identified
      // which, if used, might yield an unscaled solution that satisfies
      // the required tolerances. Can't handle cost scaling
      //
      call_status = tryToSolveUnscaledLp(highs_model_object);
      return_status =
          interpretCallStatus(call_status, return_status, "runSimplexSolver");
      if (return_status == HighsStatus::Error) return return_status;
    } else {
      // If scaling hasn't been used, then the original LP has been
      // solved to the required tolerances
      highs_model_object.unscaled_model_status_ =
          highs_model_object.scaled_model_status_;
      highs_model_object.unscaled_solution_params_ =
          highs_model_object.scaled_solution_params_;
    }
  } else {
    // If the solution isn't optimal, then clear the scaled solution
    // infeasibility parameters
    highs_model_object.unscaled_model_status_ =
        highs_model_object.scaled_model_status_;
    invalidateSolutionInfeasibilityParams(
        highs_model_object.scaled_solution_params_);
  }

#ifdef HiGHSDEV
  // Report profiling and analysis for the application of the simplex
  // method to this LP problem
  if (simplex_analysis.analyse_simplex_time)
    reportSimplexProfiling(highs_model_object);
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  if (simplex_info.report_HFactor_clock) simplex_analysis.reportFactorTimer();
  if (simplex_info.analyse_iterations) simplex_analysis.summaryReport();
  simplex_analysis.summaryReportFactor();
#endif

  // Deduce the HiGHS basis and solution from the simplex basis and solution
  HighsSimplexInterface simplex_interface(highs_model_object);
  simplex_interface.convertSimplexToHighsSolution();
  simplex_interface.convertSimplexToHighsBasis();

  copySolutionObjectiveParams(highs_model_object.scaled_solution_params_,
                              highs_model_object.unscaled_solution_params_);

  // Assess success according to the scaled model status, unless
  // something worse has happened earlier
  call_status =
      highsStatusFromHighsModelStatus(highs_model_object.scaled_model_status_);
  return_status = interpretCallStatus(call_status, return_status);
  return return_status;
}
#endif

HighsStatus solveLpEkkSimplex(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::OK;
  //  HighsStatus call_status;
  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsOptions& options = highs_model_object.options_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;

  // Reset the model status and solution parameters for the unscaled
  // LP in case of premature return
  resetModelStatusAndSolutionParams(highs_model_object);

  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = highs_model_object.lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "solveLpEkkSimplex called for LP with non-positive (%d) "
                    "number of constraints",
                    highs_model_object.lp_.numRow_);
    return HighsStatus::Error;
  }

  // If the simplex LP isn't initialised, scale and pass the current LP
  if (!simplex_lp_status.initialised) scaleAndPassLpToEkk(highs_model_object);

  // If there is no simplex basis, use the HiGHS basis
  if (!simplex_lp_status.has_basis && highs_model_object.basis_.valid_) {
    return_status = ekk_instance.setBasis(highs_model_object.basis_);
    if (return_status == HighsStatus::Error) return HighsStatus::Error;
  }

  // Solve the LP!
  return_status = ekk_instance.solve();
  if (return_status == HighsStatus::Error) return HighsStatus::Error;

  // Copy solution data into the HMO
  HighsSolutionParams& solution_params =
      highs_model_object.unscaled_solution_params_;
  highs_model_object.scaled_model_status_ = ekk_instance.scaled_model_status_;
  solution_params.objective_function_value =
      ekk_instance.simplex_info_.primal_objective_value;
  highs_model_object.iteration_counts_.simplex += ekk_instance.iteration_count_;
  highs_model_object.solution_ = ekk_instance.getSolution();
  if (highs_model_object.scale_.is_scaled_)
    unscaleSolution(highs_model_object.solution_, highs_model_object.scale_);
  highs_model_object.basis_ = ekk_instance.getHighsBasis();

  // Determine whether the unscaled LP has been solved
  double new_primal_feasibility_tolerance;
  double new_dual_feasibility_tolerance;
  //  HighsSolutionParams solution_;
  getUnscaledInfeasibilitiesAndNewTolerances(
      ekk_instance.options_, ekk_instance.simplex_lp_,
      ekk_instance.scaled_model_status_, ekk_instance.simplex_basis_,
      ekk_instance.simplex_info_, highs_model_object.scale_, solution_params,
      new_primal_feasibility_tolerance, new_dual_feasibility_tolerance);

  // Handle non-optimal status
  if (ekk_instance.scaled_model_status_ != HighsModelStatus::OPTIMAL) {
    highs_model_object.unscaled_model_status_ =
        ekk_instance.scaled_model_status_;
    return_status =
        highsStatusFromHighsModelStatus(ekk_instance.scaled_model_status_);
    return return_status;
  }

  // Now interpret the status of the unscaled solution when the scaled
  // LP is solved to optimailty
  assert(ekk_instance.scaled_model_status_ == HighsModelStatus::OPTIMAL);

  int num_unscaled_primal_infeasibilities =
      solution_params.num_primal_infeasibilities;
  int num_unscaled_dual_infeasibilities =
      solution_params.num_dual_infeasibilities;
  // Set the model and solution status according to the unscaled solution
  // parameters
  if (num_unscaled_primal_infeasibilities == 0 &&
      num_unscaled_dual_infeasibilities == 0) {
    // Optimal
    highs_model_object.unscaled_model_status_ = HighsModelStatus::OPTIMAL;
    solution_params.primal_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
    solution_params.dual_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
  } else {
    // Not optimal - should try refinement
    highs_model_object.unscaled_model_status_ = HighsModelStatus::NOTSET;
    assert(num_unscaled_primal_infeasibilities > 0 ||
           num_unscaled_dual_infeasibilities > 0);
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "Have num/max/sum primal (%d/%g/%g) and dual (%d/%g/%g) "
                    "unscaled infeasibilities",
                    num_unscaled_primal_infeasibilities,
                    solution_params.max_primal_infeasibility,
                    solution_params.sum_primal_infeasibilities,
                    num_unscaled_dual_infeasibilities,
                    solution_params.max_dual_infeasibility,
                    solution_params.sum_dual_infeasibilities);
    if (ekk_instance.scaled_model_status_ == HighsModelStatus::OPTIMAL)
      HighsLogMessage(
          highs_model_object.options_.logfile, HighsMessageType::INFO,
          "Possibly re-solve with feasibility tolerances of %g "
          "primal and %g dual",
          new_primal_feasibility_tolerance, new_dual_feasibility_tolerance);
    highs_model_object.solution_ = ekk_instance.getSolution();
    if (highs_model_object.scale_.is_scaled_)
      unscaleSolution(highs_model_object.solution_, highs_model_object.scale_);
    highs_model_object.basis_ = ekk_instance.getHighsBasis();
  }
  return return_status;
}

HighsStatus solveLpSimplex(HighsModelObject& highs_model_object) {
  HighsStatus return_status;
  const bool use_solveLpEkkSimplex = true;  // true;
  if (highs_model_object.options_.simplex_class_ekk && use_solveLpEkkSimplex) {
    return_status = solveLpEkkSimplex(highs_model_object);
  } else {
    return_status = solveLpHmoSimplex(highs_model_object);
  }
  return return_status;
}
