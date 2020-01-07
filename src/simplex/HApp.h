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
#include "lp_data/HighsStatus.h"
#include "simplex/HDual.h"
#include "simplex/HPrimal.h"
#include "simplex/HQPrimal.h"
#include "util/HighsUtils.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
//#include "HRanging.h"
#include "simplex/HSimplex.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/SimplexConst.h"
#include "simplex/SimplexTimer.h"

#ifdef HiGHSDEV
#ifdef OPENMP
#include "omp.h"
#endif
#endif

#ifdef HiGHSDEV
void reportAnalyseInvertForm(const HighsModelObject& highs_model_object) {
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  
  printf("grep_kernel,%s,%s,%d,%d,%d,",
	 highs_model_object.lp_.model_name_.c_str(),
	 highs_model_object.lp_.lp_name_.c_str(),
	 simplex_info.num_invert,
	 simplex_info.num_kernel,
	 simplex_info.num_major_kernel);
  if (simplex_info.num_kernel) printf("%g", simplex_info.sum_kernel_dim/simplex_info.num_kernel);
  printf(",%g,%g,",
	 simplex_info.running_average_kernel_dim,
	 simplex_info.max_kernel_dim);
  if (simplex_info.num_invert) printf("Fill-in,%g", simplex_info.sum_invert_fill_factor/simplex_info.num_invert);
  printf(",");
  if (simplex_info.num_kernel) printf("%g", simplex_info.sum_kernel_fill_factor/simplex_info.num_kernel);
  printf(",");
  if (simplex_info.num_major_kernel) printf("%g", simplex_info.sum_major_kernel_fill_factor/simplex_info.num_major_kernel);
  printf(",%g,%g,%g\n",
	 simplex_info.running_average_invert_fill_factor,
	 simplex_info.running_average_kernel_fill_factor,
	 simplex_info.running_average_major_kernel_fill_factor);
}
#endif

// Single function to solve the (scaled) LP according to
// options. Assumes that the LP has a positive number of rows, since
// unconstrained LPs should be solved in solveLpSimplex
//
// Also sets the solution parameters for the unscaled LP
HighsStatus runSimplexSolver(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsTimer& timer = highs_model_object.timer_;

  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = highs_model_object.lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::ERROR,
		    "runSimplexSolver called for LP with non-positive (%d) number of constraints",
		    highs_model_object.lp_.numRow_);
    return HighsStatus::Error;
  }

  // Set simplex options from HiGHS options.
  // ToDo: Should only be done when not hot-starting since strategy
  // knowledge based on run-time experience should be preserved.
  setSimplexOptions(highs_model_object);

  SimplexTimer simplex_timer;
  simplex_timer.initialiseSimplexClocks(highs_model_object);
  //
  // Transition to the best possible simplex basis and solution
  call_status = transition(highs_model_object);
  return_status = interpretCallStatus(call_status, return_status, "transition");
  if (return_status == HighsStatus::Error) return return_status;
  //
  // Given a simplex basis and solution, use the number of primal and
  // dual infeasibilities to determine whether the simplex solver is
  // needed and, if so, possibly which variant to use.
  //
  // 1. If it is "CHOOSE", in which case an approapriate stratgy is
  // used
  //
  // 2. If re-solving choose the strategy appropriate to primal or
  // dual feasibility
  //
  // Copy the simplex stratgy so that it can be modified:
  //
  const int simplex_strategy = highs_model_object.options_.simplex_strategy;
  const std::string parallel = highs_model_object.options_.parallel;
  // Set use_simplex_strategy to an illegal value to check that it is set
  const int illegal_simplex_strategy = SIMPLEX_STRATEGY_MIN-1;
  int use_simplex_strategy = illegal_simplex_strategy;
  HighsSolutionParams& scaled_solution_params = highs_model_object.scaled_solution_params_;
  if (scaled_solution_params.num_primal_infeasibilities == 0) {
    // Primal feasible
    if (scaled_solution_params.num_dual_infeasibilities == 0) {
      // Dual feasible
      // Simplex solution is optimal
      highs_model_object.scaled_model_status_ = HighsModelStatus::OPTIMAL;
      scaled_solution_params.primal_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
      scaled_solution_params.dual_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
    } else {
      // Only dual infeasible
      if (simplex_strategy == SIMPLEX_STRATEGY_CHOOSE) {
	// Choose to use primal simplex
        use_simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
      } else {
        use_simplex_strategy = simplex_strategy;
      }
    }
  } else {
    // Not primal feasible
    if (simplex_strategy == SIMPLEX_STRATEGY_CHOOSE) {
      // Choose to use serial dual simplex
      use_simplex_strategy = SIMPLEX_STRATEGY_DUAL;
    } else {
      use_simplex_strategy = simplex_strategy;
    }
  }
  // Set min/num/max_threads to correspond to serial code. If parallel
  // stratgies are used, the number of HiGHS threads used will be set
  // to the maximum of the minimum required and the maximum number of
  // HiGHS threads allowed in options_
  //
  // All this is independent of the number of OMP threads available,
  // since code with multiple HiGHS threads can be run in serial.
  const int highs_max_threads = highs_model_object.options_.highs_max_threads;
  simplex_info.max_threads = 1;
  if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL && parallel == on_string) {
  // Check for case where simplex strategy is SIMPLEX_STRATEGY_DUAL
  // but parallel option is chosen
    use_simplex_strategy = SIMPLEX_STRATEGY_DUAL_MULTI;
    simplex_info.max_threads = PARALLEL_THREADS_DEFAULT;
  } else if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS) {
    simplex_info.max_threads = max(DUAL_TASKS_MIN_THREADS, highs_max_threads);
  } else if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
    simplex_info.max_threads = max(DUAL_MULTI_MIN_THREADS, highs_max_threads);
  }
  // Set the number of HiGHS threads to be used to be the maximum
  // number to be used
  simplex_info.num_threads = simplex_info.max_threads;
  // Give a warning if the number of threads to be used is more than
  // the maximum number of HiGHS threads allowed
  if (simplex_info.num_threads > highs_max_threads) {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::WARNING,
		    "Using %d HiGHS threads for parallel strategy rather than maximum number (%d) specified in options",
		    simplex_info.num_threads, highs_max_threads);
  }
  const int highs_min_threads = highs_model_object.options_.highs_min_threads;
  // Set the minimum number of HiGHS threads to be used to be the
  // lesser of the number to be used and the minimum number of HiGHS
  // threads allowed in options_
  simplex_info.min_threads = min(simplex_info.num_threads, highs_min_threads);
  // Give a warning if the minimum number of threads to be used is
  // less than the minimum number of HiGHS threads allowed
  if (simplex_info.min_threads < highs_min_threads) {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::WARNING,
		    "Using minimum of %d HiGHS threads for parallel strategy rather than minimum number (%d) specified in options",
		    simplex_info.min_threads, highs_min_threads);
  }
#ifdef HiGHSDEV
  // reportSimplexLpStatus(simplex_lp_status, "After transition");
#endif
  if (highs_model_object.scaled_model_status_ != HighsModelStatus::OPTIMAL) {
    assert(use_simplex_strategy > illegal_simplex_strategy);
    if (use_simplex_strategy <= illegal_simplex_strategy) {
#ifdef HiGHSDEV
      printf("Failed to assign use_simplex_strategy in runSimplexSolver\n");
#endif
      return HighsStatus::Error;
    }
    simplex_info.simplex_strategy = use_simplex_strategy;
    // Official start of solver Start the solve clock - because
    // setupForSimplexSolve has simplex computations
    timer.start(timer.solve_clock);
#ifdef HiGHSDEV
    timer.start(simplex_info.clock_[SimplexTotalClock]);
#endif
    if (use_simplex_strategy == SIMPLEX_STRATEGY_PRIMAL) {
      // Use primal simplex solver
      HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
		      "Using primal simplex solver");
      HQPrimal primal_solver(highs_model_object);
      call_status = primal_solver.solve();
      return_status = interpretCallStatus(call_status, return_status, "HQPrimal::solve");
      if (return_status == HighsStatus::Error) return return_status;
    } else {
      // Use dual simplex solver
      HDual dual_solver(highs_model_object);
      dual_solver.options();
      // Solve, depending on the particular strategy
      if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS) {
        // Parallel - SIP
	HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
			"Using parallel simplex solver - SIP with %d threads", simplex_info.num_threads);
        // writePivots("tasks");
        call_status = dual_solver.solve();
	return_status = interpretCallStatus(call_status, return_status, "HDual::solve");
	if (return_status == HighsStatus::Error) return return_status;
      } else if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
        // Parallel - PAMI
	HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
			"Using parallel simplex solver - PAMI with %d threads", simplex_info.num_threads);
        // writePivots("multi");
        // if (opt.partitionFile.size() > 0)
        // {model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;}
	call_status = dual_solver.solve();
	  return_status = interpretCallStatus(call_status, return_status, "HDual::solve");
	  if (return_status == HighsStatus::Error) return return_status;
      } else {
        // Serial
	HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
			"Using dual simplex solver - serial");
        call_status = dual_solver.solve();
	return_status = interpretCallStatus(call_status, return_status, "HDual::solve");
	if (return_status == HighsStatus::Error) return return_status;
      }
    }
    if (simplex_info.dual_phase1_iteration_count +
            simplex_info.dual_phase2_iteration_count +
            simplex_info.primal_phase1_iteration_count +
            simplex_info.primal_phase2_iteration_count !=
        scaled_solution_params.simplex_iteration_count)
      printf("Iteration total error \n");

    if (highs_model_object.options_.simplex_initial_condition_check) {
      timer.start(simplex_info.clock_[BasisConditionClock]);
      double basis_condition = computeBasisCondition(highs_model_object);
      timer.stop(simplex_info.clock_[BasisConditionClock]);
      HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                      "Final basis condition estimate is %g", basis_condition);
    }

    // Official finish of solver
    timer.stop(timer.solve_clock);

    scaled_solution_params.objective_function_value = simplex_info.primal_objective_value;
#ifdef HiGHSDEV
    timer.stop(simplex_info.clock_[SimplexTotalClock]);
    reportSimplexProfiling(highs_model_object);

    if (use_simplex_strategy == SIMPLEX_STRATEGY_PRIMAL) {
      HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                      "Iterations [Ph1 %d; Ph2 %d] Total %d",
                      simplex_info.primal_phase1_iteration_count,
                      simplex_info.primal_phase2_iteration_count,
                      scaled_solution_params.simplex_iteration_count);
    } else {
      HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                      "Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d",
                      simplex_info.dual_phase1_iteration_count,
                      simplex_info.dual_phase2_iteration_count,
                      simplex_info.primal_phase2_iteration_count,
                      scaled_solution_params.simplex_iteration_count);
    }
    int omp_max_threads = 0;
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
#endif
  if (omp_max_threads <= 1 && simplex_info.report_HFactor_clock)
    highs_model_object.factor_.reportTimer();
#endif
  }

  if (simplex_info.analyse_lp_solution) {
    // Analyse the simplex basic solution, assuming that the scaled solution params are known
    const bool report = true;
    call_status = analyseSimplexBasicSolution(highs_model_object,
					      highs_model_object.scaled_solution_params_,
					      report);
    return_status = interpretCallStatus(call_status, return_status, "analyseSimplexBasicSolution");
    if (return_status == HighsStatus::Error) return return_status;

#ifdef HiGHSDEV
    // For debugging, it's good to be able to check what comes out of
    // the solver. This is done fully by analyseHighsBasicSolution,
    // which computes the primal and dual residuals so isn't cheap.
    //
    // Copy the solution and basis
    HighsSimplexInterface simplex_interface(highs_model_object);
    simplex_interface.convertSimplexToHighsSolution();
    simplex_interface.convertSimplexToHighsBasis();
    call_status = 
      analyseHighsBasicSolution(highs_model_object.options_.logfile,
				highs_model_object, "to check simplex basic solution");
    return_status = interpretCallStatus(call_status, return_status, "analyseHighsBasicSolution");
    if (return_status == HighsStatus::Error) return return_status;
    // Invalidate the basis to make sure it is set again later
    // without HiGHSDEV
    highs_model_object.basis_.valid_ = false;
#endif
  }
  return_status = highsStatusFromHighsModelStatus(highs_model_object.scaled_model_status_);
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "After running the simplex solver");
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
    // Deduce the unscaled solution parameters, and new fasibility tolerances if not primal and/or dual feasible
    call_status =
      getNewPrimalDualInfeasibilityTolerancesFromSimplexBasicSolution(highs_model_object, 
								      highs_model_object.unscaled_solution_params_,
								      new_primal_feasibility_tolerance,
								      new_dual_feasibility_tolerance);
    return_status = interpretCallStatus(call_status, return_status,
					"getNewPrimalDualInfeasibilityTolerancesFromSimplexBasicSolution");
    if (return_status == HighsStatus::Error) return return_status;
    int num_unscaled_primal_infeasibilities = highs_model_object.unscaled_solution_params_.num_primal_infeasibilities;
    int num_unscaled_dual_infeasibilities = highs_model_object.unscaled_solution_params_.num_dual_infeasibilities;
    // Set the model and solution status according to the unscaled solution parameters
    if (num_unscaled_primal_infeasibilities==0 && num_unscaled_dual_infeasibilities==0)
      highs_model_object.unscaled_model_status_ = HighsModelStatus::OPTIMAL;

    if (highs_model_object.unscaled_model_status_ == HighsModelStatus::OPTIMAL) return HighsStatus::OK;

    //Not optimal
    assert(num_unscaled_primal_infeasibilities > 0 || num_unscaled_dual_infeasibilities > 0);

    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
		    "Have %d primal and %d dual unscaled infeasibilities",
		    num_unscaled_primal_infeasibilities,
		    num_unscaled_dual_infeasibilities);
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
		    "Possibly re-solve with feasibility tolerances of %g primal and %g dual",
		    new_primal_feasibility_tolerance,
		    new_dual_feasibility_tolerance);
    const bool refinement = false;
    if (refinement) {
      HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
		      "Re-solving with refined tolerances");
      highs_model_object.scaled_solution_params_.primal_feasibility_tolerance = new_primal_feasibility_tolerance;
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance = new_dual_feasibility_tolerance;

      HighsOptions save_options = highs_model_object.options_;
      HighsOptions& options = highs_model_object.options_;
      options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
      call_status = runSimplexSolver(highs_model_object);
      options = save_options;
      return_status = interpretCallStatus(call_status, return_status, "runSimplexSolver");
      if (return_status == HighsStatus::Error) return return_status;
      // Assess success according to the scaled model status, unless
      // something worse has happened earlier
      call_status = highsStatusFromHighsModelStatus(highs_model_object.scaled_model_status_);
      return_status = interpretCallStatus(call_status, return_status);
      if (return_status == HighsStatus::Error) return return_status;
    } else {
      HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
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
HighsStatus solveLpSimplex(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::OK;
  HighsStatus call_status;
  // Reset unscaled and scaled model status and solution params - except for iteration counts
  resetModelStatusAndSolutionParams(highs_model_object);
  HighsSimplexAnalysis& simplex_analysis = highs_model_object.simplex_analysis_;
  simplex_analysis.setup(highs_model_object.lp_, highs_model_object.options_,
			 highs_model_object.scaled_solution_params_.simplex_iteration_count);
  if (!highs_model_object.lp_.numRow_) {
    // Unconstrained LP so solve directly
    call_status = solveUnconstrainedLp(highs_model_object);
    return_status = interpretCallStatus(call_status, return_status, "solveUnconstrainedLp");
    if (return_status == HighsStatus::Error) return return_status;
  }
  // (Try to) solve the scaled LP
  call_status = runSimplexSolver(highs_model_object);
  return_status = interpretCallStatus(call_status, return_status, "runSimplexSolver");
  if (return_status == HighsStatus::Error) return return_status;
#ifdef HiGHSDEV
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  if (simplex_info.analyse_invert_form) reportAnalyseInvertForm(highs_model_object);
#endif

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
      return_status = interpretCallStatus(call_status, return_status, "runSimplexSolver");
      if (return_status == HighsStatus::Error) return return_status;
    } else {
      // If scaling hasn't been used, then the original LP has been
      // solved to the required tolerances
      highs_model_object.unscaled_model_status_ = highs_model_object.scaled_model_status_;
      highs_model_object.unscaled_solution_params_ = highs_model_object.scaled_solution_params_;
    }
  } else {
      // If the solution isn't optimal, then clear the scaled solution infeasibility parameters
    highs_model_object.unscaled_model_status_ = highs_model_object.scaled_model_status_;
    invalidateSolutionInfeasibilityParams(highs_model_object.scaled_solution_params_);
  }
    
#ifdef HiGHSDEV
  if (simplex_info.analyse_iterations) simplex_analysis.summaryReport();
#endif

  // Deduce the HiGHS basis and solution from the simplex basis and solution
  HighsSimplexInterface simplex_interface(highs_model_object);
  simplex_interface.convertSimplexToHighsSolution();
  simplex_interface.convertSimplexToHighsBasis();

  copySolutionIterationCountAndObjectiveParams(highs_model_object.scaled_solution_params_,
					       highs_model_object.unscaled_solution_params_);
  
  // Assess success according to the scaled model status, unless
  // something worse has happened earlier
  call_status = highsStatusFromHighsModelStatus(highs_model_object.scaled_model_status_);
  return_status = interpretCallStatus(call_status, return_status);
  return return_status;
}
#endif
