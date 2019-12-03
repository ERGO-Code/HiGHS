/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
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
#include "util/HighsUtils.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
//#include "HRanging.h"
#include "simplex/HSimplex.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/SimplexConst.h"
#include "simplex/SimplexTimer.h"

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
// unconstrained LPs should be solved in solveModelSimplex
//
// Also sets the solution parameters for the unscaled LP
HighsStatus runSimplexSolver(HighsModelObject& highs_model_object) {
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  HighsTimer& timer = highs_model_object.timer_;
  // Set the solver_return_status to error so that not setting it
  // later is evident
  HighsStatus solver_return_status = HighsStatus::Error;

  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveModelSimplex
  bool positive_num_row = highs_model_object.lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) return HighsStatus::Error;

  // Set simplex options from HiGHS options.
  // ToDo: Should only be done when not hot-starting since strategy
  // knowledge based on run-time experience should be preserved.
  setSimplexOptions(highs_model_object);

  SimplexTimer simplex_timer;
  simplex_timer.initialiseSimplexClocks(highs_model_object);
  //
  // Transition to the best possible simplex basis and solution
  solver_return_status = transition(highs_model_object);
  if (solver_return_status != HighsStatus::OK) return solver_return_status;
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
  int use_simplex_strategy = simplex_info.simplex_strategy;
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
      // Only dual infeasible, so maybe use primal simplex
      if (use_simplex_strategy == SIMPLEX_STRATEGY_CHOOSE)
        use_simplex_strategy = SIMPLEX_STRATEGY_PRIMAL;
    }
  } else {
    // Not primal feasible, so maybe use dual simplex
    if (use_simplex_strategy == SIMPLEX_STRATEGY_CHOOSE)
      use_simplex_strategy = SIMPLEX_STRATEGY_DUAL;
  }
#ifdef HiGHSDEV
  // reportSimplexLpStatus(simplex_lp_status, "After transition");
#endif
  if (highs_model_object.scaled_model_status_ != HighsModelStatus::OPTIMAL) {
    // Official start of solver Start the solve clock - because
    // setupForSimplexSolve has simplex computations
    timer.start(timer.solve_clock);
#ifdef HiGHSDEV
    timer.start(simplex_info.clock_[SimplexTotalClock]);
#endif
    if (use_simplex_strategy == SIMPLEX_STRATEGY_PRIMAL) {
      // Use primal simplex solver
      HPrimal primal_solver(highs_model_object);
      primal_solver.solve();
    } else {
      // Use dual simplex solver
      HDual dual_solver(highs_model_object);
      dual_solver.options();
      // Solve, depending on the particular strategy
      if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL_PLAIN) {
        // Serial
        dual_solver.solve();
      } else if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL_TASKS) {
        // Parallel - SIP
        // writePivots("tasks");
        dual_solver.solve(8);
      } else if (use_simplex_strategy == SIMPLEX_STRATEGY_DUAL_MULTI) {
        // Parallel - PAMI
        // writePivots("multi");
        // if (opt.partitionFile.size() > 0)
        // {model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;}
        bool FourThreads = false;
        bool EightThreads = false;
        if (FourThreads)
          dual_solver.solve(4);
        else if (EightThreads)
          dual_solver.solve(8);
        else
          dual_solver.solve();
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
      HighsLogMessage(HighsMessageType::INFO,
                      "Final basis condition estimate is %g", basis_condition);
    }

    // Official finish of solver
    timer.stop(timer.solve_clock);

    scaled_solution_params.objective_function_value = simplex_info.primal_objective_value;
#ifdef HiGHSDEV
    timer.stop(simplex_info.clock_[SimplexTotalClock]);
    reportSimplexProfiling(highs_model_object);
    // ToDO move iterationAnalysisReport to simplex
    printf("!! Move iterationAnalysisReport() to HSimplex\n");
    //    if (simplex_info.analyseSimplexIterations) iterationAnalysisReport();

    if (use_simplex_strategy == SIMPLEX_STRATEGY_PRIMAL) {
      HighsLogMessage(HighsMessageType::INFO,
                      "Iterations [Ph1 %d; Ph2 %d] Total %d",
                      simplex_info.primal_phase1_iteration_count,
                      simplex_info.primal_phase2_iteration_count,
                      scaled_solution_params.simplex_iteration_count);
    } else {
      HighsLogMessage(HighsMessageType::INFO,
                      "Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d",
                      simplex_info.dual_phase1_iteration_count,
                      simplex_info.dual_phase2_iteration_count,
                      simplex_info.primal_phase2_iteration_count,
                      scaled_solution_params.simplex_iteration_count);
    }
#endif
  }

  if (simplex_info.analyseLpSolution) {
    // Analyse the simplex basic solution, assuming that the scaled solution params are known
    const bool report = true;
    HighsStatus return_status =
      analyseSimplexBasicSolution(highs_model_object,
				  highs_model_object.scaled_solution_params_,
				  report);
    if (return_status != HighsStatus::OK) return return_status;

#ifdef HiGHSDEV
    // For debugging, it's good to be able to check what comes out of
    // the solver. This is done fully by analyseHighsBasicSolution,
    // which computes the primal and dual residuals so isn't cheap.
    //
    // Copy the solution and basis
    HighsSimplexInterface simplex_interface(highs_model_object);
    simplex_interface.convertSimplexToHighsSolution();
    simplex_interface.convertSimplexToHighsBasis();
    return_status = 
      analyseHighsBasicSolution(highs_model_object, "to check simplex basic soulution");
    if (return_status != HighsStatus::OK) return return_status;
    // Invalidate the basis to make sure it is set again later
    // without HiGHSDEV
    highs_model_object.basis_.valid_ = false;
#endif
  }
  solver_return_status = highsStatusFromHighsModelStatus(highs_model_object.scaled_model_status_);
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "After running the simplex solver");
#endif
  return solver_return_status;
}

HighsStatus tryToSolveUnscaledLp(HighsModelObject& highs_model_object) {
  for (int pass = 0; pass < 2; pass++) {
    double new_primal_feasibility_tolerance;
    double new_dual_feasibility_tolerance;
#ifdef HiGHSDEV
    HighsLogMessage(HighsMessageType::INFO, "tryToSolveUnscaledLp pass %1d:", pass);
#endif
    // Deduce the unscaled solution parameters, and new fasibility tolerances if not primal and/or dual feasible
    HighsStatus return_status =
      getNewPrimalDualInfeasibilityTolerancesFromSimplexBasicSolution(highs_model_object, 
								      highs_model_object.unscaled_solution_params_,
								      new_primal_feasibility_tolerance,
								      new_dual_feasibility_tolerance);
    if (return_status != HighsStatus::OK) return return_status;
    int num_unscaled_primal_infeasibilities = highs_model_object.unscaled_solution_params_.num_primal_infeasibilities;
    int num_unscaled_dual_infeasibilities = highs_model_object.unscaled_solution_params_.num_dual_infeasibilities;
    // Set the model and solution status according to the unscaled solution parameters
    if (num_unscaled_primal_infeasibilities==0 && num_unscaled_dual_infeasibilities==0)
      highs_model_object.unscaled_model_status_ = HighsModelStatus::OPTIMAL;

    if (highs_model_object.unscaled_model_status_ == HighsModelStatus::OPTIMAL) return HighsStatus::OK;

    //Not optimal
    assert(num_unscaled_primal_infeasibilities > 0 || num_unscaled_dual_infeasibilities > 0);

    HighsLogMessage(HighsMessageType::INFO,
		    "Have %d primal and %d dual unscaled infeasibilities",
		    num_unscaled_primal_infeasibilities,
		    num_unscaled_dual_infeasibilities);
    HighsLogMessage(HighsMessageType::INFO,
		    "Possibly re-solve with feasibility tolerances of %g primal and %g dual",
		    new_primal_feasibility_tolerance,
		    new_dual_feasibility_tolerance);
    const bool refinement = false;
    if (refinement) {
      HighsLogMessage(HighsMessageType::INFO, "Re-solving with refined tolerances");
      HighsOptions save_options = highs_model_object.options_;
      HighsOptions& options = highs_model_object.options_;
      options.primal_feasibility_tolerance = new_primal_feasibility_tolerance;
      options.dual_feasibility_tolerance = new_dual_feasibility_tolerance;
      options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
      HighsStatus highs_status = runSimplexSolver(highs_model_object);
      options = save_options;
      if (highs_model_object.scaled_model_status_ != HighsModelStatus::OPTIMAL) return highs_status;
    } else {
      HighsLogMessage(HighsMessageType::INFO, "Not re-solving with refined tolerances");
      return HighsStatus::OK;
    }
  }
  return HighsStatus::OK;
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
HighsStatus solveModelSimplex(HighsModelObject& highs_model_object) {
  // Reset unscaled and scaled model status and solution params - except for iteration counts
  resetModelStatusAndSolutionParams(highs_model_object);

  if (!highs_model_object.lp_.numRow_) {
    // Unconstrained LP so solve directly
    HighsStatus return_status = solveUnconstrainedLp(highs_model_object);
    if (return_status != HighsStatus::OK) return HighsStatus::Error;
    return return_status;
  }

  // (Try to) solve the scaled LP
  HighsStatus highs_status = runSimplexSolver(highs_model_object);
#ifdef HiGHSDEV
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  if (simplex_info.analyse_invert_form) reportAnalyseInvertForm(highs_model_object);
#endif
  if (highs_status != HighsStatus::OK) return highs_status;

  double cost_scale = highs_model_object.scale_.cost_;
#ifdef HiGHSDEV
  if (cost_scale != 1) printf("solveModelSimplex: Can't handle cost scaling\n");
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
      highs_status = tryToSolveUnscaledLp(highs_model_object);
      if (highs_status != HighsStatus::OK) return highs_status;
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
    
  // Deduce the HiGHS basis and solution from the simplex basis and solution
  HighsSimplexInterface simplex_interface(highs_model_object);
  simplex_interface.convertSimplexToHighsSolution();
  simplex_interface.convertSimplexToHighsBasis();

  copySolutionIterationCountAndObjectiveParams(highs_model_object.scaled_solution_params_,
					       highs_model_object.unscaled_solution_params_);
  
  // Assess success according to the scaled model status
  return highsStatusFromHighsModelStatus(highs_model_object.scaled_model_status_);
}
#endif
