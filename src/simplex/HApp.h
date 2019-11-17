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
//#include "HRanging.h"
#include "simplex/HSimplex.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/SimplexConst.h"
#include "simplex/SimplexTimer.h"

// Single function to solve an lp according to options and convert
// simplex solution and basis
HighsStatus runSimplexSolver(HighsModelObject& highs_model_object) {
  HighsSimplexInterface simplex_interface(highs_model_object);
  HighsTimer& timer = highs_model_object.timer_;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;

  // Handle the case of unconstrained LPs here
  if (!highs_model_object.lp_.numRow_) {
    setSimplexOptions(highs_model_object);
    HighsStatus solver_return_status = solveUnconstrainedLp(highs_model_object);
    if (solver_return_status == HighsStatus::Error) return solver_return_status;
    int report_level = -1;
#ifdef HiGHSDEV
    report_level = 1;
#endif
    if (simplex_info.analyseLpSolution) {
      HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
      HighsSolutionParams solution_params;
      solution_params.primal_feasibility_tolerance =
	highs_model_object.options_.primal_feasibility_tolerance;
      solution_params.dual_feasibility_tolerance =
	highs_model_object.options_.dual_feasibility_tolerance;
      solution_params.iteration_count =
	highs_model_object.simplex_info_.iteration_count;
      highs_model_object.model_status_ = 
	analyseHighsSolution(highs_model_object.lp_,
			     highs_model_object.basis_,
			     highs_model_object.solution_,
			     solution_params, report_level, "after solving unconstrained LP");
      simplex_info.num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
      simplex_info.max_primal_infeasibility = solution_params.max_primal_infeasibility;
      simplex_info.sum_primal_infeasibilities = solution_params.sum_primal_infeasibilities;
      simplex_info.num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
      simplex_info.max_dual_infeasibility = solution_params.max_dual_infeasibility;
      simplex_info.sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;
    }
    return solver_return_status;
  }

  // Set simplex options from HiGHS options.
  // ToDo: Should only be done when not hot-starting since strategy
  // knowledge based on run-time experience should be preserved.
  setSimplexOptions(highs_model_object);

  SimplexTimer simplex_timer;
  simplex_timer.initialiseSimplexClocks(highs_model_object);

  //
  // Transition to the best possible simplex basis and solution
  highs_model_object.model_status_ = transition(highs_model_object);
  if (highs_model_object.model_status_ == HighsModelStatus::SOLVE_ERROR) {
    return highsStatusFromHighsModelStatus(highs_model_object.model_status_);
  }
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
  if (simplex_info.num_primal_infeasibilities == 0) {
    // Primal feasible
    if (simplex_info.num_dual_infeasibilities == 0) {
      // Dual feasible
      // Simplex solution is optimal
      highs_model_object.model_status_ = HighsModelStatus::OPTIMAL;
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
  if (highs_model_object.model_status_ != HighsModelStatus::OPTIMAL) {
    // Official start of solver Start the solve clock - because
    // setupForSimplexSolve has simplex computations
    timer.start(timer.solve_clock);
#ifdef HiGHSDEV
    timer.start(simplex_info.clock_[SimplexTotalClock]);
#endif
  }
#ifdef HiGHSDEV
  // reportSimplexLpStatus(simplex_lp_status, "After transition");
#endif
  if (highs_model_object.model_status_ != HighsModelStatus::OPTIMAL) {
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
        simplex_info.iteration_count)
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
                      simplex_info.iteration_count);
    } else {
      HighsLogMessage(HighsMessageType::INFO,
                      "Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d",
                      simplex_info.dual_phase1_iteration_count,
                      simplex_info.dual_phase2_iteration_count,
                      simplex_info.primal_phase2_iteration_count,
                      simplex_info.iteration_count);
    }
#endif
  }

  // JAJH 230819: Frig to retain highs_model_object.model_status_ 
  HighsModelStatus save_model_status = highs_model_object.model_status_;
  if (highs_model_object.model_status_ == HighsModelStatus::OPTIMAL) {
    // Optimal solution: copy the solution and basis
    simplex_interface.convertSimplexToHighsSolution();
    simplex_interface.convertSimplexToHighsBasis();
    int report_level=-1;
#ifdef HiGHSDEV
    report_level = 1;
#endif
    if (simplex_info.analyseLpSolution) {
      HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
      HighsSolutionParams solution_params;
      solution_params.primal_feasibility_tolerance =
	highs_model_object.options_.primal_feasibility_tolerance;
      solution_params.dual_feasibility_tolerance =
	highs_model_object.options_.dual_feasibility_tolerance;
      solution_params.iteration_count =	simplex_info.iteration_count;
      highs_model_object.model_status_ = 
	analyseHighsSolution(highs_model_object.lp_,
			     highs_model_object.basis_,
			     highs_model_object.solution_,
			     solution_params, report_level, "after running the simplex solver");
      simplex_info.num_primal_infeasibilities = solution_params.num_primal_infeasibilities;
      simplex_info.max_primal_infeasibility = solution_params.max_primal_infeasibility;
      simplex_info.sum_primal_infeasibilities = solution_params.sum_primal_infeasibilities;
      simplex_info.num_dual_infeasibilities = solution_params.num_dual_infeasibilities;
      simplex_info.max_dual_infeasibility = solution_params.max_dual_infeasibility;
      simplex_info.sum_dual_infeasibilities = solution_params.sum_dual_infeasibilities;
    }
  }
  // JAJH 230819: Frig to retain highs_model_object.model_status_ 
  highs_model_object.model_status_ = save_model_status;
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "After running the simplex solver");
#endif
  return highsStatusFromHighsModelStatus(highs_model_object.model_status_);
}

HighsStatus solveModelSimplex(HighsModelObject& highs_model_object) {
  const bool refinement = false;
  HighsStatus highs_status = runSimplexSolver(highs_model_object);

  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;

#ifdef HiGHSDEV
  if (simplex_info.analyse_invert_form) {
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
  HighsModelStatus model_status = highs_model_object.model_status_;
  if (model_status != HighsModelStatus::OPTIMAL) return highsStatusFromHighsModelStatus(model_status);

  HighsOptions& options = highs_model_object.options_;
  HighsLp& lp = highs_model_object.lp_;
  SimplexBasis& basis = highs_model_object.simplex_basis_;
  HighsScale& scale = highs_model_object.scale_;
  if (!scale.is_scaled_) return highs_status;
  double cost_scale = scale.cost_;
  if (cost_scale != 1) printf("solveModelSimplex: Cant't handle cost scaling\n");
  assert(cost_scale == 1);
  for (int pass = 0; pass < 2; pass++) {
    int num_scaled_dual_infeasibilities = 0;
    double max_scaled_dual_infeasibility = 0;
    double sum_scaled_dual_infeasibilities = 0;
    int num_unscaled_dual_infeasibilities = 0;
    double max_unscaled_dual_infeasibility = 0;
    double sum_unscaled_dual_infeasibilities = 0;
    int num_scaled_primal_infeasibilities = 0;
    double max_scaled_primal_infeasibility = 0;
    double sum_scaled_primal_infeasibilities = 0;
    int num_unscaled_primal_infeasibilities = 0;
    double max_unscaled_primal_infeasibility = 0;
    double sum_unscaled_primal_infeasibilities = 0;
    double new_primal_feasibility_tolerance = simplex_info.primal_feasibility_tolerance;
    double new_dual_feasibility_tolerance = simplex_info.dual_feasibility_tolerance;
    for (int iVar = 0; iVar < lp.numCol_ + lp.numRow_; iVar++) {
      // Look at the nonbasic dual infeasibilities
      if (basis.nonbasicFlag_[iVar] == NONBASIC_FLAG_FALSE) continue;
      // No dual infeasiblity for fixed rows and columns
      if (simplex_info.workLower_[iVar] == simplex_info.workUpper_[iVar]) continue;
      bool col = iVar < lp.numCol_;
      double scale_mu;
      int iCol=0;
      int iRow=0;
      if (col) {
	iCol = iVar;
	scale_mu = 1 / (scale.col_[iCol] / scale.cost_);
      } else {
	iRow = iVar - lp.numCol_;
	scale_mu = scale.row_[iRow] * scale.cost_;
      }
      double scaled_dual = simplex_info.workDual_[iVar];
      double unscaled_dual = scaled_dual * scale_mu;
      double scaled_dual_infeasibility = max(-basis.nonbasicMove_[iVar] * scaled_dual, 0.);
      double unscaled_dual_infeasibility = max(-basis.nonbasicMove_[iVar] * unscaled_dual, 0.);
      if (scaled_dual_infeasibility > simplex_info.dual_feasibility_tolerance) {
	num_scaled_dual_infeasibilities++;
      }
      max_scaled_dual_infeasibility = max(scaled_dual_infeasibility, max_scaled_dual_infeasibility);
      sum_scaled_dual_infeasibilities += scaled_dual_infeasibility;
      if (unscaled_dual_infeasibility > options.dual_feasibility_tolerance) {
	num_unscaled_dual_infeasibilities++;
	double multiplier = options.dual_feasibility_tolerance / scale_mu;
#ifdef HiGHSDEV
	double lower = simplex_info.workLower_[iVar];
	double upper = simplex_info.workUpper_[iVar];
	double value = simplex_info.workValue_[iVar];
	HighsLogMessage(HighsMessageType::INFO,
			"Var %6d (%6d, %6d): [%11.4g, %11.4g, %11.4g] %11.4g s=%11.4g %11.4g: Mu = %g",
			iVar, iCol, iRow, lower, value, upper,
			scaled_dual_infeasibility, scale_mu, unscaled_dual_infeasibility,
			multiplier);
#endif
	new_dual_feasibility_tolerance = min(multiplier, new_dual_feasibility_tolerance);
      }
      max_unscaled_dual_infeasibility = max(unscaled_dual_infeasibility, max_unscaled_dual_infeasibility);
      sum_unscaled_dual_infeasibilities += unscaled_dual_infeasibility;
    }
    for (int ix = 0; ix < lp.numRow_; ix++) {
      // Look at the basic primal infeasibilities
      int iVar = basis.basicIndex_[ix];
      // No dual infeasiblity for fixed rows and columns
      if (simplex_info.workLower_[iVar] == simplex_info.workUpper_[iVar]) continue;
      bool col = iVar < lp.numCol_;
      double scale_mu;
      int iCol=0;
      int iRow=0;
      if (col) {
	iCol = iVar;
	scale_mu = scale.col_[iCol];
      } else {
	iRow = iVar - lp.numCol_;
	scale_mu = 1 / scale.row_[iRow];
      }
      double lower = simplex_info.baseLower_[ix];
      double upper = simplex_info.baseUpper_[ix];
      double value = simplex_info.baseValue_[ix];
      double scaled_primal_infeasibility = max(max(lower-value, value-upper), 0.);
      double unscaled_primal_infeasibility = scaled_primal_infeasibility * scale_mu;
      if (scaled_primal_infeasibility > simplex_info.primal_feasibility_tolerance) {
	num_scaled_primal_infeasibilities++;
      }
      max_scaled_primal_infeasibility = max(scaled_primal_infeasibility, max_scaled_primal_infeasibility);
      sum_scaled_primal_infeasibilities += scaled_primal_infeasibility;
      if (unscaled_primal_infeasibility > options.primal_feasibility_tolerance) {
	num_unscaled_primal_infeasibilities++;
	double multiplier = options.primal_feasibility_tolerance / scale_mu;
#ifdef HiGHSDEV
	HighsLogMessage(HighsMessageType::INFO,
			"Var %6d (%6d, %6d): [%11.4g, %11.4g, %11.4g] %11.4g s=%11.4g %11.4g: Mu = %g",
			iVar, iCol, iRow, lower, value, upper,
			scaled_primal_infeasibility, scale_mu, unscaled_primal_infeasibility,
			multiplier);
#endif
	new_primal_feasibility_tolerance = min(multiplier, new_primal_feasibility_tolerance);
      }
      max_unscaled_primal_infeasibility = max(unscaled_primal_infeasibility, max_unscaled_primal_infeasibility);
      sum_unscaled_primal_infeasibilities += unscaled_primal_infeasibility;
    }
#ifdef HiGHSDEV
    HighsLogMessage(HighsMessageType::INFO, "solveModelSimplex pass %1d:", pass);
    if (num_scaled_primal_infeasibilities>0) {
      HighsLogMessage(HighsMessageType::ERROR, "  Scaled primal infeasibilities: num/max/sum = %6d/%11.4g/%11.4g",
		      num_scaled_primal_infeasibilities,
		      max_scaled_primal_infeasibility,
		      sum_scaled_primal_infeasibilities);
    }
    HighsLogMessage(HighsMessageType::INFO, "Unscaled primal infeasibilities: num/max/sum = %6d/%11.4g/%11.4g",
		    num_unscaled_primal_infeasibilities,
		    max_unscaled_primal_infeasibility,
		    sum_unscaled_primal_infeasibilities);
    if (num_scaled_dual_infeasibilities>0) {
      HighsLogMessage(HighsMessageType::ERROR, "  Scaled   dual infeasibilities: num/max/sum = %6d/%11.4g/%11.4g",
		      num_scaled_dual_infeasibilities,
		      max_scaled_dual_infeasibility,
		      sum_scaled_dual_infeasibilities);
    }
    HighsLogMessage(HighsMessageType::INFO, "Unscaled   dual infeasibilities: num/max/sum = %6d/%11.4g/%11.4g",
		    num_unscaled_dual_infeasibilities,
		    max_unscaled_dual_infeasibility,
		    sum_unscaled_dual_infeasibilities);
#endif
    if (num_unscaled_primal_infeasibilities || num_unscaled_dual_infeasibilities) {
      HighsLogMessage(HighsMessageType::INFO,
		      "Have %d primal and %d dual unscaled infeasibilities so possibly re-solve with infeasibility tolerances of %g primal and %g dual",
		      num_unscaled_primal_infeasibilities,
		      num_unscaled_dual_infeasibilities,
		      new_primal_feasibility_tolerance,
		      new_dual_feasibility_tolerance);
      if (refinement) {
	HighsLogMessage(HighsMessageType::INFO, "Re-solving with refined tolerances");
	HighsOptions save_options = highs_model_object.options_;
	HighsOptions& options = highs_model_object.options_;
	options.primal_feasibility_tolerance = new_primal_feasibility_tolerance;
	options.dual_feasibility_tolerance = new_dual_feasibility_tolerance;
	options.simplex_strategy = SIMPLEX_STRATEGY_CHOOSE;
	HighsStatus highs_status = runSimplexSolver(highs_model_object);
	if (highs_model_object.model_status_ != HighsModelStatus::OPTIMAL) return highs_status;
	options = save_options;
      } else {
	HighsLogMessage(HighsMessageType::INFO, "Not re-solving with refined tolerances");
	return highs_status;
      }
    } else {
      return highs_status;
    }
  }
  return highs_status;
}
#endif
