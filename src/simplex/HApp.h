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
#define SIMPLX_HAPP_H_

// todo: clear includes.
#include <getopt.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "HConfig.h"
#include "simplex/HPrimal.h"
#include "simplex/HDual.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsUtils.h"
//#include "HRanging.h"
#include "simplex/HSimplex.h"
#include "simplex/HighsSimplexInterface.h"
#include "simplex/SimplexConst.h"
#include "simplex/SimplexTimer.h"

// Single function to solve an lp according to options and covert simplex solution and basis
HighsStatus runSimplexSolver(const HighsOptions& opt,
                             HighsModelObject& highs_model_object) {
  HighsSimplexInterface simplex_interface(highs_model_object);
  HighsTimer &timer = highs_model_object.timer_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;

  // Set simplex options from HiGHS options
  options(highs_model_object, opt);

#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "On entry to runSimplexSolver");
#endif
  // Set the default simplex strategy from the options
  SimplexStrategy use_simplex_strategy = opt.simplex_strategy;
  if (opt.clean_up) {
    computePrimalObjectiveValueFromColumnValue(highs_model_object, &highs_model_object.solution_.col_value[0]);
    rebuildPostsolve(highs_model_object);
    // Analyse the basis and solution
    //    printf("\nOn entry to runSimplexSolver\n"); SimplexSolutionStatus lp_status = simplex_interface.analyseHighsSolutionAndSimplexBasis();
  }
  // Possibly set up the LP to be solved by the simplex method. According to options
  //
  // * Transpose the LP to be solved - deprecated since primal simplex solver is better
  //
  // * Scale the LP to be solved
  //
  // * Permute the LP to be solved - good idea to do all the time, but needs permutations to be applied to column solution
  //
  // * Tighten the bounds of LP to be solved - deprecated since presolve is better
  //
  if (!simplex_lp_status.valid) {
    setupSimplexLp(highs_model_object);
#ifdef HiGHSDEV
    //    reportSimplexLpStatus(simplex_lp_status, "After setupSimplexLp");
#endif
  }
  
  SimplexTimer simplex_timer;
  simplex_timer.initialiseSimplexClocks(highs_model_object);

  // Official start of solver
  // Start the solve closk - because setupForSimplexSolve has simplex computations
  timer.start(timer.solve_clock);
#ifdef HiGHSDEV
  timer.start(simplex_info.clock_[SimplexTotalClock]);
#endif

  // Setup the basis if not valid, taking the Highs basis if it's
  // valid, otherwise a unit basis or crash basis and perform INVERT
  // if necessary

  setupForSimplexSolve(highs_model_object);
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "After setupForSimplexSolve");
#endif

  if (use_simplex_strategy == SimplexStrategy::PRIMAL) {
    // Use primal simplex solver
    HPrimal primal_solver(highs_model_object);
    primal_solver.solve();
  } else {
    // Use dual simplex solver
    HDual dual_solver(highs_model_object);
    dual_solver.options();
    // Solve, depending on the particular strategy
    if (use_simplex_strategy == SimplexStrategy::DUAL_PLAIN) {
      // Serial
      dual_solver.solve();
    } else if (use_simplex_strategy == SimplexStrategy::DUAL_TASKS) {
      // Parallel - SIP
      // writePivots("tasks");
      dual_solver.solve(8);
    } else if (use_simplex_strategy == SimplexStrategy::DUAL_MULTI) {
      // Parallel - PAMI
      // writePivots("multi");
      // if (opt.partitionFile.size() > 0) {model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;}
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
      simplex_info.iteration_count) printf("Iteration total error \n");

  if (highs_model_object.options_.simplex_initial_condition_check) {
    timer.start(simplex_info.clock_[BasisConditionClock]);
    double basis_condition = computeBasisCondition(highs_model_object);
    timer.stop(simplex_info.clock_[BasisConditionClock]);
    printf("Optimal basis condition estimate is %g\n", basis_condition);
  }

  // Official finish of solver
  timer.stop(timer.solve_clock);

#ifdef HiGHSDEV
  timer.stop(simplex_info.clock_[SimplexTotalClock]);
  reportSimplexProfiling(highs_model_object);
  // ToDO move iterationAnalysisReport to simplex
  printf("!! Move iterationAnalysisReport() to HSimplex\n");
  //    if (simplex_info.analyseSimplexIterations) iterationAnalysisReport();

  if (use_simplex_strategy == SimplexStrategy::PRIMAL) {
    printf("Iterations [Ph1 %d; Ph2 %d] Total %d\n",
	   simplex_info.primal_phase1_iteration_count,
	   simplex_info.primal_phase2_iteration_count,
	   simplex_info.iteration_count);
  } else {
    printf("Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d\n",
	   simplex_info.dual_phase1_iteration_count,
	   simplex_info.dual_phase2_iteration_count,
	   simplex_info.primal_phase2_iteration_count,
	   simplex_info.iteration_count);
  }
#endif

#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "After solve");
  if (simplex_info.analyseLpSolution) { analyse_lp_solution(highs_model_object);}
#endif
  comparePrimalDualObjectiveValues(highs_model_object);

  HighsStatus result = simplex_interface.LpStatusToHighsStatus(simplex_lp_status.solution_status);

  // Return if not optimal
  if (result != HighsStatus::Optimal) return result;

  // Optimal solution: copy the solution and basis
  simplex_interface.convertSimplexToHighsSolution();
  simplex_interface.convertSimplexToHighsBasis();
  //  printf("\nOn leaving runSimplexSolver\n"); SimplexSolutionStatus lp_status = simplex_interface.analyseHighsSolutionAndSimplexBasis();
  //  rebuildPostsolve(highs_model_object);// Just for checking correctness of rebuildPostsolve checks
  
  return result;
}
#endif
