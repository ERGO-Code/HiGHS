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
//#include <stdexcept> // Just to hack in primal simplex solver

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

HighsStatus solveSimplex(
			 const HighsOptions& opt,
                         HighsModelObject& highs_model_object
			 ) {
  // Just solves the LP in highs_model_object.scaled_lp_
  HighsSimplexInterface simplex_interface(highs_model_object);
  HighsTimer &timer = highs_model_object.timer_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &lp = highs_model_object.lp_;

#ifdef HiGHSDEV
  bool compute_basis_condition = true;
  if (compute_basis_condition) {
    double basis_condition = computeBasisCondition(highs_model_object);
    HighsPrintMessage(ML_MINIMAL, "Initial basis condition estimate of %11.4g is", basis_condition);
    if (basis_condition > 1e12) {
      HighsPrintMessage(ML_MINIMAL, " excessive\n");
      return simplex_interface.LpStatusToHighsStatus(SimplexSolutionStatus::FAILED);
    } else {
      HighsPrintMessage(ML_MINIMAL, " OK\n");
    }
  }
  timer.start(simplex_info.clock_[SimplexTotalClock]);
#endif

  // Official start of solver
  timer.start(timer.solve_clock);

  if (opt.simplex_strategy == SimplexStrategy::PRIMAL) {
    // Use primal simplex solver
    HPrimal primal_solver(highs_model_object);
    primal_solver.solve();
  } else {
    // Use dual simplex solver
    HDual dual_solver(highs_model_object);
    dual_solver.options();
    // Solve, depending on the particular strategy
    if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_PLAIN) {
      // Serial
      dual_solver.solve();
    } else if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_TASKS) {
      // Parallel - SIP
      // writePivots("tasks");
      dual_solver.solve(8);
    } else if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_MULTI) {
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
      HighsStatus result = simplex_interface.LpStatusToHighsStatus(simplex_lp_status.solution_status);
    }
  }
  if (simplex_info.dual_phase1_iteration_count +
      simplex_info.dual_phase2_iteration_count +
      simplex_info.primal_phase1_iteration_count +
      simplex_info.primal_phase2_iteration_count !=
      simplex_info.iteration_count) printf("Iteration total error \n");

  // Official finish of solver
  timer.stop(timer.solve_clock);

#ifdef HiGHSDEV
  timer.stop(simplex_info.clock_[SimplexTotalClock]);
  reportSimplexProfiling(highs_model_object);
  if (compute_basis_condition) {
    double basis_condition = computeBasisCondition(highs_model_object);
    printf("Optimal basis condition estimate is %g\n", basis_condition);
  }
  // ToDO move iterateRpAn to simplex
  printf("!! Move iterateRpAn() to HSimplex\n");
  //    if (simplex_info.analyseSimplexIterations) iterateRpAn();

  printf("Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d\n",
	 simplex_info.dual_phase1_iteration_count,
         simplex_info.dual_phase2_iteration_count,
	 simplex_info.primal_phase2_iteration_count,
	 simplex_info.iteration_count);
#endif

#ifdef HiGHSDEV
  //  printf("report_simplex_lp_status_flags(workHMO.simplex_lp_status_)\n");cout<<flush;
  //  report_simplex_lp_status_flags(highs_model_object.simplex_lp_status_);
  if (simplex_info.analyseLpSolution) { analyse_lp_solution(highs_model_object);}
#endif

  return simplex_interface.LpStatusToHighsStatus(simplex_lp_status.solution_status);

}

// Single function to solve an lp according to options and covert simplex solution and basis
HighsStatus runSimplexSolver(const HighsOptions& opt,
                             HighsModelObject& highs_model_object) {
  HighsSimplexInterface simplex_interface(highs_model_object);

  // Set simplex options from HiGHS options
  options(highs_model_object, opt);

  // Possibly set up the LP to be solved by the simplex method. According to options
  //
  // * Transpose the LP to be solved - deprecated since primal simplex solver is better
  // * Scale the LP to be solved
  // * Permute the LP to be solved
  // * Tighten the bounds of LP to be solved - deprecated since presolve is better
  //
  if (!highs_model_object.simplex_lp_status_.valid) setupSimplexLp(highs_model_object);

  SimplexTimer simplex_timer;
  simplex_timer.initialiseSimplexClocks(highs_model_object);

  // Setup the basis if not valid, taking the Highs basis if it's
  // valid, otherwise a unit basis or crash basis and perform INVERT
  // if necessary

  setupForSimplexSolve(highs_model_object);

  // Call the simplex solver
  HighsStatus result = solveSimplex(opt, highs_model_object);

  // Return if not optimal
  if (result != HighsStatus::Optimal) return result;

  // Optimal solution: copy the solution and basis
  simplex_interface.convertSimplexToHighsSolution();
  simplex_interface.convertSimplexToHighsBasis();
  return result;
}
#endif
