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
#include "simplex/HCrash.h"
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
#include "simplex/SimplexTimer.h" // Just to hack in primal simplex solver

//using std::cout;
//using std::endl;
//using std::flush;

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

  SimplexTimer simplex_timer;
  simplex_timer.initialiseSimplexClocks(highs_model_object);
#ifdef HiGHSDEV
  timer.start(simplex_info.clock_[SimplexTotalClock]);
#endif

  bool ranging = true;
  // Initialize solver and set dual solver options from simplex options
  if (opt.simplex_strategy == SimplexStrategy::PRIMAL) {
    // Use primal simplex solver
    HPrimal primal_solver(highs_model_object);

    timer.start(timer.solve_clock);
    primal_solver.solve();
    timer.stop(timer.solve_clock);

  } else {
    // Use dual simplex solver
    HDual dual_solver(highs_model_object);
    dual_solver.options();
  
    // If after postsolve. todo: advanced basis start here.
    if (opt.clean_up) {
      timer.start(timer.solve_clock);
      dual_solver.solve();
      timer.stop(timer.solve_clock);
    } else {
      // Crash, if HighsModelObject has basis information.
      if (simplex_info.crash_strategy != SimplexCrashStrategy::OFF) {
	HCrash crash;
	timer.start(timer.crash_clock);
	crash.crash(highs_model_object, 0);
	timer.stop(timer.crash_clock);
      }
      timer.start(timer.solve_clock);
      // Solve, depending on the options.
      // Parallel.
      if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_TASKS) {
	dual_solver.solve(8);
      } else if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_MULTI) {
	//    if (opt.partitionFile.size() > 0) {model.strOption[STROPT_PARTITION_FILE] = opt.partitionFile;}
	dual_solver.solve(8);
#ifdef HiGHSDEV
	// Reinstate this once simplex::writePivots is written
	//    if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_MULTI) writePivots("multi");
	//    if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_TASKS) writePivots("tasks");
#endif
      } else {
	// Serial. Based on previous solvePlainJAJH.
	
	bool FourThreads = false;
	bool EightThreads = false;
	
	if (FourThreads)
	  dual_solver.solve(4);
	else if (EightThreads)
	  dual_solver.solve(8);
	else
	  dual_solver.solve();
	
	HighsStatus result = simplex_interface.LpStatusToHighsStatus(simplex_lp_status.solution_status);

	timer.stop(timer.solve_clock);

      }
    }
  }

  if (simplex_info.dual_phase1_iteration_count +
      simplex_info.dual_phase2_iteration_count +
      simplex_info.primal_phase1_iteration_count +
      simplex_info.primal_phase2_iteration_count !=
      simplex_info.iteration_count) {
    printf("Iteration total error \n");
  }


#ifdef HiGHSDEV
  timer.stop(simplex_info.clock_[SimplexTotalClock]);
  reportSimplexProfiling(highs_model_object);
  printf("!! Move an_bs_cond() to HSimplex\n");
  /*
    if (rp_bs_cond) {
    double bs_cond = an_bs_cond();
    printf("Optimal basis condition estimate is %g\n", bs_cond);
    }
  */
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

  // Setup the basis if not valid, taking the Highs basis if it's
  // valid, otherwise a unit basis
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
