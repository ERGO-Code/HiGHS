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

using std::runtime_error; // Just to hack in primal simplex solver
using std::cout;
using std::endl;
using std::flush;

HighsStatus LpStatusToHighsStatus(SimplexSolutionStatus simplex_solution_status) {
  switch (simplex_solution_status) {
  case SimplexSolutionStatus::OUT_OF_TIME:
      return HighsStatus::Timeout;
  case SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      return HighsStatus::ReachedDualObjectiveUpperBound;
  case SimplexSolutionStatus::FAILED:
      return HighsStatus::SolutionError;
  case SimplexSolutionStatus::SINGULAR:
      return HighsStatus::SolutionError;
  case SimplexSolutionStatus::UNBOUNDED:
      return HighsStatus::Unbounded;
  case SimplexSolutionStatus::INFEASIBLE:
      return HighsStatus::Infeasible;
  case SimplexSolutionStatus::OPTIMAL:
      return HighsStatus::Optimal;
  default:
      return HighsStatus::NotImplemented;
  }
}

HighsStatus solveSimplex(
			 const HighsOptions& opt,
                         HighsModelObject& highs_model_object
			 ) {
  // Just solves the LP in highs_model_object.scaled_lp_
  HighsTimer &timer = highs_model_object.timer_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HighsLp &lp = highs_model_object.lp_;


  SimplexTimer simplex_timer;
  simplex_timer.initialiseDualSimplexClocks(highs_model_object);
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

    // Deduce the LP basis from the simplex basis
    printf("!! Convert rather than copy from highs_model_object.simplex_basis_ to highs_model_object.basis_ !!\n");
    //    highs_model_object.basis_ = highs_model_object.simplex_basis_;


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
	
	//    reportLp(lp);
	//    reportLpSolution(highs_model_object);
	HighsStatus result = LpStatusToHighsStatus(simplex_lp_status.solution_status);

	timer.stop(timer.solve_clock);

      }
    }
  }

  if (simplex_info.dual_phase1_iteration_count +
      simplex_info.dual_phase2_iteration_count +
      simplex_info.primal_phase2_iteration_count !=
      simplex_info.iteration_count) {
    printf("Iteration total error \n");
  }


#ifdef HiGHSDEV
  timer.stop(simplex_info.clock_[SimplexTotalClock]);
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
  if (simplex_info.simplex_strategy == SimplexStrategy::PRIMAL) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportDualSimplexInnerClock(highs_model_object);
    }
  } else if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_PLAIN) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportDualSimplexInnerClock(highs_model_object);
    }
    if (simplex_info.report_simplex_outer_clock) {
      simplex_timer.reportDualSimplexIterateClock(highs_model_object);
      simplex_timer.reportDualSimplexOuterClock(highs_model_object);
    }
  }
  
  //  if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_TASKS) {
  //    int reportList[] = {
  //        HTICK_INVERT,        HTICK_CHUZR1,        HTICK_BTRAN,
  //        HTICK_PRICE,         HTICK_CHUZC1,        HTICK_CHUZC2,
  //        HTICK_CHUZC3,        HTICK_DEVEX_WT,      HTICK_FTRAN,
  //        HTICK_FTRAN_BFRT,    HTICK_FTRAN_DSE,     HTICK_UPDATE_DUAL,
  //        HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT, HTICK_UPDATE_FACTOR,
  //        HTICK_GROUP1};
  //    int reportCount = sizeof(reportList) / sizeof(int);
  //    timer.report(reportCount, reportList, 0.0);
  //  }
  
  if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_MULTI) {
    //    int reportList[] = {
    //        HTICK_INVERT,        HTICK_CHUZR1,        HTICK_BTRAN,
    //        HTICK_PRICE,         HTICK_CHUZC1,        HTICK_CHUZC2,
    //        HTICK_CHUZC3,        HTICK_DEVEX_WT,      HTICK_FTRAN,
    //        HTICK_FTRAN_BFRT,    HTICK_FTRAN_DSE,     HTICK_UPDATE_DUAL,
    //        HTICK_UPDATE_PRIMAL, HTICK_UPDATE_WEIGHT, HTICK_UPDATE_FACTOR,
    //        HTICK_UPDATE_ROW_EP};
    //    int reportCount = sizeof(reportList) / sizeof(int);
    //    timer.report(reportCount, reportList, 0.0);
    printf("PAMI   %-20s    CUTOFF  %6g    PERSISTENSE  %6g\n",
	   highs_model_object.lp_.model_name_.c_str(), simplex_info.pami_cutoff,
	   simplex_info.iteration_count / (1.0 + simplex_info.multi_iteration));
  }
  
  if (simplex_info.report_simplex_phases_clock) {
    simplex_timer.reportSimplexTotalClock(highs_model_object);
    simplex_timer.report_simplex_phases_clock(highs_model_object);
  }

  if (simplex_info.analyse_invert_time) {
    double current_run_highs_time = timer.readRunHighsClock();
    int iClock = simplex_info.clock_[InvertClock];
    simplex_info.total_inverts = timer.clock_num_call[iClock];
    simplex_info.total_invert_time = timer.clock_time[iClock];
    
    printf(
	   "Time: Total inverts =  %4d; Total invert  time = %11.4g of Total time = %11.4g",
	   simplex_info.total_inverts, simplex_info.total_invert_time, current_run_highs_time);
    if (current_run_highs_time > 0.001) {
      printf(" (%6.2f%%)\n", (100 * simplex_info.total_invert_time) / current_run_highs_time);
    } else {
      printf("\n");
    }
  }
  if (simplex_info.analyseRebuildTime) {
    double current_run_highs_time = timer.readRunHighsClock();
    HighsClockRecord totalRebuildClock;
    timer.clockInit(totalRebuildClock);
    timer.clockAdd(totalRebuildClock, simplex_info.clock_[IterateDualRebuildClock]);
    timer.clockAdd(totalRebuildClock, simplex_info.clock_[IteratePrimalRebuildClock]);
    int totalRebuilds = 0;
    double totalRebuildTime = 0;
    printf(
        "Time: Total rebuild time = %11.4g (%4d) of Total time = %11.4g",
        totalRebuildTime, totalRebuilds, current_run_highs_time);
    if (current_run_highs_time > 0.001) {
      printf(" (%6.2f%%)\n", (100 * totalRebuildTime) / current_run_highs_time);
    } else {
      printf("\n");
    }
  }
  printf("Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d\n",
	 simplex_info.dual_phase1_iteration_count,
         simplex_info.dual_phase2_iteration_count,
	 simplex_info.primal_phase2_iteration_count,
	 simplex_info.iteration_count);
#endif

#ifdef HiGHSDEV
  //  printf("report_simplex_lp_status_flags(workHMO.simplex_lp_status_)\n");cout<<flush;
  //  report_simplex_lp_status_flags(highs_model_object.simplex_lp_status_);
  if (simplex_info.analyseLpSolution) { util_analyse_lp_solution(highs_model_object);}
#endif

  return LpStatusToHighsStatus(simplex_lp_status.solution_status);

}

// Single function to solve an lp according to options and fill
// solution in solution.
HighsStatus runSimplexSolver(const HighsOptions& opt,
                             HighsModelObject& highs_model_object) {
  HighsSimplexInterface simplex_interface(highs_model_object);
  HighsTimer &timer = highs_model_object.timer_;

  // Set up aliases
  const HighsLp &lp = highs_model_object.lp_;
  HighsScale &scale = highs_model_object.scale_;
  HighsLp &simplex_lp = highs_model_object.simplex_lp_;
  HighsBasis &basis = highs_model_object.basis_;
  SimplexBasis &simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexInfo &simplex_info = highs_model_object.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = highs_model_object.simplex_lp_status_;
  HMatrix &matrix = highs_model_object.matrix_;
  HFactor &factor = highs_model_object.factor_;

  // Set simplex options from HiGHS options
  options(highs_model_object, opt);

  if (!simplex_lp_status.valid) {
    // No valid simplex LP so copy the LP to the structure to be used by the solver
    simplex_lp = lp;

    // Possibly transpose the LP to be solved. This will change the
    // numbers of rows and columns in the LP to be solved
    if (simplex_info.transpose_simplex_lp) transpose_simplex_lp(highs_model_object);

    // Now that the numbers of rows and columns in the LP to be solved
    // are fixed, initialise the real and integer random vectors
    initialise_simplex_lp_random_vectors(highs_model_object);
    //
    // Allocate memory for the basis
    // assignBasis();
    const int numTot = highs_model_object.lp_.numCol_ + highs_model_object.lp_.numRow_;
    simplex_basis.basicIndex_.resize(highs_model_object.lp_.numRow_);
    simplex_basis.nonbasicFlag_.resize(numTot);
    simplex_basis.nonbasicMove_.resize(numTot);
    //
    // Possibly scale the LP to be used by the solver
    //
    // Initialise unit scaling factors, to simplify things if no scaling
    // is performed
    scaleHighsModelInit(highs_model_object);
    if (simplex_info.scale_simplex_lp)
      scale_simplex_lp(highs_model_object);
    //
    // Possibly permute the columns of the LP to be used by the solver. 
    if (simplex_info.permute_simplex_lp)
      permute_simplex_lp(highs_model_object);
    //
    // Possibly tighten the bounds of LP to be used by the solver. 
    if (simplex_info.tighten_simplex_lp)
      tighten_simplex_lp(highs_model_object);
    //
#ifdef HiGHSDEV
    // Analyse the scaled LP
    if (simplex_info.analyseLp) {
      util_analyseLp(lp, "Unscaled");
      if (simplex_lp_status.is_scaled) {
	util_analyseVectorValues("Column scaling factors", lp.numCol_, scale.col_, false);
	util_analyseVectorValues("Row    scaling factors", lp.numRow_, scale.row_, false);
	util_analyseLp(simplex_lp, "Scaled");
      }
    }
    report_simplex_lp_status(highs_model_object.simplex_lp_status_);
#endif
  }
  if (!simplex_basis.valid_ && basis.valid_) {
    // Simplex basis is not valid, but HiGHS basis is valid so convert it to a simplex basis
      simplex_interface.convertHighsToSimplexBasis();
  }
  if (simplex_basis.valid_) {
    // Valid simplex basis so use it to initialise...
    initialise_from_nonbasic(highs_model_object); // initFromNonbasic();
    matrix.setup(simplex_lp.numCol_, simplex_lp.numRow_,
		 &simplex_lp.Astart_[0],
		 &simplex_lp.Aindex_[0],
		 &simplex_lp.Avalue_[0],
		 &simplex_basis.nonbasicFlag_[0]);
  } else {
    // ... otherwise start from a logical basis
    initialise_with_logical_basis(highs_model_object);
    matrix.setup_lgBs(simplex_lp.numCol_, simplex_lp.numRow_,
		      &simplex_lp.Astart_[0],
		      &simplex_lp.Aindex_[0],
		      &simplex_lp.Avalue_[0]);
  
  }
  factor.setup(simplex_lp.numCol_, simplex_lp.numRow_,
		&simplex_lp.Astart_[0],
		&simplex_lp.Aindex_[0],
		&simplex_lp.Avalue_[0],
		&simplex_basis.basicIndex_[0]);

  HighsStatus result = solveSimplex(opt, highs_model_object);

  if (result != HighsStatus::Optimal) return result;

  simplex_interface.convertSimplexToHighsSolution();
  simplex_interface.convertSimplexToHighsBasis();
  return result;
}
#endif
