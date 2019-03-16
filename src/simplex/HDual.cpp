/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HDual.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HDual.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>

#include "lp_data/HConst.h"
#include "simplex/HCrash.h"
#include "simplex/HPrimal.h"
#include "lp_data/HighsLp.h"
#include "io/HighsIO.h"
#include "lp_data/HighsModelObject.h"
#include "util/HighsTimer.h"
#include "simplex/SimplexTimer.h"
#include "simplex/HSimplex.h"

using std::runtime_error;
using std::cout;
using std::endl;
using std::flush;
using std::fabs;

void HDual::solve(int num_threads) {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;
  simplex_lp_status.solution_status = SimplexSolutionStatus::UNSET;
  // Cannot solve box-constrained LPs
  if (workHMO.simplex_lp_.numRow_ == 0) return;

  HighsTimer &timer = workHMO.timer_;
  invertHint = INVERT_HINT_NO;

  SimplexTimer simplex_timer;
  simplex_timer.initialiseDualSimplexClocks(workHMO);

  // Setup aspects of the model data which are needed for solve() but better
  // left until now for efficiency reasons.
#ifdef HiGHSDEV
  printf("Calling setup_for_solve(workHMO);\n");
#endif
  setup_for_solve(workHMO); //  model->setup_for_solve();

#ifdef HiGHSDEV
  timer.start(simplex_info.clock_[SimplexTotalClock]);
#endif
  // Set SolveBailout to be true if control is to be returned immediately to
  // calling function
  SolveBailout = false;

  // Initialise working environment
  // Does LOTS, including initialisation of edge weights. Should only
  // be called if model dimension changes
  init(num_threads);

  initialise_cost(workHMO, 1); //  model->initCost(1);
  if (!simplex_lp_status.has_fresh_invert) {
    int rankDeficiency = compute_factor(workHMO); // int rankDeficiency = model->computeFactor();

    if (rankDeficiency) {
      throw runtime_error("Dual initialise: singular-basis-matrix");
    }
#ifdef HiGHSDEV
    double bsCond = an_bs_cond();
    HighsPrintMessage(ML_MINIMAL, "Initial basis condition estimate of %11.4g is", bsCond);
    if (bsCond > 1e12) {
      HighsPrintMessage(ML_MINIMAL, " excessive\n");
      return;
    } else {
      HighsPrintMessage(ML_MINIMAL, " OK\n");
    }
#endif
  }
  // Consider initialising edge weights
  //
  // NB workEdWt is assigned and initialised to 1s in
  // dualRHS.setup(workHMO) so that CHUZR is well defined, even for
  // Dantzig pricing
  //
#ifdef HiGHSDEV
  //  printf("simplex_lp_status.has_dual_steepest_edge_weights 2 = %d; dual_edge_weight_mode = %d; DualEdgeWeightMode::STEEPEST_EDGE =
  //  %d\n",
  //	 simplex_lp_status.has_dual_steepest_edge_weights, dual_edge_weight_mode, DualEdgeWeightMode::STEEPEST_EDGE);cout<<flush;
  //  printf("Edge weights known? %d\n", !simplex_lp_status.has_dual_steepest_edge_weights);cout<<flush;
#endif
  if (!simplex_lp_status.has_dual_steepest_edge_weights) {
    // Edge weights are not known
    // Set up edge weights according to dual_edge_weight_mode and initialise_dual_steepest_edge_weights
    if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
      // Using dual Devex edge weights
      // Zero the number of Devex frameworks used and set up the first one
      n_dvx_fwk = 0;
      dvx_ix.assign(solver_num_tot, 0);
      iz_dvx_fwk();
    } else if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // Using dual steepest edge (DSE) weights
      int num_basic_structurals = solver_num_row - simplex_info.num_basic_logicals;
      bool computeExactDseWeights = num_basic_structurals > 0 && initialise_dual_steepest_edge_weights;
#ifdef HiGHSDEV
      n_wg_DSE_wt = 0;
      if (computeExactDseWeights) {
	printf("If (0<num_basic_structurals = %d) && %d = initialise_dual_steepest_edge_weights: Compute exact "
	       "DSE weights\n", num_basic_structurals, initialise_dual_steepest_edge_weights);
      }
#endif
      if (computeExactDseWeights) {
        // Basis is not logical and DSE weights are to be initialised
#ifdef HiGHSDEV
        printf("Compute exact DSE weights\n");  // int RpI = 1;
	int iClock = simplex_info.clock_[SimplexIzDseWtClock];
	timer.start(iClock);
#endif
        for (int i = 0; i < solver_num_row; i++) {
#ifdef HiGHSDEV
          //	  if (i==RpI) {printf("Computing exact DSE weight %d\n", i); RpI
          //= RpI*2;}
#endif
          row_ep.clear();
          row_ep.count = 1;
          row_ep.index[0] = i;
          row_ep.array[i] = 1;
          row_ep.packFlag = false;
          factor->btran(row_ep, row_epDensity);
          dualRHS.workEdWt[i] = row_ep.norm2();
          double lc_OpRsDensity = (double)row_ep.count / solver_num_row;
          uOpRsDensityRec(lc_OpRsDensity, row_epDensity);
        }
#ifdef HiGHSDEV
	timer.stop(iClock);
        double IzDseWtTT = timer.read(iClock);
        HighsPrintMessage(ML_DETAILED, "Computed %d initial DSE weights in %gs\n",
			  solver_num_row, IzDseWtTT);
#endif
      }
#ifdef HiGHSDEV
      else {
	HighsPrintMessage(ML_DETAILED, "solve:: %d basic structurals: starting from B=I so unit initial DSE weights\n",
			  num_basic_structurals);
      }
#endif
    }
    // Indicate that edge weights are known
    simplex_lp_status.has_dual_steepest_edge_weights = true;
  }

#ifdef HiGHSDEV
  bool rp_bs_cond = false;
  if (rp_bs_cond) {
    double bs_cond = an_bs_cond();
    printf("Initial basis condition estimate is %g\n", bs_cond);
  }
#endif

  compute_dual(workHMO); //  model->computeDual();
  compute_dual_infeasible_in_dual(workHMO, &dualInfeasCount);//model->computeDualInfeasInDual(&dualInfeasCount);
  solvePhase = dualInfeasCount > 0 ? 1 : 2;

  // Find largest dual. No longer adjust the dual tolerance accordingly
  double largeDual = 0;
  for (int i = 0; i < solver_num_tot; i++) {
    if (workHMO.simplex_basis_.nonbasicFlag_[i]) {
      double myDual = fabs(workDual[i] * jMove[i]);
      if (largeDual < myDual) largeDual = myDual;
    }
  }
#ifdef HiGHSDEV
  //  printf("Solve: Large dual = %g\n", largeDual);cout<<flush;
#endif

  // Check that the model is OK to solve:
  //
  // Level 0 just checks the flags
  //
  // Level 1 also checks that the basis is OK and that the necessary
  // data in work* is populated.
  //
  // Level 2 (will) checks things like the nonbasic duals and basic
  // primal values
  //
  // Level 3 (will) checks expensive things like the INVERT and
  // steepeest edge weights
  //
#ifdef HiGHSDEV
  //  if ((solvePhase != 1) && (solvePhase != 2)) {printf("In solve():
  //  solvePhase = %d\n", solvePhase);cout<<flush;}
#endif
  bool ok = ok_to_solve(workHMO, 1, solvePhase);//  model->OKtoSolve(1, solvePhase);
  if (!ok) {printf("NOT OK TO SOLVE???\n"); cout << flush;}
  assert(ok);

#ifdef HiGHSDEV
  //  Analyse the initial values of primal and dual variables
  //  an_iz_vr_v();
#endif

  // The major solving loop

  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  iterateIzAn();

  while (solvePhase) {
    int it0 = simplex_info.iteration_count;
#ifdef HiGHSDEV
    double simplexTotalTime = timer.read(simplex_info.clock_[SimplexTotalClock]);
    // printf("HDual::solve Phase %d: Iteration %d; simplexTotalTime = %g\n",
    // solvePhase, simplex_info.iteration_count, simplexTotalTime);cout<<flush;
#endif
    // When starting a new phase the (updated) dual objective function
    // value isn't known. Indicate this so that when the value
    // computed from scratch in build() isn't checked against the the
    // updated value
    simplex_lp_status.has_dual_objective_value = 0;
    switch (solvePhase) {
      case 1:
	timer.start(simplex_info.clock_[SimplexDualPhase1Clock]);
        solve_phase1();
	timer.stop(simplex_info.clock_[SimplexDualPhase1Clock]);
        simplex_info.dual_phase1_iteration_count += (simplex_info.iteration_count - it0);
        break;
      case 2:
	timer.start(simplex_info.clock_[SimplexDualPhase2Clock]);
        solve_phase2();
	timer.stop(simplex_info.clock_[SimplexDualPhase2Clock]);
        simplex_info.dual_phase2_iteration_count += (simplex_info.iteration_count - it0);
        break;
      case 4:
        break;
      default:
        solvePhase = 0;
        break;
    }
    // Jump for primal
    if (solvePhase == 4) break;
    // Possibly bail out
    if (SolveBailout) break;
  }

#ifdef HiGHSDEV
  if (simplex_info.analyseSimplexIterations) iterateRpAn();
  // Report the ticks before primal
  if (simplex_info.simplex_strategy == SimplexStrategy::DUAL_PLAIN) {
    if (simplex_info.report_simplex_inner_clock) {
      simplex_timer.reportDualSimplexInnerClock(workHMO);
    }
    if (simplex_info.report_simplex_outer_clock) {
      simplex_timer.reportDualSimplexIterateClock(workHMO);
      simplex_timer.reportDualSimplexOuterClock(workHMO);
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
             workHMO.lp_.model_name_.c_str(), pami_cutoff,
             simplex_info.iteration_count / (1.0 + multi_iteration));
    }
#endif

  if (simplex_lp_status.solution_status != SimplexSolutionStatus::OUT_OF_TIME) {
    // Use primal to clean up if not out of time
    int it0 = simplex_info.iteration_count;
    if (solvePhase == 4) {
      HPrimal hPrimal(workHMO);

      timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
      hPrimal.solvePhase2();
      timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);

    }
    simplex_info.primal_phase2_iteration_count += (simplex_info.iteration_count - it0);
  }
  // Save the solved results
#ifdef HiGHSDEV
  if (simplex_info.dual_phase1_iteration_count +
      simplex_info.dual_phase2_iteration_count +
      simplex_info.primal_phase2_iteration_count !=
      simplex_info.iteration_count) {
    printf("Iteration total error \n");
  }
  printf("Iterations [Ph1 %d; Ph2 %d; Pr %d] Total %d\n",
	 simplex_info.dual_phase1_iteration_count,
         simplex_info.dual_phase2_iteration_count,
	 simplex_info.primal_phase2_iteration_count,
	 simplex_info.iteration_count);
  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    printf("Devex: n_dvx_fwk = %d; Average n_dvx_it = %d\n", n_dvx_fwk,
           simplex_info.iteration_count / n_dvx_fwk);
  }
  if (rp_bs_cond) {
    double bs_cond = an_bs_cond();
    printf("Optimal basis condition estimate is %g\n", bs_cond);
  }
#endif
#ifdef HiGHSDEV
  //  if ((solvePhase != 1) && (solvePhase != 2)) {printf("In solve():
  //  solvePhase = %d\n", solvePhase);cout<<flush;}
#endif
  ok = ok_to_solve(workHMO, 1, solvePhase);// model->OKtoSolve(1, solvePhase);
  if (!ok) {printf("NOT OK After Solve???\n"); cout << flush;}
  assert(ok);
#ifdef HiGHSDEV
  //  printf("report_simplex_lp_status_flags(workHMO.simplex_lp_status_)\n");cout<<flush;
  //  report_simplex_lp_status_flags(workHMO.simplex_lp_status_);
  timer.stop(simplex_info.clock_[SimplexTotalClock]);
  double simplexTotalTime = timer.read(simplex_info.clock_[SimplexTotalClock]);

  if (simplex_info.report_simplex_phases_clock) {
    simplex_timer.reportSimplexTotalClock(workHMO);
    simplex_timer.report_simplex_phases_clock(workHMO);
  }
#endif

#ifdef HiGHSDEV
  if (simplex_info.analyseLpSolution) { util_analyse_lp_solution(workHMO);}
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

#endif
}

void HDual::options() {
  // Set solver options from simplex options

  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;

  interpret_dual_edge_weight_strategy(simplex_info.dual_edge_weight_strategy);
  interpret_price_strategy(simplex_info.price_strategy);

  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance = simplex_info.primal_feasibility_tolerance;
  dual_feasibility_tolerance = simplex_info.dual_feasibility_tolerance;
  dual_objective_value_upper_bound = simplex_info.dual_objective_value_upper_bound;
  //  perturb_costs = simplex_info.perturb_costs;
  //  iterationLimit = simplex_info.iterationLimit;

  // Set values of internal options
}

void HDual::init(int num_threads) {
  // Copy size, matrix and factor

  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  matrix = &workHMO.matrix_;
  factor = &workHMO.factor_;

  // Copy pointers
  jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  workDual = &workHMO.simplex_info_.workDual_[0];
  workValue = &workHMO.simplex_info_.workValue_[0];
  workRange = &workHMO.simplex_info_.workRange_[0];
  baseLower = &workHMO.simplex_info_.baseLower_[0];
  baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  baseValue = &workHMO.simplex_info_.baseValue_[0];

  // Copy tolerances
  Tp = primal_feasibility_tolerance;
  Td = dual_feasibility_tolerance;

  // Setup local vectors
  columnDSE.setup(solver_num_row);
  columnBFRT.setup(solver_num_row);
  column.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);
  //  row_ap_ultra.setup(solver_num_col);
  columnDensity = 0;
  row_epDensity = 0;
  row_apDensity = 0;
  rowdseDensity = 0;
  // Setup other buffers
  dualRow.setup();
  dualRHS.setup();

  // Initialize for tasks
  if (workHMO.simplex_info_.simplex_strategy == SimplexStrategy::DUAL_TASKS) {
    init_slice(num_threads - 2);
  }

  // Initialize for multi
  if (workHMO.simplex_info_.simplex_strategy == SimplexStrategy::DUAL_MULTI) {
    multi_num = num_threads;
    if (multi_num < 1) multi_num = 1;
    if (multi_num > HIGHS_THREAD_LIMIT) multi_num = HIGHS_THREAD_LIMIT;
    for (int i = 0; i < multi_num; i++) {
      multi_choice[i].row_ep.setup(solver_num_row);
      multi_choice[i].column.setup(solver_num_row);
      multi_choice[i].columnBFRT.setup(solver_num_row);
    }
    init_slice(multi_num - 1);
  }
  multi_iteration = 0;
  //  string partitionFile = model->strOption[STROPT_PARTITION_FILE];
  //  if (partitionFile.size())
  //  {
  //    dualRHS.setup_partition(partitionFile.c_str());
  //  }
}

void HDual::init_slice(int init_sliced_num) {
  // Number of slices
  slice_num = init_sliced_num;
  if (slice_num < 1) slice_num = 1;
  if (slice_num > HIGHS_SLICED_LIMIT) slice_num = HIGHS_SLICED_LIMIT;

  // Alias to the matrix
  const int *Astart = matrix->getAstart();
  const int *Aindex = matrix->getAindex();
  const double *Avalue = matrix->getAvalue();
  const int AcountX = Astart[solver_num_col];

  // Figure out partition weight
  double sliced_countX = AcountX / slice_num;
  slice_start[0] = 0;
  for (int i = 0; i < slice_num - 1; i++) {
    int endColumn = slice_start[i] + 1;  // At least one column
    int endX = Astart[endColumn];
    int stopX = (i + 1) * sliced_countX;
    while (endX < stopX) {
      endX = Astart[++endColumn];
    }
    slice_start[i + 1] = endColumn;
    if (endColumn >= solver_num_col) {
      slice_num = i;  // SHRINK
      break;
    }
  }
  slice_start[slice_num] = solver_num_col;

  // Partition the matrix, row_ap and related packet
  vector<int> sliced_Astart;
  for (int i = 0; i < slice_num; i++) {
    // The matrix
    int mystart = slice_start[i];
    int mycount = slice_start[i + 1] - mystart;
    int mystartX = Astart[mystart];
    sliced_Astart.resize(mycount + 1);
    for (int k = 0; k <= mycount; k++)
      sliced_Astart[k] = Astart[k + mystart] - mystartX;
    // TODO generalise this call so slice can be used with non-logical initial
    // basis
    slice_matrix[i].setup_lgBs(mycount, solver_num_row, &sliced_Astart[0],
                               Aindex + mystartX, Avalue + mystartX);

    // The row_ap and its packages
    slice_row_ap[i].setup(mycount);
    slice_dualRow[i].setupSlice(mycount);
  }
}

void HDual::solve_phase1() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;
  HighsPrintMessage(ML_DETAILED, "dual-phase-1-start\n");
  // Switch to dual phase 1 bounds
  initialise_bound(workHMO, 1); //model->initBound(1);
  initialise_value(workHMO); //  model->initValue();
  // Main solving structure
  timer.start(simplex_info.clock_[IterateClock]);
  for (;;) {
    timer.start(simplex_info.clock_[IterateDualRebuildClock]);
    rebuild();
    timer.stop(simplex_info.clock_[IterateDualRebuildClock]);
    for (;;) {
      switch (simplex_info.simplex_strategy) {
        default:
        case SimplexStrategy::DUAL_PLAIN:
          iterate();
          break;
        case SimplexStrategy::DUAL_TASKS:
          iterate_tasks();
          break;
        case SimplexStrategy::DUAL_MULTI:
          iterate_multi();
          break;
      }
      if (invertHint) break;
      // printf("HDual::solve_phase1: Iter = %d; Objective = %g\n",
      // simplex_info.iteration_count, simplex_info.dual_objective_value);
      double current_dual_objective_value = simplex_info.updatedDualObjectiveValue;
      if (current_dual_objective_value > simplex_info.dual_objective_value_upper_bound) {
#ifdef SCIP_DEV
        printf("HDual::solve_phase1: %12g = Objective > ObjectiveUB\n",
	       current_dual_objective_value, simplex_info.dual_objective_value_upper_bound);
#endif
        simplex_lp_status.solution_status = SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND;
        break;
      }
    }
    double current_run_highs_time = timer.readRunHighsClock();
    if (current_run_highs_time > simplex_info.highs_run_time_limit) {
      SolveBailout = true;
      simplex_lp_status.solution_status = SimplexSolutionStatus::OUT_OF_TIME;
      break;
    }
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) break;
  }

  timer.stop(simplex_info.clock_[IterateClock]);
  if (SolveBailout) return;

  if (rowOut == -1) {
    HighsPrintMessage(ML_DETAILED, "dual-phase-1-optimal\n");
    // Go to phase 2
    if (simplex_info.dualObjectiveValue == 0) {
      solvePhase = 2;
    } else {
      // We still have dual infeasible
      if (workHMO.simplex_info_.costs_perturbed) {
        // Clean up perturbation and go on
        cleanup();
        if (dualInfeasCount == 0) solvePhase = 2;
      } else {
        // Report dual infeasible
        solvePhase = -1;
        HighsPrintMessage(ML_MINIMAL, "dual-infeasible\n");
        simplex_lp_status.solution_status = SimplexSolutionStatus::UNBOUNDED;
      }
    }
  } else if (invertHint == INVERT_HINT_CHOOSE_COLUMN_FAIL) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = -1;
    HighsPrintMessage(ML_MINIMAL, "dual-phase-1-not-solved\n");
    simplex_lp_status.solution_status = SimplexSolutionStatus::FAILED;
  } else if (columnIn == -1) {
    // We got dual phase 1 unbounded - strange
    HighsPrintMessage(ML_MINIMAL, "dual-phase-1-unbounded\n");
    if (workHMO.simplex_info_.costs_perturbed) {
      // Clean up perturbation and go on
      cleanup();
      if (dualInfeasCount == 0) solvePhase = 2;
    } else {
      // Report strange issues
      solvePhase = -1;
      HighsPrintMessage(ML_MINIMAL, "dual-phase-1-not-solved\n");
      simplex_lp_status.solution_status = SimplexSolutionStatus::FAILED;
    }
  }

  if (solvePhase == 2) {
    initialise_bound(workHMO);//    model->initBound();
    initialise_value(workHMO);//    model->initValue();
  }
}

void HDual::solve_phase2() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;
  HighsPrintMessage(ML_DETAILED, "dual-phase-2-start\n");

  // Collect free variables
  dualRow.create_Freelist();
  // Main solving structure
  timer.start(simplex_info.clock_[IterateClock]);
  for (;;) {
    // Outer loop of solve_phase2()
    // Rebuild all values, reinverting B if updates have been performed
    timer.start(simplex_info.clock_[IterateDualRebuildClock]);
    rebuild();
    timer.stop(simplex_info.clock_[IterateDualRebuildClock]);
    if (dualInfeasCount > 0) break;
    for (;;) {
      // Inner loop of solve_phase2()
      // Performs one iteration in case SimplexStrategy::DUAL_PLAIN:
      switch (simplex_info.simplex_strategy) {
        default:
        case SimplexStrategy::DUAL_PLAIN:
          iterate();
          break;
        case SimplexStrategy::DUAL_TASKS:
          iterate_tasks();
          break;
        case SimplexStrategy::DUAL_MULTI:
          iterate_multi();
          break;
      }
      // invertHint can be true for various reasons see SimplexConst.h
      if (invertHint) break;
      double current_dual_objective_value = simplex_info.updatedDualObjectiveValue;
      if (current_dual_objective_value > simplex_info.dual_objective_value_upper_bound) {
#ifdef SCIP_DEV
        printf("HDual::solve_phase2: %12g = Objective > ObjectiveUB\n",
	       current_dual_objective_value, simplex_info.dual_objective_value_upper_bound);
#endif
        simplex_lp_status.solution_status = SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND;
        SolveBailout = true;
        break;
      }
    }
    if (simplex_lp_status.solution_status == SimplexSolutionStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND) {
      SolveBailout = true;
      break;
    }
    double current_run_highs_time = timer.readRunHighsClock();
    if (current_run_highs_time > simplex_info.highs_run_time_limit) {
      simplex_lp_status.solution_status = SimplexSolutionStatus::OUT_OF_TIME;
      SolveBailout = true;
      break;
    }
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) break;
  }
  timer.stop(simplex_info.clock_[IterateClock]);

  if (SolveBailout) {
    return;
  }
  if (dualInfeasCount > 0) {
    // There are dual infeasiblities so switch to Phase 1 and return
    HighsPrintMessage(ML_DETAILED, "dual-phase-2-found-free\n");
    solvePhase = 1;
  } else if (rowOut == -1) {
    // There is no candidate in CHUZR, even after rebuild so probably optimal
    HighsPrintMessage(ML_DETAILED, "dual-phase-2-optimal\n");
    //		printf("Rebuild: cleanup()\n");
    cleanup();
    if (dualInfeasCount > 0) {
      // There are dual infeasiblities after cleanup() so switch to primal
      // simplex
      solvePhase = 4;  // Do primal
    } else {
      // There are no dual infeasiblities after cleanup() so optimal!
      solvePhase = 0;
      HighsPrintMessage(ML_DETAILED, "problem-optimal\n");
      simplex_lp_status.solution_status = SimplexSolutionStatus::OPTIMAL;
    }
  } else if (invertHint == INVERT_HINT_CHOOSE_COLUMN_FAIL) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = -1;
    HighsPrintMessage(ML_MINIMAL, "dual-phase-2-not-solved\n");
    simplex_lp_status.solution_status = SimplexSolutionStatus::FAILED;
  } else if (columnIn == -1) {
    // There is no candidate in CHUZC, so probably dual unbounded
    HighsPrintMessage(ML_MINIMAL, "dual-phase-2-unbounded\n");
    if (workHMO.simplex_info_.costs_perturbed) {
      // If the costs have been perturbed, clean up and return
      cleanup();
    } else {
      // If the costs have not been perturbed, so dual unbounded---and hence
      // primal infeasible
      solvePhase = -1;
      HighsPrintMessage(ML_MINIMAL, "problem-infeasible\n");
      simplex_lp_status.solution_status = SimplexSolutionStatus::INFEASIBLE;
    }
  }
}

void HDual::rebuild() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus &simplex_lp_status = workHMO.simplex_lp_status_;
  // Save history information
  // Move this to Simplex class once it's created
  //  record_pivots(-1, -1, 0);  // Indicate REINVERT

  int sv_invertHint = invertHint;
  invertHint = INVERT_HINT_NO;
  // Possibly Rebuild workHMO.factor_
  bool reInvert = simplex_info.update_count > 0;
  if (!invert_if_row_out_negative) {
    // Don't reinvert if rowOut is negative [equivalently, if sv_invertHint ==
    // INVERT_HINT_POSSIBLY_OPTIMAL]
    if (sv_invertHint == INVERT_HINT_POSSIBLY_OPTIMAL) {
      assert(rowOut == -1);
      reInvert = false;
    }
  }
  if (reInvert) {
    const int *baseIndex = &workHMO.simplex_basis_.basicIndex_[0];
    // Scatter the edge weights so that, after INVERT,
    // they can be gathered according to the new

    // permutation of baseIndex
    timer.start(simplex_info.clock_[PermWtClock]);
    for (int i = 0; i < solver_num_row; i++)
      dualRHS.workEdWtFull[baseIndex[i]] = dualRHS.workEdWt[i];
    timer.stop(simplex_info.clock_[PermWtClock]);

    timer.start(simplex_info.clock_[InvertClock]);

    // Call computeFactor to perform INVERT
    int rankDeficiency = compute_factor(workHMO); //    int rankDeficiency = model->computeFactor();
    timer.stop(simplex_info.clock_[InvertClock]);

    if (rankDeficiency)
      throw runtime_error("Dual reInvert: singular-basis-matrix");
    // Gather the edge weights according to the
    // permutation of baseIndex after INVERT
    timer.start(simplex_info.clock_[PermWtClock]);
    for (int i = 0; i < solver_num_row; i++)
      dualRHS.workEdWt[i] = dualRHS.workEdWtFull[baseIndex[i]];
    timer.stop(simplex_info.clock_[PermWtClock]);

    // Possibly look at the basis condition
    //		double bsCond = an_bs_cond();
  }

  // Recompute dual solution
  timer.start(simplex_info.clock_[ComputeDualClock]);
  compute_dual(workHMO);  //  model->computeDual();
  timer.stop(simplex_info.clock_[ComputeDualClock]);

  timer.start(simplex_info.clock_[CorrectDualClock]);
  correct_dual(workHMO, &dualInfeasCount);  //model->correctDual(&dualInfeasCount);
  timer.stop(simplex_info.clock_[CorrectDualClock]);

  // Recompute primal solution
  timer.start(simplex_info.clock_[ComputePrimalClock]);
  compute_primal(workHMO);//model->computePrimal();
  timer.stop(simplex_info.clock_[ComputePrimalClock]);

  // Collect primal infeasible as a list
  timer.start(simplex_info.clock_[CollectPrIfsClock]);
  dualRHS.create_infeasArray();
  dualRHS.create_infeasList(columnDensity);
  timer.stop(simplex_info.clock_[CollectPrIfsClock]);

  // Check the objective value maintained by updating against the
  // value when computed exactly - so long as there is a value to
  // check against
  bool checkDualObjectiveValue = simplex_lp_status.has_dual_objective_value;
  // Compute the objective value
  timer.start(simplex_info.clock_[ComputeDuobjClock]);
  compute_dual_objective_value(workHMO, solvePhase);
  timer.stop(simplex_info.clock_[ComputeDuobjClock]);

  double dualObjectiveValue = simplex_info.dualObjectiveValue;
  if (checkDualObjectiveValue) {
    double absDualObjectiveError = fabs(simplex_info.updatedDualObjectiveValue - dualObjectiveValue);
    double rlvDualObjectiveError = absDualObjectiveError/max(1.0, fabs(dualObjectiveValue));
#ifdef HiGHSDEV
    /*
    // TODO Investigate these Dual objective value errors
    if (rlvDualObjectiveError >= 1e-8) {
      HighsLogMessage(HighsMessageType::WARNING, "Dual objective value error abs(rel) = %12g (%12g)",
			absDualObjectiveError, rlvDualObjectiveError);
    }
    */
#endif
  }
  simplex_info.updatedDualObjectiveValue = dualObjectiveValue;

#ifdef HiGHSDEV
  //  checkDualObjectiveValue("After computing dual objective value");
  //  printf("Checking INVERT in rebuild()\n"); workHMO.factor_.checkInvert();
#endif

  timer.start(simplex_info.clock_[ReportInvertClock]);
  iterateRpInvert(sv_invertHint);
  timer.stop(simplex_info.clock_[ReportInvertClock]);

  total_INVERT_TICK = factor->build_syntheticTick;  // Was factor->pseudoTick
  total_FT_inc_TICK = 0;
#ifdef HiGHSDEV
  total_fake = 0;
#endif
  total_syntheticTick = 0;

#ifdef HiGHSDEV
  if (simplex_info.analyseRebuildTime) {
    int iClock = simplex_info.clock_[IterateDualRebuildClock];
    int totalRebuilds = timer.clock_num_call[iClock];
    double totalRebuildTime = timer.read(iClock);
    printf(
	   "Dual  Ph%-2d rebuild %4d (%1d) on iteration %9d: Total rebuild time = %11.4g\n",
	   solvePhase, totalRebuilds, sv_invertHint, simplex_info.iteration_count, totalRebuildTime);
  }
#endif
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HDual::cleanup() {
  // Remove perturbation and recompute the dual solution
  HighsPrintMessage(ML_DETAILED, "dual-cleanup-shift\n");
  initialise_cost(workHMO);  //  model->initCost();
  initialise_bound(workHMO);  //  model->initBound();
  compute_dual(workHMO);  //  model->computeDual();
  compute_dual_objective_value(workHMO, solvePhase);
  iterateRpInvert(-1);

  compute_dual_infeasible_in_primal(workHMO, &dualInfeasCount);//model->computeDualInfeasInPrimal(&dualInfeasCount);
}

void HDual::iterate() {
  // This is the main teration loop for dual revised simplex. All the
  // methods have as their first line if (invertHint) return;, where
  // invertHint is, for example, set to 1 when CHUZR finds no
  // candidate. This causes a break from the inner loop of
  // solve_phase% and, hence, a call to rebuild()

  // Reporting:
  // Row-wise matrix after update in updateMatrix(columnIn, columnOut);
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  timer.start(simplex_info.clock_[IterateChuzrClock]);
  chooseRow();
  timer.stop(simplex_info.clock_[IterateChuzrClock]);

  timer.start(simplex_info.clock_[IterateChuzcClock]);
  chooseColumn(&row_ep);
  timer.stop(simplex_info.clock_[IterateChuzcClock]);

  timer.start(simplex_info.clock_[IterateFtranClock]);
  updateFtranBFRT();
  // updateFtran(); computes the pivotal column in the data structure "column"
  updateFtran();

  // updateFtranDSE performs the DSE FTRAN on pi_p
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) updateFtranDSE(&row_ep);
  timer.stop(simplex_info.clock_[IterateFtranClock]);

  // updateVerify() Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  timer.start(simplex_info.clock_[IterateVerifyClock]);
  updateVerify();
  timer.stop(simplex_info.clock_[IterateVerifyClock]);

  // updateDual() Updates the dual values
  timer.start(simplex_info.clock_[IterateDualClock]);
  updateDual();
  timer.stop(simplex_info.clock_[IterateDualClock]);

  // updatePrimal(&row_ep); Updates the primal values and the edge weights
  timer.start(simplex_info.clock_[IteratePrimalClock]);
  updatePrimal(&row_ep);
  timer.stop(simplex_info.clock_[IteratePrimalClock]);

  if ((dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) && (nw_dvx_fwk)) {
    timer.start(simplex_info.clock_[IterateDevexIzClock]);
    iz_dvx_fwk();
    timer.stop(simplex_info.clock_[IterateDevexIzClock]);
  }

  // Update the basis representation
  timer.start(simplex_info.clock_[IteratePivotsClock]);
  updatePivots();
  timer.stop(simplex_info.clock_[IteratePivotsClock]);

  // Analyse the iteration: possibly report; possibly switch strategy
  iterateAn();
}

void HDual::iterate_tasks() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  slice_PRICE = 1;

  // Group 1
  chooseRow();

  // Disable slice when too sparse
  if (1.0 * row_ep.count / solver_num_row < 0.01) slice_PRICE = 0;

  timer.start(simplex_info.clock_[Group1Clock]);
#pragma omp parallel
#pragma omp single
  {
#pragma omp task
    {
      columnDSE.copy(&row_ep);
      updateFtranDSE(&columnDSE);
    }
#pragma omp task
    {
      if (slice_PRICE)
        chooseColumn_slice(&row_ep);
      else
        chooseColumn(&row_ep);
#pragma omp task
      updateFtranBFRT();
#pragma omp task
      updateFtran();
#pragma omp taskwait
    }
  }
  timer.stop(simplex_info.clock_[Group1Clock]);

  updateVerify();
  updateDual();
  updatePrimal(&columnDSE);
  updatePivots();
}

void HDual::iterateIzAn() {
  HighsTimer &timer = workHMO.timer_;
  AnIterIt0 = workHMO.simplex_info_.iteration_count;
  AnIterCostlyDseFq = 0;
#ifdef HiGHSDEV
  AnIterPrevRpNumCostlyDseIt = 0;
  AnIterPrevIt = 0;
  AnIterOpRec *AnIter;
  AnIter = &AnIterOp[AnIterOpTy_Btran];
  AnIter->AnIterOpName = "Btran";
  AnIter = &AnIterOp[AnIterOpTy_Price];
  AnIter->AnIterOpName = "Price";
  AnIter = &AnIterOp[AnIterOpTy_Ftran];
  AnIter->AnIterOpName = "Ftran";
  AnIter = &AnIterOp[AnIterOpTy_FtranBFRT];
  AnIter->AnIterOpName = "FtranBFRT";
  AnIter = &AnIterOp[AnIterOpTy_FtranDSE];
  AnIter->AnIterOpName = "FtranDSE";
  for (int k = 0; k < NumAnIterOpTy; k++) {
    AnIter = &AnIterOp[k];
    AnIter->AnIterOpLog10RsDsty = 0;
    AnIter->AnIterOpSuLog10RsDsty = 0;
    if (k == AnIterOpTy_Price) {
      AnIter->AnIterOpHyperCANCEL = 1.0;
      AnIter->AnIterOpHyperTRAN = 1.0;
      AnIter->AnIterOpRsDim = solver_num_col;
    } else {
      if (k == AnIterOpTy_Btran) {
        AnIter->AnIterOpHyperCANCEL = hyperCANCEL;
        AnIter->AnIterOpHyperTRAN = hyperBTRANU;
      } else {
        AnIter->AnIterOpHyperCANCEL = hyperCANCEL;
        AnIter->AnIterOpHyperTRAN = hyperFTRANL;
      }
      AnIter->AnIterOpRsDim = solver_num_row;
    }
    AnIter->AnIterOpNumCa = 0;
    AnIter->AnIterOpNumHyperOp = 0;
    AnIter->AnIterOpNumHyperRs = 0;
    AnIter->AnIterOpRsMxNNZ = 0;
    AnIter->AnIterOpSuNumCa = 0;
    AnIter->AnIterOpSuNumHyperOp = 0;
    AnIter->AnIterOpSuNumHyperRs = 0;
  }
  int last_invert_hint = INVERT_HINT_Count-1;
  for (int k = 1; k <= last_invert_hint; k++) AnIterNumInvert[k] = 0;
  AnIterNumPrDgnIt = 0;
  AnIterNumDuDgnIt = 0;
  AnIterNumColPrice = 0;
  AnIterNumRowPrice = 0;
  AnIterNumRowPriceWSw = 0;
  AnIterNumRowPriceUltra = 0;
  int last_dual_edge_weight_mode = (int) DualEdgeWeightMode::STEEPEST_EDGE;
  for (int k = 0; k <= last_dual_edge_weight_mode; k++) {
    AnIterNumEdWtIt[k] = 0;
  }
  AnIterNumCostlyDseIt = 0;
  AnIterTraceNumRec = 0;
  AnIterTraceIterDl = 1;
  AnIterTraceRec *lcAnIter = &AnIterTrace[0];
  lcAnIter->AnIterTraceIter = AnIterIt0;
  lcAnIter->AnIterTraceTime = timer.getTime();
#endif
}

void HDual::iterateAn() {
  HighsTimer &timer = workHMO.timer_;
  // Possibly report on the iteration
  iterateRp();

  // Possibly switch from DSE to Dvx
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    double AnIterCostlyDseMeasureDen;
    //    AnIterCostlyDseMeasureDen = row_epDensity*columnDensity;
    AnIterCostlyDseMeasureDen =
        max(max(row_epDensity, columnDensity), row_apDensity);
    if (AnIterCostlyDseMeasureDen > 0) {
      AnIterCostlyDseMeasure = rowdseDensity / AnIterCostlyDseMeasureDen;
      AnIterCostlyDseMeasure = AnIterCostlyDseMeasure * AnIterCostlyDseMeasure;
    } else {
      AnIterCostlyDseMeasure = 0;
    }
    bool CostlyDseIt = AnIterCostlyDseMeasure > AnIterCostlyDseMeasureLimit &&
                       rowdseDensity > AnIterCostlyDseMnDensity;
    AnIterCostlyDseFq = (1 - runningAverageMu) * AnIterCostlyDseFq;
    if (CostlyDseIt) {
      AnIterNumCostlyDseIt++;
      AnIterCostlyDseFq += runningAverageMu * 1.0;
      int lcNumIter = workHMO.simplex_info_.iteration_count - AnIterIt0;
      if (allow_dual_steepest_edge_to_devex_switch &&
          (AnIterNumCostlyDseIt > lcNumIter * AnIterFracNumCostlyDseItbfSw) &&
          (lcNumIter > AnIterFracNumTot_ItBfSw * solver_num_tot)) {
        // At least 5% of the (at least) 0.1NumTot iterations have been costly
        // DSE so switch to Devex
#ifdef HiGHSDEV
        printf(
            "Switch from DSE to Dvx after %d costly DSE iterations of %d: "
            "Col_Dsty = %11.4g; R_Ep_Dsty = %11.4g; DSE_Dsty = %11.4g\n",
            AnIterNumCostlyDseIt, lcNumIter, rowdseDensity, row_epDensity,
            columnDensity);
#endif
        dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
        // Zero the number of Devex frameworks used and set up the first one
        n_dvx_fwk = 0;
        dvx_ix.assign(solver_num_tot, 0);
        iz_dvx_fwk();
      }
    }
  }

#ifdef HiGHSDEV
  int AnIterCuIt = workHMO.simplex_info_.iteration_count;
  bool iterLg = AnIterCuIt % 100 == 0;
  iterLg = false;
  if (iterLg) {
    int lc_NumCostlyDseIt = AnIterNumCostlyDseIt - AnIterPrevRpNumCostlyDseIt;
    AnIterPrevRpNumCostlyDseIt = AnIterNumCostlyDseIt;
    printf("Iter %10d: ", AnIterCuIt);
    iterateRpDsty(ML_MINIMAL, true);
    iterateRpDsty(ML_MINIMAL, false);
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      int lc_pct = (100 * AnIterNumCostlyDseIt) / (AnIterCuIt - AnIterIt0);
      printf("| Fq = %4.2f; Su =%5d (%3d%%)", AnIterCostlyDseFq,
             AnIterNumCostlyDseIt, lc_pct);

      if (lc_NumCostlyDseIt > 0) printf("; LcNum =%3d", lc_NumCostlyDseIt);
    }
    printf("\n");
  }

  for (int k = 0; k < NumAnIterOpTy; k++) {
    AnIterOpRec *lcAnIterOp = &AnIterOp[k];
    if (lcAnIterOp->AnIterOpNumCa) {
      lcAnIterOp->AnIterOpSuNumCa += lcAnIterOp->AnIterOpNumCa;
      lcAnIterOp->AnIterOpSuNumHyperOp += lcAnIterOp->AnIterOpNumHyperOp;
      lcAnIterOp->AnIterOpSuNumHyperRs += lcAnIterOp->AnIterOpNumHyperRs;
      lcAnIterOp->AnIterOpSuLog10RsDsty += lcAnIterOp->AnIterOpLog10RsDsty;
    }
    lcAnIterOp->AnIterOpNumCa = 0;
    lcAnIterOp->AnIterOpNumHyperOp = 0;
    lcAnIterOp->AnIterOpNumHyperRs = 0;
    lcAnIterOp->AnIterOpLog10RsDsty = 0;
  }
  if (invertHint > 0) AnIterNumInvert[invertHint]++;
  if (thetaDual <= 0) AnIterNumDuDgnIt++;
  if (thetaPrimal <= 0) AnIterNumPrDgnIt++;
  if (AnIterCuIt > AnIterPrevIt)
    AnIterNumEdWtIt[(int) dual_edge_weight_mode] += (AnIterCuIt - AnIterPrevIt);

  AnIterTraceRec *lcAnIter = &AnIterTrace[AnIterTraceNumRec];
  //  if (workHMO.simplex_info_.iteration_count ==
  //  AnIterTraceIterRec[AnIterTraceNumRec]+AnIterTraceIterDl) {
  if (workHMO.simplex_info_.iteration_count == lcAnIter->AnIterTraceIter + AnIterTraceIterDl) {
    if (AnIterTraceNumRec == AN_ITER_TRACE_MX_NUM_REC) {
      for (int rec = 1; rec <= AN_ITER_TRACE_MX_NUM_REC / 2; rec++)
        AnIterTrace[rec] = AnIterTrace[2 * rec];
      AnIterTraceNumRec = AnIterTraceNumRec / 2;
      AnIterTraceIterDl = AnIterTraceIterDl * 2;
    } else {
      AnIterTraceNumRec++;
      lcAnIter = &AnIterTrace[AnIterTraceNumRec];
      lcAnIter->AnIterTraceIter = workHMO.simplex_info_.iteration_count;
      lcAnIter->AnIterTraceTime = timer.getTime();
      lcAnIter->AnIterTraceDsty[AnIterOpTy_Btran] = row_epDensity;
      lcAnIter->AnIterTraceDsty[AnIterOpTy_Price] = row_apDensity;
      lcAnIter->AnIterTraceDsty[AnIterOpTy_Ftran] = columnDensity;
      lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranBFRT] = columnDensity;
      if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
        lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = rowdseDensity;
        lcAnIter->AnIterTraceAux0 = AnIterCostlyDseMeasure;
      } else {
        lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = 0;
        lcAnIter->AnIterTraceAux0 = 0;
      }
      lcAnIter->AnIterTrace_dual_edge_weight_mode = (int) dual_edge_weight_mode;
    }
  }
  AnIterPrevIt = AnIterCuIt;
#endif
}

void HDual::iterateRp() {
  int numIter = workHMO.simplex_info_.iteration_count;
  bool header = numIter % 10 == 1;
  //  header = true;  // JAJH10/10
  if (header) iterateRpFull(header);
  iterateRpFull(false);
}

void HDual::iterateRpFull(bool header) {
  if (header) {
    iterateRpIterPh(ML_DETAILED, true);
    iterateRpDuObj(ML_DETAILED, true);
#ifdef HiGHSDEV
    iterateRpIterDa(ML_DETAILED, true);
    iterateRpDsty(ML_DETAILED, true);
    HighsPrintMessage(ML_DETAILED, " FreeLsZ");
#endif
    HighsPrintMessage(ML_DETAILED, "\n");
  } else {
    iterateRpIterPh(ML_DETAILED, false);
    iterateRpDuObj(ML_DETAILED, false);
#ifdef HiGHSDEV
    iterateRpIterDa(ML_DETAILED, false);
    iterateRpDsty(ML_DETAILED, false);
    HighsPrintMessage(ML_DETAILED, " %7d", dualRow.freeListSize);
#endif
    HighsPrintMessage(ML_DETAILED, "\n");
  }
}

void HDual::iterateRpIterPh(int iterate_log_level, bool header) {
  if (header) {
    HighsPrintMessage(iterate_log_level, " Iteration Ph");
  } else {
    int numIter = workHMO.simplex_info_.iteration_count;
    HighsPrintMessage(iterate_log_level, " %9d %2d", numIter, solvePhase);
  }
}
void HDual::iterateRpDuObj(int iterate_log_level, bool header) {
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  if (header) {
    HighsPrintMessage(iterate_log_level, "    DualObjective    ");
  } else {
    HighsPrintMessage(iterate_log_level, " %20.10e", simplex_info.updatedDualObjectiveValue);
  }
}

void HDual::iterateRpIterDa(int iterate_log_level, bool header) {
  if (header) {
    HighsPrintMessage(iterate_log_level, " Inv       NumCk     LvR     LvC     EnC        DlPr        ThDu        ThPr          Aa");
  } else {
    HighsPrintMessage(iterate_log_level, " %3d %11.4g %7d %7d %7d %11.4g %11.4g %11.4g %11.4g", 
		      invertHint, numericalTrouble, rowOut, columnOut, columnIn, deltaPrimal,
		      thetaDual, thetaPrimal, alpha);
  }
}

void HDual::iterateRpDsty(int iterate_log_level, bool header) {
  if (header) {
    HighsPrintMessage(iterate_log_level, "  Col R_Ep R_Ap  DSE");
  } else {
    int l10ColDse = intLog10(columnDensity);
    int l10REpDse = intLog10(row_epDensity);
    int l10RapDse = intLog10(row_apDensity);
    double lc_rowdseDensity;
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      lc_rowdseDensity = rowdseDensity;
    } else {
      lc_rowdseDensity = 0;
    }
    int l10DseDse = intLog10(lc_rowdseDensity);
    HighsPrintMessage(iterate_log_level, " %4d %4d %4d %4d", l10ColDse, l10REpDse, l10RapDse, l10DseDse);
  }
}

void HDual::iterateRpInvert(int i_v) {
#ifdef HiGHSDEV
  HighsPrintMessage(ML_MINIMAL, "Iter %10d:", workHMO.simplex_info_.iteration_count);
  iterateRpDsty(ML_MINIMAL, true);
  iterateRpDsty(ML_MINIMAL, false);
  iterateRpDuObj(ML_MINIMAL, false);
  HighsPrintMessage(ML_MINIMAL, " %2d\n", i_v);
#else
  report_iteration_count_dual_objective_value(workHMO, i_v);
#endif
}

int HDual::intLog10(double v) {
  int intLog10V = -99;
  if (v > 0) intLog10V = log(v) / log(10.0);
  return intLog10V;
}

void HDual::uOpRsDensityRec(double lc_OpRsDensity, double &opRsDensity) {
  // Update an average density record for BTRAN, an FTRAN or PRICE
  opRsDensity = (1 - runningAverageMu) * (opRsDensity) +
                runningAverageMu * lc_OpRsDensity;
}

void HDual::chooseRow() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  // Choose the index of a row to leave the basis (CHUZR)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  // Choose candidates repeatedly until candidate is OK or optimality is
  // detected
  for (;;) {
    // Choose the index of a good row to leave the basis
    dualRHS.choose_normal(&rowOut);
    if (rowOut == -1) {
      // No index found so may be dual optimal. By setting
      // invertHint>0 all subsequent methods in the iteration will
      // be skipped until reinversion and rebuild have taken place
      invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
      return;
    }
    // Compute pi_p = B^{-T}e_p in row_ep
    timer.start(simplex_info.clock_[BtranClock]);
    // Set up RHS for BTRAN
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = rowOut;
    row_ep.array[rowOut] = 1;
    row_ep.packFlag = true;
#ifdef HiGHSDEV
    if (simplex_info.analyseSimplexIterations) iterateOpRecBf(AnIterOpTy_Btran, row_ep, row_epDensity);
#endif
    // Perform BTRAN
    factor->btran(row_ep, row_epDensity);
#ifdef HiGHSDEV
    if (simplex_info.analyseSimplexIterations) iterateOpRecAf(AnIterOpTy_Btran, row_ep);
#endif
    timer.stop(simplex_info.clock_[BtranClock]);
    // Verify DSE weight
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // For DSE, see how accurate the updated weight is
      // Save the updated weight
      double u_weight = dualRHS.workEdWt[rowOut];
      // Compute the weight from row_ep and over-write the updated weight
      double c_weight = dualRHS.workEdWt[rowOut] = row_ep.norm2();
      // If the weight error is acceptable then break out of the
      // loop. All we worry about is accepting rows with weights
      // which are not too small, since this can make the row look
      // unreasonably attractive
      if (u_weight >= 0.25 * c_weight) break;
#ifdef HiGHSDEV
      // Count the number of wrong DSE weights for reporting
      n_wg_DSE_wt += 1;
#endif
      // Weight error is unacceptable so look for another
      // candidate. Of course, it's possible that the same
      // candidate is chosen, but the weight will be correct (so
      // no infinite loop).
    } else {
      // If not using DSE then accept the row by breaking out of
      // the loop
      break;
    }
  }
  // Index of row to leave the basis has been found
  //
  // Assign basic info:
  //
  // Record the column (variable) associated with the leaving row
  columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  // Record the change in primal variable associated with the move to the bound
  // being violated
  if (baseValue[rowOut] < baseLower[rowOut]) {
    // Below the lower bound so set deltaPrimal = value - LB < 0
    deltaPrimal = baseValue[rowOut] - baseLower[rowOut];
  } else {
    // Above the upper bound so set deltaPrimal = value - UB > 0
    deltaPrimal = baseValue[rowOut] - baseUpper[rowOut];
  }
  // Set sourceOut to be -1 if deltaPrimal<0, otherwise +1 (since deltaPrimal>0)
  sourceOut = deltaPrimal < 0 ? -1 : 1;
  // Update the record of average row_ep (pi_p) density. This ignores
  // any BTRANs done for skipped candidates
  double lc_OpRsDensity = (double)row_ep.count / solver_num_row;
  uOpRsDensityRec(lc_OpRsDensity, row_epDensity);
}

void HDual::chooseColumn(HVector *row_ep) {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  // Compute pivot row (PRICE) and choose the index of a column to enter the
  // basis (CHUZC)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  //
  // PRICE
  //
  timer.start(simplex_info.clock_[PriceClock]);
  row_ap.clear();

#ifdef HiGHSDEV
  bool anPriceEr = false;
  bool useUltraPrice = allow_price_ultra &&
                       row_apDensity * solver_num_col * 10 < row_ap.ilP2 &&
                       row_apDensity < 1e-3;
#endif
  if (price_mode == PriceMode::COL) {
    // Column-wise PRICE
#ifdef HiGHSDEV
    if (simplex_info.analyseSimplexIterations) {
      iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
      AnIterNumColPrice++;
    }
#endif
    // Perform column-wise PRICE
    matrix->price_by_col(row_ap, *row_ep);
  } else {
    // By default, use row-wise PRICE, but possibly use column-wise
    // PRICE if the density of row_ep is too high
    double lc_dsty = (double)(*row_ep).count / solver_num_row;
    if (allow_price_by_col_switch && (lc_dsty > dstyColPriceSw)) {
      // Use column-wise PRICE due to density of row_ep
#ifdef HiGHSDEV
      if (simplex_info.analyseSimplexIterations) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
        AnIterNumColPrice++;
      }
#endif
      // Perform column-wise PRICE
      matrix->price_by_col(row_ap, *row_ep);
      // Zero the components of row_ap corresponding to basic variables
      // (nonbasicFlag[*]=0)
      const int *nonbasicFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
      for (int col = 0; col < solver_num_col; col++) {
        row_ap.array[col] = nonbasicFlag[col] * row_ap.array[col];
      }
#ifdef HiGHSDEV
      // Ultra-sparse PRICE is in development
    } else if (useUltraPrice) {
      if (simplex_info.analyseSimplexIterations) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, row_apDensity);
        AnIterNumRowPriceUltra++;
      }
      // Perform ultra-sparse row-wise PRICE
      matrix->price_by_row_ultra(row_ap, *row_ep);
      if (anPriceEr) {
        bool price_er;
        price_er = matrix->price_er_ck(row_ap, *row_ep);
        if (!price_er) HighsPrintMessage(ML_VERBOSE, "No ultra PRICE error\n");
      }
#endif
    } else if (allow_price_by_row_switch) {
      // Avoid hyper-sparse PRICE on current density of result or
      // switch if the density of row_ap becomes extreme
#ifdef HiGHSDEV
      if (simplex_info.analyseSimplexIterations) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, row_apDensity);
        AnIterNumRowPriceWSw++;
      }
#endif
      // Set the value of the density of row_ap at which the switch to
      // sparse row-wise PRICE should be made
      const double sw_dsty = matrix->price_by_row_sw_dsty;
      // Perform hyper-sparse row-wise PRICE with switching
      matrix->price_by_row_w_sw(row_ap, *row_ep, row_apDensity, 0, sw_dsty);
    } else {
      // No avoiding hyper-sparse PRICE on current density of result
      // or switch if the density of row_ap becomes extreme
#ifdef HiGHSDEV
      if (simplex_info.analyseSimplexIterations) {
        iterateOpRecBf(AnIterOpTy_Price, *row_ep, 0.0);
        AnIterNumRowPrice++;
      }
#endif
      // Perform hyper-sparse row-wise PRICE
      matrix->price_by_row(row_ap, *row_ep);
    }
  }
#ifdef HiGHSDEV
  // Possibly analyse the error in the result of PRICE
  if (anPriceEr) {
    matrix->price_er_ck(row_ap, *row_ep);
  }
#endif
  // Update the record of average row_ap density
  double lc_OpRsDensity = (double)row_ap.count / solver_num_col;
  uOpRsDensityRec(lc_OpRsDensity, row_apDensity);
#ifdef HiGHSDEV
  if (simplex_info.analyseSimplexIterations) iterateOpRecAf(AnIterOpTy_Price, row_ap);
#endif
  timer.stop(simplex_info.clock_[PriceClock]);
  //
  // CHUZC
  //
  // Section 0: Clear data and call create_Freemove to set a value of
  // nonbasicMove for all free columns to prevent their dual values
  // from being changed.
  timer.start(simplex_info.clock_[Chuzc0Clock]);
  dualRow.clear();
  dualRow.workDelta = deltaPrimal;
  dualRow.create_Freemove(row_ep);
  timer.stop(simplex_info.clock_[Chuzc0Clock]);
  //
  // Section 1: Pack row_ap and row_ep, then determine the possible
  // variables - candidates for CHUZC
  timer.start(simplex_info.clock_[Chuzc1Clock]);
  dualRow.choose_makepack(
      &row_ap, 0);  // Pack row_ap into the packIndex/Value of HDualRow
  dualRow.choose_makepack(
      row_ep, solver_num_col);  // Pack row_ep into the packIndex/Value of HDualRow
  dualRow.choose_possible();  // Determine the possible variables - candidates
                              // for CHUZC
  timer.stop(simplex_info.clock_[Chuzc1Clock]);
  //
  // Take action if the step to an expanded bound is not positive, or
  // there are no candidates for CHUZC
  columnIn = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    invertHint = INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED;  
    return;
  }
  //
  // Sections 2 and 3: Perform (bound-flipping) ratio test. This can
  // fail if the dual values are excessively large
  bool chooseColumnFail = dualRow.choose_final();
  if (chooseColumnFail) {
    invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
    return;
  }
  //
  // Section 4: Reset the nonbasicMove values for free columns
  timer.start(simplex_info.clock_[Chuzc4Clock]);
  dualRow.delete_Freemove();
  timer.stop(simplex_info.clock_[Chuzc4Clock]);
  // Record values for basis change, checking for numerical problems and update
  // of dual variables
  columnIn = dualRow.workPivot;   // Index of the column entering the basis
  alphaRow = dualRow.workAlpha;   // Pivot value computed row-wise - used for
                                  // numerical checking
  thetaDual = dualRow.workTheta;  // Dual step length

  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    timer.start(simplex_info.clock_[DevexWtClock]);
    // Determine the exact Devex weight
    double og_dvx_wt_o_rowOut = dualRHS.workEdWt[rowOut];
    double tru_dvx_wt_o_rowOut = 0;
    // Loop over [row_ap; row_ep] using the packed values
    for (int el_n = 0; el_n < dualRow.packCount; el_n++) {
      int vr_n = dualRow.packIndex[el_n];
      double pv = dvx_ix[vr_n] * dualRow.packValue[el_n];
      tru_dvx_wt_o_rowOut += pv * pv;
    }
    tru_dvx_wt_o_rowOut = max(1.0, tru_dvx_wt_o_rowOut);
    // Analyse the Devex weight to determine whether a new framework
    // should be set up
    double dvx_rao = max(og_dvx_wt_o_rowOut / tru_dvx_wt_o_rowOut,
                         tru_dvx_wt_o_rowOut / og_dvx_wt_o_rowOut);
    int i_te = solver_num_row / minRlvNumberDevexIterations;
    i_te = max(minAbsNumberDevexIterations, i_te);
    // Square maxAllowedDevexWeightRatio due to keeping squared
    // weights
    nw_dvx_fwk =
        dvx_rao > maxAllowedDevexWeightRatio * maxAllowedDevexWeightRatio ||
        n_dvx_it > i_te;
    dualRHS.workEdWt[rowOut] = tru_dvx_wt_o_rowOut;
    timer.stop(simplex_info.clock_[DevexWtClock]);
  }
  return;
}

void HDual::chooseColumn_slice(HVector *row_ep) {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  // Choose the index of a column to enter the basis (CHUZC) by
  // exploiting slices of the pivotal row - for SIP and PAMI
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  timer.start(simplex_info.clock_[Chuzr1Clock]);
  dualRow.clear();
  dualRow.workDelta = deltaPrimal;
  dualRow.create_Freemove(row_ep);

  // Row_ep:         PACK + CC1
#pragma omp task
  {
    dualRow.choose_makepack(row_ep, solver_num_col);
    dualRow.choose_possible();
  }

  // Row_ap: PRICE + PACK + CC1
  for (int i = 0; i < slice_num; i++) {
#pragma omp task
    {
      slice_row_ap[i].clear();
      slice_matrix[i].price_by_row(slice_row_ap[i], *row_ep);

      slice_dualRow[i].clear();
      slice_dualRow[i].workDelta = deltaPrimal;
      slice_dualRow[i].choose_makepack(&slice_row_ap[i], slice_start[i]);
      slice_dualRow[i].choose_possible();
    }
  }
#pragma omp taskwait

  // Join CC1 results here
  for (int i = 0; i < slice_num; i++)
    dualRow.choose_joinpack(&slice_dualRow[i]);

  // Infeasible we created before
  columnIn = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    invertHint = INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED;  
    timer.stop(simplex_info.clock_[Chuzr1Clock]);
    return;
  }
  timer.stop(simplex_info.clock_[Chuzr1Clock]);

  // Choose column 2, This only happens if didn't go out
  dualRow.choose_final();
  dualRow.delete_Freemove();
  columnIn = dualRow.workPivot;
  alphaRow = dualRow.workAlpha;
  thetaDual = dualRow.workTheta;
}

void HDual::updateFtran() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  // Compute the pivotal column (FTRAN)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  timer.start(simplex_info.clock_[FtranClock]);
  // Clear the picotal column and indicate that its values should be packed
  column.clear();
  column.packFlag = true;
  // Get the constraint matrix column by combining just one column
  // with unit multiplier
  matrix->collect_aj(column, columnIn, 1);
#ifdef HiGHSDEV
  if (simplex_info.analyseSimplexIterations) iterateOpRecBf(AnIterOpTy_Ftran, column, columnDensity);
#endif
  // Perform FTRAN
  factor->ftran(column, columnDensity);
#ifdef HiGHSDEV
  if (simplex_info.analyseSimplexIterations) iterateOpRecAf(AnIterOpTy_Ftran, column);
#endif
  // Save the pivot value computed column-wise - used for numerical checking
  alpha = column.array[rowOut];
  timer.stop(simplex_info.clock_[FtranClock]);
}

void HDual::updateFtranBFRT() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  // Compute the RHS changes corresponding to the BFRT (FTRAN-BFRT)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Only time updateFtranBFRT if dualRow.workCount > 0;
  // If dualRow.workCount = 0 then dualRow.update_flip(&columnBFRT)
  // merely clears columnBFRT so no FTRAN is performed
  bool time_updateFtranBFRT = dualRow.workCount > 0;

  if (time_updateFtranBFRT) {
    timer.start(simplex_info.clock_[FtranBfrtClock]);
  }

  dualRow.update_flip(&columnBFRT);

  if (columnBFRT.count) {
#ifdef HiGHSDEV
    if (simplex_info.analyseSimplexIterations)
      iterateOpRecBf(AnIterOpTy_FtranBFRT, columnBFRT, columnDensity);
#endif
    // Perform FTRAN BFRT
    factor->ftran(columnBFRT, columnDensity);
#ifdef HiGHSDEV
    if (simplex_info.analyseSimplexIterations) iterateOpRecAf(AnIterOpTy_FtranBFRT, columnBFRT);
#endif
  }
  if (time_updateFtranBFRT) {
    timer.stop(simplex_info.clock_[FtranBfrtClock]);
  }
}

void HDual::updateFtranDSE(HVector *DSE_Vector) {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  // Compute the vector required to update DSE weights - being FTRAN
  // applied to the pivotal column (FTRAN-DSE)
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  timer.start(simplex_info.clock_[FtranDseClock]);
#ifdef HiGHSDEV
  if (simplex_info.analyseSimplexIterations) iterateOpRecBf(AnIterOpTy_FtranDSE, *DSE_Vector, rowdseDensity);
#endif
  // Perform FTRAN DSE
  factor->ftran(*DSE_Vector, rowdseDensity);
#ifdef HiGHSDEV
  if (simplex_info.analyseSimplexIterations) iterateOpRecAf(AnIterOpTy_FtranDSE, *DSE_Vector);
#endif
  timer.stop(simplex_info.clock_[FtranDseClock]);
}

void HDual::updateVerify() {
  // Compare the pivot value computed row-wise and column-wise and
  // determine whether reinversion is advisable
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Look at the relative difference between the absolute values of the two
  // pivot values
  double aCol = fabs(alpha);
  double aRow = fabs(alphaRow);
  double aDiff = fabs(aCol - aRow);
  numericalTrouble = aDiff / min(aCol, aRow);
  // Reinvert if the relative difference is large enough, and updates hav ebeen
  // performed
  if (numericalTrouble > 1e-7 && workHMO.simplex_info_.update_count > 0) {
    invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
  }
}

void HDual::updateDual() {
  // Update the dual values
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;

  // Update - dual (shift and back)
  if (thetaDual == 0)
    // Little to do if thetaDual is zero
    shift_cost(workHMO, columnIn, -workDual[columnIn]);//model->shiftCost(columnIn, -workDual[columnIn]);
  else {
    // Update the dual values (if packCount>0)
    dualRow.update_dual(thetaDual, columnOut);
    if (workHMO.simplex_info_.simplex_strategy != SimplexStrategy::DUAL_PLAIN && slice_PRICE) {
      // Update the dual variables slice-by-slice [presumably
      // nothing is done in the previous call to
      // dualRow.update_dual. TODO: Check with Qi
      for (int i = 0; i < slice_num; i++)
        slice_dualRow[i].update_dual(thetaDual, columnOut);
    }
  }
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;
  shift_back(workHMO, columnOut);//model->shiftBack(columnOut);
}

void HDual::updatePrimal(HVector *DSE_Vector) {
  // Update the primal values and any edge weights
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  // NB DSE_Vector is only computed if dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE
  // Update - primal and weight
  dualRHS.update_primal(&columnBFRT, 1);
  dualRHS.update_infeasList(&columnBFRT);
  double x_out = baseValue[rowOut];
  double l_out = baseLower[rowOut];
  double u_out = baseUpper[rowOut];
  thetaPrimal = (x_out - (deltaPrimal < 0 ? l_out : u_out)) / alpha;
  dualRHS.update_primal(&column, thetaPrimal);
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    double thisEdWt = dualRHS.workEdWt[rowOut] / (alpha * alpha);
    dualRHS.update_weight_DSE(&column, thisEdWt, -2 / alpha,
                              &DSE_Vector->array[0]);
    dualRHS.workEdWt[rowOut] = thisEdWt;
  } else if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    // Pivotal row is for the current basis: weights are required for
    // the next basis so have to divide the current (exact) weight by
    // the pivotal value
    double thisEdWt = dualRHS.workEdWt[rowOut] / (alpha * alpha);
    double dvx_wt_o_rowOut = max(1.0, thisEdWt);
    // nw_wt is max(workEdWt[iRow], NewExactWeight*columnArray[iRow]^2);
    //
    // But NewExactWeight is dvx_wt_o_rowOut = max(1.0, dualRHS.workEdWt[rowOut] / (alpha * alpha))
    //
    // so nw_wt = max(workEdWt[iRow], dvx_wt_o_rowOut*columnArray[iRow]^2);
    //
    // Update rest of weights
    dualRHS.update_weight_Dvx(&column, dvx_wt_o_rowOut);
    dualRHS.workEdWt[rowOut] = dvx_wt_o_rowOut;
    n_dvx_it += 1;
  }
  dualRHS.update_infeasList(&column);

  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    double lc_OpRsDensity = (double)DSE_Vector->count / solver_num_row;
    uOpRsDensityRec(lc_OpRsDensity, rowdseDensity);
  }
  double lc_OpRsDensity = (double)column.count / solver_num_row;
  uOpRsDensityRec(lc_OpRsDensity, columnDensity);

  //  total_fake += column.fakeTick;
  total_syntheticTick += column.syntheticTick;
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    //      total_fake += DSE_Vector->fakeTick;
    total_syntheticTick += DSE_Vector->syntheticTick;
  }
  total_FT_inc_TICK += column.syntheticTick;  // Was .pseudoTick
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    total_FT_inc_TICK += DSE_Vector->syntheticTick;  // Was .pseudoTick
  }
}

void HDual::updatePivots() {
  // UPDATE
  //
  // If reinversion is needed then skip this method
  if (invertHint) return;
  //
  // Update the sets of indices of basic and nonbasic variables
  update_pivots(workHMO, columnIn, rowOut, sourceOut);
  //  checkDualObjectiveValue("After update_pivots");
  //
  // Update the iteration count and store the basis change if HiGHSDEV
  // is defined
  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  workHMO.simplex_info_.iteration_count++;
  //
  // Update the invertible representation of the basis matrix
  update_factor(workHMO, &column, &row_ep, &rowOut, &invertHint);
  //
  // Update the row-wise representation of the nonbasic columns
  update_matrix(workHMO, columnIn, columnOut);
  //
  // Delete Freelist entry for columnIn
  dualRow.delete_Freelist(columnIn);
  //
  // Update the primal value for the row where the basis change has
  // occurred, and set the corresponding squared primal infeasibility
  // value in dualRHS.workArray
  dualRHS.update_pivots(rowOut, workHMO.simplex_info_.workValue_[columnIn] + thetaPrimal);
  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock =
      total_syntheticTick >= factor->build_syntheticTick;
#ifdef HiGHSDEV
  //    bool reinvert_syntheticClock = total_fake >=
  //    factor->build_syntheticTick;
#endif
  if (reinvert_syntheticClock && workHMO.simplex_info_.update_count >= 50) {
    invertHint = INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT;
  }
}

void HDual::iz_dvx_fwk() {
  HighsTimer &timer = workHMO.timer_;
  HighsSimplexInfo &simplex_info = workHMO.simplex_info_;
  // Initialise the Devex framework: reference set is all basic
  // variables
  timer.start(simplex_info.clock_[DevexIzClock]);
  const int *nonbasicFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  for (int vr_n = 0; vr_n < solver_num_tot; vr_n++) {
    //      if (workHMO.simplex_basis_.nonbasicFlag_[vr_n])
    //      if (nonbasicFlag[vr_n])
    //			Nonbasic variables not in reference set
    //	dvx_ix[vr_n] = dvx_not_in_R;
    //      else
    //			Basic variables in reference set
    //	dvx_ix[vr_n] = dvx_in_R;
    //
    // Assume all nonbasic variables have |nonbasicFlag[vr_n]|=1
    //
    // nonbasicFlag[vr_n]*nonbasicFlag[vr_n] evaluates faster than
    // abs(nonbasicFlag[vr_n])
    //
    // int dvx_ix_o_vr = 1-abs(nonbasicFlag[vr_n]);
    //
    int dvx_ix_o_vr = 1 - nonbasicFlag[vr_n] * nonbasicFlag[vr_n];
    dvx_ix[vr_n] = dvx_ix_o_vr;
  }
  dualRHS.workEdWt.assign(solver_num_row, 1.0);  // Set all initial weights to 1
  n_dvx_it = 0;    // Zero the count of iterations with this Devex framework
  n_dvx_fwk += 1;  // Increment the number of Devex frameworks
  nw_dvx_fwk =
      false;  // Indicate that there's no need for a new Devex framework
  timer.stop(simplex_info.clock_[DevexIzClock]);
}

void HDual::interpret_dual_edge_weight_strategy(SimplexDualEdgeWeightStrategy dual_edge_weight_strategy) {
  if (dual_edge_weight_strategy == SimplexDualEdgeWeightStrategy::DANTZIG) {
    dual_edge_weight_mode = DualEdgeWeightMode::DANTZIG;
  } else if (dual_edge_weight_strategy == SimplexDualEdgeWeightStrategy::DEVEX) {
    dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
  } else if (dual_edge_weight_strategy == SimplexDualEdgeWeightStrategy::STEEPEST_EDGE) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = false;
  } else if (dual_edge_weight_strategy == SimplexDualEdgeWeightStrategy::STEEPEST_EDGE_UNIT_INITIAL) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = false;
    allow_dual_steepest_edge_to_devex_switch = false;
  } else if (dual_edge_weight_strategy == SimplexDualEdgeWeightStrategy::STEEPEST_EDGE_TO_DEVEX_SWITCH) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = true;
  } else {
   HighsPrintMessage(ML_MINIMAL, "HDual::interpret_dual_edge_weight_strategy: unrecognised dual_edge_weight_strategy = %d - using dual steepest edge with possible switch to Devex\n", dual_edge_weight_strategy);
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = true;
  }
}

void HDual::interpret_price_strategy(SimplexPriceStrategy price_strategy) {
  allow_price_by_col_switch = false;
  allow_price_by_row_switch = false;
  allow_price_ultra = false;
  if (price_strategy == SimplexPriceStrategy::COL) {
    price_mode = PriceMode::COL;
  } else if (price_strategy == SimplexPriceStrategy::ROW) {
    price_mode = PriceMode::ROW;
  } else if (price_strategy == SimplexPriceStrategy::ROW_SWITCH) {
    price_mode = PriceMode::ROW;
    allow_price_by_row_switch = true;
  } else if (price_strategy == SimplexPriceStrategy::ROW_SWITCH_COL_SWITCH) {
    price_mode = PriceMode::ROW;
    allow_price_by_col_switch = true;
    allow_price_by_row_switch = true;
  } else if (price_strategy == SimplexPriceStrategy::ROW_ULTRA) {
    price_mode = PriceMode::ROW;
    allow_price_by_col_switch = true;
    allow_price_by_row_switch = true;
    allow_price_ultra = true;
  } else {
    HighsPrintMessage(ML_MINIMAL, "HDual::interpret_price_strategy: unrecognised price_strategy = %d - using row Price with switch or colump price switch\n", price_strategy);
    price_mode = PriceMode::ROW;
    allow_price_by_col_switch = true;
    allow_price_by_row_switch = true;
  }
}

#ifdef HiGHSDEV
double HDual::checkDualObjectiveValue(const char *message, int phase) {
  static double previousUpdatedDualObjectiveValue = 0;
  static double previousDualObjectiveValue = 0;
  compute_dual_objective_value(workHMO, phase);
  double updatedDualObjectiveValue = workHMO.simplex_info_.updatedDualObjectiveValue;
  double dualObjectiveValue = workHMO.simplex_info_.dualObjectiveValue;
  double changeInUpdatedDualObjectiveValue = updatedDualObjectiveValue - previousUpdatedDualObjectiveValue;
  double changeInDualObjectiveValue = dualObjectiveValue - previousDualObjectiveValue;
  double updatedDualObjectiveError = dualObjectiveValue - updatedDualObjectiveValue;
  double rlvUpdatedDualObjectiveError = fabs(updatedDualObjectiveError)/max(1.0, fabs(dualObjectiveValue));
  bool erFd = rlvUpdatedDualObjectiveError > 1e-8;
  if (erFd)
    printf("Phase %1d: duObjV = %11.4g (%11.4g); updated duObjV = %11.4g (%11.4g); Error(|Rel|) = %11.4g (%11.4g) |%s\n",
	   phase,
	   dualObjectiveValue, changeInDualObjectiveValue,
	   updatedDualObjectiveValue, changeInUpdatedDualObjectiveValue,
	   updatedDualObjectiveError, rlvUpdatedDualObjectiveError,
	   message);
  previousDualObjectiveValue = dualObjectiveValue;
  previousUpdatedDualObjectiveValue = dualObjectiveValue;
  workHMO.simplex_info_.updatedDualObjectiveValue = dualObjectiveValue;
  // Now have dual objective value
  workHMO.simplex_lp_status_.has_dual_objective_value = true;
  return updatedDualObjectiveError;
}
#endif

// Utility to get a row of the inverse of B for SCIP
int HDual::util_getBasisInvRow(int r, double *coef, int *inds, int *ninds) {
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = r;
  row_ep.array[r] = 1;
  row_ep.packFlag = true;
  factor->btran(row_ep, row_epDensity);
  //  printf("util_getBasisInvRow: nnz = %4d/%4d\n", row_ep.count, solver_num_row);
  for (int row = 0; row < solver_num_row; row++) {
    //    printf("BasisInvRow(%4d) = %11g\n", row,  row_ep.array[row]);
    coef[row] = row_ep.array[row];
  }
  if (0 <= row_ep.count && row_ep.count <= solver_num_row) {
    for (int ix = 0; ix < row_ep.count; ix++) inds[ix] = row_ep.index[ix];
    ninds[0] = row_ep.count;
  } else {
    printf(
        "util_getBasisInvRow: row_ep.count < 0 or row_ep.count > solver_num_row: %4d; "
        "%4d\n",
        row_ep.count, solver_num_row);
    ninds[0] = -1;
  }
  cout << flush;
  return 0;
}

double HDual::an_bs_cond() {
  // Alias to the matrix
  matrix = &workHMO.matrix_;
  const int *Astart = matrix->getAstart();
  const double *Avalue = matrix->getAvalue();
  // Compute the Hager condition number estimate for the basis matrix
  double NoDensity = 1;
  bs_cond_x.resize(solver_num_row);
  bs_cond_y.resize(solver_num_row);
  bs_cond_z.resize(solver_num_row);
  bs_cond_w.resize(solver_num_row);
  // x = ones(n,1)/n;
  // y = A\x;
  double mu = 1.0 / solver_num_row;
  double norm_Binv;
  for (int r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = mu;
  row_ep.clear();
  row_ep.count = solver_num_row;
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    row_ep.index[r_n] = r_n;
    row_ep.array[r_n] = bs_cond_x[r_n];
  }
  for (int ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    factor->ftran(row_ep, NoDensity);
    // zeta = sign(y);
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_y[r_n] = row_ep.array[r_n];
      if (bs_cond_y[r_n] > 0)
        bs_cond_w[r_n] = 1.0;
      else if (bs_cond_y[r_n] < 0)
        bs_cond_w[r_n] = -1.0;
      else
        bs_cond_w[r_n] = 0.0;
    }
    // z=A'\zeta;
    row_ep.clear();
    row_ep.count = solver_num_row;
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      row_ep.index[r_n] = r_n;
      row_ep.array[r_n] = bs_cond_w[r_n];
    }
    row_ep.packFlag = false;
    factor->btran(row_ep, NoDensity);
    // norm_z = norm(z,'inf');
    // ztx = z'*x ;
    // NormEst = norm(y,1);
    // fd_i = 0;
    // for i=1:n
    //    if abs(z(i)) == norm_z
    //        fd_i = i;
    //        break
    //    end
    // end
    double norm_z = 0.0;
    double ztx = 0.0;
    norm_Binv = 0.0;
    int argmax_z = -1;
    for (int r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_z[r_n] = row_ep.array[r_n];
      double abs_z_v = abs(bs_cond_z[r_n]);
      if (abs_z_v > norm_z) {
        norm_z = abs_z_v;
        argmax_z = r_n;
      }
      ztx += bs_cond_z[r_n] * bs_cond_x[r_n];
      norm_Binv += abs(bs_cond_y[r_n]);
    }
    // printf("%2d: ||z||_inf = %8.2g; z^T*x = %8.2g; ||y||_1 = %g\n", ps_n,
    // norm_z, ztx, norm_Binv);
    if (norm_z <= ztx) break;
    // x = zeros(n,1);
    // x(fd_i) = 1;
    for (int r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = 0.0;
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = argmax_z;
    row_ep.array[argmax_z] = 1.0;
    bs_cond_x[argmax_z] = 1.0;
  }
  double norm_B = 0.0;
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    int vr_n = workHMO.simplex_basis_.basicIndex_[r_n];
    double c_norm = 0.0;
    if (vr_n < solver_num_col)
      for (int el_n = Astart[vr_n]; el_n < Astart[vr_n + 1]; el_n++)
        c_norm += abs(Avalue[el_n]);
    else
      c_norm += 1.0;
    norm_B = max(c_norm, norm_B);
  }
  double cond_B = norm_Binv * norm_B;
  printf("Hager estimate of ||B^{-1}||_1 = %g; ||B||_1 = %g so cond_1(B) estimate is %g\n",
	 norm_Binv, norm_B, cond_B);
  return cond_B;
}

#ifdef HiGHSDEV
void HDual::iterateOpRecBf(int opTy, HVector &vector, double hist_dsty) {
  AnIterOpRec *AnIter = &AnIterOp[opTy];
  AnIter->AnIterOpNumCa++;
  double curr_dsty = 1.0 * vector.count / solver_num_row;
  //  printf("%10s: %g<= %g;  %g<= %g\n", AnIter->AnIterOpName.c_str(),
  //	 curr_dsty, AnIter->AnIterOpHyperCANCEL,
  //	 hist_dsty, AnIter->AnIterOpHyperTRAN);
  if (curr_dsty <= AnIter->AnIterOpHyperCANCEL &&
      hist_dsty <= AnIter->AnIterOpHyperTRAN)
    AnIter->AnIterOpNumHyperOp++;
}

void HDual::iterateOpRecAf(int opTy, HVector &vector) {
  AnIterOpRec *AnIter = &AnIterOp[opTy];
  double rsDsty = 1.0 * vector.count / AnIter->AnIterOpRsDim;
  if (rsDsty <= hyperRESULT) AnIter->AnIterOpNumHyperRs++;
  AnIter->AnIterOpRsMxNNZ = max(vector.count, AnIter->AnIterOpRsMxNNZ);
  if (rsDsty > 0) {
    AnIter->AnIterOpLog10RsDsty += log(rsDsty) / log(10.0);
  } else {
    /*
    // TODO Investigate these zero norms
    double vectorNorm = 0;

    for (int index = 0; index < AnIter->AnIterOpRsDim; index++) {
      double vectorValue = vector.array[index];
      vectorNorm += vectorValue * vectorValue;
    }
    vectorNorm = sqrt(vectorNorm);
    printf("Strange: operation %s has result density = %g: ||vector|| = %g\n", AnIter->AnIterOpName.c_str(), rsDsty, vectorNorm);
    */
  }
}

void HDual::iterateRpAn() {
  HighsTimer &timer = workHMO.timer_;
  int AnIterNumIter = workHMO.simplex_info_.iteration_count - AnIterIt0;
  printf("\nAnalysis of %d iterations (%d to %d)\n", AnIterNumIter,
         AnIterIt0 + 1, workHMO.simplex_info_.iteration_count);
  if (AnIterNumIter <= 0) return;
  int lc_EdWtNumIter;
  lc_EdWtNumIter = AnIterNumEdWtIt[(int) DualEdgeWeightMode::STEEPEST_EDGE];
  if (lc_EdWtNumIter > 0)
    printf("DSE for %12d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[(int) DualEdgeWeightMode::DEVEX];
  if (lc_EdWtNumIter > 0)
    printf("Dvx for %12d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  lc_EdWtNumIter = AnIterNumEdWtIt[(int) DualEdgeWeightMode::DANTZIG];
  if (lc_EdWtNumIter > 0)
    printf("Dan for %12d (%3d%%) iterations\n", lc_EdWtNumIter,
           (100 * lc_EdWtNumIter) / AnIterNumIter);
  printf("\n");
  for (int k = 0; k < NumAnIterOpTy; k++) {
    AnIterOpRec *AnIter = &AnIterOp[k];
    int lcNumCa = AnIter->AnIterOpSuNumCa;
    printf("\n%-9s performed %d times\n", AnIter->AnIterOpName.c_str(),
           AnIter->AnIterOpSuNumCa);
    if (lcNumCa > 0) {
      int lcHyperOp = AnIter->AnIterOpSuNumHyperOp;
      int lcHyperRs = AnIter->AnIterOpSuNumHyperRs;
      int pctHyperOp = (100 * lcHyperOp) / lcNumCa;
      int pctHyperRs = (100 * lcHyperRs) / lcNumCa;
      double lcRsDsty = pow(10.0, AnIter->AnIterOpSuLog10RsDsty / lcNumCa);
      int lcAnIterOpRsDim = AnIter->AnIterOpRsDim;
      int lcNumNNz = lcRsDsty * lcAnIterOpRsDim;
      int lcMxNNz = AnIter->AnIterOpRsMxNNZ;
      double lcMxNNzDsty = (1.0 * lcMxNNz) / AnIter->AnIterOpRsDim;
      printf("%12d hyper-sparse operations (%3d%%)\n", lcHyperOp,
             pctHyperOp);
      printf("%12d hyper-sparse results    (%3d%%)\n", lcHyperRs,
             pctHyperRs);
      printf("%12g density of result (%d / %d nonzeros)\n", lcRsDsty, lcNumNNz, lcAnIterOpRsDim);
      printf("%12g density of result with max (%d / %d) nonzeros\n",
             lcMxNNzDsty, lcMxNNz, lcAnIterOpRsDim);
    }
  }
  int NumInvert = 0;

  int last_invert_hint = INVERT_HINT_Count-1;
  for (int k = 1; k <= last_invert_hint; k++) NumInvert += AnIterNumInvert[k];
  if (NumInvert > 0) {
    int lcNumInvert = 0;
    printf("\nInvert    performed %d times: average frequency = %d\n",
           NumInvert, AnIterNumIter / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_UPDATE_LIMIT_REACHED];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to update limit reached\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_SYNTHETIC_CLOCK_SAYS_INVERT];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to pseudo-clock\n", lcNumInvert,
             (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_OPTIMAL];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to possibly optimal\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to possibly primal unbounded\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_DUAL_UNBOUNDED];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to possibly dual unbounded\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_POSSIBLY_SINGULAR_BASIS];
    if (lcNumInvert > 0)
      printf("%12d (%3d%%) Invert operations due to possibly singular basis\n",
             lcNumInvert, (100 * lcNumInvert) / NumInvert);
    lcNumInvert = AnIterNumInvert[INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX];
    if (lcNumInvert > 0)
      printf(
          "%12d (%3d%%) Invert operations due to primal infeasible in primal "
          "simplex\n",
          lcNumInvert, (100 * lcNumInvert) / NumInvert);
  }
  printf("\n%12d (%3d%%) primal degenerate iterations\n", AnIterNumPrDgnIt,
         (100 * AnIterNumPrDgnIt) / AnIterNumIter);
  printf("%12d (%3d%%)   dual degenerate iterations\n", AnIterNumDuDgnIt,
         (100 * AnIterNumDuDgnIt) / AnIterNumIter);
  int suPrice = AnIterNumColPrice + AnIterNumRowPrice + AnIterNumRowPriceWSw +
                AnIterNumRowPriceUltra;
  if (suPrice > 0) {
    printf("\n%12d Price operations:\n", suPrice);
    printf("%12d Col Price      (%3d%%)\n", AnIterNumColPrice,
           (100 * AnIterNumColPrice) / suPrice);
    printf("%12d Row Price      (%3d%%)\n", AnIterNumRowPrice,
           (100 * AnIterNumRowPrice) / suPrice);
    printf("%12d Row PriceWSw   (%3d%%)\n", AnIterNumRowPriceWSw,
           (100 * AnIterNumRowPriceWSw / suPrice));
    printf("%12d Row PriceUltra (%3d%%)\n", AnIterNumRowPriceUltra,
           (100 * AnIterNumRowPriceUltra / suPrice));
  }
  printf("\n%12d (%3d%%) costly DSE        iterations\n", AnIterNumCostlyDseIt,
         (100 * AnIterNumCostlyDseIt) / AnIterNumIter);

  //
  // Add a record for the final iterations: may end up with one more
  // than AN_ITER_TRACE_MX_NUM_REC records, so ensure that there is enough
  // space in the arrays
  //
  AnIterTraceNumRec++;
  AnIterTraceRec *lcAnIter;
  lcAnIter = &AnIterTrace[AnIterTraceNumRec];
  lcAnIter->AnIterTraceIter = workHMO.simplex_info_.iteration_count;
  lcAnIter->AnIterTraceTime = timer.getTime();
  lcAnIter->AnIterTraceDsty[AnIterOpTy_Btran] = row_epDensity;
  lcAnIter->AnIterTraceDsty[AnIterOpTy_Price] = row_apDensity;
  lcAnIter->AnIterTraceDsty[AnIterOpTy_Ftran] = columnDensity;
  lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranBFRT] = columnDensity;
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = rowdseDensity;
    lcAnIter->AnIterTraceAux0 = AnIterCostlyDseMeasure;
  } else {
    lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE] = 0;
    lcAnIter->AnIterTraceAux0 = 0;
  }
  lcAnIter->AnIterTrace_dual_edge_weight_mode = (int) dual_edge_weight_mode;

  if (AnIterTraceIterDl >= 100) {
    printf("\n Iteration speed analysis\n");
    lcAnIter = &AnIterTrace[0];
    int fmIter = lcAnIter->AnIterTraceIter;
    double fmTime = lcAnIter->AnIterTraceTime;
    printf(
        "        Iter (      FmIter:      ToIter)     Time      Iter/sec |  Col R_Ep R_Ap  DSE | "
        "EdWt | Aux0\n");
    for (int rec = 1; rec <= AnIterTraceNumRec; rec++) {
      lcAnIter = &AnIterTrace[rec];
      int toIter = lcAnIter->AnIterTraceIter;
      double toTime = lcAnIter->AnIterTraceTime;
      int dlIter = toIter - fmIter;
      if (rec < AnIterTraceNumRec && dlIter != AnIterTraceIterDl)
        printf("STRANGE: %d = dlIter != AnIterTraceIterDl = %d\n", dlIter,
               AnIterTraceIterDl);
      double dlTime = toTime - fmTime;
      int iterSpeed = 0;
      if (dlTime > 0) iterSpeed = dlIter / dlTime;
      int lc_dual_edge_weight_mode = lcAnIter->AnIterTrace_dual_edge_weight_mode;
      int l10ColDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Ftran]);
      int l10REpDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Btran]);
      int l10RapDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_Price]);
      int l10DseDse = intLog10(lcAnIter->AnIterTraceDsty[AnIterOpTy_FtranDSE]);
      int l10Aux0 = intLog10(lcAnIter->AnIterTraceAux0);
      std::string str_dual_edge_weight_mode;
      if (lc_dual_edge_weight_mode == (int) DualEdgeWeightMode::STEEPEST_EDGE)
        str_dual_edge_weight_mode = "DSE";
      else if (lc_dual_edge_weight_mode == (int) DualEdgeWeightMode::DEVEX)
        str_dual_edge_weight_mode = "Dvx";
      else if (lc_dual_edge_weight_mode == (int) DualEdgeWeightMode::DANTZIG)
        str_dual_edge_weight_mode = "Dan";
      else
        str_dual_edge_weight_mode = "XXX";
      printf("%12d (%12d:%12d) %8.4f  %12d | %4d %4d %4d %4d |  %3s | %4d\n",
             dlIter, fmIter, toIter, dlTime, iterSpeed, l10ColDse, l10REpDse,
             l10RapDse, l10DseDse, str_dual_edge_weight_mode.c_str(), l10Aux0);
      fmIter = toIter;
      fmTime = toTime;
    }
    printf("\n");
  }
}

// TODO Put this in the right place - try to identify when workValue is set
// [nonbasic primals are set to the right bound for dual feasibility?]
void HDual::an_iz_vr_v() {
  double norm_bc_pr_vr = 0;
  double norm_bc_du_vr = 0;
  for (int r_n = 0; r_n < solver_num_row; r_n++) {
    int vr_n = workHMO.simplex_basis_.basicIndex_[r_n];
    norm_bc_pr_vr += baseValue[r_n] * baseValue[r_n];
    norm_bc_du_vr += workDual[vr_n] * workDual[vr_n];
  }
  double norm_nonbc_pr_vr = 0;
  double norm_nonbc_du_vr = 0;
  for (int vr_n = 0; vr_n < solver_num_tot; vr_n++) {
    if (workHMO.simplex_basis_.nonbasicFlag_[vr_n]) {
      double pr_act_v = workHMO.simplex_info_.workValue_[vr_n];
      norm_nonbc_pr_vr += pr_act_v * pr_act_v;
      norm_nonbc_du_vr += workDual[vr_n] * workDual[vr_n];
    }
  }
  // printf("Initial point has %d dual infeasibilities\n", dualInfeasCount);
  // printf("Norm of the basic    primal activites is %g\n",
  // sqrt(norm_bc_pr_vr)); printf("Norm of the basic    dual   activites is
  // %g\n", sqrt(norm_bc_du_vr)); printf("Norm of the nonbasic primal activites
  // is %g\n", sqrt(norm_nonbc_pr_vr)); printf("Norm of the nonbasic dual
  // activites is %g\n", sqrt(norm_nonbc_du_vr));
}
#endif
