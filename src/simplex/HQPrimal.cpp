/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HQPrimal.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HQPrimal.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsRandom.h"
#include "util/HighsUtils.h"

#include <cassert>
#include <cstdio>
#include <iostream>

using std::runtime_error;

HighsStatus HQPrimal::solve() {
  HighsOptions& options = workHMO.options_;
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  workHMO.scaled_model_status_ = HighsModelStatus::NOTSET;
  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = workHMO.simplex_lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
		    "HPrimal::solve called for LP with non-positive (%d) number of constraints",
		    workHMO.simplex_lp_.numRow_);
    return HighsStatus::Error;
  }

  HighsTimer& timer = workHMO.timer_;
  invertHint = INVERT_HINT_NO;

  // Setup aspects of the model data which are needed for solve() but better
  // left until now for efficiency reasons.
  // ToDo primal simplex version
  // setup_for_solve(workHMO);

  // Set SolveBailout to be true if control is to be returned immediately to
  // calling function
  // ToDo Move to Simplex
  //  SolveBailout = false;

  // Initialise working environment
  // Does LOTS, including initialisation of edge weights. Should only
  // be called if model dimension changes
  // ToDo primal simplex version
  // init(num_threads);

  // ToDo primal simplex version
  // initialise_cost(workHMO, 1); //  model->initCost(1);
  assert(simplex_lp_status.has_fresh_invert);
  if (!simplex_lp_status.has_fresh_invert) {
    printf(
        "ERROR: Should enter with fresh INVERT - unless no_invert_on_optimal "
        "is set\n");
  }
  // Consider initialising edge weights - create Primal variants
  //
#ifdef HiGHSDEV
  //  printf("simplex_lp_status.has_dual_steepest_edge_weights 2 = %d;
  //  dual_edge_weight_mode = %d; DualEdgeWeightMode::STEEPEST_EDGE = %d\n",
  //	 simplex_lp_status.has_dual_steepest_edge_weights,
  //dual_edge_weight_mode, DualEdgeWeightMode::STEEPEST_EDGE);cout<<flush;
  //  printf("Edge weights known? %d\n",
  //  !simplex_lp_status.has_dual_steepest_edge_weights);cout<<flush;
#endif
  /*
  if (!simplex_lp_status.has_dual_steepest_edge_weights) {
    // Edge weights are not known
    // Set up edge weights according to dual_edge_weight_mode and
  initialise_dual_steepest_edge_weights
    // Using dual Devex edge weights
    // Zero the number of Devex frameworks used and set up the first one
    n_dvx_fwk = 0;
    dvx_ix.assign(solver_num_tot, 0);
    iz_dvx_fwk();
    // Indicate that edge weights are known
    simplex_lp_status.has_dual_steepest_edge_weights = true;
  }
  */

  // ToDo Determine primal simplex phase from initial primal values
  //
  /*
  compute_primal(workHMO);
  compute_primal_infeasible_in_??(workHMO, &dualInfeasCount);
  solvePhase = ??InfeasCount > 0 ? 1 : 2;
  */
  solvePhase = 0;  // Frig to skip while (solvePhase) {*}

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
  // ToDo Write primal simplex equivalent
  /*
  bool ok = ok_to_solve(workHMO, 1, solvePhase);
  if (!ok) {printf("NOT OK TO SOLVE???\n"); cout << flush;}
  assert(ok);
  */
#ifdef HiGHSDEV
  //  reportSimplexLpStatus(simplex_lp_status, "Before HQPrimal major solving
  //  loop");
#endif
  // The major solving loop

  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  // ToDo Move to simplex and adapt so it's OK for primal and dual
  //  iterationAnalysisInitialise();

  while (solvePhase) {
    /*
    int it0 = scaled_solution_params.simplex_iteration_count;
    switch (solvePhase) {
      case 1:
        timer.start(simplex_info.clock_[SimplexPrimalPhase1Clock]);
        solvePhase1();
        timer.stop(simplex_info.clock_[SimplexPrimalPhase1Clock]);
        simplex_info.primal_phase1_iteration_count +=
    (scaled_solution_params.simplex_iteration_count - it0); break; case 2:
        timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
        solvePhase2();
        timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);
        simplex_info.primal_phase2_iteration_count +=
    (scaled_solution_params.simplex_iteration_count - it0); break; case 4: break; default:
        solvePhase = 0;
        break;
    }
    // Jump for primal
    if (solvePhase == 4) break;
    // Possibly bail out
    if (SolveBailout) break;
    */
  }
  solvePhase = 2;
  if (workHMO.scaled_model_status_ != HighsModelStatus::REACHED_TIME_LIMIT) {
    if (solvePhase == 2) {
      int it0 = scaled_solution_params.simplex_iteration_count;

      timer.start(simplex_info.clock_[SimplexPrimalPhase2Clock]);
      solvePhase2();
      timer.stop(simplex_info.clock_[SimplexPrimalPhase2Clock]);

      simplex_info.primal_phase2_iteration_count +=
          (scaled_solution_params.simplex_iteration_count - it0);
    }
  }
#ifdef HiGHSDEV
  /*
  if (primal_edge_weight_mode == PrimalEdgeWeightMode::DEVEX) {
    printf("Devex: n_dvx_fwk = %d; Average n_dvx_it = %d\n", n_dvx_fwk,
           scaled_solution_params.simplex_iteration_count / n_dvx_fwk);
  }
  */
#endif
  /*
  // ToDo Adapt ok_to_solve to be used by primal
  bool ok = ok_to_solve(workHMO, 1, solvePhase);// model->OKtoSolve(1,
  solvePhase); if (!ok) {printf("NOT OK After Solve???\n"); cout << flush;}
  assert(ok);
  */
  return HighsStatus::OK;
}

void HQPrimal::solvePhase2() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
  printf("HQPrimal::solvePhase2\n");
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = 0;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 2;
  // Set up local copies of model dimensions
  solver_num_col = workHMO.simplex_lp_.numCol_;
  solver_num_row = workHMO.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  analysis = &workHMO.simplex_analysis_;

  // Setup update limits
  simplex_info.update_limit =
      min(100 + solver_num_row / 100,
          1000);  // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);
  columnDensity = 0;
  row_epDensity = 0;

  ph1SorterR.reserve(solver_num_row);
  ph1SorterT.reserve(solver_num_row);

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < solver_num_tot; iCol++) {
    if (highs_isInfinity(-workHMO.simplex_info_.workLower_[iCol])) {
      if (highs_isInfinity(workHMO.simplex_info_.workUpper_[iCol])) {
        // Free column
        no_free_columns = false;
        break;
      }
    }
  }
#ifdef HiGHSDEV
  if (no_free_columns) {
    printf("Model has no free columns\n");
  } else {
    printf("Model has free columns\n");
  }
#endif

  // Setup other buffers

  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, 
		    "primal-phase2-start\n");
  // Main solving structure
  for (;;) {
    timer.start(simplex_info.clock_[IteratePrimalRebuildClock]);
    primalRebuild();
    timer.stop(simplex_info.clock_[IteratePrimalRebuildClock]);

    if (isPrimalPhase1) {
      for (;;) {
        /* Primal phase 1 choose column */
        phase1ChooseColumn();
        if (columnIn == -1) {
          printf("==> Primal phase 1 choose column failed.\n");
          invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
          break;
        }

        /* Primal phsae 1 choose row */
        phase1ChooseRow();
        if (rowOut == -1) {
          printf("Primal phase 1 choose row failed.\n");
          exit(0);
        }

        /* Primal phase 1 update */
        phase1Update();
        if (invertHint) {
          break;
        }
      }
      /* Go to the next rebuild */
      if (invertHint) {
        /* Stop when the invert is new */
        if (simplex_lp_status.has_fresh_rebuild) {
          break;
        }
        continue;
      }
    }

    for (;;) {
      primalChooseColumn();
      if (columnIn == -1) {
        invertHint = INVERT_HINT_POSSIBLY_OPTIMAL;
        break;
      }
      primalChooseRow();
      if (rowOut == -1) {
        invertHint = INVERT_HINT_POSSIBLY_PRIMAL_UNBOUNDED;
        break;
      }
      primalUpdate();
      if (invertHint) {
        break;
      }
    }

    double currentRunHighsTime = timer.readRunHighsClock();
    if (currentRunHighsTime > workHMO.options_.time_limit) {
      workHMO.scaled_model_status_ = HighsModelStatus::REACHED_TIME_LIMIT;
      break;
    }
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) {
#ifdef HiGHSDEV
      if (num_flip_since_rebuild)
        printf("Consider doing a primal rebuild if flips have occurred\n");
#endif
      //      if (num_flip_since_rebuild == 0)
      break;
    }
  }

  if (workHMO.scaled_model_status_ == HighsModelStatus::REACHED_TIME_LIMIT) {
    return;
  }

  if (columnIn == -1) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, 
		      "primal-optimal\n");
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, 
		      "problem-optimal\n");
    workHMO.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL, 
		      "primal-unbounded\n");
    workHMO.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  computeDualObjectiveValue(workHMO);
}

void HQPrimal::primalRebuild() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = workHMO.simplex_lp_status_;
  HighsTimer& timer = workHMO.timer_;
  // Move this to Simplex class once it's created
  //  simplex_method.record_pivots(-1, -1, 0);  // Indicate REINVERT

  // Rebuild workHMO.factor_ - only if we got updates
  int sv_invertHint = invertHint;
  invertHint = INVERT_HINT_NO;
  // Possibly Rebuild workHMO.factor_
  bool reInvert = simplex_info.update_count > 0;
  if (!invert_if_row_out_negative) {
    // Don't reinvert if columnIn is negative [equivalently, if sv_invertHint ==
    // INVERT_HINT_POSSIBLY_OPTIMAL]
    if (sv_invertHint == INVERT_HINT_POSSIBLY_OPTIMAL) {
      assert(columnIn == -1);
      reInvert = false;
    }
  }
  if (reInvert) {
    timer.start(simplex_info.clock_[InvertClock]);
    int rankDeficiency = compute_factor(workHMO);
    timer.stop(simplex_info.clock_[InvertClock]);
    if (rankDeficiency) {
      throw runtime_error("Primal reInvert: singular-basis-matrix");
    }
    simplex_info.update_count = 0;
  }
  timer.start(simplex_info.clock_[ComputeDualClock]);
  compute_dual(workHMO);
  timer.stop(simplex_info.clock_[ComputeDualClock]);

  timer.start(simplex_info.clock_[ComputePrimalClock]);
  compute_primal(workHMO);
  timer.stop(simplex_info.clock_[ComputePrimalClock]);

  // Primal objective section
  timer.start(simplex_info.clock_[ComputePrObjClock]);
  computePrimalObjectiveValue(workHMO);
  timer.stop(simplex_info.clock_[ComputePrObjClock]);

  double primal_objective_value = simplex_info.primal_objective_value;
#ifdef HiGHSDEV
  // Check the objective value maintained by updating against the
  // value when computed exactly - so long as there is a value to
  // check against
  if (simplex_lp_status.has_primal_objective_value) {
    double absPrimalObjectiveError = fabs(
        simplex_info.updated_primal_objective_value - primal_objective_value);
    double rlvPrimalObjectiveError =
        absPrimalObjectiveError / max(1.0, fabs(primal_objective_value));
    // TODO Investigate these Primal objective value errors
    if (rlvPrimalObjectiveError >= 1e-8) {
      HighsLogMessage(workHMO.options_.logfile, HighsMessageType::WARNING,
                      "Primal objective value error |rel| = %12g (%12g)",
                      absPrimalObjectiveError, rlvPrimalObjectiveError);
    }
  }
#endif
  simplex_info.updated_primal_objective_value = primal_objective_value;

  timer.start(simplex_info.clock_[ComputePrIfsClock]);
  computePrimalInfeasible(workHMO);
  timer.stop(simplex_info.clock_[ComputePrIfsClock]);

  timer.start(simplex_info.clock_[ComputeDuIfsClock]);
  computeDualInfeasible(workHMO);
  timer.stop(simplex_info.clock_[ComputeDuIfsClock]);

  /* Whether to switch to primal phase 1 */
  isPrimalPhase1 = 0;
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  if (scaled_solution_params.num_primal_infeasibilities > 0) {
    isPrimalPhase1 = 1;
    phase1ComputeDual();
  }

  timer.start(simplex_info.clock_[ReportRebuildClock]);
  iterationReportRebuild(sv_invertHint);
  timer.stop(simplex_info.clock_[ReportRebuildClock]);
  // Indicate that a header must be printed before the next iteration log
  previous_iteration_report_header_iteration_count = -1;
#ifdef HiGHSDEV
  if (simplex_info.analyseRebuildTime) {
    int iClock = simplex_info.clock_[IteratePrimalRebuildClock];
    int totalRebuilds = timer.clock_num_call[iClock];
    double totalRebuildTime = timer.read(iClock);
    printf(
        "Primal     rebuild %d (%1d) on iteration %9d: Total rebuild time %g\n",
        totalRebuilds, sv_invertHint, scaled_solution_params.simplex_iteration_count,
        totalRebuildTime);
  }
#endif
  num_flip_since_rebuild = 0;
  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HQPrimal::primalChooseColumn() {
  HighsRandom& random = workHMO.random_;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  const int* jFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  const int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double dualTolerance = workHMO.scaled_solution_params_.dual_feasibility_tolerance;

  timer.start(simplex_info.clock_[ChuzcPrimalClock]);
  columnIn = -1;
  double bestInfeas = 0;
  if (no_free_columns) {
    const int numSection = 1;
    int startSection = random.integer() % numSection;
    int deltaCol = (solver_num_tot + numSection - 1) / numSection;
    int fromCol = startSection * deltaCol;
    int toCol = min(fromCol + deltaCol, solver_num_tot);
    int numPass = 1;
    //    printf("\nstartSection = %1d; deltaCol = %d\n", startSection,
    //    deltaCol);
    for (;;) {
      //      printf("CHUZC: %1d [%6d, %6d] %6d\n", numPass, fromCol, toCol,
      //      solver_num_tot);
      for (int iCol = fromCol; iCol < toCol; iCol++) {
        // Then look at dual infeasible
        if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
          if (bestInfeas * devexWeight[iCol] < fabs(workDual[iCol])) {
            bestInfeas = fabs(workDual[iCol]) / devexWeight[iCol];
            columnIn = iCol;
          }
        }
      }
      if (columnIn >= 0 || numPass == numSection) {
        //	printf("Break from CHUZC after %d passes\n", numPass);
        break;
      }
      if (toCol == solver_num_tot) {
        fromCol = 0;
        toCol = deltaCol;
      } else {
        fromCol = toCol;
        toCol = min(fromCol + deltaCol, solver_num_tot);
      }
      numPass++;
    }
  } else {
    for (int iCol = 0; iCol < solver_num_tot; iCol++) {
      if (jFlag[iCol] && fabs(workDual[iCol]) > dualTolerance) {
        // Always take free
        // TODO: if we found free,
        // Then deal with it in dual phase 1
        if (workLower[iCol] == -HIGHS_CONST_INF &&
            workUpper[iCol] == HIGHS_CONST_INF) {
          columnIn = iCol;
          break;
        }
        // Then look at dual infeasible
        if (jMove[iCol] * workDual[iCol] < -dualTolerance) {
          if (bestInfeas * devexWeight[iCol]  < fabs(workDual[iCol])) {
            bestInfeas = fabs(workDual[iCol]) / devexWeight[iCol];
            columnIn = iCol;
          }
        }
      }
    }
  }
  timer.stop(simplex_info.clock_[ChuzcPrimalClock]);
}

void HQPrimal::primalChooseRow() {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  HighsTimer& timer = workHMO.timer_;
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance =
      workHMO.scaled_solution_params_.primal_feasibility_tolerance;

  // Compute pivot column
  timer.start(simplex_info.clock_[FtranClock]);
  col_aq.clear();
  col_aq.packFlag = true;
  workHMO.matrix_.collect_aj(col_aq, columnIn, 1);
  analysis->equalDensity(columnDensity, analysis->col_aq_density);
  workHMO.factor_.ftran(col_aq, analysis->col_aq_density);
  timer.stop(simplex_info.clock_[FtranClock]);

  analysis->equalDensity(columnDensity, analysis->col_aq_density);
  columnDensity = 0.95 * columnDensity + 0.05 * col_aq.count / solver_num_row;
  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density, analysis->col_aq_density);
  analysis->equalDensity(columnDensity, analysis->col_aq_density);

  const bool check_dual = false;
  if (check_dual) {
    const double* workCost = &workHMO.simplex_info_.workCost_[0];
    const double* workDual = &workHMO.simplex_info_.workDual_[0];
    const int* basicIndex = &workHMO.simplex_basis_.basicIndex_[0];
    double check_dual_value = workCost[columnIn];
    for (int i = 0; i < col_aq.count; i++) {
      int row = col_aq.index[i];
      int col = basicIndex[row];
      double value = col_aq.array[row];
      double cost = workCost[col];
      check_dual_value -= value * cost;
      //    printf("Entry %2d: [%2d, %12g] Cost = %12g; check_dual_value =
      //    %12g\n", i, row, value, cost, check_dual_value);
    }
    thetaDual = workDual[columnIn];
    double dual_error =
        fabs(check_dual_value - thetaDual) / max(1.0, fabs(thetaDual));
    if (dual_error > 1e-8)
      printf("Checking dual: updated = %12g; direct = %12g; error = %12g\n",
             thetaDual, check_dual_value, dual_error);
  }

  timer.start(simplex_info.clock_[Chuzr1Clock]);
  // Initialize
  rowOut = -1;

  // Choose row pass 1
  double alphaTol = workHMO.simplex_info_.update_count < 10
                        ? 1e-9
                        : workHMO.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  const int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  int moveIn = jMove[columnIn];
  if (moveIn == 0) {
    // If there's still free in the N
    // We would report not-solved
    // Need to handle free
  }
  double relaxTheta = 1e100;
  double relaxSpace;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    alpha = col_aq.array[index] * moveIn;
    if (alpha > alphaTol) {
      relaxSpace = baseValue[index] - baseLower[index] + primalTolerance;
      if (relaxSpace < relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    } else if (alpha < -alphaTol) {
      relaxSpace = baseValue[index] - baseUpper[index] - primalTolerance;
      if (relaxSpace > relaxTheta * alpha) relaxTheta = relaxSpace / alpha;
    }
  }
  timer.stop(simplex_info.clock_[Chuzr1Clock]);

  timer.start(simplex_info.clock_[Chuzr2Clock]);
  double bestAlpha = 0;
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    alpha = col_aq.array[index] * moveIn;
    if (alpha > alphaTol) {
      // Positive pivotal column entry
      double tightSpace = baseValue[index] - baseLower[index];
      if (tightSpace < relaxTheta * alpha) {
        if (bestAlpha < alpha) {
          bestAlpha = alpha;
          rowOut = index;
        }
      }
    } else if (alpha < -alphaTol) {
      // Negative pivotal column entry
      double tightSpace = baseValue[index] - baseUpper[index];
      if (tightSpace > relaxTheta * alpha) {
        if (bestAlpha < -alpha) {
          bestAlpha = -alpha;
          rowOut = index;
        }
      }
    }
  }
  timer.stop(simplex_info.clock_[Chuzr2Clock]);
}

void HQPrimal::primalUpdate() {
  HighsTimer& timer = workHMO.timer_;
  int* jMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* workValue = &workHMO.simplex_info_.workValue_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const double primalTolerance =
      workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;

  // Compute thetaPrimal
  int moveIn = jMove[columnIn];
  //  int
  columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  //  double
  alpha = col_aq.array[rowOut];
  //  double
  thetaPrimal = 0;
  if (alpha * moveIn > 0) {
    // Lower bound
    thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alpha;
  } else {
    // Upper bound
    thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alpha;
  }

  // 1. Make sure it is inside bounds or just flip bound
  double lowerIn = workLower[columnIn];
  double upperIn = workUpper[columnIn];
  double valueIn = workValue[columnIn] + thetaPrimal;
  bool flipped = false;
  if (jMove[columnIn] == 1) {
    if (valueIn > upperIn + primalTolerance) {
      // Flip to upper
      workValue[columnIn] = upperIn;
      thetaPrimal = upperIn - lowerIn;
      flipped = true;
      jMove[columnIn] = -1;
    }
  } else if (jMove[columnIn] == -1) {
    if (valueIn < lowerIn - primalTolerance) {
      // Flip to lower
      workValue[columnIn] = lowerIn;
      thetaPrimal = lowerIn - upperIn;
      flipped = true;
      jMove[columnIn] = 1;
    }
  }

  timer.start(simplex_info.clock_[UpdatePrimalClock]);
  for (int i = 0; i < col_aq.count; i++) {
    int index = col_aq.index[i];
    baseValue[index] -= thetaPrimal * col_aq.array[index];
  }
  timer.stop(simplex_info.clock_[UpdatePrimalClock]);

  simplex_info.updated_primal_objective_value +=
      workDual[columnIn] * thetaPrimal;

  computePrimalInfeasible(workHMO);

  // If flipped, then no need touch the pivots
  if (flipped) {
    rowOut = -1;
    numericalTrouble = 0;
    thetaDual = workDual[columnIn];
    iterationReport();
    num_flip_since_rebuild++;
    return;
  }

  // Pivot in
  int sourceOut = alpha * moveIn > 0 ? -1 : 1;
  timer.start(simplex_info.clock_[IteratePivotsClock]);
  update_pivots(workHMO, columnIn, rowOut, sourceOut);
  timer.stop(simplex_info.clock_[IteratePivotsClock]);

  baseValue[rowOut] = valueIn;

  timer.start(simplex_info.clock_[CollectPrIfsClock]);
  // Check for any possible infeasible
  for (int iRow = 0; iRow < solver_num_row; iRow++) {
    if (baseValue[iRow] < baseLower[iRow] - primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    } else if (baseValue[iRow] > baseUpper[iRow] + primalTolerance) {
      invertHint = INVERT_HINT_PRIMAL_INFEASIBLE_IN_PRIMAL_SIMPLEX;
    }
  }
  timer.stop(simplex_info.clock_[CollectPrIfsClock]);

  // 2. Now we can update the dual

  timer.start(simplex_info.clock_[BtranClock]);
  row_ep.clear();
  row_ap.clear();
  row_ep.count = 1;
  row_ep.index[0] = rowOut;
  row_ep.array[rowOut] = 1;
  row_ep.packFlag = true;
  analysis->equalDensity(row_epDensity, analysis->row_ep_density);
#ifdef HiGHSDEV
  //  if (simplex_info.analyseSimplexIterations)
  //  iterateOpRecBf(AnIterOpTy_Btran, row_ep, row_epDensity);
#endif
  workHMO.factor_.btran(row_ep, analysis->row_ep_density);
#ifdef HiGHSDEV
  //  if (simplex_info.analyseSimplexIterations)
  //  iterateOpRecAf(AnIterOpTy_Btran, row_ep);
#endif
  timer.stop(simplex_info.clock_[BtranClock]);

  timer.start(simplex_info.clock_[PriceClock]);
  workHMO.matrix_.price_by_row(row_ap, row_ep);
  timer.stop(simplex_info.clock_[PriceClock]);
  analysis->equalDensity(row_epDensity, analysis->row_ep_density);
  const double local_row_ep_density = (double)row_ep.count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_ep_density, analysis->row_ep_density);
  row_epDensity = 0.95 * row_epDensity + 0.05 * row_ep.count / solver_num_row;
  analysis->equalDensity(row_epDensity, analysis->row_ep_density);

  timer.start(simplex_info.clock_[UpdateDualClock]);
  //  double
  thetaDual = workDual[columnIn] / alpha;
  for (int i = 0; i < row_ap.count; i++) {
    int iCol = row_ap.index[i];
    workDual[iCol] -= thetaDual * row_ap.array[iCol];
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iGet = row_ep.index[i];
    int iCol = iGet + solver_num_col;
    workDual[iCol] -= thetaDual * row_ep.array[iGet];
  }
  timer.stop(simplex_info.clock_[UpdateDualClock]);

  /* Update the devex weight */
  devexUpdate();

  // After dual update in primal simplex the dual objective value is not known
  workHMO.simplex_lp_status_.has_dual_objective_value = false;

  // updateVerify for primal
  numericalTrouble = 0;
  /*
  double aCol = fabs(alpha);
  double alphaRow;
  if (columnIn < workHMO.simplex_lp_.numCol_) {
    alphaRow = row_ap.array[columnIn];
  } else {
    alphaRow = row_ep.array[rowOut];
  }
  double aRow = fabs(alphaRow);
  double aDiff = fabs(aCol - aRow);
  numericalTrouble = aDiff / min(aCol, aRow);
  if (numericalTrouble > 1e-7)
    printf("Numerical check: alphaCol = %12g, alphaRow = a%12g, aDiff = a%12g:
  measure = %12g\n", alpha, alphaRow, aDiff, numericalTrouble);
  // Reinvert if the relative difference is large enough, and updates have been
  performed
  //  if (numericalTrouble > 1e-7 && workHMO.simplex_info_.update_count > 0)
  invertHint = INVERT_HINT_POSSIBLY_SINGULAR_BASIS;
  */
  // Dual for the pivot
  workDual[columnIn] = 0;
  workDual[columnOut] = -thetaDual;

  // Update workHMO.factor_ basis
  update_factor(workHMO, &col_aq, &row_ep, &rowOut, &invertHint);
  update_matrix(workHMO, columnIn, columnOut);
  if (simplex_info.update_count >= simplex_info.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }
  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  workHMO.scaled_solution_params_.simplex_iteration_count++;

  /* Reset the devex when there are too many errors */
  if(nBadDevexWeight > 3) {
    devexReset();
  }

  // Report on the iteration
  iterationReport();
}

/* Compute the reduced cost for primal phase 1 with artificial cost. */
void HQPrimal::phase1ComputeDual() {
  /* Alias to problem size, tolerance and work arrays */
  const int nRow = workHMO.lp_.numRow_;
  const int nCol = workHMO.lp_.numCol_;
  const double dFeasTol = workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  const double *baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double *baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  const double *baseValue = &workHMO.simplex_info_.baseValue_[0];

  /* Setup artificial cost and compute pi with BTran */
  HVector buffer;
  buffer.setup(nRow);
  buffer.clear();
  for (int iRow = 0; iRow < nRow; iRow++) {
    buffer.index[iRow] = iRow;
    if (baseValue[iRow] <  baseLower[iRow] - dFeasTol) {
      buffer.array[iRow] = -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
      buffer.array[iRow] = 1.0;
    } else {
      buffer.array[iRow] = 0.0;
    }
  }
  workHMO.factor_.btran(buffer, 1);

  /* The compute the phase 1 reduced cost for all variables by SpMV */
  HVector bufferLong;
  bufferLong.setup(nCol);
  bufferLong.clear();
  workHMO.matrix_.price_by_col(bufferLong, buffer);
  const int* nbFlag = &workHMO.simplex_basis_.nonbasicFlag_[0];
  double *workDual = &workHMO.simplex_info_.workDual_[0];
  for (int iSeq = 0; iSeq < nCol + nRow; iSeq++) {
    workDual[iSeq] = 0.0;
  }
  for (int iSeq = 0; iSeq < nCol; iSeq++) {
    if (nbFlag[iSeq])
      workDual[iSeq] = -bufferLong.array[iSeq];
  }
  for (int iRow = 0, iSeq = nCol; iRow < nRow; iRow++, iSeq++) {
    if (nbFlag[iSeq])
      workDual[iSeq] = -buffer.array[iRow];
  }

  /* Recompute number of dual infeasible variables with the phase 1 cost */
  computeDualInfeasible(workHMO);
}

/* Choose a pivot column for the phase 1 primal simplex method */
void HQPrimal::phase1ChooseColumn() {
  const int nSeq = workHMO.lp_.numCol_ + workHMO.lp_.numRow_;
  const int* nbMove = &workHMO.simplex_basis_.nonbasicMove_[0];
  const double* workDual = &workHMO.simplex_info_.workDual_[0];
  const double dDualTol = workHMO.scaled_solution_params_.dual_feasibility_tolerance;
  double dBestScore = 0;
  columnIn = -1;
  for (int iSeq = 0; iSeq < nSeq; iSeq++) {
    double dMyDual = nbMove[iSeq] * workDual[iSeq];
    double dMyScore = dMyDual / devexWeight[iSeq];
    if (dMyDual < -dDualTol && dMyScore < dBestScore) {
      dBestScore = dMyScore;
      columnIn = iSeq;
    }
  }
}

/* Choose a pivot row for the phase 1 primal simplex method */
void HQPrimal::phase1ChooseRow() {
  /* Alias to work arrays */
  const double dFeasTol = workHMO.scaled_solution_params_.primal_feasibility_tolerance;
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  const double* baseValue = &workHMO.simplex_info_.baseValue_[0];

  /* Compute the transformed pivot column and update its density */
  col_aq.clear();
  col_aq.packFlag = true;
  workHMO.matrix_.collect_aj(col_aq, columnIn, 1);
  analysis->equalDensity(columnDensity, analysis->col_aq_density);
  workHMO.factor_.ftran(col_aq, analysis->col_aq_density);

  analysis->equalDensity(columnDensity, analysis->col_aq_density);
  columnDensity = 0.95 * columnDensity + 0.05 * col_aq.count / solver_num_row;
  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density, analysis->col_aq_density);
  analysis->equalDensity(columnDensity, analysis->col_aq_density);

  /* Compute the reducedc cost for the pivot column and compare it with the kept value */
  double dCompDual = 0.0;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    if (baseValue[iRow] < baseLower[iRow] - dFeasTol) {
      dCompDual -= col_aq.array[iRow] * -1.0;
    } else if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
      dCompDual -= col_aq.array[iRow] * +1.0;
    }
  }
  if (fabs(workHMO.simplex_info_.workDual_[columnIn] - dCompDual) > (fabs(dCompDual) + 1.0) * 1e-9) {
    printf("==> Phase 1 reduced cost. Updated %g, Computed %g\n", workHMO.simplex_info_.workDual_[columnIn], dCompDual);
  }

  /* Collect phase 1 theta lists */
  int nRow = workHMO.lp_.numRow_;
  const int iMoveIn = workHMO.simplex_basis_.nonbasicMove_[columnIn];
  const double dPivotTol = workHMO.simplex_info_.update_count < 10 ? 1e-9 :
                           workHMO.simplex_info_.update_count < 20 ? 1e-8 : 1e-7;
  ph1SorterR.clear();
  ph1SorterT.clear();
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    double dAlpha = col_aq.array[iRow] * iMoveIn;

    /* When the basic variable x[i] decrease */
    if (dAlpha > +dPivotTol) {
      /* Whether it can become feasible by going below its upper bound */
      if (baseValue[iRow] > baseUpper[iRow] + dFeasTol) {
        double dFeasTheta = (baseValue[iRow] - baseUpper[iRow] - dFeasTol) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow));
      }
      /* Whether it can become infeasible (again) by going below its lower bound */
      if (baseValue[iRow] > baseLower[iRow] - dFeasTol && baseLower[iRow] > -HIGHS_CONST_INF) {
        double dRelaxTheta = (baseValue[iRow] - baseLower[iRow] + dFeasTol) / dAlpha;
        double dTightTheta = (baseValue[iRow] - baseLower[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow - nRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow - nRow));
      }
    }

    /* When the basic variable x[i] increase */
    if (dAlpha < -dPivotTol) {
      /* Whether it can become feasible by going above its lower bound */
      if (baseValue[iRow] < baseLower[iRow] - dFeasTol) {
        double dFeasTheta = (baseValue[iRow] - baseLower[iRow] + dFeasTol) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dFeasTheta, iRow - nRow));
        ph1SorterT.push_back(std::make_pair(dFeasTheta, iRow - nRow));
      }

      /* Whether it can become infeasible (again) by going above its upper bound */
      if (baseValue[iRow] < baseUpper[iRow] + dFeasTol && baseUpper[iRow] < +HIGHS_CONST_INF) {
        double dRelaxTheta = (baseValue[iRow] - baseUpper[iRow] - dFeasTol) / dAlpha;
        double dTightTheta = (baseValue[iRow] - baseUpper[iRow]) / dAlpha;
        ph1SorterR.push_back(std::make_pair(dRelaxTheta, iRow));
        ph1SorterT.push_back(std::make_pair(dTightTheta, iRow));
      }
    }
  }

  /* When there is no candidates at all, we can leave it here */
  if (ph1SorterR.empty()) {
    rowOut = -1;
    columnOut = -1;
    return;
  }

  /*
   * Now sort the relaxed theta to find the final break point.
   * TODO: Consider partial sort.
   *       Or heapify [O(n)] and then pop k points [kO(log(n))].
   */
  std::sort(ph1SorterR.begin(), ph1SorterR.end());
  double dMaxTheta = ph1SorterR.at(0).first;
  double dGradient = fabs(workHMO.simplex_info_.workDual_[columnIn]);
  for (unsigned int i = 0; i < ph1SorterR.size(); i++) {
    double dMyTheta = ph1SorterR.at(i).first;
    int index = ph1SorterR.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    dGradient -= fabs(col_aq.array[iRow]);
    /* Stop when the gradient start to decrease */
    if (dGradient <= 0) {
      break;
    }
    dMaxTheta = dMyTheta;
  }

  /* Find out the biggest possible alpha for pivot */
  std::sort(ph1SorterT.begin(), ph1SorterT.end());
  double dMaxAlpha = 0.0;
  unsigned int iLast = ph1SorterT.size();
  for (unsigned int i = 0; i < ph1SorterT.size(); i++) {
    double dMyTheta = ph1SorterT.at(i).first;
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    /* Stop when the theta is too large */
    if (dMyTheta > dMaxTheta) {
      iLast = i;
      break;
    }
    /* Update the maximal possible alpha */
    if (dMaxAlpha < dAbsAlpha) {
      dMaxAlpha = dAbsAlpha;
    }
  }

  /* Finally choose a pivot with good enough alpha, backwardly */
  rowOut = -1;
  columnOut = -1;
  phase1OutBnd = 0;
  for (int i = iLast - 1; i >= 0; i--) {
    int index = ph1SorterT.at(i).second;
    int iRow = index >= 0 ? index : index + nRow;
    double dAbsAlpha = fabs(col_aq.array[iRow]);
    if (dAbsAlpha > dMaxAlpha * 0.1) {
      rowOut = iRow;
      phase1OutBnd = index > 0 ? 1 : -1;
      break;
    }
  }
  if(rowOut != -1) {
    columnOut = workHMO.simplex_basis_.basicIndex_[rowOut];
  }
}

/* Update the primal and dual solutions for the primal phase 1 */
void HQPrimal::phase1Update() {
  /* Alias to bounds arrays */
  const double* workLower = &workHMO.simplex_info_.workLower_[0];
  const double* workUpper = &workHMO.simplex_info_.workUpper_[0];
  const double* baseLower = &workHMO.simplex_info_.baseLower_[0];
  const double* baseUpper = &workHMO.simplex_info_.baseUpper_[0];
  double* workValue = &workHMO.simplex_info_.workValue_[0];
  double* baseValue = &workHMO.simplex_info_.baseValue_[0];
  const int iMoveIn = workHMO.simplex_basis_.nonbasicMove_[columnIn];
  const double dFeasTol = workHMO.scaled_solution_params_.primal_feasibility_tolerance;

  /* Compute the primal theta and see if we should have do bound flip instead */
  alpha = col_aq.array[rowOut];
  thetaPrimal = 0.0;
  if(phase1OutBnd == 1) {
    thetaPrimal = (baseValue[rowOut] - baseUpper[rowOut]) / alpha;
  } else {
    thetaPrimal = (baseValue[rowOut] - baseLower[rowOut]) / alpha;
  }

  double lowerIn = workLower[columnIn];
  double upperIn = workUpper[columnIn];
  double valueIn = workValue[columnIn] + thetaPrimal;
  int ifFlip = 0;
  if (iMoveIn == +1 && valueIn > upperIn + dFeasTol) {
    workValue[columnIn] = upperIn;
    thetaPrimal = upperIn - lowerIn;
    ifFlip = 1;
    workHMO.simplex_basis_.nonbasicMove_[columnIn] = -1;
  }
  if (iMoveIn == -1 && valueIn < lowerIn - dFeasTol) {
    workValue[columnIn] = lowerIn;
    thetaPrimal = lowerIn - upperIn;
    ifFlip = 1;
    workHMO.simplex_basis_.nonbasicMove_[columnIn] = +1;
  }

  /* Update for the flip case */
  if(ifFlip) {
    /* Recompute things on flip */
    if (invertHint == 0) {
      compute_primal(workHMO);
      computePrimalInfeasible(workHMO);
      if (workHMO.scaled_solution_params_.num_primal_infeasibilities > 0) {
        isPrimalPhase1 = 1;
        phase1ComputeDual();
      } else {
        invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
      }
    }
    return;
  }

  /* Compute BTran for update LU */
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = rowOut;
  row_ep.array[rowOut] = 1;
  row_ep.packFlag = true;
  analysis->equalDensity(row_epDensity, analysis->row_ep_density);
  workHMO.factor_.btran(row_ep, analysis->row_ep_density);

  analysis->equalDensity(row_epDensity, analysis->row_ep_density);
  const double local_row_ep_density = (double)row_ep.count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_ep_density, analysis->row_ep_density);
  row_epDensity = 0.95 * row_epDensity + 0.05 * row_ep.count / solver_num_row;
  analysis->equalDensity(row_epDensity, analysis->row_ep_density);

  /* Compute the whole pivot row for updating the devex weight */
  row_ap.clear();
  workHMO.matrix_.price_by_row(row_ap, row_ep);


  /* Update the devex weight */
  devexUpdate();

   /* Update other things */
  update_pivots(workHMO, columnIn, rowOut, phase1OutBnd);
  update_factor(workHMO, &col_aq, &row_ep, &rowOut, &invertHint);
  update_matrix(workHMO, columnIn, columnOut);
  if (workHMO.simplex_info_.update_count >= workHMO.simplex_info_.update_limit) {
    invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
  }


  /* Recompute dual and primal */
  if (invertHint == 0) {
    compute_primal(workHMO);
    computePrimalInfeasible(workHMO);
    if (workHMO.scaled_solution_params_.num_primal_infeasibilities > 0) {
      isPrimalPhase1 = 1;
      phase1ComputeDual();
    } else {
      invertHint = INVERT_HINT_UPDATE_LIMIT_REACHED;
    }
  }

  /* Reset the devex framework when necessary */
  if(nBadDevexWeight > 3) {
    devexReset();
  }


  // Move this to Simplex class once it's created
  // simplex_method.record_pivots(columnIn, columnOut, alpha);
  workHMO.scaled_solution_params_.simplex_iteration_count++;
}

/* Reset the devex weight */
void HQPrimal::devexReset() {
  const int nSeq = workHMO.lp_.numCol_ + workHMO.lp_.numRow_;
  devexWeight.assign(nSeq, 1.0);
  devexRefSet.assign(nSeq, 0);
  for (int iSeq = 0; iSeq < nSeq; iSeq++) {
    if (workHMO.simplex_basis_.nonbasicFlag_[iSeq]) {
      devexRefSet[iSeq] = 1;
    }
  }
  nBadDevexWeight = 0;
}

void HQPrimal::devexUpdate() {
  /* Compute the pivot weight from the reference set */
  double dPivotWeight = 0.0;
  for (int i = 0; i < col_aq.count; i++) {
    int iRow = col_aq.index[i];
    int iSeq = workHMO.simplex_basis_.basicIndex_[iRow];
    if (devexRefSet[iSeq]) {
      double dAlpha = col_aq.array[iRow];
      dPivotWeight += dAlpha * dAlpha;
    }
  }
  if (devexRefSet[columnIn]) {
    dPivotWeight += 1.0;
  }
  dPivotWeight = sqrt(dPivotWeight);

  /* Check if the saved weight is too large */
  if (devexWeight[columnIn] > 3.0 * dPivotWeight) {
    nBadDevexWeight++;
  }

  /* Update the devex weight for all */
  double dPivot = col_aq.array[rowOut];
  dPivotWeight /= fabs(dPivot);

  for (int i = 0; i < row_ap.count; i++) {
    int iSeq = row_ap.index[i];
    double alpha = row_ap.array[iSeq];
    double devex = dPivotWeight * fabs(alpha);
    if (devexRefSet[iSeq]) {
      devex += 1.0;
    }
    if (devexWeight[iSeq] < devex) {
      devexWeight[iSeq] = devex;
    }
  }
  for (int i = 0; i < row_ep.count; i++) {
    int iPtr = row_ep.index[i];
    int iSeq = row_ep.index[i] + solver_num_col;
    double alpha = row_ep.array[iPtr];
    double devex = dPivotWeight * fabs(alpha);
    if (devexRefSet[iSeq]) {
      devex += 1.0;
    }
    if (devexWeight[iSeq] < devex) {
      devexWeight[iSeq] = devex;
    }
  }

  /* Update devex weight for the pivots */
  devexWeight[columnOut] = max(1.0, dPivotWeight);
  devexWeight[columnIn] = 1.0;
}

void HQPrimal::iterationReport() {
  int iteration_count = workHMO.scaled_solution_params_.simplex_iteration_count;
  int iteration_count_difference = iteration_count -
    previous_iteration_report_header_iteration_count;
  bool header = (previous_iteration_report_header_iteration_count < 0)
    || (iteration_count_difference > 10);
  if (header) {
    iterationReportFull(header);
    previous_iteration_report_header_iteration_count = iteration_count;
  }
  iterationReportFull(false);
}

void HQPrimal::iterationReportFull(bool header) {
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  bool report_condition = simplex_info.analyse_invert_condition;
#endif
  if (header) {
    iterationReportIterationAndPhase(ML_DETAILED, true);
    iterationReportPrimalObjective(ML_DETAILED, true);
#ifdef HiGHSDEV
    iterationReportIterationData(ML_DETAILED, true);
    //    iterationReportDsty(ML_DETAILED, true);
    //    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, " FreeLsZ");
    if (report_condition)
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, 
			"   Condition");
#endif
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, "\n");
  } else {
    iterationReportIterationAndPhase(ML_DETAILED, false);
    iterationReportPrimalObjective(ML_DETAILED, false);
#ifdef HiGHSDEV
    iterationReportIterationData(ML_DETAILED, false);
    //    iterationReportDsty(ML_DETAILED, false);
    //    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, " %7d", dualRow.freeListSize);
    if (report_condition)
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, 
			" %11.4g", simplex_info.invert_condition);
#endif
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED, "\n");
  }
}

void HQPrimal::iterationReportIterationAndPhase(int iterate_log_level,
                                               bool header) {
  int solvePhase = 2;
  if (header) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
		      " Iteration Ph");
  } else {
    int iteration_count = workHMO.scaled_solution_params_.simplex_iteration_count;
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
		      " %9d %2d", iteration_count, solvePhase);
  }
}

void HQPrimal::iterationReportPrimalObjective(int iterate_log_level,
                                             bool header) {
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  if (header) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
		      "      PrimalObjective");
  } else {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
		      " %20.10e", simplex_info.updated_primal_objective_value);
  }
}

void HQPrimal::iterationReportIterationData(int iterate_log_level, bool header) {
  if (header) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level,
                      " Inv       NumCk     EnC     LvR     LvC        ThDu    "
                      "    ThPr          Aa");
  } else {
    bool flipped = rowOut < 0;
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
		      " %3d %11.4g %7d", invertHint, numericalTrouble, columnIn);
    if (flipped) {
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level,
                        "                 %11.4g %11.4g            ", thetaDual, thetaPrimal);
    } else {
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
			" %7d %7d %11.4g %11.4g %11.4g", rowOut, columnOut, thetaDual, thetaPrimal, alpha);
    }
  }
}

/*
void HQPrimal::iterationReportDsty(int iterate_log_level, bool header) {
  bool rp_dual_steepest_edge = dual_edge_weight_mode ==
DualEdgeWeightMode::STEEPEST_EDGE; if (header) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
    "  Col R_Ep R_Ap");
    if (rp_dual_steepest_edge) {
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
      "  DSE");
    } else {
      HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
      "     ");
    }
  } else {
    int l10ColDse = intLog10(analysis->col_aq_density);
    int l10REpDse = intLog10(analysis->row_ep_density);
    int l10RapDse = intLog10(analysis->row_ap_density);
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
    " %4d %4d %4d", l10ColDse, l10REpDse,
l10RapDse); if (rp_dual_steepest_edge) { int l10DseDse =
intLog10(analysis->row_DSE_density); HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
" %4d",
l10DseDse); } else { HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, iterate_log_level, 
"     ");
    }
  }
}
int HQPrimal::intLog10(double v) {
  int intLog10V = -99;
  if (v > 0) intLog10V = log(v) / log(10.0);
  return intLog10V;
}

*/
void HQPrimal::iterationReportRebuild(const int i_v) {
#ifdef HiGHSDEV
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  bool report_condition = simplex_info.analyse_invert_condition;
  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL,
                    "Iter %10d:", workHMO.scaled_solution_params_.simplex_iteration_count);
  //  iterationReportDsty(ML_MINIMAL, true);
  //  iterationReportDsty(ML_MINIMAL, false);
  iterationReportPrimalObjective(ML_MINIMAL, false);
  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL, 
		    " PrPh%1d(%2d)", solvePhase, i_v);
  if (report_condition)
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL, 
		      " k(B)%10.4g", simplex_info.invert_condition);
  if (solvePhase == 2) reportInfeasibility();
  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL, "\n");
#else
  logRebuild(workHMO, true, solvePhase);
#endif
}

void HQPrimal::reportInfeasibility() {
  HighsSolutionParams& scaled_solution_params = workHMO.scaled_solution_params_;
  if (scaled_solution_params.sum_primal_infeasibilities > 0) {
    HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL, 
		      " Pr: %d(%g);",
                      scaled_solution_params.num_primal_infeasibilities,
                      scaled_solution_params.sum_primal_infeasibilities);
  }
  HighsPrintMessage(workHMO.options_.output, workHMO.options_.message_level, ML_MINIMAL, 
		    " Du: %d(%g)",
                    scaled_solution_params.num_dual_infeasibilities,
                    scaled_solution_params.sum_dual_infeasibilities);
}
