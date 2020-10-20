/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkPrimal.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HEkkPrimal.h"

//#include <cassert>
//#include <cstdio>
//#include <iostream>

//#include "io/HighsIO.h"
//#include "lp_data/HConst.h"
//#include "simplex/HSimplex.h"
//#include "simplex/HSimplexDebug.h"
//#include "simplex/SimplexTimer.h"
//#include "util/HighsRandom.h"
//#include "util/HighsUtils.h"

using std::runtime_error;

HighsStatus HEkkPrimal::solve() {

  HighsOptions& options = ekk_instance.options_;
  HighsLp& simplex_lp = ekk_instance.simplex_lp_;
  HighsLogMessage(options.logfile, HighsMessageType::INFO,
                    "HEkkPrimal::solve called for LP with %d columns, %d rows and %d entries",
      simplex_lp.numCol_, simplex_lp.numRow_,
      simplex_lp.Astart_[simplex_lp.numCol_]);
  /* TEMP
  HighsSimplexInfo& simplex_info = ekk_instance.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;
  ekk_instance.scaled_model_status_ = HighsModelStatus::NOTSET;
  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = ekk_instance.simplex_lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HPrimal::solve called for LP with non-positive (%d) "
                    "number of constraints",
                    ekk_instance.simplex_lp_.numRow_);
    return HighsStatus::Error;
  }
  invertHint = INVERT_HINT_NO;

  // Setup aspects of the model data which are needed for solve() but better
  // left until now for efficiency reasons.
  // ToDo primal simplex version
  // setup_for_solve(ekk_instance);

  // Set SolveBailout to be true if control is to be returned immediately to
  // calling function
  // ToDo Move to Simplex
  //  SolveBailout = false;

  // Initialise working environment
  // Does LOTS, including initialisation of edge weights. Should only
  // be called if model dimension changes
  // ToDo primal simplex version
  // init();
  // initParallel();

  // ToDo primal simplex version
  // initialiseCost(ekk_instance, 1); //  model->initCost(1);
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
  // dual_edge_weight_mode, DualEdgeWeightMode::STEEPEST_EDGE);cout<<flush;
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
    devex_index.assign(solver_num_tot, 0);
    initialiseDevexFramework();
    // Indicate that edge weights are known
    simplex_lp_status.has_dual_steepest_edge_weights = true;
  }
  */

  // ToDo Determine primal simplex phase from initial primal values
  //
  /*
  computePrimal(ekk_instance);
  compute_primal_infeasible_in_??(ekk_instance, &dualInfeasCount);
  solvePhase = ??InfeasCount > 0 ? 1 : 2;
  */

  /* TEMP


  solvePhase = 0;  // Frig to skip while (solvePhase) {*}
  assert(simplex_lp_status.has_primal_objective_value);
  simplex_info.updated_primal_objective_value =
      simplex_info.primal_objective_value;
  solve_bailout = false;
  // Possibly bail out immediately if iteration limit is current value
  if (bailout()) return HighsStatus::Warning;
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
  // ToDo Adapt debugOkForSolve to be used by primal
    if (debugOkForSolve(ekk_instance, solvePhase) == HighsDebugStatus::LOGICAL_ERROR)
    return HighsStatus::Error;
    */
#ifdef HiGHSDEV
    //  reportSimplexLpStatus(simplex_lp_status, "Before HEkkPrimal major solving
    //  loop");
#endif
  // The major solving loop

  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  // ToDo Move to simplex and adapt so it's OK for primal and dual
  //  iterationAnalysisInitialise();

  /* TEMP
  while (solvePhase) {
    /*
    int it0 = ekk_instance.iteration_counts_.simplex;
    switch (solvePhase) {
      case 1:
        analysis->simplexTimerStart(SimplexPrimalPhase1Clock);
        solvePhase1();
        analysis->simplexTimerStop(SimplexPrimalPhase1Clock);
        simplex_info.primal_phase1_iteration_count +=
    (ekk_instance.iteration_counts_.simplex - it0); break; case 2:
        analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
        solvePhase2();
        analysis->simplexTimerStop(SimplexPrimalPhase2Clock);
        simplex_info.primal_phase2_iteration_count +=
    (ekk_instance.iteration_counts_.simplex - it0); break; case 4:
    break; default: solvePhase = 0; break;
    }
    // Jump for primal
    if (solvePhase == 4) break;
    // Possibly bail out
    if (SolveBailout) break;
    */
  /* TEMP
  }
  solvePhase = 2;
  assert(solve_bailout == false);
  analysis = &ekk_instance.simplex_analysis_;
  if (solvePhase == 2) {
    int it0 = ekk_instance.iteration_counts_.simplex;

    analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
    solvePhase2();
    analysis->simplexTimerStop(SimplexPrimalPhase2Clock);

    simplex_info.primal_phase2_iteration_count +=
        (ekk_instance.iteration_counts_.simplex - it0);
    if (bailout()) return HighsStatus::Warning;
  }
  /*
  // ToDo Adapt debugOkForSolve to be used by primal
  if (debugOkForSolve(ekk_instance, solvePhase) == HighsDebugStatus::LOGICAL_ERROR)
    return HighsStatus::Error;
  */
  return HighsStatus::OK;
}

void HEkkPrimal::solvePhase2() {
  /*
  HighsSimplexInfo& simplex_info = ekk_instance.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;
  // When starting a new phase the (updated) primal objective function
  // value isn't known. Indicate this so that when the value
  // computed from scratch in build() isn't checked against the the
  // updated value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Set invertHint so that it's assigned when first tested
  invertHint = INVERT_HINT_NO;
  // Set solvePhase=2 so it's set if solvePhase2() is called directly
  solvePhase = 2;
  solve_bailout = false;
  // Possibly bail out immediately if iteration limit is current value
  if (bailout()) return;
  // Set up local copies of model dimensions
  solver_num_col = ekk_instance.simplex_lp_.numCol_;
  solver_num_row = ekk_instance.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  analysis = &ekk_instance.simplex_analysis_;

  // Setup update limits
  simplex_info.update_limit =
      min(100 + solver_num_row / 100,
          1000);  // TODO: Consider allowing the dual limit to be used
  simplex_info.update_count = 0;

  // Setup local vectors
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);

  ph1SorterR.reserve(solver_num_row);
  ph1SorterT.reserve(solver_num_row);

#ifdef HiGHSDEV
  printf(
      "HEkkPrimal::solvePhase2 - WARNING: Not setting analysis->col_aq_density = "
      "0\n");
  printf(
      "HEkkPrimal::solvePhase2 - WARNING: Not setting analysis->row_ep_density = "
      "0\n");
#endif
  //  analysis->col_aq_density = 0;
  //  analysis->row_ep_density = 0;

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < solver_num_tot; iCol++) {
    if (highs_isInfinity(-ekk_instance.simplex_info_.workLower_[iCol])) {
      if (highs_isInfinity(ekk_instance.simplex_info_.workUpper_[iCol])) {
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

  HighsPrintMessage(ekk_instance.options_.output, ekk_instance.options_.message_level,
                    ML_DETAILED, "primal-phase2-start\n");
  // Main solving structure
  for (;;) {
    analysis->simplexTimerStart(IteratePrimalRebuildClock);
    primalRebuild();
    analysis->simplexTimerStop(IteratePrimalRebuildClock);

    if (isPrimalPhase1) {
      for (;;) {
      // Primal phase 1 choose column
        phase1ChooseColumn();
        if (columnIn == -1) {
          invertHint = INVERT_HINT_CHOOSE_COLUMN_FAIL;
          break;
        }

        // Primal phase 1 choose row
        phase1ChooseRow();
        if (rowOut == -1) {
          HighsLogMessage(ekk_instance.options_.logfile, HighsMessageType::ERROR,
                          "Primal phase 1 choose row failed");
          exit(0);
        }

        // Primal phase 1 update 
        phase1Update();
        if (invertHint) {
          break;
        }
        if (bailout()) return;
      }
      // Go to the next rebuild
      if (invertHint) {
        // Stop when the invert is new
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
      if (bailout()) return;
      if (invertHint) {
        break;
      }
    }
    // If the data are fresh from rebuild() and no flips have occurred, break
    // out of the outer loop to see what's ocurred
    if (simplex_lp_status.has_fresh_rebuild && num_flip_since_rebuild == 0)
      break;
  }
  // If bailing out, should have returned already
  assert(!solve_bailout);

  if (columnIn == -1) {
    HighsPrintMessage(ekk_instance.options_.output, ekk_instance.options_.message_level,
                      ML_DETAILED, "primal-optimal\n");
    HighsPrintMessage(ekk_instance.options_.output, ekk_instance.options_.message_level,
                      ML_DETAILED, "problem-optimal\n");
    ekk_instance.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(ekk_instance.options_.output, ekk_instance.options_.message_level,
                      ML_MINIMAL, "primal-unbounded\n");
    ekk_instance.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  computeDualObjectiveValue(ekk_instance);
  */
}

