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

#include "simplex/HEkkDebug.h"

//#include <cassert>
//#include <cstdio>
//#include <iostream>

//#include "io/HighsIO.h"
//#include "lp_data/HConst.h"
//#include "simplex/HSimplex.h"
//#include "simplex/HSimplexDebug.h"
#include "simplex/SimplexTimer.h"
//#include "util/HighsRandom.h"
//#include "util/HighsUtils.h"

// using std::runtime_error;

HighsStatus HEkkPrimal::solve() {
  HighsOptions& options = ekk_instance_.options_;
  HighsLp& simplex_lp = ekk_instance_.simplex_lp_;
  HighsLogMessage(
      options.logfile, HighsMessageType::INFO,
      "HEkkPrimal::solve called for LP with %d columns, %d rows and %d entries",
      simplex_lp.numCol_, simplex_lp.numRow_,
      simplex_lp.Astart_[simplex_lp.numCol_]);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = ekk_instance_.simplex_lp_.numRow_ > 0;
  if (!positive_num_row) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HEkkPrimal::solve called for LP with non-positive (%d) "
                    "number of constraints",
                    ekk_instance_.simplex_lp_.numRow_);
    assert(positive_num_row);
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }
  if (ekk_instance_.bailoutOnTimeIterations())
    return ekk_instance_.returnFromSolve(HighsStatus::Warning);

  // Set up bound perturbation as cost perturbation in HDual
  if (!use_bound_perturbation)
    HighsLogMessage(options.logfile, HighsMessageType::INFO,
                    "HEkkPrimal::solve not using bound perturbation");

  if (!simplex_lp_status.has_invert) {
    HighsLogMessage(options.logfile, HighsMessageType::ERROR,
                    "HEkkPrimal::solve called without INVERT");
    assert(simplex_lp_status.has_fresh_invert);
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }

  if (use_bound_perturbation) {
    ekk_instance_.computePrimal();
    ekk_instance_.computeSimplexPrimalInfeasible();
  }
  solvePhase =
      ekk_instance_.simplex_info_.num_primal_infeasibilities > 0 ? 1 : 2;

  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         use_bound_perturbation) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);

  //  if (ekkDebugSimplex(ekk_instance_, algorithm, solvePhase) ==
  //  HighsDebugStatus::LOGICAL_ERROR) return
  //  ekk_instance_.returnFromSolve(HighsStatus::Error);

  // The major solving loop
  // Initialise the iteration analysis. Necessary for strategy, but
  // much is for development and only switched on with HiGHSDEV
  // ToDo Move to simplex and adapt so it's OK for primal and dual
  //  iterationAnalysisInitialise();

  //  while (solvePhase) {
  //    int it0 = ekk_instance_.iteration_count_;
  //    switch (solvePhase) {
  //      case 1:
  //        analysis->simplexTimerStart(SimplexPrimalPhase1Clock);
  //        solvePhase1();
  //        analysis->simplexTimerStop(SimplexPrimalPhase1Clock);
  //        simplex_info.primal_phase1_iteration_count +=
  //        (ekk_instance_.iteration_count_ - it0);
  //	break;
  //    case 2:
  //      analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
  //      solvePhase2();
  //      analysis->simplexTimerStop(SimplexPrimalPhase2Clock);
  //      simplex_info.primal_phase2_iteration_count +=
  //      (ekk_instance_.iteration_count_ - it0); break;
  //    case 4:
  //      break;
  //    default:
  //      solvePhase = 0;
  //      break;
  //    }
  //    // Jump for primal
  //    if (solvePhase == 4) break;
  //    // Possibly bail out
  //    if (SolveBailout) break;
  //  }

  solvePhase = 2;
  assert(ekk_instance_.solve_bailout_ == false);
  //  HighsSimplexAnalysis& analysis = ekk_instance_.simplex_analysis_;
  if (solvePhase == 2) {
    int it0 = ekk_instance_.iteration_count_;

    analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
    solvePhase2();
    analysis->simplexTimerStop(SimplexPrimalPhase2Clock);

    simplex_info.primal_phase2_iteration_count +=
        (ekk_instance_.iteration_count_ - it0);
    if (ekk_instance_.bailoutReturn()) return HighsStatus::Warning;
  }
  if (ekkDebugOkForSolve(ekk_instance_, algorithm, solvePhase,
                         use_bound_perturbation) ==
      HighsDebugStatus::LOGICAL_ERROR)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  return HighsStatus::OK;
}

void HEkkPrimal::solvePhase2() {
  //  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
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
  // Possibly bail out immediately if iteration limit is current value
  if (ekk_instance_.bailoutReturn()) return;
  // Setup update limits
  //  simplex_info.update_limit =
  //      min(100 + solver_num_row / 100,
  //          1000);  // TODO: Consider allowing the dual limit to be used
  //  simplex_info.update_count = 0;

  HighsPrintMessage(ekk_instance_.options_.output,
                    ekk_instance_.options_.message_level, ML_DETAILED,
                    "primal-phase2-start\n");
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
          HighsLogMessage(ekk_instance_.options_.logfile,
                          HighsMessageType::ERROR,
                          "Primal phase 1 choose row failed");
          exit(0);
        }

        // Primal phase 1 update
        phase1Update();
        if (invertHint) {
          break;
        }
        if (ekk_instance_.bailoutReturn()) return;
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
      if (ekk_instance_.bailoutReturn()) return;
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
  assert(!ekk_instance_.solve_bailout_);

  if (columnIn == -1) {
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_DETAILED,
                      "primal-optimal\n");
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_DETAILED,
                      "problem-optimal\n");
    ekk_instance_.scaled_model_status_ = HighsModelStatus::OPTIMAL;
  } else {
    HighsPrintMessage(ekk_instance_.options_.output,
                      ekk_instance_.options_.message_level, ML_MINIMAL,
                      "primal-unbounded\n");
    ekk_instance_.scaled_model_status_ = HighsModelStatus::PRIMAL_UNBOUNDED;
  }
  ekk_instance_.computeDualObjectiveValue();
}

void HEkkPrimal::initialise() {
  analysis = &ekk_instance_.analysis_;

  num_col = ekk_instance_.simplex_lp_.numCol_;
  num_row = ekk_instance_.simplex_lp_.numRow_;
  num_tot = num_col + num_row;

  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance =
      ekk_instance_.scaled_solution_params_.primal_feasibility_tolerance;
  dual_feasibility_tolerance =
      ekk_instance_.scaled_solution_params_.dual_feasibility_tolerance;

  ekk_instance_.simplex_lp_status_.has_primal_objective_value = false;
  ekk_instance_.simplex_lp_status_.has_dual_objective_value = false;
  ekk_instance_.scaled_model_status_ = HighsModelStatus::NOTSET;
  ekk_instance_.solve_bailout_ = false;

  // Setup local vectors
  col_aq.setup(num_row);
  row_ep.setup(num_row);
  row_ap.setup(num_col);

  ph1SorterR.reserve(num_row);
  ph1SorterT.reserve(num_row);

  devexReset();

  no_free_columns = true;
  for (int iCol = 0; iCol < num_tot; iCol++) {
    if (ekk_instance_.simplex_info_.workLower_[iCol] == -HIGHS_CONST_INF &&
        ekk_instance_.simplex_info_.workUpper_[iCol] == HIGHS_CONST_INF) {
      // Free column
      no_free_columns = false;
      break;
    }
  }
  if (!no_free_columns)
    HighsLogMessage(ekk_instance_.options_.logfile, HighsMessageType::INFO,
                    "HEkkPrimal:: LP has free columns");
}

void HEkkPrimal::primalRebuild() {
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  isPrimalPhase1 = 0;
  simplex_lp_status.has_fresh_rebuild = true;
  num_flip_since_rebuild = 0;
}

void HEkkPrimal::primalChooseColumn() { columnIn = -1; }

void HEkkPrimal::primalChooseRow() {}

void HEkkPrimal::primalUpdate() {}

void HEkkPrimal::phase1ComputeDual() {}

void HEkkPrimal::phase1ChooseColumn() {}

void HEkkPrimal::phase1ChooseRow() {}

void HEkkPrimal::phase1Update() {}

void HEkkPrimal::devexReset() {
  devex_weight.assign(num_tot, 1.0);
  devex_index.assign(num_tot, 0);
  for (int iVar = 0; iVar < num_tot; iVar++) {
    const int nonbasicFlag = ekk_instance_.simplex_basis_.nonbasicFlag_[iVar];
    devex_index[iVar] = nonbasicFlag * nonbasicFlag;
  }
  num_devex_iterations = 0;
  num_bad_devex_weight = 0;
}

void HEkkPrimal::devexUpdate() {}

void HEkkPrimal::iterationAnalysisData() {}

void HEkkPrimal::iterationAnalysis() {}

void HEkkPrimal::reportRebuild(const int rebuild_invert_hint) {}
