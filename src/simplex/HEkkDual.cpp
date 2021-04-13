/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkDual.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "HEkkDual.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HCrash.h"
#include "simplex/HEkkPrimal.h"
//#include "simplex/HSimplex.h"
#include "simplex/HEkkDebug.h"
#include "simplex/HSimplexReport.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsTimer.h"

#ifdef OPENMP
#include "omp.h"
#endif

using std::cout;
using std::endl;
using std::fabs;
using std::flush;
using std::runtime_error;

HighsStatus HEkkDual::solve() {
  assert(SOLVE_PHASE_ERROR == -3);
  assert(SOLVE_PHASE_EXIT == -2);
  assert(SOLVE_PHASE_UNKNOWN == -1);
  assert(SOLVE_PHASE_OPTIMAL == 0);
  assert(SOLVE_PHASE_1 == 1);
  assert(SOLVE_PHASE_2 == 2);
  assert(SOLVE_PHASE_CLEANUP == 4);
  HighsOptions& options = ekk_instance_.options_;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  HighsModelStatus& scaled_model_status = ekk_instance_.scaled_model_status_;
  scaled_model_status = HighsModelStatus::NOTSET;
  if (debugDualSimplex("Initialise", true) == HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  // Assumes that the LP has a positive number of rows
  bool positive_num_row = ekk_instance_.simplex_lp_.numRow_ > 0;
  if (!positive_num_row) {
    highsLogUser(
        options.log_options, HighsLogType::ERROR,
        "HEkkDual::solve called for LP with non-positive (%" HIGHSINT_FORMAT
        ") "
        "number of constraints\n",
        ekk_instance_.simplex_lp_.numRow_);
    assert(positive_num_row);
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }
  rebuild_reason = REBUILD_REASON_NO;

  // Set solve_bailout to be true if control is to be returned immediately to
  // calling function
  solve_bailout = false;
  if (bailoutOnTimeIterations())
    return ekk_instance_.returnFromSolve(HighsStatus::Warning);

  // Initialise working environment.
  init();
  initParallel();

  bool dual_info_ok = dualInfoOk(ekk_instance_.simplex_lp_);
  if (!dual_info_ok) {
    highsLogUser(options.log_options, HighsLogType::ERROR,
                 "HPrimalDual::solve has error in dual information\n");
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }

  // Decide whether to use LiDSE by not storing squared primal infeasibilities
  simplex_info.store_squared_primal_infeasibility = true;
  if (options.less_infeasible_DSE_check) {
    if (isLessInfeasibleDSECandidate(options.log_options,
                                     ekk_instance_.simplex_lp_)) {
      // LP is a candidate for LiDSE
      if (options.less_infeasible_DSE_choose_row)
        // Use LiDSE
        simplex_info.store_squared_primal_infeasibility = false;
    }
  }

  // Determine whether the solution is near-optimal. Value 1 is
  // unimportant, as the sum of primal infeasiblilities for
  // near-optimal solutions is typically many orders of magnitude
  // smaller than 1, and the sum of primal infeasiblilities will be
  // very much larger for non-trivial LPs that are dual feasible for a
  // logical or crash basis.
  const bool near_optimal = simplex_info.num_dual_infeasibility == 0 &&
                            simplex_info.sum_primal_infeasibility < 1;
  if (near_optimal)
    highsLogDev(options.log_options, HighsLogType::DETAILED,
                "Dual feasible and num / max / sum primal infeasibilities are "
                "%" HIGHSINT_FORMAT
                " / %g "
                "/ %g, so near-optimal\n",
                simplex_info.num_primal_infeasibility,
                simplex_info.max_primal_infeasibility,
                simplex_info.sum_primal_infeasibility);

  // Perturb costs according to whether the solution is near-optimnal
  const bool perturb_costs = !near_optimal;
  if (!perturb_costs)
    highsLogDev(options.log_options, HighsLogType::DETAILED,
                "Near-optimal, so don't use cost perturbation\n");
  ekk_instance_.initialiseCost(SimplexAlgorithm::DUAL, SOLVE_PHASE_UNKNOWN,
                               perturb_costs);
  assert(simplex_lp_status.has_invert);
  if (!simplex_lp_status.has_invert) {
    highsLogUser(options.log_options, HighsLogType::ERROR,
                 "HPrimalDual:: Should enter solve with INVERT\n");
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }
  // Consider initialising edge weights
  //
  // NB workEdWt is assigned and initialised to 1s in
  // dualRHS.setup(ekk_instance_) so that CHUZR is well defined, even for
  // Dantzig pricing
  //
  if (!simplex_lp_status.has_dual_steepest_edge_weights) {
    // Edge weights are not known
    // Set up edge weights according to dual_edge_weight_mode and
    // initialise_dual_steepest_edge_weights

    if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
      // Using dual Devex edge weights, so set up the first framework
      simplex_info.devex_index_.assign(solver_num_tot, 0);
      initialiseDevexFramework();
    } else if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // Intending to using dual steepest edge (DSE) weights
      if (initialise_dual_steepest_edge_weights) {
        // Exact DSE weights need to be computed if the basis contains
        // structurals
        bool logical_basis = true;
        for (HighsInt iRow = 0; iRow < solver_num_row; iRow++) {
          if (ekk_instance_.simplex_basis_.basicIndex_[iRow] < solver_num_col) {
            logical_basis = false;
            break;
          }
        }
        if (!logical_basis) {
          if (near_optimal) {
            // Basis is not logical but near optimal, so use Devex
            // rather than initialise DSE weights
            highsLogDev(
                options.log_options, HighsLogType::DETAILED,
                "Basis is not logical, but near-optimal so use Devex rather "
                "than compute exact weights for DSE\n");
            dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
            simplex_info.devex_index_.assign(solver_num_tot, 0);
            initialiseDevexFramework();
          } else {
            // Basis is not logical and DSE weights are to be initialised
            highsLogDev(options.log_options, HighsLogType::DETAILED,
                        "Basis is not logical, so compute exact DSE weights\n");
            if (ekk_instance_.analysis_.analyse_simplex_time) {
              analysis->simplexTimerStart(SimplexIzDseWtClock);
              analysis->simplexTimerStart(DseIzClock);
            }
            for (HighsInt i = 0; i < solver_num_row; i++) {
              row_ep.clear();
              row_ep.count = 1;
              row_ep.index[0] = i;
              row_ep.array[i] = 1;
              row_ep.packFlag = false;
              factor->btran(row_ep, analysis->row_ep_density,
                            analysis->pointer_serial_factor_clocks);
              dualRHS.workEdWt[i] = row_ep.norm2();
              const double local_row_ep_density =
                  (double)row_ep.count / solver_num_row;
              analysis->updateOperationResultDensity(local_row_ep_density,
                                                     analysis->row_ep_density);
              ekk_instance_.updateOperationResultDensity(
                  local_row_ep_density,
                  ekk_instance_.simplex_info_.row_ep_density);
            }
            if (ekk_instance_.analysis_.analyse_simplex_time) {
              analysis->simplexTimerStop(SimplexIzDseWtClock);
              analysis->simplexTimerStop(DseIzClock);
              double IzDseWtTT =
                  analysis->simplexTimerRead(SimplexIzDseWtClock);
              highsLogDev(options.log_options, HighsLogType::DETAILED,
                          "Computed %" HIGHSINT_FORMAT
                          " initial DSE weights in %gs\n",
                          solver_num_row, IzDseWtTT);
            }
          }
        } else {
          highsLogDev(
              options.log_options, HighsLogType::DETAILED,
              "solve:: Starting from B=I so unit initial DSE weights\n");
        }
      }
    }
    // Indicate that edge weights are known
    simplex_lp_status.has_dual_steepest_edge_weights = true;
  }
  // Resize the copy of scattered edge weights for backtracking
  simplex_info.backtracking_basis_edge_weights_.resize(solver_num_tot);

  // Compute the dual values
  ekk_instance_.computeDual();
  // Determine the number of dual infeasibilities, and hence the solve phase
  ekk_instance_.computeDualInfeasibleWithFlips();
  dualInfeasCount = simplex_info.num_dual_infeasibility;
  solvePhase = dualInfeasCount > 0 ? SOLVE_PHASE_1 : SOLVE_PHASE_2;
  if (ekkDebugOkForSolve(ekk_instance_, SimplexAlgorithm::DUAL, solvePhase,
                         ekk_instance_.scaled_model_status_) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  //
  // The major solving loop
  //
  while (solvePhase) {
    HighsInt it0 = ekk_instance_.iteration_count_;
    // When starting a new phase the (updated) dual objective function
    // value isn't known. Indicate this so that when the value
    // computed from scratch in rebuild() isn't checked against the
    // the updated value
    simplex_lp_status.has_dual_objective_value = false;
    if (solvePhase == SOLVE_PHASE_UNKNOWN) {
      // Reset the phase 2 bounds so that true number of dual
      // infeasibilities can be determined
      ekk_instance_.initialiseBound(SimplexAlgorithm::DUAL,
                                    SOLVE_PHASE_UNKNOWN);
      ekk_instance_.initialiseNonbasicValueAndMove();
      // Determine the number of dual infeasibilities, and hence the solve phase
      ekk_instance_.computeDualInfeasibleWithFlips();
      dualInfeasCount = simplex_info.num_dual_infeasibility;
      solvePhase = dualInfeasCount > 0 ? SOLVE_PHASE_1 : SOLVE_PHASE_2;
      if (simplex_info.backtracking_) {
        // Backtracking, so set the bounds and primal values
        ekk_instance_.initialiseBound(SimplexAlgorithm::DUAL, solvePhase);
        ekk_instance_.initialiseNonbasicValueAndMove();
        // Can now forget that we might have been backtracking
        simplex_info.backtracking_ = false;
      }
    }
    assert(solvePhase == SOLVE_PHASE_1 || solvePhase == SOLVE_PHASE_2);
    if (solvePhase == SOLVE_PHASE_1) {
      // Phase 1
      analysis->simplexTimerStart(SimplexDualPhase1Clock);
      solvePhase1();
      analysis->simplexTimerStop(SimplexDualPhase1Clock);
      simplex_info.dual_phase1_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else if (solvePhase == SOLVE_PHASE_2) {
      // Phase 2
      analysis->simplexTimerStart(SimplexDualPhase2Clock);
      solvePhase2();
      analysis->simplexTimerStop(SimplexDualPhase2Clock);
      simplex_info.dual_phase2_iteration_count +=
          (ekk_instance_.iteration_count_ - it0);
    } else {
      // Should only be SOLVE_PHASE_1 or SOLVE_PHASE_2
      scaled_model_status = HighsModelStatus::SOLVE_ERROR;
      return ekk_instance_.returnFromSolve(HighsStatus::Error);
    }
    // Return if bailing out from solve
    if (solve_bailout)
      return ekk_instance_.returnFromSolve(HighsStatus::Warning);
    // Can have all possible cases of solvePhase
    assert(solvePhase >= SOLVE_PHASE_MIN && solvePhase <= SOLVE_PHASE_MAX);
    // Look for scenarios when the major solving loop ends
    if (solvePhase == SOLVE_PHASE_ERROR) {
      // Solver error so return HighsStatus::Error
      assert(scaled_model_status == HighsModelStatus::SOLVE_ERROR);
      return ekk_instance_.returnFromSolve(HighsStatus::Error);
    }
    if (solvePhase == SOLVE_PHASE_EXIT) {
      // LP identified as not having an optimal solution
      assert(scaled_model_status == HighsModelStatus::PRIMAL_DUAL_INFEASIBLE ||
             scaled_model_status == HighsModelStatus::PRIMAL_INFEASIBLE);
      break;
    }
    if (solvePhase == SOLVE_PHASE_1 &&
        scaled_model_status == HighsModelStatus::DUAL_INFEASIBLE) {
      // Dual infeasibilities after phase 2 for a problem known to be dual
      // infeasible.
      break;
    }
    if (solvePhase == SOLVE_PHASE_CLEANUP) {
      // Dual infeasibilities after phase 2 for a problem not known to
      // be dual infeasible. Primal feasible with dual infeasibilities
      // so use primal simplex to clean up
      break;
    }
    // If solvePhase == SOLVE_PHASE_OPTIMAL == 0 then major solving
    // loop ends naturally since solvePhase is false
  }
  // If bailing out, should have returned already
  assert(!solve_bailout);
  // Should only have these cases
  assert(solvePhase == SOLVE_PHASE_EXIT || solvePhase == SOLVE_PHASE_UNKNOWN ||
         solvePhase == SOLVE_PHASE_OPTIMAL || solvePhase == SOLVE_PHASE_1 ||
         solvePhase == SOLVE_PHASE_CLEANUP);
  if (solvePhase == SOLVE_PHASE_1) {
    assert(scaled_model_status == HighsModelStatus::DUAL_INFEASIBLE);
    // Resolve case of LP that is dual infeasible (and not primal
    // feasible since that would yield solvePhase ==
    // SOLVE_PHASE_CLEANUP Looking to identify primal infeasiblilty or
    // primal unboundedness Cleanup with phase 1 for new primal code
    ekk_instance_.computePrimalObjectiveValue();

    highsLogUser(options.log_options, HighsLogType::ERROR,
                 "Primal simplex solver unavailable\n");
    //    if (scaled_model_status == HighsModelStatus::OPTIMAL) {
    //      if
    //      (ekk_instance_.simplex_info_..num_primal_infeasibility) {
    //        // Optimal with primal infeasibilities => primal infeasible
    //        assert(ekk_instance_.simplex_info_.num_primal_infeasibility
    //        > 0); scaled_model_status =
    //        HighsModelStatus::PRIMAL_DUAL_INFEASIBLE;
    //      }
    //    } else {
    //      // Should only be primal unbounded
    //      assert(scaled_model_status == HighsModelStatus::PRIMAL_UNBOUNDED);
    //    }
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  }
  if (solvePhase == SOLVE_PHASE_CLEANUP) {
    ekk_instance_.computePrimalObjectiveValue();
    if (options.dual_simplex_cleanup_strategy ==
        DUAL_SIMPLEX_CLEANUP_STRATEGY_NONE) {
      // No clean up. Dual simplex was optimal with perturbed costs,
      // so say that the scaled LP has been solved
      // optimally. Optimality (unlikely) for the unscaled LP will
      // still be assessed honestly, so leave it to the user to
      // deceide whether the solution can be accepted.
      scaled_model_status = HighsModelStatus::OPTIMAL;
    } else {
      // Use primal to clean up
      if (options.dual_simplex_cleanup_strategy) {
        analysis->simplexTimerStart(SimplexPrimalPhase2Clock);
        // Cleanup with primal code
        // Switch off any bound perturbation
        double save_primal_simplex_bound_perturbation_multiplier =
            simplex_info.primal_simplex_bound_perturbation_multiplier;
        simplex_info.primal_simplex_bound_perturbation_multiplier = 0;
        HEkkPrimal hEkkPrimal(ekk_instance_);
        hEkkPrimal.solve();
        // Restore any bound perturbation
        simplex_info.primal_simplex_bound_perturbation_multiplier =
            save_primal_simplex_bound_perturbation_multiplier;
        analysis->simplexTimerStop(SimplexPrimalPhase2Clock);
      }
    }
  }
  if (ekkDebugOkForSolve(ekk_instance_, SimplexAlgorithm::DUAL, solvePhase,
                         ekk_instance_.scaled_model_status_) ==
      HighsDebugStatus::kLogicalError)
    return ekk_instance_.returnFromSolve(HighsStatus::Error);
  ekk_instance_.computePrimalObjectiveValue();
  return ekk_instance_.returnFromSolve(HighsStatus::OK);
}

void HEkkDual::options() {
  // Set solver options from simplex options

  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;

  interpretDualEdgeWeightStrategy(simplex_info.dual_edge_weight_strategy);

  // Copy values of simplex solver options to dual simplex options
  primal_feasibility_tolerance =
      ekk_instance_.options_.primal_feasibility_tolerance;
  dual_feasibility_tolerance =
      ekk_instance_.options_.dual_feasibility_tolerance;
  dual_objective_value_upper_bound =
      ekk_instance_.options_.dual_objective_value_upper_bound;
  //  perturb_costs = simplex_info.perturb_costs;
  //  iterationLimit = simplex_info.iterationLimit;

  // Set values of internal options
}

void HEkkDual::init() {
  // Copy size, matrix and factor

  solver_num_col = ekk_instance_.simplex_lp_.numCol_;
  solver_num_row = ekk_instance_.simplex_lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;

  matrix = &ekk_instance_.matrix_;
  factor = &ekk_instance_.factor_;
  analysis = &ekk_instance_.analysis_;

  // Copy pointers
  jMove = &ekk_instance_.simplex_basis_.nonbasicMove_[0];
  workDual = &ekk_instance_.simplex_info_.workDual_[0];
  workValue = &ekk_instance_.simplex_info_.workValue_[0];
  workRange = &ekk_instance_.simplex_info_.workRange_[0];
  baseLower = &ekk_instance_.simplex_info_.baseLower_[0];
  baseUpper = &ekk_instance_.simplex_info_.baseUpper_[0];
  baseValue = &ekk_instance_.simplex_info_.baseValue_[0];

  // Copy tolerances
  Tp = primal_feasibility_tolerance;
  Td = dual_feasibility_tolerance;

  // Setup local vectors
  col_DSE.setup(solver_num_row);
  col_BFRT.setup(solver_num_row);
  col_aq.setup(solver_num_row);
  row_ep.setup(solver_num_row);
  row_ap.setup(solver_num_col);
  // Setup other buffers
  dualRow.setup();
  dualRHS.setup();
}

void HEkkDual::initParallel() {
  // Identify the (current) number of HiGHS tasks to be used
  const HighsInt num_threads = ekk_instance_.simplex_info_.num_threads;

  // Initialize for tasks
  if (ekk_instance_.simplex_info_.simplex_strategy ==
      SIMPLEX_STRATEGY_DUAL_TASKS) {
    const HighsInt pass_num_slice = num_threads - 2;
    assert(pass_num_slice > 0);
    if (pass_num_slice <= 0) {
      highsLogUser(ekk_instance_.options_.log_options, HighsLogType::WARNING,
                   "SIP trying to use using %" HIGHSINT_FORMAT
                   " slices due to number of "
                   "threads (%" HIGHSINT_FORMAT
                   ") being too small: results unpredictable\n",
                   pass_num_slice, num_threads);
    }
    initSlice(pass_num_slice);
  }

  // Initialize for multi
  if (ekk_instance_.simplex_info_.simplex_strategy ==
      SIMPLEX_STRATEGY_DUAL_MULTI) {
    multi_num = num_threads;
    if (multi_num < 1) multi_num = 1;
    if (multi_num > kHighsThreadLimit) multi_num = kHighsThreadLimit;
    for (HighsInt i = 0; i < multi_num; i++) {
      multi_choice[i].row_ep.setup(solver_num_row);
      multi_choice[i].col_aq.setup(solver_num_row);
      multi_choice[i].col_BFRT.setup(solver_num_row);
    }
    const HighsInt pass_num_slice = max(multi_num - 1, HighsInt{1});
    assert(pass_num_slice > 0);
    if (pass_num_slice <= 0) {
      highsLogUser(ekk_instance_.options_.log_options, HighsLogType::WARNING,
                   "PAMI trying to use using %" HIGHSINT_FORMAT
                   " slices due to number of "
                   "threads (%" HIGHSINT_FORMAT
                   ") being too small: results unpredictable\n",
                   pass_num_slice, num_threads);
    }
    initSlice(pass_num_slice);
  }
  multi_iteration = 0;
  //  string partitionFile = model->strOption[STROPT_PARTITION_FILE];
  //  if (partitionFile.size())
  //  {
  //    dualRHS.setup_partition(partitionFile.c_str());
  //  }
}

void HEkkDual::initSlice(const HighsInt initial_num_slice) {
  // Number of slices
  slice_num = initial_num_slice;
  if (slice_num < 1) slice_num = 1;
  assert(slice_num <= HIGHS_SLICED_LIMIT);
  if (slice_num > HIGHS_SLICED_LIMIT) {
    highsLogUser(ekk_instance_.options_.log_options, HighsLogType::WARNING,
                 "WARNING: %" HIGHSINT_FORMAT
                 " = slice_num > HIGHS_SLICED_LIMIT = %" HIGHSINT_FORMAT
                 " so truncating "
                 "slice_num\n",
                 slice_num, HIGHS_SLICED_LIMIT);
    slice_num = HIGHS_SLICED_LIMIT;
  }

  // Alias to the matrix
  const HighsInt* Astart = matrix->getAstart();
  const HighsInt* Aindex = matrix->getAindex();
  const double* Avalue = matrix->getAvalue();
  const HighsInt AcountX = Astart[solver_num_col];

  // Figure out partition weight
  double sliced_countX = AcountX / slice_num;
  slice_start[0] = 0;
  for (HighsInt i = 0; i < slice_num - 1; i++) {
    HighsInt endColumn = slice_start[i] + 1;  // At least one column
    HighsInt endX = Astart[endColumn];
    HighsInt stopX = (i + 1) * sliced_countX;
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
  vector<HighsInt> sliced_Astart;
  for (HighsInt i = 0; i < slice_num; i++) {
    // The matrix
    HighsInt mystart = slice_start[i];
    HighsInt mycount = slice_start[i + 1] - mystart;
    HighsInt mystartX = Astart[mystart];
    sliced_Astart.resize(mycount + 1);
    for (HighsInt k = 0; k <= mycount; k++)
      sliced_Astart[k] = Astart[k + mystart] - mystartX;
    slice_matrix[i].setup_lgBs(mycount, solver_num_row, &sliced_Astart[0],
                               Aindex + mystartX, Avalue + mystartX);

    // The row_ap and its packages
    slice_row_ap[i].setup(mycount);
    slice_dualRow[i].setupSlice(mycount);
  }
}

void HEkkDual::solvePhase1() {
  // Performs dual phase 1 iterations. Returns solvePhase with value
  //
  // SOLVE_PHASE_ERROR => Solver error
  //
  // SOLVE_PHASE_UNKNOWN => Back-tracking due to singularity
  //
  // SOLVE_PHASE_1 => Dual infeasibility suspected, but have to go out
  // and back in to solvePhase1 to perform fresh rebuild. Also if
  // bailing out due to reaching time/iteration limit.
  //
  // SOLVE_PHASE_2 => Continue with dual phase 2 iterations

  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  HighsModelStatus& scaled_model_status = ekk_instance_.scaled_model_status_;
  // When starting a new phase the (updated) dual objective function
  // value isn't known. Indicate this so that when the value computed
  // from scratch in build() isn't checked against the the updated
  // value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Set rebuild_reason so that it's assigned when first tested
  rebuild_reason = REBUILD_REASON_NO;
  // Set solvePhase = SOLVE_PHASE_1 and solve_bailout = false so they are set if
  // solvePhase1() is called directly
  solvePhase = SOLVE_PHASE_1;
  solve_bailout = false;
  if (bailoutOnTimeIterations()) return;
  // Report the phase start
  highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
              "dual-phase-1-start\n");
  // Switch to dual phase 1 bounds
  ekk_instance_.initialiseBound(SimplexAlgorithm::DUAL, solvePhase);
  ekk_instance_.initialiseNonbasicValueAndMove();

  // If there's no backtracking basis, save the initial basis in case of
  // backtracking
  if (!simplex_info.valid_backtracking_basis_)
    ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  analysis->simplexTimerStart(IterateClock);
  for (;;) {
    analysis->simplexTimerStart(IterateDualRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IterateDualRebuildClock);
    if (solvePhase == SOLVE_PHASE_ERROR) {
      scaled_model_status = HighsModelStatus::SOLVE_ERROR;
      return;
    }
    if (solvePhase == SOLVE_PHASE_UNKNOWN) {
      // If backtracking, may change phase, so drop out
      analysis->simplexTimerStop(IterateClock);
      return;
    }
    if (bailoutOnTimeIterations()) break;
    for (;;) {
      if (debugDualSimplex("Before iteration") ==
          HighsDebugStatus::kLogicalError) {
        solvePhase = SOLVE_PHASE_ERROR;
        return;
      }
      switch (simplex_info.simplex_strategy) {
        default:
        case SIMPLEX_STRATEGY_DUAL_PLAIN:
          iterate();
          break;
        case SIMPLEX_STRATEGY_DUAL_TASKS:
          iterateTasks();
          break;
        case SIMPLEX_STRATEGY_DUAL_MULTI:
          iterateMulti();
          break;
      }
      if (bailoutOnTimeIterations()) break;
      if (rebuild_reason) break;
    }
    if (solve_bailout) break;
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) break;
  }

  analysis->simplexTimerStop(IterateClock);
  // Possibly return due to bailing out, having now stopped
  // IterateClock
  if (bailoutReturn()) return;

  // If bailing out, should have done so already
  assert(!solve_bailout);
  // Assess outcome of dual phase 1
  if (row_out == -1) {
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
                "dual-phase-1-optimal\n");
    // Optimal in phase 1
    if (simplex_info.dual_objective_value == 0) {
      // Zero phase 1 objective so go to phase 2
      //
      // This is the usual way to exit phase 1. Although the test
      // looks ambitious, the dual objective is the sum of products of
      // primal and dual values for nonbasic variables. For dual
      // simplex phase 1, the primal bounds are set so that when the
      // dual value is feasible, the primal value is set to
      // zero. Otherwise the value is +1/-1 according to the required
      // sign of the dual, except for free variables, where the bounds
      // are [-1000, 1000].
      //
      // OK if costs are perturbed, since they remain perturbed in phase 2 until
      // the final clean-up
      solvePhase = SOLVE_PHASE_2;
    } else {
      // A negative dual objective value at an optimal solution of
      // phase 1 means that there are dual infeasibilities. If the
      // objective value is very negative, then it's clear
      // enough. However, if it's small, it could be the sum of
      // values, all of which are smaller than the dual deasibility
      // tolerance. Plus there may be cost perturbations to remove
      assessPhase1Optimality();
    }
  } else if (rebuild_reason == REBUILD_REASON_CHOOSE_COLUMN_FAIL) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = SOLVE_PHASE_ERROR;
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                "dual-phase-1-not-solved\n");
    scaled_model_status = HighsModelStatus::SOLVE_ERROR;
  } else if (variable_in == -1) {
    // We got dual phase 1 unbounded - strange
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                "dual-phase-1-unbounded\n");
    if (ekk_instance_.simplex_info_.costs_perturbed) {
      // Clean up perturbation
      cleanup();
      highsLogUser(ekk_instance_.options_.log_options, HighsLogType::WARNING,
                   "Cleaning up cost perturbation when unbounded in phase 1\n");
      if (dualInfeasCount == 0) {
        // No dual infeasibilities and (since unbounded) at least zero
        // phase 1 objective so go to phase 2
        solvePhase = SOLVE_PHASE_2;
      }
    } else {
      // Report strange issues
      solvePhase = SOLVE_PHASE_ERROR;
      highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                  "dual-phase-1-not-solved\n");
      scaled_model_status = HighsModelStatus::SOLVE_ERROR;
    }
  }

  // Debug here since not all simplex data will be correct until after
  // rebuild() when switching to Phase 2
  //
  // Also have to avoid debug when the model status is not set and
  // there are dual infeasibilities, since this happens legitimately
  // when the LP is dual infeasible. However, the model status can't
  // be set to dual infeasible until perturbations have been removed.
  //
  const bool no_debug =
      ekk_instance_.simplex_info_.num_dual_infeasibility > 0 &&
      scaled_model_status == HighsModelStatus::NOTSET;
  if (!no_debug) {
    if (debugDualSimplex("End of solvePhase1") ==
        HighsDebugStatus::kLogicalError) {
      solvePhase = SOLVE_PHASE_ERROR;
      return;
    }
  }

  if (solvePhase == SOLVE_PHASE_2) {
    // Moving to phase 2 so allow cost perturbation. It may have been
    // prevented to avoid cleanup-perturbation loops when optimal in
    // phase 1
    simplex_info.allow_cost_perturbation = true;
    ekk_instance_.initialiseBound(SimplexAlgorithm::DUAL, solvePhase, true);
    ekk_instance_.initialiseNonbasicValueAndMove();
  }
  return;
}

void HEkkDual::solvePhase2() {
  // Performs dual phase 2 iterations. Returns solvePhase with value
  //
  // SOLVE_PHASE_ERROR => Solver error
  //
  // SOLVE_PHASE_EXIT => LP identified as not having an optimal solution
  //
  // SOLVE_PHASE_UNKNOWN => Back-tracking due to singularity
  //
  // SOLVE_PHASE_OPTIMAL => Primal feasible and no dual infeasibilities =>
  // Optimal
  //
  // SOLVE_PHASE_1 => Primal feasible and dual infeasibilities for a
  // problem known to be dual infeasible => Use primal phase 2 to
  // determine primal unboundedness.
  //
  // SOLVE_PHASE_2 => Dual unboundedness suspected, but have to go out
  // and back in to solvePhase2 to perform fresh rebuild. Also if
  // bailing out due to reaching time/iteration limit or dual
  // objective
  //
  // SOLVE_PHASE_CLEANUP => Contrinue with primal phase 2 iterations to clean up
  // dual infeasibilities
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  HighsModelStatus& scaled_model_status = ekk_instance_.scaled_model_status_;
  // When starting a new phase the (updated) dual objective function
  // value isn't known. Indicate this so that when the value computed
  // from scratch in build() isn't checked against the the updated
  // value
  simplex_lp_status.has_primal_objective_value = false;
  simplex_lp_status.has_dual_objective_value = false;
  // Set rebuild_reason so that it's assigned when first tested
  rebuild_reason = REBUILD_REASON_NO;
  // Set solvePhase = SOLVE_PHASE_2 and solve_bailout = false so they are set if
  // solvePhase2() is called directly
  solvePhase = SOLVE_PHASE_2;
  solve_bailout = false;
  if (bailoutOnTimeIterations()) return;
  // Report the phase start
  highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
              "dual-phase-2-start\n");
  // Collect free variables
  dualRow.createFreelist();

  // If there's no backtracking basis Save the initial basis in case of
  // backtracking
  if (!simplex_info.valid_backtracking_basis_)
    ekk_instance_.putBacktrackingBasis();

  // Main solving structure
  analysis->simplexTimerStart(IterateClock);
  for (;;) {
    // Outer loop of solvePhase2()
    // Rebuild all values, reinverting B if updates have been performed
    analysis->simplexTimerStart(IterateDualRebuildClock);
    rebuild();
    analysis->simplexTimerStop(IterateDualRebuildClock);
    if (solvePhase == SOLVE_PHASE_ERROR) {
      scaled_model_status = HighsModelStatus::SOLVE_ERROR;
      return;
    }
    if (solvePhase == SOLVE_PHASE_UNKNOWN) {
      // If backtracking, may change phase, so drop out
      analysis->simplexTimerStop(IterateClock);
      return;
    }
    if (bailoutOnTimeIterations()) break;
    if (bailoutOnDualObjective()) break;
    if (dualInfeasCount > 0) break;
    for (;;) {
      // Inner loop of solvePhase2()
      // Performs one iteration in case SIMPLEX_STRATEGY_DUAL_PLAIN:
      if (debugDualSimplex("Before iteration") ==
          HighsDebugStatus::kLogicalError) {
        solvePhase = SOLVE_PHASE_ERROR;
        return;
      }
      switch (simplex_info.simplex_strategy) {
        default:
        case SIMPLEX_STRATEGY_DUAL_PLAIN:
          iterate();
          break;
        case SIMPLEX_STRATEGY_DUAL_TASKS:
          iterateTasks();
          break;
        case SIMPLEX_STRATEGY_DUAL_MULTI:
          iterateMulti();
          break;
      }
      if (bailoutOnTimeIterations()) break;
      if (bailoutOnDualObjective()) break;
      if (rebuild_reason) break;
    }
    if (solve_bailout) break;
    // If the data are fresh from rebuild(), break out of
    // the outer loop to see what's ocurred
    // Was:	if (simplex_info.update_count == 0) break;
    if (simplex_lp_status.has_fresh_rebuild) break;
  }
  analysis->simplexTimerStop(IterateClock);
  // Possibly return due to bailing out, having now stopped
  // IterateClock
  if (bailoutReturn()) return;

  // If bailing out, should have done so already
  assert(!solve_bailout);
  // Assess outcome of dual phase 2
  if (dualInfeasCount > 0) {
    // There are dual infeasiblities so possibly switch to Phase 1 and
    // return. "Possibly" because, if dual infeasibility has already
    // been shown, primal simplex is used to distinguish primal
    // unboundedness from primal infeasibility
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
                "dual-phase-2-found-free\n");
    solvePhase = SOLVE_PHASE_1;
  } else if (row_out == -1) {
    // There is no candidate in CHUZR, even after rebuild so probably optimal
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
                "dual-phase-2-optimal\n");
    // Remove any cost perturbations and see if basis is still dual feasible
    cleanup();
    if (dualInfeasCount > 0) {
      // There are dual infeasiblities, so consider performing primal
      // simplex iterations to get dual feasibility
      solvePhase = SOLVE_PHASE_CLEANUP;
    } else {
      // There are no dual infeasiblities so optimal!
      solvePhase = SOLVE_PHASE_OPTIMAL;
      highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
                  "problem-optimal\n");
      scaled_model_status = HighsModelStatus::OPTIMAL;
    }
  } else if (rebuild_reason == REBUILD_REASON_CHOOSE_COLUMN_FAIL) {
    // chooseColumn has failed
    // Behave as "Report strange issues" below
    solvePhase = SOLVE_PHASE_ERROR;
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                "dual-phase-2-not-solved\n");
    scaled_model_status = HighsModelStatus::SOLVE_ERROR;
  } else if (variable_in == -1) {
    // There is no candidate in CHUZC, so probably dual unbounded
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                "dual-phase-2-unbounded\n");
    if (ekk_instance_.simplex_info_.costs_perturbed) {
      // If the costs have been perturbed, clean up and return
      cleanup();
    } else {
      // The costs have not been perturbed, so dual unbounded---and
      // hence primal infeasible.
      solvePhase = SOLVE_PHASE_EXIT;
      if (scaled_model_status == HighsModelStatus::DUAL_INFEASIBLE) {
        // LP has already been shown to be dual infeasible, so a cost
        // shift was necessary to achieve dual feasibility in phase
        // 2. In this case dual unboundedness (hopefully) implies
        // primal and dual infeasiblility. Without a primal simplex
        // solver, this is the only way to identify primal
        // infeasiblilty for and LP that's also dual infeasible. With
        // a primal simplex solver, once dual infeasiblility has been
        // identified, use primal phase 1 to identify primal
        // infeasiblilty.
        highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                    "problem-primal-dual-infeasible\n");
        scaled_model_status = HighsModelStatus::PRIMAL_DUAL_INFEASIBLE;
      } else {
        // Dual feasible and dual unbounded, so save dual ray
        saveDualRay();
        // Model status should be unset?
        assert(scaled_model_status == HighsModelStatus::NOTSET);
        highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                    "problem-primal-infeasible\n");
        scaled_model_status = HighsModelStatus::PRIMAL_INFEASIBLE;
      }
    }
  }
  // Before primal simplex clean-up there will be dual infeasibilities
  if (solvePhase != SOLVE_PHASE_CLEANUP) {
    if (debugDualSimplex("End of solvePhase2") ==
        HighsDebugStatus::kLogicalError) {
      solvePhase = SOLVE_PHASE_ERROR;
      return;
    }
  }
  return;
}

void HEkkDual::rebuild() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance_.simplex_lp_status_;
  // Save history information
  // Move this to Simplex class once it's created
  //  record_pivots(-1, -1, 0);  // Indicate REINVERT

  const HighsInt reason_for_rebuild = rebuild_reason;
  rebuild_reason = REBUILD_REASON_NO;
  // Possibly Rebuild ekk_instance_.factor_
  bool reInvert = simplex_info.update_count > 0;
  if (reInvert) {
    // Get a nonsingular inverse if possible. One of three things
    // happens: Current basis is nonsingular; Current basis is
    // singular and last nonsingular basis is refactorized as
    // nonsingular - or found singular. Latter is code failure.
    if (!ekk_instance_.getNonsingularInverse(solvePhase)) {
      solvePhase = SOLVE_PHASE_ERROR;
      return;
    }
  }

  if (!ekk_instance_.simplex_lp_status_.has_matrix) {
    // Don't have the matrix either row-wise or col-wise, so
    // reinitialise it
    assert(simplex_info.backtracking_);
    HighsLp& simplex_lp = ekk_instance_.simplex_lp_;
    analysis->simplexTimerStart(matrixSetupClock);
    ekk_instance_.matrix_.setup(simplex_lp.numCol_, simplex_lp.numRow_,
                                &simplex_lp.Astart_[0], &simplex_lp.Aindex_[0],
                                &simplex_lp.Avalue_[0],
                                &ekk_instance_.simplex_basis_.nonbasicFlag_[0]);
    simplex_lp_status.has_matrix = true;
    analysis->simplexTimerStop(matrixSetupClock);
  }
  // Record whether the update objective value should be tested. If
  // the objective value is known, then the updated objective value
  // should be correct - once the correction due to recomputing the
  // dual values has been applied.
  //
  // Note that computePrimalObjectiveValue sets
  // has_primal_objective_value
  const bool check_updated_objective_value =
      simplex_lp_status.has_dual_objective_value;
  double previous_dual_objective_value;
  if (check_updated_objective_value) {
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase,
    //    "Before computeDual");
    previous_dual_objective_value = simplex_info.updated_dual_objective_value;
  } else {
    // Reset the knowledge of previous objective values
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, -1, "");
  }
  // Recompute dual solution
  ekk_instance_.computeDual();

  if (simplex_info.backtracking_) {
    // If backtracking, may change phase, so drop out
    solvePhase = SOLVE_PHASE_UNKNOWN;
    return;
  }
  analysis->simplexTimerStart(CorrectDualClock);
  ekk_instance_.correctDual(&dualInfeasCount);
  analysis->simplexTimerStop(CorrectDualClock);

  // Recompute primal solution
  ekk_instance_.computePrimal();

  // Collect primal infeasible as a list
  analysis->simplexTimerStart(CollectPrIfsClock);
  dualRHS.createArrayOfPrimalInfeasibilities();
  dualRHS.createInfeasList(analysis->col_aq_density);
  analysis->simplexTimerStop(CollectPrIfsClock);

  // Dual objective section
  //
  ekk_instance_.computeDualObjectiveValue(solvePhase);

  if (check_updated_objective_value) {
    // Apply the objective value correction due to computing duals
    // from scratch.
    const double dual_objective_value_correction =
        simplex_info.dual_objective_value - previous_dual_objective_value;
    simplex_info.updated_dual_objective_value +=
        dual_objective_value_correction;
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm);
  }
  // Now that there's a new dual_objective_value, reset the updated
  // value
  simplex_info.updated_dual_objective_value = simplex_info.dual_objective_value;

  if (!simplex_info.run_quiet) {
    // Report the primal infeasiblities
    ekk_instance_.computeSimplexPrimalInfeasible();
    if (solvePhase == SOLVE_PHASE_1) {
      // In phase 1, report the simplex LP dual infeasiblities
      ekk_instance_.computeSimplexLpDualInfeasible();
    } else {
      // In phase 2, report the simplex dual infeasiblities
      ekk_instance_.computeSimplexDualInfeasible();
    }
    reportRebuild(reason_for_rebuild);
  }

  ekk_instance_.build_syntheticTick_ = factor->build_syntheticTick;
  ekk_instance_.total_syntheticTick_ = 0;

  // Dual simplex doesn't maintain the number of primal
  // infeasiblities, so set it to an illegal value now
  ekk_instance_.invalidatePrimalInfeasibilityRecord();
  // Although dual simplex should always be dual feasible,
  // infeasiblilities are only corrected in rebuild
  ekk_instance_.invalidateDualInfeasibilityRecord();

  // Data are fresh from rebuild
  simplex_lp_status.has_fresh_rebuild = true;
}

void HEkkDual::cleanup() {
  highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
              "dual-cleanup-shift\n");
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  // Remove perturbation and don't permit further perturbation
  ekk_instance_.initialiseCost(SimplexAlgorithm::DUAL, SOLVE_PHASE_UNKNOWN);
  simplex_info.allow_cost_perturbation = false;
  // No solvePhase term in initialiseBound is surely an omission -
  // when cleanup called in phase 1
  ekk_instance_.initialiseBound(SimplexAlgorithm::DUAL, solvePhase);
  // Possibly take a copy of the original duals before recomputing them
  vector<double> original_workDual;
  if (ekk_instance_.options_.highs_debug_level > kHighsDebugLevelCheap)
    original_workDual = simplex_info.workDual_;
  // Compute the dual values
  ekk_instance_.computeDual();
  // Possibly analyse the change in duals
  //  debugCleanup(ekk_instance_, original_workDual);
  // Compute the dual infeasibilities
  ekk_instance_.computeSimplexDualInfeasible();
  dualInfeasCount = ekk_instance_.simplex_info_.num_dual_infeasibility;

  // Compute the dual objective value
  ekk_instance_.computeDualObjectiveValue(solvePhase);
  // Now that there's a new dual_objective_value, reset the updated
  // value
  simplex_info.updated_dual_objective_value = simplex_info.dual_objective_value;

  if (!simplex_info.run_quiet) {
    // Report the primal infeasiblities
    ekk_instance_.computeSimplexPrimalInfeasible();
    // In phase 1, report the simplex LP dual infeasiblities
    // In phase 2, report the simplex dual infeasiblities (known)
    if (solvePhase == SOLVE_PHASE_1)
      ekk_instance_.computeSimplexLpDualInfeasible();
    reportRebuild();
  }
}

void HEkkDual::iterate() {
  // This is the main teration loop for dual revised simplex. All the
  // methods have as their first line if (rebuild_reason) return;, where
  // rebuild_reason is, for example, set to 1 when CHUZR finds no
  // candidate. This causes a break from the inner loop of
  // solve_phase% and, hence, a call to rebuild()

  // Reporting:
  // Row-wise matrix after update in updateMatrix(variable_in, variable_out);
  analysis->simplexTimerStart(IterateChuzrClock);
  chooseRow();
  analysis->simplexTimerStop(IterateChuzrClock);

  analysis->simplexTimerStart(IterateChuzcClock);
  chooseColumn(&row_ep);
  analysis->simplexTimerStop(IterateChuzcClock);

  analysis->simplexTimerStart(IterateFtranClock);
  updateFtranBFRT();

  // updateFtran(); computes the pivotal column in the data structure "column"
  updateFtran();

  // updateFtranDSE performs the DSE FTRAN on pi_p
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE)
    updateFtranDSE(&row_ep);
  analysis->simplexTimerStop(IterateFtranClock);

  // updateVerify() Checks row-wise pivot against column-wise pivot for
  // numerical trouble
  analysis->simplexTimerStart(IterateVerifyClock);
  updateVerify();
  analysis->simplexTimerStop(IterateVerifyClock);

  // updateDual() Updates the dual values
  analysis->simplexTimerStart(IterateDualClock);
  updateDual();
  analysis->simplexTimerStop(IterateDualClock);

  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "Before
  //  updatePrimal");
  // updatePrimal(&row_ep); Updates the primal values and the edge weights
  analysis->simplexTimerStart(IteratePrimalClock);
  updatePrimal(&row_ep);
  analysis->simplexTimerStop(IteratePrimalClock);
  // After primal update in dual simplex the primal objective value is not known
  ekk_instance_.simplex_lp_status_.has_primal_objective_value = false;
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "After
  //  updatePrimal");

  // Update the basis representation
  analysis->simplexTimerStart(IteratePivotsClock);
  updatePivots();
  analysis->simplexTimerStop(IteratePivotsClock);

  if (new_devex_framework) {
    // Initialise new Devex framework
    analysis->simplexTimerStart(IterateDevexIzClock);
    initialiseDevexFramework();
    analysis->simplexTimerStop(IterateDevexIzClock);
  }

  // Analyse the iteration: possibly report; possibly switch strategy
  iterationAnalysis();
}

void HEkkDual::iterateTasks() {
  slice_PRICE = 1;

  // Group 1
  chooseRow();

  // Disable slice when too sparse
  if (1.0 * row_ep.count / solver_num_row < 0.01) slice_PRICE = 0;

  analysis->simplexTimerStart(Group1Clock);
#pragma omp parallel
#pragma omp single
  {
#pragma omp task
    {
      col_DSE.copy(&row_ep);
      updateFtranDSE(&col_DSE);
    }
#pragma omp task
    {
      if (slice_PRICE)
        chooseColumnSlice(&row_ep);
      else
        chooseColumn(&row_ep);
#pragma omp task
      updateFtranBFRT();
#pragma omp task
      updateFtran();
#pragma omp taskwait
    }
  }
  analysis->simplexTimerStop(Group1Clock);

  updateVerify();
  updateDual();
  updatePrimal(&col_DSE);
  updatePivots();
}

void HEkkDual::iterationAnalysisData() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  analysis->simplex_strategy = simplex_info.simplex_strategy;
  analysis->edge_weight_mode = dual_edge_weight_mode;
  analysis->solve_phase = solvePhase;
  analysis->simplex_iteration_count = ekk_instance_.iteration_count_;
  analysis->devex_iteration_count = num_devex_iterations;
  analysis->pivotal_row_index = row_out;
  analysis->leaving_variable = variable_out;
  analysis->entering_variable = variable_in;
  analysis->rebuild_reason = rebuild_reason;
  analysis->reduced_rhs_value = 0;
  analysis->reduced_cost_value = 0;
  analysis->edge_weight = 0;
  analysis->primal_delta = delta_primal;
  analysis->primal_step = theta_primal;
  analysis->dual_step = theta_dual;
  analysis->pivot_value_from_column = alpha_col;
  analysis->pivot_value_from_row = alpha_row;
  analysis->factor_pivot_threshold = simplex_info.factor_pivot_threshold;
  analysis->numerical_trouble = numericalTrouble;
  analysis->objective_value = simplex_info.updated_dual_objective_value;
  // Since maximization is achieved by minimizing the LP with negated
  // costs, in phase 2 the dual objective value is negated, so flip
  // its sign according to the LP sense
  if (solvePhase == SOLVE_PHASE_2)
    analysis->objective_value *= (HighsInt)ekk_instance_.simplex_lp_.sense_;
  analysis->num_primal_infeasibility = simplex_info.num_primal_infeasibility;
  analysis->sum_primal_infeasibility = simplex_info.sum_primal_infeasibility;
  if (solvePhase == SOLVE_PHASE_1) {
    analysis->num_dual_infeasibility =
        analysis->num_dual_phase_1_lp_dual_infeasibility;
    analysis->sum_dual_infeasibility =
        analysis->sum_dual_phase_1_lp_dual_infeasibility;
  } else {
    analysis->num_dual_infeasibility = simplex_info.num_dual_infeasibility;
    analysis->sum_dual_infeasibility = simplex_info.sum_dual_infeasibility;
  }
  if ((dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) &&
      (num_devex_iterations == 0))
    analysis->num_devex_framework++;
}

void HEkkDual::iterationAnalysis() {
  // Possibly report on the iteration
  iterationAnalysisData();
  analysis->iterationReport();

  // Possibly switch from DSE to Devex
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    const bool switch_to_devex = ekk_instance_.switchToDevex();
    if (switch_to_devex) {
      dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
      // Using dual Devex edge weights, so set up the first framework
      ekk_instance_.simplex_info_.devex_index_.assign(solver_num_tot, 0);
      initialiseDevexFramework();
    }
  }
  if (analysis->analyse_simplex_data) analysis->iterationRecord();
}

void HEkkDual::reportRebuild(const HighsInt reason_for_rebuild) {
  analysis->simplexTimerStart(ReportRebuildClock);
  iterationAnalysisData();
  analysis->rebuild_reason = reason_for_rebuild;
  analysis->invertReport();
  analysis->simplexTimerStop(ReportRebuildClock);
}

void HEkkDual::chooseRow() {
  // Choose the index of a row to leave the basis (CHUZR)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  // Choose candidates repeatedly until candidate is OK or optimality is
  // detected
  for (;;) {
    // Choose the index of a good row to leave the basis
    dualRHS.chooseNormal(&row_out);
    if (row_out == -1) {
      // No index found so may be dual optimal.
      rebuild_reason = REBUILD_REASON_POSSIBLY_OPTIMAL;
      return;
    }
    // Compute pi_p = B^{-T}e_p in row_ep
    analysis->simplexTimerStart(BtranClock);
    // Set up RHS for BTRAN
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = row_out;
    row_ep.array[row_out] = 1;
    row_ep.packFlag = true;
    if (analysis->analyse_simplex_data)
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep,
                                      analysis->row_ep_density);
    // Perform BTRAN
    factor->btran(row_ep, analysis->row_ep_density,
                  analysis->pointer_serial_factor_clocks);
    if (analysis->analyse_simplex_data)
      analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_EP, row_ep);
    analysis->simplexTimerStop(BtranClock);
    // Verify DSE weight
    if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
      // For DSE, see how accurate the updated weight is
      // Save the updated weight
      double updated_edge_weight = dualRHS.workEdWt[row_out];
      // Compute the weight from row_ep and over-write the updated weight
      computed_edge_weight = dualRHS.workEdWt[row_out] = row_ep.norm2();
      // If the weight error is acceptable then break out of the
      // loop. All we worry about is accepting rows with weights
      // which are not too small, since this can make the row look
      // unreasonably attractive
      if (acceptDualSteepestEdgeWeight(updated_edge_weight)) break;
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
  variable_out = ekk_instance_.simplex_basis_.basicIndex_[row_out];
  // Record the change in primal variable associated with the move to the bound
  // being violated
  if (baseValue[row_out] < baseLower[row_out]) {
    // Below the lower bound so set delta_primal = value - LB < 0
    delta_primal = baseValue[row_out] - baseLower[row_out];
  } else {
    // Above the upper bound so set delta_primal = value - UB > 0
    delta_primal = baseValue[row_out] - baseUpper[row_out];
  }
  // Set move_out to be -1 if delta_primal<0, otherwise +1 (since
  // delta_primal>0)
  move_out = delta_primal < 0 ? -1 : 1;
  // Update the record of average row_ep (pi_p) density. This ignores
  // any BTRANs done for skipped candidates
  const double local_row_ep_density = (double)row_ep.count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_ep_density,
                                         analysis->row_ep_density);
  ekk_instance_.updateOperationResultDensity(
      local_row_ep_density, ekk_instance_.simplex_info_.row_ep_density);
}

bool HEkkDual::acceptDualSteepestEdgeWeight(const double updated_edge_weight) {
  // Accept the updated weight if it is at least a quarter of the
  // computed weight. Excessively large updated weights don't matter!
  const double accept_weight_threshold = 0.25;
  const bool accept_weight =
      updated_edge_weight >= accept_weight_threshold * computed_edge_weight;
  //  if (analysis->analyse_simplex_data)
  ekk_instance_.assessDSEWeightError(computed_edge_weight, updated_edge_weight);
  analysis->dualSteepestEdgeWeightError(computed_edge_weight,
                                        updated_edge_weight);
  return accept_weight;
}

bool HEkkDual::newDevexFramework(const double updated_edge_weight) {
  // Analyse the Devex weight to determine whether a new framework
  // should be set up
  double devex_ratio = max(updated_edge_weight / computed_edge_weight,
                           computed_edge_weight / updated_edge_weight);
  HighsInt i_te = solver_num_row / minRlvNumberDevexIterations;
  i_te = max(minAbsNumberDevexIterations, i_te);
  // Square maxAllowedDevexWeightRatio due to keeping squared
  // weights
  const double accept_ratio_threshold =
      maxAllowedDevexWeightRatio * maxAllowedDevexWeightRatio;
  const bool accept_ratio = devex_ratio <= accept_ratio_threshold;
  const bool accept_it = num_devex_iterations <= i_te;
  bool return_new_devex_framework;
  return_new_devex_framework = !accept_ratio || !accept_it;
  /*
  if (return_new_devex_framework) {
    printf("New Devex framework: (Iter %" HIGHSINT_FORMAT ") updated weight =
  %11.4g; computed weight = %11.4g; Devex ratio = %11.4g\n",
           ekk_instance_.iteration_count_,
           updated_edge_weight, computed_edge_weight, devex_ratio);
    return true;
  }
  */
  return return_new_devex_framework;
}

void HEkkDual::chooseColumn(HVector* row_ep) {
  // Compute pivot row (PRICE) and choose the index of a column to enter the
  // basis (CHUZC)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  //
  // PRICE
  //
  ekk_instance_.tableauRowPrice(*row_ep, row_ap);
  //
  // CHUZC
  //
  // Section 0: Clear data and call createFreemove to set a value of
  // nonbasicMove for all free columns to prevent their dual values
  // from being changed.
  analysis->simplexTimerStart(Chuzc0Clock);
  dualRow.clear();
  dualRow.workDelta = delta_primal;
  dualRow.createFreemove(row_ep);
  analysis->simplexTimerStop(Chuzc0Clock);
  //
  // Section 1: Pack row_ap and row_ep, then determine the possible
  // variables - candidates for CHUZC
  analysis->simplexTimerStart(Chuzc1Clock);
  // Pack row_ap into the packIndex/Value of HEkkDualRow
  dualRow.chooseMakepack(&row_ap, 0);
  // Pack row_ep into the packIndex/Value of HEkkDualRow
  dualRow.chooseMakepack(row_ep, solver_num_col);
  // Determine the possible variables - candidates for CHUZC
  dualRow.choosePossible();
  analysis->simplexTimerStop(Chuzc1Clock);
  //
  // Take action if the step to an expanded bound is not positive, or
  // there are no candidates for CHUZC
  variable_in = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    rebuild_reason = REBUILD_REASON_POSSIBLY_DUAL_UNBOUNDED;
    return;
  }
  //
  // Sections 2 and 3: Perform (bound-flipping) ratio test. This can
  // fail if the dual values are excessively large
  bool chooseColumnFail = dualRow.chooseFinal();
  if (chooseColumnFail) {
    rebuild_reason = REBUILD_REASON_CHOOSE_COLUMN_FAIL;
    return;
  }
  //
  // Section 4: Reset the nonbasicMove values for free columns
  analysis->simplexTimerStart(Chuzc4Clock);
  dualRow.deleteFreemove();
  analysis->simplexTimerStop(Chuzc4Clock);
  // Record values for basis change, checking for numerical problems and update
  // of dual variables
  variable_in = dualRow.workPivot;  // Index of the column entering the basis
  alpha_row = dualRow.workAlpha;    // Pivot value computed row-wise - used for
                                    // numerical checking
  theta_dual = dualRow.workTheta;   // Dual step length

  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX &&
      !new_devex_framework) {
    // When using Devex, unless a new framework is to be used, get the
    // exact weight for the pivotal row and, based on its accuracy,
    // determine that a new framework is to be used. In serial
    // new_devex_framework should only ever be false at this point in
    // this method, but in PAMI, this method may be called multiple
    // times in minor iterations and the new framework is set up in
    // majorUpdate.
    analysis->simplexTimerStart(DevexWtClock);
    // Determine the exact Devex weight
    dualRow.computeDevexWeight();
    computed_edge_weight = dualRow.computed_edge_weight;
    computed_edge_weight = max(1.0, computed_edge_weight);
    analysis->simplexTimerStop(DevexWtClock);
  }
  return;
}

void HEkkDual::chooseColumnSlice(HVector* row_ep) {
  // Choose the index of a column to enter the basis (CHUZC) by
  // exploiting slices of the pivotal row - for SIP and PAMI
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  analysis->simplexTimerStart(Chuzc0Clock);
  dualRow.clear();
  dualRow.workDelta = delta_primal;
  dualRow.createFreemove(row_ep);
  analysis->simplexTimerStop(Chuzc0Clock);

  //  const HighsInt solver_num_row = ekk_instance_.simplex_lp_.numRow_;
  const double local_density = 1.0 * row_ep->count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  choosePriceTechnique(simplex_info.price_strategy, local_density,
                       use_col_price, use_row_price_w_switch);

  if (analysis->analyse_simplex_data) {
    const HighsInt row_ep_count = row_ep->count;
    if (use_col_price) {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                      row_ep_count, 0.0);
      analysis->num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                      row_ep_count, analysis->row_ep_density);
      analysis->num_row_price_with_switch++;
    } else {
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                      row_ep_count, analysis->row_ep_density);
      analysis->num_row_price++;
    }
  }
  analysis->simplexTimerStart(PriceChuzc1Clock);
  // Row_ep:         PACK + CC1

  /*
  HighsInt row_ep_thread_id = 0;
  vector<HighsInt> row_ap_thread_id;
  row_ap_thread_id.resize(slice_num);
  */

#pragma omp task
  {
    dualRow.chooseMakepack(row_ep, solver_num_col);
    dualRow.choosePossible();
#ifdef OPENMP
    //    HighsInt row_ep_thread_id = omp_get_thread_num();
    //    printf("Hello world from Row_ep:         PACK + CC1 thread %"
    //    HIGHSINT_FORMAT "\n", row_ep_thread_id);
#endif
  }

  // Row_ap: PRICE + PACK + CC1
  for (HighsInt i = 0; i < slice_num; i++) {
#pragma omp task
    {
#ifdef OPENMP
      //      HighsInt row_ap_thread_id = omp_get_thread_num();
      //      printf("Hello world from omp Row_ap: PRICE + PACK + CC1 [%1"
      //      HIGHSINT_FORMAT "] thread %" HIGHSINT_FORMAT "\n", i,
      //      row_ap_thread_id);
#endif
      slice_row_ap[i].clear();

      //      slice_matrix[i].priceByRowSparseResult(slice_row_ap[i], *row_ep);

      if (use_col_price) {
        // Perform column-wise PRICE
        slice_matrix[i].priceByColumn(slice_row_ap[i], *row_ep);
      } else if (use_row_price_w_switch) {
        // Perform hyper-sparse row-wise PRICE, but switch if the density of
        // row_ap becomes extreme
        slice_matrix[i].priceByRowSparseResultWithSwitch(
            slice_row_ap[i], *row_ep, analysis->row_ap_density, 0,
            slice_matrix[i].hyperPRICE);
      } else {
        // Perform hyper-sparse row-wise PRICE
        slice_matrix[i].priceByRowSparseResult(slice_row_ap[i], *row_ep);
      }

      slice_dualRow[i].clear();
      slice_dualRow[i].workDelta = delta_primal;
      slice_dualRow[i].chooseMakepack(&slice_row_ap[i], slice_start[i]);
      slice_dualRow[i].choosePossible();
    }
  }
#pragma omp taskwait

  if (analysis->analyse_simplex_data) {
    // Determine the nonzero count of the whole row
    HighsInt row_ap_count = 0;
    for (HighsInt i = 0; i < slice_num; i++)
      row_ap_count += slice_row_ap[i].count;
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_AP,
                                   row_ap_count);
  }

  // Join CC1 results here
  for (HighsInt i = 0; i < slice_num; i++) {
    dualRow.chooseJoinpack(&slice_dualRow[i]);
  }

  analysis->simplexTimerStop(PriceChuzc1Clock);

  // Infeasible we created before
  variable_in = -1;
  if (dualRow.workTheta <= 0 || dualRow.workCount == 0) {
    rebuild_reason = REBUILD_REASON_POSSIBLY_DUAL_UNBOUNDED;
    return;
  }

  // Choose column 2, This only happens if didn't go out
  HighsInt return_code = dualRow.chooseFinal();
  if (return_code) {
    if (return_code < 0) {
      rebuild_reason = REBUILD_REASON_CHOOSE_COLUMN_FAIL;
    } else {
      rebuild_reason = REBUILD_REASON_POSSIBLY_DUAL_UNBOUNDED;
    }
    return;
  }

  analysis->simplexTimerStart(Chuzc4Clock);
  dualRow.deleteFreemove();
  analysis->simplexTimerStop(Chuzc4Clock);

  variable_in = dualRow.workPivot;
  alpha_row = dualRow.workAlpha;
  theta_dual = dualRow.workTheta;

  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX &&
      !new_devex_framework) {
    // When using Devex, unless a new framework is to be used, get the
    // exact weight for the pivotal row and, based on its accuracy,
    // determine that a new framework is to be used. In serial
    // new_devex_framework should only ever be false at this point in
    // this method, but in PAMI, this method may be called multiple
    // times in minor iterations and the new framework is set up in
    // majorUpdate.
    analysis->simplexTimerStart(DevexWtClock);
    // Determine the partial sums of the exact Devex weight
    // First the partial sum for row_ep
    dualRow.computeDevexWeight();
    // Second the partial sums for the slices of row_ap
    for (HighsInt i = 0; i < slice_num; i++)
      slice_dualRow[i].computeDevexWeight(i);
    // Accumulate the partial sums
    // Initialse with the partial sum for row_ep
    computed_edge_weight = dualRow.computed_edge_weight;
    // Update with the partial sum for row_ep
    for (HighsInt i = 0; i < slice_num; i++)
      computed_edge_weight += slice_dualRow[i].computed_edge_weight;
    computed_edge_weight = max(1.0, computed_edge_weight);
    analysis->simplexTimerStop(DevexWtClock);
  }
}

void HEkkDual::updateFtran() {
  // Compute the pivotal column (FTRAN)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  analysis->simplexTimerStart(FtranClock);
  // Clear the picotal column and indicate that its values should be packed
  col_aq.clear();
  col_aq.packFlag = true;
  // Get the constraint matrix column by combining just one column
  // with unit multiplier
  matrix->collect_aj(col_aq, variable_in, 1);
  if (analysis->analyse_simplex_data)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq,
                                    analysis->col_aq_density);
  // Perform FTRAN
  factor->ftran(col_aq, analysis->col_aq_density,
                analysis->pointer_serial_factor_clocks);
  if (analysis->analyse_simplex_data)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN, col_aq);
  const double local_col_aq_density = (double)col_aq.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_aq_density,
                                         analysis->col_aq_density);
  ekk_instance_.updateOperationResultDensity(
      local_col_aq_density, ekk_instance_.simplex_info_.col_aq_density);
  // Save the pivot value computed column-wise - used for numerical checking
  alpha_col = col_aq.array[row_out];
  analysis->simplexTimerStop(FtranClock);
}

void HEkkDual::updateFtranBFRT() {
  // Compute the RHS changes corresponding to the BFRT (FTRAN-BFRT)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  // Only time updateFtranBFRT if dualRow.workCount > 0;
  // If dualRow.workCount = 0 then dualRow.updateFlip(&col_BFRT)
  // merely clears col_BFRT so no FTRAN is performed
  bool time_updateFtranBFRT = dualRow.workCount > 0;

  if (time_updateFtranBFRT) {
    analysis->simplexTimerStart(FtranBfrtClock);
  }

  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "Before
  //  update_flip");
  dualRow.updateFlip(&col_BFRT);
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "After
  //  update_flip");

  if (col_BFRT.count) {
    if (analysis->analyse_simplex_data)
      analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN_BFRT,
                                      col_BFRT, analysis->col_BFRT_density);
    // Perform FTRAN BFRT
    factor->ftran(col_BFRT, analysis->col_BFRT_density,
                  analysis->pointer_serial_factor_clocks);
    if (analysis->analyse_simplex_data)
      analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN_BFRT,
                                     col_BFRT);
  }
  if (time_updateFtranBFRT) {
    analysis->simplexTimerStop(FtranBfrtClock);
  }
  const double local_col_BFRT_density = (double)col_BFRT.count / solver_num_row;
  analysis->updateOperationResultDensity(local_col_BFRT_density,
                                         analysis->col_BFRT_density);
  ekk_instance_.updateOperationResultDensity(
      local_col_BFRT_density, ekk_instance_.simplex_info_.col_BFRT_density);
}

void HEkkDual::updateFtranDSE(HVector* DSE_Vector) {
  // Compute the vector required to update DSE weights - being FTRAN
  // applied to the pivotal column (FTRAN-DSE)
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  analysis->simplexTimerStart(FtranDseClock);
  if (analysis->analyse_simplex_data)
    analysis->operationRecordBefore(ANALYSIS_OPERATION_TYPE_FTRAN_DSE,
                                    *DSE_Vector, analysis->row_DSE_density);
  // Perform FTRAN DSE
  factor->ftran(*DSE_Vector, analysis->row_DSE_density,
                analysis->pointer_serial_factor_clocks);
  if (analysis->analyse_simplex_data)
    analysis->operationRecordAfter(ANALYSIS_OPERATION_TYPE_FTRAN_DSE,
                                   *DSE_Vector);
  analysis->simplexTimerStop(FtranDseClock);
  const double local_row_DSE_density =
      (double)DSE_Vector->count / solver_num_row;
  analysis->updateOperationResultDensity(local_row_DSE_density,
                                         analysis->row_DSE_density);
  ekk_instance_.updateOperationResultDensity(
      local_row_DSE_density, ekk_instance_.simplex_info_.row_DSE_density);
}

void HEkkDual::updateVerify() {
  // Compare the pivot value computed row-wise and column-wise and
  // determine whether reinversion is advisable
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  // Use the two pivot values to identify numerical trouble
  if (ekk_instance_.reinvertOnNumericalTrouble(
          "HEkkDual::updateVerify", numericalTrouble, alpha_col, alpha_row,
          numerical_trouble_tolerance)) {
    rebuild_reason = REBUILD_REASON_POSSIBLY_SINGULAR_BASIS;
  }
}

void HEkkDual::updateDual() {
  // Update the dual values
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;

  // Update - dual (shift and back)
  if (theta_dual == 0) {
    // Little to do if theta_dual is zero
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase,
    //    "Before shift_cost");
    shiftCost(variable_in, -workDual[variable_in]);
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase,
    //    "After shift_cost");
  } else {
    // Update the whole vector of dual values
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase,
    //    "Before calling dualRow.updateDual");
    dualRow.updateDual(theta_dual);
    if (ekk_instance_.simplex_info_.simplex_strategy !=
            SIMPLEX_STRATEGY_DUAL_PLAIN &&
        slice_PRICE) {
      // Update the slice-by-slice copy of dual variables
      for (HighsInt i = 0; i < slice_num; i++)
        slice_dualRow[i].updateDual(theta_dual);
    }
    //    debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase,
    //    "After calling dualRow.updateDual");
  }
  // Identify the changes in the dual objective
  double dual_objective_value_change;
  const double variable_in_delta_dual = workDual[variable_in];
  const double variable_in_value = workValue[variable_in];
  const HighsInt variable_in_nonbasicFlag =
      ekk_instance_.simplex_basis_.nonbasicFlag_[variable_in];
  dual_objective_value_change =
      variable_in_nonbasicFlag * (-variable_in_value * variable_in_delta_dual);
  dual_objective_value_change *= ekk_instance_.cost_scale_;
  ekk_instance_.simplex_info_.updated_dual_objective_value +=
      dual_objective_value_change;
  // Surely variable_out_nonbasicFlag is always 0 since it's basic - so there's
  // no dual objective change
  const HighsInt variable_out_nonbasicFlag =
      ekk_instance_.simplex_basis_.nonbasicFlag_[variable_out];
  assert(variable_out_nonbasicFlag == 0);
  if (variable_out_nonbasicFlag) {
    const double variable_out_delta_dual = workDual[variable_out] - theta_dual;
    const double variable_out_value = workValue[variable_out];
    dual_objective_value_change =
        variable_out_nonbasicFlag *
        (-variable_out_value * variable_out_delta_dual);
    dual_objective_value_change *= ekk_instance_.cost_scale_;
    ekk_instance_.simplex_info_.updated_dual_objective_value +=
        dual_objective_value_change;
  }
  workDual[variable_in] = 0;
  workDual[variable_out] = -theta_dual;

  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "Before
  //  shift_back");
  shiftBack(variable_out);
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "After
  //  shift_back");
}

void HEkkDual::updatePrimal(HVector* DSE_Vector) {
  // Update the primal values and any edge weights
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    const double updated_edge_weight = dualRHS.workEdWt[row_out];
    dualRHS.workEdWt[row_out] = computed_edge_weight;
    new_devex_framework = newDevexFramework(updated_edge_weight);
  }
  // DSE_Vector is either col_DSE = B^{-1}B^{-T}e_p (if using dual
  // steepest edge weights) or row_ep = B^{-T}e_p.
  //
  // Update - primal and weight
  dualRHS.updatePrimal(&col_BFRT, 1);
  dualRHS.updateInfeasList(&col_BFRT);
  double x_out = baseValue[row_out];
  double l_out = baseLower[row_out];
  double u_out = baseUpper[row_out];
  theta_primal = (x_out - (delta_primal < 0 ? l_out : u_out)) / alpha_col;
  dualRHS.updatePrimal(&col_aq, theta_primal);
  if (dual_edge_weight_mode == DualEdgeWeightMode::STEEPEST_EDGE) {
    const double new_pivotal_edge_weight =
        dualRHS.workEdWt[row_out] / (alpha_col * alpha_col);
    const double Kai = -2 / alpha_col;
    dualRHS.updateWeightDualSteepestEdge(&col_aq, new_pivotal_edge_weight, Kai,
                                         &DSE_Vector->array[0]);
    dualRHS.workEdWt[row_out] = new_pivotal_edge_weight;
  } else if (dual_edge_weight_mode == DualEdgeWeightMode::DEVEX) {
    // Pivotal row is for the current basis: weights are required for
    // the next basis so have to divide the current (exact) weight by
    // the pivotal value
    double new_pivotal_edge_weight =
        dualRHS.workEdWt[row_out] / (alpha_col * alpha_col);
    new_pivotal_edge_weight = max(1.0, new_pivotal_edge_weight);
    // nw_wt is max(workEdWt[iRow], NewExactWeight*columnArray[iRow]^2);
    //
    // But NewExactWeight is new_pivotal_edge_weight = max(1.0,
    // dualRHS.workEdWt[row_out] / (alpha * alpha))
    //
    // so nw_wt = max(workEdWt[iRow],
    // new_pivotal_edge_weight*columnArray[iRow]^2);
    //
    // Update rest of weights
    dualRHS.updateWeightDevex(&col_aq, new_pivotal_edge_weight);
    dualRHS.workEdWt[row_out] = new_pivotal_edge_weight;
    num_devex_iterations++;
  }
  dualRHS.updateInfeasList(&col_aq);

  // Whether or not dual steepest edge weights are being used, have to
  // add in DSE_Vector->syntheticTick_ since this contains the
  // contribution from forming row_ep = B^{-T}e_p.
  ekk_instance_.total_syntheticTick_ += col_aq.syntheticTick;
  ekk_instance_.total_syntheticTick_ += DSE_Vector->syntheticTick;
}

// Record the shift in the cost of a particular column
void HEkkDual::shiftCost(const HighsInt iCol, const double amount) {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  simplex_info.costs_perturbed = true;
  if (simplex_info.workShift_[iCol] != 0) {
    printf("Column %" HIGHSINT_FORMAT " already has nonzero shift of %g\n",
           iCol, simplex_info.workShift_[iCol]);
  }
  assert(simplex_info.workShift_[iCol] == 0);
  simplex_info.workShift_[iCol] = amount;
}

// Undo the shift in the cost of a particular column
void HEkkDual::shiftBack(const HighsInt iCol) {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  simplex_info.workDual_[iCol] -= simplex_info.workShift_[iCol];
  simplex_info.workShift_[iCol] = 0;
}

void HEkkDual::updatePivots() {
  // UPDATE
  //
  // If reinversion is needed then skip this method
  if (rebuild_reason) return;
  //
  // Update the sets of indices of basic and nonbasic variables
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "Before
  //  update_pivots");
  ekk_instance_.updatePivots(variable_in, row_out, move_out);
  //  debugUpdatedObjectiveValue(ekk_instance_, algorithm, solvePhase, "After
  //  update_pivots");
  //
  ekk_instance_.iteration_count_++;
  //
  // Update the invertible representation of the basis matrix
  ekk_instance_.updateFactor(&col_aq, &row_ep, &row_out, &rebuild_reason);
  //
  // Update the row-wise representation of the nonbasic columns
  ekk_instance_.updateMatrix(variable_in, variable_out);
  //
  // Delete Freelist entry for variable_in
  dualRow.deleteFreelist(variable_in);
  //
  // Update the primal value for the row where the basis change has
  // occurred, and set the corresponding primal infeasibility value in
  // dualRHS.work_infeasibility
  dualRHS.updatePivots(
      row_out,
      ekk_instance_.simplex_info_.workValue_[variable_in] + theta_primal);

  /*
  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock = total_syntheticTick >= build_syntheticTick;
  const bool performed_min_updates =
      ekk_instance_.simplex_info_.update_count >=
      synthetic_tick_reinversion_min_update_count;
  if (reinvert_syntheticClock && performed_min_updates)
    rebuild_reason = REBUILD_REASON_SYNTHETIC_CLOCK_SAYS_INVERT;
  */
}

void HEkkDual::initialiseDevexFramework(const bool parallel) {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  // Initialise the Devex framework: reference set is all basic
  // variables
  analysis->simplexTimerStart(DevexIzClock);
  const vector<int8_t>& nonbasicFlag =
      ekk_instance_.simplex_basis_.nonbasicFlag_;
  // Initialise the devex framework. The devex reference set is
  // initialise to be the current set of basic variables - and never
  // changes until a new framework is set up. In a simplex iteration,
  // to compute the exact Devex weight for the pivotal row requires
  // summing the squares of the its entries over the indices in the
  // reference set. This is achieved by summing over all indices, but
  // multiplying the entry by the value in devex_index before
  // equaring. Thus devex_index contains 1 for indices in the
  // reference set, and 0 otherwise. This is achieved by setting the
  // values of devex_index to be 1-nonbasicFlag^2, ASSUMING
  // |nonbasicFlag|=1 iff the corresponding variable is nonbasic
  for (HighsInt vr_n = 0; vr_n < solver_num_tot; vr_n++)
    simplex_info.devex_index_[vr_n] =
        1 - nonbasicFlag[vr_n] * nonbasicFlag[vr_n];
  // Set all initial weights to 1, zero the count of iterations with
  // this Devex framework, increment the number of Devex frameworks
  // and indicate that there's no need for a new Devex framework
  dualRHS.workEdWt.assign(solver_num_row, 1.0);
  num_devex_iterations = 0;
  new_devex_framework = false;
  minor_new_devex_framework = false;
  analysis->simplexTimerStop(DevexIzClock);
}

void HEkkDual::interpretDualEdgeWeightStrategy(
    const HighsInt dual_edge_weight_strategy) {
  if (dual_edge_weight_strategy == SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_CHOOSE) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = true;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DANTZIG) {
    dual_edge_weight_mode = DualEdgeWeightMode::DANTZIG;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_DEVEX) {
    dual_edge_weight_mode = DualEdgeWeightMode::DEVEX;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = false;
  } else if (dual_edge_weight_strategy ==
             SIMPLEX_DUAL_EDGE_WEIGHT_STRATEGY_STEEPEST_EDGE_UNIT_INITIAL) {
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = false;
    allow_dual_steepest_edge_to_devex_switch = false;
  } else {
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                "HEkkDual::interpretDualEdgeWeightStrategy: "
                "unrecognised dual_edge_weight_strategy = %" HIGHSINT_FORMAT
                " - using "
                "dual steepest edge with possible switch to Devex\n",
                dual_edge_weight_strategy);
    dual_edge_weight_mode = DualEdgeWeightMode::STEEPEST_EDGE;
    initialise_dual_steepest_edge_weights = true;
    allow_dual_steepest_edge_to_devex_switch = true;
  }
}

void HEkkDual::saveDualRay() {
  ekk_instance_.simplex_lp_status_.has_dual_ray = true;
  ekk_instance_.simplex_info_.dual_ray_row_ = row_out;
  ekk_instance_.simplex_info_.dual_ray_sign_ = move_out;
}

void HEkkDual::assessPhase1Optimality() {
  // Should only be called when optimal in phase 1 (row_out == -1) with negative
  // dual activity
  assert(solvePhase == SOLVE_PHASE_1);
  assert(row_out == -1);
  //  assert(ekk_instance_.simplex_info_.dual_objective_value < 0);
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsModelStatus& scaled_model_status = ekk_instance_.scaled_model_status_;
  // We still have dual infeasibilities, so clean up any perturbations
  // before concluding dual infeasibility
  //
  // What if the dual objective is positive but tiny?
  if (fabs(simplex_info.dual_objective_value) <= primal_feasibility_tolerance)
    highsLogUser(ekk_instance_.options_.log_options, HighsLogType::INFO,
                 "Optimal in phase 1 but not jumping to phase 2 since "
                 "dual objective is %10.4g: Costs perturbed = %" HIGHSINT_FORMAT
                 "\n",
                 simplex_info.dual_objective_value,
                 ekk_instance_.simplex_info_.costs_perturbed);
  if (ekk_instance_.simplex_info_.costs_perturbed) {
    // Clean up perturbation
    cleanup();
    // If there are now dual infeasibilities with respect to phase 1
    // bounds, have to got back to rebuild()
    if (dualInfeasCount == 0) {
      // No dual infeasibilities with respect to phase 1 bounds.
      if (simplex_info.dual_objective_value == 0) {
        // No dual infeasibilities (with respect to phase 2 bounds) so
        // go to phase 2
        highsLogUser(ekk_instance_.options_.log_options, HighsLogType::INFO,
                     "LP is dual feasible after removing cost perturbations "
                     "so go to phase 2\n");
      } else {
        // LP is dual infeasible if the dual objective is sufficiently
        // positive, so no conclusions on the primal LP can be deduced
        // - could be primal unbounded or primal infeasible.
        //
        // Shift any dual infeasibilities and go to dual phase 2. If a
        // primal feasible point is found then the shifts are removed
        // and primal phase 2 will identify whether the LP is primal
        // unbounded. If dual unboundedness is found, then no
        // conclusion can be drawn. Have to use primal phase 1 (and
        // possibly phase 2) to determine whether the LP is primal
        // infeasible or unbounded.
        //
        // What's important is that the solver doesn't go back to dual
        // phase 1, otherwise it can fail to terminate
        //
        // Indicate the conclusion of dual infeasiblility by setting
        // the scaled model status
        reportOnPossibleLpDualInfeasibility();
        scaled_model_status = HighsModelStatus::DUAL_INFEASIBLE;
      }
      solvePhase = SOLVE_PHASE_2;
    }
  } else {
    // Phase 1 problem is optimal with original costs and negative
    // dual objective. In this case, hsol deduces dual infeasibility
    // and returns UNBOUNDED as a status, but this is wrong if the LP
    // is primal infeasible. As discussed above, this can only be
    // determined by going to dual phase 2, and then primal phase 1,
    // if necessary.
    //
    reportOnPossibleLpDualInfeasibility();
    scaled_model_status = HighsModelStatus::DUAL_INFEASIBLE;
    solvePhase = SOLVE_PHASE_2;
  }
  if (dualInfeasCount > 0) {
    // Must still be solvePhase = SOLVE_PHASE_1 since dual infeasibilities with
    // respect to phase 1 bounds mean that primal values must
    // change, so primal feasibility is unknown
    assert(solvePhase == SOLVE_PHASE_1);
  } else {
    // Optimal in dual phase 1, so must go to phase 2
    assert(solvePhase == SOLVE_PHASE_2);
    // Reset the duals, if necessary shifting costs of free variable
    // so that their duals are zero
    exitPhase1ResetDuals();
  }
}

void HEkkDual::exitPhase1ResetDuals() {
  const HighsLp& simplex_lp = ekk_instance_.simplex_lp_;
  const SimplexBasis& simplex_basis = ekk_instance_.simplex_basis_;
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;

  const bool reperturb_costs = true;
  if (reperturb_costs) {
    if (simplex_info.costs_perturbed) {
      highsLogDev(ekk_instance_.options_.log_options, HighsLogType::INFO,
                  "Costs are already perturbed in exitPhase1ResetDuals\n");
    } else {
      highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
                  "Re-perturbing costs when optimal in phase 1\n");
      ekk_instance_.initialiseCost(SimplexAlgorithm::DUAL, SOLVE_PHASE_2, true);
      ekk_instance_.computeDual();
    }
  }

  const HighsInt numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  HighsInt num_shift = 0;
  double sum_shift = 0;
  for (HighsInt iVar = 0; iVar < numTot; iVar++) {
    if (simplex_basis.nonbasicFlag_[iVar]) {
      double lp_lower;
      double lp_upper;
      if (iVar < simplex_lp.numCol_) {
        lp_lower = simplex_lp.colLower_[iVar];
        lp_upper = simplex_lp.colUpper_[iVar];
      } else {
        HighsInt iRow = iVar - simplex_lp.numCol_;
        lp_lower = simplex_lp.rowLower_[iRow];
        lp_upper = simplex_lp.rowUpper_[iRow];
      }
      if (lp_lower <= -kHighsInf && lp_upper >= kHighsInf) {
        const double shift = -simplex_info.workDual_[iVar];
        simplex_info.workDual_[iVar] = 0;
        simplex_info.workCost_[iVar] = simplex_info.workCost_[iVar] + shift;
        num_shift++;
        sum_shift += fabs(shift);
        highsLogDev(ekk_instance_.options_.log_options, HighsLogType::VERBOSE,
                    "Variable %" HIGHSINT_FORMAT
                    " is free: shift cost to zero dual of %g\n",
                    iVar, shift);
      }
    }
  }
  if (num_shift)
    highsLogDev(ekk_instance_.options_.log_options, HighsLogType::DETAILED,
                "Performed %" HIGHSINT_FORMAT
                " cost shift(s) for free variables to zero "
                "dual values: total = %g\n",
                num_shift, sum_shift);
}

void HEkkDual::reportOnPossibleLpDualInfeasibility() {
  HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HighsSimplexAnalysis& analysis = ekk_instance_.analysis_;
  assert(solvePhase == SOLVE_PHASE_1);
  assert(row_out == -1);
  //  assert(simplex_info.dual_objective_value < 0);
  assert(!simplex_info.costs_perturbed);
  std::string lp_dual_status;
  if (analysis.num_dual_phase_1_lp_dual_infeasibility) {
    lp_dual_status = "infeasible";
  } else {
    lp_dual_status = "feasible";
  }
  highsLogUser(ekk_instance_.options_.log_options, HighsLogType::INFO,
               "LP is dual %s with dual phase 1 objective %10.4g and num / "
               "max / sum dual infeasibilities = %" HIGHSINT_FORMAT
               " / %9.4g / %9.4g\n",
               lp_dual_status.c_str(), simplex_info.dual_objective_value,
               analysis.num_dual_phase_1_lp_dual_infeasibility,
               analysis.max_dual_phase_1_lp_dual_infeasibility,
               analysis.sum_dual_phase_1_lp_dual_infeasibility);
}

bool HEkkDual::dualInfoOk(const HighsLp& lp) {
  HighsInt lp_numCol = lp.numCol_;
  HighsInt lp_numRow = lp.numRow_;
  bool dimensions_ok;
  dimensions_ok = lp_numCol == solver_num_col && lp_numRow == solver_num_row;
  assert(dimensions_ok);
  if (!dimensions_ok) {
    printf("LP-Solver dimension incompatibility (%" HIGHSINT_FORMAT
           ", %" HIGHSINT_FORMAT ") != (%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
           ")\n",
           lp_numCol, solver_num_col, lp_numRow, solver_num_row);
    return false;
  }
  dimensions_ok = lp_numCol == factor->numCol && lp_numRow == factor->numRow;
  assert(dimensions_ok);
  if (!dimensions_ok) {
    printf("LP-Factor dimension incompatibility (%" HIGHSINT_FORMAT
           ", %" HIGHSINT_FORMAT ") != (%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
           ")\n",
           lp_numCol, factor->numCol, lp_numRow, factor->numRow);
    return false;
  }
  return true;
}

bool HEkkDual::bailoutReturn() {
  if (solve_bailout) {
    // If bailout has already been decided: check that it's for one of
    // these reasons
    assert(ekk_instance_.scaled_model_status_ ==
               HighsModelStatus::REACHED_TIME_LIMIT ||
           ekk_instance_.scaled_model_status_ ==
               HighsModelStatus::REACHED_ITERATION_LIMIT ||
           ekk_instance_.scaled_model_status_ ==
               HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);
  }
  return solve_bailout;
}

bool HEkkDual::bailoutOnTimeIterations() {
  HighsModelStatus& scaled_model_status = ekk_instance_.scaled_model_status_;
  if (solve_bailout) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(scaled_model_status == HighsModelStatus::REACHED_TIME_LIMIT ||
           scaled_model_status == HighsModelStatus::REACHED_ITERATION_LIMIT ||
           scaled_model_status ==
               HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);
  } else if (ekk_instance_.timer_.readRunHighsClock() >
             ekk_instance_.options_.time_limit) {
    solve_bailout = true;
    scaled_model_status = HighsModelStatus::REACHED_TIME_LIMIT;
  } else if (ekk_instance_.iteration_count_ >=
             ekk_instance_.options_.simplex_iteration_limit) {
    solve_bailout = true;
    scaled_model_status = HighsModelStatus::REACHED_ITERATION_LIMIT;
  }
  return solve_bailout;
}

bool HEkkDual::bailoutOnDualObjective() {
  if (solve_bailout) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(ekk_instance_.scaled_model_status_ ==
               HighsModelStatus::REACHED_TIME_LIMIT ||
           ekk_instance_.scaled_model_status_ ==
               HighsModelStatus::REACHED_ITERATION_LIMIT ||
           ekk_instance_.scaled_model_status_ ==
               HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);
  } else if (ekk_instance_.simplex_lp_.sense_ == ObjSense::MINIMIZE &&
             solvePhase == SOLVE_PHASE_2) {
    if (ekk_instance_.simplex_info_.updated_dual_objective_value >
        ekk_instance_.options_.dual_objective_value_upper_bound)
      solve_bailout = reachedExactDualObjectiveValueUpperBound();
  }
  return solve_bailout;
}

bool HEkkDual::reachedExactDualObjectiveValueUpperBound() {
  // Solving a minimization in dual simplex phase 2, and dual
  // objective exceeds the prescribed upper bound. However, costs
  // will be perturbed, so need to check whether exact dual
  // objective value exceeds the prescribed upper bound. This can be
  // a relatively expensive calculation, so determine whether to do
  // it according to the sparsity of the pivotal row
  bool reached_exact_dual_objective_value_upper_bound = false;
  double use_row_ap_density =
      std::min(std::max(analysis->row_ap_density, 0.01), 1.0);
  HighsInt check_frequency = 1.0 / use_row_ap_density;
  assert(check_frequency > 0);

  bool check_exact_dual_objective_value =
      ekk_instance_.simplex_info_.update_count % check_frequency == 0;

  if (check_exact_dual_objective_value) {
    const double dual_objective_value_upper_bound =
        ekk_instance_.options_.dual_objective_value_upper_bound;
    const double perturbed_dual_objective_value =
        ekk_instance_.simplex_info_.updated_dual_objective_value;
    const double perturbed_value_residual =
        perturbed_dual_objective_value - dual_objective_value_upper_bound;
    const double exact_dual_objective_value = computeExactDualObjectiveValue();
    const double exact_value_residual =
        exact_dual_objective_value - dual_objective_value_upper_bound;
    std::string action;
    if (exact_dual_objective_value > dual_objective_value_upper_bound) {
#ifdef SCIP_DEV
      printf("HEkkDual::solvePhase2: %12g = Objective > ObjectiveUB\n",
             ekk_instance_.simplex_info_.updated_dual_objective_value,
             dual_objective_value_upper_bound);
#endif
      action = "Have DualUB bailout";
      reached_exact_dual_objective_value_upper_bound = true;
      ekk_instance_.scaled_model_status_ =
          HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND;
    } else {
      action = "No   DualUB bailout";
    }
    highsLogUser(ekk_instance_.options_.log_options, HighsLogType::INFO,
                 "%s on iteration %" HIGHSINT_FORMAT
                 ": Density %11.4g; Frequency %" HIGHSINT_FORMAT
                 ": "
                 "Residual(Perturbed = %g; Exact = %g)\n",
                 action.c_str(), ekk_instance_.iteration_count_,
                 use_row_ap_density, check_frequency, perturbed_value_residual,
                 exact_value_residual);
  }
  return reached_exact_dual_objective_value_upper_bound;
}

double HEkkDual::computeExactDualObjectiveValue() {
  const HighsLp& simplex_lp = ekk_instance_.simplex_lp_;
  const SimplexBasis& simplex_basis = ekk_instance_.simplex_basis_;
  const HighsSimplexInfo& simplex_info = ekk_instance_.simplex_info_;
  HMatrix& matrix = ekk_instance_.matrix_;
  HFactor& factor = ekk_instance_.factor_;
  // Create a local buffer for the pi vector
  HVector dual_col;
  dual_col.setup(simplex_lp.numRow_);
  dual_col.clear();
  for (HighsInt iRow = 0; iRow < simplex_lp.numRow_; iRow++) {
    HighsInt iVar = simplex_basis.basicIndex_[iRow];
    if (iVar < simplex_lp.numCol_) {
      const double value = simplex_lp.colCost_[iVar];
      if (value) {
        dual_col.array[iRow] = value;
        dual_col.index[dual_col.count++] = iRow;
      }
    }
  }
  // Create a local buffer for the dual vector
  const HighsInt numTot = simplex_lp.numCol_ + simplex_lp.numRow_;
  HVector dual_row;
  dual_row.setup(simplex_lp.numCol_);
  dual_row.clear();
  if (dual_col.count) {
    const double historical_density_for_non_hypersparse_operation = 1;
    factor.btran(dual_col, historical_density_for_non_hypersparse_operation);
    matrix.priceByColumn(dual_row, dual_col);
  }
  double dual_objective = simplex_lp.offset_;
  double norm_dual = 0;
  double norm_delta_dual = 0;
  for (HighsInt iCol = 0; iCol < simplex_lp.numCol_; iCol++) {
    if (!simplex_basis.nonbasicFlag_[iCol]) continue;
    double exact_dual = simplex_lp.colCost_[iCol] - dual_row.array[iCol];
    double residual = fabs(exact_dual - simplex_info.workDual_[iCol]);
    norm_dual += fabs(exact_dual);
    norm_delta_dual += residual;
    if (residual > 1e10)
      highsLogUser(
          ekk_instance_.options_.log_options, HighsLogType::WARNING,
          "Col %4" HIGHSINT_FORMAT
          ": ExactDual = %11.4g; WorkDual = %11.4g; Residual = %11.4g\n",
          iCol, exact_dual, simplex_info.workDual_[iCol], residual);
    dual_objective += simplex_info.workValue_[iCol] * exact_dual;
  }
  for (HighsInt iVar = simplex_lp.numCol_; iVar < numTot; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    HighsInt iRow = iVar - simplex_lp.numCol_;
    double exact_dual = -dual_col.array[iRow];
    double residual = fabs(exact_dual - simplex_info.workDual_[iVar]);
    norm_dual += fabs(exact_dual);
    norm_delta_dual += residual;
    if (residual > 1e10)
      highsLogUser(
          ekk_instance_.options_.log_options, HighsLogType::WARNING,
          "Row %4" HIGHSINT_FORMAT
          ": ExactDual = %11.4g; WorkDual = %11.4g; Residual = %11.4g\n",
          iRow, exact_dual, simplex_info.workDual_[iVar], residual);
    dual_objective += simplex_info.workValue_[iVar] * exact_dual;
  }
  double relative_delta = norm_delta_dual / std::max(norm_dual, 1.0);
  if (relative_delta > 1e-3)
    highsLogUser(
        ekk_instance_.options_.log_options, HighsLogType::WARNING,
        "||exact dual vector|| = %g; ||delta dual vector|| = %g: ratio = %g\n",
        norm_dual, norm_delta_dual, relative_delta);
  return dual_objective;
}

HighsDebugStatus HEkkDual::debugDualSimplex(const std::string message,
                                            const bool initialise) {
  HighsDebugStatus return_status = ekkDebugSimplex(
      message, ekk_instance_, algorithm, solvePhase, initialise);
  if (return_status == HighsDebugStatus::kLogicalError) return return_status;
  if (initialise) return return_status;
  return HighsDebugStatus::kOk;
}
