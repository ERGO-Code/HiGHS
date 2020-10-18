/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HEkk.h"

#include <cassert>
#include <iostream>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
//#include "simplex/HSimplex.h"
//#include "simplex/HSimplexDebug.h"
#include "simplex/HFactorDebug.h"
//#include "simplex/SimplexTimer.h"
//#include "util/HighsRandom.h"
//#include "util/HighsUtils.h"

using std::cout;
using std::endl;

HighsStatus HEkk::init() {
  solver_num_col = lp_.numCol_;
  solver_num_row = lp_.numRow_;
  solver_num_tot = solver_num_col + solver_num_row;
  return HighsStatus::OK;
}
HighsStatus HEkk::solve() {
  init();
  HighsLogMessage(options_.logfile, HighsMessageType::INFO,
		  "HEkk::solve called for LP with %d columns, %d rows and %d entries",
		  solver_num_col, solver_num_row, lp_.Astart_[solver_num_col]);
  bool positive_num_row = solver_num_row > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    HighsLogMessage(options_.logfile, HighsMessageType::ERROR,
                    "HPrimal::solve called for LP with non-positive (%d) "
                    "number of constraints",
                    solver_num_row);
    return HighsStatus::Error;
  }
  invertHint = INVERT_HINT_NO;

  // First setup the factor arrays if they don't exist
  if (!simplex_lp_status.has_factor_arrays) {
    factor.setup(solver_num_col, solver_num_row, &lp_.Astart_[0],
                 &lp_.Aindex_[0], &lp_.Avalue_[0],
                 &simplex_basis_.basicIndex_[0], options_.highs_debug_level,
                 options_.logfile, options_.output, options_.message_level,
                 options_.factor_pivot_threshold,
                 options_.factor_pivot_tolerance);
    simplex_lp_status.has_factor_arrays = true;
  }
  // Reinvert if there isn't a fresh INVERT. ToDo Override this for MIP hot
  // start
  //  bool reinvert = !simplex_lp_status.has_fresh_invert;
  bool reinvert = !simplex_lp_status.has_invert;
  if (reinvert) {
    //    analysis.simplexTimerStart(InvertClock);
    const int rank_deficiency = computeFactor();
    //    analysis.simplexTimerStop(InvertClock);
    if (rank_deficiency) return HighsStatus::Error;
  }
  assert(simplex_lp_status.has_fresh_invert);

  // Possibly set up the HMatrix column-wise and row-wise copies of the matrix
  if (!simplex_lp_status.has_matrix_col_wise ||
      !simplex_lp_status.has_matrix_row_wise) {
    //    analysis.simplexTimerStart(matrixSetupClock);
    matrix.setup(solver_num_col, solver_num_row, &lp_.Astart_[0],
                 &lp_.Aindex_[0], &lp_.Avalue_[0],
                 &simplex_basis_.nonbasicFlag_[0]);
    simplex_lp_status.has_matrix_col_wise = true;
    simplex_lp_status.has_matrix_row_wise = true;
    //    analysis.simplexTimerStop(matrixSetupClock);
  }

  // Set up the simplex work and base arrays
  //  analysis.simplexTimerStart(allocateSimplexArraysClock);
  allocateWorkAndBaseArrays();
  //  analysis.simplexTimerStop(allocateSimplexArraysClock);

  //  analysis.simplexTimerStart(initialiseSimplexCostBoundsClock);
  initialiseCost();
  initialiseBound();
  //  analysis.simplexTimerStop(initialiseSimplexCostBoundsClock);

  // Possibly solve for the basic primal and nonbasic dual values to determine
  // which simplex solver to use, unless it's forced
  //  if (simplex_lp_status.has_basic_primal_values) {
  initialiseNonbasicWorkValue();
  //  analysis.simplexTimerStart(ComputePrimalClock);
  computePrimal();
  //  analysis.simplexTimerStop(ComputePrimalClock);
  simplex_lp_status.has_basic_primal_values = true;
  //}
  //  if (simplex_lp_status.has_basic_dual_values) {
  //  analysis.simplexTimerStart(ComputeDualClock);
  computeDual();
  //  analysis.simplexTimerStop(ComputeDualClock);
  simplex_lp_status.has_nonbasic_dual_values = true;
  //}

  solvePhase = 0;  // Frig to skip while (solvePhase) {*}
  assert(simplex_lp_status.has_primal_objective_value);
  updated_primal_objective_value = primal_objective_value;

  return HighsStatus::Error;
}

// Private methods
int HEkk::computeFactor() {
  const int rank_deficiency = factor.build();
  const bool force = rank_deficiency;
  debugCheckInvert(options_, factor, force);
  if (rank_deficiency) {
    // Have an invertible representation, but of B with column(s)
    // replacements due to singularity. So no (fresh) representation of
    // B^{-1}
    simplex_lp_status.has_invert = false;
    simplex_lp_status.has_fresh_invert = false;
  } else {
    // Now have a representation of B^{-1}, and it is fresh!
    simplex_lp_status.has_invert = true;
    simplex_lp_status.has_fresh_invert = true;
  }
  // Set the update count to zero since the corrected invertible
  // representation may be used for an initial basis. In any case the
  // number of updates shouldn't be positive
  update_count = 0;
  return rank_deficiency;
}

void HEkk::computePrimal() {
  //  HighsSimplexAnalysis* analysis = &highs_model_object.simplex_analysis_;
  // Setup a local buffer for the values of basic variables
  HVector primal_col;
  primal_col.setup(solver_num_row);
  primal_col.clear();
  for (int i = 0; i < solver_num_tot; i++) {
    if (simplex_basis_.nonbasicFlag_[i] && workValue_[i] != 0)
      matrix.collect_aj(primal_col, i, workValue_[i]);
  }
  // If debugging, take a copy of the RHS
  vector<double> debug_primal_rhs;
  if (options_.highs_debug_level >= HIGHS_DEBUG_LEVEL_COSTLY)
    debug_primal_rhs = primal_col.array;

  // It's possible that the buffer has no nonzeros, so performing
  // FTRAN is unnecessary. Not much of a saving, but the zero density
  // looks odd in the analysis!
  if (primal_col.count) {
    factor.ftran(primal_col, analysis->primal_col_density);
    const double local_primal_col_density =
        (double)primal_col.count / solver_num_row;
    analysis->updateOperationResultDensity(local_primal_col_density,
                                           analysis->primal_col_density);
  }
  for (int i = 0; i < solver_num_row; i++) {
    int iCol = simplex_basis_.basicIndex_[i];
    baseValue_[i] = -primal_col.array[i];
    baseLower_[i] = workLower_[iCol];
    baseUpper_[i] = workUpper_[iCol];
  }
  //  debugComputePrimal(highs_model_object, debug_primal_rhs);
  // Now have basic primals
  simplex_lp_status.has_basic_primal_values = true;
}

void HEkk::computeDual() {
  // Create a local buffer for the pi vector
  HVector dual_col;
  dual_col.setup(solver_num_row);
  dual_col.clear();
  for (int iRow = 0; iRow < solver_num_row; iRow++) {
    const double value =
        workCost_[simplex_basis_.basicIndex_[iRow]] +
        workShift_[simplex_basis_.basicIndex_[iRow]];
    if (value) {
      dual_col.count++;
      dual_col.index[iRow] = iRow;
      dual_col.array[iRow] = value;
    }
  }
  // If debugging, take a copy of the basic costs and any previous duals
  vector<double> debug_previous_workDual;
  vector<double> debug_basic_costs;
  if (options_.highs_debug_level >=
      HIGHS_DEBUG_LEVEL_COSTLY) {
    debug_basic_costs = dual_col.array;
    if (simplex_lp_status.has_nonbasic_dual_values)
      debug_previous_workDual = workDual_;
  }
  // Copy the costs in case the basic costs are all zero
  for (int i = 0; i < solver_num_tot; i++)
    workDual_[i] = workCost_[i];
  if (dual_col.count) {
    // RHS of row dual calculation is nonzero
#ifdef HiGHSDEV
    if (analyse_iterations)
      analysis.operationRecordBefore(ANALYSIS_OPERATION_TYPE_BTRAN_FULL,
                                     dual_col, analysis.dual_col_density);
#endif
    factor.btran(dual_col, analysis.dual_col_density,
                 analysis.pointer_serial_factor_clocks);
#ifdef HiGHSDEV
    if (analyse_iterations)
      analysis.operationRecordAfter(ANALYSIS_OPERATION_TYPE_BTRAN_FULL,
                                    dual_col);
#endif
    const double local_dual_col_density =
        (double)dual_col.count / solver_num_row;
    analysis.updateOperationResultDensity(local_dual_col_density,
                                          analysis.dual_col_density);
    // Create a local buffer for the values of reduced costs
    HVector dual_row;
    dual_row.setup(solver_num_col);
    dual_row.clear();
#ifdef HiGHSDEV
    double price_full_historical_density = 1;
    if (analyse_iterations)
      analysis.operationRecordBefore(ANALYSIS_OPERATION_TYPE_PRICE_FULL,
                                     dual_row, price_full_historical_density);
#endif
    matrix.priceByColumn(dual_row, dual_col);
#ifdef HiGHSDEV
    if (analyse_iterations)
      analysis.operationRecordAfter(ANALYSIS_OPERATION_TYPE_PRICE_FULL,
                                    dual_row);
#endif
    for (int i = 0; i < solver_num_col; i++)
      workDual_[i] -= dual_row.array[i];
    for (int i = solver_num_col; i < solver_num_tot; i++)
      workDual_[i] -= dual_col.array[i - solver_num_col];
    // Possibly analyse the computed dual values
    debugComputeDual(highs_model_object, debug_previous_workDual,
                     debug_basic_costs, dual_col.array);
  }
  // Now have nonbasic duals
  simplex_lp_status.has_nonbasic_dual_values = true;
}

void HEkk::allocateWorkAndBaseArrays() {
  // Allocate bounds and solution spaces
  workCost_.resize(solver_num_tot);
  workDual_.resize(solver_num_tot);
  workShift_.resize(solver_num_tot);

  workLower_.resize(solver_num_tot);
  workUpper_.resize(solver_num_tot);
  workRange_.resize(solver_num_tot);
  workValue_.resize(solver_num_tot);

  baseLower_.resize(solver_num_row);
  baseUpper_.resize(solver_num_row);
  baseValue_.resize(solver_num_row);
}

void HEkk::initialisePhase2ColCost() {
  // Copy the Phase 2 cost and zero the shift
  for (int col = 0; col < solver_num_col; col++) {
    int var = col;
    workCost_[var] = (int)lp_.sense_ * lp_.colCost_[col];
    workShift_[var] = 0;
  }
}

void HEkk::initialisePhase2RowCost() {
  // Zero the cost and shift
  for (int iVar = solver_num_col; iVar < solver_num_col + solver_num_row; iVar++) {
    workCost_[iVar] = 0;
    workShift_[iVar] = 0;
  }
}

void HEkk::initialiseCost(int perturb) {
  // Copy the cost
  initialisePhase2ColCost();
  initialisePhase2RowCost();
  // See if we want to skip perturbation
  costs_perturbed = 0;
  if (perturb == 0 || dual_simplex_cost_perturbation_multiplier == 0) return;
  costs_perturbed = 1;

  // Perturb the original costs, scale down if is too big
  int num_original_nonzero_cost = 0;
  double bigc = 0;
  for (int i = 0; i < solver_num_col; i++) {
    const double abs_cost = fabs(workCost_[i]);
    bigc = max(bigc, abs_cost);
    if (abs_cost) num_original_nonzero_cost++;
  }
  const int pct0 = (100 * num_original_nonzero_cost) / solver_num_col;
  double average_cost = 0;
  if (num_original_nonzero_cost) average_cost = bigc / num_original_nonzero_cost;
  HighsPrintMessage(options_.output, options_.message_level, ML_DETAILED,
		    "Cost perturbation for %s: Have %d nonzero costs (%3d%%) with bigc = %g and average = %g",
		    lp_.model_name_.c_str(), num_original_nonzero_cost, pct0, bigc, average_cost);
  if (bigc > 100) {
    bigc = sqrt(sqrt(bigc));
    HighsPrintMessage(options_.output, options_.message_level, ML_DETAILED,
		      " is large so set bigc = sqrt(bigc) = %g\n", bigc);
  } else {
    HighsPrintMessage(options_.output, options_.message_level, ML_DETAILED, "\n");
  }

  // If there's few boxed variables, we will just use simple perturbation
  double boxedRate = 0;
  const int solver_num_tot = solver_num_col + solver_num_row;
  for (int i = 0; i < solver_num_tot; i++)
    boxedRate += (workRange_[i] < 1e30);
  boxedRate /= solver_num_tot;
  if (boxedRate < 0.01) {
    bigc = min(bigc, 1.0);
    HighsPrintMessage(options_.output, options_.message_level, ML_DETAILED,
		      "small boxedRate (%g) so set bigc = min(bigc, 1.0) = %g\n",
		      boxedRate, bigc);
  
  }
  // Determine the perturbation base
  double base = 5e-7 * bigc;
  // Now do the perturbation
  //  const vector<double> previous_workCost = workCost;
  for (int i = 0; i < solver_num_col; i++) {
    double lower = lp_.colLower_[i];
    double upper = lp_.colUpper_[i];
    double xpert = (fabs(workCost_[i]) + 1) * base *
                   dual_simplex_cost_perturbation_multiplier *
                   (1 + numTotRandomValue_[i]);
    if (lower <= -HIGHS_CONST_INF && upper >= HIGHS_CONST_INF) {
      // Free - no perturb
    } else if (upper >= HIGHS_CONST_INF) {  // Lower
      workCost_[i] += xpert;
    } else if (lower <= -HIGHS_CONST_INF) {  // Upper
      workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      workCost_[i] +=
          (workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
  }
    // Pass workCost_ and previous_workCost; to analysis code
    //    perturbation1[i] = fabs(workCost_[i] - previous_cost[i]);
    //    if (perturbation1) updateValueDistribution(perturbation1, analysis->cost_perturbation1_distribution);

  vector<double> perturbation2(solver_num_col);
  for (int i = solver_num_col; i < solver_num_tot; i++) {
    perturbation2[i] =
        (0.5 - numTotRandomValue_[i]) *
        dual_simplex_cost_perturbation_multiplier * 1e-12;
    workCost_[i] += perturbation2[i];
    
  }
    // Pass perturbation2 to analysis code
    //    perturbation2 = fabs(perturbation2);
    //    updateValueDistribution(perturbation2, analysis->cost_perturbation2_distribution);
}

void HEkk::initialisePhase2ColBound() {
  // Copy bounds and compute ranges
  for (int iCol = 0; iCol < solver_num_col; iCol++) {
    workLower_[iCol] = lp_.colLower_[iCol];
    workUpper_[iCol] = lp_.colUpper_[iCol];
    workRange_[iCol] = workUpper_[iCol] - workLower_[iCol];
  }
}

void HEkk::initialisePhase2RowBound() {
  // Copy bounds and compute ranges
  for (int row = 0; row < solver_num_row; row++) {
    int var = solver_num_col + row;
    workLower_[var] = -lp_.rowUpper_[row];
    workUpper_[var] = -lp_.rowLower_[row];
    workRange_[var] = workUpper_[var] - workLower_[var];
  }
}

void HEkk::initialiseBound(int phase) {
  // Initialise the Phase 2 bounds (and ranges). NB Phase 2 bounds
  // necessary to compute Phase 1 bounds
  initialisePhase2ColBound();
  initialisePhase2RowBound();
  if (phase == 2) return;

  // The dual objective is the sum of products of primal and dual
  // values for nonbasic variables. For dual simplex phase 1, the
  // primal bounds are set so that when the dual value is feasible, the
  // primal value is set to zero. Otherwise the value is +1/-1
  // according to the required sign of the dual, except for free
  // variables, where the bounds are [-1000, 1000]. Hence the dual
  // objective is the negation of the sum of infeasibilities, unless there are
  // free In Phase 1: change to dual phase 1 bound.
  const double inf = HIGHS_CONST_INF;
  for (int i = 0; i < solver_num_tot; i++) {
    if (workLower_[i] == -inf &&
        workUpper_[i] == inf) {
      // Don't change for row variables: they should never become
      // nonbasic when starting from a logical basis, and no crash
      // should make a free row nonbasic, but could an advanced basis
      // make a free row nonbasic.
      // But what it it happened?
      if (i >= solver_num_col) continue;
      workLower_[i] = -1000,
      workUpper_[i] = 1000;  // FREE
    } else if (workLower_[i] == -inf) {
      workLower_[i] = -1, workUpper_[i] = 0;  // UPPER
    } else if (workUpper_[i] == inf) {
      workLower_[i] = 0, workUpper_[i] = 1;  // LOWER
    } else {
      workLower_[i] = 0,
      workUpper_[i] = 0;  // BOXED or FIXED
    }
    workRange_[i] = workUpper_[i] - workLower_[i];
  }
}

void HEkk::initialiseNonbasicWorkValue() {
  for (int iVar = 0; iVar < solver_num_tot; iVar++) {
    if (!simplex_basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic variable
    const double lower = workLower_[iVar];
    const double upper = workUpper_[iVar];
    double value;
    if (lower == upper) {
      value = lower;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_UP) {
      value = lower;
    } else if (simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_DN) {
      value = upper;
    } else {
      assert(simplex_basis_.nonbasicMove_[iVar] == NONBASIC_MOVE_ZE);
      value = 0;
    }
    workValue_[iVar] = value;
  }
}
