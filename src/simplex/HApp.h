/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef SIMPLEX_HAPP_H_
#define SIMPLEX_HAPP_H_

// todo: clear includes.
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
#include "lp_data/HighsStatus.h"
#include "simplex/HEkk.h"
#include "simplex/HEkkPrimal.h"
#include "simplex/HSimplex.h"
#include "simplex/HSimplexNlaDebug.h"
#include "simplex/HSimplexReport.h"
#include "simplex/SimplexConst.h"
#include "simplex/SimplexTimer.h"
#include "util/HighsUtils.h"

// Single method to solve an LP with the simplex method. Solves the
// scaled LP then analyses the unscaled solution. If it doesn't satisfy
// the required tolerances, tolerances for the scaled LP are
// identified which, if used, might yield an unscaled solution that
// satisfies the required tolerances.
//
// If possible, it sets the HiGHS basis and solution

HighsStatus returnFromSolveLpSimplex(HighsLpSolverObject& solver_object,
                                     HighsStatus return_status) {
  HighsOptions& options = solver_object.options_;
  HEkk& ekk_instance = solver_object.ekk_instance_;
  HighsLp& incumbent_lp = solver_object.lp_;
  HSimplexNla& simplex_nla = ekk_instance.simplex_nla_;
  // Copy the simplex iteration count to highs_info_ from ekk_instance
  solver_object.highs_info_.simplex_iteration_count =
      ekk_instance.iteration_count_;
  // Ensure that the incumbent LP is neither moved, nor scaled
  assert(!incumbent_lp.is_moved_);
  assert(!incumbent_lp.is_scaled_);
  // Cannot expect any more with an error return
  if (return_status == HighsStatus::kError) return return_status;
  //
  // Ensure that there is an invert for the current LP
  assert(ekk_instance.status_.has_invert);
  // Ensure that simplex NLA is set up and has the right scaling
  assert(simplex_nla.is_setup_);
  // Set the simplex NLA scaling
  simplex_nla.setLpAndScalePointers(&incumbent_lp);
  if (incumbent_lp.scale_.has_scaling) {
    // The LP has scaling, so ensure that the simplex NLA has its scaling
    void* nla_scale = (void*)simplex_nla.scale_;
    void* lp_scale = (void*)(&incumbent_lp.scale_);
    assert(nla_scale == lp_scale);
  } else {
    // The LP has no scaling, so ensure that the simplex NLA has no scaling
    assert(!simplex_nla.scale_);
  }
  // ToDo Need to switch off this forced debug
  const bool force_debug = true;
  if (simplex_nla.options_->highs_debug_level >= kHighsDebugLevelCostly ||
      force_debug) {
    simplex_nla.lp_ = &incumbent_lp;
    if (debugCheckInvert(simplex_nla, true) == HighsDebugStatus::kError) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Error in basis matrix inverse after solving the LP\n");
      return_status = HighsStatus::kError;
    }
  }
  return return_status;
}

HighsStatus solveLpSimplex(HighsLpSolverObject& solver_object) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsOptions& options = solver_object.options_;
  HighsLp& incumbent_lp = solver_object.lp_;
  HighsSolution& solution = solver_object.solution_;
  HighsModelStatus& unscaled_model_status =
      solver_object.unscaled_model_status_;
  HighsModelStatus& scaled_model_status = solver_object.scaled_model_status_;
  HighsInfo& highs_info = solver_object.highs_info_;
  HighsBasis& basis = solver_object.basis_;

  HEkk& ekk_instance = solver_object.ekk_instance_;
  HighsLp& ekk_lp = ekk_instance.lp_;
  HighsSimplexInfo& ekk_info = ekk_instance.info_;
  SimplexBasis& ekk_basis = ekk_instance.basis_;
  HighsSimplexStatus& status = ekk_instance.status_;
  HSimplexNla& simplex_nla = ekk_instance.simplex_nla_;

  // Copy the simplex iteration count from highs_info_ to ekk_instance, just for
  // convenience
  ekk_instance.iteration_count_ = highs_info.simplex_iteration_count;

  // Reset the model status and HighsInfo values in case of premature
  // return
  resetModelStatusAndHighsInfo(solver_object);

  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLp
  bool positive_num_row = solver_object.lp_.num_row_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "solveLpSimplex called for LP with non-positive (%" HIGHSINT_FORMAT
        ") "
        "number of constraints\n",
        incumbent_lp.num_row_);
    return returnFromSolveLpSimplex(solver_object, HighsStatus::kError);
  }
  // On entry to solveLpSimplex, the incumbent LP is assumed to be
  // unscaled and not moved
  assert(!incumbent_lp.is_scaled_);
  assert(!incumbent_lp.is_moved_);
  // Consider scaling the LP - either with any existing scaling, or by
  // considering computing scaling factors if there are none - and
  // then move to EKK
  considerScaling(options, incumbent_lp);
  // Move the LP to EKK, updating other EKK pointers and any simplex
  // NLA pointers, since they may have moved if the LP has been
  // modified
  ekk_instance.moveLp(std::move(incumbent_lp), &solver_object);
  incumbent_lp.is_moved_ = true;
  if (!status.initialised) {
    // The simplex instance isn't initialised
    call_status = ekk_instance.setup();
    if (call_status == HighsStatus::kError) {
      incumbent_lp.moveBackLpAndUnapplyScaling(ekk_lp);
      return returnFromSolveLpSimplex(solver_object, call_status);
    }
  }
  // If there is no simplex basis, use the HiGHS basis
  if (!status.has_basis && basis.valid) {
    call_status = ekk_instance.setBasis(basis);
    if (call_status == HighsStatus::kError) {
      incumbent_lp.moveBackLpAndUnapplyScaling(ekk_lp);
      return returnFromSolveLpSimplex(solver_object, call_status);
    }
  }
  // These local illegal values are over-written with correct values
  // if the scaled LP is solved in order to take correct algorithmic
  // decisions if the unscaled LP is solved later
  HighsInt num_unscaled_primal_infeasibilities =
      kHighsIllegalInfeasibilityCount;
  HighsInt num_unscaled_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  if (!incumbent_lp.scale_.has_scaling) {
    //
    // Solve the unscaled LP with unscaled NLA
    //
    return_status = ekk_instance.solve();
    //
  } else {
    // Indicate that there is no (current) need to refine the solution
    // by solving the unscaled LP with scaled NLA
    bool refine_solution = false;
    if (options.simplex_unscaled_solution_strategy ==
            kSimplexUnscaledSolutionStrategyNone ||
        options.simplex_unscaled_solution_strategy ==
            kSimplexUnscaledSolutionStrategyRefine) {
      //
      // Solve the scaled LP!
      //
      return_status = ekk_instance.solve();
      //
      if (return_status == HighsStatus::kError) {
        incumbent_lp.moveBackLpAndUnapplyScaling(ekk_lp);
        return returnFromSolveLpSimplex(solver_object, return_status);
      }
      // Copy solution data from the EKK instance
      scaled_model_status = ekk_instance.model_status_;
      highs_info.objective_function_value = ekk_info.primal_objective_value;
      highs_info.simplex_iteration_count = ekk_instance.iteration_count_;
      solution = ekk_instance.getSolution();
      basis = ekk_instance.getHighsBasis();
      assert(basis.valid);
      highs_info.basis_validity = kBasisValidityValid;
      incumbent_lp.moveBackLpAndUnapplyScaling(ekk_lp);
      // Now that the incumbent LP is unscaled, to use the simplex NLA
      // requires scaling to be applied
      simplex_nla.setLpAndScalePointers(&incumbent_lp);
      unscaleSolution(solution, incumbent_lp.scale_);
      // Determine whether the unscaled LP has been solved
      getUnscaledInfeasibilities(options, incumbent_lp.scale_, ekk_basis,
                                 ekk_info, highs_info);
      num_unscaled_primal_infeasibilities =
          highs_info.num_primal_infeasibilities;
      num_unscaled_dual_infeasibilities = highs_info.num_dual_infeasibilities;
      // Determine whether the unscaled solution has infeasibilities
      // after the scaled LP has been solved to optimality
      const bool scaled_optimality_but_unscaled_infeasibilities =
          scaled_model_status == HighsModelStatus::kOptimal &&
          (num_unscaled_primal_infeasibilities ||
           num_unscaled_dual_infeasibilities);
      if (scaled_optimality_but_unscaled_infeasibilities)
        highsLogDev(options.log_options, HighsLogType::kInfo,
                    "Have num/max/sum primal (%" HIGHSINT_FORMAT
                    "/%g/%g) and dual (%" HIGHSINT_FORMAT
                    "/%g/%g) "
                    "unscaled infeasibilities\n",
                    highs_info.num_primal_infeasibilities,
                    highs_info.max_primal_infeasibility,
                    highs_info.sum_primal_infeasibilities,
                    highs_info.num_dual_infeasibilities,
                    highs_info.max_dual_infeasibility,
                    highs_info.sum_dual_infeasibilities);
      // Determine whether refinement will take place
      refine_solution =
          options.simplex_unscaled_solution_strategy ==
              kSimplexUnscaledSolutionStrategyRefine &&
          (scaled_optimality_but_unscaled_infeasibilities ||
           scaled_model_status == HighsModelStatus::kInfeasible ||
           scaled_model_status == HighsModelStatus::kUnboundedOrInfeasible ||
           scaled_model_status == HighsModelStatus::kUnbounded ||
           scaled_model_status == HighsModelStatus::kObjectiveBound ||
           scaled_model_status == HighsModelStatus::kObjectiveTarget);
      // Handle the case when refinement will not take place
      if (!refine_solution) {
        unscaled_model_status = scaled_model_status;
        return_status = highsStatusFromHighsModelStatus(unscaled_model_status);
        return returnFromSolveLpSimplex(solver_object, return_status);
      }
    }
    assert(options.simplex_unscaled_solution_strategy ==
               kSimplexUnscaledSolutionStrategyDirect ||
           refine_solution);
    // Solve the unscaled LP using scaled NLA. This requires pointers of
    // a scaled matrix to be passed to the HFactor instance. Use the
    // incumbent LP for ths.
    //
    // Check that the incumbent LP has been moved back and is unscaled
    assert(!incumbent_lp.is_moved_);
    assert(!incumbent_lp.is_scaled_);
    // Move the incumbent LP
    ekk_instance.moveLp(std::move(incumbent_lp), &solver_object);
    incumbent_lp.is_moved_ = true;
    // Save options/strategies that may be changed
    HighsInt simplex_strategy = options.simplex_strategy;
    double dual_simplex_cost_perturbation_multiplier =
        options.dual_simplex_cost_perturbation_multiplier;
    HighsInt simplex_dual_edge_weight_strategy =
        ekk_info.dual_edge_weight_strategy;
    if (num_unscaled_primal_infeasibilities == 0) {
      // Only dual infeasibilities, so use primal simplex
      options.simplex_strategy = kSimplexStrategyPrimal;
    } else {
      // Using dual simplex, so force Devex if starting from an advanced
      // basis with no steepest edge weights
      //    if (status.has_basis || basis.valid) {
      // ToDo Track whether steepest edge weights are known &&
      // !status.has_dual_steepest_edge_weights) {
      ekk_info.dual_edge_weight_strategy = kSimplexDualEdgeWeightStrategyDevex;
      // options.dual_simplex_cost_perturbation_multiplier = 0;
    }
    //
    // Solve the unscaled LP with scaled NLA
    //
    return_status = ekk_instance.solve();
    //
    // Restore the options/strategies that may have been changed
    options.simplex_strategy = simplex_strategy;
    options.dual_simplex_cost_perturbation_multiplier =
        dual_simplex_cost_perturbation_multiplier;
    ekk_info.dual_edge_weight_strategy = simplex_dual_edge_weight_strategy;
  }
  // Copy solution data from the EKK istance
  //
  scaled_model_status = ekk_instance.model_status_;
  highs_info.objective_function_value = ekk_info.primal_objective_value;
  highs_info.simplex_iteration_count = ekk_instance.iteration_count_;
  solution = ekk_instance.getSolution();
  basis = ekk_instance.getHighsBasis();
  assert(basis.valid);
  highs_info.basis_validity = kBasisValidityValid;
  // Move the incumbent LP back from Ekk
  incumbent_lp = std::move(ekk_lp);
  incumbent_lp.is_moved_ = false;
  simplex_nla.setLpAndScalePointers(&incumbent_lp);
  if (return_status == HighsStatus::kError) {
    return returnFromSolveLpSimplex(solver_object, HighsStatus::kError);
  }
  // The unscaled LP has been solved - either directly, or because
  // there was no scaling. Copy values into the HighsInfo
  // that are set (above) by the call to getUnscaledInfeasibilities
  highs_info.num_primal_infeasibilities = ekk_info.num_primal_infeasibilities;
  highs_info.max_primal_infeasibility = ekk_info.max_primal_infeasibility;
  highs_info.sum_primal_infeasibilities = ekk_info.sum_primal_infeasibilities;
  highs_info.num_dual_infeasibilities = ekk_info.num_dual_infeasibilities;
  highs_info.max_dual_infeasibility = ekk_info.max_dual_infeasibility;
  highs_info.sum_dual_infeasibilities = ekk_info.sum_dual_infeasibilities;
  setSolutionStatus(highs_info);
  unscaled_model_status = scaled_model_status;
  return_status = highsStatusFromHighsModelStatus(unscaled_model_status);
  return returnFromSolveLpSimplex(solver_object, return_status);
}
#endif
