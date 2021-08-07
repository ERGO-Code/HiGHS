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
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsSolve.h"
#include "lp_data/HighsStatus.h"
#include "simplex/HEkk.h"
#include "simplex/HEkkPrimal.h"
#include "simplex/HSimplex.h"
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

void recoverIncumbentAndSimplexLp(bool& incumbent_lp_moved,
                                  bool& incumbent_lp_scaled,
                                  HighsLp& incumbent_lp, HighsLp& ekk_lp,
                                  const SimplexScale& scale) {
  if (incumbent_lp_moved) incumbent_lp = std::move(ekk_lp);
  incumbent_lp_moved = false;
  ekk_lp = incumbent_lp;
  if (incumbent_lp_scaled) {
    const bool force_unscale = true;
    unscaleSimplexLp(incumbent_lp, scale, force_unscale);
    incumbent_lp_scaled = false;
  } else {
    const bool force_scale = true;
    scaleSimplexLp(ekk_lp, scale, force_scale);
  }
}

HighsStatus solveLpSimplex(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsOptions& options = highs_model_object.options_;
  HighsLp& lp = highs_model_object.lp_;
  SimplexScale& scale = highs_model_object.scale_;
  HighsSolution& solution = highs_model_object.solution_;
  HighsModelStatus& unscaled_model_status =
      highs_model_object.unscaled_model_status_;
  HighsModelStatus& scaled_model_status =
      highs_model_object.scaled_model_status_;
  HighsSolutionParams& solution_params = highs_model_object.solution_params_;
  HighsBasis& basis = highs_model_object.basis_;

  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsSimplexStatus& status = ekk_instance.status_;
  HighsLp& ekk_lp = ekk_instance.lp_;
  HSimplexNla& simplex_nla = ekk_instance.simplex_nla_;

  // Reset the model status and solution parameters for the unscaled
  // LP in case of premature return
  resetModelStatusAndSolutionParams(highs_model_object);

  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLp
  bool positive_num_row = highs_model_object.lp_.num_row_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "solveLpSimplex called for LP with non-positive (%" HIGHSINT_FORMAT
        ") "
        "number of constraints\n",
        lp.num_row_);
    return HighsStatus::kError;
  }
  // On entry to solveLpSimplex, the incumbent LP is assumed to be
  // unscaled and, if the simplex instance is initialised, it is
  // assumed to contain the scaled LP.
  //
  // Have to keep track of whether the incumbent LP is scaled, and
  // whether it has been moved to the simplex instance. Currently it
  // is unscaled and has not been moved.
  bool incumbent_lp_scaled = false;
  bool incumbent_lp_moved = false;
  //
  // On return from solveLpSimplex, the incumbent LP must be unscaled
  // and the simplex instance must contain the scaled LP.
  const bool changed_scaling_strategy =
      options.simplex_scale_strategy != kSimplexScaleStrategyChoose &&
      options.simplex_scale_strategy != scale.strategy;
  if (!status.initialised || changed_scaling_strategy) {
    // The simplex instance isn't initialised or the scaling factors
    // for the incumbent LP correspond to a different strategy, so
    // consider finding scaling factors
    getSimplexScaling(options, lp, scale);
    incumbent_lp_scaled = true;
    // Any scaling has been applied to the incumbent LP in the course
    // of finding the scaling factors, so move the LP to Ekk
    ekk_instance.moveNewLp(std::move(lp));
    incumbent_lp_moved = true;
  } else {
    // Update the pointers to the HFactor matrix since they may have
    // moved if the LP has been modified
    ekk_instance.updateFactorMatrixPointers();
  }
  // If there is no simplex basis, use the HiGHS basis
  if (!status.has_basis && basis.valid) {
    return_status = ekk_instance.setBasis(basis);
    if (return_status == HighsStatus::kError) {
      recoverIncumbentAndSimplexLp(incumbent_lp_moved, incumbent_lp_scaled, lp,
                                   ekk_lp, scale);
      assert(!incumbent_lp_moved && !incumbent_lp_scaled);
      return return_status;
    }
  }
  HighsInt num_unscaled_primal_infeasibility = kHighsIllegalInfeasibilityCount;
  HighsInt num_unscaled_dual_infeasibility = kHighsIllegalInfeasibilityCount;
  if (options.simplex_unscaled_solution_strategy ==
          kSimplexUnscaledSolutionStrategyNone ||
      options.simplex_unscaled_solution_strategy ==
          kSimplexUnscaledSolutionStrategyRefine) {
    // Solve the scaled LP!
    return_status = ekk_instance.solve();
    if (return_status == HighsStatus::kError) {
      recoverIncumbentAndSimplexLp(incumbent_lp_moved, incumbent_lp_scaled, lp,
                                   ekk_lp, scale);
      assert(!incumbent_lp_moved && !incumbent_lp_scaled);
      return return_status;
    }
    // Copy solution data into the HMO
    scaled_model_status = ekk_instance.model_status_;
    solution_params.objective_function_value =
        ekk_instance.info_.primal_objective_value;
    highs_model_object.iteration_counts_.simplex +=
        ekk_instance.iteration_count_;
    solution = ekk_instance.getSolution();
    if (scale.is_scaled) unscaleSolution(solution, scale);
    basis = ekk_instance.getHighsBasis();

    // Determine whether the unscaled LP has been solved
    getUnscaledInfeasibilities(options, ekk_lp, ekk_instance.basis_,
                               ekk_instance.info_, scale, solution_params);

    num_unscaled_primal_infeasibility =
        solution_params.num_primal_infeasibility;
    num_unscaled_dual_infeasibility = solution_params.num_dual_infeasibility;

    assert(solution_params.num_primal_infeasibility >= 0);
    assert(solution_params.num_dual_infeasibility >= 0);
    if (num_unscaled_primal_infeasibility) {
      solution_params.primal_solution_status = kSolutionStatusInfeasible;
    } else {
      solution_params.primal_solution_status = kSolutionStatusFeasible;
    }
    if (num_unscaled_dual_infeasibility) {
      solution_params.dual_solution_status = kSolutionStatusInfeasible;
    } else {
      solution_params.dual_solution_status = kSolutionStatusFeasible;
    }
    // Determine whether the unscaled solution has infeasibilities
    // after the scaled LP has been solved to optimality
    const bool scaled_optimality_but_unscaled_infeasibilities =
        scaled_model_status == HighsModelStatus::kOptimal &&
        (num_unscaled_primal_infeasibility || num_unscaled_dual_infeasibility);
    if (scaled_optimality_but_unscaled_infeasibilities)
      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Have num/max/sum primal (%" HIGHSINT_FORMAT
                  "/%g/%g) and dual (%" HIGHSINT_FORMAT
                  "/%g/%g) "
                  "unscaled infeasibilities\n",
                  num_unscaled_primal_infeasibility,
                  solution_params.max_primal_infeasibility,
                  solution_params.sum_primal_infeasibility,
                  num_unscaled_dual_infeasibility,
                  solution_params.max_dual_infeasibility,
                  solution_params.sum_dual_infeasibility);
    // Determine whether refinement will take place
    const bool refine_solution =
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
      recoverIncumbentAndSimplexLp(incumbent_lp_moved, incumbent_lp_scaled, lp,
                                   ekk_lp, scale);
      assert(!incumbent_lp_moved && !incumbent_lp_scaled);
      return return_status;
    }
  }
  assert(options.simplex_unscaled_solution_strategy ==
             kSimplexUnscaledSolutionStrategyDirect ||
         options.simplex_unscaled_solution_strategy ==
             kSimplexUnscaledSolutionStrategyRefine);
  // Solve the unscaled LP using scaled NLA. This requires pointers of
  // a scaled matrix to be passed to the HFactor instance. Use the
  // incumbent LP for ths.
  //
  // Possibly move back the incumbent LP
  if (incumbent_lp_moved) {
    lp = std::move(ekk_lp);
    incumbent_lp_moved = false;
  }

  // Create a scaled LP for the HFactor instance, by taking a copy of
  // the incumbent LP, and ensuring that the incumbent LP is unscaled
  HighsLp scaled_lp = lp;
  if (incumbent_lp_scaled) {
    // The incumbent LP is scaled, so scaled_lp is scaled. Unscale the
    // incumbent LP
    const bool force_unscale = true;
    unscaleSimplexLp(lp, scale, force_unscale);
    incumbent_lp_scaled = false;
  } else {
    // The incumbent LP is unscaled, so scaled_lp is unscaled: scale
    // it
    const bool force_scale = true;
    scaleSimplexLp(scaled_lp, scale, force_scale);
  }
  // Now the incumbent LP is unscaled, and the LP for the HFactor
  // instance is scaled
  //
  // Move the incumbent LP and pass pointers to the scaling factors
  // and scaled matrix for the HFactor instance
  ekk_instance.moveUnscaledLp(std::move(lp), &scale, &scaled_lp.a_matrix_,
			      &scaled_lp.a_start_[0],
                              &scaled_lp.a_index_[0], &scaled_lp.a_value_[0]);
  incumbent_lp_moved = true;
  // Save options/strategies that may be changed
  HighsInt simplex_strategy = options.simplex_strategy;
  double dual_simplex_cost_perturbation_multiplier =
      options.dual_simplex_cost_perturbation_multiplier;
  HighsInt simplex_dual_edge_weight_strategy =
      ekk_instance.info_.dual_edge_weight_strategy;
  if (num_unscaled_primal_infeasibility == 0) {
    // Only dual infeasibilities, so use primal simplex
    options.simplex_strategy = kSimplexStrategyPrimal;
  } else {
    // Using dual simplex, so force Devex if starting from an advanced
    // basis with no steepest edge weights
    //    if (status.has_basis || basis.valid) {
    // ToDo Track whether steepest edge weights are known &&
    // !status.has_dual_steepest_edge_weights) {
    ekk_instance.info_.dual_edge_weight_strategy =
        kSimplexDualEdgeWeightStrategyDevex;
  }
  //    }
  //    options.dual_simplex_cost_perturbation_multiplier = 0;
  //  }
  // Solve the unscaled LP!
  return_status = ekk_instance.solve();
  // Restore the options/strategies that may have been changed
  options.simplex_strategy = simplex_strategy;
  options.dual_simplex_cost_perturbation_multiplier =
      dual_simplex_cost_perturbation_multiplier;
  ekk_instance.info_.dual_edge_weight_strategy =
      simplex_dual_edge_weight_strategy;
  // Move the incumbent LP back from Ekk
  lp = std::move(ekk_lp);
  incumbent_lp_moved = false;
  ekk_instance.passScaledLp(scaled_lp);
  if (return_status == HighsStatus::kError) return HighsStatus::kError;

  // Copy solution data into the HMO
  //
  scaled_model_status = ekk_instance.model_status_;
  solution_params.objective_function_value =
      ekk_instance.info_.primal_objective_value;
  highs_model_object.iteration_counts_.simplex += ekk_instance.iteration_count_;
  solution = ekk_instance.getSolution();
  basis = ekk_instance.getHighsBasis();
  // Copy values into the HighsSolutionParams that are set (above) by
  // the call to getUnscaledInfeasibilities
  solution_params.num_primal_infeasibility =
      ekk_instance.info_.num_primal_infeasibility;
  solution_params.max_primal_infeasibility =
      ekk_instance.info_.max_primal_infeasibility;
  solution_params.sum_primal_infeasibility =
      ekk_instance.info_.sum_primal_infeasibility;
  solution_params.num_dual_infeasibility =
      ekk_instance.info_.num_dual_infeasibility;
  solution_params.max_dual_infeasibility =
      ekk_instance.info_.max_dual_infeasibility;
  solution_params.sum_dual_infeasibility =
      ekk_instance.info_.sum_dual_infeasibility;

  num_unscaled_primal_infeasibility = solution_params.num_primal_infeasibility;
  num_unscaled_dual_infeasibility = solution_params.num_dual_infeasibility;

  assert(solution_params.num_primal_infeasibility >= 0);
  assert(solution_params.num_dual_infeasibility >= 0);
  if (num_unscaled_primal_infeasibility) {
    solution_params.primal_solution_status = kSolutionStatusInfeasible;
  } else {
    solution_params.primal_solution_status = kSolutionStatusFeasible;
  }
  if (num_unscaled_dual_infeasibility) {
    solution_params.dual_solution_status = kSolutionStatusInfeasible;
  } else {
    solution_params.dual_solution_status = kSolutionStatusFeasible;
  }
  unscaled_model_status = scaled_model_status;
  return_status = highsStatusFromHighsModelStatus(unscaled_model_status);
  return return_status;
}
#endif
