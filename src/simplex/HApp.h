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
// It sets the HiGHS basis within highs_model_object and, if optimal,
// the HiGHS solution, too
HighsStatus solveLpSimplex(HighsModelObject& highs_model_object) {
  HighsStatus return_status = HighsStatus::kOk;
  //  HighsStatus call_status;
  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsOptions& options = highs_model_object.options_;
  HighsSimplexStatus& status = ekk_instance.status_;

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
        highs_model_object.lp_.num_row_);
    return HighsStatus::kError;
  }

  // If the simplex LP isn't initialised, scale and pass the current LP
  if (!status.initialised) scaleAndPassLpToEkk(highs_model_object);

  // If there is no simplex basis, use the HiGHS basis
  if (!status.has_basis && highs_model_object.basis_.valid) {
    return_status = ekk_instance.setBasis(highs_model_object.basis_);
    if (return_status == HighsStatus::kError) return HighsStatus::kError;
  }
  // Solve the LP!
  return_status = ekk_instance.solve();
  if (return_status == HighsStatus::kError) return HighsStatus::kError;

  // Copy solution data into the HMO
  HighsSolutionParams& solution_params = highs_model_object.solution_params_;
  highs_model_object.scaled_model_status_ = ekk_instance.model_status_;
  solution_params.objective_function_value =
      ekk_instance.info_.primal_objective_value;
  highs_model_object.iteration_counts_.simplex += ekk_instance.iteration_count_;
  highs_model_object.solution_ = ekk_instance.getSolution();
  if (highs_model_object.scale_.is_scaled)
    unscaleSolution(highs_model_object.solution_, highs_model_object.scale_);
  highs_model_object.basis_ = ekk_instance.getHighsBasis();

  // Determine whether the unscaled LP has been solved
  getUnscaledInfeasibilities(ekk_instance.options_, ekk_instance.lp_,
                             ekk_instance.model_status_, ekk_instance.basis_,
                             ekk_instance.info_, highs_model_object.scale_,
                             solution_params);

  HighsInt& num_unscaled_primal_infeasibility =
      solution_params.num_primal_infeasibility;
  HighsInt& num_unscaled_dual_infeasibility =
      solution_params.num_dual_infeasibility;

  if (num_unscaled_primal_infeasibility > 0) {
    solution_params.primal_solution_status = kSolutionStatusInfeasible;
  } else {
    solution_params.primal_solution_status = kSolutionStatusFeasible;
  }
  if (num_unscaled_dual_infeasibility > 0) {
    solution_params.dual_solution_status = kSolutionStatusInfeasible;
  } else {
    solution_params.dual_solution_status = kSolutionStatusFeasible;
  }

  // Handle non-optimal status
  if (ekk_instance.model_status_ != HighsModelStatus::kOptimal) {
    highs_model_object.unscaled_model_status_ = ekk_instance.model_status_;
    return_status = highsStatusFromHighsModelStatus(ekk_instance.model_status_);
    return return_status;
  }

  // Now interpret the status of the unscaled solution when the scaled
  // LP is solved to optimailty
  assert(ekk_instance.model_status_ == HighsModelStatus::kOptimal);

  // Set the model and solution status according to the unscaled solution
  // parameters
  if (num_unscaled_primal_infeasibility == 0 &&
      num_unscaled_dual_infeasibility == 0) {
    // Optimal
    highs_model_object.unscaled_model_status_ = HighsModelStatus::kOptimal;
  } else {
    // Not optimal - should try refinement
    assert(num_unscaled_primal_infeasibility > 0 ||
           num_unscaled_dual_infeasibility > 0);
    highs_model_object.unscaled_model_status_ = HighsModelStatus::kNotset;
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
    HighsLp scaled_lp;
    if (options.simplex_unscaled_solution_strategy ==
        kSimplexUnscaledSolutionStrategyRefine) {
      highsLogUser(options.log_options, HighsLogType::kInfo,
                   "Removing unscaled primal/dual infeasibilities\n");
      // Move back the scaled LP
      HighsLp& lp = highs_model_object.lp_;
      HighsLp& ekk_lp = ekk_instance.lp_;
      HSimplexNla& simplex_nla = ekk_instance.simplex_nla_;
      lp = std::move(ekk_lp);
      // Take a copy of the scaled LP
      scaled_lp = lp;
      // Unscale the LP
      unscaleSimplexLp(lp, highs_model_object.scale_);
      // Move the unscaled LP to Ekk
      ekk_lp = std::move(lp);
      // Pass scaling factors and scaled matrix pointers to the
      // simplex NLA
      simplex_nla.passScaleAndFactorMatrixPointers(
          &highs_model_object.scale_, &scaled_lp.a_start_[0],
          &scaled_lp.a_index_[0], &scaled_lp.a_value_[0]);
      // Force reinitialisation of the matrix for the simplex solver
      // now that the LP is unscaled
      ekk_instance.initialiseMatrix(true);
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
        // Using dual simplex, so force Devex
        options.dual_simplex_cost_perturbation_multiplier = 0;
        ekk_instance.info_.dual_edge_weight_strategy =
            kSimplexDualEdgeWeightStrategyDevex;
      }
      // Solve the unscaled LP!
      return_status = ekk_instance.solve();
      // Restore the options/strategies that may have been changed
      options.simplex_strategy = simplex_strategy;
      options.dual_simplex_cost_perturbation_multiplier =
          dual_simplex_cost_perturbation_multiplier;
      ekk_instance.info_.dual_edge_weight_strategy =
          simplex_dual_edge_weight_strategy;
      // Restore the unscaled LP
      lp = std::move(ekk_lp);
      // Move the scaled LP to Ekk to leave HiGHS and the simplex
      // solver in their traditional state
      ekk_lp = std::move(scaled_lp);

      if (return_status == HighsStatus::kError) return HighsStatus::kError;

      // Tidy up as above - until refinement is directed from calling method
      //
      // Copy solution data into the HMO
      highs_model_object.scaled_model_status_ = ekk_instance.model_status_;
      solution_params.objective_function_value =
          ekk_instance.info_.primal_objective_value;
      highs_model_object.iteration_counts_.simplex +=
          ekk_instance.iteration_count_;
      highs_model_object.solution_ = ekk_instance.getSolution();
      highs_model_object.basis_ = ekk_instance.getHighsBasis();

      // Determine whether the unscaled LP has been solved
      const bool scaled_simplex_lp = false;
      getUnscaledInfeasibilities(
          ekk_instance.options_, ekk_instance.lp_, ekk_instance.model_status_,
          ekk_instance.basis_, ekk_instance.info_, highs_model_object.scale_,
          solution_params, scaled_simplex_lp);
      if (num_unscaled_primal_infeasibility > 0) {
        solution_params.primal_solution_status = kSolutionStatusInfeasible;
      } else {
        solution_params.primal_solution_status = kSolutionStatusFeasible;
      }
      if (num_unscaled_dual_infeasibility > 0) {
        solution_params.dual_solution_status = kSolutionStatusInfeasible;
      } else {
        solution_params.dual_solution_status = kSolutionStatusFeasible;
      }

      // Handle non-optimal status
      if (ekk_instance.model_status_ != HighsModelStatus::kOptimal) {
        highs_model_object.unscaled_model_status_ = ekk_instance.model_status_;
        return_status =
            highsStatusFromHighsModelStatus(ekk_instance.model_status_);
        return return_status;
      }

      // Now interpret the status of the unscaled solution when the scaled
      // LP is solved to optimailty
      assert(ekk_instance.model_status_ == HighsModelStatus::kOptimal);

      highsLogDev(options.log_options, HighsLogType::kInfo,
                  "Now have num/max/sum primal (%" HIGHSINT_FORMAT
                  "/%g/%g) and dual (%" HIGHSINT_FORMAT
                  "/%g/%g) "
                  "unscaled infeasibilities\n",
                  num_unscaled_primal_infeasibility,
                  solution_params.max_primal_infeasibility,
                  solution_params.sum_primal_infeasibility,
                  num_unscaled_dual_infeasibility,
                  solution_params.max_dual_infeasibility,
                  solution_params.sum_dual_infeasibility);

      // Set the model and solution status according to the unscaled solution
      // parameters
      if (num_unscaled_primal_infeasibility == 0 &&
          num_unscaled_dual_infeasibility == 0) {
        // Optimal
        highs_model_object.unscaled_model_status_ = HighsModelStatus::kOptimal;
      } else {
        // Should be optimal!
        highsLogDev(options.log_options, HighsLogType::kError,
                    "Still have num/max/sum primal (%" HIGHSINT_FORMAT
                    "/%g/%g) and dual (%" HIGHSINT_FORMAT
                    "/%g/%g) "
                    "unscaled infeasibilities\n",
                    num_unscaled_primal_infeasibility,
                    solution_params.max_primal_infeasibility,
                    solution_params.sum_primal_infeasibility,
                    num_unscaled_dual_infeasibility,
                    solution_params.max_dual_infeasibility,
                    solution_params.sum_dual_infeasibility);
        assert(num_unscaled_primal_infeasibility == 0 &&
               num_unscaled_dual_infeasibility == 0);
        return HighsStatus::kError;
      }
    }
    //    highs_model_object.solution_ = ekk_instance.getSolution();
    //    if (highs_model_object.scale_.is_scaled)
    //      unscaleSolution(highs_model_object.solution_,
    //      highs_model_object.scale_);
    highs_model_object.basis_ = ekk_instance.getHighsBasis();
  }
  return return_status;
}
#endif
