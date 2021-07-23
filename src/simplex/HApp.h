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

#ifdef OPENMP
#include "omp.h"
#endif

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
  bool positive_num_row = highs_model_object.lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    highsLogUser(
        options.log_options, HighsLogType::kError,
        "solveLpSimplex called for LP with non-positive (%" HIGHSINT_FORMAT
        ") "
        "number of constraints\n",
        highs_model_object.lp_.numRow_);
    return HighsStatus::kError;
  }

  // If the simplex LP isn't initialised, scale and pass the current LP
  if (!status.initialised) scaleAndPassLpToEkk(highs_model_object);

  // If there is no simplex basis, use the HiGHS basis
  if (!status.has_basis && highs_model_object.basis_.valid) {
    return_status = ekk_instance.setBasis(highs_model_object.basis_);
    if (return_status == HighsStatus::kError) return HighsStatus::kError;
  }
  HighsInt simplex_iteration_limit = options.simplex_iteration_limit;
  //  if (kRefineSimplex) options.simplex_iteration_limit = 3;
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
  double new_primal_feasibility_tolerance;
  double new_dual_feasibility_tolerance;
  getUnscaledInfeasibilitiesAndNewTolerances(
      ekk_instance.options_, ekk_instance.lp_, ekk_instance.model_status_,
      ekk_instance.basis_, ekk_instance.info_, highs_model_object.scale_,
      solution_params, new_primal_feasibility_tolerance,
      new_dual_feasibility_tolerance);

  HighsInt num_unscaled_primal_infeasibility =
      solution_params.num_primal_infeasibility;
  HighsInt num_unscaled_dual_infeasibility =
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
    if (kRefineSimplex) {
      printf("Ignoring non-optimal status %s\n", utilModelStatusToString(ekk_instance.model_status_).c_str());
    } else {
      highs_model_object.unscaled_model_status_ = ekk_instance.model_status_;
      return_status = highsStatusFromHighsModelStatus(ekk_instance.model_status_);
      return return_status;
    }
  }

  // Now interpret the status of the unscaled solution when the scaled
  // LP is solved to optimailty
  if (kRefineSimplex) {
    printf("Ignoring assert(ekk_instance.model_status_ == HighsModelStatus::kOptimal)\n");
  } else {
    assert(ekk_instance.model_status_ == HighsModelStatus::kOptimal);
  }

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
    highsLogDev(highs_model_object.options_.log_options, HighsLogType::kInfo,
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
    if (kRefineSimplex) {
      const bool unscale = true;
      HighsLp scaled_lp;
      if (unscale) {
	// Move back the scaled LP
	HighsLp& lp = highs_model_object.lp_;
	HighsLp& ekk_lp = highs_model_object.ekk_instance_.lp_;
	HSimplexNla& simplex_nla = highs_model_object.ekk_instance_.simplex_nla_;
	lp = std::move(ekk_lp);
	// Take a copy of the scaled LP
        scaled_lp = lp;
	// Unscale the LP
	unscaleSimplexLp(lp, highs_model_object.scale_);
	// Move the unscaled LP to Ekk
	ekk_lp = std::move(lp);
	// Pass scaling factors and scaled matrix pointers to the
	// simplex NLA
	simplex_nla.passScaleAndMatrixPointers(&highs_model_object.scale_,
					       &scaled_lp.Astart_[0], 
					       &scaled_lp.Aindex_[0], 
					       &scaled_lp.Avalue_[0]);
	// Reinitialise the matrix for the simplex solver now that the
	// LP is unscaled
	ekk_instance.status_.has_matrix = false;
	highs_model_object.ekk_instance_.initialiseMatrix();
      }
      // Restore the iteration limit (temp)
      options.simplex_iteration_limit = simplex_iteration_limit;
      // Save the dual simplex cost perturbation multiplier and set
      // the option to zero
      double dual_simplex_cost_perturbation_multiplier =
	options.dual_simplex_cost_perturbation_multiplier;
      options.dual_simplex_cost_perturbation_multiplier = 0;
      // Save the simplex dual edge weight strategy and set
      // the option to Devex
      HighsInt simplex_dual_edge_weight_strategy = ekk_instance.info_.dual_edge_weight_strategy;
      ekk_instance.info_.dual_edge_weight_strategy = kSimplexDualEdgeWeightStrategyDevex;
      return_status = ekk_instance.solve();
      // Restore the dual simplex cost perturbation multiplier and
      // simplex dual edge weight strategy
      options.dual_simplex_cost_perturbation_multiplier = dual_simplex_cost_perturbation_multiplier;
      ekk_instance.info_.dual_edge_weight_strategy = simplex_dual_edge_weight_strategy;

      printf("Abort in solveLpSimplex\n"); fflush(stdout);
      abort();
    } else {
      if (ekk_instance.model_status_ == HighsModelStatus::kOptimal)
	highsLogDev(highs_model_object.options_.log_options, HighsLogType::kInfo,
		    "Possibly re-solve with feasibility tolerances of %g "
		    "primal and %g dual\n",
		    new_primal_feasibility_tolerance,
		    new_dual_feasibility_tolerance);
    }
    highs_model_object.solution_ = ekk_instance.getSolution();
    if (highs_model_object.scale_.is_scaled)
      unscaleSolution(highs_model_object.solution_, highs_model_object.scale_);
    highs_model_object.basis_ = ekk_instance.getHighsBasis();
  }
  return return_status;
}
#endif
