/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
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
  HighsStatus return_status = HighsStatus::OK;
  //  HighsStatus call_status;
  HEkk& ekk_instance = highs_model_object.ekk_instance_;
  HighsOptions& options = highs_model_object.options_;
  HighsSimplexLpStatus& simplex_lp_status = ekk_instance.simplex_lp_status_;

  // Reset the model status and solution parameters for the unscaled
  // LP in case of premature return
  resetModelStatusAndSolutionParams(highs_model_object);

  // Assumes that the LP has a positive number of rows, since
  // unconstrained LPs should be solved in solveLpSimplex
  bool positive_num_row = highs_model_object.lp_.numRow_ > 0;
  assert(positive_num_row);
  if (!positive_num_row) {
    highsLogUser(options.log_options, HighsLogType::ERROR,
                 "solveLpEkkSimplex called for LP with non-positive (%d) "
                 "number of constraints\n",
                 highs_model_object.lp_.numRow_);
    return HighsStatus::Error;
  }

  // If the simplex LP isn't initialised, scale and pass the current LP
  if (!simplex_lp_status.initialised) scaleAndPassLpToEkk(highs_model_object);

  // If there is no simplex basis, use the HiGHS basis
  if (!simplex_lp_status.has_basis && highs_model_object.basis_.valid_) {
    return_status = ekk_instance.setBasis(highs_model_object.basis_);
    if (return_status == HighsStatus::Error) return HighsStatus::Error;
  }

  // Solve the LP!
  return_status = ekk_instance.solve();
  if (return_status == HighsStatus::Error) return HighsStatus::Error;

  // Copy solution data into the HMO
  HighsSolutionParams& solution_params = highs_model_object.solution_params_;
  highs_model_object.scaled_model_status_ = ekk_instance.scaled_model_status_;
  solution_params.objective_function_value =
      ekk_instance.simplex_info_.primal_objective_value;
  highs_model_object.iteration_counts_.simplex += ekk_instance.iteration_count_;
  highs_model_object.solution_ = ekk_instance.getSolution();
  if (highs_model_object.scale_.is_scaled_)
    unscaleSolution(highs_model_object.solution_, highs_model_object.scale_);
  highs_model_object.basis_ = ekk_instance.getHighsBasis();

  // Determine whether the unscaled LP has been solved
  double new_primal_feasibility_tolerance;
  double new_dual_feasibility_tolerance;
  //  HighsSolutionParams solution_;
  getUnscaledInfeasibilitiesAndNewTolerances(
      ekk_instance.options_, ekk_instance.simplex_lp_,
      ekk_instance.scaled_model_status_, ekk_instance.simplex_basis_,
      ekk_instance.simplex_info_, highs_model_object.scale_, solution_params,
      new_primal_feasibility_tolerance, new_dual_feasibility_tolerance);

  // Handle non-optimal status
  if (ekk_instance.scaled_model_status_ != HighsModelStatus::OPTIMAL) {
    highs_model_object.unscaled_model_status_ =
        ekk_instance.scaled_model_status_;
    return_status =
        highsStatusFromHighsModelStatus(ekk_instance.scaled_model_status_);
    return return_status;
  }

  // Now interpret the status of the unscaled solution when the scaled
  // LP is solved to optimailty
  assert(ekk_instance.scaled_model_status_ == HighsModelStatus::OPTIMAL);

  HighsInt num_unscaled_primal_infeasibility =
      solution_params.num_primal_infeasibility;
  HighsInt num_unscaled_dual_infeasibility = solution_params.num_dual_infeasibility;
  // Set the model and solution status according to the unscaled solution
  // parameters
  if (num_unscaled_primal_infeasibility == 0 &&
      num_unscaled_dual_infeasibility == 0) {
    // Optimal
    highs_model_object.unscaled_model_status_ = HighsModelStatus::OPTIMAL;
    solution_params.primal_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
    solution_params.dual_status = PrimalDualStatus::STATUS_FEASIBLE_POINT;
  } else {
    // Not optimal - should try refinement
    highs_model_object.unscaled_model_status_ = HighsModelStatus::NOTSET;
    assert(num_unscaled_primal_infeasibility > 0 ||
           num_unscaled_dual_infeasibility > 0);
    highsLogUser(highs_model_object.options_.log_options, HighsLogType::INFO,
                 "Have num/max/sum primal (%d/%g/%g) and dual (%d/%g/%g) "
                 "unscaled infeasibilities\n",
                 num_unscaled_primal_infeasibility,
                 solution_params.max_primal_infeasibility,
                 solution_params.sum_primal_infeasibility,
                 num_unscaled_dual_infeasibility,
                 solution_params.max_dual_infeasibility,
                 solution_params.sum_dual_infeasibility);
    if (ekk_instance.scaled_model_status_ == HighsModelStatus::OPTIMAL)
      highsLogUser(highs_model_object.options_.log_options, HighsLogType::INFO,
                   "Possibly re-solve with feasibility tolerances of %g "
                   "primal and %g dual\n",
                   new_primal_feasibility_tolerance,
                   new_dual_feasibility_tolerance);
    highs_model_object.solution_ = ekk_instance.getSolution();
    if (highs_model_object.scale_.is_scaled_)
      unscaleSolution(highs_model_object.solution_, highs_model_object.scale_);
    highs_model_object.basis_ = ekk_instance.getHighsBasis();
  }
  return return_status;
}
#endif
