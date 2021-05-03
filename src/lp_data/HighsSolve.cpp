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
/**@file lp_data/HighsSolve.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "lp_data/HighsInfo.h"
#include "lp_data/HighsModelObject.h"
#include "lp_data/HighsSolution.h"
#include "simplex/HApp.h"
#include "util/HighsUtils.h"
#ifdef IPX_ON
#include "ipm/IpxWrapper.h"
#else
#include "ipm/IpxWrapperEmpty.h"
#endif

// The method below runs simplex or ipx solver on the lp.
HighsStatus solveLp(HighsModelObject& model, const string message) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsOptions& options = model.options_;
  // Reset unscaled model status and solution params - except for
  // iteration counts
  resetModelStatusAndSolutionParams(model);
  highsLogUser(options.log_options, HighsLogType::kInfo,
               (message + "\n").c_str());
#ifdef HIGHSDEV
  // Shouldn't have to check validity of the LP since this is done when it is
  // loaded or modified
  call_status = assessLp(model.lp_, options_);
  // If any errors have been found or normalisation carried out,
  // call_status will be ERROR or WARNING, so only valid return is OK.
  assert(call_status == HighsStatus::kOk);
  return_status = interpretCallStatus(call_status, return_status, "assessLp");
  if (return_status == HighsStatus::kError) return return_status;
#endif
  if (!model.lp_.numRow_) {
    // Unconstrained LP so solve directly
    call_status = solveUnconstrainedLp(model);
    return_status =
        interpretCallStatus(call_status, return_status, "solveUnconstrainedLp");
    if (return_status == HighsStatus::kError) return return_status;
    // Set the scaled model status for completeness
    model.scaled_model_status_ = model.unscaled_model_status_;
  } else if (options.solver == kIpmString) {
    // Use IPM
#ifdef IPX_ON
    bool imprecise_solution;
    // Use IPX to solve the LP
    try {
      call_status = solveLpIpx(options, model.timer_, model.lp_, imprecise_solution, model.basis_,
			       model.solution_, model.iteration_counts_, model.unscaled_model_status_,
			       model.solution_params_);
    } catch (const std::runtime_error& exception) {
      highsLogDev(options.log_options, HighsLogType::kError,
		  "Runtime error %s in solveLpIpx\n", exception.what());
      call_status = HighsStatus::kError;
    }
    return_status =
        interpretCallStatus(call_status, return_status, "solveLpIpx");
    if (return_status == HighsStatus::kError) return return_status;
    // Non-error return requires a primal solution
    assert(model.solution_.value_valid);
    // Set the scaled model status for completeness
    model.scaled_model_status_ = model.unscaled_model_status_;
    // Get the infeasibilities and objective value
    // ToDo: This should take model.basis_ and use it if it's valid
    getPrimalDualInfeasibilities(model.lp_, model.solution_, model.solution_params_);
    model.solution_params_.objective_function_value =
      computeObjectiveValue(model.lp_, model.solution_);
    if ((model.unscaled_model_status_ == HighsModelStatus::kUnknown ||
	 (model.unscaled_model_status_ ==
	  HighsModelStatus::kUnboundedOrInfeasible && !options.allow_unbounded_or_infeasible)) &&
	options.run_crossover) {
      // IPX has returned a model status that HiGHS would rather
      // avoid, so perform simplex clean-up if crossover was allowed.
      //
      // This is an unusual situation, and the cost will usually be
      // acceptable. Worst case is if crossover wasn't run, in which
      // case there's no basis to start simplex
      //
      // ToDo: Check whether simplex can exploit the primal solution returned by IPX
      highsLogUser(options.log_options, HighsLogType::kWarning,
		   "Imprecise solution returned from IPX, so use simplex to "
		   "clean up\n");
      // Reset the return status since it will now be determined by
      // the outcome of the simplex solve
      return_status = HighsStatus::kOk;
      call_status = solveLpSimplex(model);
      return_status =
	interpretCallStatus(call_status, return_status, "solveLpSimplex");
      if (return_status == HighsStatus::kError) return return_status;
      if (!isSolutionRightSize(model.lp_, model.solution_)) {
	highsLogUser(options.log_options, HighsLogType::kError,
		     "Inconsistent solution returned from solver\n");
	return HighsStatus::kError;
      }
    }
#else
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Model cannot be solved with IPM\n");
    return HighsStatus::kError;
#endif
  } else {
    // Use Simplex
    call_status = solveLpSimplex(model);
    return_status =
        interpretCallStatus(call_status, return_status, "solveLpSimplex");
    if (return_status == HighsStatus::kError) return return_status;
    if (!isSolutionRightSize(model.lp_, model.solution_)) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Inconsistent solution returned from solver\n");
      return HighsStatus::kError;
    }
  }
  // ToDo: Possibly analyse the HiGHS (basic) solution - generalise to
  // handle solution without basis? But doesn't returnFromRun do this?
  //
  // NB IPX may not yield a basic solution
  if (model.basis_.valid) debugHighsBasicSolution(message, model);

  return return_status;
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsSolutionParams
HighsStatus solveUnconstrainedLp(HighsModelObject& highs_model_object) {
  return (solveUnconstrainedLp(
      highs_model_object.options_, highs_model_object.lp_,
      highs_model_object.unscaled_model_status_,
      highs_model_object.solution_params_, highs_model_object.solution_,
      highs_model_object.basis_));
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsSolutionParams
HighsStatus solveUnconstrainedLp(const HighsOptions& options, const HighsLp& lp,
                                 HighsModelStatus& model_status,
                                 HighsSolutionParams& solution_params,
                                 HighsSolution& solution, HighsBasis& basis) {
  // Aliase to model status and solution parameters
  resetModelStatusAndSolutionParams(model_status, solution_params, options);

  // Check that the LP really is unconstrained!
  assert(lp.numRow_ == 0);
  if (lp.numRow_ != 0) return HighsStatus::kError;

  highsLogUser(options.log_options, HighsLogType::kInfo,
               "Solving an unconstrained LP with %" HIGHSINT_FORMAT
               " columns\n",
               lp.numCol_);

  solution.col_value.assign(lp.numCol_, 0);
  solution.col_dual.assign(lp.numCol_, 0);
  basis.col_status.assign(lp.numCol_, HighsBasisStatus::kNonbasic);

  double primal_feasibility_tolerance =
      solution_params.primal_feasibility_tolerance;
  double dual_feasibility_tolerance =
      solution_params.dual_feasibility_tolerance;

  // Initialise the objective value calculation. Done using
  // HighsSolution so offset is vanilla
  double objective = lp.offset_;
  bool infeasible = false;
  bool unbounded = false;

  solution_params.num_primal_infeasibility = 0;
  solution_params.max_primal_infeasibility = 0;
  solution_params.sum_primal_infeasibility = 0;
  solution_params.num_dual_infeasibility = 0;
  solution_params.max_dual_infeasibility = 0;
  solution_params.sum_dual_infeasibility = 0;

  for (HighsInt iCol = 0; iCol < lp.numCol_; iCol++) {
    double cost = lp.colCost_[iCol];
    double dual = (HighsInt)lp.sense_ * cost;
    double lower = lp.colLower_[iCol];
    double upper = lp.colUpper_[iCol];
    double value;
    double primal_infeasibility = 0;
    HighsBasisStatus status = HighsBasisStatus::kNonbasic;
    if (lower > upper) {
      // Inconsistent bounds, so set the variable to lower bound,
      // unless it's infinite. Otherwise set the variable to upper
      // bound, unless it's infinite. Otherwise set the variable to
      // zero.
      if (highs_isInfinity(lower)) {
        // Lower bound of +inf
        if (highs_isInfinity(-upper)) {
          // Unite upper bound of -inf
          value = 0;
          status = HighsBasisStatus::kZero;
          primal_infeasibility = kHighsInf;
        } else {
          value = upper;
          status = HighsBasisStatus::kUpper;
          primal_infeasibility = lower - value;
        }
      } else {
        value = lower;
        status = HighsBasisStatus::kLower;
        primal_infeasibility = value - upper;
      }
    } else if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free column: must have zero cost
      value = 0;
      status = HighsBasisStatus::kZero;
      if (fabs(dual) > dual_feasibility_tolerance) unbounded = true;
    } else if (dual >= dual_feasibility_tolerance) {
      // Column with sufficiently positive dual: set to lower bound
      // and check for unboundedness
      if (highs_isInfinity(-lower)) unbounded = true;
      value = lower;
      status = HighsBasisStatus::kLower;
    } else if (dual <= -dual_feasibility_tolerance) {
      // Column with sufficiently negative dual: set to upper bound
      // and check for unboundedness
      if (highs_isInfinity(upper)) unbounded = true;
      value = upper;
      status = HighsBasisStatus::kUpper;
    } else {
      // Column with sufficiently small dual: set to lower bound (if
      // finite) otherwise upper bound
      if (highs_isInfinity(-lower)) {
        value = upper;
        status = HighsBasisStatus::kUpper;
      } else {
        value = lower;
        status = HighsBasisStatus::kLower;
      }
    }
    assert(status != HighsBasisStatus::kNonbasic);
    solution.col_value[iCol] = value;
    solution.col_dual[iCol] = (HighsInt)lp.sense_ * dual;
    basis.col_status[iCol] = status;
    objective += value * cost;
    solution_params.sum_primal_infeasibility += primal_infeasibility;
    if (primal_infeasibility > primal_feasibility_tolerance) {
      infeasible = true;
      solution_params.num_primal_infeasibility++;
      solution_params.max_primal_infeasibility =
          max(primal_infeasibility, solution_params.max_primal_infeasibility);
    }
  }
  solution_params.objective_function_value = objective;
  basis.valid = true;

  if (infeasible) {
    model_status = HighsModelStatus::kInfeasible;
    solution_params.primal_status = kHighsPrimalDualStatusInfeasiblePoint;
    solution_params.dual_status = kHighsPrimalDualStatusUnknown;
  } else {
    solution_params.primal_status = kHighsPrimalDualStatusFeasiblePoint;
    if (unbounded) {
      model_status = HighsModelStatus::kUnbounded;
      solution_params.dual_status = kHighsPrimalDualStatusInfeasiblePoint;
    } else {
      model_status = HighsModelStatus::kOptimal;
      solution_params.dual_status = kHighsPrimalDualStatusFeasiblePoint;
    }
  }
  return HighsStatus::kOk;
}
