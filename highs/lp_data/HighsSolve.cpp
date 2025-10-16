/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolve.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include <sstream>

#include "ipm/IpxWrapper.h"
#include "lp_data/HighsSolutionDebug.h"
#include "pdlp/CupdlpWrapper.h"
#include "simplex/HApp.h"

// The method below runs simplex, ipx or pdlp solver on the lp.
HighsStatus solveLp(HighsLpSolverObject& solver_object, const string message) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsOptions& options = solver_object.options_;
  // Reset unscaled model status and solution params - except for
  // iteration counts
  resetModelStatusAndHighsInfo(solver_object);
  highsLogUser(options.log_options, HighsLogType::kInfo,
               (message + "\n").c_str());
  if (options.highs_debug_level > kHighsDebugLevelMin) {
    // Shouldn't have to check validity of the LP since this is done when it is
    // loaded or modified
    call_status = assessLp(solver_object.lp_, options);
    // If any errors have been found or normalisation carried out,
    // call_status will be ERROR or WARNING, so only valid return is OK.
    assert(call_status == HighsStatus::kOk);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "assessLp");
    if (return_status == HighsStatus::kError) return return_status;
  }
  if (!solver_object.lp_.num_row_ || solver_object.lp_.a_matrix_.numNz() == 0) {
    // LP is unconstrained due to having no rows or a zero constraint
    // matrix, so solve directly
    call_status = solveUnconstrainedLp(solver_object);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "solveUnconstrainedLp");
    if (return_status == HighsStatus::kError) return return_status;
  } else if (options.solver == kIpmString || options.run_centring ||
             options.solver == kPdlpString) {
    // Use IPM or PDLP
    if (options.solver == kIpmString || options.run_centring) {
      // Use IPX to solve the LP
      try {
        call_status = solveLpIpx(solver_object);
      } catch (const std::exception& exception) {
        highsLogDev(options.log_options, HighsLogType::kError,
                    "Exception %s in solveLpIpx\n", exception.what());
        call_status = HighsStatus::kError;
      }
      return_status = interpretCallStatus(options.log_options, call_status,
                                          return_status, "solveLpIpx");
    } else {
      // Use cuPDLP-C to solve the LP
      try {
        call_status = solveLpCupdlp(solver_object);
      } catch (const std::exception& exception) {
        highsLogDev(options.log_options, HighsLogType::kError,
                    "Exception %s in solveLpCupdlp\n", exception.what());
        call_status = HighsStatus::kError;
      }
      return_status = interpretCallStatus(options.log_options, call_status,
                                          return_status, "solveLpCupdlp");
    }
    // Check for error return
    if (return_status == HighsStatus::kError) return return_status;

    // Non-error return requires a primal solution
    assert(solver_object.solution_.value_valid);

    if (options.solver == kIpmString || options.run_centring) {
      // Setting the IPM-specific values of (highs_)info_ has been done in
      // solveLpIpx
      const bool unwelcome_ipx_status =
          solver_object.model_status_ == HighsModelStatus::kUnknown ||
          (solver_object.model_status_ ==
               HighsModelStatus::kUnboundedOrInfeasible &&
           !options.allow_unbounded_or_infeasible);
      if (unwelcome_ipx_status) {
        // When performing an analytic centre calculation, the setting
        // of options.run_crossover is ignored, so simplex clean-up is
        // not possible - or desirable, anyway!
        highsLogUser(
            options.log_options, HighsLogType::kWarning,
            "Unwelcome IPX status of %s: basis is %svalid; solution is "
            "%svalid; run_crossover is \"%s\"\n",
            utilModelStatusToString(solver_object.model_status_).c_str(),
            solver_object.basis_.valid ? "" : "not ",
            solver_object.solution_.value_valid ? "" : "not ",
            options.run_centring ? kHighsOffString.c_str()
                                 : options.run_crossover.c_str());
        const bool allow_simplex_cleanup =
            options.run_crossover != kHighsOffString && !options.run_centring;
        if (allow_simplex_cleanup) {
          // IPX has returned a model status that HiGHS would rather
          // avoid, so perform simplex clean-up if crossover was allowed.
          //
          // This is an unusual situation, and the cost will usually be
          // acceptable. Worst case is if crossover wasn't run, in which
          // case there's no basis to start simplex
          //
          // ToDo: Check whether simplex can exploit the primal solution
          // returned by IPX
          highsLogUser(options.log_options, HighsLogType::kWarning,
                       "IPX solution is imprecise, so clean up with simplex\n");
          // Reset the return status since it will now be determined by
          // the outcome of the simplex solve
          return_status = HighsStatus::kOk;
          call_status = solveLpSimplex(solver_object);
          return_status = interpretCallStatus(options.log_options, call_status,
                                              return_status, "solveLpSimplex");
          if (return_status == HighsStatus::kError) return return_status;
          if (!isSolutionRightSize(solver_object.lp_,
                                   solver_object.solution_)) {
            highsLogUser(options.log_options, HighsLogType::kError,
                         "Inconsistent solution returned from solver\n");
            return HighsStatus::kError;
          }
        }  // options.run_crossover == kHighsOnString
           // clang-format off
      }  // unwelcome_ipx_status
      // clang-format on
    }
  } else {
    // Use Simplex
    call_status = solveLpSimplex(solver_object);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "solveLpSimplex");
    if (return_status == HighsStatus::kError) return return_status;
    if (!isSolutionRightSize(solver_object.lp_, solver_object.solution_)) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Inconsistent solution returned from solver\n");
      return HighsStatus::kError;
    }
  }
  // Analyse the HiGHS (basic) solution
  if (debugHighsLpSolution(message, solver_object) ==
      HighsDebugStatus::kLogicalError)
    return_status = HighsStatus::kError;
  return return_status;
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsInfo
HighsStatus solveUnconstrainedLp(HighsLpSolverObject& solver_object) {
  return (solveUnconstrainedLp(solver_object.options_, solver_object.lp_,
                               solver_object.model_status_,
                               solver_object.highs_info_,
                               solver_object.solution_, solver_object.basis_));
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsInfo
HighsStatus solveUnconstrainedLp(const HighsOptions& options, const HighsLp& lp,
                                 HighsModelStatus& model_status,
                                 HighsInfo& highs_info, HighsSolution& solution,
                                 HighsBasis& basis) {
  // Aliase to model status and solution parameters
  resetModelStatusAndHighsInfo(model_status, highs_info);

  // Check that the LP really is unconstrained!
  assert(lp.num_row_ == 0 || lp.a_matrix_.numNz() == 0);
  if (lp.num_row_ > 0) {
    // LP has rows, but should only be here if the constraint matrix
    // is zero
    if (lp.a_matrix_.numNz() > 0) return HighsStatus::kError;
  }

  highsLogUser(options.log_options, HighsLogType::kInfo,
               "Solving an unconstrained LP with %" HIGHSINT_FORMAT
               " columns\n",
               lp.num_col_);
  solution.col_value.assign(lp.num_col_, 0);
  solution.col_dual.assign(lp.num_col_, 0);
  basis.col_status.assign(lp.num_col_, HighsBasisStatus::kNonbasic);
  // No rows for primal solution, dual solution or basis
  solution.row_value.clear();
  solution.row_dual.clear();
  basis.row_status.clear();

  double primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  double dual_feasibility_tolerance = options.dual_feasibility_tolerance;

  // Initialise the objective value calculation. Done using
  // HighsSolution so offset is vanilla
  double objective = lp.offset_;

  highs_info.num_primal_infeasibilities = 0;
  highs_info.max_primal_infeasibility = 0;
  highs_info.sum_primal_infeasibilities = 0;
  highs_info.num_dual_infeasibilities = 0;
  highs_info.max_dual_infeasibility = 0;
  highs_info.sum_dual_infeasibilities = 0;

  if (lp.num_row_ > 0) {
    // Assign primal, dual and basis status for rows, checking for
    // infeasibility
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
      double primal_infeasibility = 0;
      double lower = lp.row_lower_[iRow];
      double upper = lp.row_upper_[iRow];
      if (lower > primal_feasibility_tolerance) {
        // Lower bound too large for zero activity
        primal_infeasibility = lower;
      } else if (upper < -primal_feasibility_tolerance) {
        // Upper bound too small for zero activity
        primal_infeasibility = -upper;
      }
      solution.row_value.push_back(0);
      solution.row_dual.push_back(0);
      basis.row_status.push_back(HighsBasisStatus::kBasic);
      if (primal_infeasibility > primal_feasibility_tolerance)
        highs_info.num_primal_infeasibilities++;
      highs_info.sum_primal_infeasibilities += primal_infeasibility;
      highs_info.max_primal_infeasibility =
          std::max(primal_infeasibility, highs_info.max_primal_infeasibility);
    }
  }

  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double cost = lp.col_cost_[iCol];
    double dual = (HighsInt)lp.sense_ * cost;
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    double value;
    double primal_infeasibility = 0;
    double dual_infeasibility = -1;
    HighsBasisStatus status = HighsBasisStatus::kNonbasic;
    if (lower > upper) {
      // Inconsistent bounds, so set the variable to lower bound,
      // unless it's infinite. Otherwise set the variable to upper
      // bound, unless it's infinite. Otherwise set the variable to
      // zero.
      if (highs_isInfinity(lower)) {
        // Lower bound of +inf
        if (highs_isInfinity(-upper)) {
          // Upper bound of -inf
          value = 0;
          status = HighsBasisStatus::kZero;
          primal_infeasibility = kHighsInf;
          dual_infeasibility = std::fabs(dual);
        } else {
          // Finite upper bound - since lower exceeds it
          value = upper;
          status = HighsBasisStatus::kUpper;
          primal_infeasibility = lower - value;
          dual_infeasibility = std::max(dual, 0.);
        }
      } else {
        // Finite lower bound
        value = lower;
        status = HighsBasisStatus::kLower;
        primal_infeasibility = value - upper;
        dual_infeasibility = std::max(-dual, 0.);
      }
    } else if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free column: set to zero and record dual infeasibility
      value = 0;
      status = HighsBasisStatus::kZero;
      dual_infeasibility = std::fabs(dual);
    } else if (dual >= dual_feasibility_tolerance) {
      // Column with sufficiently positive dual
      if (!highs_isInfinity(-lower)) {
        // Set to this finite lower bound
        value = lower;
        status = HighsBasisStatus::kLower;
        dual_infeasibility = 0;
      } else {
        // Infinite lower bound so set to upper bound and record dual
        // infeasibility
        value = upper;
        status = HighsBasisStatus::kUpper;
        dual_infeasibility = dual;
      }
    } else if (dual <= -dual_feasibility_tolerance) {
      // Column with sufficiently negative dual
      if (!highs_isInfinity(upper)) {
        // Set to this finite upper bound
        value = upper;
        status = HighsBasisStatus::kUpper;
        dual_infeasibility = 0;
      } else {
        // Infinite upper bound so set to lower bound and record dual
        // infeasibility
        value = lower;
        status = HighsBasisStatus::kLower;
        dual_infeasibility = -dual;
      }
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
      dual_infeasibility = std::fabs(dual);
    }
    assert(status != HighsBasisStatus::kNonbasic);
    assert(dual_infeasibility >= 0);
    solution.col_value[iCol] = value;
    solution.col_dual[iCol] = (HighsInt)lp.sense_ * dual;
    basis.col_status[iCol] = status;
    objective += value * cost;
    if (primal_infeasibility > primal_feasibility_tolerance)
      highs_info.num_primal_infeasibilities++;
    highs_info.sum_primal_infeasibilities += primal_infeasibility;
    highs_info.max_primal_infeasibility =
        std::max(primal_infeasibility, highs_info.max_primal_infeasibility);
    if (dual_infeasibility > dual_feasibility_tolerance)
      highs_info.num_dual_infeasibilities++;
    highs_info.sum_dual_infeasibilities += dual_infeasibility;
    highs_info.max_dual_infeasibility =
        std::max(dual_infeasibility, highs_info.max_dual_infeasibility);
  }
  highs_info.objective_function_value = objective;
  solution.value_valid = true;
  solution.dual_valid = true;
  basis.valid = true;
  basis.useful = true;
  highs_info.basis_validity = kBasisValidityValid;
  setSolutionStatus(highs_info);
  if (highs_info.num_primal_infeasibilities) {
    // Primal infeasible
    model_status = HighsModelStatus::kInfeasible;
  } else if (highs_info.num_dual_infeasibilities) {
    // Dual infeasible => primal unbounded for unconstrained LP
    model_status = HighsModelStatus::kUnbounded;
  } else {
    model_status = HighsModelStatus::kOptimal;
  }

  return HighsStatus::kOk;
}

// Assuming that any user scaling in user_scale_data has been applied,
// determine the model coefficient ranges, assess it for values
// outside the [small, large] range, and give appropriate scaling
// recommendations
void assessExcessiveObjectiveBoundScaling(const HighsLogOptions log_options,
                                          const HighsModel& model,
                                          HighsUserScaleData& user_scale_data) {
  const HighsLp& lp = model.lp_;
  if (lp.num_col_ == 0 || lp.num_row_ == 0) return;
  const bool user_cost_or_bound_scale =
      user_scale_data.user_objective_scale || user_scale_data.user_bound_scale;
  const double small_objective_coefficient =
      kExcessivelySmallObjectiveCoefficient;
  const double large_objective_coefficient =
      kExcessivelyLargeObjectiveCoefficient;
  const double small_bound = kExcessivelySmallBoundValue;
  const double large_bound = kExcessivelyLargeBoundValue;
  std::stringstream message;
  if (user_cost_or_bound_scale) {
    if (user_scale_data.user_objective_scale)
      message << highsFormatToString(" user_cost_scale option value of %d",
                                     user_scale_data.user_objective_scale);
    if (user_scale_data.user_bound_scale) {
      if (user_scale_data.user_objective_scale) message << " and";
      message << highsFormatToString(" user_bound_scale option value of %d",
                                     user_scale_data.user_bound_scale);
    }
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Assessing costs and bounds after applying%s\n",
                 message.str().c_str());
  }
  // Lambda for assessing a finite nonzero
  auto assessFiniteNonzero = [&](const double value, double& min_value,
                                 double& max_value) {
    double abs_value = std::abs(value);
    if (abs_value > 0 && abs_value < kHighsInf) {
      min_value = std::min(abs_value, min_value);
      max_value = std::max(abs_value, max_value);
    }
  };
  double min_continuous_col_cost = kHighsInf;
  double min_noncontinuous_col_cost = kHighsInf;
  double max_continuous_col_cost = -kHighsInf;
  double max_noncontinuous_col_cost = -kHighsInf;
  double min_continuous_col_bound = kHighsInf;
  double min_noncontinuous_col_bound = kHighsInf;
  double max_continuous_col_bound = -kHighsInf;
  double max_noncontinuous_col_bound = -kHighsInf;
  double min_continuous_matrix_value = kHighsInf;
  double min_noncontinuous_matrix_value = kHighsInf;
  double max_continuous_matrix_value = -kHighsInf;
  double max_noncontinuous_matrix_value = -kHighsInf;
  const bool is_mip = lp.integrality_.size();
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (is_mip && lp.integrality_[iCol] != HighsVarType::kContinuous) {
      assessFiniteNonzero(lp.col_cost_[iCol], min_noncontinuous_col_cost,
                          max_noncontinuous_col_cost);
      assessFiniteNonzero(lp.col_lower_[iCol], min_noncontinuous_col_bound,
                          max_noncontinuous_col_bound);
      assessFiniteNonzero(lp.col_upper_[iCol], min_noncontinuous_col_bound,
                          max_noncontinuous_col_bound);
    } else {
      assessFiniteNonzero(lp.col_cost_[iCol], min_continuous_col_cost,
                          max_continuous_col_cost);
      assessFiniteNonzero(lp.col_lower_[iCol], min_continuous_col_bound,
                          max_continuous_col_bound);
      assessFiniteNonzero(lp.col_upper_[iCol], min_continuous_col_bound,
                          max_continuous_col_bound);
    }
  }
  double min_col_cost =
      std::min(min_continuous_col_cost, min_noncontinuous_col_cost);
  double max_col_cost =
      std::max(max_continuous_col_cost, max_noncontinuous_col_cost);
  double min_col_bound =
      std::min(min_continuous_col_bound, min_noncontinuous_col_bound);
  double max_col_bound =
      std::max(max_continuous_col_bound, max_noncontinuous_col_bound);

  double min_matrix_value = kHighsInf;
  double max_matrix_value = -kHighsInf;
  const HighsInt num_matrix_nz = lp.a_matrix_.numNz();
  for (HighsInt iEl = 0; iEl < num_matrix_nz; iEl++)
    assessFiniteNonzero(lp.a_matrix_.value_[iEl], min_matrix_value,
                        max_matrix_value);

  double min_row_bound = kHighsInf;
  double max_row_bound = -kHighsInf;
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    assessFiniteNonzero(lp.row_lower_[iRow], min_row_bound, max_row_bound);
    assessFiniteNonzero(lp.row_upper_[iRow], min_row_bound, max_row_bound);
  }

  double min_continuous_hessian_value = kHighsInf;
  double max_continuous_hessian_value = -kHighsInf;
  const HighsInt num_hessian_nz = model.hessian_.numNz();
  for (HighsInt iEl = 0; iEl < num_hessian_nz; iEl++)
    assessFiniteNonzero(model.hessian_.value_[iEl],
                        min_continuous_hessian_value,
                        max_continuous_hessian_value);

  // Determine the minimum and maximum overall bounds that can be
  // scaled with user_bound_scale before zeroing extrema due to
  // absence of finite nonzero bounds

  double min_scalable_bound = std::min(min_continuous_col_bound, min_row_bound);
  double max_scalable_bound = std::max(max_continuous_col_bound, max_row_bound);
  if (min_scalable_bound == kHighsInf) min_scalable_bound = 0;
  if (max_scalable_bound == -kHighsInf) max_scalable_bound = 0;

  if (min_col_cost == kHighsInf) min_col_cost = 0;
  if (max_col_cost == -kHighsInf) max_col_cost = 0;
  if (min_col_bound == kHighsInf) min_col_bound = 0;
  if (max_col_bound == -kHighsInf) max_col_bound = 0;
  if (min_row_bound == kHighsInf) min_row_bound = 0;
  if (max_row_bound == -kHighsInf) max_row_bound = 0;

  double min_hessian_value = min_continuous_hessian_value;
  double max_hessian_value = max_continuous_hessian_value;
  if (min_hessian_value == kHighsInf) min_hessian_value = 0;
  if (max_hessian_value == -kHighsInf) max_hessian_value = 0;

  // Report on the coefficient ranges
  highsLogUser(log_options, HighsLogType::kInfo, "Coefficient ranges:\n");
  if (num_matrix_nz)
    highsLogUser(log_options, HighsLogType::kInfo, "  Matrix  [%5.0e, %5.0e]\n",
                 min_matrix_value, max_matrix_value);
  if (lp.num_col_) {
    highsLogUser(log_options, HighsLogType::kInfo, "  Cost    [%5.0e, %5.0e]\n",
                 min_col_cost, max_col_cost);
    if (num_hessian_nz)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "  Hessian [%5.0e, %5.0e]\n", min_hessian_value,
                   max_hessian_value);
    highsLogUser(log_options, HighsLogType::kInfo, "  Bound   [%5.0e, %5.0e]\n",
                 min_col_bound, max_col_bound);
  }
  if (lp.num_row_)
    highsLogUser(log_options, HighsLogType::kInfo, "  RHS     [%5.0e, %5.0e]\n",
                 min_row_bound, max_row_bound);

  // LPs with no columns or no finite nonzero costs will have
  // max_col_cost = 0
  assert(max_col_cost >= 0);
  // LPs with no columns or no finite nonzero bounds will have
  // max_col_bound = 0
  assert(max_col_bound >= 0);
  // LPs with no rows or no finite nonzero bounds will have
  // max_row_bound = 0
  assert(max_row_bound >= 0);

  const std::string problem =
      user_cost_or_bound_scale ? "User-scaled problem" : "Problem";

  if (0 < min_col_cost && min_col_cost < small_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small costs\n", problem.c_str());
  if (max_col_cost > large_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large costs\n", problem.c_str());
  if (0 < min_hessian_value && min_hessian_value < small_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small Hessian values\n",
                 problem.c_str());
  if (max_hessian_value > large_objective_coefficient)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large Hessian values\n",
                 problem.c_str());
  if (0 < min_col_bound && min_col_bound < small_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small column bounds\n",
                 problem.c_str());
  if (max_col_bound > large_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large column bounds\n",
                 problem.c_str());
  if (0 < min_row_bound && min_row_bound < small_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively small row bounds\n", problem.c_str());
  if (max_row_bound > large_bound)
    highsLogUser(log_options, HighsLogType::kWarning,
                 "%s has some excessively large row bounds\n", problem.c_str());

  // Determine recommended user scaling values
  auto suggestScaling = [&](double min_value, double max_value,
                            double small_value, double large_value) {
    double ratio = 1;
    if (max_value > large_value) {
      // Max scalable value is large, so suggest scaling values down
      // so that the max value is large_value
      ratio = large_value / max_value;
      assert(0 == 11);
    } else if (0 < max_value && max_value < small_value) {
      // All scalable values are small, so suggest scaling them up so
      // the max value is small_value
      ratio = small_value / max_value;
    }
    assert(ratio);
    return ratio;
  };

  double suggested_bound_scaling = suggestScaling(
      min_scalable_bound, max_scalable_bound, small_bound, large_bound);
  // Determine the suggested (new) value for user_bound_scale,
  // allowing for the fact that the current value has been applied
  HighsInt dl_user_bound_scale = std::ceil(std::log2(suggested_bound_scaling));
  user_scale_data.suggested_user_bound_scale =
      user_scale_data.user_bound_scale + dl_user_bound_scale;
  // Determine the order of magnitude of the suggested bound scaling -
  // just for logging
  HighsInt suggested_bound_scale_order_of_magnitude =
      std::ceil(std::log10(suggested_bound_scaling));
  // Applying the suggested bound scaling requires the costs and
  // matrix columns of non-continuous variables to be scaled, and any
  // Hessian entries are also scaled
  //
  // Determine the corresponding extreme non-continuous costs and
  // update the extreme costs so that objective scaling can be
  // suggested
  double suggested_user_bound_scale_value =
      pow(2.0, user_scale_data.suggested_user_bound_scale);
  min_noncontinuous_col_cost *= suggested_user_bound_scale_value;
  max_noncontinuous_col_cost *= suggested_user_bound_scale_value;
  min_hessian_value /= suggested_user_bound_scale_value;
  max_hessian_value /= suggested_user_bound_scale_value;

  min_col_cost = std::min(min_continuous_col_cost, min_noncontinuous_col_cost);
  max_col_cost = std::max(max_continuous_col_cost, max_noncontinuous_col_cost);
  double min_objective_coefficient =
      std::min(min_col_cost, min_continuous_hessian_value);
  double max_objective_coefficient =
      std::max(max_col_cost, max_continuous_hessian_value);
  if (min_objective_coefficient == kHighsInf) min_objective_coefficient = 0;
  if (max_objective_coefficient == -kHighsInf) max_objective_coefficient = 0;

  double suggested_objective_scaling =
      suggestScaling(min_objective_coefficient, max_objective_coefficient,
                     small_objective_coefficient, large_objective_coefficient);
  // Determine the suggested (new) value for user_objective_scale,
  // allowing for the fact that the current value has been applied
  HighsInt dl_user_objective_scale =
      std::ceil(std::log2(suggested_objective_scaling));
  user_scale_data.suggested_user_objective_scale =
      user_scale_data.user_objective_scale + dl_user_objective_scale;
  // Determine the order of magnitude of the suggested objective scaling -
  // just for logging
  HighsInt suggested_objective_scale_order_of_magnitude =
      std::ceil(std::log10(suggested_objective_scaling));

  // Only report the order of magnitude scaling if there is no user
  // scaling
  bool order_of_magnitude_message =
      suggested_objective_scale_order_of_magnitude &&
      !user_scale_data.user_objective_scale;
  message.str(std::string());
  if (order_of_magnitude_message)
    message << highsFormatToString(
        "   Consider scaling the objective by 1e%+1d",
        int(suggested_objective_scale_order_of_magnitude));
  if (dl_user_objective_scale) {
    if (!order_of_magnitude_message) {
      message << "   Consider";
    } else {
      message << ", or";
    }
    message << highsFormatToString(
        " setting the user_cost_scale option to %d",
        int(user_scale_data.suggested_user_objective_scale));
  }
  if (order_of_magnitude_message || dl_user_objective_scale)
    highsLogUser(log_options, HighsLogType::kWarning, "%s\n",
                 message.str().c_str());

  message.str(std::string());
  order_of_magnitude_message = suggested_bound_scale_order_of_magnitude &&
                               !user_scale_data.user_bound_scale;
  message.str(std::string());
  if (order_of_magnitude_message)
    message << highsFormatToString(
        "   Consider scaling the    bounds by 1e%+1d",
        int(suggested_bound_scale_order_of_magnitude));
  if (dl_user_bound_scale) {
    if (!order_of_magnitude_message) {
      message << "   Consider";
    } else {
      message << ", or";
    }
    message << highsFormatToString(
        " setting the user_bound_scale option to %d",
        int(user_scale_data.suggested_user_bound_scale));
  }
  if (order_of_magnitude_message || dl_user_bound_scale)
    highsLogUser(log_options, HighsLogType::kWarning, "%s\n",
                 message.str().c_str());
}
