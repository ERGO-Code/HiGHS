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
/**@file lp_data/HighsSolutionDebug.cpp
 * @brief
 */
#include "lp_data/HighsSolutionDebug.h"

#include <math.h>

#include <vector>

#include "lp_data/HighsDebug.h"
#include "lp_data/HighsModelUtils.h"
#include "util/HighsUtils.h"

const double large_relative_solution_param_error = 1e-12;
const double excessive_relative_solution_param_error =
    sqrt(large_relative_solution_param_error);

const double large_residual_error = 1e-12;
const double excessive_residual_error = sqrt(large_residual_error);

HighsDebugStatus debugHighsSolution(const string message,
                                    const HighsOptions& options,
                                    const HighsLp& lp,
                                    const HighsSolution& solution,
                                    const HighsBasis& basis) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to report on KKT errors after solving a model when only
  // the solution (possibly only primal) and (possibly) basis are
  // known
  //
  // Set up a HighsModelStatus and HighsSolutionParams just to
  // complete the parameter list.By setting
  // check_model_status_and_solution_params to be false they waren't
  // used.
  HighsModelStatus dummy_model_status;
  HighsSolutionParams dummy_solution_params;
  // Call resetModelStatusAndSolutionParams to side-step compiler
  // warning.
  resetModelStatusAndSolutionParams(dummy_model_status, dummy_solution_params,
                                    options);
  const bool check_model_status_and_solution_params = false;
  return debugHighsSolution(message, options, lp, solution, basis,
                            dummy_model_status, dummy_solution_params,
                            check_model_status_and_solution_params);
}

HighsDebugStatus debugHighsSolution(const std::string message,
                                    const HighsModelObject& model) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to check the unscaled model status and solution params
  const bool check_model_status_and_solution_params = true;
  return debugHighsSolution(message, model.options_, model.lp_, model.solution_,
                            model.basis_, model.unscaled_model_status_,
                            model.solution_params_,
                            check_model_status_and_solution_params);
}

HighsDebugStatus debugHighsSolution(
    const string message, const HighsOptions& options, const HighsLp& lp,
    const HighsSolution& solution, const HighsBasis& basis,
    const HighsModelStatus model_status, const HighsInfo& info) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to check the HiGHS model_status and info
  //
  // Copy the data from info to solution_params so general method can be used
  //
  HighsSolutionParams solution_params;
  copyFromInfo(solution_params, info);
  solution_params.primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  solution_params.dual_feasibility_tolerance =
      options.dual_feasibility_tolerance;
  const bool check_model_status_and_solution_params = true;
  return debugHighsSolution(message, options, lp, solution, basis, model_status,
                            solution_params,
                            check_model_status_and_solution_params);
}

HighsDebugStatus debugHighsSolution(
    const std::string message, const HighsOptions& options, const HighsLp& lp,
    const HighsSolution& solution, const HighsBasis& basis,
    const HighsModelStatus model_status,
    const HighsSolutionParams& solution_params,
    const bool check_model_status_and_solution_params) {
  // Non-trivially expensive analysis of a solution to a model
  //
  // Called to possibly check the model_status and solution_params,
  // and then report on the KKT and model status, plus any errors
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status;
  // Use local_model_status to for checking - or if it's not known
  HighsModelStatus local_model_status = HighsModelStatus::kNotset;
  // Use local_solution_params to determine solution_params for
  // checking - or if it's not known
  HighsSolutionParams local_solution_params;
  local_solution_params.primal_feasibility_tolerance =
      options.primal_feasibility_tolerance;
  local_solution_params.dual_feasibility_tolerance =
      options.dual_feasibility_tolerance;
  if (check_model_status_and_solution_params) {
    double local_objective_function_value = 0;
    if (solution.value_valid)
      local_objective_function_value = computeObjectiveValue(lp, solution);
    local_solution_params.objective_function_value =
        local_objective_function_value;
  }
  HighsPrimalDualErrors primal_dual_errors;
  // Determine the extent to which KKT conditions are not satisfied,
  // accumulating data on primal/dual errors relating to any basis
  // implications and excessive residuals
  const bool get_residuals =
      true;  // options.highs_debug_level >= kHighsDebugLevelCostly;
  getKktFailures(lp, solution, basis, local_solution_params, primal_dual_errors,
                 get_residuals);
  HighsInt& num_primal_infeasibility =
      local_solution_params.num_primal_infeasibility;
  HighsInt& num_dual_infeasibility =
      local_solution_params.num_dual_infeasibility;
  if (check_model_status_and_solution_params) {
    // Can assume that model_status and solution_params are known, so should be
    // checked
    local_model_status = model_status;
    // Check that solution_params is the same as when computed from scratch
    return_status = debugCompareSolutionParams(options, solution_params,
                                               local_solution_params);
    if (return_status != HighsDebugStatus::kOk) return return_status;
  } else {
    // Determine whether optimality can be reported
    if (num_primal_infeasibility == 0 && num_dual_infeasibility == 0)
      local_model_status = HighsModelStatus::kOptimal;
  }
  if (check_model_status_and_solution_params &&
      model_status == HighsModelStatus::kOptimal) {
    bool error_found = false;
    if (num_primal_infeasibility > 0) {
      error_found = true;
      highsLogDev(options.log_options, HighsLogType::kError,
                  "debugHighsSolution: %" HIGHSINT_FORMAT
                  " primal infeasiblilities but model status is %s\n",
                  num_primal_infeasibility,
                  utilModelStatusToString(model_status).c_str());
    }
    if (num_dual_infeasibility > 0) {
      error_found = true;
      highsLogDev(options.log_options, HighsLogType::kError,
                  "debugHighsSolution: %" HIGHSINT_FORMAT
                  " dual infeasiblilities but model status is %s\n",
                  num_dual_infeasibility,
                  utilModelStatusToString(model_status).c_str());
    }
    if (error_found) return HighsDebugStatus::kLogicalError;
  }
  // Report on the solution
  debugReportHighsSolution(message, options.log_options, local_solution_params,
                           local_model_status);
  // Analyse the primal and dual errors
  return_status = debugAnalysePrimalDualErrors(options, primal_dual_errors);
  return return_status;
}

void debugReportHighsSolution(const string message,
                              const HighsLogOptions& log_options,
                              const HighsSolutionParams& solution_params,
                              const HighsModelStatus model_status) {
  highsLogDev(log_options, HighsLogType::kInfo, "\nHiGHS solution: %s\n",
              message.c_str());
  highsLogDev(log_options, HighsLogType::kInfo,
              "Infeas:                Pr %" HIGHSINT_FORMAT
              "(Max %.4g, Sum %.4g); Du %" HIGHSINT_FORMAT
              "(Max %.4g, "
              "Sum %.4g); Status: %s\n",
              solution_params.num_primal_infeasibility,
              solution_params.max_primal_infeasibility,
              solution_params.sum_primal_infeasibility,
              solution_params.num_dual_infeasibility,
              solution_params.max_dual_infeasibility,
              solution_params.sum_dual_infeasibility,
              utilModelStatusToString(model_status).c_str());
}

HighsDebugStatus debugBasisConsistent(const HighsOptions& options,
                                      const HighsLp lp,
                                      const HighsBasis& basis) {
  // Cheap analysis of a HiGHS basis, checking vector sizes, numbers
  // of basic/nonbasic variables
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (!basis.valid) return return_status;
  bool consistent = isBasisConsistent(lp, basis);
  if (!consistent) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS basis inconsistency\n");
    assert(consistent);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugBasisRightSize(const HighsOptions& options,
                                     const HighsLp lp,
                                     const HighsBasis& basis) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool right_size = isBasisRightSize(lp, basis);
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS basis size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugPrimalSolutionRightSize(const HighsOptions& options,
                                              const HighsLp lp,
                                              const HighsSolution& solution) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool right_size = isPrimalSolutionRightSize(lp, solution);
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS primal solution size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugDualSolutionRightSize(const HighsOptions& options,
                                            const HighsLp lp,
                                            const HighsSolution& solution) {
  if (options.highs_debug_level < kHighsDebugLevelCheap)
    return HighsDebugStatus::kNotChecked;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  bool right_size = isDualSolutionRightSize(lp, solution);
  if (!right_size) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "HiGHS dual solution size error\n");
    assert(right_size);
    return_status = HighsDebugStatus::kLogicalError;
  }
  return return_status;
}

HighsDebugStatus debugSolutionRightSize(const HighsOptions& options,
                                        const HighsLp lp,
                                        const HighsSolution& solution) {
  HighsDebugStatus return_status;
  return_status = debugPrimalSolutionRightSize(options, lp, solution);
  if (return_status != HighsDebugStatus::kOk) return return_status;
  return debugDualSolutionRightSize(options, lp, solution);
}

// Methods below are not called externally

HighsDebugStatus debugHaveBasisAndSolutionData(const HighsLp& lp,
                                               const HighsBasis& basis,
                                               const HighsSolution& solution) {
  if (!isSolutionRightSize(lp, solution))
    return HighsDebugStatus::kLogicalError;
  if (!isBasisRightSize(lp, basis) && basis.valid)
    return HighsDebugStatus::kLogicalError;
  return HighsDebugStatus::kOk;
}

HighsDebugStatus debugAnalysePrimalDualErrors(
    const HighsOptions& options, HighsPrimalDualErrors& primal_dual_errors) {
  std::string value_adjective;
  HighsLogType report_level;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  const bool force_report = options.highs_debug_level >= kHighsDebugLevelCostly;
  if (primal_dual_errors.num_nonzero_basic_duals) {
    value_adjective = "Error";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kLogicalError;
  } else {
    value_adjective = "";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  if (force_report) report_level = HighsLogType::kInfo;
  highsLogDev(
      options.log_options, report_level,
      "PrDuErrors : %-9s Nonzero basic duals:       num = %2" HIGHSINT_FORMAT
      "; "
      "max = %9.4g; sum = %9.4g\n",
      value_adjective.c_str(), primal_dual_errors.num_nonzero_basic_duals,
      primal_dual_errors.max_nonzero_basic_dual,
      primal_dual_errors.sum_nonzero_basic_duals);

  if (primal_dual_errors.num_off_bound_nonbasic) {
    value_adjective = "Error";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kLogicalError;
  } else {
    value_adjective = "";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  if (force_report) report_level = HighsLogType::kInfo;
  highsLogDev(
      options.log_options, report_level,
      "PrDuErrors : %-9s Off-bound nonbasic values: num = %2" HIGHSINT_FORMAT
      "; "
      "max = %9.4g; sum = %9.4g\n",
      value_adjective.c_str(), primal_dual_errors.num_off_bound_nonbasic,
      primal_dual_errors.max_off_bound_nonbasic,
      primal_dual_errors.sum_off_bound_nonbasic);

  if (primal_dual_errors.max_primal_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (primal_dual_errors.max_primal_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  if (force_report) report_level = HighsLogType::kInfo;
  highsLogDev(
      options.log_options, report_level,
      "PrDuErrors : %-9s Primal residual:           num = %2" HIGHSINT_FORMAT
      "; "
      "max = %9.4g; sum = %9.4g\n",
      value_adjective.c_str(), primal_dual_errors.num_primal_residual,
      primal_dual_errors.max_primal_residual,
      primal_dual_errors.sum_primal_residual);

  if (primal_dual_errors.max_dual_residual > excessive_residual_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (primal_dual_errors.max_dual_residual > large_residual_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "";
    report_level = HighsLogType::kVerbose;
    return_status = HighsDebugStatus::kOk;
  }
  if (force_report) report_level = HighsLogType::kInfo;
  highsLogDev(
      options.log_options, report_level,
      "PrDuErrors : %-9s Dual residual:             num = %2" HIGHSINT_FORMAT
      "; "
      "max = %9.4g; sum = %9.4g\n",
      value_adjective.c_str(), primal_dual_errors.num_dual_residual,
      primal_dual_errors.max_dual_residual,
      primal_dual_errors.sum_dual_residual);

  return return_status;
}

HighsDebugStatus debugCompareSolutionParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  return_status =
      debugWorseStatus(debugCompareSolutionObjectiveParams(
                           options, solution_params0, solution_params1),
                       return_status);
  return_status =
      debugWorseStatus(debugCompareSolutionStatusParams(
                           options, solution_params0, solution_params1),
                       return_status);
  return_status =
      debugWorseStatus(debugCompareSolutionInfeasibilityParams(
                           options, solution_params0, solution_params1),
                       return_status);
  return return_status;
}

HighsDebugStatus debugCompareSolutionObjectiveParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  return debugCompareSolutionParamValue(
      "objective_function_value", options,
      solution_params0.objective_function_value,
      solution_params1.objective_function_value);
}

HighsDebugStatus debugCompareSolutionStatusParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  return_status = debugWorseStatus(
      debugCompareSolutionParamInteger("primal_status", options,
                                       solution_params0.primal_solution_status,
                                       solution_params1.primal_solution_status),
      return_status);
  return_status =
      debugWorseStatus(debugCompareSolutionParamInteger(
                           "dual_status", options, solution_params0.dual_solution_status,
                           solution_params1.dual_solution_status),
                       return_status);
  return return_status;
}

HighsDebugStatus debugCompareSolutionInfeasibilityParams(
    const HighsOptions& options, const HighsSolutionParams& solution_params0,
    const HighsSolutionParams& solution_params1) {
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  return_status =
      debugWorseStatus(debugCompareSolutionParamInteger(
                           "num_primal_infeasibility", options,
                           solution_params0.num_primal_infeasibility,
                           solution_params1.num_primal_infeasibility),
                       return_status);
  return_status = debugWorseStatus(
      debugCompareSolutionParamValue("sum_primal_infeasibility", options,
                                     solution_params0.sum_primal_infeasibility,
                                     solution_params1.sum_primal_infeasibility),
      return_status);
  return_status = debugWorseStatus(
      debugCompareSolutionParamValue("max_primal_infeasibility", options,
                                     solution_params0.max_primal_infeasibility,
                                     solution_params1.max_primal_infeasibility),
      return_status);

  return_status = debugWorseStatus(
      debugCompareSolutionParamInteger("num_dual_infeasibility", options,
                                       solution_params0.num_dual_infeasibility,
                                       solution_params1.num_dual_infeasibility),
      return_status);
  return_status = debugWorseStatus(
      debugCompareSolutionParamValue("sum_dual_infeasibility", options,
                                     solution_params0.sum_dual_infeasibility,
                                     solution_params1.sum_dual_infeasibility),
      return_status);
  return_status = debugWorseStatus(
      debugCompareSolutionParamValue("max_dual_infeasibility", options,
                                     solution_params0.max_dual_infeasibility,
                                     solution_params1.max_dual_infeasibility),
      return_status);
  return return_status;
}

HighsDebugStatus debugCompareSolutionParamValue(const string name,
                                                const HighsOptions& options,
                                                const double v0,
                                                const double v1) {
  if (v0 == v1) return HighsDebugStatus::kOk;
  double delta = highsRelativeDifference(v0, v1);
  std::string value_adjective;
  HighsLogType report_level;
  HighsDebugStatus return_status = HighsDebugStatus::kOk;
  if (delta > excessive_relative_solution_param_error) {
    value_adjective = "Excessive";
    report_level = HighsLogType::kError;
    return_status = HighsDebugStatus::kError;
  } else if (delta > large_relative_solution_param_error) {
    value_adjective = "Large";
    report_level = HighsLogType::kDetailed;
    return_status = HighsDebugStatus::kWarning;
  } else {
    value_adjective = "OK";
    report_level = HighsLogType::kVerbose;
  }
  highsLogDev(options.log_options, report_level,
              "SolutionPar:  %-9s relative difference of %9.4g for %s\n",
              value_adjective.c_str(), delta, name.c_str());
  return return_status;
}

HighsDebugStatus debugCompareSolutionParamInteger(const string name,
                                                  const HighsOptions& options,
                                                  const HighsInt v0,
                                                  const HighsInt v1) {
  if (v0 == v1) return HighsDebugStatus::kOk;
  highsLogDev(options.log_options, HighsLogType::kError,
              "SolutionPar:  difference of %" HIGHSINT_FORMAT " for %s\n",
              v1 - v0, name.c_str());
  return HighsDebugStatus::kLogicalError;
}

void debugReportHighsBasicSolution(const string message,
                                   const HighsOptions& options,
                                   const HighsSolutionParams& solution_params,
                                   const HighsModelStatus model_status) {
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "\nHiGHS basic solution: %s\n", message.c_str());
  highsLogDev(options.log_options, HighsLogType::kInfo,
              "Infeas:                Pr %" HIGHSINT_FORMAT
              "(Max %.4g, Sum %.4g); Du %" HIGHSINT_FORMAT
              "(Max %.4g, "
              "Sum %.4g); Status: %s\n",
              solution_params.num_primal_infeasibility,
              solution_params.max_primal_infeasibility,
              solution_params.sum_primal_infeasibility,
              solution_params.num_dual_infeasibility,
              solution_params.max_dual_infeasibility,
              solution_params.sum_dual_infeasibility,
              utilModelStatusToString(model_status).c_str());
}
