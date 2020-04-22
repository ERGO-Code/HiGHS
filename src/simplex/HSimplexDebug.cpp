/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HSimplexDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "simplex/HSimplexDebug.h"

#include "simplex/HSimplex.h"
#include "simplex/SimplexTimer.h"

const double computed_primal_excessive_absolute_norm = 1e12;
const double computed_primal_large_absolute_norm =
    sqrt(computed_primal_excessive_absolute_norm);
const double computed_primal_excessive_relative_norm = 1e6;
const double computed_primal_large_relative_norm =
    sqrt(computed_primal_excessive_relative_norm);
const double computed_dual_excessive_absolute_norm = 1e12;
const double computed_dual_large_absolute_norm =
    sqrt(computed_dual_excessive_absolute_norm);
const double computed_dual_excessive_relative_norm = 1e6;
const double computed_dual_large_relative_norm =
    sqrt(computed_dual_excessive_relative_norm);
const double computed_dual_small_relative_basic_dual_norm = 1e-12;
const double computed_dual_large_relative_basic_dual_norm =
    sqrt(computed_dual_small_relative_basic_dual_norm);
const double computed_dual_small_relative_change = 1e-12;
const double computed_dual_large_relative_change =
    sqrt(computed_dual_small_relative_change);
const double computed_dual_small_absolute_change = 1e-6;
const double computed_dual_large_absolute_change =
    sqrt(computed_dual_small_absolute_change);
const double updated_objective_small_relative_error = 1e-12;
const double updated_objective_large_relative_error =
    sqrt(updated_objective_small_relative_error);
const double updated_objective_small_absolute_error = 1e-6;
const double updated_objective_large_absolute_error =
    sqrt(updated_objective_small_absolute_error);
const double excessive_basis_condition = 1e16;
const double large_basis_condition = sqrt(excessive_basis_condition);
const double fair_basis_condition = sqrt(large_basis_condition);

HighsDebugStatus debugComputePrimal(const HighsModelObject& highs_model_object,
                                    const std::vector<double>& primal_rhs) {
  // Non-trivially expensive analysis of computed primal values.
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;
  const std::vector<double>& primal_value =
      highs_model_object.simplex_info_.baseValue_;

  int num_row = highs_model_object.simplex_lp_.numRow_;

  // Use the size of the RHS to determine whether to use it
  const bool have_primal_rhs = (int)primal_rhs.size() == num_row;

  double primal_rhs_norm = 0;
  if (have_primal_rhs) {
    for (int iRow = 0; iRow < num_row; iRow++)
      primal_rhs_norm += fabs(primal_rhs[iRow]);
  }
  double computed_primal_absolute_norm = 0;
  for (int iRow = 0; iRow < num_row; iRow++)
    computed_primal_absolute_norm += fabs(primal_value[iRow]);

  std::string value_adjective;
  int report_level;
  return_status = HighsDebugStatus::OK;
  double computed_primal_relative_norm;
  if (primal_rhs_norm) {
    computed_primal_relative_norm =
        computed_primal_absolute_norm / primal_rhs_norm;
  } else {
    computed_primal_relative_norm = -1;
  }
  if (computed_primal_relative_norm > computed_primal_excessive_relative_norm ||
      computed_primal_absolute_norm > computed_primal_excessive_absolute_norm) {
    value_adjective = "Excessive";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else if (computed_primal_relative_norm >
                 computed_primal_large_relative_norm ||
             computed_primal_absolute_norm >
                 computed_primal_large_absolute_norm) {
    value_adjective = "Large";
    report_level = ML_DETAILED;
    return_status = HighsDebugStatus::WARNING;
  } else {
    value_adjective = "SMALL";
    report_level = ML_VERBOSE;
  }
  HighsPrintMessage(highs_model_object.options_.output,
                    highs_model_object.options_.message_level, report_level,
                    "ComputePrimal: %s absolute (%g) or relative (%g) norm of "
                    "primal values\n",
                    value_adjective.c_str(), computed_primal_absolute_norm,
                    computed_primal_relative_norm);
  if (have_primal_rhs && !primal_rhs_norm) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "ComputePrimal: |PrimalRHS| = %g", primal_rhs_norm);
    return_status = HighsDebugStatus::WARNING;
  }
  return return_status;
}
HighsDebugStatus debugComputedDual(const HighsModelObject& highs_model_object,
                                   const std::vector<double>& previous_dual,
                                   const std::vector<double>& basic_costs,
                                   const std::vector<double>& row_dual) {
  // Non-trivially expensive analysis of computed dual values.
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsDebugStatus return_status = HighsDebugStatus::NOT_CHECKED;
  const std::vector<double>& new_dual =
      highs_model_object.simplex_info_.workDual_;

  int num_row = highs_model_object.simplex_lp_.numRow_;
  int num_col = highs_model_object.simplex_lp_.numCol_;

  // Make sure that previous_dual and basic_costs have the right size
  assert((int)basic_costs.size() == num_row);
  assert((int)row_dual.size() == num_row);

  double basic_costs_norm = 0;
  double row_dual_norm = 0;
  for (int iRow = 0; iRow < num_row; iRow++) {
    basic_costs_norm += fabs(basic_costs[iRow]);
    row_dual_norm += fabs(row_dual[iRow]);
  }
  double basic_dual_norm = 0;
  double nonbasic_dual_norm = 0;
  for (int iVar = 0; iVar < num_row + num_col; iVar++) {
    if (!highs_model_object.simplex_basis_.nonbasicFlag_[iVar]) {
      basic_dual_norm += fabs(new_dual[iVar]);
      continue;
    }
    nonbasic_dual_norm += fabs(new_dual[iVar]);
  }
  std::string value_adjective;
  int report_level;
  return_status = HighsDebugStatus::OK;
  // Comment on the norm of basic duals (relative to |c_B|) which, as
  // c_B-BB^{-1}c_B, should be zero
  if (basic_costs_norm) {
    double computed_dual_relative_basic_dual_norm =
        basic_dual_norm / basic_costs_norm;
    if (computed_dual_relative_basic_dual_norm >
        computed_dual_large_relative_basic_dual_norm) {
      value_adjective = "Large";
      report_level = ML_DETAILED;
      return_status = HighsDebugStatus::WARNING;
    } else if (computed_dual_relative_basic_dual_norm >
               computed_dual_small_relative_basic_dual_norm) {
      value_adjective = "Small";
      report_level = ML_DETAILED;
      return_status = HighsDebugStatus::WARNING;
    } else {
      value_adjective = "OK";
      report_level = ML_VERBOSE;
    }
    HighsPrintMessage(highs_model_object.options_.output,
                      highs_model_object.options_.message_level, report_level,
                      "ComputedDual: %s norm of basic duals (%g) relative to "
                      "basic costs (%g)\n",
                      value_adjective.c_str(),
                      computed_dual_relative_basic_dual_norm, basic_costs_norm);
  } else {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "ComputedDual: basic cost norm is = %g", basic_costs_norm);
    return_status = HighsDebugStatus::WARNING;
  }
  std::string change_adjective;
  double computed_dual_absolute_change = 0;
  double computed_dual_relative_change = 0;
  if (previous_dual.size()) {
    for (int iVar = 0; iVar < num_row + num_col; iVar++) {
      if (!highs_model_object.simplex_basis_.nonbasicFlag_[iVar]) continue;
      computed_dual_absolute_change +=
          fabs(new_dual[iVar] - previous_dual[iVar]);
    }
    if (nonbasic_dual_norm)
      computed_dual_relative_change =
          computed_dual_absolute_change / nonbasic_dual_norm;
  }
  if (nonbasic_dual_norm) {
    if (computed_dual_relative_change > computed_dual_large_relative_change ||
        computed_dual_absolute_change > computed_dual_large_absolute_change) {
      change_adjective = "Large";
      report_level = ML_DETAILED;
      return_status = HighsDebugStatus::WARNING;
    } else if (computed_dual_relative_change >
                   computed_dual_small_relative_change ||
               computed_dual_absolute_change >
                   computed_dual_small_absolute_change) {
      change_adjective = "Small";
      report_level = ML_DETAILED;
      return_status = HighsDebugStatus::WARNING;
    } else {
      change_adjective = "OK";
      report_level = ML_VERBOSE;
    }
    HighsPrintMessage(highs_model_object.options_.output,
                      highs_model_object.options_.message_level, report_level,
                      "ComputedDual: B.pi=c_B has |c_B|=%g; |pi|=%g; "
                      "|pi^TA-c|=[Nonbasic %g; Basic %g]\n",
                      basic_costs_norm, row_dual_norm, nonbasic_dual_norm,
                      basic_dual_norm);
    HighsPrintMessage(
        highs_model_object.options_.output,
        highs_model_object.options_.message_level, report_level,
        "ComputedDual: %s absolute (%g) or relative (%g) change\n",
        change_adjective.c_str(), computed_dual_absolute_change,
        computed_dual_relative_change);
  } else {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "ComputedDual: |NonbasicDual| = %g", nonbasic_dual_norm);
    return_status = HighsDebugStatus::WARNING;
  }
  return return_status;
}

HighsDebugStatus debugUpdatedObjectiveValue(
    HighsModelObject& highs_model_object, const SimplexAlgorithm algorithm,
    const int phase, const std::string message) {
  // Non-trivially expensive check of updated objective value. Computes the
  // exact objective value
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;

  static bool have_previous_primal_objective_value;
  static double previous_primal_objective_value;
  static double previous_updated_primal_objective_value;
  static double updated_primal_objective_correction;

  static bool have_previous_dual_objective_value;
  static double previous_dual_objective_value;
  static double previous_updated_dual_objective_value;
  static double updated_dual_objective_correction;
  if (phase < 0) {
    if (algorithm == SimplexAlgorithm::PRIMAL) {
      have_previous_primal_objective_value = false;
    } else {
      have_previous_dual_objective_value = false;
    }
    return HighsDebugStatus::OK;
  }
  double objective_value;
  double updated_objective_value;
  bool have_previous_objective_value;
  double previous_objective_value;
  double previous_updated_objective_value;
  double updated_objective_correction;
  std::string algorithm_name;
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    algorithm_name = "primal";
    have_previous_objective_value = have_previous_primal_objective_value;
    if (have_previous_objective_value) {
      previous_objective_value = previous_primal_objective_value;
      previous_updated_objective_value =
          previous_updated_primal_objective_value;
      updated_objective_correction = updated_primal_objective_correction;
    }
    updated_objective_value = simplex_info.updated_primal_objective_value;
    // Save the current objective value so that it can be recovered
    // after calling computePrimalObjectiveValue
    double save_objective_value = simplex_info.primal_objective_value;
    computePrimalObjectiveValue(highs_model_object);
    objective_value = simplex_info.primal_objective_value;
    simplex_info.primal_objective_value = save_objective_value;
  } else {
    algorithm_name = "dual";
    have_previous_objective_value = have_previous_dual_objective_value;
    if (have_previous_objective_value) {
      previous_objective_value = previous_dual_objective_value;
      previous_updated_objective_value = previous_updated_dual_objective_value;
      updated_objective_correction = updated_dual_objective_correction;
    }
    updated_objective_value = simplex_info.updated_dual_objective_value;
    // Save the current objective value so that it can be recovered
    // after calling computeDualObjectiveValue
    double save_objective_value = simplex_info.dual_objective_value;
    computeDualObjectiveValue(highs_model_object, phase);
    objective_value = simplex_info.dual_objective_value;
    simplex_info.dual_objective_value = save_objective_value;
  }
  double change_in_objective_value = 0;
  double change_in_updated_objective_value = 0;
  if (have_previous_objective_value) {
    change_in_objective_value = objective_value - previous_objective_value;
    change_in_updated_objective_value =
        updated_objective_value - previous_updated_objective_value;
    updated_objective_value += updated_objective_correction;
  } else {
    updated_objective_correction = 0;
  }
  const double updated_objective_error =
      objective_value - updated_objective_value;
  const double updated_objective_absolute_error = fabs(updated_objective_error);
  const double updated_objective_relative_error =
      updated_objective_absolute_error / max(1.0, fabs(objective_value));
  updated_objective_correction += updated_objective_error;

  // Now update the records of previous objective value
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    have_previous_primal_objective_value = true;
    previous_primal_objective_value = objective_value;
    previous_updated_primal_objective_value = updated_objective_value;
    updated_primal_objective_correction = updated_objective_correction;
  } else {
    have_previous_dual_objective_value = true;
    previous_dual_objective_value = objective_value;
    previous_updated_dual_objective_value = updated_objective_value;
    updated_dual_objective_correction = updated_objective_correction;
  }

  // Now analyse the error
  HighsDebugStatus return_status = HighsDebugStatus::OK;
  if (updated_objective_relative_error >
          updated_objective_small_relative_error ||
      updated_objective_absolute_error >
          updated_objective_small_absolute_error) {
    if (updated_objective_relative_error >
            updated_objective_large_relative_error ||
        updated_objective_absolute_error >
            updated_objective_large_absolute_error) {
      HighsPrintMessage(
          highs_model_object.options_.output,
          highs_model_object.options_.message_level, ML_ALWAYS,
          "Updated %s objective value: large absolute (%g) or relative (%g) "
          "error"
          " - objective change - exact (%g) updated (%g) | %s\n",
          algorithm_name.c_str(), updated_objective_error,
          updated_objective_relative_error, change_in_objective_value,
          change_in_updated_objective_value, message.c_str());
      return_status = HighsDebugStatus::LARGE_ERROR;
    } else {
      HighsPrintMessage(
          highs_model_object.options_.output,
          highs_model_object.options_.message_level, ML_DETAILED,
          "Updated %s objective value: small absolute (%g) or relative (%g) "
          "error"
          " - objective change - exact (%g) updated (%g) | %s\n",
          algorithm_name.c_str(), updated_objective_error,
          updated_objective_relative_error, change_in_objective_value,
          change_in_updated_objective_value, message.c_str());
      return_status = HighsDebugStatus::SMALL_ERROR;
    }
    HighsPrintMessage(
        highs_model_object.options_.output,
        highs_model_object.options_.message_level, ML_VERBOSE,
        "Updated %s objective value: OK absolute (%g) or relative (%g) error\n",
        algorithm_name.c_str(), updated_objective_error,
        updated_objective_relative_error);
    return_status = HighsDebugStatus::OK;
  }
  return return_status;
}

HighsDebugStatus debugUpdatedObjectiveValue(
    const HighsModelObject& highs_model_object,
    const SimplexAlgorithm algorithm) {
  // Cheap check of updated objective value - assumes that the
  // objective value computed directly is correct, so only call after
  // this has been done
  if (highs_model_object.options_.highs_debug_level == HIGHS_DEBUG_LEVEL_NONE)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  std::string algorithm_name = "dual";
  if (algorithm == SimplexAlgorithm::PRIMAL) algorithm_name = "primal";
  double exact_objective;
  double updated_objective;
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    assert(highs_model_object.simplex_lp_status_.has_primal_objective_value);
    exact_objective = simplex_info.primal_objective_value;
    updated_objective = simplex_info.updated_primal_objective_value;
  } else {
    assert(highs_model_object.simplex_lp_status_.has_dual_objective_value);
    exact_objective = simplex_info.dual_objective_value;
    updated_objective = simplex_info.updated_dual_objective_value;
  }
  const double updated_objective_absolute_error =
      fabs(updated_objective - exact_objective);
  const double updated_objective_relative_error =
      updated_objective_absolute_error / max(1.0, fabs(exact_objective));

  HighsDebugStatus return_status = HighsDebugStatus::OK;
  // Now analyse the error
  if (updated_objective_relative_error >
          updated_objective_small_relative_error ||
      updated_objective_absolute_error >
          updated_objective_small_absolute_error) {
    if (updated_objective_relative_error >
            updated_objective_large_relative_error ||
        updated_objective_absolute_error >
            updated_objective_large_absolute_error) {
      HighsPrintMessage(
          highs_model_object.options_.output,
          highs_model_object.options_.message_level, ML_ALWAYS,
          "Updated %s objective value: large absolute (%g) or relative (%g) "
          "error\n",
          algorithm_name.c_str(), updated_objective_absolute_error,
          updated_objective_relative_error);
      return_status = HighsDebugStatus::LARGE_ERROR;
    } else {
      HighsPrintMessage(
          highs_model_object.options_.output,
          highs_model_object.options_.message_level, ML_DETAILED,
          "Updated %s objective value: small absolute (%g) or relative (%g) "
          "error\n",
          algorithm_name.c_str(), updated_objective_absolute_error,
          updated_objective_relative_error);
      return_status = HighsDebugStatus::SMALL_ERROR;
    }
    HighsPrintMessage(
        highs_model_object.options_.output,
        highs_model_object.options_.message_level, ML_VERBOSE,
        "Updated %s objective value: OK absolute (%g) or relative (%g) error\n",
        algorithm_name.c_str(), updated_objective_absolute_error,
        updated_objective_relative_error);
    return_status = HighsDebugStatus::OK;
  }
  return return_status;
}

HighsDebugStatus debugFixedNonbasicMove(
    const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check of nonbasicMove for fixed variables
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  int num_fixed_variable_move_errors = 0;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    if (simplex_info.workLower_[iVar] == simplex_info.workUpper_[iVar] &&
        simplex_basis.nonbasicMove_[iVar])
      num_fixed_variable_move_errors++;
  }
  assert(num_fixed_variable_move_errors == 0);
  if (num_fixed_variable_move_errors) {
    HighsPrintMessage(highs_model_object.options_.output,
                      highs_model_object.options_.message_level, ML_ALWAYS,
                      "There are %d fixed nonbasicMove errors",
                      num_fixed_variable_move_errors);
    return HighsDebugStatus::LOGICAL_ERROR;
  }
  return HighsDebugStatus::OK;
}

HighsDebugStatus debugNonbasicMove(const HighsModelObject& highs_model_object) {
  // Non-trivially expensive check of NonbasicMove
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  int num_free_variable_move_errors = 0;
  int num_lower_bounded_variable_move_errors = 0;
  int num_upper_bounded_variable_move_errors = 0;
  int num_boxed_variable_move_errors = 0;
  int num_fixed_variable_move_errors = 0;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = simplex_info.workLower_[iVar];
    const double upper = simplex_info.workUpper_[iVar];

    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free
        if (simplex_basis.nonbasicMove_[iVar]) {
          num_free_variable_move_errors++;
        }
      } else {
        // Only lower bounded
        if (simplex_basis.nonbasicMove_[iVar] != NONBASIC_MOVE_UP) {
          num_lower_bounded_variable_move_errors++;
        }
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded
        if (simplex_basis.nonbasicMove_[iVar] != NONBASIC_MOVE_DN) {
          num_upper_bounded_variable_move_errors++;
        }
      } else {
        // Boxed or fixed
        if (lower != upper) {
          // Boxed
          if (!simplex_basis.nonbasicMove_[iVar]) {
            num_boxed_variable_move_errors++;
          }
        } else {
          // Fixed
          if (simplex_basis.nonbasicMove_[iVar]) {
            num_fixed_variable_move_errors++;
          }
        }
      }
    }
  }
  int num_errors =
      num_free_variable_move_errors + num_lower_bounded_variable_move_errors +
      num_upper_bounded_variable_move_errors + num_boxed_variable_move_errors +
      num_fixed_variable_move_errors;

  if (num_errors) {
    HighsPrintMessage(
        highs_model_object.options_.output,
        highs_model_object.options_.message_level, ML_ALWAYS,
        "There are %d nonbasicMove errors: %d free; %d lower; %d upper; %d "
        "boxed; %d fixed",
        num_errors, num_free_variable_move_errors,
        num_lower_bounded_variable_move_errors,
        num_upper_bounded_variable_move_errors, num_boxed_variable_move_errors,
        num_fixed_variable_move_errors);
  }
  assert(num_errors == 0);
  if (num_errors) return HighsDebugStatus::LOGICAL_ERROR;
  return HighsDebugStatus::OK;
}

HighsDebugStatus debugBasisCondition(const HighsModelObject& highs_model_object,
                                     const std::string message) {
  // Non-trivially expensive assessment of basis condition
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  double basis_condition = computeBasisCondition(highs_model_object);
  if (basis_condition > excessive_basis_condition) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "%s basis condition estimate of %g is excessive",
                    message.c_str(), basis_condition);
    return HighsDebugStatus::WARNING;
  } else if (basis_condition > large_basis_condition) {
    HighsLogMessage(highs_model_object.options_.logfile,
                    HighsMessageType::WARNING,
                    "%s basis condition estimate of %g is large",
                    message.c_str(), basis_condition);
    return HighsDebugStatus::WARNING;
  } else if (basis_condition > fair_basis_condition) {
    HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                    "%s basis condition estimate of %g is fair",
                    message.c_str(), basis_condition);
    return HighsDebugStatus::OK;
  }
  HighsLogMessage(highs_model_object.options_.logfile, HighsMessageType::INFO,
                  "%s basis condition estimate of %g is small", message.c_str(),
                  basis_condition);
  return HighsDebugStatus::OK;
}

HighsDebugStatus debugCleanup(HighsModelObject& highs_model_object,
                              const std::vector<double>& original_dual) {
  if (highs_model_object.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsLp& simplex_lp = highs_model_object.simplex_lp_;
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const SimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  HighsSimplexAnalysis& analysis = highs_model_object.simplex_analysis_;
  // Make sure that the original_dual has been set up
  assert((int)original_dual.size() == simplex_lp.numCol_ + simplex_lp.numRow_);
  const std::vector<double>& new_dual = simplex_info.workDual_;

  const double dual_feasibility_tolerance =
      highs_model_object.scaled_solution_params_.dual_feasibility_tolerance;
  int num_dual_sign_change = 0;
  double norm_nonbasic_dual_change = 0;
  for (int iVar = 0; iVar < simplex_lp.numCol_ + simplex_lp.numRow_; iVar++) {
    if (!simplex_basis.nonbasicFlag_[iVar]) continue;
    const double dual_change = std::fabs(new_dual[iVar] - original_dual[iVar]);
    updateValueDistribution(dual_change,
                            analysis.cleanup_dual_change_distribution);
    norm_nonbasic_dual_change += dual_change;
    const double max_dual =
        std::max(std::fabs(new_dual[iVar]), std::fabs(original_dual[iVar]));
    if (max_dual > dual_feasibility_tolerance &&
        new_dual[iVar] * original_dual[iVar] < 0)
      num_dual_sign_change++;
  }
  printf(
      "Dual cleanup for %s has |dual change| = %g, with %d meaningful sign "
      "change(s)\n",
      simplex_lp.model_name_.c_str(), norm_nonbasic_dual_change,
      num_dual_sign_change);
  return HighsDebugStatus::OK;
}
