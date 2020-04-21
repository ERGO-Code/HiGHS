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

HighsDebugStatus debugComputedDual(const HighsModelObject& workHMO,
                                   const std::vector<double>& previous_dual,
                                   const std::vector<double>& basic_costs,
                                   const std::vector<double>& row_dual) {
  // Non-trivially expensive analysis of computed dual values.
  if (workHMO.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  const std::vector<double>& new_daul = workHMO.simplex_info_.workDual_;

  int num_row = workHMO.simplex_lp_.numRow_;
  int num_col = workHMO.simplex_lp_.numCol_;
  double basic_costs_norm = 0;
  double row_dual_norm = 0;
  for (int iRow = 0; iRow < num_row; iRow++) {
    basic_costs_norm += fabs(basic_costs[iRow]);
    row_dual_norm += fabs(row_dual[iRow]);
  }
  double new_dual_norm = 0;
  for (int iVar = 0; iVar < num_row + num_col; iVar++)
    new_dual_norm += fabs(new_daul[iVar]);
  double computed_dual_absolute_change = 0;
  double computed_dual_relative_change = 0;
  if (previous_dual.size()) {
    for (int iVar = 0; iVar < num_row + num_col; iVar++)
      computed_dual_absolute_change +=
          fabs(new_daul[iVar] - previous_dual[iVar]);
    if (new_dual_norm)
      computed_dual_relative_change =
          computed_dual_absolute_change / new_dual_norm;
  }
  if (new_dual_norm) {
    if (computed_dual_relative_change > computed_dual_small_relative_change ||
        computed_dual_absolute_change > computed_dual_small_absolute_change) {
      if (computed_dual_relative_change > computed_dual_large_relative_change ||
          computed_dual_absolute_change > computed_dual_large_absolute_change) {
        HighsPrintMessage(
            workHMO.options_.output, workHMO.options_.message_level,
            ML_DETAILED,
            "ComputedDual: B.pi=c_B has |c_B|=%g; |pi|=%g; |pi^TA-c|=%g\n",
            basic_costs_norm, row_dual_norm, new_dual_norm);
        HighsPrintMessage(
            workHMO.options_.output, workHMO.options_.message_level,
            ML_DETAILED,
            "ComputedDual: Large absolute (%g) or relative (%g) change\n",
            computed_dual_absolute_change, computed_dual_relative_change);
      } else {
        HighsPrintMessage(
            workHMO.options_.output, workHMO.options_.message_level,
            ML_DETAILED,
            "ComputedDual: B.pi=c_B has |c_B|=%g; |pi|=%g; |pi^TA-c|=%g\n",
            basic_costs_norm, row_dual_norm, new_dual_norm);
        HighsPrintMessage(
            workHMO.options_.output, workHMO.options_.message_level,
            ML_DETAILED,
            "ComputedDual: Small absolute (%g) or relative (%g) change\n",
            computed_dual_absolute_change, computed_dual_relative_change);
      }
      HighsPrintMessage(
          workHMO.options_.output, workHMO.options_.message_level, ML_VERBOSE,
          "ComputedDual: B.pi=c_B has |c_B|=%g; |pi|=%g; |pi^TA-c|=%g\n",
          basic_costs_norm, row_dual_norm, new_dual_norm);
      HighsPrintMessage(
          workHMO.options_.output, workHMO.options_.message_level, ML_VERBOSE,
          "ComputedDual: OK absolute (%g) or relative (%g) change\n",
          computed_dual_absolute_change, computed_dual_relative_change);
    }
  } else {
    HighsLogMessage(workHMO.options_.logfile, HighsMessageType::WARNING,
                    "ComputedDual: |Dual| = %g", new_dual_norm);
    return HighsDebugStatus::WARNING;
  }
  return HighsDebugStatus::OK;
}

HighsDebugStatus debugUpdatedObjectiveValue(HighsModelObject& workHMO,
                                            const SimplexAlgorithm algorithm,
                                            const int phase,
                                            const std::string message) {
  // Non-trivially expensive check of updated objective value. Computes the
  // exact objective value
  if (workHMO.options_.highs_debug_level < HIGHS_DEBUG_LEVEL_COSTLY)
    return HighsDebugStatus::NOT_CHECKED;
  HighsSimplexInfo& simplex_info = workHMO.simplex_info_;

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
    computePrimalObjectiveValue(workHMO);
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
    computeDualObjectiveValue(workHMO, phase);
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
  const double updated_objective_error = objective_value - updated_objective_value;
  const double updated_objective_absolute_error =
      fabs(updated_objective_error);
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
          workHMO.options_.output, workHMO.options_.message_level, ML_ALWAYS,
          "Updated %s objective value: large absolute (%g) or relative (%g) "
          "error"
          " - objective change - exact (%g) updated (%g) | %s\n",
          algorithm_name.c_str(), updated_objective_error,
          updated_objective_relative_error, change_in_objective_value,
          change_in_updated_objective_value, message.c_str());
      return_status = HighsDebugStatus::LARGE_ERROR;
    } else {
      HighsPrintMessage(
          workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
          "Updated %s objective value: small absolute (%g) or relative (%g) "
          "error"
          " - objective change - exact (%g) updated (%g) | %s\n",
          algorithm_name.c_str(), updated_objective_error,
          updated_objective_relative_error, change_in_objective_value,
          change_in_updated_objective_value, message.c_str());
      return_status = HighsDebugStatus::SMALL_ERROR;
    }
    HighsPrintMessage(
        workHMO.options_.output, workHMO.options_.message_level, ML_VERBOSE,
        "Updated %s objective value: OK absolute (%g) or relative (%g) error\n",
        algorithm_name.c_str(), updated_objective_error,
        updated_objective_relative_error);
    return_status = HighsDebugStatus::OK;
  }
  return return_status;
}

HighsDebugStatus debugUpdatedObjectiveValue(const HighsModelObject& workHMO,
                                            const SimplexAlgorithm algorithm) {
  // Cheap check of updated objective value - assumes that the
  // objective value computed directly is correct, so only call after
  // this has been done
  if (workHMO.options_.highs_debug_level == HIGHS_DEBUG_LEVEL_NONE)
    return HighsDebugStatus::NOT_CHECKED;
  const HighsSimplexInfo& simplex_info = workHMO.simplex_info_;
  std::string algorithm_name = "dual";
  if (algorithm == SimplexAlgorithm::PRIMAL) algorithm_name = "primal";
  double exact_objective;
  double updated_objective;
  if (algorithm == SimplexAlgorithm::PRIMAL) {
    assert(workHMO.simplex_lp_status_.has_primal_objective_value);
    exact_objective = simplex_info.primal_objective_value;
    updated_objective = simplex_info.updated_primal_objective_value;
  } else {
    assert(workHMO.simplex_lp_status_.has_dual_objective_value);
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
          workHMO.options_.output, workHMO.options_.message_level, ML_ALWAYS,
          "Updated %s objective value: large absolute (%g) or relative (%g) "
          "error\n",
          algorithm_name.c_str(), updated_objective_absolute_error,
          updated_objective_relative_error);
      return_status = HighsDebugStatus::LARGE_ERROR;
    } else {
      HighsPrintMessage(
          workHMO.options_.output, workHMO.options_.message_level, ML_DETAILED,
          "Updated %s objective value: small absolute (%g) or relative (%g) "
          "error\n",
          algorithm_name.c_str(), updated_objective_absolute_error,
          updated_objective_relative_error);
      return_status = HighsDebugStatus::SMALL_ERROR;
    }
    HighsPrintMessage(
        workHMO.options_.output, workHMO.options_.message_level, ML_VERBOSE,
        "Updated %s objective value: OK absolute (%g) or relative (%g) error\n",
        algorithm_name.c_str(), updated_objective_absolute_error,
        updated_objective_relative_error);
    return_status = HighsDebugStatus::OK;
  }
  return return_status;
}
