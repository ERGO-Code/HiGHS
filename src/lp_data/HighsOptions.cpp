/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsOptions.cpp
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsOptions.h"
#include "io/HighsIO.h"

OptionStatus setOptionValue(HighsOptions& options, const std::string& option, const std::string& value) {
  if (option == presolve_string)
    return setPresolveValue(options, value);

  else if (option == crash_string)
    return setCrashValue(options, value);

  else if (option == parallel_string)
    return setParallelValue(options, value);

  else if (option == simplex_string)
    return setSimplexValue(options, value);

  else if (option == ipm_string)
    return setIpmValue(options, value);

  else if (option == highs_run_time_limit_string) 
    return setHighsRunTimeLimitValue(options, atof(value.c_str()));

  else if (option == simplex_iteration_limit_string) 
    return setSimplexIterationLimitValue(options, atoi(value.c_str()));

  else if (option == parser_type_string) 
    return setParserTypeValue(options, value);

  else if (option == mip_string)
    return setMipValue(options, value);

  else if (option == find_feasibility_string)
    return setFindFeasibilityValue(options, value);

  else if (option == find_feasibility_strategy_string)
    return setFindFeasibilityStrategyValue(options, value);

  else if (option == find_feasibility_dualize_string)
    return setFindFeasibilityDualizeValue(options, value);

  else if (option == infinite_cost_string) 
    return setInfiniteCostValue(options, atof(value.c_str()));

   else if (option == infinite_bound_string) 
    return setInfiniteBoundValue(options, atof(value.c_str()));

   else if (option == small_matrix_value_string) 
    return setSmallMatrixValueValue(options, atof(value.c_str()));

   else if (option == large_matrix_value_string) 
    return setLargeMatrixValueValue(options, atof(value.c_str()));

   else if (option == primal_feasibility_tolerance_string) 
    return setPrimalFeasibilityToleranceValue(options, atof(value.c_str()));

   else if (option == dual_feasibility_tolerance_string) 
    return setDualFeasibilityToleranceValue(options, atof(value.c_str()));

   else if (option == dual_objective_value_upper_bound_string) 
    return setDualObjectiveValueUpperBoundValue(options, atof(value.c_str()));

   else {
    HighsLogMessage(HighsMessageType::WARNING, "Unknown option: %s.", option.c_str());
    return OptionStatus::UNKNOWN_OPTION;
  }

  return OptionStatus::OK;
}


OptionStatus checkOptionsValue(HighsOptions& options) {
  if (options.simplex_option == SimplexOption::ON && options.ipx == true) {
    HighsPrintMessage(ML_MINIMAL, "Options error: both simplex and ipm set to on.\n");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setPresolveValue(HighsOptions& options, const std::string& value) {
  if (value == on_string)
    options.presolve_option = PresolveOption::ON;
  else if (value == off_string)
    options.presolve_option = PresolveOption::OFF;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "presolve option value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), on_string.c_str(), off_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setCrashValue(HighsOptions& options, const std::string& value) {
  if (value == on_string)
    options.crash_option = CrashOption::ON;
  else if (value == off_string)
    options.crash_option = CrashOption::OFF;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "crash option value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), on_string.c_str(), off_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setParallelValue(HighsOptions& options, const std::string& value) {
  if (value == on_string)
    options.pami = true;
  else if (value == off_string)
    options.pami = false;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "parallel option value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), on_string.c_str(), off_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setSimplexValue(HighsOptions& options, const std::string& value) {
  if (value == on_string)
    options.simplex_option = SimplexOption::ON;
  else if (value == off_string)
    options.simplex_option = SimplexOption::OFF;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "simplex option value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), on_string.c_str(), off_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setIpmValue(HighsOptions& options, const std::string& value) {
  if (value == on_string)
    options.ipx = true;
  else if (value == off_string)
    options.ipx = false;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "ipm option value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), on_string.c_str(), off_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setHighsRunTimeLimitValue(HighsOptions& options, const double& value) {
  if (value >= 0)
    options.highs_run_time_limit = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "HiGHS run time limit value \"%s\" is not permitted: legal values are no less than %d\n",
		    value, 0);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setSimplexIterationLimitValue(HighsOptions& options, const int& value) {
  if (value >= 0)
    options.simplex_iteration_limit = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "Simplex iteration limit value \"%s\" is not permitted: legal values are no less than %d\n",
		    value, 0);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setMipValue(HighsOptions& options, const std::string& value) {
  if (value == "on" || value == "true")
    options.mip = true;
  else if (value == "off" || value == "false")
    options.mip = false;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "mip value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setFindFeasibilityValue(HighsOptions& options, const std::string& value) {
  if (value == "on" || value == "true")
    options.find_feasibility = true;
  else if (value == "off" || value == "false")
    options.find_feasibility = false;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "find_feasibility value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setFindFeasibilityStrategyValue(HighsOptions& options, const std::string& value) {
  if (value == "approx_component")
    options.feasibility_strategy = FeasibilityStrategy::kApproxComponentWise;
  else if (value == "approx_exact")
    options.feasibility_strategy = FeasibilityStrategy::kApproxExact;
  else if (value == "direct") 
    options.feasibility_strategy = FeasibilityStrategy::kDirectSolve;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "feasibility component-wise value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setFindFeasibilityDualizeValue(HighsOptions& options, const std::string& value) {
  if (value == "on" || value == "true")
    options.feasibility_strategy_dualize = true;
  else if (value == "off" || value == "false")
    options.feasibility_strategy_dualize = false;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "feasibility dualize value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setParserTypeValue(HighsOptions& options, const std::string& value) {
  if (value == "fixed")
    options.parser_type = HighsMpsParserType::fixed;
  else if (value == "free")
    options.parser_type = HighsMpsParserType::free;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "parser type value \"%s\" is not permitted: legal values are \"%s\" and \"%s\"\n",
		    value.c_str(), fixed_string.c_str(), free_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setInfiniteCostValue(HighsOptions& options, const double& value) {
  if (value >= INFINITE_COST_MIN && value <= INFINITE_COST_MAX)
    options.infinite_cost = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "infinite cost value \"%s\" is not permitted: legal values are between %d and %d\n",
		    value, INFINITE_COST_MIN, INFINITE_COST_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setInfiniteBoundValue(HighsOptions& options, const double& value) {
  if (value >= INFINITE_BOUND_MIN && value <= INFINITE_BOUND_MAX)
    options.infinite_bound = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "infinite bound value \"%s\" is not permitted: legal values are between %d and %d\n",
		    value, INFINITE_BOUND_MIN, INFINITE_BOUND_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setSmallMatrixValueValue(HighsOptions& options, const double& value) {
  if (value >= SMALL_MATRIX_VALUE_MIN && value <= SMALL_MATRIX_VALUE_MAX)
    options.small_matrix_value = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "small matrix value \"%s\" is not permitted: legal values are between %d and %d\n",
		    value, SMALL_MATRIX_VALUE_MIN, SMALL_MATRIX_VALUE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setLargeMatrixValueValue(HighsOptions& options, const double& value) {
  if (value >= LARGE_MATRIX_VALUE_MIN && value <= LARGE_MATRIX_VALUE_MAX)
    options.large_matrix_value = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "large matrix value \"%s\" is not permitted: legal values are between %d and %d\n",
		    value, LARGE_MATRIX_VALUE_MIN, LARGE_MATRIX_VALUE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setPrimalFeasibilityToleranceValue(HighsOptions& options, const double& value) {
  if (value >= PRIMAL_FEASIBILITY_TOLERANCE_MIN && value <= PRIMAL_FEASIBILITY_TOLERANCE_MAX)
    options.primal_feasibility_tolerance = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "primal feasibility tolerance value \"%s\" is not permitted: legal values are between %d and %d\n",
		    value, PRIMAL_FEASIBILITY_TOLERANCE_MIN, PRIMAL_FEASIBILITY_TOLERANCE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setDualFeasibilityToleranceValue(HighsOptions& options, const double& value) {
  if (value >= DUAL_FEASIBILITY_TOLERANCE_MIN && value <= DUAL_FEASIBILITY_TOLERANCE_MAX)
    options.dual_feasibility_tolerance = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
		    "dual feasibility tolerance value \"%s\" is not permitted: legal values are between %d and %d\n",
		    value, DUAL_FEASIBILITY_TOLERANCE_MIN, DUAL_FEASIBILITY_TOLERANCE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setDualObjectiveValueUpperBoundValue(HighsOptions& options, const double& value) {
  options.dual_objective_value_upper_bound = value;
  return OptionStatus::OK;
}

