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

OptionStatus setOptionValue(HighsOptions& options, const std::string& option,
                            const std::string& value) {
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

  else if (option == mps_parser_type_string)
    return setParserTypeValue(options, value);

  else if (option == mip_string)
    return setMipValue(options, value);

  else if (option == find_feasibility_string)
    return setFindFeasibilityValue(options, value);

  else if (option == find_feasibility_strategy_string)
    return setFindFeasibilityStrategyValue(options, value);

  else if (option == find_feasibility_dualize_string)
    return setFindFeasibilityDualizeValue(options, value);

  else if (option == run_as_hsol_string)
    return setRunAsHsolValue(options, atoi(value.c_str()));

  else if (option == infinite_cost_string)
    return setInfiniteCostValue(options, atof(value.c_str()));

  else if (option == infinite_bound_string)
    return setInfiniteBoundValue(options, atof(value.c_str()));

  else if (option == small_matrix_value_string)
    return setSmallMatrixValueValue(options, atof(value.c_str()));

  else if (option == large_matrix_value_string)
    return setLargeMatrixValueValue(options, atof(value.c_str()));

  else if (option == allowed_simplex_scale_factor_string)
    return setAllowedSimplexScaleFactorValue(options, atoi(value.c_str()));

  else if (option == primal_feasibility_tolerance_string)
    return setPrimalFeasibilityToleranceValue(options, atof(value.c_str()));

  else if (option == dual_feasibility_tolerance_string)
    return setDualFeasibilityToleranceValue(options, atof(value.c_str()));

  else if (option == dual_objective_value_upper_bound_string)
    return setDualObjectiveValueUpperBoundValue(options, atof(value.c_str()));

  else if (option == simplex_strategy_string)
    return setSimplexStrategyValue(options, atoi(value.c_str()));

  else if (option == simplex_dualise_strategy_string)
    return setSimplexDualiseStrategyValue(options, atoi(value.c_str()));

  else if (option == simplex_permute_strategy_string)
    return setSimplexPermuteStrategyValue(options, atoi(value.c_str()));

  else if (option == simplex_scale_strategy_string)
    return setSimplexScaleStrategyValue(options, atoi(value.c_str()));

  else if (option == simplex_crash_strategy_string)
    return setSimplexCrashStrategyValue(options, atoi(value.c_str()));

  else if (option == simplex_dual_edge_weight_strategy_string)
    return setSimplexDualEdgeWeightStrategyValue(options, atoi(value.c_str()));

  else if (option == simplex_primal_edge_weight_strategy_string)
    return setSimplexPrimalEdgeWeightStrategyValue(options,
                                                   atoi(value.c_str()));

  else if (option == simplex_price_strategy_string)
    return setSimplexPriceStrategyValue(options, atoi(value.c_str()));

  else if (option == simplex_initial_condition_check_string)
    return setSimplexInitialConditionCheckValue(options, atoi(value.c_str()));

  else if (option == simplex_initial_condition_tolerance_string)
    return setSimplexInitialConditionToleranceValue(options,
                                                    atof(value.c_str()));

  else if (option == message_level_string)
    return setMessageLevelValue(options, atoi(value.c_str()));

  else {
    HighsLogMessage(HighsMessageType::WARNING, "Unknown option: %s.",
                    option.c_str());
    return OptionStatus::UNKNOWN_OPTION;
  }

  return OptionStatus::OK;
}

OptionStatus checkOptionsValue(HighsOptions& options) {
  if (options.simplex_option == SimplexOption::ON && options.ipx == true) {
    HighsPrintMessage(ML_MINIMAL,
                      "Options error: both simplex and ipm set to on.\n");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

// Set values of options so that HiGHS runs as Hsol
void setHsolOptions(HighsOptions& options) {
  // Set command line options to their hsol values
  options.presolve_option = PresolveOption::OFF;
  options.crash_option = CrashOption::OFF;
  options.parallel_option = ParallelOption::OFF;
  options.simplex_option = SimplexOption::ON;
  options.highs_run_time_limit = HIGHS_CONST_INF;
  options.simplex_iteration_limit = HIGHS_CONST_I_INF;
  options.mps_parser_type = HighsMpsParserType::fixed;
  options.infinite_cost = HIGHS_CONST_INF;
  options.infinite_bound = HIGHS_CONST_INF;
  options.small_matrix_value = 0;
  options.large_matrix_value = HIGHS_CONST_INF;
  options.allowed_simplex_scale_factor = HIGHS_CONST_I_INF;
  options.primal_feasibility_tolerance = 1e-7;
  options.dual_feasibility_tolerance = 1e-7;
  options.dual_objective_value_upper_bound = HIGHS_CONST_INF;
  options.simplex_strategy = SimplexStrategy::DUAL_PLAIN;
  options.simplex_dualise_strategy = SimplexDualiseStrategy::OFF;
  options.simplex_permute_strategy = SimplexPermuteStrategy::OFF;
  options.simplex_scale_strategy = SimplexScaleStrategy::HSOL;
  options.simplex_crash_strategy = SimplexCrashStrategy::OFF;
  options.simplex_dual_edge_weight_strategy = SimplexDualEdgeWeightStrategy::STEEPEST_EDGE;
  options.simplex_primal_edge_weight_strategy = SimplexPrimalEdgeWeightStrategy::DANTZIG;
  options.simplex_price_strategy = SimplexPriceStrategy::ROW;
}

void reportStringOptionValue(const int report_level, const string option_string,
                             const string option_value,
                             const string option_default) {
  string default_space;
  if (report_level) {
    default_space = "       ";
  } else {
    default_space = "";
  }
  bool is_default = option_value == option_default;
  if (!is_default || report_level) {
    if (is_default) {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has default value \"%s\"\n",
		      option_string.c_str(), option_default.c_str());
    } else {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s value \"%s\": default value is \"%s\"",
		      option_string.c_str(), default_space.c_str(), option_value.c_str(), option_default.c_str());
    }
  }
}

void reportIntOptionValue(const int report_level, const string option_string,
                          const int option_value, const int option_default,
                          const int* option_min, const int* option_max) {
  string default_space;
  if (report_level) {
    default_space = "       ";
  } else {
    default_space = "";
  }
  char value_char [100];
  char range_char [100];
  int value_num_char;
  int range_num_char;
  bool is_default = option_value == option_default;
  if (!is_default || report_level > 0) {
    if (is_default) {
      value_num_char = sprintf(value_char, " default value %12d", option_default);
    } else {
      value_num_char = sprintf(value_char, "%s value %12d: default value is %12d", default_space.c_str(), option_value, option_default);
    }
    if (option_min == NULL) {
      if (option_max == NULL) {
        range_num_char = sprintf(range_char, " ");
      } else {
        range_num_char = sprintf(range_char, ": valid range is [-Inf, %6d]", *option_max);
      }
    } else {
      if (option_max == NULL) {
        range_num_char = sprintf(range_char, ": valid range is [%6d, Inf]", *option_min);
      } else {
        range_num_char = sprintf(range_char, ": valid range is [%6d, %6d]", *option_min, *option_max);
      }
    }
    HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s%s", option_string.c_str(), value_char, range_char);
  }
}

void reportDoubleOptionValue(const int report_level, const string option_string,
                             const double option_value,
                             const double option_default,
                             const double* option_min,
                             const double* option_max) {
  string default_space;
  if (report_level) {
    default_space = "       ";
  } else {
    default_space = "";
  }
  char value_char [100];
  char range_char [100];
  int value_num_char;
  int range_num_char;
  bool is_default = option_value == option_default;
  if (!is_default || report_level > 0) {
    if (is_default) {
      value_num_char = sprintf(value_char, " default value %12g", option_default);
    } else {
      value_num_char = sprintf(value_char,
	      "%s value %12g: default value is %12g", default_space.c_str(), option_value, option_default);
    }
    if (option_min == NULL) {
      if (option_max == NULL) {
        range_num_char = sprintf(range_char, " ");
      } else {
        range_num_char = sprintf(range_char, ": valid range is [-Inf, %12g]", *option_max);
      }
    } else {
      if (option_max == NULL) {
        range_num_char = sprintf(range_char, ": valid range is [%12g, Inf]", *option_min);
      } else {
        range_num_char = sprintf(range_char, ": valid range is [%12g, %12g]", *option_min, *option_max);
      }
    }
    HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s%s", option_string.c_str(), value_char, range_char);
  }
}

void reportOptionsValue(const HighsOptions& options, const int report_level) {
  // Report on the command line options
  bool is_default;
  string default_space;
  if (report_level) {
    default_space = "       ";
  } else {
    default_space = "";
  }
  // Model file name
  //
  // Don't report this since hiGHS returns an error if it's not
  // changed from the default.
  //
  //  reportStringOptionValue(report_level, file_string, options.filename,
  //                          FILENAME_DEFAULT);
  // Options file name
  reportStringOptionValue(report_level, options_file_string,
                          options.options_file, OPTIONS_FILE_DEFAULT);
  // Presolve option
  is_default = options.presolve_option == PresolveOption::DEFAULT;
  if (!is_default || report_level) {
    if (is_default) {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has default value \"off\"",
             presolve_string.c_str());
    } else {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s value \"on\": default value is \"off\"",
	     presolve_string.c_str(), default_space.c_str());
    }
  }
  // Crash option
  is_default = options.crash_option == CrashOption::DEFAULT;
  if (!is_default || report_level) {
    if (is_default) {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has default value \"off\"", crash_string.c_str());
    } else {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s value \"on\": default value is \"off\"",
	     crash_string.c_str(), default_space.c_str());
    }
  }
  // Parallel option
  is_default = options.parallel_option == ParallelOption::DEFAULT;
  if (!is_default || report_level) {
    if (is_default) {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has default value \"off\"",
             parallel_string.c_str());
    } else {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s value \"on\": default value is \"off\"",
	     parallel_string.c_str(), default_space.c_str());
    }
  }
  // Simplex option
  is_default = options.simplex_option == SimplexOption::DEFAULT;
  if (!is_default || report_level) {
    if (is_default) {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has default value \"off\"",
             simplex_string.c_str());
    } else {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s value \"on\": default value is \"off\"",
	     simplex_string.c_str(), default_space.c_str());
    }
  }
  // Ipx option ToDo This is a mess name-wise and should use
  // IpmOption::OFF/ON/DEFAULT
  is_default = options.ipx == false;
  if (!is_default || report_level) {
    if (is_default) {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has default value \"false\"", ipm_string.c_str());
    } else {
      HighsLogMessage(HighsMessageType::INFO, "Option: %-32s has%s value \"true\": default value is "
	     "\"false\"",
	     ipm_string.c_str(), default_space.c_str());
    }
  }
  // HiGHS run time limit
  reportDoubleOptionValue(report_level, highs_run_time_limit_string,
                          options.highs_run_time_limit,
                          HIGHS_RUN_TIME_LIMIT_DEFAULT, NULL, NULL);
  // Simplex iteration limit
  reportIntOptionValue(report_level, simplex_iteration_limit_string,
                       options.simplex_iteration_limit,
                       SIMPLEX_ITERATION_LIMIT_DEFAULT, NULL, NULL);
  /*

const string mps_parser_type_string = "mps_parser_type";
const string mip_string = "mip";
const string find_feasibility_string = "find_feasibility";
const string find_feasibility_strategy_string = "feasibility_strategy";
const string find_feasibility_dualize_string = "feasibility_dualize";
  */

  // run_as_hsol
  reportIntOptionValue(report_level, run_as_hsol_string,
		       options.run_as_hsol, RUN_AS_HSOL_DEFAULT,
		       &RUN_AS_HSOL_MIN, &RUN_AS_HSOL_MAX);

  // infinite_cost
  reportDoubleOptionValue(report_level, infinite_cost_string,
                          options.infinite_cost, INFINITE_COST_DEFAULT,
                          &INFINITE_COST_MIN, &INFINITE_COST_MAX);

  // infinite_bound
  reportDoubleOptionValue(report_level, infinite_bound_string,
                          options.infinite_bound, INFINITE_BOUND_DEFAULT,
                          &INFINITE_BOUND_MIN, &INFINITE_BOUND_MAX);

  // small_matrix_value
  reportDoubleOptionValue(report_level, small_matrix_value_string,
                          options.small_matrix_value,
                          SMALL_MATRIX_VALUE_DEFAULT, &SMALL_MATRIX_VALUE_MIN,
                          &SMALL_MATRIX_VALUE_MAX);

  // large_matrix_value
  reportDoubleOptionValue(report_level, large_matrix_value_string,
                          options.large_matrix_value,
                          LARGE_MATRIX_VALUE_DEFAULT, &LARGE_MATRIX_VALUE_MIN,
                          &LARGE_MATRIX_VALUE_MAX);

  // Allowed simplex scale factor
  reportIntOptionValue(report_level, allowed_simplex_scale_factor_string,
                       options.allowed_simplex_scale_factor,
                       ALLOWED_SIMPLEX_SCALE_FACTOR_DEFAULT,
                       &ALLOWED_SIMPLEX_SCALE_FACTOR_MIN,
                       &ALLOWED_SIMPLEX_SCALE_FACTOR_MAX);

  // primal_feasibility_tolerance
  reportDoubleOptionValue(report_level, primal_feasibility_tolerance_string,
                          options.primal_feasibility_tolerance,
                          PRIMAL_FEASIBILITY_TOLERANCE_DEFAULT,
                          &PRIMAL_FEASIBILITY_TOLERANCE_MIN,
                          &PRIMAL_FEASIBILITY_TOLERANCE_MAX);

  // dual_feasibility_tolerance
  reportDoubleOptionValue(
      report_level, dual_feasibility_tolerance_string,
      options.dual_feasibility_tolerance, DUAL_FEASIBILITY_TOLERANCE_DEFAULT,
      &DUAL_FEASIBILITY_TOLERANCE_MIN, &DUAL_FEASIBILITY_TOLERANCE_MAX);

  // dual_objective_value_upper_bound
  reportDoubleOptionValue(report_level, dual_objective_value_upper_bound_string,
                          options.dual_objective_value_upper_bound,
                          DUAL_OBJECTIVE_VALUE_UPPER_BOUND_DEFAULT, NULL, NULL);

  // Simplex strategy
  reportIntOptionValue(report_level, simplex_strategy_string,
                       (int)options.simplex_strategy,
                       (int)SimplexStrategy::DEFAULT, NULL, NULL);

  // Simplex dualise strategy
  reportIntOptionValue(report_level, simplex_dualise_strategy_string,
                       (int)options.simplex_dualise_strategy,
                       (int)SimplexDualiseStrategy::DEFAULT, NULL, NULL);

  // Simplex permute strategy
  reportIntOptionValue(report_level, simplex_permute_strategy_string,
                       (int)options.simplex_permute_strategy,
                       (int)SimplexPermuteStrategy::DEFAULT, NULL, NULL);

  // Simplex scale strategy
  reportIntOptionValue(report_level, simplex_scale_strategy_string,
                       (int)options.simplex_scale_strategy,
                       (int)SimplexScaleStrategy::DEFAULT, NULL, NULL);

  // Simplex crash strategy
  reportIntOptionValue(report_level, simplex_crash_strategy_string,
                       (int)options.simplex_crash_strategy,
                       (int)SimplexCrashStrategy::DEFAULT, NULL, NULL);
}

OptionStatus setPresolveValue(HighsOptions& options, const std::string& value) {
  if (value == on_string)
    options.presolve_option = PresolveOption::ON;
  else if (value == off_string)
    options.presolve_option = PresolveOption::OFF;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "presolve option value \"%s\" is not permitted: legal "
                    "values are \"%s\" and \"%s\"\n",
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
                    "crash option value \"%s\" is not permitted: legal values "
                    "are \"%s\" and \"%s\"\n",
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
                    "parallel option value \"%s\" is not permitted: legal "
                    "values are \"%s\" and \"%s\"\n",
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
                    "simplex option value \"%s\" is not permitted: legal "
                    "values are \"%s\" and \"%s\"\n",
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
                    "ipm option value \"%s\" is not permitted: legal values "
                    "are \"%s\" and \"%s\"\n",
                    value.c_str(), on_string.c_str(), off_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setHighsRunTimeLimitValue(HighsOptions& options,
                                       const double& value) {
  if (value >= 0)
    options.highs_run_time_limit = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "HiGHS run time limit value \"%s\" is not permitted: legal "
                    "values are no less than %d\n",
                    value, 0);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setSimplexIterationLimitValue(HighsOptions& options,
                                           const int& value) {
  if (value >= 0)
    options.simplex_iteration_limit = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "Simplex iteration limit value \"%s\" is not permitted: "
                    "legal values are no less than %d\n",
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
                    "mip value \"%s\" is not permitted: legal values are "
                    "\"%s\" and \"%s\"\n",
                    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setFindFeasibilityValue(HighsOptions& options,
                                     const std::string& value) {
  if (value == "on" || value == "true")
    options.find_feasibility = true;
  else if (value == "off" || value == "false")
    options.find_feasibility = false;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "find_feasibility value \"%s\" is not permitted: legal "
                    "values are \"%s\" and \"%s\"\n",
                    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setFindFeasibilityStrategyValue(HighsOptions& options,
                                             const std::string& value) {
  if (value == "approx_component")
    options.feasibility_strategy = FeasibilityStrategy::kApproxComponentWise;
  else if (value == "approx_exact")
    options.feasibility_strategy = FeasibilityStrategy::kApproxExact;
  else if (value == "direct")
    options.feasibility_strategy = FeasibilityStrategy::kDirectSolve;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "feasibility component-wise value \"%s\" is not permitted: "
                    "legal values are \"%s\" and \"%s\"\n",
                    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setFindFeasibilityDualizeValue(HighsOptions& options,
                                            const std::string& value) {
  if (value == "on" || value == "true")
    options.feasibility_strategy_dualize = true;
  else if (value == "off" || value == "false")
    options.feasibility_strategy_dualize = false;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "feasibility dualize value \"%s\" is not permitted: legal "
                    "values are \"%s\" and \"%s\"\n",
                    value.c_str(), "on", "off");
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setParserTypeValue(HighsOptions& options,
                                const std::string& value) {
  if (value == "fixed")
    options.mps_parser_type = HighsMpsParserType::fixed;
  else if (value == "free")
    options.mps_parser_type = HighsMpsParserType::free;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "parser type value \"%s\" is not permitted: legal values "
                    "are \"%s\" and \"%s\"\n",
                    value.c_str(), fixed_string.c_str(), free_string.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setRunAsHsolValue(HighsOptions& options, const int& value) {
  if (value >= RUN_AS_HSOL_MIN && value <= RUN_AS_HSOL_MAX)
    options.run_as_hsol = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "infinite cost value \"%s\" is not permitted: legal values "
                    "are between %d and %d\n",
                    value, RUN_AS_HSOL_MIN, RUN_AS_HSOL_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setInfiniteCostValue(HighsOptions& options, const double& value) {
  if (value >= INFINITE_COST_MIN && value <= INFINITE_COST_MAX)
    options.infinite_cost = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "infinite cost value \"%s\" is not permitted: legal values "
                    "are between %d and %d\n",
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
                    "infinite bound value \"%s\" is not permitted: legal "
                    "values are between %d and %d\n",
                    value, INFINITE_BOUND_MIN, INFINITE_BOUND_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setSmallMatrixValueValue(HighsOptions& options,
                                      const double& value) {
  if (value >= SMALL_MATRIX_VALUE_MIN && value <= SMALL_MATRIX_VALUE_MAX)
    options.small_matrix_value = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "small matrix value \"%s\" is not permitted: legal values "
                    "are between %d and %d\n",
                    value, SMALL_MATRIX_VALUE_MIN, SMALL_MATRIX_VALUE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setLargeMatrixValueValue(HighsOptions& options,
                                      const double& value) {
  if (value >= LARGE_MATRIX_VALUE_MIN && value <= LARGE_MATRIX_VALUE_MAX)
    options.large_matrix_value = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "large matrix value \"%s\" is not permitted: legal values "
                    "are between %d and %d\n",
                    value, LARGE_MATRIX_VALUE_MIN, LARGE_MATRIX_VALUE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setAllowedSimplexScaleFactorValue(HighsOptions& options,
                                               const int& value) {
  if (value >= ALLOWED_SIMPLEX_SCALE_FACTOR_MIN &&
      value <= ALLOWED_SIMPLEX_SCALE_FACTOR_MAX)
    options.allowed_simplex_scale_factor = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "allowed simplex scale factor value \"%s\" is not "
                    "permitted: legal values are between %d and %d\n",
                    value, ALLOWED_SIMPLEX_SCALE_FACTOR_MIN,
                    ALLOWED_SIMPLEX_SCALE_FACTOR_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setPrimalFeasibilityToleranceValue(HighsOptions& options,
                                                const double& value) {
  if (value >= PRIMAL_FEASIBILITY_TOLERANCE_MIN &&
      value <= PRIMAL_FEASIBILITY_TOLERANCE_MAX)
    options.primal_feasibility_tolerance = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "primal feasibility tolerance value \"%s\" is not "
                    "permitted: legal values are between %d and %d\n",
                    value, PRIMAL_FEASIBILITY_TOLERANCE_MIN,
                    PRIMAL_FEASIBILITY_TOLERANCE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setDualFeasibilityToleranceValue(HighsOptions& options,
                                              const double& value) {
  if (value >= DUAL_FEASIBILITY_TOLERANCE_MIN &&
      value <= DUAL_FEASIBILITY_TOLERANCE_MAX)
    options.dual_feasibility_tolerance = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "dual feasibility tolerance value \"%s\" is not permitted: "
                    "legal values are between %d and %d\n",
                    value, DUAL_FEASIBILITY_TOLERANCE_MIN,
                    DUAL_FEASIBILITY_TOLERANCE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setDualObjectiveValueUpperBoundValue(HighsOptions& options,
                                                  const double& value) {
  options.dual_objective_value_upper_bound = value;
  return OptionStatus::OK;
}

OptionStatus setSimplexStrategyValue(HighsOptions& options, const int& value) {
  options.simplex_strategy = intToSimplexStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexDualiseStrategyValue(HighsOptions& options,
                                            const int& value) {
  options.simplex_dualise_strategy = intToSimplexDualiseStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexPermuteStrategyValue(HighsOptions& options,
                                            const int& value) {
  options.simplex_permute_strategy = intToSimplexPermuteStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexScaleStrategyValue(HighsOptions& options,
                                          const int& value) {
  options.simplex_scale_strategy = intToSimplexScaleStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexCrashStrategyValue(HighsOptions& options,
                                          const int& value) {
  options.simplex_crash_strategy = intToSimplexCrashStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexDualEdgeWeightStrategyValue(HighsOptions& options,
                                                   const int& value) {
  options.simplex_dual_edge_weight_strategy =
      intToSimplexDualEdgeWeightStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexPrimalEdgeWeightStrategyValue(HighsOptions& options,
                                                     const int& value) {
  options.simplex_primal_edge_weight_strategy =
      intToSimplexPrimalEdgeWeightStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexPriceStrategyValue(HighsOptions& options,
                                          const int& value) {
  options.simplex_price_strategy = intToSimplexPriceStrategy(value);
  return OptionStatus::OK;
}

OptionStatus setSimplexInitialConditionCheckValue(HighsOptions& options,
                                                  const int& value) {
  options.simplex_initial_condition_check = value;
  return OptionStatus::OK;
}

OptionStatus setSimplexInitialConditionToleranceValue(HighsOptions& options,
                                                      const double& value) {
  if (value >= SIMPLEX_INITIAL_CONDITION_TOLERANCE_MIN &&
      value <= SIMPLEX_INITIAL_CONDITION_TOLERANCE_MAX)
    options.simplex_initial_condition_tolerance = value;
  else {
    HighsLogMessage(HighsMessageType::ERROR,
                    "simplex initial condition tolerance value \"%s\" is not "
                    "permitted: legal values are between %d and %d\n",
                    value, SIMPLEX_INITIAL_CONDITION_TOLERANCE_MIN,
                    SIMPLEX_INITIAL_CONDITION_TOLERANCE_MAX);
    return OptionStatus::ILLEGAL_VALUE;
  }
  return OptionStatus::OK;
}

OptionStatus setMessageLevelValue(HighsOptions& options, const int& value) {
  HighsSetMessagelevel(value);
  options.messageLevel = value;
  return OptionStatus::OK;
}

SimplexStrategy intToSimplexStrategy(const int& value) {
  if (value == (int)SimplexStrategy::CHOOSE) return SimplexStrategy::CHOOSE;
  if (value == (int)SimplexStrategy::DUAL_PLAIN)
    return SimplexStrategy::DUAL_PLAIN;
  if (value == (int)SimplexStrategy::DUAL_MULTI)
    return SimplexStrategy::DUAL_MULTI;
  if (value == (int)SimplexStrategy::PRIMAL) return SimplexStrategy::PRIMAL;
  return SimplexStrategy::DEFAULT;
}
SimplexDualiseStrategy intToSimplexDualiseStrategy(const int& value) {
  if (value == (int)SimplexDualiseStrategy::OFF)
    return SimplexDualiseStrategy::OFF;
  if (value == (int)SimplexDualiseStrategy::CHOOSE)
    return SimplexDualiseStrategy::CHOOSE;
  if (value == (int)SimplexDualiseStrategy::ON)
    return SimplexDualiseStrategy::ON;
  return SimplexDualiseStrategy::DEFAULT;
}
SimplexPermuteStrategy intToSimplexPermuteStrategy(const int& value) {
  if (value == (int)SimplexPermuteStrategy::OFF)
    return SimplexPermuteStrategy::OFF;
  if (value == (int)SimplexPermuteStrategy::CHOOSE)
    return SimplexPermuteStrategy::CHOOSE;
  if (value == (int)SimplexPermuteStrategy::ON)
    return SimplexPermuteStrategy::ON;
  return SimplexPermuteStrategy::DEFAULT;
}
SimplexScaleStrategy intToSimplexScaleStrategy(const int& value) {
  if (value == (int)SimplexScaleStrategy::OFF) return SimplexScaleStrategy::OFF;
  if (value == (int)SimplexScaleStrategy::CHOOSE)
    return SimplexScaleStrategy::CHOOSE;
  if (value == (int)SimplexScaleStrategy::HSOL)
    return SimplexScaleStrategy::HSOL;
  if (value == (int)SimplexScaleStrategy::HIGHS)
    return SimplexScaleStrategy::HIGHS;
  return SimplexScaleStrategy::DEFAULT;
}
SimplexCrashStrategy intToSimplexCrashStrategy(const int& value) {
  if (value == (int)SimplexCrashStrategy::LTSSF_K)
    return SimplexCrashStrategy::LTSSF_K;
  if (value == (int)SimplexCrashStrategy::LTSSF_PRI)
    return SimplexCrashStrategy::LTSSF_PRI;
  if (value == (int)SimplexCrashStrategy::BIXBY)
    return SimplexCrashStrategy::BIXBY;
  printf("Crash strategy option %d cannot be parsed yet!\n", value);
  return SimplexCrashStrategy::DEFAULT;
}
SimplexDualEdgeWeightStrategy intToSimplexDualEdgeWeightStrategy(
    const int& value) {
  if (value == (int)SimplexDualEdgeWeightStrategy::DANTZIG)
    return SimplexDualEdgeWeightStrategy::DANTZIG;
  if (value == (int)SimplexDualEdgeWeightStrategy::DEVEX)
    return SimplexDualEdgeWeightStrategy::DEVEX;
  if (value == (int)SimplexDualEdgeWeightStrategy::STEEPEST_EDGE)
    return SimplexDualEdgeWeightStrategy::STEEPEST_EDGE;
  if (value ==
      (int)SimplexDualEdgeWeightStrategy::STEEPEST_EDGE_TO_DEVEX_SWITCH)
    return SimplexDualEdgeWeightStrategy::STEEPEST_EDGE_TO_DEVEX_SWITCH;
  return SimplexDualEdgeWeightStrategy::DEFAULT;
}
SimplexPrimalEdgeWeightStrategy intToSimplexPrimalEdgeWeightStrategy(
    const int& value) {
  if (value == (int)SimplexPrimalEdgeWeightStrategy::DANTZIG)
    return SimplexPrimalEdgeWeightStrategy::DANTZIG;
  if (value == (int)SimplexPrimalEdgeWeightStrategy::DEVEX)
    return SimplexPrimalEdgeWeightStrategy::DEVEX;
  return SimplexPrimalEdgeWeightStrategy::DEFAULT;
}
SimplexPriceStrategy intToSimplexPriceStrategy(const int& value) {
  if (value == (int)SimplexPriceStrategy::COL) return SimplexPriceStrategy::COL;
  if (value == (int)SimplexPriceStrategy::ROW) return SimplexPriceStrategy::ROW;
  if (value == (int)SimplexPriceStrategy::ROW_SWITCH)
    return SimplexPriceStrategy::ROW_SWITCH;
  if (value == (int)SimplexPriceStrategy::ROW_SWITCH_COL_SWITCH)
    return SimplexPriceStrategy::ROW_SWITCH_COL_SWITCH;
  return SimplexPriceStrategy::DEFAULT;
}
