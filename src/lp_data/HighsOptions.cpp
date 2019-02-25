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

OptionStatus setUserOptionValue(HighsOptions& options, const std::string& option, const std::string& value) {
  if (option == "presolve")
    return setPresolveValue(options, value);
  else if (option == "crash")
    return setCrashValue(options, value);
  else if (option == "parallel")
    return setParallelValue(options, value);
  else if (option == "ipm")
    return setIpmValue(options, value);
  else if (option == "simplex")
    return setSimplexValue(options, value);
  else
    return OptionStatus::UNKNOWN_OPTION;

  return OptionStatus::OK;
}

// Used for options read from file or set by the user from another code.
OptionStatus setOptionValue(HighsOptions& options, const std::string& option, const std::string& value) {
  // Return if one of the main run-time options is set or is given an illegal value.
  OptionStatus return_status = setUserOptionValue(options, option, value);
  if (return_status == OptionStatus::OK || return_status == OptionStatus::ILLEGAL_VALUE) return return_status;

  // Value of option was unknown, so see if it's one of the other options
  assert(return_status == OptionStatus::UNKNOWN_OPTION);
											   
  if (option == "infinite_cost") {
      options.infinite_cost = atof(value.c_str());
  } else if (option == "infinite_bound") {
      options.infinite_bound = atof(value.c_str());
  } else if (option == "small_matrix_value") {
      options.small_matrix_value = atof(value.c_str());
  } else if (option == "large_matrix_value") {
      options.large_matrix_value = atof(value.c_str());
  } else if (option == "dual_objective_value_upper_bound") {
      options.dual_objective_value_upper_bound = atof(value.c_str());
      /*
 } else if (option == "parser_type") {
    if (value == "free")
      options.parser_type = HighsMpsParserType::free;
    else if (value == "fixed")
      options.parser_type = HighsMpsParserType::fixed;
    else
      return OptionStatus::ILLEGAL_VALUE;
      */
  } else {
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
