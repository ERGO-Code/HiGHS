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
#include "HighsOptions.h"
#include "HighsIO.h"

bool setUserOptionValue(HighsOptions& options, const std::string& option, const std::string& value) {
  if (option == "presolve") {
    if (value == "on")
      options.presolve_option = PresolveOption::ON;
    else if (value == "off")
      options.presolve_option = PresolveOption::OFF;
    else
      return false;
  } else if (option == "crash") {
    if (value == "on")
      options.crash_option = CrashOption::ON;
    else if (value == "off")
      options.crash_option = CrashOption::OFF;
    else
      return false;
  } else if (option == "parallel") {
    if (value == "on")
      options.pami = true;
    else if (value == "off")
      options.pami = false;
    else
      return false;
  } else if (option == "ipm") {
    if (value == "on")
      options.ipx = true;
    else if (value == "off")
      options.ipx = false;
    else
      return false;
  } else if (option == "simplex") {
    if (value == "on")
      options.simplex_option = SimplexOption::ON;
    else if (value == "off")
      options.simplex_option = SimplexOption::OFF;
    else
      return false;
  } else {
    return false;
  }

  return true;
}

// Used for options read from file or set by the user from another code.
bool setOptionValue(HighsOptions& options, const std::string& option, const std::string& value) {
 if (option == "presolve") {
    if (value == "on")
      options.presolve_option = PresolveOption::ON;
    else if (value == "off")
      options.presolve_option = PresolveOption::OFF;
    else
      return false;
  } else if (option == "crash") {
    if (value == "on")
      options.crash_option = CrashOption::ON;
    else if (value == "off")
      options.crash_option = CrashOption::OFF;
    else
      return false;
  } else if (option == "parallel") {
    if (value == "on")
      options.pami = true;
    else if (value == "off")
      options.pami = false;
    else
      return false;
  } else if (option == "ipm") {
    if (value == "on")
      options.ipx = true;
    else if (value == "off")
      options.ipx = false;
    else
      return false;
  } else if (option == "simplex") {
    if (value == "on")
      options.simplex_option = SimplexOption::ON;
    else if (value == "off")
      options.simplex_option = SimplexOption::OFF;
    else
      return false;
  } else if (option == "small_matrix_value") {
      options.small_matrix_value = atof(value.c_str());
  } else {
    HighsLogMessage(HighsMessageType::WARNING, "Unknown option: %s.", option.c_str());
    return false;
  }

  return true;
}


bool checkOptionsValue(HighsOptions& options) {
  if (options.simplex_option == SimplexOption::ON && options.ipx == true) {
    HighsPrintMessage(ML_MINIMAL, "Options error: both simplex and ipm set to on.\n");
    return false;
  }
  return true;
}