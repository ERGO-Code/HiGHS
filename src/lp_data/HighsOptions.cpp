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
  }
  // todo: else if (option == "ipm")
  // todo: else if (option == "simplex")

  return true;
}