/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsDebug.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsDebug.h"

#include <algorithm>  // For std::max

HighsStatus debugDebugToHighsStatus(const HighsDebugStatus debug_status) {
  switch (debug_status) {
    case HighsDebugStatus::NOT_CHECKED:
    case HighsDebugStatus::OK:
    case HighsDebugStatus::SMALL_ERROR:
      return HighsStatus::OK;
    case HighsDebugStatus::WARNING:
      return HighsStatus::Warning;
    case HighsDebugStatus::LARGE_ERROR:
    case HighsDebugStatus::LOGICAL_ERROR:
      return HighsStatus::Error;
  }
  return HighsStatus::OK;
}

HighsDebugStatus debugWorseStatus(const HighsDebugStatus status0,
                                  const HighsDebugStatus status1) {
  return static_cast<HighsDebugStatus>(std::max((int)status0, (int)status1));
}

bool rightSizeVector(FILE* logfile, const std::string name0,
                     const std::string name1, const std::vector<double> v,
                     const int right_size) {
  const int v_size = v.size();
  if (v_size != right_size) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "%s: %s size is %d, not numCol = %d", name0.c_str(),
                    name1.c_str(), v_size, right_size);
    return false;
  }
  return true;
}

bool rightSizeVector(FILE* logfile, const std::string name0,
                     const std::string name1, const std::vector<int> v,
                     const int right_size) {
  const int v_size = v.size();
  if (v_size != right_size) {
    HighsLogMessage(logfile, HighsMessageType::ERROR,
                    "%s: %s size is %d, not numCol = %d", name0.c_str(),
                    name1.c_str(), v_size, right_size);
    return false;
  }
  return true;
}
