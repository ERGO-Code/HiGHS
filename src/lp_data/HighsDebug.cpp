/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
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
#include <cassert>    // For std::max

HighsStatus debugDebugToHighsStatus(const HighsDebugStatus debug_status) {
  switch (debug_status) {
    case HighsDebugStatus::kNotChecked:
    case HighsDebugStatus::kOk:
    case HighsDebugStatus::kSmallError:
      return HighsStatus::OK;
    case HighsDebugStatus::kWarning:
    case HighsDebugStatus::kLargeError:
      return HighsStatus::Warning;
    case HighsDebugStatus::kError:
    case HighsDebugStatus::kExcessiveError:
    case HighsDebugStatus::kLogicalError:
      return HighsStatus::Error;
    default:
      return HighsStatus::OK;
  }
}

HighsDebugStatus debugWorseStatus(const HighsDebugStatus status0,
                                  const HighsDebugStatus status1) {
  return static_cast<HighsDebugStatus>(
      std::max((HighsInt)status0, (HighsInt)status1));
}

bool debugVectorRightSize(const std::vector<double> v,
                          const HighsInt right_size) {
  const HighsInt v_size = v.size();
  const bool is_right_size = v_size == right_size;
  assert(is_right_size);
  return is_right_size;
}

bool debugVectorRightSize(const std::vector<HighsInt> v,
                          const HighsInt right_size) {
  const HighsInt v_size = v.size();
  const bool is_right_size = v_size == right_size;
  assert(is_right_size);
  return is_right_size;
}
