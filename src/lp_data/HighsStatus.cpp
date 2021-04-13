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
#include "lp_data/HighsStatus.h"

#include "io/HighsIO.h"

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status) {
  switch (status) {
    case HighsStatus::kOk:
      return "OK";
      break;
    case HighsStatus::kWarning:
      return "Warning";
      break;
    case HighsStatus::kError:
      return "Error";
      break;
    default:
#ifdef HiGHSDEV
      printf("HiGHS status %" HIGHSINT_FORMAT " not recognised\n",
             (HighsInt)status);
#endif
      return "Unrecognised HiGHS status";
      break;
  }
  return "";
}

HighsStatus interpretCallStatus(const HighsStatus call_status,
                                const HighsStatus from_return_status,
                                const std::string& message) {
  HighsStatus to_return_status;
  to_return_status = worseStatus(call_status, from_return_status);
#ifdef HiGHSDEV
  if (call_status != HighsStatus::kOk) {
    if (message != "") {
      printf("HighsStatus::%s return from %s\n",
             HighsStatusToString(call_status).c_str(), message.c_str());
    } else {
      printf("HighsStatus::%s return\n",
             HighsStatusToString(call_status).c_str());
    }
  }
#endif
  return to_return_status;
}

HighsStatus worseStatus(const HighsStatus status0, const HighsStatus status1) {
  HighsStatus return_status = HighsStatus::kError;
  if (status0 == HighsStatus::kError || status1 == HighsStatus::kError)
    return_status = HighsStatus::kError;
  else if (status0 == HighsStatus::kWarning || status1 == HighsStatus::kWarning)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;
  return return_status;
}
