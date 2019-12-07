#include "lp_data/HighsStatus.h"
#include "io/HighsIO.h"

// Report a HighsStatus.
void HighsStatusReport(const char* message, HighsStatus status) {
  HighsLogMessage(HighsMessageType::INFO, "%s: HighsStatus = %d - %s\n",
                  message, (int)status, HighsStatusToString(status).c_str());
}

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status) {
  switch (status) {
    case HighsStatus::OK:
      return "OK";
      break;
    case HighsStatus::Warning:
      return "Warning";
      break;
    case HighsStatus::Error:
      return "Error";
      break;
    default:
#ifdef HiGHSDEV
      printf("HiGHS status %d not recognised\n", (int)status);
#endif
      return "Unrecognised HiGHS status";
      break;
  }
  return "";
}

HighsStatus interpretCallStatus(const HighsStatus call_status,
				const HighsStatus from_return_status,
				const std::string message) {
  HighsStatus to_return_status;
  to_return_status = worseStatus(call_status, from_return_status);
#ifdef HiGHSDEV
  if (call_status != HighsStatus::OK) {
    if (message != "") {
      printf("HighsStatus::%s return from %s\n",
	     HighsStatusToString(call_status).c_str(),
	     message.c_str());
    } else {
      printf("HighsStatus::%s return\n",
	     HighsStatusToString(call_status).c_str());
    }
  }
#endif
  return to_return_status;
}

HighsStatus worseStatus(const HighsStatus status0, const HighsStatus status1) {
  HighsStatus return_status = HighsStatus::Error;
  if (status0 == HighsStatus::Error || status1 == HighsStatus::Error)
    return_status = HighsStatus::Error;
  else if (status0 == HighsStatus::Warning || status1 == HighsStatus::Warning)
    return_status = HighsStatus::Warning;
  else
    return_status = HighsStatus::OK;
  return return_status;
}
