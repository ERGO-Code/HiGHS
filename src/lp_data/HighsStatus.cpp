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
      return "Status toString() not implemented.";
      break;
  }
  return "";
}

HighsStatus worseStatus(HighsStatus status0, HighsStatus status1) {
  HighsStatus return_status = HighsStatus::Error;
  if (status0 == HighsStatus::Error || status1 == HighsStatus::Error)
    return_status = HighsStatus::Error;
  else if (status0 == HighsStatus::Warning || status1 == HighsStatus::Warning)
    return_status = HighsStatus::Warning;
  else
    return_status = HighsStatus::OK;
  return return_status;
}
