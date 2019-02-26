#include "lp_data/HighsStatus.h"
#include "io/HighsIO.h"

// Report a HighsStatus.
void HighsStatusReport(const char* message, HighsStatus status) {
  HighsLogMessage(HighsMessageType::INFO, "%s: HighsStatus = %d - %s\n", message, (int)status, HighsStatusToString(status).c_str());
}

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status) {
  switch (status) {
    case HighsStatus::NotSet:
      return "Not Set";
      break;
    case HighsStatus::OK:
      return "OK";
      break;
    case HighsStatus::Info:
      return "Info";
      break;
    case HighsStatus::Warning:
      return "Warning";
      break;
    case HighsStatus::Error:
      return "Error";
      break;
    case HighsStatus::NotImplemented:
      return "Not implemented";
      break;
    case HighsStatus::Init:
      return "Init";
      break;
    case HighsStatus::LpError:
      return "LP Error";
      break;
    case HighsStatus::OptionsError:
      return "Options Error";
      break;
    case HighsStatus::PresolveError:
      return "Presolve Error";
      break;
    case HighsStatus::SolutionError:
      return "Solution Error";
      break;
    case HighsStatus::PostsolveError:
      return "Postsolve Error";
      break;
    case HighsStatus::LpEmpty:
      return "LP is empty";
      break;
    case HighsStatus::Unbounded:
      return "Unbounded";
      break;
    case HighsStatus::Infeasible:
      return "Infeasible";
      break;
    case HighsStatus::Feasible:
      return "Feasible";
      break;
    case HighsStatus::Optimal:
      return "Optimal";
      break;
    case HighsStatus::Timeout:
      return "Timeout";
      break;
    default:
      return "Status toString() not implemented.";
      break;
  }
  return "";
}

HighsStatus worseStatus(HighsStatus status0, HighsStatus status1) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (status0 == HighsStatus::Error || status1 == HighsStatus::Error)
    return_status = HighsStatus::Error;
  else if (status0 == HighsStatus::Warning || status1 == HighsStatus::Warning)
    return_status = HighsStatus::Warning;
  else if (status0 == HighsStatus::Info || status1 == HighsStatus::Info)
    return_status = HighsStatus::Info;
  else 
    return_status = HighsStatus::OK;
  return return_status;
}
