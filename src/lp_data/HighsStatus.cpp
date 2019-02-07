#include "HighsStatus.h"

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status) {
  switch (status) {
    case HighsStatus::OK:
      return "OK";
      break;
    case HighsStatus::Init:
      return "Init";
      break;
    case HighsStatus::LpError:
      return "Lp Error";
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
    case HighsStatus::NotImplemented:
      return "Not implemented";
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
  }
  return "";
}