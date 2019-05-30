#ifndef LP_DATA_HIGHS_STATUS_H_
#define LP_DATA_HIGHS_STATUS_H_

#include <string>

// HiGHS status
enum class HighsStatus
{
  NotSet,
  OK,
  Info,
  Warning,
  Error,
  NotImplemented,
  Init,
  LpError,
  OptionsError,
  PresolveError,
  SolutionError,
  PostsolveError,
  LpEmpty,
  ReachedDualObjectiveUpperBound,
  Unbounded,
  Infeasible,
  Feasible,
  Optimal,
  Timeout,
  ReachedIterationLimit,
  NumericalDifficulties
};

// Report a HighsStatus.
void HighsStatusReport(const char* message, HighsStatus status);

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status);

// Return the maximum of two HighsStatus
HighsStatus worseStatus(HighsStatus status0, HighsStatus status1);
#endif
