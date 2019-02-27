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
  Init,
  LpError,
  OptionsError,
  PresolveError,
  SolutionError,
  PostsolveError,
  NotImplemented,
  ReachedDualObjectiveUpperBound,
  Unbounded,
  Infeasible,
  Feasible,
  Optimal,
  Timeout,
  ReachedIterationLimit,
  NumericalDifficulties
};

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status);

// Return the maximum of two HighsStatus
HighsStatus worse_status(HighsStatus status0, HighsStatus status1);
#endif
