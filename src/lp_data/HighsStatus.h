#ifndef LP_DATA_HIGHS_STATUS_H_
#define LP_DATA_HIGHS_STATUS_H_

#include <string>

// HiGHS status
enum class HighsStatus {
  OK = 0,
  Warning,
  Error
};

// Report a HighsStatus.
void HighsStatusReport(const char* message, HighsStatus status);

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status);

// Return the maximum of two HighsStatus
HighsStatus worseStatus(HighsStatus status0, HighsStatus status1);
#endif
