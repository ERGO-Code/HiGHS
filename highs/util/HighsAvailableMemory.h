/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsAvailableMemory.h
 * @brief Memory available to the process
 */
#ifndef UTIL_HIGHSAVAILABLEMEMORY_H_
#define UTIL_HIGHSAVAILABLEMEMORY_H_

#include <new>

// Memory that can be used without swapping, in bytes, or negative if unknown.
// On Linux this respects the cgroup limits of a container or systemd unit. A
// process limit (ulimit) needs no check, since allocation then fails
double availableMemoryBytes();

// Thrown when a solver stops because it's predicted not to fit in memory.
// This is planned, so it's logged as a warning rather than an error
struct HighsPredictedOutOfMemory : std::bad_alloc {
  const char* what() const noexcept override {
    return "predicted to run out of memory";
  }
};

#endif  // UTIL_HIGHSAVAILABLEMEMORY_H_
