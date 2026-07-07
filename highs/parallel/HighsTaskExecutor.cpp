/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "parallel/HighsTaskExecutor.h"

#include "parallel/HighsParallel.h"
#include "util/HighsHash.h"

#if defined(__linux__)
#include <sched.h>
#elif defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

using namespace highs;

// Like std::thread::hardware_concurrency(), but respects the process
// CPU affinity mask set by the OS (e.g., taskset or start /affinity).
unsigned int highs::parallel::available_concurrency() {
#if defined(__linux__)
  // Query the set of CPUs this process is allowed to run on
  cpu_set_t set;
  if (sched_getaffinity(0, sizeof(set), &set) == 0) {
    int count = CPU_COUNT(&set);
    if (count > 0) return static_cast<unsigned int>(count);
  }
#elif defined(_WIN32) || defined(_WIN64)
  // Query the process affinity bitmask and count set bits
  DWORD_PTR process_mask, system_mask;
  if (GetProcessAffinityMask(GetCurrentProcess(), &process_mask,
                             &system_mask)) {
    int count = HighsHashHelpers::popcnt(static_cast<uint64_t>(process_mask));
    if (count > 0) return static_cast<unsigned int>(count);
  }
#endif
  // Fallback when affinity query is unavailable or fails
  return std::thread::hardware_concurrency();
}

#ifdef _WIN32
static thread_local HighsSplitDeque* threadLocalWorkerDequePtr{nullptr};
HighsSplitDeque*& HighsTaskExecutor::threadLocalWorkerDeque() {
  return threadLocalWorkerDequePtr;
}

static thread_local HighsTaskExecutor::ExecutorHandle globalExecutorHandle{};

HighsTaskExecutor::ExecutorHandle&
HighsTaskExecutor::threadLocalExecutorHandle() {
  return globalExecutorHandle;
}
#else
thread_local HighsSplitDeque* HighsTaskExecutor::threadLocalWorkerDequePtr{
    nullptr};
thread_local HighsTaskExecutor::ExecutorHandle
    HighsTaskExecutor::globalExecutorHandle{};
#endif

void HighsTaskExecutor::ExecutorHandle::dispose() {
  if (ptr == nullptr) return;
  if (isMain) {
    ptr->stopWorkerThreads(false);
  }

  // check to see if we are the last handle and if so, delete the executor
  if (--ptr->referenceCount == 0) {
    cache_aligned::Deleter<HighsTaskExecutor>()(ptr);
  }

  ptr = nullptr;
}
