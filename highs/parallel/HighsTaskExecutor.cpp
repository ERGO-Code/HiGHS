/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "parallel/HighsTaskExecutor.h"

#include "parallel/HighsParallel.h"

#if defined(__linux__)
#include <sched.h>

#include <fstream>
#include <set>
#include <string>
#elif defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>

#include <vector>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#endif

using namespace highs;

// Fallback: assume 2 vthreads per core
static unsigned int fallback_core_count() {
  return (std::thread::hardware_concurrency() + 1) / 2;
}

// Returns the number of physical cores available to this process,
// respecting the OS affinity mask.
unsigned int highs::parallel::available_core_count() {
#if defined(__linux__)
  // Get the affinity mask, then count unique physical core IDs
  cpu_set_t affinity;
  if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0) {
    std::set<int> physical_cores;
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
      if (!CPU_ISSET(cpu, &affinity)) continue;
      std::ifstream topology_file("/sys/devices/system/cpu/cpu" +
                                  std::to_string(cpu) + "/topology/core_id");
      if (topology_file.is_open()) {
        int core_id;
        topology_file >> core_id;
        physical_cores.insert(core_id);
      }
    }
    if (!physical_cores.empty())
      return static_cast<unsigned int>(physical_cores.size());
  }
#elif defined(_WIN32) || defined(_WIN64)
  // Count physical cores whose logical processors overlap with process affinity
  DWORD_PTR process_mask, system_mask;
  if (!GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask))
    return fallback_core_count();

  // First call with nullptr to query the required buffer size
  DWORD length = 0;
  GetLogicalProcessorInformation(nullptr, &length);
  if (GetLastError() != ERROR_INSUFFICIENT_BUFFER) return fallback_core_count();

  std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION> info(
      length / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
  if (!GetLogicalProcessorInformation(info.data(), &length))
    return fallback_core_count();

  unsigned int physical_cores = 0;
  for (const auto& entry : info) {
    if (entry.Relationship != RelationProcessorCore) continue;
    if (entry.ProcessorMask & process_mask) physical_cores++;
  }
  if (physical_cores > 0) return physical_cores;
#elif defined(__APPLE__)
  // macOS does not support affinity; return physical core count
  int n;
  size_t size = sizeof(n);
  if (sysctlbyname("hw.physicalcpu", &n, &size, nullptr, 0) == 0 && n > 0)
    return static_cast<unsigned int>(n);
#endif
  return fallback_core_count();
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
