/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_PARALLEL_H_
#define HIGHS_PARALLEL_H_

#include <thread>

#if defined(__linux__)
#include <sched.h>
#elif defined(_WIN32)
#include <windows.h>
#endif

#include "parallel/HighsMutex.h"
#include "parallel/HighsTaskExecutor.h"

namespace highs {

namespace parallel {

using mutex = HighsMutex;

inline unsigned int available_concurrency() {
#if defined(__linux__)
  cpu_set_t set;
  if (sched_getaffinity(0, sizeof(set), &set) == 0) {
    int count = CPU_COUNT(&set);
    if (count > 0) return static_cast<unsigned int>(count);
  }
#elif defined(_WIN32)
  DWORD_PTR process_mask, system_mask;
  if (GetProcessAffinityMask(GetCurrentProcess(), &process_mask,
                             &system_mask)) {
    int count = static_cast<int>(__popcnt64(process_mask));
    if (count > 0) return static_cast<unsigned int>(count);
  }
#endif
  return std::thread::hardware_concurrency();
}

inline void initialize_scheduler(int numThreads = 0) {
  if (numThreads == 0) {
#ifdef HIGHS_NO_DEFAULT_THREADS
    numThreads = 1;
#else
    numThreads = (available_concurrency() + 1) / 2;
#endif
  }
  HighsTaskExecutor::initialize(numThreads);
}

inline int num_threads() {
  return HighsTaskExecutor::getThisWorkerDeque()->getNumWorkers();
}

inline int thread_num() {
  return HighsTaskExecutor::getThisWorkerDeque()->getOwnerId();
}

template <typename F>
void spawn(HighsSplitDeque* localDeque, F&& f) {
  localDeque->push(std::forward<F>(f));
}

template <typename F>
void spawn(F&& f) {
  spawn(HighsTaskExecutor::getThisWorkerDeque(), std::forward<F>(f));
}

inline void sync(HighsSplitDeque* localDeque) {
  std::pair<HighsSplitDeque::Status, HighsTask*> popResult = localDeque->pop();
  switch (popResult.first) {
    case HighsSplitDeque::Status::kEmpty:
      assert(false);
      // fall through
    case HighsSplitDeque::Status::kOverflown:
      // when the local deque is overflown the task has been executed during
      // spawn already
      break;
    case HighsSplitDeque::Status::kStolen:
      HighsTaskExecutor::sync_stolen_task(localDeque, popResult.second);
      break;
    case HighsSplitDeque::Status::kWork:
      popResult.second->run();
  }
}

inline void sync() { sync(HighsTaskExecutor::getThisWorkerDeque()); }
class TaskGroup {
  HighsSplitDeque* workerDeque;
  int dequeHead;

 public:
  TaskGroup() {
    workerDeque = HighsTaskExecutor::getThisWorkerDeque();
    dequeHead = workerDeque->getCurrentHead();
  }

  template <typename F>
  void spawn(F&& f) const {
    highs::parallel::spawn(workerDeque, std::forward<F>(f));
  }

  void sync() const {
    assert(workerDeque->getCurrentHead() > dequeHead);
    highs::parallel::sync(workerDeque);
  }

  void taskWait() const {
    while (workerDeque->getCurrentHead() > dequeHead)
      highs::parallel::sync(workerDeque);
  }

  void cancel() {
    for (int i = dequeHead; i < workerDeque->getCurrentHead(); ++i)
      workerDeque->cancelTask(i);
  }

  ~TaskGroup() noexcept(false) {
    cancel();

    // Setting the rootTask to null ensures that taskWait does not throw
    // exceptions, so that all tasks can be successfully completed before
    // destroying the object.
    HighsTask* prevRoot = workerDeque->setRootTask(nullptr);

    taskWait();

    // Reinstate the original rootTask and potentially throw exception.
    // dtor needs to be noexcept(false) for this to work.
    workerDeque->setRootTask(prevRoot);
    workerDeque->checkInterrupt();
  }
};

template <typename F>
void for_each(HighsInt start, HighsInt end, F&& f, HighsInt grainSize = 1) {
  if (end - start <= grainSize) {
    f(start, end);
  } else {
    TaskGroup tg;

    do {
      HighsInt split = (start + end) >> 1;
      tg.spawn([split, end, grainSize, &f]() {
        for_each(split, end, f, grainSize);
      });
      end = split;
    } while (end - start > grainSize);

    f(start, end);
    tg.taskWait();
  }
}

}  // namespace parallel

}  // namespace highs

#endif
