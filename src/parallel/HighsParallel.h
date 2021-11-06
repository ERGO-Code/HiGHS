/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_PARALLEL_H_
#define HIGHS_PARALLEL_H_

#include "parallel/HighsMutex.h"
#include "parallel/HighsTaskExecutor.h"

namespace highs {

namespace parallel {

using mutex = HighsMutex;

inline void initialize_scheduler(int numThreads = 0) {
  if (numThreads == 0)
    numThreads = (std::thread::hardware_concurrency() + 1) / 2;
  HighsTaskExecutor::initialize(numThreads);
}

inline int num_threads() { return HighsTaskExecutor::getNumWorkerThreads(); }

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
      HighsTaskExecutor::getGlobalTaskExecutor()->sync_stolen_task(
          localDeque, popResult.second);
      break;
    case HighsSplitDeque::Status::kWork:
      popResult.second->run();
  }
}

inline void sync() { sync(HighsTaskExecutor::getThisWorkerDeque()); }

template <typename F>
void for_each(HighsInt start, HighsInt end, F&& f, HighsInt grainSize = 1) {
  if (end - start <= grainSize) {
    f(start, end);
    return;
  }

  HighsSplitDeque* workerDeque = HighsTaskExecutor::getThisWorkerDeque();

  HighsInt split = (start + end) >> 1;
  spawn(workerDeque,
        [split, end, grainSize, &f]() { for_each(split, end, f, grainSize); });

  for_each(start, split, f, grainSize);

  sync(workerDeque);
}

}  // namespace parallel

}  // namespace highs

#endif