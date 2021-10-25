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

#include "parallel/HighsTaskExecutor.h"

namespace highs {

namespace parallel {

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
void spawn(HighsTaskExecutor::WorkerDeque* localDeque, F&& f) {
  localDeque->push(std::forward<F>(f));
}

template <typename F>
void spawn(F&& f) {
  spawn(HighsTaskExecutor::getThisWorkerDeque(), std::forward<F>(f));
}

inline void sync(HighsTaskExecutor::WorkerDeque* localDeque) {
  std::pair<HighsTaskExecutor::WorkerDeque::Status, HighsTaskExecutor::Task*>
      popResult = localDeque->pop();
  switch (popResult.first) {
    case HighsTaskExecutor::WorkerDeque::Status::kEmpty:
      assert(false);
      // fall through
    case HighsTaskExecutor::WorkerDeque::Status::kOverflown:
      // when the local deque is overflown the task has been executed during
      // spawn already
      break;
    case HighsTaskExecutor::WorkerDeque::Status::kStolen:
      HighsTaskExecutor::getGlobalTaskExecutor()->sync_stolen_task(
          localDeque, popResult.second);
      break;
    case HighsTaskExecutor::WorkerDeque::Status::kWork:
      popResult.second->run();
  }
}

inline void sync() { sync(HighsTaskExecutor::getThisWorkerDeque()); }

template <typename F>
void for_each(HighsInt start, HighsInt end, F&& f, HighsInt grainSize = 1) {
  HighsTaskExecutor::WorkerDeque* workerDeque =
      HighsTaskExecutor::getThisWorkerDeque();

  HighsInt numTasks = 0;
  HighsInt numThreads = 2 * HighsTaskExecutor::getNumWorkerThreads();
  while (end - start > grainSize) {
    HighsInt split = (start + end) >> 1;
    ++numTasks;
    workerDeque->pushPrivate(
        [split, end, &f, grainSize]() { for_each(split, end, f, grainSize); });
    end = split;

    if ((numThreads >> numTasks) == 1) workerDeque->publishTasks(numTasks);
  }

  if ((numThreads >> numTasks) > 1) workerDeque->publishTasks(numTasks);

  f(start, end);

  for (HighsInt i = 0; i < numTasks; ++i) sync(workerDeque);
}

template <typename F>
void for_each_static(HighsInt start, HighsInt end, F&& f) {
  HighsTaskExecutor::WorkerDeque* workerDeque =
      HighsTaskExecutor::getThisWorkerDeque();

  HighsInt totalSize = end - start;
  HighsInt numTasks =
      std::min(HighsTaskExecutor::getNumWorkerThreads(), totalSize);
  double chunkMult = double(totalSize) / numTasks;

  for (HighsInt i = numTasks - 1; i > 0; --i) {
    HighsInt chunkStart = i * chunkMult;
    workerDeque->pushPrivate([chunkStart, end, &f] { f(chunkStart, end); });
    end = chunkStart;
  }

  workerDeque->publishTasks(numTasks);

  f(start, end);

  for (HighsInt j = 0; j < numTasks; ++j) sync(workerDeque);
}

}  // namespace parallel

}  // namespace highs

#endif