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
#ifndef HIGHS_TASKEXECUTOR_H_
#define HIGHS_TASKEXECUTOR_H_

#include <cassert>
#include <chrono>
#include <condition_variable>
#include <memory>
#include <thread>
#include <vector>

#include "parallel/HighsCacheAlign.h"
#include "parallel/HighsSplitDeque.h"
#include "util/HighsInt.h"
#include "util/HighsRandom.h"

class HighsTaskExecutor {
  static constexpr int kTaskArraySize = 8192;
  static constexpr int kTaskSize = 64;
  static constexpr int kNumTryFac = 16;
  static constexpr int kMicroSecsBeforeSleep = 5000;
  static constexpr int kMicroSecsBeforeGlobalSync = 1000;

  using cache_aligned = highs::cache_aligned;

 public:
  using WorkerDeque = HighsSplitDeque<kTaskArraySize, kTaskSize>;
  using Task = WorkerDeque::Task;

 private:
  static thread_local WorkerDeque* threadLocalWorkerDeque;
  static cache_aligned::shared_ptr<HighsTaskExecutor> globalExecutor;

  std::vector<cache_aligned::unique_ptr<WorkerDeque>> workerDeques;
  cache_aligned::shared_ptr<WorkerDeque::GlobalQueue> globalQueue;

  Task* random_steal_loop(WorkerDeque* localDeque) {
    const int numWorkers = workerDeques.size();

    int numTries = 16 * (numWorkers - 1);

    auto tStart = std::chrono::high_resolution_clock::now();

    while (true) {
      for (int s = 0; s < numTries; ++s) {
        Task* task = localDeque->randomSteal(workerDeques.data(), numWorkers);
        if (task) return task;
      }

      auto numMicroSecs =
          std::chrono::duration_cast<std::chrono::microseconds>(
              std::chrono::high_resolution_clock::now() - tStart)
              .count();

      if (numMicroSecs < kMicroSecsBeforeGlobalSync)
        numTries *= 2;
      else
        break;
    }

    return nullptr;
  }

  void run_worker(int workerId) {
    WorkerDeque* localDeque = workerDeques[workerId].get();
    threadLocalWorkerDeque = localDeque;
    Task* currentTask =
        globalQueue->waitForNewTask(localDeque, kMicroSecsBeforeSleep);
    while (true) {
      assert(currentTask != nullptr);

      currentTask->run(localDeque);

      currentTask = random_steal_loop(localDeque);
      if (currentTask != nullptr) continue;

      currentTask =
          globalQueue->waitForNewTask(localDeque, kMicroSecsBeforeSleep);
    }
  }

 public:
  HighsTaskExecutor(int numThreads) {
    assert(numThreads > 0);
    workerDeques.resize(numThreads);
    globalQueue = cache_aligned::make_shared<WorkerDeque::GlobalQueue>(
        workerDeques.data());

    for (int i = 0; i < numThreads; ++i)
      workerDeques[i] = cache_aligned::make_unique<WorkerDeque>(globalQueue, i);

    threadLocalWorkerDeque = workerDeques[0].get();
    for (int i = 1; i < numThreads; ++i)
      std::thread([&](int id) { run_worker(id); }, i).detach();
  }

  static WorkerDeque* getThisWorkerDeque() { return threadLocalWorkerDeque; }

  static HighsTaskExecutor* getGlobalTaskExecutor() {
    return globalExecutor.get();
  }

  static int getNumWorkerThreads() {
    return globalExecutor->workerDeques.size();
  }

  static void initialize(int numThreads) {
    if (!globalExecutor)
      globalExecutor =
          cache_aligned::make_shared<HighsTaskExecutor>(numThreads);
  }

  void sync_stolen_task(WorkerDeque* localDeque, Task* stolenTask) {
    if (!stolenTask->leapfrog(localDeque)) {
      const int numWorkers = workerDeques.size();
      int numTries = kNumTryFac * (numWorkers - 1);

      auto tStart = std::chrono::high_resolution_clock::now();

      while (true) {
        for (int s = 0; s < numTries; ++s) {
          if (stolenTask->isFinished()) {
            localDeque->popStolen();
            return;
          }
          Task* task = localDeque->randomSteal(workerDeques.data(), numWorkers);
          if (task) task->run(localDeque);
        }

        auto numMicroSecs =
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - tStart)
                .count();

        if (numMicroSecs < kMicroSecsBeforeSleep)
          numTries *= 2;
        else
          break;
      }

      if (!stolenTask->isFinished())
        stolenTask->waitForTaskToFinish(localDeque);
    }

    localDeque->popStolen();
  }
};

#endif
