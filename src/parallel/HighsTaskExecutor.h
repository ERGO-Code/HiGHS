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

#include "libhighs_export.h"
#include "parallel/HighsCacheAlign.h"
#include "parallel/HighsSplitDeque.h"
#include "util/HighsInt.h"
#include "util/HighsRandom.h"
class HighsTaskExecutor {
 public:
  static constexpr int kNumTryFac = 16;
  static constexpr int kMicroSecsBeforeSleep = 5000;
  static constexpr int kMicroSecsBeforeGlobalSync = 1000;

 private:
  using cache_aligned = highs::cache_aligned;

#ifdef _WIN32
  static HighsSplitDeque*& threadLocalWorkerDeque();
#else
  static thread_local HighsSplitDeque* threadLocalWorkerDequePtr;
  static HighsSplitDeque*& threadLocalWorkerDeque() {
    return threadLocalWorkerDequePtr;
  }
#endif

  static LIBHIGHS_EXPORT cache_aligned::shared_ptr<HighsTaskExecutor>
      globalExecutor;
  static LIBHIGHS_EXPORT HighsSpinMutex initMutex;

  std::vector<cache_aligned::unique_ptr<HighsSplitDeque>> workerDeques;
  cache_aligned::shared_ptr<HighsSplitDeque::WorkerBunk> workerBunk;
  std::atomic<bool> active;

  HighsTask* random_steal_loop(HighsSplitDeque* localDeque) {
    const int numWorkers = workerDeques.size();

    int numTries = 16 * (numWorkers - 1);

    auto tStart = std::chrono::high_resolution_clock::now();

    while (true) {
      for (int s = 0; s < numTries; ++s) {
        HighsTask* task = localDeque->randomSteal();
        if (task) return task;
      }

      if (!workerBunk->haveJobs.load(std::memory_order_relaxed)) break;

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
    // spin until the global executor pointer is set up
    while (!active.load(std::memory_order_acquire))
      HighsSpinMutex::yieldProcessor();
    // now acquire a reference count of the global executor
    auto executor = globalExecutor;
    HighsSplitDeque* localDeque = workerDeques[workerId].get();
    threadLocalWorkerDeque() = localDeque;
    HighsTask* currentTask = workerBunk->waitForNewTask(localDeque);
    while (currentTask != nullptr) {
      localDeque->runStolenTask(currentTask);

      currentTask = random_steal_loop(localDeque);
      if (currentTask != nullptr) continue;

      currentTask = workerBunk->waitForNewTask(localDeque);
    }
  }

 public:
  HighsTaskExecutor(int numThreads) {
    assert(numThreads > 0);
    active.store(false, std::memory_order_relaxed);
    workerDeques.resize(numThreads);
    workerBunk = cache_aligned::make_shared<HighsSplitDeque::WorkerBunk>();
    for (int i = 0; i < numThreads; ++i)
      workerDeques[i] = cache_aligned::make_unique<HighsSplitDeque>(
          workerBunk, workerDeques.data(), i, numThreads);

    threadLocalWorkerDeque() = workerDeques[0].get();
    for (int i = 1; i < numThreads; ++i)
      std::thread([&](int id) { run_worker(id); }, i).detach();
  }

  static HighsSplitDeque* getThisWorkerDeque() {
    return threadLocalWorkerDeque();
  }

  static HighsTaskExecutor* getGlobalTaskExecutor() {
    return globalExecutor.get();
  }

  static int getNumWorkerThreads() {
    return globalExecutor->workerDeques.size();
  }

  static void initialize(int numThreads) {
    if (!globalExecutor) {
      std::lock_guard<HighsSpinMutex> lg{initMutex};
      if (!globalExecutor) {
        globalExecutor =
            cache_aligned::make_shared<HighsTaskExecutor>(numThreads);
        globalExecutor->active.store(true, std::memory_order_release);
      }
    }
  }

  static void shutdown(bool blocking = false) {
    if (globalExecutor) {
      // first spin until every worker has acquired its executor reference
      while (globalExecutor.use_count() != globalExecutor->workerDeques.size())
        HighsSpinMutex::yieldProcessor();
      // set the active flag to false first with release ordering
      globalExecutor->active.store(false, std::memory_order_release);
      // now inject the null task as termination signal to every worker
      for (auto& workerDeque : globalExecutor->workerDeques)
        workerDeque->injectTaskAndNotify(nullptr);
      // finally release the global executor reference
      if (blocking) {
        while (globalExecutor.use_count() != 1)
          HighsSpinMutex::yieldProcessor();
      }

      globalExecutor.reset();
    }
  }

  void sync_stolen_task(HighsSplitDeque* localDeque, HighsTask* stolenTask) {
    if (!localDeque->leapfrogStolenTask(stolenTask)) {
      const int numWorkers = workerDeques.size();
      int numTries = kNumTryFac * (numWorkers - 1);

      auto tStart = std::chrono::high_resolution_clock::now();

      while (true) {
        for (int s = 0; s < numTries; ++s) {
          if (stolenTask->isFinished()) {
            localDeque->popStolen();
            return;
          }
          localDeque->yield();
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
        localDeque->waitForTaskToFinish(stolenTask);
    }

    localDeque->popStolen();
  }
};

#endif
