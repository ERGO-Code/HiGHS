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

#include <iostream>

#include <cassert>
#include <chrono>
#include <condition_variable>
#include <memory>
#include <thread>
#include <vector>

#include "parallel/HighsCacheAlign.h"
#include "parallel/HighsSchedulerConstants.h"
#include "parallel/HighsSplitDeque.h"
#include "util/HighsInt.h"
#include "util/HighsRandom.h"

class HighsTaskExecutor {
 public:
  using cache_aligned = highs::cache_aligned;
  struct ExecutorHandle {
    HighsTaskExecutor* ptr{nullptr};
    bool isMain{false};

    void dispose();
    ~ExecutorHandle() { dispose(); }
  };

 private:
#ifdef _WIN32
  static HighsSplitDeque*& threadLocalWorkerDeque();
  static ExecutorHandle& threadLocalExecutorHandle();
#else
  static thread_local HighsSplitDeque* threadLocalWorkerDequePtr;
  static thread_local ExecutorHandle globalExecutorHandle;

  static HighsSplitDeque*& threadLocalWorkerDeque() {
    return threadLocalWorkerDequePtr;
  }

  static ExecutorHandle& threadLocalExecutorHandle() {
    return globalExecutorHandle;
  }
#endif

  std::atomic<int> referenceCount;

  cache_aligned::shared_ptr<HighsSplitDeque::WorkerBunk> workerBunk;
  std::vector<cache_aligned::unique_ptr<HighsSplitDeque>> workerDeques;
  std::vector<std::thread> workerThreads;

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

      if (numMicroSecs < HighsSchedulerConstants::kMicroSecsBeforeGlobalSync)
        numTries *= 2;
      else
        break;
    }

    return nullptr;
  }

  static void run_worker(int workerId, HighsTaskExecutor* ptr) {
    auto &executorHandle = threadLocalExecutorHandle();
    executorHandle.ptr = ptr;

    // check if main thread has shutdown before thread has started
    // if (ptr->mainWorkerId.load() != std::thread::id()) {

    // ???
    if (!(executorHandle.isMain)) {

      HighsSplitDeque* localDeque = ptr->workerDeques[workerId].get();
      threadLocalWorkerDeque() = localDeque;

      HighsTask* currentTask = ptr->workerBunk->waitForNewTask(localDeque);
      while (currentTask != nullptr) {
        localDeque->runStolenTask(currentTask);

        currentTask = ptr->random_steal_loop(localDeque);
        if (currentTask != nullptr) continue;

        currentTask = ptr->workerBunk->waitForNewTask(localDeque);
      }
    } else {
      std::cout << " WHAT NO"<< std::endl;
    }

    threadLocalExecutorHandle().dispose();
  }

 public:
  HighsTaskExecutor(int numThreads) {
    assert(numThreads > 0);

    // mainWorkerId.store(std::this_thread::get_id());
    // isMainThread = true;

    workerDeques.resize(numThreads);
    workerBunk = cache_aligned::make_shared<HighsSplitDeque::WorkerBunk>();
    for (int i = 0; i < numThreads; ++i)
      workerDeques[i] = cache_aligned::make_unique<HighsSplitDeque>(
          workerBunk, workerDeques.data(), i, numThreads);

    threadLocalWorkerDeque() = workerDeques[0].get();
    workerThreads.reserve(numThreads - 1);
    referenceCount.store(numThreads);

    for (int i = 1, numThreads = workerDeques.size(); i < numThreads; ++i) {
      workerThreads.emplace_back(
          std::move(std::thread(&HighsTaskExecutor::run_worker, i, this)));
    }
  }

  void stopWorkerThreads(bool blocking = false) {
    // auto id = mainWorkerId.exchange(std::thread::id());
    // if (id == std::thread::id()) return;  // already been called
    // ??? 

    // now inject the null task as termination signal to every worker
    for (auto& workerDeque : workerDeques) {
      workerDeque->injectTaskAndNotify(nullptr);
    }

    // only block if called on main thread, otherwise deadlock may occur
    // if (blocking && std::this_thread::get_id() == id) {
    // threadLocalExecutorHandle();
    if (blocking && threadLocalExecutorHandle().isMain) {
      for (auto& workerThread : workerThreads) {
        workerThread.join();
      }
    } 
    else {
      for (auto& workerThread : workerThreads) {
        // do not detach if shutdown() will be called from the driver code
        // after we are done with the highs library
        //workerThread.detach();
        workerThread.detach();
      }
    }
  }

  static HighsSplitDeque* getThisWorkerDeque() {
    return threadLocalWorkerDeque();
  }

  static int getNumWorkerThreads() {
    return threadLocalWorkerDeque()->getNumWorkers();
  }

  static void initialize(int numThreads) {
    auto& executorHandle = threadLocalExecutorHandle();
    if (executorHandle.ptr == nullptr) {
      executorHandle.ptr = new (cache_aligned::alloc(sizeof(HighsTaskExecutor))) HighsTaskExecutor(numThreads);
      executorHandle.isMain = true;
    }
  }

  // can be called on main or worker threads
  // blocking ignored unless called on main thread
  static void shutdown(bool blocking = false) {
    auto& executorHandle = threadLocalExecutorHandle();

    if (executorHandle.ptr != nullptr) {
      executorHandle.ptr->stopWorkerThreads(blocking);
      executorHandle.dispose();
    }
  }

  static void sync_stolen_task(HighsSplitDeque* localDeque,
                               HighsTask* stolenTask) {
    HighsSplitDeque* stealer;
    if (!localDeque->leapfrogStolenTask(stolenTask, stealer)) {
      const int numWorkers = localDeque->getNumWorkers();
      int numTries = HighsSchedulerConstants::kNumTryFac * (numWorkers - 1);

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

        if (numMicroSecs < HighsSchedulerConstants::kMicroSecsBeforeSleep)
          numTries *= 2;
        else
          break;
      }

      if (!stolenTask->isFinished())
        localDeque->waitForTaskToFinish(stolenTask, stealer);
    }

    localDeque->popStolen();
  }
};

#endif
