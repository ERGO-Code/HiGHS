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
#ifndef HIGHS_SPLIT_DEQUE_H_
#define HIGHS_SPLIT_DEQUE_H_

#include <array>
#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <mutex>
#include <random>
#include <thread>

#include "parallel/HighsCacheAlign.h"
#include "util/HighsInt.h"
#include "util/HighsRandom.h"

template <uint32_t TaskArraySize, uint32_t MaxTaskSize>
class HighsSplitDeque {
  using cache_aligned = highs::cache_aligned;

  class CallableBase {
   public:
    virtual void operator()() = 0;
  };

  template <typename F>
  class Callable : public CallableBase {
    F functor;

   public:
    Callable(F&& f) : functor(std::forward<F>(f)) {}

    virtual void operator()() {
      F callFunctor = std::move(functor);
      callFunctor();
    }
  };

 public:
  class Task;
  struct GlobalQueue;

 private:
  struct OwnerData {
    cache_aligned::shared_ptr<GlobalQueue> globalQueue = nullptr;
    HighsRandom randgen;
    uint32_t head = 0;
    uint32_t splitCopy = 0;
    int ownerId = -1;
    bool allStolenCopy = true;
  };

  struct WaitForTaskData {
    std::mutex waitMutex;
    std::condition_variable waitCondition;
    std::atomic<Task*> injectedTask{nullptr};

    void injectTaskAndNotify(Task* t) {
      Task* prev = injectedTask.exchange(t, std::memory_order_relaxed);
      if (prev == reinterpret_cast<Task*>(uintptr_t{1})) {
        std::unique_lock<std::mutex> lg(waitMutex);
        waitCondition.notify_one();
      }
    }
  };

  struct StealerData {
    cache_aligned::unique_ptr<WaitForTaskData> waitForTaskData{nullptr};
    std::atomic<uint64_t> ts{0};
    std::atomic_bool allStolen{true};
  };

  using SplitDeque = HighsSplitDeque<TaskArraySize, MaxTaskSize>;

  struct TaskMetadata {
    std::atomic<SplitDeque*> stealer;
  };

  static constexpr uint64_t makeTailSplit(uint32_t tail, uint32_t split) {
    return (uint64_t(tail) << 32) | split;
  }

  static constexpr uint32_t tail(uint64_t tailSplit) { return tailSplit >> 32; }

  static constexpr uint32_t split(uint64_t tailSplit) {
    return static_cast<uint32_t>(tailSplit);
  }

 public:
  class Task {
    TaskMetadata metadata;
    char taskData[MaxTaskSize - sizeof(TaskMetadata)];

   public:
    template <typename F>
    void setTaskData(F&& f) {
      static_assert(sizeof(F) <= sizeof(taskData),
                    "given task type exceeds maximum size allowed for deque\n");
      static_assert(std::is_trivially_destructible<F>::value,
                    "given task type must be trivially destructible\n");
      metadata.stealer.store(nullptr, std::memory_order_relaxed);
      new (taskData) Callable<F>(std::forward<F>(f));
    }

    void run() { reinterpret_cast<CallableBase*>(taskData)->operator()(); }

    void run(SplitDeque* stealer) {
      metadata.stealer.store(stealer, std::memory_order_relaxed);
      reinterpret_cast<CallableBase*>(taskData)->operator()();
      SplitDeque* waitingOwner =
          metadata.stealer.exchange(reinterpret_cast<SplitDeque*>(uintptr_t{1}),
                                    std::memory_order_release);

      if (waitingOwner != stealer) {
        std::unique_lock<std::mutex> lg(
            waitingOwner->stealerData.waitForTaskData->waitMutex);

        waitingOwner->stealerData.waitForTaskData->waitCondition.notify_one();
      }
    }

    bool isFinished() const {
      SplitDeque* stealer = metadata.stealer.load(std::memory_order_acquire);
      return reinterpret_cast<uintptr_t>(stealer) == 1;
    }

    // steal from the stealer until this task is finished or no more work can be
    // stolen from the stealer. If the task is finished then true is returned,
    // otherwise false is returned
    bool leapfrog(SplitDeque* owner) const {
      SplitDeque* stealer = metadata.stealer.load(std::memory_order_acquire);
      if (reinterpret_cast<uintptr_t>(stealer) == 1)
        return true;
      else {
        while (stealer == nullptr) {
          // the task has been stolen, but the stealer has not yet started
          // executing the task in this case, yield and check again in a spin
          // loop until the stealer executes the task and becomes visible to
          // this thread
          std::this_thread::yield();
          stealer = metadata.stealer.load(std::memory_order_acquire);
        }
      }

      while (true) {
        if (reinterpret_cast<uintptr_t>(stealer) == 1) return true;

        Task* t = stealer->stealWithRetryLoop();
        if (t != nullptr) {
          t->run(owner);
          stealer = metadata.stealer.load(std::memory_order_acquire);
          continue;
        }

        stealer = metadata.stealer.load(std::memory_order_acquire);
        break;
      }

      return reinterpret_cast<uintptr_t>(stealer) == 1;
    }

    void waitForTaskToFinish(SplitDeque* owner) {
      std::unique_lock<std::mutex> lg(
          owner->stealerData.waitForTaskData->waitMutex);
      // exchange the value stored in stealer with the pointer to owner
      // so that the stealer will see this pointer instead of nullptr
      // when it stores whether the task is finished. In that case the
      // stealer will additionally acquire the wait mutex and signal the owner
      // thread that the task is finished
      SplitDeque* stealer =
          metadata.stealer.exchange(owner, std::memory_order_relaxed);

      assert(stealer != nullptr);
      // wait on the condition variable until the stealer is finished with the
      // task
      while (reinterpret_cast<uintptr_t>(stealer) != 1) {
        owner->stealerData.waitForTaskData->waitCondition.wait(lg);
        stealer = metadata.stealer.load(std::memory_order_relaxed);
      }
    }
  };

  struct GlobalQueueData {
    std::atomic<SplitDeque*> nextSleeper{nullptr};
    std::atomic<SplitDeque*> nextUnsignaled{nullptr};
    int ownerId;
  };

  struct GlobalQueue {
    static constexpr uint64_t kAbaTagShift = 20;
    static constexpr uint64_t kIndexMask = (uint64_t{1} << kAbaTagShift) - 1;
    cache_aligned::unique_ptr<SplitDeque>* workers;
    alignas(64) std::atomic_uint64_t sleeperStack;
    alignas(64) std::atomic_uint64_t unsignaledStack;

    GlobalQueue(cache_aligned::unique_ptr<SplitDeque>* workers)
        : workers(workers), sleeperStack(0), unsignaledStack(0) {}

    void pushSleeper(SplitDeque* deque) {
      uint64_t stackState = sleeperStack.load(std::memory_order_relaxed);
      uint64_t newStackState;

      do {
        SplitDeque* head = stackState & kIndexMask
                               ? workers[(stackState & kIndexMask) - 1].get()
                               : nullptr;
        deque->globalQueueData.nextSleeper.store(head,
                                                 std::memory_order_relaxed);

        newStackState = (stackState >> kAbaTagShift) + 1;
        newStackState = (newStackState << kAbaTagShift) |
                        uint64_t(deque->readOwnersId() + 1);
      } while (!sleeperStack.compare_exchange_weak(stackState, newStackState,
                                                   std::memory_order_release,
                                                   std::memory_order_relaxed));
    }

    SplitDeque* popSleeper() {
      uint64_t stackState = sleeperStack.load(std::memory_order_relaxed);
      SplitDeque* head;
      uint64_t newStackState;

      do {
        if ((stackState & kIndexMask) == 0) return nullptr;
        head = workers[(stackState & kIndexMask) - 1].get();
        SplitDeque* newHead =
            head->globalQueueData.nextSleeper.load(std::memory_order_relaxed);
        int newHeadId = newHead != nullptr ? newHead->readOwnersId() + 1 : 0;
        newStackState = (stackState >> kAbaTagShift) + 1;
        newStackState = (newStackState << kAbaTagShift) | uint64_t(newHeadId);
      } while (!sleeperStack.compare_exchange_weak(stackState, newStackState,
                                                   std::memory_order_relaxed,
                                                   std::memory_order_acquire));

      head->globalQueueData.nextSleeper.store(nullptr,
                                              std::memory_order_relaxed);

      return head;
    }

    void pushUnsignaled(SplitDeque* deque) {
      uint64_t stackState = unsignaledStack.load(std::memory_order_relaxed);
      uint64_t newStackState;

      do {
        SplitDeque* head = stackState & kIndexMask
                               ? workers[(stackState & kIndexMask) - 1].get()
                               : nullptr;
        deque->globalQueueData.nextUnsignaled.store(head,
                                                    std::memory_order_relaxed);
        newStackState = (stackState >> kAbaTagShift) + 1;
        newStackState = (newStackState << kAbaTagShift) |
                        uint64_t(deque->readOwnersId() + 1);
      } while (!unsignaledStack.compare_exchange_weak(
          stackState, newStackState, std::memory_order_release,
          std::memory_order_relaxed));
    }

    SplitDeque* popUnsignaled() {
      uint64_t stackState = unsignaledStack.load(std::memory_order_relaxed);
      SplitDeque* head;
      uint64_t newStackState;

      do {
        if ((stackState & kIndexMask) == 0) return nullptr;
        head = workers[(stackState & kIndexMask) - 1].get();
        SplitDeque* newHead = head->globalQueueData.nextUnsignaled.load(
            std::memory_order_relaxed);
        int newHeadId = newHead != nullptr ? newHead->readOwnersId() + 1 : 0;
        newStackState = (stackState >> kAbaTagShift) + 1;
        newStackState = (newStackState << kAbaTagShift) | uint64_t(newHeadId);
      } while (!unsignaledStack.compare_exchange_weak(
          stackState, newStackState, std::memory_order_relaxed,
          std::memory_order_acquire));

      head->globalQueueData.nextUnsignaled.store(nullptr,
                                                 std::memory_order_relaxed);

      return head;
    }

    void publishWork(SplitDeque* localDeque) {
      SplitDeque* sleeper = popSleeper();
      bool resetSignal = true;
      while (sleeper) {
        uint32_t t = localDeque->stealWithRetryLoopAndGetTail(false);
        if (t >= localDeque->ownerData.splitCopy) {
          if (localDeque->ownerData.head == localDeque->ownerData.splitCopy) {
            localDeque->ownerData.allStolenCopy = true;
            localDeque->stealerData.allStolen.store(true,
                                                    std::memory_order_relaxed);
          }
          pushSleeper(sleeper);

          SplitDeque* unsignaled;
          while ((unsignaled = popUnsignaled()) != nullptr)
            unsignaled->requestPublishGlobalQueue();

          return;
        } else {
          sleeper->stealerData.waitForTaskData->injectTaskAndNotify(
              &localDeque->taskArray[t]);
        }

        if (t == localDeque->ownerData.splitCopy - 1) {
          if (localDeque->ownerData.head == localDeque->ownerData.splitCopy) {
            localDeque->ownerData.allStolenCopy = true;
            localDeque->stealerData.allStolen.store(true,
                                                    std::memory_order_relaxed);
          }
          return;
        }
        sleeper = popSleeper();
        if (resetSignal && !sleeper) {
          resetSignal = false;
          // reset the publish global request
          localDeque->splitRequest.fetch_and(~kPublishGlobal,
                                             std::memory_order_relaxed);
          // then push this worker as unsignaled
          pushUnsignaled(localDeque);
          // now read again whether there is a new sleeper that might have not
          // seen the flag being reset and try to feed it more work or at least
          // signal it
          sleeper = popSleeper();
        }
      }
    }

    Task* waitForNewTask(SplitDeque* localDeque, int numMicroSecsBeforeSleep) {
      pushSleeper(localDeque);
      SplitDeque* unsignaled;
      while ((unsignaled = popUnsignaled()) != nullptr)
        unsignaled->requestPublishGlobalQueue();

      auto tStart = std::chrono::high_resolution_clock::now();
      int spinIters = 10;
      do {
        for (int i = 0; i < spinIters; ++i) {
          Task* t = localDeque->stealerData.waitForTaskData->injectedTask.load(
              std::memory_order_relaxed);
          if (t != nullptr) {
            localDeque->stealerData.waitForTaskData->injectedTask.store(
                nullptr, std::memory_order_relaxed);
            return t;
          }
          std::this_thread::yield();
        }

        auto numMicroSecs =
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - tStart)
                .count();

        if (numMicroSecs < numMicroSecsBeforeSleep)
          spinIters *= 2;
        else
          break;
      } while (true);

      {
        std::unique_lock<std::mutex> lg(
            localDeque->stealerData.waitForTaskData->waitMutex);

        Task* t =
            localDeque->stealerData.waitForTaskData->injectedTask.exchange(
                reinterpret_cast<Task*>(1), std::memory_order_relaxed);

        if (t != nullptr) {
          localDeque->stealerData.waitForTaskData->injectedTask.store(
              nullptr, std::memory_order_relaxed);
          return t;
        }

        while (true) {
          localDeque->stealerData.waitForTaskData->waitCondition.wait(lg);
          t = localDeque->stealerData.waitForTaskData->injectedTask.load(
              std::memory_order_relaxed);
          if (t != reinterpret_cast<Task*>(1)) {
            localDeque->stealerData.waitForTaskData->injectedTask.store(
                nullptr, std::memory_order_relaxed);
            return t;
          }
        }
      }
    }
  };

 private:
  static_assert(sizeof(OwnerData) <= 64,
                "sizeof(OwnerData) exceeds cache line size");
  static_assert(sizeof(StealerData) <= 64,
                "sizeof(StealerData) exceeds cache line size");
  static_assert(sizeof(GlobalQueueData) <= 64,
                "sizeof(GlobalQueueData) exceeds cache line size");

  alignas(64) OwnerData ownerData;
  alignas(64) std::atomic_uint splitRequest;
  alignas(64) StealerData stealerData;
  alignas(64) GlobalQueueData globalQueueData;
  alignas(64) std::array<Task, TaskArraySize> taskArray;

  enum Request : unsigned int {
    kNone = 0u,
    kSplitRequest = 1u,
    kPublishGlobal = 2u,
  };

  void growShared(bool publishAllTasks) {
    uint32_t newSplit;

    if (publishAllTasks)
      newSplit = std::min(TaskArraySize, ownerData.head);
    else
      newSplit = std::min(TaskArraySize,
                          (ownerData.splitCopy + ownerData.head + 1) / 2);
    assert(newSplit > ownerData.splitCopy);

    // we want to replace the old split point with the new splitPoint
    // but not alter the upper 32 bits of tail.
    // Hence with xor we can xor or the copy of the current split point
    // to set the lower bits to zero and then xor the bits of the new split
    // point to the lower bits that are then zero. First doing the xor of the
    // old and new split point and then doing the xor with the stealerData
    // will not alter the result. Also the upper 32 bits of the xor mask are
    // zero and will therefore not alter the value of tail.
    uint64_t xorMask = ownerData.splitCopy ^ newSplit;
    // since we publish the task data here, we need to use
    // std::memory_order_release which ensures all writes to set up the task
    // are done
    assert((xorMask >> 32) == 0);

    stealerData.ts.fetch_xor(xorMask, std::memory_order_release);
    ownerData.splitCopy = newSplit;
    unsigned int splitRq =
        splitRequest.fetch_and(~kSplitRequest, std::memory_order_relaxed);
    if (splitRq & kPublishGlobal) ownerData.globalQueue->publishWork(this);
  }

  bool shrinkShared() {
    uint64_t ts = stealerData.ts.load(std::memory_order_relaxed);
    uint32_t t = tail(ts);
    uint32_t s = split(ts);

    if (t != s) {
      uint32_t newSplit = (t + s) / 2;
      uint64_t xorMask = newSplit ^ ownerData.splitCopy;
      ownerData.splitCopy = newSplit;
      t = tail(stealerData.ts.fetch_xor(xorMask, std::memory_order_acq_rel));
      if (t != s) {
        if (t > newSplit) {
          newSplit = (t + s) / 2;
          xorMask = newSplit ^ ownerData.splitCopy;
          ownerData.splitCopy = newSplit;
          stealerData.ts.fetch_xor(xorMask, std::memory_order_relaxed);
        }

        return false;
      }
    }

    stealerData.allStolen.store(true, std::memory_order_relaxed);
    ownerData.allStolenCopy = true;
    return true;
  }

 public:
  HighsSplitDeque(const cache_aligned::shared_ptr<GlobalQueue>& globalQueue,
                  int ownerId) {
    ownerData.ownerId = ownerId;
    globalQueueData.ownerId = ownerId;
    std::random_device rd;
    ownerData.randgen.initialise(rd());
    ownerData.globalQueue = globalQueue;

    stealerData.waitForTaskData = cache_aligned::make_unique<WaitForTaskData>();
    splitRequest.store(kNone, std::memory_order_relaxed);
    globalQueue->pushUnsignaled(this);

    assert((reinterpret_cast<uintptr_t>(this) & 63u) == 0);
    static_assert(offsetof(SplitDeque, splitRequest) == 64,
                  "alignas failed to guarantee 64 byte alignment");
    static_assert(offsetof(SplitDeque, stealerData) == 128,
                  "alignas failed to guarantee 64 byte alignment");
    static_assert(offsetof(SplitDeque, globalQueueData) == 192,
                  "alignas failed to guarantee 64 byte alignment");
    static_assert(offsetof(SplitDeque, taskArray) == 256,
                  "alignas failed to guarantee 64 byte alignment");
  }

  template <typename F>
  void push(F&& f) {
    if (ownerData.head >= TaskArraySize) {
      // task queue is full, execute task directly
      if (ownerData.splitCopy < TaskArraySize && !ownerData.allStolenCopy) {
        unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
        if (splitRq) growShared(splitRq & kPublishGlobal);
      }

      ownerData.head += 1;
      f();
      return;
    }

    taskArray[ownerData.head++].setTaskData(std::forward<F>(f));
    if (ownerData.allStolenCopy) {
      assert(ownerData.head > 0);
      stealerData.ts.store(makeTailSplit(ownerData.head - 1, ownerData.head),
                           std::memory_order_release);
      stealerData.allStolen.store(false, std::memory_order_relaxed);
      ownerData.splitCopy = ownerData.head;
      ownerData.allStolenCopy = false;

      unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
      if (splitRq) {
        splitRq =
            splitRequest.fetch_and(~kSplitRequest, std::memory_order_relaxed);
        if (splitRq & kPublishGlobal) ownerData.globalQueue->publishWork(this);
      }
    } else {
      unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
      if (splitRq) growShared(splitRq & kPublishGlobal);
    }
  }

  enum class Status {
    kEmpty,
    kStolen,
    kWork,
    kOverflown,
  };

  std::pair<Status, Task*> pop() {
    if (ownerData.head == 0) return std::make_pair(Status::kEmpty, nullptr);

    if (ownerData.head > TaskArraySize) {
      // task queue was full and the overflown tasks have
      // been directly executed
      ownerData.head -= 1;
      return std::make_pair(Status::kOverflown, nullptr);
    }

    if (ownerData.allStolenCopy)
      return std::make_pair(Status::kStolen, &taskArray[ownerData.head - 1]);

    if (ownerData.splitCopy == ownerData.head) {
      if (shrinkShared())
        return std::make_pair(Status::kStolen, &taskArray[ownerData.head - 1]);
    }

    ownerData.head -= 1;

    if (ownerData.head == 0) {
      if (!ownerData.allStolenCopy) {
        ownerData.allStolenCopy = true;
        stealerData.allStolen.store(true, std::memory_order_relaxed);
      }
    } else if (ownerData.head != ownerData.splitCopy) {
      unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
      if (splitRq) growShared(splitRq & kPublishGlobal);
    }

    return std::make_pair(Status::kWork, &taskArray[ownerData.head]);
  }

  void popStolen() {
    ownerData.head -= 1;
    if (!ownerData.allStolenCopy) {
      ownerData.allStolenCopy = true;
      stealerData.allStolen.store(true, std::memory_order_relaxed);
    }
  }

  Task* steal() {
    if (stealerData.allStolen.load(std::memory_order_relaxed)) return nullptr;

    uint64_t ts = stealerData.ts.load(std::memory_order_relaxed);
    uint32_t t = tail(ts);
    uint32_t s = split(ts);
    if (t < s) {
      if (stealerData.ts.compare_exchange_weak(ts, makeTailSplit(t + 1, s),
                                               std::memory_order_acquire,
                                               std::memory_order_relaxed))
        return &taskArray[t];

      t = tail(ts);
      s = split(ts);
      if (t < s) return nullptr;
    }

    if (t < TaskArraySize && !splitRequest.load(std::memory_order_relaxed))
      splitRequest.fetch_or(kSplitRequest, std::memory_order_relaxed);

    return nullptr;
  }

  Task* stealWithRetryLoop(bool setSplitRq = true) {
    if (stealerData.allStolen.load(std::memory_order_relaxed)) return nullptr;

    uint64_t ts = stealerData.ts.load(std::memory_order_relaxed);
    uint32_t t = tail(ts);
    uint32_t s = split(ts);

    while (t < s) {
      if (stealerData.ts.compare_exchange_weak(ts, makeTailSplit(t + 1, s),
                                               std::memory_order_acquire,
                                               std::memory_order_relaxed))
        return &taskArray[t];

      t = tail(ts);
      s = split(ts);
    }

    if (setSplitRq && t < TaskArraySize &&
        !splitRequest.load(std::memory_order_relaxed))
      splitRequest.fetch_or(kSplitRequest, std::memory_order_relaxed);

    return nullptr;
  }

  uint32_t stealWithRetryLoopAndGetTail(bool setSplitRq = true) {
    if (stealerData.allStolen.load(std::memory_order_relaxed))
      return std::numeric_limits<uint32_t>::max();

    uint64_t ts = stealerData.ts.load(std::memory_order_relaxed);
    uint32_t t = tail(ts);
    uint32_t s = split(ts);

    while (t < s) {
      if (stealerData.ts.compare_exchange_weak(ts, makeTailSplit(t + 1, s),
                                               std::memory_order_acquire,
                                               std::memory_order_relaxed))
        return t;

      t = tail(ts);
      s = split(ts);
    }

    if (setSplitRq && t < TaskArraySize &&
        !splitRequest.load(std::memory_order_relaxed))
      splitRequest.fetch_or(kSplitRequest, std::memory_order_relaxed);

    return s;
  }

  void requestPublishGlobalQueue() {
    // if (!(splitRequest.load(std::memory_order_relaxed) & kPublishGlobal)) {
    splitRequest.store(kSplitRequest | kPublishGlobal,
                       std::memory_order_relaxed);
    //}
  }

  Task* randomSteal(cache_aligned::unique_ptr<SplitDeque>* workers,
                    int numWorkers) {
    HighsInt next = ownerData.randgen.integer(numWorkers - 1);
    next += next >= ownerData.ownerId;
    assert(next != ownerData.ownerId);
    assert(next >= 0);
    assert(next < numWorkers);

    return workers[next]->steal();
  }

  int getOwnerId() const { return ownerData.ownerId; }

  int readOwnersId() const { return globalQueueData.ownerId; }
};

#endif
