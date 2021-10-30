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
#include <thread>

#include "parallel/HighsCacheAlign.h"
#include "parallel/HighsTask.h"
#include "util/HighsInt.h"
#include "util/HighsRandom.h"
#include "parallel/HighsSpinMutex.h"

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
  enum Constants {
    kTaskArraySize = 8192,
  };
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
    std::atomic<HighsTask*> injectedTask{nullptr};
  };

  struct StealerData {
    cache_aligned::unique_ptr<WaitForTaskData> waitForTaskData{nullptr};
    std::atomic<uint64_t> ts{0};
    std::atomic_bool allStolen{true};
  };

  struct TaskMetadata {
    std::atomic<HighsSplitDeque*> stealer;
  };

  static constexpr uint64_t makeTailSplit(uint32_t tail, uint32_t split) {
    return (uint64_t(tail) << 32) | split;
  }

  static constexpr uint32_t tail(uint64_t tailSplit) { return tailSplit >> 32; }

  static constexpr uint32_t split(uint64_t tailSplit) {
    return static_cast<uint32_t>(tailSplit);
  }

  static constexpr int abaTagShift() { return 20; }
  static constexpr uint64_t stackIndexMask() {
    return (uint64_t{1} << abaTagShift()) - 1;
  }

 public:
  struct ThreadGroupData {
    std::atomic<HighsSplitDeque*> nextSleeper{nullptr};
    std::atomic<HighsSplitDeque*> nextUnsignaled{nullptr};
    int workerIndex;
  };

  class SleeperStack {
    std::atomic_uint64_t stackState{0};

   public:
    void push(HighsSplitDeque** const workerArray, HighsSplitDeque* deque) {
      uint64_t state = stackState.load(std::memory_order_relaxed);
      uint64_t newState;

      do {
        HighsSplitDeque* head =
            state & stackIndexMask()
                ? workerArray[(state & stackIndexMask()) - 1]
                : nullptr;
        deque->globalQueueData.nextSleeper.store(head,
                                                 std::memory_order_relaxed);

        newState = (state >> abaTagShift()) + 1;
        newState = (newState << abaTagShift()) |
                   uint64_t(deque->readThreadGroupIndex() + 1);
      } while (!stackState.compare_exchange_weak(state, newState,
                                                 std::memory_order_release,
                                                 std::memory_order_relaxed));
    }

    HighsSplitDeque* pop(HighsSplitDeque** const workerArray) {
      uint64_t state = stackState.load(std::memory_order_relaxed);
      HighsSplitDeque* head;
      uint64_t newState;

      do {
        if ((state & stackIndexMask()) == 0) return nullptr;
        head = workerArray[(state & stackIndexMask()) - 1];
        HighsSplitDeque* newHead =
            head->globalQueueData.nextSleeper.load(std::memory_order_relaxed);
        int newHeadId =
            newHead != nullptr ? newHead->readThreadGroupIndex() + 1 : 0;
        newState = (state >> abaTagShift()) + 1;
        newState = (newState << abaTagShift()) | uint64_t(newHeadId);
      } while (!stackState.compare_exchange_weak(state, newState,
                                                 std::memory_order_relaxed,
                                                 std::memory_order_acquire));

      head->globalQueueData.nextSleeper.store(nullptr,
                                              std::memory_order_relaxed);

      return head;
    }
  };

  class UnsginaledStack {
    std::atomic_uint64_t stackState{0};

   public:
    void push(HighsSplitDeque** const workerArray, HighsSplitDeque* deque) {
      uint64_t state = stackState.load(std::memory_order_relaxed);
      uint64_t newState;

      do {
        HighsSplitDeque* head =
            state & stackIndexMask()
                ? workerArray[(state & stackIndexMask()) - 1]
                : nullptr;
        deque->globalQueueData.nextUnsignaled.store(head,
                                                    std::memory_order_relaxed);

        newState = (state >> abaTagShift()) + 1;
        newState = (newState << abaTagShift()) |
                   uint64_t(deque->readThreadGroupIndex() + 1);
      } while (!stackState.compare_exchange_weak(state, newState,
                                                 std::memory_order_release,
                                                 std::memory_order_relaxed));
    }

    HighsSplitDeque* pop(HighsSplitDeque** const workerArray) {
      uint64_t state = stackState.load(std::memory_order_relaxed);
      HighsSplitDeque* head;
      uint64_t newState;

      do {
        if ((state & stackIndexMask()) == 0) return nullptr;
        head = workerArray[(state & stackIndexMask()) - 1];
        HighsSplitDeque* newHead = head->globalQueueData.nextUnsignaled.load(
            std::memory_order_relaxed);
        int newHeadId = newHead != nullptr ? newHead->readThreadGroupIndex() + 1 : 0;
        newState = (state >> abaTagShift()) + 1;
        newState = (newState << abaTagShift()) | uint64_t(newHeadId);
      } while (!stackState.compare_exchange_weak(state, newState,
                                                 std::memory_order_relaxed,
                                                 std::memory_order_acquire));

      head->globalQueueData.nextUnsignaled.store(nullptr,
                                                 std::memory_order_relaxed);

      return head;
    }
  };

  struct GlobalQueueData {
    std::atomic<HighsSplitDeque*> nextSleeper{nullptr};
    std::atomic<HighsSplitDeque*> nextUnsignaled{nullptr};
    int ownerId;
  };

  struct GlobalQueue {
    static constexpr uint64_t kAbaTagShift = 20;
    static constexpr uint64_t kIndexMask = (uint64_t{1} << kAbaTagShift) - 1;
    cache_aligned::unique_ptr<HighsSplitDeque>* workers;
    alignas(64) std::atomic_uint64_t sleeperStack;
    alignas(64) std::atomic_uint64_t unsignaledStack;

    GlobalQueue(cache_aligned::unique_ptr<HighsSplitDeque>* workers)
        : workers(workers), sleeperStack(0), unsignaledStack(0) {}

    void pushSleeper(HighsSplitDeque* deque) {
      uint64_t stackState = sleeperStack.load(std::memory_order_relaxed);
      uint64_t newStackState;

      do {
        HighsSplitDeque* head =
            stackState & kIndexMask
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

    HighsSplitDeque* popSleeper() {
      uint64_t stackState = sleeperStack.load(std::memory_order_relaxed);
      HighsSplitDeque* head;
      uint64_t newStackState;

      do {
        if ((stackState & kIndexMask) == 0) return nullptr;
        head = workers[(stackState & kIndexMask) - 1].get();
        HighsSplitDeque* newHead =
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

    void pushUnsignaled(HighsSplitDeque* deque) {
      uint64_t stackState = unsignaledStack.load(std::memory_order_relaxed);
      uint64_t newStackState;

      do {
        HighsSplitDeque* head =
            stackState & kIndexMask
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

    HighsSplitDeque* popUnsignaled() {
      uint64_t stackState = unsignaledStack.load(std::memory_order_relaxed);
      HighsSplitDeque* head;
      uint64_t newStackState;

      do {
        if ((stackState & kIndexMask) == 0) return nullptr;
        head = workers[(stackState & kIndexMask) - 1].get();
        HighsSplitDeque* newHead = head->globalQueueData.nextUnsignaled.load(
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

    void publishWork(HighsSplitDeque* localDeque) {
      HighsSplitDeque* sleeper = popSleeper();
      bool resetSignal = true;
      while (sleeper) {
        uint32_t t = localDeque->selfStealAndGetTail();
        if (t == localDeque->ownerData.splitCopy) {
          if (localDeque->ownerData.head == localDeque->ownerData.splitCopy) {
            localDeque->ownerData.allStolenCopy = true;
            localDeque->stealerData.allStolen.store(true,
                                                    std::memory_order_relaxed);
          }
          pushSleeper(sleeper);

          HighsSplitDeque* unsignaled;
          while ((unsignaled = popUnsignaled()) != nullptr)
            unsignaled->requestPublishGlobalQueue();

          return;
        } else {
          sleeper->injectTaskAndNotify(&localDeque->taskArray[t]);
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

    HighsTask* waitForNewTask(HighsSplitDeque* localDeque,
                              int numMicroSecsBeforeSleep) {
      pushSleeper(localDeque);
      HighsSplitDeque* unsignaled;
      while ((unsignaled = popUnsignaled()) != nullptr)
        unsignaled->requestPublishGlobalQueue();

      auto tStart = std::chrono::high_resolution_clock::now();
      int spinIters = 10;
      do {
        for (int i = 0; i < spinIters; ++i) {
          HighsTask* t =
              localDeque->stealerData.waitForTaskData->injectedTask.load(
                  std::memory_order_relaxed);
          if (t != nullptr) {
            localDeque->stealerData.waitForTaskData->injectedTask.store(
                nullptr, std::memory_order_relaxed);
            return t;
          }

          HighsSpinMutex::yieldProcessor();
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

        HighsTask* t =
            localDeque->stealerData.waitForTaskData->injectedTask.exchange(
                reinterpret_cast<HighsTask*>(1), std::memory_order_relaxed);

        if (t != nullptr) {
          localDeque->stealerData.waitForTaskData->injectedTask.store(
              nullptr, std::memory_order_relaxed);
          return t;
        }

        while (true) {
          localDeque->stealerData.waitForTaskData->waitCondition.wait(lg);
          t = localDeque->stealerData.waitForTaskData->injectedTask.load(
              std::memory_order_relaxed);
          if (t != reinterpret_cast<HighsTask*>(1)) {
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

  OwnerData ownerData;
  alignas(64) std::atomic_uint splitRequest;
  alignas(64) StealerData stealerData;
  alignas(64) GlobalQueueData globalQueueData;
  alignas(64) std::array<HighsTask, kTaskArraySize> taskArray;

  enum Request : unsigned int {
    kNone = 0u,
    kSplitRequest = 0x00ffu,
    kPublishGlobal = 0xff00u,
  };

  void growShared() {
    unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
    if (!splitRq) return;

    uint32_t newSplit;

    if (splitRq & kPublishGlobal)
      newSplit = std::min(uint32_t{kTaskArraySize}, ownerData.head);
    else
      newSplit = std::min(uint32_t{kTaskArraySize},
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
    splitRq = splitRequest.fetch_and(~kSplitRequest, std::memory_order_relaxed);
    if (splitRq & kPublishGlobal) ownerData.globalQueue->publishWork(this);
  }

  bool shrinkShared() {
    uint32_t t = tail(stealerData.ts.load(std::memory_order_relaxed));
    uint32_t s = ownerData.splitCopy;

    if (t != s) {
      ownerData.splitCopy = (t + s) / 2;
      t = tail(stealerData.ts.fetch_add(uint64_t{ownerData.splitCopy} - s,
                                        std::memory_order_acq_rel));
      if (t != s) {
        if (t > ownerData.splitCopy) {
          ownerData.splitCopy = (t + s) / 2;
          stealerData.ts.store(makeTailSplit(t, ownerData.splitCopy),
                               std::memory_order_relaxed);
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
    ownerData.randgen.initialise(ownerId);
    ownerData.globalQueue = globalQueue;

    stealerData.waitForTaskData = cache_aligned::make_unique<WaitForTaskData>();
    splitRequest.store(kNone, std::memory_order_relaxed);
    globalQueue->pushUnsignaled(this);

    assert((reinterpret_cast<uintptr_t>(this) & 63u) == 0);
    static_assert(offsetof(HighsSplitDeque, splitRequest) == 64,
                  "alignas failed to guarantee 64 byte alignment");
    static_assert(offsetof(HighsSplitDeque, stealerData) == 128,
                  "alignas failed to guarantee 64 byte alignment");
    static_assert(offsetof(HighsSplitDeque, globalQueueData) == 192,
                  "alignas failed to guarantee 64 byte alignment");
    static_assert(offsetof(HighsSplitDeque, taskArray) == 256,
                  "alignas failed to guarantee 64 byte alignment");
  }

  template <typename F>
  void push(F&& f) {
    if (ownerData.head >= kTaskArraySize) {
      // task queue is full, execute task directly
      if (ownerData.splitCopy < kTaskArraySize && !ownerData.allStolenCopy)
        growShared();

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
    } else
      growShared();
  }

  enum class Status {
    kEmpty,
    kStolen,
    kWork,
    kOverflown,
  };

  std::pair<Status, HighsTask*> pop() {
    if (ownerData.head == 0) return std::make_pair(Status::kEmpty, nullptr);

    if (ownerData.head > kTaskArraySize) {
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
    } else if (ownerData.head != ownerData.splitCopy)
      growShared();

    return std::make_pair(Status::kWork, &taskArray[ownerData.head]);
  }

  void popStolen() {
    ownerData.head -= 1;
    if (!ownerData.allStolenCopy) {
      ownerData.allStolenCopy = true;
      stealerData.allStolen.store(true, std::memory_order_relaxed);
    }
  }

  HighsTask* steal() {
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

    if (t < kTaskArraySize && !splitRequest.load(std::memory_order_relaxed))
      splitRequest.fetch_or(kSplitRequest, std::memory_order_relaxed);

    return nullptr;
  }

  HighsTask* stealWithRetryLoop() {
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

    if (t < kTaskArraySize && !splitRequest.load(std::memory_order_relaxed))
      splitRequest.fetch_or(kSplitRequest, std::memory_order_relaxed);

    return nullptr;
  }

  uint32_t selfStealAndGetTail() {
    if (ownerData.allStolenCopy) return ownerData.splitCopy;

    // when we steal from ourself we can simply do a fetch_add predictively
    // instead of a cas loop. If the tail we read like this ends up to be
    // above already equal to the splitPoint then we correct it with a simple
    // store. When tail > split instead of tail == split no wrong result can
    // occur as long as we know that the task at taskArray[split] is not
    // actually considered to be stolen and tail is corrected before the owner
    // enters shrinkShared.

    uint32_t t = tail(stealerData.ts.fetch_add(makeTailSplit(1, 0),
                                               std::memory_order_relaxed));

    if (t == ownerData.splitCopy)
      stealerData.ts.store(makeTailSplit(t, ownerData.splitCopy),
                           std::memory_order_relaxed);

    return t;
  }

  void requestPublishGlobalQueue() {
    splitRequest.store(kSplitRequest | kPublishGlobal,
                       std::memory_order_relaxed);
  }

  HighsTask* randomSteal(cache_aligned::unique_ptr<HighsSplitDeque>* workers,
                         int numWorkers) {
    HighsInt next = ownerData.randgen.integer(numWorkers - 1);
    next += next >= ownerData.ownerId;
    assert(next != ownerData.ownerId);
    assert(next >= 0);
    assert(next < numWorkers);

    return workers[next]->steal();
  }

  void injectTaskAndNotify(HighsTask* t) {
    HighsTask* prev = stealerData.waitForTaskData->injectedTask.exchange(
        t, std::memory_order_relaxed);
    if (prev == reinterpret_cast<HighsTask*>(uintptr_t{1})) {
      std::unique_lock<std::mutex> lg(stealerData.waitForTaskData->waitMutex);
      stealerData.waitForTaskData->waitCondition.notify_one();
    }
  }

  void notify() {
    std::unique_lock<std::mutex> lg(stealerData.waitForTaskData->waitMutex);
    stealerData.waitForTaskData->waitCondition.notify_one();
  }

  void runStolenTask(HighsTask* task) {
    HighsSplitDeque* owner = task->run(this);
    if (owner) owner->notify();
  }

  // steal from the stealer until this task is finished or no more work can be
  // stolen from the stealer. If the task is finished then true is returned,
  // otherwise false is returned
  bool leapfrogStolenTask(HighsTask* task) {
    HighsSplitDeque* stealer = task->getStealerIfUnfinished();

    if (stealer == nullptr) return true;

    do {
      HighsTask* t = stealer->stealWithRetryLoop();
      if (t == nullptr) break;
      runStolenTask(t);
    } while (!task->isFinished());

    return task->isFinished();
  }

  void waitForTaskToFinish(HighsTask* t) {
    std::unique_lock<std::mutex> lg(stealerData.waitForTaskData->waitMutex);
    // exchange the value stored in stealer with the pointer to owner
    // so that the stealer will see this pointer instead of nullptr
    // when it stores whether the task is finished. In that case the
    // stealer will additionally acquire the wait mutex and signal the owner
    // thread that the task is finished

    if (!t->requestNotifyWhenFinished(this)) return;

    // wait on the condition variable until the stealer is finished with the
    // task
    do {
      stealerData.waitForTaskData->waitCondition.wait(lg);
    } while (!t->isFinished());
  }

  int getOwnerId() const { return ownerData.ownerId; }

  int readOwnersId() const { return globalQueueData.ownerId; }

  int readThreadGroupIndex() const { return 0; } //return threadGroupData.workerIndex; }
};

#endif
