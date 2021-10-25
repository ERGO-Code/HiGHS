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

  struct OwnerData {
    HighsRandom randgen;
    uint32_t head = 0;
    uint32_t splitCopy = 0;
    int ownerId = -1;
    bool allStolenCopy = true;
  };

  struct WaitForTaskData {
    std::mutex waitMutex;
    std::condition_variable waitCondition;
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
      SplitDeque* waitingOwner = metadata.stealer.exchange(
          reinterpret_cast<SplitDeque*>(1), std::memory_order_release);

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

  struct GlobalQueue;

  struct GlobalQueueData {
    SplitDeque* next = nullptr;
    SplitDeque* prev = nullptr;
    SplitDeque** headPtr = nullptr;
    SplitDeque** tailPtr = nullptr;
    cache_aligned::shared_ptr<GlobalQueue> globalQueue = nullptr;
  };

  struct GlobalQueue {
    std::mutex mtx;
    std::condition_variable queueNonEmpty;

    SplitDeque* queueHead = nullptr;
    SplitDeque* queueTail = nullptr;

    SplitDeque* activeHead = nullptr;
    SplitDeque* activeTail = nullptr;

    void push(SplitDeque* localDeque) {
      std::unique_lock<std::mutex> lg(mtx);
      localDeque->unlinkGlobalQueue();

      // since we needed to go to the global queue we know that there are
      // sleeping workers that do not have work yet. Hence, we directly store a
      // split request for ourself so that it gets less likely we need to
      // acquire the lock of the global queue again.
      localDeque->splitRequest.store(kSplitRequest, std::memory_order_relaxed);
      localDeque->linkGlobalQueue(queueHead, queueTail);
      queueNonEmpty.notify_one();
    }

    Task* waitForNewTask(SplitDeque* localDeque) {
      std::unique_lock<std::mutex> lg(mtx);
      localDeque->unlinkGlobalQueue();
      localDeque->splitRequest.store(kNone, std::memory_order_relaxed);
      while (true) {
        while (queueTail != nullptr) {
          SplitDeque* q = queueTail;
          if (q != localDeque) {
            Task* t = q->stealWithRetryLoop(false);
            if (t != nullptr) {
              localDeque->linkGlobalQueue(queueHead, queueTail);
              queueNonEmpty.notify_one();
              return t;
            }

            q->unlinkGlobalQueue();
            q->linkGlobalQueue(activeHead, activeTail);
          } else
            q->unlinkGlobalQueue();
        }

        // set request on all active workers to publish their work to the global
        // queue
        for (SplitDeque* q = activeHead; q != nullptr;
             q = q->globalQueueData.next)
          q->requestPublishGlobalQueue();

        // now wait until the global queue becomes non-empty
        while (queueTail == nullptr) queueNonEmpty.wait(lg);
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
    // old and new split point and then doing the xor with the stealerData will
    // not alter the result. Also the upper 32 bits of the xor mask are zero and
    // will therefore not alter the value of tail.
    uint64_t xorMask = ownerData.splitCopy ^ newSplit;
    // since we publish the task data here, we need to use
    // std::memory_order_release which ensures all writes to set up the task are
    // done
    assert((xorMask >> 32) == 0);

    stealerData.ts.fetch_xor(xorMask, std::memory_order_release);
    ownerData.splitCopy = newSplit;
    unsigned int splitRq =
        splitRequest.fetch_and(~kSplitRequest, std::memory_order_relaxed);
    if (splitRq & kPublishGlobal) globalQueueData.globalQueue->push(this);
  }

  void unlinkGlobalQueue() {
    assert(globalQueueData.tailPtr != nullptr);
    assert(globalQueueData.headPtr != nullptr);
    SplitDeque* tmpNext = globalQueueData.next;
    SplitDeque* tmpPrev = globalQueueData.prev;

    globalQueueData.next = nullptr;
    globalQueueData.prev = nullptr;

    if (tmpNext != nullptr)
      tmpNext->globalQueueData.prev = tmpPrev;
    else
      *globalQueueData.tailPtr = tmpPrev;

    if (tmpPrev != nullptr)
      tmpPrev->globalQueueData.next = tmpNext;
    else
      *globalQueueData.headPtr = tmpNext;
  }

  void linkGlobalQueue(SplitDeque*& head, SplitDeque*& tail) {
    globalQueueData.headPtr = &head;
    globalQueueData.tailPtr = &tail;

    globalQueueData.next = head;
    globalQueueData.prev = nullptr;
    head = this;

    if (globalQueueData.next == nullptr)
      tail = this;
    else
      globalQueueData.next->globalQueueData.prev = this;
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
    std::random_device rd;
    ownerData.randgen.initialise(rd());

    stealerData.waitForTaskData = cache_aligned::make_unique<WaitForTaskData>();
    globalQueueData.globalQueue = globalQueue;
    splitRequest.store(kNone, std::memory_order_relaxed);
    linkGlobalQueue(globalQueue->activeHead, globalQueue->activeTail);

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

      unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
      if (splitRq) {
        splitRq =
            splitRequest.fetch_and(~kSplitRequest, std::memory_order_relaxed);
        if (splitRq & kPublishGlobal) globalQueueData.globalQueue->push(this);
      }
      ownerData.splitCopy = ownerData.head;
      ownerData.allStolenCopy = false;
    } else {
      unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
      if (splitRq) growShared(splitRq & kPublishGlobal);
    }
  }

  template <typename F>
  void pushPrivate(F&& f) {
    if (ownerData.head >= TaskArraySize) {
      ownerData.head += 1;
      f();
      return;
    }

    taskArray[ownerData.head++].setTaskData(std::forward<F>(f));
  }

  void publishTasks(HighsInt numTasks) {
    if (ownerData.allStolenCopy) {
      assert(ownerData.head >= numTasks);
      stealerData.ts.store(
          makeTailSplit(ownerData.head - numTasks, ownerData.head),
          std::memory_order_release);
      stealerData.allStolen.store(false, std::memory_order_relaxed);

      unsigned int splitRq = splitRequest.load(std::memory_order_relaxed);
      if (splitRq) {
        splitRq =
            splitRequest.fetch_and(~kSplitRequest, std::memory_order_relaxed);
        if (splitRq & kPublishGlobal) globalQueueData.globalQueue->push(this);
      }
      ownerData.splitCopy = ownerData.head;
      ownerData.allStolenCopy = false;
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

    if (ownerData.head != ownerData.splitCopy) {
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

  void requestPublishGlobalQueue() {
    if (!(splitRequest.load(std::memory_order_relaxed) & kPublishGlobal)) {
      splitRequest.store(kSplitRequest | kPublishGlobal,
                         std::memory_order_relaxed);
    }
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
};

#endif
