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
#ifndef HIGHS_TASK_H_
#define HIGHS_TASK_H_

#include <atomic>
#include <cassert>
#include <cstring>
#include <type_traits>

#include "parallel/HighsSpinMutex.h"

class HighsSplitDeque;

class HighsTask {
  friend class HighsSplitDeque;

 public:
  enum Constants {
    kMaxTaskSize = 64,
  };

 private:
  struct Metadata {
    std::atomic<HighsSplitDeque*> stealer;
  };

  class CallableBase {
   public:
    virtual void operator()() = 0;
  };

  template <typename F>
  class Callable : public CallableBase {
    F functor;

   public:
    Callable(F&& functor) : functor(std::forward<F>(functor)) {}

    virtual void operator()() override {
      F callFunctor = std::move(functor);
      callFunctor();
    }
  };

  char taskData[kMaxTaskSize - sizeof(Metadata)];
  Metadata metadata;

  CallableBase& getCallable() {
    union {
      CallableBase* callablePtr;
      char* storagePtr;
    } u;

    u.storagePtr = this->taskData;
    return *u.callablePtr;
  }

  /// run task as stealer and return the owner's split deque if the owner is
  /// waiting and needs to be signaled
  HighsSplitDeque* run(HighsSplitDeque* stealer) {
    metadata.stealer.store(stealer, std::memory_order_relaxed);
    getCallable()();
    HighsSplitDeque* waitingOwner = metadata.stealer.exchange(
        reinterpret_cast<HighsSplitDeque*>(uintptr_t{1}),
        std::memory_order_release);

    if (waitingOwner != stealer) return waitingOwner;

    return nullptr;
  }

 public:
  /// initialize the task with given callable type. Task is considered
  /// unfinished after setTaskData
  template <typename F>
  void setTaskData(F&& f) {
    static_assert(sizeof(F) <= sizeof(taskData),
                  "given task type exceeds maximum size allowed for deque\n");
    static_assert(std::is_trivially_destructible<F>::value,
                  "given task type must be trivially destructible\n");
    metadata.stealer.store(nullptr, std::memory_order_relaxed);

    new (taskData) Callable<F>(std::forward<F>(f));

    assert(static_cast<CallableBase*>(reinterpret_cast<Callable<F>*>(
               taskData)) == reinterpret_cast<CallableBase*>(taskData));
  }

  /// run task as owner
  void run() { getCallable()(); }

  /// request notification of the owner when the is finished.
  /// Should be called while the owner holds its wait mutex
  /// and only after getStealerIfUnfinished() has been called.
  /// Returns true if the notification was set and false if it was not set
  /// because the task was finished in the meantime.
  bool requestNotifyWhenFinished(HighsSplitDeque* owner) {
    HighsSplitDeque* stealer =
        metadata.stealer.exchange(owner, std::memory_order_relaxed);

    assert(stealer != nullptr);

    return reinterpret_cast<uintptr_t>(stealer) != 1;
  }

  /// check if task is finished
  bool isFinished() const {
    HighsSplitDeque* stealer = metadata.stealer.load(std::memory_order_acquire);
    return reinterpret_cast<uintptr_t>(stealer) == 1;
  }

  /// get the stealer of a stolen task, or nullptr if the stealer finished
  /// executing the task. Spin waits for the stealer to have started executing
  /// the task if necessary.
  HighsSplitDeque* getStealerIfUnfinished() {
    HighsSplitDeque* stealer = metadata.stealer.load(std::memory_order_acquire);
    if (reinterpret_cast<uintptr_t>(stealer) == 1)
      return nullptr;
    else {
      do {
        // the task has been stolen, but the stealer has not yet started
        // executing the task in this case, yield and check again in a spin
        // loop until the stealer executes the task and becomes visible to
        // this thread
        HighsSpinMutex::yieldProcessor();
        stealer = metadata.stealer.load(std::memory_order_acquire);
      } while (stealer == nullptr);
    }

    if (reinterpret_cast<uintptr_t>(stealer) == 1) return nullptr;

    return stealer;
  }
};

#endif
