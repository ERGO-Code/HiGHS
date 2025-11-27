/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HIGHS_SPIN_MUTEX_H_
#define HIGHS_SPIN_MUTEX_H_

#include <atomic>

#include "HConfig.h"

#ifdef HIGHS_HAVE_MM_PAUSE
#include <immintrin.h>
#else
#include <thread>
#endif

class HighsSpinMutex {
  std::atomic<bool> flag{false};

 public:
  static void yieldProcessor() {
#ifdef HIGHS_HAVE_MM_PAUSE
    _mm_pause();
#else
    // ToDo: See if this is OK on Mac M1
    std::this_thread::yield();
#endif
  }

  bool try_lock() { return !flag.exchange(true, std::memory_order_acquire); }

  void lock() {
    while (true) {
      if (!flag.exchange(true, std::memory_order_acquire)) return;

      while (flag.load(std::memory_order_relaxed)) yieldProcessor();
    }
  }

  void unlock() { flag.store(false, std::memory_order_release); }
};

#endif
