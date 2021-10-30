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

#ifndef HIGHS_SPIN_MUTEX_H_
#define HIGHS_SPIN_MUTEX_H_

#include <atomic>

#include "HConfig.h"

#ifdef HIGHS_MM_PAUSE_REQUIRES_IMMINTRIN
#include "immintrin.h"
#elif !defined(HIGHS_HAVE_MM_PAUSE)
#include <thread>
#endif

class HighsSpinMutex {
  std::atomic_bool flag{true};

 public:
  static void yieldProcessor() {
#ifdef HIGHS_HAVE_MM_PAUSE
    _mm_pause();
#else
    std::this_thread::yield();
#endif
  }

  void lock() noexcept {
    while (true) {
      if (!flag.exchange(true, std::memory_order_acquire)) return;

      while (flag.load(std::memory_order_relaxed)) yieldProcessor();
    }
  }

  void unlock() noexcept { flag.store(false, std::memory_order_release); }
};

#endif