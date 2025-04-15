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

#ifndef HIGHS_RACE_TIMER_H_
#define HIGHS_RACE_TIMER_H_

#include <atomic>
#include <limits>

template <typename T>
class HighsRaceTimer {
  std::atomic<T> limit;

 public:
  HighsRaceTimer() : limit(std::numeric_limits<T>::max()) {}

  void decreaseLimit(T newLimit) {
    T current = limit.load(std::memory_order_relaxed);

    while (current > newLimit) {
      if (limit.compare_exchange_weak(current, newLimit,
                                      std::memory_order_relaxed,
                                      std::memory_order_relaxed))
        break;
    }
  }

  bool limitReached(const T currentTime) const {
    return currentTime > limit.load(std::memory_order_relaxed);
  }
};

#endif
