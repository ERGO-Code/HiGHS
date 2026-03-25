/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveTimer.h
 * @brief Indices of presolve iClocks
 */
#ifndef PRESOLVE_PRESOLVETIMER_H_
#define PRESOLVE_PRESOLVETIMER_H_

// Clocks for profiling presolve
enum iClockPresolve {
  kPresolveClockTotal = 0,
  kPresolveClockPresolve,
  kNumPresolveClock  //!< Number of PRESOLVE clocks
};

static const double kPresolveClockTolerancePercentReport = 0.1;

class PresolveTimer {
 public:
  void initialisePresolveClocks(HighsTimerClock& presolve_timer_clock) {
    HighsTimer* timer_pointer = presolve_timer_clock.timer_pointer_;
    std::vector<HighsInt>& clock = presolve_timer_clock.clock_;

    clock.resize(kNumPresolveClock);
    clock[kPresolveClockTotal] = 0;
    clock[kPresolveClockPresolve] = timer_pointer->clock_def("Presolve");
  };

  bool reportPresolveClockList(const char* grepStamp,
                          const std::vector<HighsInt> presolve_clock_list,
                          const HighsTimerClock& presolve_timer_clock,
                          const HighsInt kPresolveClockIdeal = kPresolveClockTotal,
                          const double tolerance_percent_report_ = -1) {
    HighsTimer* timer_pointer = presolve_timer_clock.timer_pointer_;
    if (!timer_pointer->printf_flag) return false;
    const std::vector<HighsInt>& clock = presolve_timer_clock.clock_;
    HighsInt presolve_clock_list_size = presolve_clock_list.size();
    std::vector<HighsInt> clockList;
    clockList.resize(presolve_clock_list_size);
    for (HighsInt en = 0; en < presolve_clock_list_size; en++) {
      clockList[en] = clock[presolve_clock_list[en]];
    }
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kPresolveClockIdeal]];
    const double tolerance_percent_report =
        tolerance_percent_report_ >= 0 ? tolerance_percent_report_ : 1e-8;
    return timer_pointer->reportOnTolerance(
        grepStamp, clockList, ideal_sum_time, tolerance_percent_report);
  };

  void reportPresolveCoreClock(const HighsTimerClock& presolve_timer_clock) {
    const std::vector<HighsInt> presolve_clock_list{
        kPresolveClockPresolve};
    reportPresolveClockList("PresolveCore_", presolve_clock_list, presolve_timer_clock,
                       kPresolveClockTotal);
  };
};

#endif /* PRESOLVE_PRESOLVETIMER_H_ */
