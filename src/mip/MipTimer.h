/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/MipTimer.h
 * @brief Indices of mip iClocks
 */
#ifndef MIP_MIPTIMER_H_
#define MIP_MIPTIMER_H_

// Clocks for profiling the MIP dual mip solver
enum iClockMip {
  kMipClockTotal = 0,
  kMipClockPresolve,
  kMipClockSolve,
  kMipClockPostsolve,
  kMipClockInit,
  kMipClockRunPresolve,
  kMipClockRunSetup,
  kMipClockEvaluateRootNode,
  kMipClockPerformAging,
  kMipClockEvaluateNode,
  kMipClockCurrentNodePruned,
  kNumMipClock  //!< Number of MIP clocks
};

class MipTimer {
 public:
  void initialiseMipClocks(HighsTimerClock& mip_timer_clock) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    clock.resize(kNumMipClock);
    clock[kMipClockTotal] = timer_pointer->total_clock;
    clock[kMipClockPresolve] = timer_pointer->presolve_clock;
    clock[kMipClockSolve] = timer_pointer->solve_clock;
    clock[kMipClockPostsolve] = timer_pointer->postsolve_clock;
    clock[kMipClockInit] = timer_pointer->clock_def("Initialise");
    clock[kMipClockRunPresolve] = timer_pointer->clock_def("Run presolve");
    clock[kMipClockRunSetup] = timer_pointer->clock_def("Run setup");
    clock[kMipClockEvaluateRootNode] = timer_pointer->clock_def("Evaluate root node");
    clock[kMipClockPerformAging] = timer_pointer->clock_def("Perform aging");
    clock[kMipClockEvaluateNode] = timer_pointer->clock_def("Evaluate node");
    clock[kMipClockCurrentNodePruned] = timer_pointer->clock_def("Current node pruned");
    //    clock[] = timer_pointer->clock_def("");
  }

  bool reportMipClockList(const char* grepStamp,
                          const std::vector<HighsInt> mip_clock_list,
                          const HighsTimerClock& mip_timer_clock,
                          const double tolerance_percent_report_ = -1) {
    HighsTimer* timer_pointer = mip_timer_clock.timer_pointer_;
    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    HighsInt mip_clock_list_size = mip_clock_list.size();
    std::vector<HighsInt> clockList;
    clockList.resize(mip_clock_list_size);
    for (HighsInt en = 0; en < mip_clock_list_size; en++) {
      clockList[en] = clock[mip_clock_list[en]];
    }
    const double ideal_sum_time =
        timer_pointer->clock_time[clock[kMipClockTotal]];
    const double tolerance_percent_report =
        tolerance_percent_report_ >= 0 ? tolerance_percent_report_ : 1e-8;
    return timer_pointer->reportOnTolerance(
        grepStamp, clockList, ideal_sum_time, tolerance_percent_report);
  };

  void reportMipCoreClock(const HighsTimerClock& mip_timer_clock) {
    //    const std::vector<HighsInt>& clock = mip_timer_clock.clock_;
    const std::vector<HighsInt> mip_clock_list{
        kMipClockPresolve, kMipClockSolve, kMipClockPostsolve};
    reportMipClockList("MipCore", mip_clock_list, mip_timer_clock);
  };

  void reportMipLevel1Clock(const HighsTimerClock& mip_timer_clock) {
    const std::vector<HighsInt> mip_clock_list{
      kMipClockInit,
      kMipClockRunPresolve,
      kMipClockRunSetup,
      kMipClockEvaluateRootNode,
      kMipClockPerformAging,
      kMipClockEvaluateNode,
      kMipClockCurrentNodePruned};
    reportMipClockList("MipLevel1", mip_clock_list, mip_timer_clock);
  };
};

#endif /* MIP_MIPTIMER_H_ */
