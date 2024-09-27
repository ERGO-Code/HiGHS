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
/**@file simplex/HighsMipAnalysis.cpp
 * @brief
 */
#include "mip/HighsMipAnalysis.h"

#include <cmath>

#include "mip/MipTimer.h"

const HighsInt check_mip_clock = -4;

void HighsMipAnalysis::setup(const HighsLp& lp, const HighsOptions& options) {
  setupMipTime(options);
}

void HighsMipAnalysis::setupMipTime(const HighsOptions& options) {
  analyse_mip_time = kHighsAnalysisLevelMipTime & options.highs_analysis_level;
  assert(analyse_mip_time);
  if (analyse_mip_time) {
    HighsTimerClock clock;
    clock.timer_pointer_ = timer_;
    MipTimer mip_timer;
    mip_timer.initialiseMipClocks(clock);
    mip_clocks = clock;
  }
}

void HighsMipAnalysis::mipTimerStart(const HighsInt mip_clock
                                     // , const HighsInt thread_id
) {
  if (!analyse_mip_time) return;
  HighsInt highs_timer_clock = mip_clocks.clock_[mip_clock];
  if (highs_timer_clock == check_mip_clock) {
    std::string clock_name = mip_clocks.timer_pointer_->clock_names[check_mip_clock];
    printf("MipTimer: starting clock %d: %s\n", int(check_mip_clock), clock_name.c_str());
  }
  mip_clocks.timer_pointer_->start(highs_timer_clock);
}

void HighsMipAnalysis::mipTimerStop(const HighsInt mip_clock
                                    // , const HighsInt thread_id
) {
  if (!analyse_mip_time) return;
  HighsInt highs_timer_clock = mip_clocks.clock_[mip_clock];
  if (highs_timer_clock == check_mip_clock) {
    std::string clock_name = mip_clocks.timer_pointer_->clock_names[check_mip_clock];
    printf("MipTimer: stopping clock %d: %s\n", int(check_mip_clock), clock_name.c_str());
  }
  mip_clocks.timer_pointer_->stop(highs_timer_clock);
}

bool HighsMipAnalysis::mipTimerRunning(const HighsInt mip_clock
                                     // , const HighsInt thread_id
) {
  if (!analyse_mip_time) return;
  HighsInt highs_timer_clock = mip_clocks.clock_[mip_clock];
  return mip_clocks.timer_pointer_->running(highs_timer_clock);
}

double HighsMipAnalysis::mipTimerRead(const HighsInt mip_clock
                                     // , const HighsInt thread_id
) {
  if (!analyse_mip_time) return;
  HighsInt highs_timer_clock = mip_clocks.clock_[mip_clock];
  return mip_clocks.timer_pointer_->read(highs_timer_clock);
}

void HighsMipAnalysis::reportMipTimer() {
  if (!analyse_mip_time) return;
  //  assert(analyse_mip_time);
  MipTimer mip_timer;
  mip_timer.reportMipCoreClock(mip_clocks);
  mip_timer.reportMipLevel1Clock(mip_clocks);
  mip_timer.reportMipSearchClock(mip_clocks);
  mip_timer.reportMipDiveClock(mip_clocks);
  mip_timer.reportMipPrimalHeuristicsClock(mip_clocks);
}
