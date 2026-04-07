/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsMipAnalysis.cpp
 * @brief
 */
#include "mip/HighsMipAnalysis.h"

#include <cmath>

#include "mip/HighsSeparator.h"
#include "mip/MipTimer.h"
#include "parallel/HighsParallel.h"
#include "util/HighsUtils.h"

const HighsInt check_mip_clock = -4;

void HighsMipAnalysis::setupMipTime(const HighsOptions& options) {
  analyse_mip_time = kHighsAnalysisLevelMipTime & options.highs_analysis_level;
  if (analyse_mip_time) {
    // Set up the thread clocks
    HighsInt max_threads = highs::parallel::num_threads();
    thread_mip_clocks.clear();
    for (HighsInt i = 0; i < max_threads; i++) {
      HighsTimerClock clock;
      clock.timer_pointer_ = timer_;
      thread_mip_clocks.push_back(clock);
      thread_mip_clocks_.push_back(clock);
      thread_submip_clocks_.push_back(clock);
    }
    MipTimer mip_timer;
    // Some sub-solver timings are extracted from the MIP clocks, and
    // are assumed to be specific global clock IDs, but this no longer
    // happens with mult-threaded clocks, as clock IDs for each thread
    // have an offset due to clocks defined for earlier threads.
    HighsInt thread_mip_clock_offset = 0;
    for (HighsTimerClock& clock : thread_mip_clocks) {
      mip_timer.initialiseMipClocks(clock, thread_mip_clock_offset);
      thread_mip_clock_offset += kNumThreadMipClock;
    }
    for (HighsTimerClock& clock : thread_mip_clocks_) {
      mip_timer.initialiseMipClocks(clock, thread_mip_clock_offset);
      thread_mip_clock_offset += kNumThreadMipClock;
    }
    for (HighsTimerClock& clock : thread_submip_clocks_) {
      mip_timer.initialiseMipClocks(clock, thread_mip_clock_offset);
      thread_mip_clock_offset += kNumThreadMipClock;
    }
    mip_clocks = thread_mip_clocks[0];
    sepa_name_clock.push_back(
        std::make_pair(kImplboundSepaString, kMipClockImplboundSepa));
    sepa_name_clock.push_back(
        std::make_pair(kCliqueSepaString, kMipClockCliqueSepa));
    sepa_name_clock.push_back(
        std::make_pair(kTableauSepaString, kMipClockTableauSepa));
    sepa_name_clock.push_back(
        std::make_pair(kPathAggrSepaString, kMipClockPathAggrSepa));
    sepa_name_clock.push_back(
        std::make_pair(kModKSepaString, kMipClockModKSepa));
  }
}

void HighsMipAnalysis::mipTimerStart(const HighsInt mip_clock
                                     // , const HighsInt thread_id
) const {
  if (!analyse_mip_time) return;
  HighsInt local_thread_id = 0;  // highs::parallel::thread_num();
  HighsInt highs_timer_clock =
      thread_mip_clocks[local_thread_id].clock_[mip_clock];
  if (local_thread_id > 0) {
    printf(
        "mipTimerStart with MIP clock %2d and thread %2d for HiGHS clock %4d "
        "(%s)\n",
        int(mip_clock), int(local_thread_id), int(highs_timer_clock),
        thread_mip_clocks[local_thread_id]
            .timer_pointer_->clock_names[highs_timer_clock]
            .c_str());
  }

  thread_mip_clocks[local_thread_id].timer_pointer_->start(highs_timer_clock);
}

void HighsMipAnalysis::mipTimerStop(const HighsInt mip_clock
                                    // , const HighsInt thread_id
) const {
  if (!analyse_mip_time) return;
  HighsInt local_thread_id = 0;  // highs::parallel::thread_num();
  HighsInt highs_timer_clock =
      thread_mip_clocks[local_thread_id].clock_[mip_clock];
  thread_mip_clocks[local_thread_id].timer_pointer_->stop(highs_timer_clock);
}

bool HighsMipAnalysis::mipTimerRunning(const HighsInt mip_clock
                                       // , const HighsInt thread_id
) const {
  if (!analyse_mip_time) return false;
  HighsInt local_thread_id = 0;  // highs::parallel::thread_num();
  HighsInt highs_timer_clock =
      thread_mip_clocks[local_thread_id].clock_[mip_clock];
  return thread_mip_clocks[local_thread_id].timer_pointer_->running(
      highs_timer_clock);
}

double HighsMipAnalysis::mipTimerRead(const HighsInt mip_clock
                                      // , const HighsInt thread_id
) const {
  if (!analyse_mip_time) return 0;
  HighsInt local_thread_id = 0;  // highs::parallel::thread_num();
  HighsInt highs_timer_clock =
      thread_mip_clocks[local_thread_id].clock_[mip_clock];
  return thread_mip_clocks[local_thread_id].timer_pointer_->read(
      highs_timer_clock);
}

HighsInt HighsMipAnalysis::mipTimerNumCall(const HighsInt mip_clock
                                           // , const HighsInt thread_id
) const {
  if (!analyse_mip_time) return 0;
  HighsInt local_thread_id = 0;  // highs::parallel::thread_num();
  HighsInt highs_timer_clock =
      thread_mip_clocks[local_thread_id].clock_[mip_clock];
  return thread_mip_clocks[local_thread_id].timer_pointer_->numCall(
      highs_timer_clock);
}

void HighsMipAnalysis::reportMipSolveLpClock(const bool header) {
  if (header) {
    printf(
        ",simplex time,IPM time,#simplex,#IPM,simplex/total time,IPM/total "
        "time,#No basis solve,simplex/#Basis solve,simplex/#No basis solve\n");
    return;
  }
  if (!analyse_mip_time) return;
  double total_time = mip_clocks.timer_pointer_->read(0);
  if (total_time < 0.01) return;
  HighsInt simplex_basis_solve_iclock =
      mip_clocks.clock_[kMipClockDuSimplexBasisSolveLp];
  HighsInt simplex_no_basis_solve_iclock =
      mip_clocks.clock_[kMipClockDuSimplexNoBasisSolveLp];
  HighsInt ipm_solve_iclock = mip_clocks.clock_[kMipClockIpxSolveLp];
  //  HighsInt num_no_basis_solve =
  //  mip_clocks.timer_pointer_->clock_num_call[no_basis_solve_iclock]; HighsInt
  //  num_basis_solve =
  //  mip_clocks.timer_pointer_->clock_num_call[basis_solve_iclock];
  HighsInt num_simplex_basis_solve =
      mip_clocks.timer_pointer_->clock_num_call[simplex_basis_solve_iclock];
  HighsInt num_simplex_no_basis_solve =
      mip_clocks.timer_pointer_->clock_num_call[simplex_no_basis_solve_iclock];
  HighsInt num_ipm_solve =
      mip_clocks.timer_pointer_->clock_num_call[ipm_solve_iclock];
  HighsInt num_simplex_solve =
      num_simplex_basis_solve + num_simplex_no_basis_solve;
  //  assert(num_no_basis_solve+num_basis_solve == num_simplex_solve);
  double simplex_basis_solve_time =
      mip_clocks.timer_pointer_->read(simplex_basis_solve_iclock);
  double simplex_no_basis_solve_time =
      mip_clocks.timer_pointer_->read(simplex_no_basis_solve_iclock);
  double simplex_solve_time =
      simplex_basis_solve_time + simplex_no_basis_solve_time;
  double ipm_solve_time = mip_clocks.timer_pointer_->read(ipm_solve_iclock);
  double frac_simplex_solve_time = simplex_solve_time / total_time;
  double frac_ipm_solve_time = ipm_solve_time / total_time;
  double average_simplex_basis_solve_time =
      num_simplex_basis_solve > 0
          ? simplex_basis_solve_time / int(num_simplex_basis_solve)
          : 0.0;
  double average_simplex_no_basis_solve_time =
      num_simplex_no_basis_solve > 0
          ? simplex_no_basis_solve_time / int(num_simplex_no_basis_solve)
          : 0.0;
  printf(",%11.2g,%11.2g,%d,%d,%11.2g,%11.2g,%d,%11.2g,%11.2g\n",
         simplex_solve_time, ipm_solve_time, int(num_simplex_solve),
         int(num_ipm_solve), frac_simplex_solve_time, frac_ipm_solve_time,
         int(num_simplex_no_basis_solve), average_simplex_basis_solve_time,
         average_simplex_no_basis_solve_time);
  printf(
      "LP solver analysis: %d LP with %d simplex (%11.2g CPU), %d IPM (%11.2g "
      "CPU) and %d solved without basis; average simplex solve time "
      "(basis/no_basis) = (%11.2g, %11.2g)\n",
      int(num_simplex_solve + num_ipm_solve), int(num_simplex_solve),
      simplex_solve_time, int(num_ipm_solve), ipm_solve_time,
      int(num_simplex_no_basis_solve), average_simplex_basis_solve_time,
      average_simplex_no_basis_solve_time);
};

void HighsMipAnalysis::reportMipTimer() {
  if (!analyse_mip_time) return;
  MipTimer mip_timer;
  mip_timer.reportMipCoreClock(mip_clocks);
  mip_timer.reportMipLevel1Clock(mip_clocks);
  mip_timer.reportMipEvaluateRootNodeClock(mip_clocks);
  //  mip_timer.reportAltEvaluateRootNodeClock(mip_clocks);
  //  mip_timer.reportMipPresolveClock(mip_clocks);
  //  mip_timer.reportMipRootSeparationClock(mip_clocks);
  //  mip_timer.reportMipSearchClock(mip_clocks);
  //  mip_timer.reportMipDiveClock(mip_clocks);
  //  mip_timer.reportMipNodeSearchClock(mip_clocks);
  //  mip_timer.reportMipDivePrimalHeuristicsClock(mip_clocks);
  mip_timer.reportMipSubMipSolveClock(mip_clocks);
  mip_timer.reportMipSeparationClock(mip_clocks);
  mip_timer.reportMipSolveLpClock(mip_clocks);
  //  mip_timer.csvMipClock(this->model_name, mip_clocks, true, false);
  //  reportMipSolveLpClock(true);
  //
  //  mip_timer.csvMipClock(this->model_name, mip_clocks, false, false);
  //  reportMipSolveLpClock(false);
  //
  //  mip_timer.csvEvaluateRootNodeClock(this->model_name, mip_clocks, true,
  //  true);
  //
  //  mip_timer.csvEvaluateRootNodeClock(this->model_name, mip_clocks, false,
  //  true);
  //
  //  analyseVectorValues(nullptr, "Node search time",
  //                      HighsInt(node_search_time.size()), node_search_time);
  //
  //  analyseVectorValues(nullptr, "Dive time", HighsInt(dive_time.size()),
  //                      dive_time);
  // mip_timer.reportFjClock(this->model_name, mip_clocks);
}

HighsInt HighsMipAnalysis::getSepaClockIndex(const std::string& name) const {
  HighsInt num_sepa_clock = this->sepa_name_clock.size();
  assert(num_sepa_clock > 0);
  for (HighsInt iSepaClock = 0; iSepaClock < num_sepa_clock; iSepaClock++) {
    if (this->sepa_name_clock[iSepaClock].first == name)
      return this->sepa_name_clock[iSepaClock].second;
  }
  return -1;
}
