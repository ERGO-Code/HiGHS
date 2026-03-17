/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsMipAnalysis.h
 * @brief Analyse MIP iterations, both for run-time control and data
 * gathering
 */
#ifndef MIP_HIGHSMIPANALYSIS_H_
#define MIP_HIGHSMIPANALYSIS_H_

#include "lp_data/HighsAnalysis.h"
#include "lp_data/HighsLp.h"
#include "util/HighsTimer.h"

class HighsMipAnalysis {
 public:
  HighsMipAnalysis()
      : timer_(nullptr),
        sub_solver_call_time_(nullptr),
        analyse_mip_time(false) {}

  HighsTimer* timer_;
  HighsSubSolverCallTime* sub_solver_call_time_;
  void setup(const HighsLp& lp, const HighsOptions& options);

  void setupMipTime(const HighsOptions& options);
  void mipTimerStart(const HighsInt mip_clock,
                     const HighsInt thread_id = 0) const;
  void mipTimerStop(const HighsInt mip_clock,
                    const HighsInt thread_id = 0) const;
  bool mipTimerRunning(const HighsInt mip_clock,
                       const HighsInt thread_id = 0) const;
  double mipTimerRead(const HighsInt mip_clock,
                      const HighsInt thread_id = 0) const;
  HighsInt mipTimerNumCall(const HighsInt mip_clock,
                           const HighsInt thread_id = 0) const;
  void mipTimerAdd(const HighsInt mip_clock, const HighsInt num_call,
                   const double time, const HighsInt thread_id = 0) const;
  void mipTimerUpdate(const HighsSubSolverCallTime& sub_solver_call_time,
                      const bool valid_basis, const bool presolve,
                      const bool analytic_centre,
                      const HighsInt thread_id = 0) const;
  void reportMipSolveLpClock(const bool header);
  void reportMipTimer();

  HighsTimerClock* getThreadMipTimerClockPointer();

  const std::vector<HighsTimerClock>& getThreadMipTimerClocks() {
    return thread_mip_clocks;
  }
  HighsTimerClock* getThreadMipTimerClockPtr(HighsInt i) {
    assert(i >= 0 && i < (HighsInt)thread_mip_clocks.size());
    return &thread_mip_clocks[i];
  }

  HighsInt getSepaClockIndex(const std::string& name) const;
  void addSubSolverCallTime(const HighsSubSolverCallTime& sub_solver_call_time,
                            const bool analytic_centre = false) const;
  void checkSubSolverCallTime(
      const HighsSubSolverCallTime& sub_solver_call_time);
  std::string model_name;
  HighsTimerClock mip_clocks;
  std::vector<HighsTimerClock> thread_mip_clocks;
  HighsTimerClock* pointer_serial_mip_clocks; // No longer used

  bool analyse_mip_time;
  std::vector<double> dive_time;
  std::vector<double> node_search_time;
  std::vector<std::pair<std::string, HighsInt>> sepa_name_clock;
};

#endif /* MIP_HIGHSMIPANALYSIS_H_ */
