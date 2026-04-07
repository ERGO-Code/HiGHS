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

struct HighsMipTimerClock {
  HighsTimer* timer_pointer_;
  std::vector<HighsInt> mip_clock_;
  std::vector<HighsInt> submip_clock_;
};

class HighsMipAnalysis {
 public:
  HighsMipAnalysis()
      : timer_(nullptr),
        sub_solver_call_time_(nullptr),
        analyse_mip_time(false) {}

  HighsTimer* timer_;
  HighsSubSolverCallTime* sub_solver_call_time_;
  void setupMipTime(const HighsOptions& options);
  void mipTimerStart(const HighsInt mip_clock = 0) const;
  void mipTimerStop(const HighsInt mip_clock = 0) const;
  bool mipTimerRunning(const HighsInt mip_clock = 0) const;
  double mipTimerRead(const HighsInt mip_clock = 0) const;
  HighsInt mipTimerNumCall(const HighsInt mip_clock = 0) const;
  void mipTimerAdd(const HighsInt mip_clock, const HighsInt num_call,
                   const double time) const;
  void mipTimerUpdate(const HighsSubSolverCallTime& sub_solver_call_time,
                      const bool valid_basis, const bool presolve,
                      const bool analytic_centre = false) const;
  void reportMipSolveLpClock(const bool header);
  void reportMipTimer();

  HighsInt getSepaClockIndex(const std::string& name) const;
  void addSubSolverCallTime(const HighsSubSolverCallTime& sub_solver_call_time,
                            const bool analytic_centre = false) const;
  void checkSubSolverCallTime(
      const HighsSubSolverCallTime& sub_solver_call_time);
  std::string model_name;
  HighsTimerClock mip_clocks;
  std::vector<HighsTimerClock> thread_mip_clocks;

  bool submip_;
  std::vector<HighsTimerClock> thread_mip_clocks_;
  std::vector<HighsTimerClock> thread_submip_clocks_;

  bool analyse_mip_time;
  std::vector<double> dive_time;
  std::vector<double> node_search_time;
  std::vector<std::pair<std::string, HighsInt>> sepa_name_clock;
};

#endif /* MIP_HIGHSMIPANALYSIS_H_ */
