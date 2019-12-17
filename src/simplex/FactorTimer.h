/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/FactorTimer.h
 * @brief Indices of factor iClocks
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_FACTORTIMER_H_
#define SIMPLEX_FACTORTIMER_H_

#include "HConfig.h"
#include "util/HighsTimer.h"

// Clocks for profiling the dual simplex solver
enum iClockFactor {
  FactorFtran = 0,   //!< FTRAN
  FactorFtranLower,  //!< FTRAN Lower part
  FactorFtranUpper,  //!< FTRAN Upper part
  FactorNumClock     //!< Number of factor clocks
};

class FactorTimer {
 public:
  void initialiseFactorClocks(HighsTimer& timer, std::vector<int>& clock) {
    clock.resize(FactorNumClock);
    clock[FactorFtran] =
      timer.clock_def("FTRAN", "FTR");
    clock[FactorFtranLower] =
      timer.clock_def("FTRAN Lower", "FTL");
    clock[FactorFtranUpper] =
      timer.clock_def("FTRAN Upper", "FTU");
  };

  void reportFactorClockList(const char* grepStamp,
			     HighsTimer& timer, std::vector<int>& clock,
			     std::vector<int> factor_clock_list) {
    int factor_clock_list_size = factor_clock_list.size();
    std::vector<int> clockList;
    clockList.resize(factor_clock_list_size);
    for (int en = 0; en < factor_clock_list_size; en++) {
      clockList[en] = clock[factor_clock_list[en]];
    }
    timer.report_tl(grepStamp, clockList, 1e-8);
  };
  
  void reportFactorLevel0Clock(HighsTimer& timer, std::vector<int>& clock) {
    std::vector<int> factor_clock_list{
      FactorFtran
	};
    reportFactorClockList("FactorLevel0", timer, clock, factor_clock_list);
  };
  
  void reportFactorLevel1Clock(HighsTimer& timer, std::vector<int>& clock) {
    std::vector<int> factor_clock_list{
      FactorFtranLower, FactorFtranUpper
	};
    reportFactorClockList("FactorLevel1", timer, clock, factor_clock_list);
  };
  
};
#endif /* SIMPLEX_FACTORTIMER_H_ */
