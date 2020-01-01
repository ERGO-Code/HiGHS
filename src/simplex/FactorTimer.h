/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
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
  FactorInvert = 0,    //!< INVERT
  FactorInvertSimple,  //!< INVERT simple
  FactorInvertKernel,  //!< INVERT kernel
  FactorInvertDeficient, //!< INVERT deficient
  FactorInvertFinish,  //!< INVERT finish
  FactorFtran,         //!< FTRAN
  FactorFtranLower,    //!< FTRAN Lower part
  FactorFtranLowerAPF, //!< FTRAN Lower part APF
  FactorFtranLowerSps, //!< FTRAN Lower part sparse
  FactorFtranLowerHys, //!< FTRAN Lower part hyper-sparse
  FactorFtranUpper,    //!< FTRAN Upper part 
  FactorFtranUpperFT,  //!< FTRAN Upper part FT
  FactorFtranUpperMPF, //!< FTRAN Upper part MPF
  FactorFtranUpperSps, //!< FTRAN Upper part sparse
  FactorFtranUpperHys, //!< FTRAN Upper part hyper-sparse
  FactorFtranUpperPF,  //!< FTRAN Upper part PF
  FactorBtran,         //!< BTRAN
  FactorBtranLower,    //!< BTRAN Lower part
  FactorBtranLowerSps, //!< BTRAN Lower part sparse
  FactorBtranLowerHys, //!< BTRAN Lower part hyper-sparse
  FactorBtranLowerAPF, //!< BTRAN Lower part APF
  FactorBtranUpper,    //!< BTRAN Upper part
  FactorBtranUpperPF,  //!< BTRAN Upper part PF
  FactorBtranUpperSps, //!< BTRAN Upper part sparse
  FactorBtranUpperHys, //!< BTRAN Upper part hyper-sparse
  FactorBtranUpperFT,  //!< BTRAN Upper part FT
  FactorBtranUpperMPF, //!< BTRAN Upper part MPF
  FactorNumClock     //!< Number of factor clocks
};

class FactorTimer {
 public:
  void initialiseFactorClocks(HighsTimer& timer, std::vector<int>& clock) {
    clock.resize(FactorNumClock);
    clock[FactorInvert] = timer.clock_def("INVERT", "");
    clock[FactorInvertSimple] =  timer.clock_def("INVERT Simple", "INV");
    clock[FactorInvertKernel] =  timer.clock_def("INVERT Kernel", "IVK");
    clock[FactorInvertDeficient] = timer.clock_def("INVERT Deficient", "IVD");
    clock[FactorInvertFinish] =  timer.clock_def("INVERT Finish", "IVF");
    clock[FactorFtran] =         timer.clock_def("FTRAN", "FTR");
    clock[FactorFtranLower] =    timer.clock_def("FTRAN Lower", "FTL");
    clock[FactorFtranLowerAPF] = timer.clock_def("FTRAN Lower APF", "FLA");
    clock[FactorFtranLowerSps] = timer.clock_def("FTRAN Lower Sparse", "FLS");
    clock[FactorFtranLowerHys] = timer.clock_def("FTRAN Lower Hyper", "FLH");
    clock[FactorFtranUpper] =    timer.clock_def("FTRAN Upper", "FTU");
    clock[FactorFtranUpperFT] =  timer.clock_def("FTRAN Upper FT", "FUF");
    clock[FactorFtranUpperMPF] = timer.clock_def("FTRAN Upper MPF", "FUM");
    clock[FactorFtranUpperSps] = timer.clock_def("FTRAN Upper Sparse", "FUS");
    clock[FactorFtranUpperHys] = timer.clock_def("FTRAN Upper Hyper", "FUH");
    clock[FactorFtranUpperPF] =  timer.clock_def("FTRAN Upper PF", "FUP");
    clock[FactorBtran] =         timer.clock_def("BTRAN", "BTR");
    clock[FactorBtranLower] =    timer.clock_def("BTRAN Lower", "BTL");
    clock[FactorBtranLowerSps] = timer.clock_def("BTRAN Lower Sparse", "BLS");
    clock[FactorBtranLowerHys] = timer.clock_def("BTRAN Lower Hyper", "BLH");
    clock[FactorBtranLowerAPF] = timer.clock_def("BTRAN Lower APF", "BLA");
    clock[FactorBtranUpper] =    timer.clock_def("BTRAN Upper", "BTU");
    clock[FactorBtranUpperPF] =  timer.clock_def("BTRAN Upper PF", "BUP");
    clock[FactorBtranUpperSps] = timer.clock_def("BTRAN Upper Sparse", "BUS");
    clock[FactorBtranUpperHys] = timer.clock_def("BTRAN Upper Hyper", "BUH");
    clock[FactorBtranUpperFT] =  timer.clock_def("BTRAN Upper FT", "BUF");
    clock[FactorBtranUpperMPF] = timer.clock_def("BTRAN Upper MPS", "BUM");
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
      FactorInvert, FactorFtran, FactorBtran
	};
    reportFactorClockList("FactorLevel0", timer, clock, factor_clock_list);
  };
  
  void reportFactorLevel1Clock(HighsTimer& timer, std::vector<int>& clock) {
    std::vector<int> factor_clock_list{
      FactorInvertSimple, FactorInvertKernel,
	FactorInvertDeficient, FactorInvertFinish,
	FactorFtranLower, FactorFtranUpper,
	FactorBtranLower, FactorBtranUpper
	};
    reportFactorClockList("FactorLevel1", timer, clock, factor_clock_list);
  };
  
  void reportFactorLevel2Clock(HighsTimer& timer, std::vector<int>& clock) {
    std::vector<int> factor_clock_list{
      FactorInvertSimple, FactorInvertKernel,
	FactorInvertDeficient, FactorInvertFinish,
      FactorFtranLowerAPF, FactorFtranLowerSps, FactorFtranLowerHys,
	FactorFtranUpperFT, FactorFtranUpperMPF, FactorFtranUpperSps,
	FactorFtranUpperHys, FactorFtranUpperPF,
	FactorBtranLowerSps, FactorBtranLowerHys, FactorBtranLowerAPF,
	FactorBtranUpperPF, FactorBtranUpperSps, FactorBtranUpperHys,
	FactorBtranUpperFT, FactorBtranUpperMPF
	};
    reportFactorClockList("FactorLevel2", timer, clock, factor_clock_list);
  };
  
};
#endif /* SIMPLEX_FACTORTIMER_H_ */
