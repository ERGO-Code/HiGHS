/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/SimplexTimer.h
 * @brief Indices of simplex iClocks
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_SIMPLEXTIMER_H_
#define SIMPLEX_SIMPLEXTIMER_H_

#include "HighsModelObject.h"

// Clocks for profiling the dual simplex solver
enum iClockSimplex {
  SimplexTotalClock = 0, //!< Total time for simplex
  SimplexIzDseWtClock, //!< Total time to initialise DSE weights
  SimplexDualPhase1Clock, //!< Total time for dual simplex phase 1
  SimplexDualPhase2Clock, //!< Total time for dual simplex phase 2
  SimplexPrimalPhase2Clock, //!< Total time for primal simplex phase 2
  Group1Clock,            //!< Group for SIP

  IterateClock,           //!< Top level timing of HDual::solve_phase1() and HDual::solve_phase2()
  IterateDualRebuildClock,   //!< Second level timing of dual rebuild()
  IterateChuzrClock,     //!< Second level timing of CHUZR
  IterateChuzcClock,     //!< Second level timing of CHUZC
  IterateFtranClock,     //!< Second level timing of FTRAN
  IterateVerifyClock,    //!< Second level timing of numerical check
  IterateDualClock,      //!< Second level timing of dual update
  IteratePrimalClock,    //!< Second level timing of primal update
  IterateDevexIzClock,  //!< Second level timing of initialise Devex
  IteratePivotsClock,    //!< Second level timing of pivoting

  InvertClock,          //!< Invert in dual rebuild()
  PermWtClock,         //!< Permutation of SED weights each side of INVERT in dual rebuild()
  ComputeDualClock,    //!< Computation of dual values in dual rebuild()
  CorrectDualClock,    //!< Correction of dual values in dual rebuild()
  CollectPrIfsClock,  //!< Identification of primal infeasibilities in dual rebuild()
  ComputePrimalClock,  //!< Computation of primal values in dual rebuild()
  ComputeDuobjClock,   //!< Computation of dual objective value in dual rebuild()
  ReportInvertClock,   //!< Reporting of log line in dual rebuild()
  Chuzr1Clock,          //!< CHUZR
  Chuzc0Clock,          //!< CHUZC - stage 0
  Chuzc1Clock,          //!< CHUZC - stage 1
  Chuzc2Clock,          //!< CHUZC - stage 2
  Chuzc3Clock,          //!< CHUZC - stage 3
  Chuzc4Clock,          //!< CHUZC - stage 4
  DevexWtClock,        //!< Calculation of Devex weight of entering variable
  FtranClock,           //!< FTRAN - pivotal column
  BtranClock,           //!< BTRAN
  PriceClock,           //!< PRICE
  FtranDseClock,       //!< FTRAN for DSE weights
  FtranMixClock,       //!< FTRAN for PAMI
  FtranBfrtClock,      //!< FTRAN for BFRT
  UpdateDualClock,     //!< Update of dual values
  UpdatePrimalClock,   //!< Update of primal values
  DevexIzClock,        //!< Initialisation of new Devex framework
  UpdateWeightClock,   //!< Update of DSE or Devex weights
  UpdatePivotsClock,   //!< Update indices of basic and nonbasic after basis change
  UpdateFactorClock,   //!< Update the representation of \f$B^{-1}\f$
  UpdateMatrixClock,   //!< Update the row-wise copy of the constraint matrix for nonbasic columns
  UpdateRowEpClock,   //!< Update the tableau rows in PAMI

  IteratePrimalRebuildClock,   //!< Second level timing of primal rebuild()

  SimplexNumClock   //!< Number of simplex clocks
};

class SimplexTimer {
 public:
  void initialiseDualSimplexClocks(HighsModelObject & model_object) {
    HighsTimer & timer = model_object.timer_;
    HighsSimplexInfo & simplex_info = model_object.simplex_info_;
    simplex_info.clock_.resize(SimplexNumClock);
    simplex_info.clock_[SimplexTotalClock] = timer.clockDef("Simplex total", "STT");
    simplex_info.clock_[SimplexIzDseWtClock] = timer.clockDef("Iz DSE Wt", "IWT");
    simplex_info.clock_[SimplexDualPhase1Clock] = timer.clockDef("Dual Phase 1", "DP1");
    simplex_info.clock_[SimplexDualPhase2Clock] = timer.clockDef("Dual Phase 2", "DP2");
    simplex_info.clock_[SimplexPrimalPhase2Clock] = timer.clockDef("Primal Phase 2", "PP2");
    simplex_info.clock_[Group1Clock] = timer.clockDef("GROUP1", "GP1");
    simplex_info.clock_[IterateClock] = timer.clockDef("ITERATE", "ITR");
    simplex_info.clock_[IterateDualRebuildClock] = timer.clockDef("DUAL REBUILD", "DRB");
    simplex_info.clock_[IterateChuzrClock] = timer.clockDef("CHUZR", "CZR");
    simplex_info.clock_[IterateChuzcClock] = timer.clockDef("CHUZC", "CZC");
    simplex_info.clock_[IterateFtranClock] = timer.clockDef("FTRAN", "FTR");
    simplex_info.clock_[IterateVerifyClock] = timer.clockDef("VERIFY", "VRF");
    simplex_info.clock_[IterateDualClock] = timer.clockDef("DUAL", "UDU");
    simplex_info.clock_[IteratePrimalClock] = timer.clockDef("PRIMAL", "UPR");
    simplex_info.clock_[IterateDevexIzClock] = timer.clockDef("DEVEX_IZ", "DIZ");
    simplex_info.clock_[IteratePivotsClock] = timer.clockDef("PIVOTS", "PIV");
    simplex_info.clock_[InvertClock] = timer.clockDef("INVERT", "INV");
    simplex_info.clock_[PermWtClock] = timer.clockDef("PERM_WT", "PWT");
    simplex_info.clock_[ComputeDualClock] = timer.clockDef("COMPUTE_DUAL", "CPD");
    simplex_info.clock_[CorrectDualClock] = timer.clockDef("CORRECT_DUAL", "CRD");
    simplex_info.clock_[ComputePrimalClock] = timer.clockDef("COMPUTE_PRIMAL", "CPP");
    simplex_info.clock_[CollectPrIfsClock] = timer.clockDef("COLLECT_PR_IFS", "IFS");
    simplex_info.clock_[ComputeDuobjClock] = timer.clockDef("COMPUTE_DUOBJ", "DOB");
    simplex_info.clock_[ReportInvertClock] = timer.clockDef("REPORT_INVERT", "RPI");
    simplex_info.clock_[Chuzr1Clock] = timer.clockDef("CHUZR1", "CR1");
    simplex_info.clock_[Chuzc0Clock] = timer.clockDef("CHUZC0", "CC0");
    simplex_info.clock_[Chuzc1Clock] = timer.clockDef("CHUZC1", "CC1");
    simplex_info.clock_[Chuzc2Clock] = timer.clockDef("CHUZC2", "CC2");
    simplex_info.clock_[Chuzc3Clock] = timer.clockDef("CHUZC3", "CC3");
    simplex_info.clock_[Chuzc4Clock] = timer.clockDef("CHUZC4", "CC4");
    simplex_info.clock_[DevexWtClock] = timer.clockDef("DEVEX_WT", "DWT");
    simplex_info.clock_[FtranClock] = timer.clockDef("FTRAN", "COL");
    simplex_info.clock_[BtranClock] = timer.clockDef("BTRAN", "REP");
    simplex_info.clock_[PriceClock] = timer.clockDef("PRICE", "RAP");
    simplex_info.clock_[FtranDseClock] = timer.clockDef("FTRAN_DSE", "DSE");
    simplex_info.clock_[FtranMixClock] = timer.clockDef("FTRAN_MIX", "MIX");
    simplex_info.clock_[FtranBfrtClock] = timer.clockDef("FTRAN_BFRT", "BFR");
    simplex_info.clock_[UpdateDualClock] = timer.clockDef("UPDATE_DUAL", "UPD");
    simplex_info.clock_[UpdatePrimalClock] = timer.clockDef("UPDATE_PRIMAL", "UPP");
    simplex_info.clock_[DevexIzClock] = timer.clockDef("DEVEX_IZ", "DIZ");
    simplex_info.clock_[UpdateWeightClock] = timer.clockDef("UPDATE_WEIGHT", "UPW");
    simplex_info.clock_[UpdatePivotsClock] = timer.clockDef("UPDATE_PIVOTS", "UPP");
    simplex_info.clock_[UpdateFactorClock] = timer.clockDef("UPDATE_FACTOR", "UPF");
    simplex_info.clock_[UpdateMatrixClock] = timer.clockDef("UPDATE_MATRIX", "UPM");
    simplex_info.clock_[UpdateRowEpClock] = timer.clockDef("UPDATE_ROW_EP", "UPR");
    simplex_info.clock_[IteratePrimalRebuildClock] = timer.clockDef("PRIMAL REBUILD", "PRB");
  };
  
  void reportSimplexClockList(const char *grepStamp, std::vector<int> simplexClockList, HighsModelObject & model_object) {
    HighsTimer & timer = model_object.timer_;
    HighsSimplexInfo & simplex_info = model_object.simplex_info_;
    int simplexClockListSize = simplexClockList.size();
    std::vector<int> clockList;
    clockList.resize(simplexClockListSize);
    for (int en=0; en<simplexClockListSize; en++) {
      clockList[en] = simplex_info.clock_[simplexClockList[en]];
    }
    timer.report_tl(grepStamp, clockList, 0.0);
  };
  
  void reportSimplexTotalClock(HighsModelObject & model_object) {
    std::vector<int> simplexClockList{
      SimplexTotalClock
	};
    reportSimplexClockList("SimplexTotal", simplexClockList, model_object);
  };
  
  void reportSimplexPhasesClock(HighsModelObject & model_object) {
    std::vector<int> simplexClockList{
      SimplexIzDseWtClock,
	SimplexDualPhase1Clock, SimplexDualPhase2Clock, SimplexPrimalPhase2Clock
	};
    reportSimplexClockList("SimplexPhases", simplexClockList, model_object);
  };
  
  void reportDualSimplexIterateClock(HighsModelObject & model_object) {
    std::vector<int> simplexClockList{
      IterateClock
	};
    reportSimplexClockList("SimplexIterate", simplexClockList, model_object);
  };
  
  void reportDualSimplexOuterClock(HighsModelObject & model_object) {
    std::vector<int> simplexClockList{
      IterateDualRebuildClock, IterateChuzrClock, IterateChuzcClock,
	IterateFtranClock, IterateVerifyClock, IterateDualClock,
	IteratePrimalClock, IterateDevexIzClock, IteratePivotsClock
	};
    reportSimplexClockList("SimplexOuter", simplexClockList, model_object);
  };

  void reportDualSimplexInnerClock(HighsModelObject & model_object) {
    std::vector<int> simplexClockList{
      InvertClock, PermWtClock, ComputeDualClock, 
	CorrectDualClock, ComputePrimalClock, CollectPrIfsClock, 
	ComputeDuobjClock, ReportInvertClock, Chuzr1Clock, 
	BtranClock, PriceClock, Chuzc0Clock, 
	Chuzc1Clock, Chuzc2Clock, Chuzc3Clock, 
	Chuzc4Clock, DevexWtClock, FtranClock, 
	FtranBfrtClock, FtranDseClock, UpdateDualClock, 
	UpdatePrimalClock, UpdateWeightClock, DevexIzClock, 
	UpdatePivotsClock, UpdateFactorClock, UpdateMatrixClock
	};
    reportSimplexClockList("SimplexInner", simplexClockList, model_object);
  };
  
};
#endif /* SIMPLEX_SIMPLEXTIMER_H_ */
