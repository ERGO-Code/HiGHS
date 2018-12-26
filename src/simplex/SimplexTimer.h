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
  Group1Clock = 0,            //!< Group for SIP
  IterateClock,           //!< Top level timing of HDual::solve_phase1() and HDual::solve_phase2()
  IterateRebuildClock,   //!< Second level timing of rebuild()
  IterateChuzrClock,     //!< Second level timing of CHUZR
  IterateChuzcClock,     //!< Second level timing of CHUZC
  IterateFtranClock,     //!< Second level timing of FTRAN
  IterateVerifyClock,    //!< Second level timing of numerical check
  IterateDualClock,      //!< Second level timing of dual update
  IteratePrimalClock,    //!< Second level timing of primal update
  IterateDevexIzClock,  //!< Second level timing of initialise Devex
  IteratePivotsClock,    //!< Second level timing of pivoting
  InvertClock,          //!< Invert in rebuild()
  PermWtClock,         //!< Permutation of SED weights each side of INVERT in rebuild()
  ComputeDualClock,    //!< Computation of dual values in rebuild()
  CorrectDualClock,    //!< Correction of dual values in rebuild()
  CollectPrIfsClock,  //!< Identification of primal infeasibilities in rebuild()
  ComputePrimalClock,  //!< Computation of primal values in rebuild()
  ComputeDuobjClock,   //!< Computation of dual objective value in rebuild()
  ReportInvertClock,   //!< Reporting of log line in rebuild()
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
  SimplexNumClock   //!< Number of simplex clocks
};

class SimplexTimer {
 public:
  void initialiseDualSimplexClocks(HighsModelObject & model_object) {
  HighsTimer & timer = model_object.timer_;
  HighsSimplexInfo & simplex = model_object.simplex_;
  simplex.clock_.resize(SimplexNumClock);
  simplex.clock_[Group1Clock] = timer.clockDef("GROUP1", "GP1");
  simplex.clock_[IterateClock] = timer.clockDef("ITERATE", "ITR");
  simplex.clock_[IterateRebuildClock] = timer.clockDef("REBUILD", "INV");
  simplex.clock_[IterateChuzrClock] = timer.clockDef("CHUZR", "CZR");
  simplex.clock_[IterateChuzcClock] = timer.clockDef("CHUZC", "CZC");
  simplex.clock_[IterateFtranClock] = timer.clockDef("FTRAN", "FTR");
  simplex.clock_[IterateVerifyClock] = timer.clockDef("VERIFY", "VRF");
  simplex.clock_[IterateDualClock] = timer.clockDef("DUAL", "UDU");
  simplex.clock_[IteratePrimalClock] = timer.clockDef("PRIMAL", "UPR");
  simplex.clock_[IterateDevexIzClock] = timer.clockDef("DEVEX_IZ", "DIZ");
  simplex.clock_[IteratePivotsClock] = timer.clockDef("PIVOTS", "PIV");
  simplex.clock_[InvertClock] = timer.clockDef("INVERT", "INV");
  simplex.clock_[PermWtClock] = timer.clockDef("PERM_WT", "PWT");
  simplex.clock_[ComputeDualClock] = timer.clockDef("COMPUTE_DUAL", "CPD");
  simplex.clock_[CorrectDualClock] = timer.clockDef("CORRECT_DUAL", "CRD");
  simplex.clock_[ComputePrimalClock] = timer.clockDef("COMPUTE_PRIMAL", "CPP");
  simplex.clock_[CollectPrIfsClock] = timer.clockDef("COLLECT_PR_IFS", "IFS");
  simplex.clock_[ComputeDuobjClock] = timer.clockDef("COMPUTE_DUOBJ", "DOB");
  simplex.clock_[ReportInvertClock] = timer.clockDef("REPORT_INVERT", "RPI");
  simplex.clock_[Chuzr1Clock] = timer.clockDef("CHUZR1", "CR1");
  simplex.clock_[Chuzc0Clock] = timer.clockDef("CHUZC0", "CC0");
  simplex.clock_[Chuzc1Clock] = timer.clockDef("CHUZC1", "CC1");
  simplex.clock_[Chuzc2Clock] = timer.clockDef("CHUZC2", "CC2");
  simplex.clock_[Chuzc3Clock] = timer.clockDef("CHUZC3", "CC3");
  simplex.clock_[Chuzc4Clock] = timer.clockDef("CHUZC4", "CC4");
  simplex.clock_[DevexWtClock] = timer.clockDef("DEVEX_WT", "DWT");
  simplex.clock_[FtranClock] = timer.clockDef("FTRAN", "COL");
  simplex.clock_[BtranClock] = timer.clockDef("BTRAN", "REP");
  simplex.clock_[PriceClock] = timer.clockDef("PRICE", "RAP");
  simplex.clock_[FtranDseClock] = timer.clockDef("FTRAN_DSE", "DSE");
  simplex.clock_[FtranMixClock] = timer.clockDef("FTRAN_MIX", "MIX");
  simplex.clock_[FtranBfrtClock] = timer.clockDef("FTRAN_BFRT", "BFR");
  simplex.clock_[UpdateDualClock] = timer.clockDef("UPDATE_DUAL", "UPD");
  simplex.clock_[UpdatePrimalClock] = timer.clockDef("UPDATE_PRIMAL", "UPP");
  simplex.clock_[DevexIzClock] = timer.clockDef("DEVEX_IZ", "DIZ");
  simplex.clock_[UpdateWeightClock] = timer.clockDef("UPDATE_WEIGHT", "UPW");
  simplex.clock_[UpdatePivotsClock] = timer.clockDef("UPDATE_PIVOTS", "UPP");
  simplex.clock_[UpdateFactorClock] = timer.clockDef("UPDATE_FACTOR", "UPF");
  simplex.clock_[UpdateMatrixClock] = timer.clockDef("UPDATE_MATRIX", "UPM");
  simplex.clock_[UpdateRowEpClock] = timer.clockDef("UPDATE_ROW_EP", "UPR");
};

void reportSimplexClockList(std::vector<int> simplexClockList, HighsModelObject & model_object) {
  HighsTimer & timer = model_object.timer_;
  HighsSimplexInfo & simplex = model_object.simplex_;
  int simplexClockListSize = simplexClockList.size();
  std::vector<int> clockList;
  clockList.resize(simplexClockListSize);
  for (int en=0; en<simplexClockListSize; en++) {
    clockList[en] = simplex.clock_[simplexClockList[en]];
  }
  timer.report_tl(clockList, 0.0);
};

void reportDualSimplexIterateClock(HighsModelObject & model_object) {
  std::vector<int> simplexClockList{
    IterateClock
      };
  reportSimplexClockList(simplexClockList, model_object);
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
  reportSimplexClockList(simplexClockList, model_object);
};

void reportDualSimplexOuterClock(HighsModelObject & model_object) {
  std::vector<int> simplexClockList{
    IterateRebuildClock, IterateChuzrClock, IterateChuzcClock,
      IterateFtranClock, IterateVerifyClock, IterateDualClock,
      IteratePrimalClock, IterateDevexIzClock, IteratePivotsClock
      };
  reportSimplexClockList(simplexClockList, model_object);
};
};
#endif /* SIMPLEX_SIMPLEXTIMER_H_ */
