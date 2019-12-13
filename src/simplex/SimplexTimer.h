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

#include "lp_data/HighsModelObject.h"

// Clocks for profiling the dual simplex solver
enum iClockSimplex {
  SimplexTotalClock = 0,     //!< Total time for simplex
  SimplexIzDseWtClock,       //!< Total time to initialise DSE weights
  SimplexDualPhase1Clock,    //!< Total time for dual simplex phase 1
  SimplexDualPhase2Clock,    //!< Total time for dual simplex phase 2
  SimplexPrimalPhase1Clock,  //!< Total time for primal simplex phase 1
  SimplexPrimalPhase2Clock,  //!< Total time for primal simplex phase 2
  Group1Clock,               //!< Group for SIP

  IterateClock,               //!< Top level timing of HDual::solve_phase1() and
                              //!< HDual::solve_phase2()
  IterateDualRebuildClock,    //!< Second level timing of dual rebuild()
  IteratePrimalRebuildClock,  //!< Second level timing of primal rebuild()
  IterateChuzrClock,          //!< Second level timing of CHUZR
  IterateChuzcClock,          //!< Second level timing of CHUZC
  IterateFtranClock,          //!< Second level timing of FTRAN
  IterateVerifyClock,         //!< Second level timing of numerical check
  IterateDualClock,           //!< Second level timing of dual update
  IteratePrimalClock,         //!< Second level timing of primal update
  IterateDevexIzClock,        //!< Second level timing of initialise Devex
  IteratePivotsClock,         //!< Second level timing of pivoting

  ScaleClock,           //!< Scale
  CrashClock,           //!< Crash
  BasisConditionClock,  //!< Basis condition estimation
  DseIzClock,           //!< DSE weight initialisation
  InvertClock,          //!< Invert in dual rebuild()
  PermWtClock,       //!< Permutation of SED weights each side of INVERT in dual
                     //!< rebuild()
  ComputeDualClock,  //!< Computation of dual values in dual rebuild()
  CorrectDualClock,  //!< Correction of dual values in dual rebuild()
  CollectPrIfsClock,   //!< Identification of primal infeasibilities in dual
                       //!< rebuild()
  ComputePrIfsClock,   //!< Computation of num/max/sum of primal infeasibilities
  ComputeDuIfsClock,   //!< Computation of num/max/sum of dual infeasibilities
  ComputePrimalClock,  //!< Computation of primal values in dual rebuild()
  ComputeDuObjClock,  //!< Computation of dual objective value in dual rebuild()
  ComputePrObjClock,  //!< Computation of primalal objective value in primal
                      //!< rebuild()
  ReportRebuildClock,  //!< Reporting of log line in dual rebuild()
  ChuzrDualClock,      //!< CHUZR - Dual
  Chuzr1Clock,         //!< CHUZR - Primal stage 1
  Chuzr2Clock,         //!< CHUZR - Primal stage 2
  ChuzcPrimalClock,    //!< CHUZC - Primal
  Chuzc0Clock,         //!< CHUZC - Dual stage 0
  PriceChuzc1Clock,    //!< PRICE + CHUZC - Dual stage 1: parallel
  Chuzc1Clock,         //!< CHUZC - Dual stage 1
  Chuzc2Clock,         //!< CHUZC - Dual stage 2
  Chuzc3Clock,         //!< CHUZC - Dual stage 3
  Chuzc4Clock,         //!< CHUZC - Dual stage 4
  DevexWtClock,        //!< Calculation of Devex weight of entering variable
  FtranClock,          //!< FTRAN - pivotal column
  BtranClock,          //!< BTRAN
  PriceClock,          //!< PRICE
  FtranDseClock,       //!< FTRAN for DSE weights
  FtranMixParClock,    //!< FTRAN for PAMI - parallel
  FtranMixFinalClock,  //!< FTRAN for PAMI - final
  FtranBfrtClock,      //!< FTRAN for BFRT
  UpdateRowClock,     //!< Update of dual values
  UpdateDualClock,     //!< Update of dual values
  UpdatePrimalClock,   //!< Update of primal values
  DevexIzClock,        //!< Initialisation of new Devex framework
  UpdateWeightClock,   //!< Update of DSE or Devex weights
  UpdatePivotsClock,   //!< Update indices of basic and nonbasic after basis
                       //!< change
  UpdateFactorClock,   //!< Update the representation of \f$B^{-1}\f$
  UpdateMatrixClock,  //!< Update the row-wise copy of the constraint matrix for
                      //!< nonbasic columns
  UpdateRowEpClock,   //!< Update the tableau rows in PAMI

  SimplexNumClock  //!< Number of simplex clocks
};

class SimplexTimer {
 public:
  void initialiseSimplexClocks(HighsModelObject& model_object) {
    HighsTimer& timer = model_object.timer_;
    HighsSimplexInfo& simplex_info = model_object.simplex_info_;
    simplex_info.clock_.resize(SimplexNumClock);
    simplex_info.clock_[SimplexTotalClock] =
        timer.clock_def("Simplex total", "STT");
    simplex_info.clock_[SimplexIzDseWtClock] =
        timer.clock_def("Iz DSE Wt", "IWT");
    simplex_info.clock_[SimplexDualPhase1Clock] =
        timer.clock_def("Dual Phase 1", "DP1");
    simplex_info.clock_[SimplexDualPhase2Clock] =
        timer.clock_def("Dual Phase 2", "DP2");
    simplex_info.clock_[SimplexPrimalPhase1Clock] =
        timer.clock_def("Primal Phase 1", "PP1");
    simplex_info.clock_[SimplexPrimalPhase2Clock] =
        timer.clock_def("Primal Phase 2", "PP2");
    simplex_info.clock_[Group1Clock] = timer.clock_def("GROUP1", "GP1");
    simplex_info.clock_[IterateClock] = timer.clock_def("ITERATE", "ITR");
    simplex_info.clock_[IterateDualRebuildClock] =
        timer.clock_def("DUAL REBUILD", "DRB");
    simplex_info.clock_[IteratePrimalRebuildClock] =
        timer.clock_def("PRIMAL REBUILD", "PRB");
    simplex_info.clock_[IterateChuzrClock] = timer.clock_def("CHUZR", "CZR");
    simplex_info.clock_[IterateChuzcClock] = timer.clock_def("CHUZC", "CZC");
    simplex_info.clock_[IterateFtranClock] = timer.clock_def("FTRAN", "FTR");
    simplex_info.clock_[IterateVerifyClock] = timer.clock_def("VERIFY", "VRF");
    simplex_info.clock_[IterateDualClock] = timer.clock_def("DUAL", "UDU");
    simplex_info.clock_[IteratePrimalClock] = timer.clock_def("PRIMAL", "UPR");
    simplex_info.clock_[IterateDevexIzClock] =
        timer.clock_def("DEVEX_IZ", "DVI");
    simplex_info.clock_[IteratePivotsClock] = timer.clock_def("PIVOTS", "PIV");
    simplex_info.clock_[ScaleClock] = timer.clock_def("SCALE", "SCL");
    simplex_info.clock_[CrashClock] = timer.clock_def("CRASH", "CSH");
    simplex_info.clock_[BasisConditionClock] =
        timer.clock_def("BASIS_CONDITION", "CON");
    simplex_info.clock_[DseIzClock] = timer.clock_def("DSE_IZ", "DEI");
    simplex_info.clock_[InvertClock] = timer.clock_def("INVERT", "INV");
    simplex_info.clock_[PermWtClock] = timer.clock_def("PERM_WT", "PWT");
    simplex_info.clock_[ComputeDualClock] =
        timer.clock_def("COMPUTE_DUAL", "CPD");
    simplex_info.clock_[CorrectDualClock] =
        timer.clock_def("CORRECT_DUAL", "CRD");
    simplex_info.clock_[ComputePrimalClock] =
        timer.clock_def("COMPUTE_PRIMAL", "CPP");
    simplex_info.clock_[CollectPrIfsClock] =
        timer.clock_def("COLLECT_PR_IFS", "IFS");
    simplex_info.clock_[ComputePrIfsClock] =
        timer.clock_def("COMPUTE_PR_IFS", "PIF");
    simplex_info.clock_[ComputeDuIfsClock] =
        timer.clock_def("COMPUTE_DU_IFS", "DIF");
    simplex_info.clock_[ComputeDuObjClock] =
        timer.clock_def("COMPUTE_DUOBJ", "DOB");
    simplex_info.clock_[ComputePrObjClock] =
        timer.clock_def("COMPUTE_PROBJ", "POB");
    simplex_info.clock_[ReportRebuildClock] =
        timer.clock_def("REPORT_REBUILD", "RPR");
    simplex_info.clock_[ChuzrDualClock] = timer.clock_def("CHUZR_DUAL", "CRD");
    simplex_info.clock_[Chuzr1Clock] = timer.clock_def("CHUZR1", "CR1");
    simplex_info.clock_[Chuzr2Clock] = timer.clock_def("CHUZR2", "CR2");
    simplex_info.clock_[ChuzcPrimalClock] =
        timer.clock_def("CHUZC_PRIMAL", "CCP");
    simplex_info.clock_[Chuzc0Clock] = timer.clock_def("CHUZC0", "CC0");
    simplex_info.clock_[PriceChuzc1Clock] = timer.clock_def("PRICE_CHUZC1", "PC1");
    simplex_info.clock_[Chuzc1Clock] = timer.clock_def("CHUZC1", "CC1");
    simplex_info.clock_[Chuzc2Clock] = timer.clock_def("CHUZC2", "CC2");
    simplex_info.clock_[Chuzc3Clock] = timer.clock_def("CHUZC3", "CC3");
    simplex_info.clock_[Chuzc4Clock] = timer.clock_def("CHUZC4", "CC4");
    simplex_info.clock_[DevexWtClock] = timer.clock_def("DEVEX_WT", "DWT");
    simplex_info.clock_[FtranClock] = timer.clock_def("FTRAN", "COL");
    simplex_info.clock_[BtranClock] = timer.clock_def("BTRAN", "REP");
    simplex_info.clock_[PriceClock] = timer.clock_def("PRICE", "RAP");
    simplex_info.clock_[FtranDseClock] = timer.clock_def("FTRAN_DSE", "DSE");
    simplex_info.clock_[FtranMixParClock] = timer.clock_def("FTRAN_MIX_PAR", "FMP");
    simplex_info.clock_[FtranMixFinalClock] = timer.clock_def("FTRAN_MIX_FINAL", "FMF");
    simplex_info.clock_[FtranBfrtClock] = timer.clock_def("FTRAN_BFRT", "BFR");
    simplex_info.clock_[UpdateRowClock] =
        timer.clock_def("UPDATE_ROW", "UPR");
    simplex_info.clock_[UpdateDualClock] =
        timer.clock_def("UPDATE_DUAL", "UPD");
    simplex_info.clock_[UpdatePrimalClock] =
        timer.clock_def("UPDATE_PRIMAL", "UPP");
    simplex_info.clock_[DevexIzClock] = timer.clock_def("DEVEX_IZ", "DIZ");
    simplex_info.clock_[UpdateWeightClock] =
        timer.clock_def("UPDATE_WEIGHT", "UPW");
    simplex_info.clock_[UpdatePivotsClock] =
        timer.clock_def("UPDATE_PIVOTS", "UPP");
    simplex_info.clock_[UpdateFactorClock] =
        timer.clock_def("UPDATE_FACTOR", "UPF");
    simplex_info.clock_[UpdateMatrixClock] =
        timer.clock_def("UPDATE_MATRIX", "UPM");
    simplex_info.clock_[UpdateRowEpClock] =
        timer.clock_def("UPDATE_ROW_EP", "UPR");
  };

  void report_simplex_clock_list(const char* grepStamp,
                                 std::vector<int> simplex_clock_list,
                                 HighsModelObject& model_object) {
    HighsTimer& timer = model_object.timer_;
    HighsSimplexInfo& simplex_info = model_object.simplex_info_;
    int simplex_clock_list_size = simplex_clock_list.size();
    std::vector<int> clockList;
    clockList.resize(simplex_clock_list_size);
    for (int en = 0; en < simplex_clock_list_size; en++) {
      clockList[en] = simplex_info.clock_[simplex_clock_list[en]];
    }
    timer.report_tl(grepStamp, clockList, 1e-8);
  };

  void reportSimplexTotalClock(HighsModelObject& model_object) {
    std::vector<int> simplex_clock_list{SimplexTotalClock};
    report_simplex_clock_list("SimplexTotal", simplex_clock_list, model_object);
  };

  void reportSimplexPhasesClock(HighsModelObject& model_object) {
    std::vector<int> simplex_clock_list{
        SimplexIzDseWtClock, SimplexDualPhase1Clock, SimplexDualPhase2Clock,
        SimplexPrimalPhase2Clock};
    report_simplex_clock_list("SimplexPhases", simplex_clock_list,
                              model_object);
  };

  void reportDualSimplexIterateClock(HighsModelObject& model_object) {
    std::vector<int> simplex_clock_list{IterateClock};
    report_simplex_clock_list("SimplexIterate", simplex_clock_list,
                              model_object);
  };

  void reportDualSimplexOuterClock(HighsModelObject& model_object) {
    std::vector<int> simplex_clock_list{
        IterateDualRebuildClock, IterateChuzrClock,   IterateChuzcClock,
        IterateFtranClock,       IterateVerifyClock,  IterateDualClock,
        IteratePrimalClock,      IterateDevexIzClock, IteratePivotsClock};
    report_simplex_clock_list("SimplexOuter", simplex_clock_list, model_object);
  };

  void reportSimplexInnerClock(HighsModelObject& model_object) {
    std::vector<int> simplex_clock_list{
        ScaleClock,
        CrashClock,        BasisConditionClock, DseIzClock,
        InvertClock,       PermWtClock,         ComputeDualClock,
        CorrectDualClock,  ComputePrimalClock,  CollectPrIfsClock,
        ComputePrIfsClock, ComputeDuIfsClock,   ComputeDuObjClock,
        ComputePrObjClock, ReportRebuildClock,  ChuzrDualClock,
        Chuzr1Clock,       Chuzr2Clock,         BtranClock,
        PriceClock,        ChuzcPrimalClock,    Chuzc0Clock,
        Chuzc1Clock,       Chuzc2Clock,         Chuzc3Clock,
        Chuzc4Clock,       DevexWtClock,        FtranClock,
        FtranBfrtClock,    FtranDseClock,       UpdateDualClock,
        UpdatePrimalClock, UpdateWeightClock,   DevexIzClock,
        UpdatePivotsClock, UpdateFactorClock,   UpdateMatrixClock};
    report_simplex_clock_list("SimplexInner", simplex_clock_list, model_object);
  };

  void reportSimplexMultiInnerClock(HighsModelObject& model_object) {
    std::vector<int> simplex_clock_list{
        ScaleClock,
        CrashClock,        BasisConditionClock, DseIzClock,
        InvertClock,       PermWtClock,         ComputeDualClock,
        CorrectDualClock,  ComputePrimalClock,  CollectPrIfsClock,
        ComputePrIfsClock, ComputeDuIfsClock,   ComputeDuObjClock,
        ComputePrObjClock, ReportRebuildClock,  ChuzrDualClock,
        Chuzr1Clock,       Chuzr2Clock,         BtranClock,
        PriceClock,        ChuzcPrimalClock,    Chuzc0Clock,
        PriceChuzc1Clock,
        Chuzc1Clock,       Chuzc2Clock,         Chuzc3Clock,
        Chuzc4Clock,       DevexWtClock,        FtranClock,
        FtranBfrtClock,    FtranDseClock,       FtranMixParClock,
	FtranMixFinalClock,UpdateRowClock,      UpdateDualClock,
        UpdatePrimalClock, UpdateWeightClock,   DevexIzClock,
        UpdatePivotsClock, UpdateFactorClock,   UpdateMatrixClock};
    report_simplex_clock_list("SimplexMultiInner", simplex_clock_list, model_object);
  };
};
#endif /* SIMPLEX_SIMPLEXTIMER_H_ */
