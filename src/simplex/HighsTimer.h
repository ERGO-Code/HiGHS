/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HighsTimer.h
 * @brief Profiling facility for computational components in HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HIGHSTIMER_H_
#define SIMPLEX_HIGHSTIMER_H_

#include <sys/time.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>
/**
 * @brief Class for profiling facility for computational components in HiGHS
 */
class HighsTimer {
 public:

  HighsTimer() {
    startTime = getWallTime();
    startTick = getWallTick();
    numClock = 0;
    initialiseDualSimplexClocks();
  }

  /**
   * @brief Define a clock
   */
  int clockDef(
	       const char *name,   //!< Full-length name (<=16 characters) for the clock
	       const char *ch3name //!< 3-character name for the clock
	       ) {
    int iClock = numClock;
    clockNumCall.push_back(0);
    clockStart.push_back(initialClockStart);
    clockTicks.push_back(0);
    clockNames.push_back(name);
    clockCh3Names.push_back(ch3name);
    numClock++;
    return iClock;
  }

  /**
   * @brief Reset the data for all clocks
   */
  void reset() {
    for (int i = 0; i < numClock; i++) {
      clockNumCall[i] = 0;
      clockStart[i] = initialClockStart;
      clockTicks[i] = 0;
    }
    startTime = getWallTime();
    startTick = getWallTick();
  }

  /**
   * @brief Start a clock
   */
  void start(
	     int iClock  //!< Index of the clock to be started
  ) {
    assert(iClock >= 0);
    assert(iClock < numClock);
    // Check that the clock's been stopped. It should be set to
    // getWallTick() >= 0 (or initialised to initialClockStart > 0)
#ifdef HiGHSDEV
    if (clockStart[iClock] < 0) {
      printf("recordStart [%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
             iClock, clockStart[iClock], clockTicks[iClock],
             clockNumCall[iClock]);
      fflush(stdout);
    }
#endif
    assert(clockStart[iClock] > 0);
    // Set the start to be the negation of the WallTick to check that
    // the clock's been started in recordFinish
    clockStart[iClock] = -getWallTick();
  }

  /**
   * @brief Stop a clock
   */
  void stop(
	    int iClock  //!< Index of the clock to be stopped
  ) {
    assert(iClock >= 0);
    assert(iClock < numClock);
    // Check that the clock's been started. It should be set to
    // -getWallTick() <= 0
#ifdef HiGHSDEV
    if (clockStart[iClock] > 0)
      printf("recordFinish[%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
             iClock, clockStart[iClock], clockTicks[iClock],
             clockNumCall[iClock]);
#endif
    assert(clockStart[iClock] < 0);
    double wallTick = getWallTick();
    clockTicks[iClock] += (wallTick + clockStart[iClock]);
    clockNumCall[iClock]++;
    // Set the start to be the WallTick to check that the clock's been
    // stopped in recordStart
    clockStart[iClock] = wallTick;
  }

  /**
   * @brief Report timing information for the clock indices in the list
   */
  void report(
	      std::vector<int>&clockList          //!< List of indices to report
  ) {
    double tlPerCentReport = 1.0;
    report_tl(clockList, tlPerCentReport);
  }

  void report_tl(
		 std::vector<int>&clockList,          //!< List of indices to report
		 double tlPerCentReport
  ) {
    const bool reportForExcel = false;
    int numClockListEntries = clockList.size();

    // Report in one line the per-mille contribution from each clock
    double totalTick = getTick();
    printf("TXT-PROFILE-name  ");
    for (int i = 0; i < numClockListEntries; i++) {
      int iClock = clockList[i];
      assert(iClock >= 0);
      assert(iClock < numClock);
      printf(" %-3s", clockCh3Names[iClock].c_str());
    }
    printf("\n");
    printf("TXT-PROFILE-clock ");
    double suPerMille = 0;
    for (int i = 0; i < numClockListEntries; i++) {
      int iClock = clockList[i];
      double perMille = 1000.0 * clockTicks[iClock] / totalTick;
      int int_PerMille = (perMille + 0.5);  // Forcing proper rounding
      if (int_PerMille>0) {
	printf("%4d", int_PerMille); // Just in case one time is 1000!
      } else {
	printf("    "); // Just in case one time is 1000!
      }
      suPerMille += perMille;
    }
    int int_suPerMille = (suPerMille + 0.5);  // Forcing proper rounding
    printf(" per mille: Sum = %d", int_suPerMille);
    printf("\n");
    // Report one line per clock, the time, number of calls and time per call
    printf("TXT-PROFILE-time ");
#ifdef HiGHSDEV
    printf("ID: ");
#endif
    printf("Operation       :    Time             :   Calls   Time/Call\n");
    // Convert approximate seconds
    double tick2sec = 3.6e-10;
    double suTick = 0;
    double suTi = 0;
    for (int i = 0; i < numClockListEntries; i++) {
      int iClock = clockList[i];
      double tick = clockTicks[iClock];
      double ti = tick2sec * tick;
      double perCent = 100.0 * tick / totalTick;
      double tiPerCall = 0;
      if (clockNumCall[iClock] > 0) {
	tiPerCall = ti / clockNumCall[iClock];
	if (perCent >= tlPerCentReport) {
	  printf("TXT-PROFILE-time ");
#ifdef HiGHSDEV
	  printf("%2d: ", iClock);
#endif
	  printf("%-16s: %11.4e (%5.1f%%): %7d %11.4e\n",
		 clockNames[iClock].c_str(), ti, perCent, clockNumCall[iClock],
		 tiPerCall);
	}
      }
      suTi += ti;
      suTick += tick;
    }
    double perCent = 100.0 * suTick / totalTick;
    printf("TXT-PROFILE-time ");
#ifdef HiGHSDEV
    printf("    ");
#endif
    printf("SUM             : %11.4e (%5.1f%%)\n", suTi, perCent);
    printf("TXT-PROFILE-time ");
#ifdef HiGHSDEV
    printf("    ");
#endif
    printf("TOTAL           : %11.4e\n", tick2sec * totalTick);
    if (reportForExcel) {
      // Repeat reporting for Excel
      printf("grep_excel-profile-name");
      for (int i = 0; i < numClockListEntries; i++) {
	int iClock = clockList[i];
	printf(",%s", clockNames[iClock].c_str());
      }
      printf(",SumTime");
      printf(",TotalTime\n");
      printf("grep_excel-profile-time");
      for (int i = 0; i < numClockListEntries; i++) {
	int iClock = clockList[i];
	printf(",%e", tick2sec * clockTicks[iClock]);
      }
      printf(",%e", suTi);
      printf(",%e\n", tick2sec * totalTick);
      printf("grep_excel-profile-calls");
      for (int i = 0; i < numClockListEntries; i++) {
	int iClock = clockList[i];
	printf(",%d", clockNumCall[iClock]);
      }
      printf("\n");
    }
  }

  void reportDualSimplexInnerClock() {report_tl(dualSimplexInnerClockList, 0.0);}
  void reportDualSimplexOuterClock() {report_tl(dualSimplexOuterClockList, 0.0);}
  void reportDualSimplexIterateClock() {report_tl(dualSimplexIterateClockList, 0.0);}
  /**
   * @brief Return the wall-clock time since the clocks were reset
   */
  double getTime() { return getWallTime() - startTime; }

  /**
   * @brief Return the CPU ticks since the clocks were reset
   */
  double getTick() { return getWallTick() - startTick; }

  /**
   * @brief Return the current wall-clock time
   */
  double getWallTime() {
    double walltime;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    walltime = tv.tv_sec;
    walltime += (double)tv.tv_usec / 1000000.0;
    return walltime;
  }

  /**
   * @brief Return the current CPU ticks
   */
  double getWallTick() {
    unsigned a, d;
    asm volatile("rdtsc" : "=a"(a), "=d"(d));
    return ((unsigned long long)a) | (((unsigned long long)d) << 32);
  }

  // Clocks for profiling the dual simplex solver
  int Group1Clock;            //!< Group for SIP
  int IterateClock;           //!< Top level timing of HDual::solve_phase1() and HDual::solve_phase2()
  int IterateRebuildClock;   //!< Second level timing of rebuild()
  int IterateChuzrClock;     //!< Second level timing of CHUZR
  int IterateChuzcClock;     //!< Second level timing of CHUZC
  int IterateFtranClock;     //!< Second level timing of FTRAN
  int IterateVerifyClock;    //!< Second level timing of numerical check
  int IterateDualClock;      //!< Second level timing of dual update
  int IteratePrimalClock;    //!< Second level timing of primal update
  int IterateDevexIzClock;  //!< Second level timing of initialise Devex
  int IteratePivotsClock;    //!< Second level timing of pivoting
  int InvertClock;          //!< Invert in rebuild()
  int PermWtClock;         //!< Permutation of SED weights each side of INVERT in rebuild()
  int ComputeDualClock;    //!< Computation of dual values in rebuild()
  int CorrectDualClock;    //!< Correction of dual values in rebuild()
  int CollectPrIfsClock;  //!< Identification of primal infeasibilities in rebuild()
  int ComputePrimalClock;  //!< Computation of primal values in rebuild()
  int ComputeDuobjClock;   //!< Computation of dual objective value in rebuild()
  int ReportInvertClock;   //!< Reporting of log line in rebuild()
  int Chuzr1Clock;          //!< CHUZR
  int Chuzc0Clock;          //!< CHUZC - stage 0
  int Chuzc1Clock;          //!< CHUZC - stage 1
  int Chuzc2Clock;          //!< CHUZC - stage 2
  int Chuzc3Clock;          //!< CHUZC - stage 3
  int Chuzc4Clock;          //!< CHUZC - stage 4
  int DevexWtClock;        //!< Calculation of Devex weight of entering variable
  int FtranClock;           //!< FTRAN - pivotal column
  int BtranClock;           //!< BTRAN
  int PriceClock;           //!< PRICE
  int FtranDseClock;       //!< FTRAN for DSE weights
  int FtranMixClock;       //!< FTRAN for PAMI
  int FtranBfrtClock;      //!< FTRAN for BFRT
  int UpdateDualClock;     //!< Update of dual values
  int UpdatePrimalClock;   //!< Update of primal values
  int DevexIzClock;        //!< Initialisation of new Devex framework
  int UpdateWeightClock;   //!< Update of DSE or Devex weights
  int UpdatePivotsClock;   //!< Update indices of basic and nonbasic after basis change
  int UpdateFactorClock;   //!< Update the representation of \f$B^{-1}\f$
  int UpdateMatrixClock;   //!< Update the row-wise copy of the constraint matrix for nonbasic columns
  int UpdateRowEpClock;   //!< Update the tableau rows in PAMI

  // private: 
  double startTime;  //!< Elapsed time when the clocks were reset
  double startTick;  //!< CPU ticks when the clocks were reset

  const double initialClockStart = 1.0; //!< Dummy positive start time
					//!for clocks - so they can be
					//!checked as having been
					//!stopped
  int numClock;
  std::vector<int> clockNumCall;
  std::vector<double> clockStart;
  std::vector<double> clockTicks;
  std::vector<std::string> clockNames;
  std::vector<std::string> clockCh3Names;
  std::vector<int> dualSimplexInnerClockList;
  std::vector<int> dualSimplexOuterClockList;
  std::vector<int> dualSimplexIterateClockList;

  /**
   * @brief Return the current CPU ticks
   */
  double initialiseDualSimplexClocks() {
    Group1Clock = clockDef("GROUP1", "GP1");
    IterateClock = clockDef("ITERATE", "ITR");
    IterateRebuildClock = clockDef("REBUILD", "INV");
    IterateChuzrClock = clockDef("CHUZR", "CZR");
    IterateChuzcClock = clockDef("CHUZC", "CZC");
    IterateFtranClock = clockDef("FTRAN", "FTR");
    IterateVerifyClock = clockDef("VERIFY", "VRF");
    IterateDualClock = clockDef("DUAL", "UDU");
    IteratePrimalClock = clockDef("PRIMAL", "UPR");
    IterateDevexIzClock = clockDef("DEVEX_IZ", "DIZ");
    IteratePivotsClock = clockDef("PIVOTS", "PIV");
    InvertClock = clockDef("INVERT", "INV");
    PermWtClock = clockDef("PERM_WT", "PWT");
    ComputeDualClock = clockDef("COMPUTE_DUAL", "CPD");
    CorrectDualClock = clockDef("CORRECT_DUAL", "CRD");
    CollectPrIfsClock = clockDef("COMPUTE_PRIMAL", "CPP");
    ComputePrimalClock = clockDef("COLLECT_PR_IFS", "IFS");
    ComputeDuobjClock = clockDef("COMPUTE_DUOBJ", "DOB");
    ReportInvertClock = clockDef("REPORT_INVERT", "RPI");
    Chuzr1Clock = clockDef("CHUZR1", "CR1");
    Chuzc0Clock = clockDef("CHUZC0", "CC0");
    Chuzc1Clock = clockDef("CHUZC1", "CC1");
    Chuzc2Clock = clockDef("CHUZC2", "CC2");
    Chuzc3Clock = clockDef("CHUZC3", "CC3");
    Chuzc4Clock = clockDef("CHUZC4", "CC4");
    DevexWtClock = clockDef("DEVEX_WT", "DWT");
    FtranClock = clockDef("FTRAN", "COL");
    BtranClock = clockDef("BTRAN", "REP");
    PriceClock = clockDef("PRICE", "RAP");
    FtranDseClock = clockDef("FTRAN_DSE", "DSE");
    FtranMixClock = clockDef("FTRAN_MIX", "MIX");
    FtranBfrtClock = clockDef("FTRAN_BFRT", "BFR");
    UpdateDualClock = clockDef("UPDATE_DUAL", "UPD");
    UpdatePrimalClock = clockDef("UPDATE_PRIMAL", "UPP");
    DevexIzClock = clockDef("UPDATE_WEIGHT", "UPW");
    UpdateWeightClock = clockDef("DEVEX_IZ", "DIZ");
    UpdatePivotsClock = clockDef("UPDATE_PIVOTS", "UPP");
    UpdateFactorClock = clockDef("UPDATE_FACTOR", "UPF");
    UpdateMatrixClock = clockDef("UPDATE_MATRIX", "UPM");
    UpdateRowEpClock = clockDef("UPDATE_ROW_EP", "UPR");

    dualSimplexInnerClockList = {
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

    dualSimplexOuterClockList = {
      IterateRebuildClock, IterateChuzrClock, IterateChuzcClock,
	IterateFtranClock, IterateVerifyClock, IterateDualClock,
	IteratePrimalClock, IterateDevexIzClock, IteratePivotsClock};

    dualSimplexIterateClockList = {
      IterateClock
	};

  }
};
#endif /* SIMPLEX_HIGHSTIMER_H_ */
