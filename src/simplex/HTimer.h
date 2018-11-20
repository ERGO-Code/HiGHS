/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HTimer.h
 * @brief Profiling facility for computational components in HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HTIMER_H_
#define SIMPLEX_HTIMER_H_

#include <sys/time.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>
using namespace std;

/**
 * Define the indices of clocks used to time computation
 */
enum HTickItem {
  HTICK_GROUP1 = 0,        //!< Group for SIP
  HTICK_ITERATE,           //!< Top level timing of HDual::solve_phase1() and
                           //!< HDual::solve_phase2()
  HTICK_ITERATE_REBUILD,   //!< Second level timing of rebuild()
  HTICK_ITERATE_CHUZR,     //!< Second level timing of CHUZR
  HTICK_ITERATE_CHUZC,     //!< Second level timing of CHUZC
  HTICK_ITERATE_FTRAN,     //!< Second level timing of FTRAN
  HTICK_ITERATE_VERIFY,    //!< Second level timing of numerical check
  HTICK_ITERATE_DUAL,      //!< Second level timing of dual update
  HTICK_ITERATE_PRIMAL,    //!< Second level timing of primal update
  HTICK_ITERATE_DEVEX_IZ,  //!< Second level timing of initialise Devex
  HTICK_ITERATE_PIVOTS,    //!< Second level timing of pivoting

  HTICK_INVERT,          //!< INVERT in rebuild()
  HTICK_PERM_WT,         //!< Permutation of SED weights each side of INVERT in
                         //!< rebuild()
  HTICK_COMPUTE_DUAL,    //!< Computation of dual values in rebuild()
  HTICK_CORRECT_DUAL,    //!< Correction of dual values in rebuild()
  HTICK_COLLECT_PR_IFS,  //!< Identification of primal infeasibilities in
                         //!< rebuild()
  HTICK_COMPUTE_PRIMAL,  //!< Computation of primal values in rebuild()
  HTICK_COMPUTE_DUOBJ,   //!< Computation of dual objective value in rebuild()
  HTICK_REPORT_INVERT,   //!< Reporting of log line in rebuild()
  HTICK_CHUZR1,          //!< CHUZR
  HTICK_CHUZC0,          //!< CHUZC - stage 0
  HTICK_CHUZC1,          //!< CHUZC - stage 1
  HTICK_CHUZC2,          //!< CHUZC - stage 2
  HTICK_CHUZC3,          //!< CHUZC - stage 3
  HTICK_CHUZC4,          //!< CHUZC - stage 4
  HTICK_DEVEX_WT,        //!< Calculation of Devex weight of entering variable
  HTICK_FTRAN,           //!< FTRAN - pivotal column
  HTICK_BTRAN,           //!< BTRAN
  HTICK_PRICE,           //!< PRICE
  HTICK_FTRAN_DSE,       //!< FTRAN for DSE weights
  HTICK_FTRAN_MIX,       //!< FTRAN for PAMI
  HTICK_FTRAN_BFRT,      //!< FTRAN for BFRT
  HTICK_UPDATE_DUAL,     //!< Update of dual values
  HTICK_UPDATE_PRIMAL,   //!< Update of primal values
  HTICK_DEVEX_IZ,        //!< Initialisation of new Devex framework
  HTICK_UPDATE_WEIGHT,   //!< Update of DSE or Devex weights
  HTICK_UPDATE_PIVOTS,   //!< Update indices of basic and nonbasic after basis
                         //!< change
  HTICK_UPDATE_FACTOR,   //!< Update the representation of \f$B^{-1}\f$
  HTICK_UPDATE_MATRIX,   //!< Update the row-wise copy of the constraint matrix
                         //!< for nonbasic columns
  HTICK_UPDATE_ROW_EP,   //!< Update the tableau rows in PAMI
  HTICK_ITEMS_COUNT      //!< Total number of clocks
};

/**
 * @brief Class for profiling facility for computational components in HiGHS
 */
class HTimer {
 public:
  /**
   * @brief Set up full and 3-character names for each clock
   */
  HTimer() {
    itemNames[HTICK_GROUP1] = "GROUP1";
    itemCh3Names[HTICK_GROUP1] = "GP1";

    itemNames[HTICK_ITERATE] = "ITERATE";
    itemCh3Names[HTICK_ITERATE] = "ITR";
    itemNames[HTICK_ITERATE_REBUILD] = "REBUILD";
    itemCh3Names[HTICK_ITERATE_REBUILD] = "INV";
    itemNames[HTICK_ITERATE_CHUZR] = "CHUZR";
    itemCh3Names[HTICK_ITERATE_CHUZR] = "CZR";
    itemNames[HTICK_ITERATE_CHUZC] = "CHUZC";
    itemCh3Names[HTICK_ITERATE_CHUZC] = "CZC";
    itemNames[HTICK_ITERATE_FTRAN] = "FTRAN";
    itemCh3Names[HTICK_ITERATE_FTRAN] = "FTR";
    itemNames[HTICK_ITERATE_VERIFY] = "VERIFY";
    itemCh3Names[HTICK_ITERATE_VERIFY] = "VRF";
    itemNames[HTICK_ITERATE_DUAL] = "DUAL";
    itemCh3Names[HTICK_ITERATE_DUAL] = "UDU";
    itemNames[HTICK_ITERATE_PRIMAL] = "PRIMAL";
    itemCh3Names[HTICK_ITERATE_PRIMAL] = "UPR";
    itemNames[HTICK_ITERATE_DEVEX_IZ] = "DEVEX_IZ";
    itemCh3Names[HTICK_ITERATE_DEVEX_IZ] = "DIZ";
    itemNames[HTICK_ITERATE_PIVOTS] = "PIVOTS";
    itemCh3Names[HTICK_ITERATE_PIVOTS] = "PIV";

    itemNames[HTICK_INVERT] = "INVERT";
    itemCh3Names[HTICK_INVERT] = "INV";
    itemNames[HTICK_PERM_WT] = "PERM_WT";
    itemCh3Names[HTICK_PERM_WT] = "PWT";
    itemNames[HTICK_COMPUTE_DUAL] = "COMPUTE_DUAL";
    itemCh3Names[HTICK_COMPUTE_DUAL] = "CPD";
    itemNames[HTICK_CORRECT_DUAL] = "CORRECT_DUAL";
    itemCh3Names[HTICK_CORRECT_DUAL] = "CRD";
    itemNames[HTICK_COMPUTE_PRIMAL] = "COMPUTE_PRIMAL";
    itemCh3Names[HTICK_COMPUTE_PRIMAL] = "CPP";
    itemNames[HTICK_COLLECT_PR_IFS] = "COLLECT_PR_IFS";
    itemCh3Names[HTICK_COLLECT_PR_IFS] = "IFS";
    itemNames[HTICK_COMPUTE_DUOBJ] = "COMPUTE_DUOBJ";
    itemCh3Names[HTICK_COMPUTE_DUOBJ] = "DOB";
    itemNames[HTICK_REPORT_INVERT] = "REPORT_INVERT";
    itemCh3Names[HTICK_REPORT_INVERT] = "RPI";
    itemNames[HTICK_CHUZR1] = "CHUZR1";
    itemCh3Names[HTICK_CHUZR1] = "CR1";
    itemNames[HTICK_CHUZC0] = "CHUZC0";
    itemCh3Names[HTICK_CHUZC0] = "CC0";
    itemNames[HTICK_CHUZC1] = "CHUZC1";
    itemCh3Names[HTICK_CHUZC1] = "CC1";
    itemNames[HTICK_CHUZC2] = "CHUZC2";
    itemCh3Names[HTICK_CHUZC2] = "CC2";
    itemNames[HTICK_CHUZC3] = "CHUZC3";
    itemCh3Names[HTICK_CHUZC3] = "CC3";
    itemNames[HTICK_CHUZC4] = "CHUZC4";
    itemCh3Names[HTICK_CHUZC4] = "CC4";
    itemNames[HTICK_DEVEX_WT] = "DEVEX_WT";
    itemCh3Names[HTICK_DEVEX_WT] = "DWT";
    itemNames[HTICK_FTRAN] = "FTRAN";
    itemCh3Names[HTICK_FTRAN] = "COL";
    itemNames[HTICK_BTRAN] = "BTRAN";
    itemCh3Names[HTICK_BTRAN] = "REP";
    itemNames[HTICK_PRICE] = "PRICE";
    itemCh3Names[HTICK_PRICE] = "RAP";
    itemNames[HTICK_FTRAN_DSE] = "FTRAN_DSE";
    itemCh3Names[HTICK_FTRAN_DSE] = "DSE";
    itemNames[HTICK_FTRAN_MIX] = "FTRAN_MIX";
    itemCh3Names[HTICK_FTRAN_MIX] = "MIX";
    itemNames[HTICK_FTRAN_BFRT] = "FTRAN_BFRT";
    itemCh3Names[HTICK_FTRAN_BFRT] = "BFR";
    itemNames[HTICK_UPDATE_DUAL] = "UPDATE_DUAL";
    itemCh3Names[HTICK_UPDATE_DUAL] = "UPD";
    itemNames[HTICK_UPDATE_PRIMAL] = "UPDATE_PRIMAL";
    itemCh3Names[HTICK_UPDATE_PRIMAL] = "UPP";
    itemNames[HTICK_UPDATE_WEIGHT] = "UPDATE_WEIGHT";
    itemCh3Names[HTICK_UPDATE_WEIGHT] = "UPW";
    itemNames[HTICK_DEVEX_IZ] = "DEVEX_IZ";
    itemCh3Names[HTICK_DEVEX_IZ] = "DIZ";
    itemNames[HTICK_UPDATE_PIVOTS] = "UPDATE_PIVOTS";
    itemCh3Names[HTICK_UPDATE_PIVOTS] = "UPP";
    itemNames[HTICK_UPDATE_FACTOR] = "UPDATE_FACTOR";
    itemCh3Names[HTICK_UPDATE_FACTOR] = "UPF";
    itemNames[HTICK_UPDATE_MATRIX] = "UPDATE_MATRIX";
    itemCh3Names[HTICK_UPDATE_MATRIX] = "UPM";
    itemNames[HTICK_UPDATE_ROW_EP] = "UPDATE_ROW_EP";
    itemCh3Names[HTICK_UPDATE_ROW_EP] = "UPR";

    reset();
  }

  /**
   * @brief Reset the data for all clocks
   */
  void reset() {
    for (int i = 0; i < HTICK_ITEMS_COUNT; i++) {
      itemNumCall[i] = 0;
      itemStart[i] = initialItemStart;
      itemTicks[i] = 0;
    }
    startTime = getWallTime();
    startTick = getWallTick();
  }

  /**
   * @brief Start the clock for an item
   */
  void recordStart(int itemKey  //!< Index of the clock to be started
  ) {
    // Check that the clock's been stopped. It should be set to
    // getWallTick() >= 0 (or initialised to initialItemStart > 0)
#ifdef HiGHSDEV
    if (itemStart[itemKey] < 0) {
      printf("recordStart [%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
             itemKey, itemStart[itemKey], itemTicks[itemKey],
             itemNumCall[itemKey]);
      fflush(stdout);
    }
#endif
    assert(itemStart[itemKey] > 0);
    // Set the start to be the negation of the WallTick to check that
    // the clock's been started in recordFinish
    itemStart[itemKey] = -getWallTick();
  }

  /**
   * @brief Stop the clock for an item
   */
  void recordFinish(int itemKey  //!< Index of the clock to be stopped
  ) {
    // Check that the clock's been started. It should be set to
    // -getWallTick() <= 0
#ifdef HiGHSDEV
    if (itemStart[itemKey] > 0)
      printf("recordFinish[%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
             itemKey, itemStart[itemKey], itemTicks[itemKey],
             itemNumCall[itemKey]);
#endif
    assert(itemStart[itemKey] < 0);
    double wallTick = getWallTick();
    itemTicks[itemKey] += (wallTick + itemStart[itemKey]);
    itemNumCall[itemKey]++;
    // Set the start to be the WallTick to check that the clock's been
    // stopped in recordStart
    itemStart[itemKey] = wallTick;
  }

  /**
   * @brief Report timing information for the clock indices in the list
   */
  void report(int itemCount,          //!< Number of indices to report
              int *itemList,          //!< List of indices to report
              double tlPerCentReport  //!< Threshhold of percentage of
                                      //! overall time for a clock to be
                                      //! reported in detail
  ) {
    // Report in one line the per-mille contribution from each clock
    double totalTick = getTick();
    printf("txt-profile-name ");
    for (int i = 0; i < itemCount; i++)
      printf(" %-3s", itemCh3Names[itemList[i]].c_str());
    printf("\n");
    printf("txt-profile-item ");
    double suPerMille = 0;
    for (int i = 0; i < itemCount; i++) {
      double perMille = 1000.0 * itemTicks[itemList[i]] / totalTick;
      int int_PerMille = (perMille + 0.5);  // Forcing proper rounding
      printf(" %3d", int_PerMille);
      suPerMille += perMille;
    }
    int int_suPerMille = suPerMille;
    printf(" per mille: Sum = %d", int_suPerMille);
    printf("\n");
    // Report one line per clock, the time, number of calls and time per call
    printf(
        "txt-profile-time ID: Operation       :    Time            :   Calls   "
        "Time/Call\n");
    // Convert approximate seconds
    double tick2sec = 3.6e-10;
    double suTick = 0;
    double suTi = 0;
    for (int i = 0; i < itemCount; i++) {
      int item = itemList[i];
      double tick = itemTicks[item];
      double ti = tick2sec * tick;
      double perCent = 100.0 * tick / totalTick;
      double tiPerCall = 0;
      if (itemNumCall[item] > 0) tiPerCall = ti / itemNumCall[item];
      if (perCent >= tlPerCentReport) {
        printf("txt-profile-time %2d: %-16s: %11.4e (%4.1f%%): %7d %11.4e\n",
               item, itemNames[item].c_str(), ti, perCent, itemNumCall[item],
               tiPerCall);
      }
      suTi += ti;
      suTick += tick;
    }
    double perCent = 100.0 * suTick / totalTick;
    printf("txt-profile-time   : SUM             : %11.4e (%4.1f%%)\n", suTi,
           perCent);
    printf("txt-profile-time   : TOTAL           : %11.4e\n",
           tick2sec * totalTick);
    // Repeat reporting for Excel
    printf("grep_excel-profile-name");
    for (int i = 0; i < itemCount; i++)
      printf(",%s", itemNames[itemList[i]].c_str());
    printf(",SumTime");
    printf(",TotalTime\n");
    printf("grep_excel-profile-time");
    for (int i = 0; i < itemCount; i++)
      printf(",%e", tick2sec * itemTicks[itemList[i]]);
    printf(",%e", suTi);
    printf(",%e\n", tick2sec * totalTick);
    printf("grep_excel-profile-calls");
    for (int i = 0; i < itemCount; i++) printf(",%d", itemNumCall[itemList[i]]);
    printf("\n");
  }

  /**
   * @brief Return the elapsed time since the clocks were reset
   */
  double getTime() { return getWallTime() - startTime; }

  /**
   * @brief Return the CPU ticks since the clocks were reset
   */
  double getTick() { return getWallTick() - startTick; }

  double startTime;  //!< Elapsed time when the clocks were reset
  double startTick;  //!< CPU ticks when the clocks were reset

  const double initialItemStart =
      1.0;  //!< Dummy positive start time for clocks - so they can be checked
            //!< as having been stopped
  int itemNumCall[HTICK_ITEMS_COUNT];   //!< Number of times each clock has been
                                        //!< run
  double itemStart[HTICK_ITEMS_COUNT];  //!< When running: negation of start
                                        //!< time. When not running: WallTick
                                        //!< when stopped
  double itemTicks[HTICK_ITEMS_COUNT];  //!< Total number of CPU ticks for each
                                        //!< clock
  string itemNames[HTICK_ITEMS_COUNT];  //!< Full name of each clock
  string itemCh3Names[HTICK_ITEMS_COUNT];  //!< 3-character name for each clock

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
};

#endif /* SIMPLEX_HTIMER_H_ */
