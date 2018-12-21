/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HTimerPre.h
 * @brief Profiling facility for presolve in HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef SIMPLEX_HTIMERPRE_H_
#define SIMPLEX_HTIMERPRE_H_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <sys/time.h>

using std::string;

enum HTickItemPre {
  // Presolve items
  EMPTY_ROW,
  FIXED_COL,
  SING_ROW,
  DOUBLETON_EQUATION,
  FORCING_ROW,
  REDUNDANT_ROW,
  DOMINATED_ROW_BOUNDS,
  FREE_SING_COL,
  SING_COL_DOUBLETON_INEQ,
  IMPLIED_FREE_SING_COL,
  DOMINATED_COLS,
  WEAKLY_DOMINATED_COLS,
  DOMINATED_COL_BOUNDS,
  EMPTY_COL,
  // HTICK_PRE_DUPLICATE_ROWS,
  // HTICK_PRE_DUPLICATE_COLUMNS,

  // The total count so far:
  HTICK_ITEMS_COUNT_PRE,

  // Items required by postsolve
  DOUBLETON_EQUATION_ROW_BOUNDS_UPDATE,
  DOUBLETON_EQUATION_X_ZERO_INITIALLY,
  DOUBLETON_EQUATION_NEW_X_NONZERO,
  DOUBLETON_EQUATION_NEW_X_ZERO_AR_UPDATE,
  DOUBLETON_EQUATION_NEW_X_ZERO_A_UPDATE,
  SING_COL_DOUBLETON_INEQ_SECOND_SING_COL,
  FORCING_ROW_VARIABLE,
};

class HTimerPre {
 public:
  HTimerPre() {
    itemNames.resize(HTICK_ITEMS_COUNT_PRE);

    itemNames[EMPTY_ROW] = "emptR";
    itemNames[FIXED_COL] = "fixCo";
    itemNames[DOUBLETON_EQUATION] = "douEq";
    itemNames[SING_ROW] = "singR";
    itemNames[FORCING_ROW] = "forcR";
    itemNames[REDUNDANT_ROW] = "reduR";
    itemNames[FREE_SING_COL] = "colSi";
    itemNames[SING_COL_DOUBLETON_INEQ] = "colSD";
    itemNames[IMPLIED_FREE_SING_COL] = "colSI";
    itemNames[DOMINATED_COLS] = "domCo";
    itemNames[WEAKLY_DOMINATED_COLS] = "wdomC";
    itemNames[EMPTY_COL] = "emptC";

    // itemNames[HTICK_PRE_DUPLICATE_ROWS] = "duplR";
    // itemNames[HTICK_PRE_DUPLICATE_COLUMNS] = "duplC";

    reset();
  }

  // Reset the timer
  void reset() {
    for (int i = 0; i < HTICK_ITEMS_COUNT_PRE; i++) itemTicks[i] = 0;
    startTime = getWallTime();
    startTick = getWallTick();
  }

  // Record tick of one item
  void recordStart(int itemKey) { itemStart[itemKey] = getWallTime(); }
  void recordFinish(int itemKey) {
    itemTicks[itemKey] += (getWallTime() - itemStart[itemKey]);
  }

  // Report specific items in the list
  void report(int itemCount, int *itemList) {
    double totalTick = getTick();
    printf("Presolve times ");
    for (int i = 0; i < itemCount; i++)
      printf(" %s", itemNames[itemList[i]].c_str());
    printf("\n");
    printf("profile-item ");
    for (int i = 0; i < itemCount; i++) {
      int percent = 1000.0 * itemTicks[itemList[i]] / totalTick;
      printf(" %5d", percent);
    }
    printf("\n");
  }

  // Current eclipsed time
  double getTime() { return getWallTime() - startTime; }

  // Current eclipsed CPU tick
  double getTick() { return getWallTick() - startTick; }

  double itemTicks[HTICK_ITEMS_COUNT_PRE];
  vector<string> itemNames;

 private:
  // The start tick and time
  double startTime;
  double startTick;

  // The items tick and name
  double itemStart[HTICK_ITEMS_COUNT_PRE];

  // Current wall time
  double getWallTime() {
    double walltime;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    walltime = tv.tv_sec;
    walltime += (double)tv.tv_usec / 1000000.0;
    return walltime;
  }

  // Current wall tick
  double getWallTick() {
    unsigned a, d;
    asm volatile("rdtsc" : "=a"(a), "=d"(d));
    return ((unsigned long long)a) | (((unsigned long long)d) << 32);
  }
};

#endif /* SIMPLEX_HTIMERPRE_H_ */
