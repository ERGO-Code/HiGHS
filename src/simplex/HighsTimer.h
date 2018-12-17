/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
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
	      int *clockList          //!< List of indices to report
  ) {
    const bool reportForExcel = false;
    int numClockListEntries = sizeof(clockList) / sizeof(int);
    double tlPerCentReport = 0.1;

    printf("report: clockList[] = {"); for (int i = 0; i < numClockListEntries; i++) {printf(" %d", clockList[i]);} printf("}\n");

    // Report in one line the per-mille contribution from each clock
    double totalTick = getTick();
    printf("txt-profile-name  ");
    for (int i = 0; i < numClockListEntries; i++) {
      int iClock = clockList[i];
      assert(iClock >= 0);
      assert(iClock < numClock);
      printf(" %-3s", clockCh3Names[iClock].c_str());
    }
    printf("\n");
    printf("txt-profile-clock ");
    double suPerMille = 0;
    for (int i = 0; i < numClockListEntries; i++) {
      int iClock = clockList[i];
      double perMille = 1000.0 * clockTicks[iClock] / totalTick;
      int int_PerMille = (perMille + 0.5);  // Forcing proper rounding
      printf(" %3d", int_PerMille);
      suPerMille += perMille;
    }
    int int_suPerMille = suPerMille;
    printf(" per mille: Sum = %d", int_suPerMille);
    printf("\n");
    // Report one line per clock, the time, number of calls and time per call
    printf(
	   "txt-profile-time ID: Operation       :    Time             :   Calls   "
	   "Time/Call\n");
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
      if (clockNumCall[iClock] > 0) tiPerCall = ti / clockNumCall[iClock];
      if (perCent >= tlPerCentReport) {
        printf("txt-profile-time %2d: %-16s: %11.4e (%5.1f%%): %7d %11.4e\n",
               iClock, clockNames[iClock].c_str(), ti, perCent, clockNumCall[iClock],
               tiPerCall);
      }
      suTi += ti;
      suTick += tick;
    }
    double perCent = 100.0 * suTick / totalTick;
    printf("txt-profile-time   : SUM             : %11.4e (%5.1f%%)\n", suTi,
           perCent);
    printf("txt-profile-time   : TOTAL           : %11.4e\n",
           tick2sec * totalTick);
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

 private: 
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
};

#endif /* SIMPLEX_HIGHSTIMER_H_ */
