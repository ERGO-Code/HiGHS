/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsTimer.h
 * @brief Profiling facility for computational components in HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSTIMER_H_
#define UTIL_HIGHSTIMER_H_

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
    int iClock = clockDef("Run HiGHS","RnH");
    assert(iClock == 0);
    runHighsClock = iClock;
    runHighsClockTime = 0;
    runHighsClockStartTime = initialClockStart;

    presolveClock = clockDef("Presolve", "Pre");
    scaleClock = clockDef("Scale", "Scl");
    crashClock = clockDef("Crash", "Csh");
    solveClock = clockDef("Solve", "Slv");
    postsolveClock = clockDef("Postsolve", "Pst");

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
	     int iClock //!< Index of the clock to be started
  ) {
    assert(iClock >= 0);
    assert(iClock < numClock);
    // Check that the clock's been stopped. It should be set to
    // getWallTick() >= 0 (or initialised to initialClockStart > 0)
#ifdef HiGHSDEV
    if (clockStart[iClock] <= 0) {
      printf("recordStart [%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
             iClock, clockStart[iClock], clockTicks[iClock],
             clockNumCall[iClock]);
      fflush(stdout);
    }
#endif
    assert(clockStart[iClock] > 0);
    // Set the start to be the negation of the WallTick to check that
    // the clock's been started when it's next stopped
    clockStart[iClock] = -getWallTick();
  }

  /**
   * @brief Stop a clock
   */
  void stop(
	    int iClock //!< Index of the clock to be stopped
  ) {
    assert(iClock >= 0);
    assert(iClock < numClock);
    // Check that the clock's been started. It should be set to
    // -getWallTick() <= 0
#ifdef HiGHSDEV
    if (clockStart[iClock] > 0) {
      printf("recordFinish[%2d] is %11.4g: Ticks = %11.4g: NumCall = %d\n",
             iClock, clockStart[iClock], clockTicks[iClock],
             clockNumCall[iClock]);
      fflush(stdout);
    }
#endif
    assert(clockStart[iClock] < 0);
    double wallTick = getWallTick();
    clockTicks[iClock] += (wallTick + clockStart[iClock]);
    clockNumCall[iClock]++;
    // Set the start to be the WallTick to check that the clock's been
    // stopped when it's next started
    clockStart[iClock] = wallTick;
  }

  /**
   * @brief Read the time of a clock
   */
  double read(
	    int iClock //!< Index of the clock to be read
  ) {
    assert(iClock >= 0);
    assert(iClock < numClock);
    double readTick;
    double wallTick;
    if (clockStart[iClock] < 0) {
      // The clock's been started, so find the current time
      wallTick = getWallTick();
      readTick = wallTick + clockStart[iClock];
    } else {
      // The clock is currently stopped, so read the current time
      readTick = clockTicks[iClock];
    }
    double readTime = readTick*tick2sec;
    return readTime;
  }

  /**
   * @brief Start the RunHighs clock
   */
  void startRunHighsClock() {
    start(runHighsClock);
    assert(runHighsClockStartTime > 0);
    double wallTime = getWallTime();
    // Set the clock start to be the negation of WallTime to check that the clock's been
    // started when it's next stopped
    runHighsClockStartTime = -wallTime;
    //    printf("Set runHighsClockStartTime = %g\n", runHighsClockStartTime);
    //    printf("startRunHighsClock() clockTicks = %g; clockStart = %g, runHighsClockStartTime = %g\n",
    //	   clockTicks[runHighsClock], clockStart[runHighsClock], runHighsClockStartTime);
  }

  /**
   * @brief Stop the RunHighs clock
   */
  void stopRunHighsClock() {
    stop(runHighsClock);
    // Get the wall time to update tick2sec
    double wallTime = getWallTime();
    runHighsClockTime += (wallTime + runHighsClockStartTime);
    if (runHighsClockTime > 1e-2) {
      double NWtick2sec = runHighsClockTime/clockTicks[runHighsClock];
      //      printf("Updating tick2sec = %12g to %12g\n", tick2sec, NWtick2sec);
      tick2sec = NWtick2sec;
    }
    // Set the clock start to be the WallTime to check that the clock's been
    // stopped when it's next started
    runHighsClockStartTime = wallTime;
    //    printf("stopRunHighsClock() clockTicks = %g; clockStart = %g, runHighsClockStartTime = %g\n",
    //	   clockTicks[runHighsClock], clockStart[runHighsClock], runHighsClockStartTime);
  }

  /**
   * @brief Read the RunHighs clock
   */
  double readRunHighsClock() {
    int iClock = runHighsClock;
    double readTick;
    double wallTick;
    if (clockStart[iClock] < 0) {
      // The clock's been started, so find the current time
      wallTick = getWallTick();
      readTick = wallTick + clockStart[iClock];

      // Get the wall time to update tick2sec
      double wallTime = getWallTime();
      double currentRunClockTime = runHighsClockTime + (wallTime + runHighsClockStartTime);
      if (currentRunClockTime > 1e-2) {
	double NWtick2sec = currentRunClockTime/readTick;
	//	printf("Updating tick2sec = %12g to %12g/%12g = %12g\n", tick2sec, currentRunClockTime, readTick, NWtick2sec);
	tick2sec = NWtick2sec;
      }
    } else {
      // The clock is currently stopped, so read the current time
      readTick = clockTicks[iClock];
    }
    double readTime = readTick*tick2sec;
    //    printf("readRunHighsClock() clockTicks = %g; clockStart = %g, runHighsClockStartTime = %g\n",
    //	   clockTicks[runHighsClock], clockStart[runHighsClock], runHighsClockStartTime);
    return readTime;
  }

  /**
   * @brief Test whether the RunHighs clock is running
   */
  bool runningRunHighsClock() {return clockStart[runHighsClock] < 0;}

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

    // Check validity of the clock list and check no clocks are still running
    for (int i = 0; i < numClockListEntries; i++) {
      int iClock = clockList[i];
      assert(iClock >= 0);
      assert(iClock < numClock);
      // Check that the clock's not still running. It should be set to
      // getWallTick() >= 0 (or initialised to initialClockStart > 0)
#ifdef HiGHSDEV
      if (clockStart[iClock] <= 0) {
	printf("Clock %2d is still running: Start = %11.4g: Ticks = %11.4g: NumCall = %d\n",
	       iClock, clockStart[iClock], clockTicks[iClock],
	       clockNumCall[iClock]);
	fflush(stdout);
      }
#endif
      assert(clockStart[iClock] > 0);
    }

    // Report in one line the per-mille contribution from each clock
    // First give the 3-character clock names as column headers
    printf("txt-profile-name  ");
    for (int i = 0; i < numClockListEntries; i++) {
      printf(" %-3s", clockCh3Names[clockList[i]].c_str());
    }
    printf("\n");


    // Then give the per-mille contribution relative to the total
    // HiGHS run time, and then relative to the sum of ticks for these
    // clocks
    double currentRunHighsTime = readRunHighsClock();
    double currentRunHighsTick = currentRunHighsTime / tick2sec;
    double suClockTicks = 0;
    for (int passNum = 0; passNum < 2; passNum++) {
      double suPerMille = 0;
      printf("txt-profile-clock ");
      for (int i = 0; i < numClockListEntries; i++) {
	int iClock = clockList[i];
	double perMille;
	if (passNum == 0) {
	  perMille = 1000.0 * clockTicks[iClock] / currentRunHighsTick;
	} else {
	  perMille = 1000.0 * clockTicks[iClock] / suClockTicks;
	}
	int int_PerMille = (perMille + 0.5);  // Forcing proper rounding
	if (int_PerMille>0) {
	  printf("%4d", int_PerMille); // Just in case one time is 1000!
	} else {
	  printf("    "); // Just in case one time is 1000!
	}
	suPerMille += perMille;
	if (passNum == 0) {
	  suClockTicks += clockTicks[iClock];
	}
      }
      int int_suPerMille = (suPerMille + 0.5);  // Forcing proper rounding
      printf(" per mille: Sum = %4d", int_suPerMille);
      printf("\n");
    }

    // Report one line per clock, the time, number of calls and time per call
    printf("txt-profile-time ");
#ifdef HiGHSDEV
    printf("ID: ");
#endif
    printf("Operation       :    Time                     :   Calls   Time/Call\n");
    // Convert approximate seconds
    double suTick = 0;
    double suTi = 0;
    for (int i = 0; i < numClockListEntries; i++) {
      int iClock = clockList[i];
      double tick = clockTicks[iClock];
      double ti = tick2sec * tick;
      double perCentRunHighs = 100.0 * tick / currentRunHighsTick;
      double perCentSumClockTicks = 100.0 * tick / suClockTicks;
      double tiPerCall = 0;
      if (clockNumCall[iClock] > 0) {
	tiPerCall = ti / clockNumCall[iClock];
	if (perCentSumClockTicks >= tlPerCentReport) {
	  printf("txt-profile-time ");
#ifdef HiGHSDEV
	  printf("%2d: ", iClock);
#endif
	  printf("%-16s: %11.4e (%5.1f%%; %5.1f%%): %7d %11.4e\n",
		 clockNames[iClock].c_str(), ti, perCentSumClockTicks, perCentRunHighs, clockNumCall[iClock],
		 tiPerCall);
	}
      }
      suTi += ti;
      suTick += tick;
    }
    double perCentRunHighs = 100.0 * suTick / currentRunHighsTick;
    double perCentSumClockTicks = 100.0;
    printf("txt-profile-time ");
#ifdef HiGHSDEV
    printf("    ");
#endif
    printf("SUM             : %11.4e (%5.1f%%; %5.1f%%)\n", suTi, perCentSumClockTicks, perCentRunHighs);
    printf("txt-profile-time ");
#ifdef HiGHSDEV
    printf("    ");
#endif
    printf("TOTAL           : %11.4e\n", tick2sec * currentRunHighsTick);
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
      printf(",%e\n", tick2sec * currentRunHighsTick);
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

  // private: 
  double startTime;  //!< Elapsed time when the clocks were reset
  double startTick;  //!< CPU ticks when the clocks were reset
  int runHighsClock; //!< The index of the RunHighsClock - should always be 0
  double runHighsClockTime; //!< HiGHS run time - used to scale ticks to time
  double runHighsClockStartTime; //!< HiGHS run start time - used to compute HiGHS run time

  const double initialClockStart = 1.0; //!< Dummy positive start ticks for clocks - so they can be
					//!checked as having been stopped
  int numClock;
  std::vector<int> clockNumCall;
  std::vector<double> clockStart;
  std::vector<double> clockTicks;
  std::vector<std::string> clockNames;
  std::vector<std::string> clockCh3Names;
  double tick2sec = 3.6e-10;

  // Fundamental Highs clocks
  int presolveClock;
  int scaleClock;
  int crashClock;
  int solveClock;
  int postsolveClock;
};



#endif /* UTIL_HIGHSTIMER_H_ */
