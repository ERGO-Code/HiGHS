/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SCHEDULER_CONSTANTS_H_
#define HIGHS_SCHEDULER_CONSTANTS_H_

struct HighsSchedulerConstants {
  enum Constants {
    kNumTryFac = 16,
    kMicroSecsBeforeSleep = 5000,
    kMicroSecsBeforeGlobalSync = 1000,
  };
};

#endif
