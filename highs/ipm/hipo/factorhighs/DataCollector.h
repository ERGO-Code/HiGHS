#ifndef FACTORHIGHS_DATA_COLLECTOR_H
#define FACTORHIGHS_DATA_COLLECTOR_H

#include <mutex>
#include <vector>

#include "FactorHighsSettings.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/Logger.h"

namespace hipo {

// DataCollector collects times and numbers of BLAS calls during analyse,
// factorise and solve, for debugging. They are collected only if
// HIPO_TIMING_LEVEL > 0 at compile time; otherwise the functions are empty and
// there is no performance hit. Collection uses a mutex, so that it can be done
// concurrently.

class DataCollector {
  // Record of times and BLAS calls
  std::vector<double> times_{};
  std::vector<Int> blas_calls_{};

  // Mutex for concurrent access
  std::mutex mutex_;

 public:
  DataCollector();

  void sumTime(TimeItems i, double t);
  void printTimes(const Logger& logger) const;
};

}  // namespace hipo

#endif