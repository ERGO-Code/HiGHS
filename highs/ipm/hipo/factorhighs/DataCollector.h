#ifndef FACTORHIGHS_DATA_COLLECTOR_H
#define FACTORHIGHS_DATA_COLLECTOR_H

#include <atomic>
#include <limits>
#include <mutex>
#include <vector>

#include "FactorHiGHSSettings.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

struct IterData {
// data of a given ipm iteration

// factorisation data
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  double minD = std::numeric_limits<double>::max();
  double maxD = 0.0;
  double minL = std::numeric_limits<double>::max();
  double maxL = 0.0;
  double max_reg = 0.0;
  Int n_reg_piv = 0;
  Int n_swap = 0;
  Int n_2x2 = 0;
  Int n_wrong_sign = 0;
  double max_wrong_sign = 0.0;
  double nw_back_err{};
  double cw_back_err{};
  Int large_components_cw{};
#endif

  // ipm data
  double sigma_aff;
  double sigma;
  Int correctors;
  Int num_solves = 0;
  double min_theta;
  double max_theta;
  double min_prod;
  double max_prod;
  Int num_small_prod;
  Int num_large_prod;
  double omega{};
  double M_norm1;
  double M_maxdiag;
};

struct CounterData {
  // Record of times and BLAS calls
  std::vector<double> times{};
  std::vector<Int> blas_calls{};
};

// DataCollector is used to collect debug data during the ipm and factorisation.
// DataCollector is a singleton object. Only one copy of it can exist and it
// does not have a public constructor or destructor. Use:
// - DataCollector::initialise() to allocate the DataCollector.
// - DataCollector::terminate() to deallocate the DataCollector.
// - DataCollector::get()->... to access any non-static member function.
//
// Expensive data related to the factorisation is only collected if
// HIPO_COLLECT_EXPENSIVE_DATA is defined. "Cheaper" data related to the ipm is
// always collected and printed only if log_dev_level is high enough.
//
// HIPO_COLLECT_EXPENSIVE_DATA is selected at compile time, because some data is
// collected many times, during many steps of the factorisation, using locks; if
// HIPO_COLLECT_EXPENSIVE_DATA is off at compile time, then there is no
// performance hit, since the functions are empty and are completely ignored by
// the compiler.
//
// Times are only collected if HIPO_TIMING_LEVEL is defined at compile time, for
// the same reason.

class DataCollector {
  // Record of times and BLAS calls
  CounterData counter_data_;

  // record of data of ipm iterations
  std::vector<IterData> iter_data_record_{};

  // Mutexes for concurrent access
  std::mutex times_mutex_;
  std::mutex iter_data_mutex_;

  // Instance of DataCollector
  static DataCollector* ptr_;

  // Private ctor and dtor
  DataCollector();
  ~DataCollector() = default;

 public:
  // Access to the object
  static DataCollector* get();
  static void initialise();
  static void terminate();

  // ========================================================================
  // The functions below can only be accessed via DataCollector::get()->
  // ========================================================================

  IterData& back();
  void append();

  // Functions with lock, they can be accessed simultaneously
  void sumTime(TimeItems i, double t);
  void countRegPiv();
  void countSwap();
  void count2x2();
  void setWrongSign(double p);
  void setMaxReg(double new_reg);
  void setExtremeEntries(double minD, double maxD, double minoffD,
                         double maxoffD);

  // Functions without lock
  void countSolves();
  void setOmega(double omega);
  void setNorms(double norm1, double maxdiag);
  void setSigma(double sigma, bool affinescaling = false);
  void setCorrectors(Int correctors);
  void setBackError(double nw, double cw, Int large_components);
  void setExtremeTheta(const std::vector<double>& scaling);
  void setProducts(double min_prod, double max_prod, Int num_small,
                   Int num_large);

  // Const functions
  void printTimes() const;
  void printIter() const;
};

}  // namespace hipo

#endif