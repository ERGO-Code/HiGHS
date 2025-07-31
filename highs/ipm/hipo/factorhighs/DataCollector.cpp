#include "DataCollector.h"

#include "FactorHiGHSSettings.h"
#include "ipm/hipo/auxiliary/Log.h"

namespace hipo {

// instance of DataCollector
DataCollector* DataCollector::ptr_ = nullptr;

// Functions to manage DataCollector

DataCollector::DataCollector() {
  times.resize(kTimeSize);
  blas_calls.resize(kTimeBlasEnd - kTimeBlasStart + 1);
}
DataCollector* DataCollector::get() { return ptr_; }
void DataCollector::initialise() {
  if (!ptr_) ptr_ = new DataCollector();
}
void DataCollector::terminate() {
  delete ptr_;
  ptr_ = nullptr;
}
void DataCollector::append() {
  // add an empty IterData object to the record
  iter_data_record_.push_back(IterData());
}
IterData& DataCollector::back() {
  // access most recent record of data
  return iter_data_record_.back();
}

//
// Data-collecting functions
//

void DataCollector::sumTime(TimeItems i, double t) {
#if HIPO_TIMING_LEVEL > 0
  // Keep track of times and blas calls.
  std::lock_guard<std::mutex> lock(mutex_);
  times[i] += t;
#if HIPO_TIMING_LEVEL >= 3
  if (i >= kTimeBlasStart && i <= kTimeBlasEnd)
    ++blas_calls[i - kTimeBlasStart];
#endif
#endif
}
void DataCollector::setExtremeEntries(double minD, double maxD, double minoffD,
                                      double maxoffD) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  // Store max and min entries of D and L.
  std::lock_guard<std::mutex> lock(mutex_);
  back().minD = std::min(back().minD, minD);
  back().maxD = std::max(back().maxD, maxD);
  back().minL = std::min(back().minL, minoffD);
  back().maxL = std::max(back().maxL, maxoffD);
#endif
}
void DataCollector::countRegPiv() {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  // Increase the number of dynamically regularised pivots.
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_reg_piv;
#endif
}
void DataCollector::countSwap() {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_swap;
#endif
}
void DataCollector::count2x2() {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_2x2;
#endif
}
void DataCollector::setWrongSign(double p) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  ++back().n_wrong_sign;
  back().max_wrong_sign = std::max(back().max_wrong_sign, std::abs(p));
#endif
}
void DataCollector::setMaxReg(double new_reg) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  // Keep track of maximum regularisation used.
  std::lock_guard<std::mutex> lock(mutex_);
  back().max_reg = std::max(back().max_reg, new_reg);
#endif
}

void DataCollector::setNorms(double norm1, double maxdiag) {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(mutex_);
  back().M_norm1 = norm1;
  back().M_maxdiag = maxdiag;
#endif
}

//
// Printing
//

void DataCollector::printTimes(const Log* log) const {
#if HIPO_TIMING_LEVEL >= 1

  std::stringstream log_stream;

  log_stream << "----------------------------------------------------\n";
  log_stream << "Analyse time            \t" << fix(times[kTimeAnalyse], 8, 4)
             << '\n';

#if HIPO_TIMING_LEVEL >= 2

  log_stream << "\tMetis:                  "
             << fix(times[kTimeAnalyseMetis], 8, 4) << " ("
             << fix(times[kTimeAnalyseMetis] / times[kTimeAnalyse] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tTree:                   "
             << fix(times[kTimeAnalyseTree], 8, 4) << " ("
             << fix(times[kTimeAnalyseTree] / times[kTimeAnalyse] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tCounts:                 "
             << fix(times[kTimeAnalyseCount], 8, 4) << " ("
             << fix(times[kTimeAnalyseCount] / times[kTimeAnalyse] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tSupernodes:             " << fix(times[kTimeAnalyseSn], 8, 4)
             << " ("
             << fix(times[kTimeAnalyseSn] / times[kTimeAnalyse] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tReorder:                "
             << fix(times[kTimeAnalyseReorder], 8, 4) << " ("
             << fix(times[kTimeAnalyseReorder] / times[kTimeAnalyse] * 100, 4,
                    1)
             << "%)\n";
  log_stream << "\tSn sparsity pattern:    "
             << fix(times[kTimeAnalysePattern], 8, 4) << " ("
             << fix(times[kTimeAnalysePattern] / times[kTimeAnalyse] * 100, 4,
                    1)
             << "%)\n";
  log_stream << "\tRelative indices:       "
             << fix(times[kTimeAnalyseRelInd], 8, 4) << " ("
             << fix(times[kTimeAnalyseRelInd] / times[kTimeAnalyse] * 100, 4, 1)
             << "%)\n";
#endif

  log_stream << "----------------------------------------------------\n";
  log_stream << "Factorise time          \t" << fix(times[kTimeFactorise], 8, 4)
             << "\n";

#if HIPO_TIMING_LEVEL >= 2
  log_stream << "\tPrepare:                "
             << fix(times[kTimeFactorisePrepare], 8, 4) << " ("
             << fix(times[kTimeFactorisePrepare] / times[kTimeFactorise] * 100,
                    4, 1)
             << "%)\n";
  log_stream << "\tAssemble children in F: "
             << fix(times[kTimeFactoriseAssembleOriginal], 8, 4) << " ("
             << fix(times[kTimeFactoriseAssembleOriginal] /
                        times[kTimeFactorise] * 100,
                    4, 1)
             << "%)\n";
  log_stream << "\tAssemble children in C: "
             << fix(times[kTimeFactoriseAssembleChildrenFrontal], 8, 4) << " ("
             << fix(times[kTimeFactoriseAssembleChildrenFrontal] /
                        times[kTimeFactorise] * 100,
                    4, 1)
             << "%)\n";
  log_stream << "\tDense factorisation:    "
             << fix(times[kTimeFactoriseDenseFact], 8, 4) << " ("
             << fix(times[kTimeFactoriseDenseFact] / times[kTimeFactorise] *
                        100,
                    4, 1)
             << "%)\n";

  log_stream << "\tTerminate:              "
             << fix(times[kTimeFactoriseTerminate], 8, 4) << " ("
             << fix(times[kTimeFactoriseTerminate] / times[kTimeFactorise] *
                        100,
                    4, 1)
             << "%)\n";

  log_stream << "\t\tmain:           " << fix(times[kTimeDenseFact_main], 8, 4)
             << "\n";
  log_stream << "\t\tSchur:          " << fix(times[kTimeDenseFact_schur], 8, 4)
             << "\n";
  log_stream << "\t\tkernel:         "
             << fix(times[kTimeDenseFact_kernel], 8, 4) << "\n";
  log_stream << "\t\tconvert:        "
             << fix(times[kTimeDenseFact_convert], 8, 4) << "\n";
  log_stream << "\t\tpivoting:       "
             << fix(times[kTimeDenseFact_pivoting], 8, 4) << "\n";
#endif

  log_stream << "----------------------------------------------------\n";
  log_stream << "Solve time              \t" << fix(times[kTimeSolve], 8, 4)
             << "\n";

#if HIPO_TIMING_LEVEL >= 2
  log_stream << "\tPrepare:                "
             << fix(times[kTimeSolvePrepare], 8, 4) << " ("
             << fix(times[kTimeSolvePrepare] / times[kTimeSolve] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tSolve:                  "
             << fix(times[kTimeSolveSolve], 8, 4) << " ("
             << fix(times[kTimeSolveSolve] / times[kTimeSolve] * 100, 4, 1)
             << "%)\n";
  log_stream << "\t\tdense:          "
             << fix(times[kTimeSolveSolve_dense], 8, 4) << "\n";
  log_stream << "\t\tsparse:         "
             << fix(times[kTimeSolveSolve_sparse], 8, 4) << "\n";
  log_stream << "\t\tswap:           " << fix(times[kTimeSolveSolve_swap], 8, 4)
             << "\n";
  log_stream << "\tResidual:               "
             << fix(times[kTimeSolveResidual], 8, 4) << " ("
             << fix(times[kTimeSolveResidual] / times[kTimeSolve] * 100, 4, 1)
             << "%)\n";
  log_stream << "\tOmega:                  "
             << fix(times[kTimeSolveOmega], 8, 4) << " ("
             << fix(times[kTimeSolveOmega] / times[kTimeSolve] * 100, 4, 1)
             << "%)\n";
#endif
  log_stream << "----------------------------------------------------\n";

#if HIPO_TIMING_LEVEL >= 3

  double total_blas_time =
      times[kTimeBlas_copy] + times[kTimeBlas_axpy] + times[kTimeBlas_scal] +
      times[kTimeBlas_swap] + times[kTimeBlas_gemv] + times[kTimeBlas_trsv] +
      times[kTimeBlas_tpsv] + times[kTimeBlas_ger] + times[kTimeBlas_trsm] +
      times[kTimeBlas_syrk] + times[kTimeBlas_gemm];

  log_stream << "BLAS time               \t" << fix(total_blas_time, 8, 4)
             << '\n';
  log_stream << "\tcopy:           \t" << fix(times[kTimeBlas_copy], 8, 4)
             << " ("
             << fix(times[kTimeBlas_copy] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_copy - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\taxpy:           \t" << fix(times[kTimeBlas_axpy], 8, 4)
             << " ("
             << fix(times[kTimeBlas_axpy] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_axpy - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tscal:           \t" << fix(times[kTimeBlas_scal], 8, 4)
             << " ("
             << fix(times[kTimeBlas_scal] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_scal - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tswap:           \t" << fix(times[kTimeBlas_swap], 8, 4)
             << " ("
             << fix(times[kTimeBlas_swap] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_swap - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tgemv:           \t" << fix(times[kTimeBlas_gemv], 8, 4)
             << " ("
             << fix(times[kTimeBlas_gemv] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_gemv - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\ttrsv:           \t" << fix(times[kTimeBlas_trsv], 8, 4)
             << " ("
             << fix(times[kTimeBlas_trsv] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_trsv - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\ttpsv:           \t" << fix(times[kTimeBlas_tpsv], 8, 4)
             << " ("
             << fix(times[kTimeBlas_tpsv] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_tpsv - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tger:            \t" << fix(times[kTimeBlas_ger], 8, 4)
             << " ("
             << fix(times[kTimeBlas_ger] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_ger - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\ttrsm:           \t" << fix(times[kTimeBlas_trsm], 8, 4)
             << " ("
             << fix(times[kTimeBlas_trsm] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_trsm - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tsyrk:           \t" << fix(times[kTimeBlas_syrk], 8, 4)
             << " ("
             << fix(times[kTimeBlas_syrk] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_syrk - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "\tgemm:           \t" << fix(times[kTimeBlas_gemm], 8, 4)
             << " ("
             << fix(times[kTimeBlas_gemm] / times[kTimeSolve] * 100, 4, 1)
             << "%) in "
             << integer(blas_calls[kTimeBlas_gemm - kTimeBlasStart], 10)
             << " calls\n";
  log_stream << "----------------------------------------------------\n";
#endif

  if (log) log->print(log_stream);
#endif
}

void DataCollector::printIter(const Log* log) const {
#ifdef HIPO_COLLECT_EXPENSIVE_DATA

  std::stringstream log_stream;

  log_stream << "\niter |   min D     max D     min L     max L  |"
                "    reg   swap    2x2    ws | "
                " max_reg   max_ws    norm1     maxdiag|\n";

  for (Int i = 0; i < iter_data_record_.size(); ++i) {
    const IterData& iter = iter_data_record_[i];

    log_stream << integer(i, 3) << "  |";
    log_stream << sci(iter.minD, 9, 1) << " ";
    log_stream << sci(iter.maxD, 9, 1) << " ";
    log_stream << sci(iter.minL, 9, 1) << " ";
    log_stream << sci(iter.maxL, 9, 1) << " |";
    log_stream << integer(iter.n_reg_piv, 6) << " ";
    log_stream << integer(iter.n_swap, 6) << " ";
    log_stream << integer(iter.n_2x2, 6) << " ";
    log_stream << integer(iter.n_wrong_sign, 6) << " |";
    log_stream << sci(iter.max_reg, 9, 1) << " ";
    log_stream << sci(iter.max_wrong_sign, 9, 1) << " ";
    log_stream << sci(iter.M_norm1, 9, 1) << " ";
    log_stream << sci(iter.M_maxdiag, 9, 1) << "|\n";
  }
  if (log) log->print(log_stream);
#endif
}

}  // namespace hipo