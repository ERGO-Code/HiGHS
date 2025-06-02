#include "DataCollector.h"

#include "FactorHiGHSSettings.h"
#include "ipm/hpm/auxiliary/HpmLog.h"

namespace highspm {

// instance of DataCollector
DataCollector* DataCollector::ptr_ = nullptr;

// Functions to manage DataCollector

DataCollector::DataCollector() {
  counter_data_.times.resize(kTimeSize);
  counter_data_.blas_calls.resize(kTimeBlasEnd - kTimeBlasStart + 1);
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

// Expensive functions

void DataCollector::sumTime(TimeItems i, double t) {
#if HPM_TIMING_LEVEL > 0
  // Keep track of times and blas calls.
  std::lock_guard<std::mutex> lock(times_mutex_);
  counter_data_.times[i] += t;
#if HPM_TIMING_LEVEL >= 3
  if (i >= kTimeBlasStart && i <= kTimeBlasEnd)
    ++counter_data_.blas_calls[i - kTimeBlasStart];
#endif
#endif
}
void DataCollector::setExtremeEntries(double minD, double maxD, double minoffD,
                                      double maxoffD) {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  // Store max and min entries of D and L.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  back().minD = std::min(back().minD, minD);
  back().maxD = std::max(back().maxD, maxD);
  back().minL = std::min(back().minL, minoffD);
  back().maxL = std::max(back().maxL, maxoffD);
#endif
}
void DataCollector::countRegPiv() {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  // Increase the number of dynamically regularised pivots.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_reg_piv;
#endif
}
void DataCollector::countSwap() {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_swap;
#endif
}
void DataCollector::count2x2() {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_2x2;
#endif
}
void DataCollector::setWrongSign(double p) {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  ++back().n_wrong_sign;
  back().max_wrong_sign = std::max(back().max_wrong_sign, std::abs(p));
#endif
}
void DataCollector::setMaxReg(double new_reg) {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  // Keep track of maximum regularisation used.
  std::lock_guard<std::mutex> lock(iter_data_mutex_);
  back().max_reg = std::max(back().max_reg, new_reg);
#endif
}
void DataCollector::setBackError(double nw, double cw, Int large_components) {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  back().nw_back_err = std::max(back().nw_back_err, nw);
  back().cw_back_err = std::max(back().cw_back_err, cw);
  back().large_components_cw =
      std::max(back().large_components_cw, large_components);
#endif
}

// Cheap functions

void DataCollector::countSolves() { ++back().num_solves; }
void DataCollector::setOmega(double omega) {
  back().omega = std::max(back().omega, omega);
}
void DataCollector::setNorms(double norm1, double maxdiag) {
  back().M_norm1 = norm1;
  back().M_maxdiag = maxdiag;
}
void DataCollector::setSigma(double sigma, bool affinescaling) {
  if (affinescaling)
    back().sigma_aff = sigma;
  else
    back().sigma = sigma;
}
void DataCollector::setCorrectors(Int correctors) {
  back().correctors = correctors;
}
void DataCollector::setExtremeTheta(const std::vector<double>& scaling) {
  back().min_theta = std::numeric_limits<double>::infinity();
  back().max_theta = 0.0;
  for (double d : scaling) {
    if (d != 0.0) {
      back().min_theta = std::min(back().min_theta, 1.0 / d);
      back().max_theta = std::max(back().max_theta, 1.0 / d);
    }
  }
}
void DataCollector::setProducts(double min_prod, double max_prod, Int num_small,
                                Int num_large) {
  back().min_prod = min_prod;
  back().max_prod = max_prod;
  back().num_small_prod = num_small;
  back().num_large_prod = num_large;
}

void DataCollector::printTimes() const {
#if HPM_TIMING_LEVEL >= 1

  const std::vector<double>& times = counter_data_.times;

  Log::printf("----------------------------------------------------\n");
  Log::printf("Analyse time            \t%8.4f\n", times[kTimeAnalyse]);

#if HPM_TIMING_LEVEL >= 2
  Log::printf("\tMetis:                  %8.4f (%4.1f%%)\n",
              times[kTimeAnalyseMetis],
              times[kTimeAnalyseMetis] / times[kTimeAnalyse] * 100);
  Log::printf("\tTree:                   %8.4f (%4.1f%%)\n",
              times[kTimeAnalyseTree],
              times[kTimeAnalyseTree] / times[kTimeAnalyse] * 100);
  Log::printf("\tCounts:                 %8.4f (%4.1f%%)\n",
              times[kTimeAnalyseCount],
              times[kTimeAnalyseCount] / times[kTimeAnalyse] * 100);
  Log::printf("\tSupernodes:             %8.4f (%4.1f%%)\n",
              times[kTimeAnalyseSn],
              times[kTimeAnalyseSn] / times[kTimeAnalyse] * 100);
  Log::printf("\tReorder:                %8.4f (%4.1f%%)\n",
              times[kTimeAnalyseReorder],
              times[kTimeAnalyseReorder] / times[kTimeAnalyse] * 100);
  Log::printf("\tSn sparsity pattern:    %8.4f (%4.1f%%)\n",
              times[kTimeAnalysePattern],
              times[kTimeAnalysePattern] / times[kTimeAnalyse] * 100);
  Log::printf("\tRelative indices:       %8.4f (%4.1f%%)\n",
              times[kTimeAnalyseRelInd],
              times[kTimeAnalyseRelInd] / times[kTimeAnalyse] * 100);
#endif

  Log::printf("----------------------------------------------------\n");
  Log::printf("Factorise time          \t%8.4f\n", times[kTimeFactorise]);

#if HPM_TIMING_LEVEL >= 2
  Log::printf("\tPrepare:                %8.4f (%4.1f%%)\n",
              times[kTimeFactorisePrepare],
              times[kTimeFactorisePrepare] / times[kTimeFactorise] * 100);
  Log::printf(
      "\tAssembly original:      %8.4f (%4.1f%%)\n",
      times[kTimeFactoriseAssembleOriginal],
      times[kTimeFactoriseAssembleOriginal] / times[kTimeFactorise] * 100);
  Log::printf("\tAssemble children in F: %8.4f (%4.1f%%)\n",
              times[kTimeFactoriseAssembleChildrenFrontal],
              times[kTimeFactoriseAssembleChildrenFrontal] /
                  times[kTimeFactorise] * 100);
  Log::printf("\tAssemble children in C: %8.4f (%4.1f%%)\n",
              times[kTimeFactoriseAssembleChildrenClique],
              times[kTimeFactoriseAssembleChildrenClique] /
                  times[kTimeFactorise] * 100);
  Log::printf("\tDense factorisation:    %8.4f (%4.1f%%)\n",
              times[kTimeFactoriseDenseFact],
              times[kTimeFactoriseDenseFact] / times[kTimeFactorise] * 100);
  Log::printf("\t\tmain:           %8.4f\n", times[kTimeDenseFact_main]);
  Log::printf("\t\tSchur:          %8.4f\n", times[kTimeDenseFact_schur]);
  Log::printf("\t\tkernel:         %8.4f\n", times[kTimeDenseFact_kernel]);
  Log::printf("\t\tconvert:        %8.4f\n", times[kTimeDenseFact_convert]);
  Log::printf("\t\tpivoting:       %8.4f\n", times[kTimeDenseFact_pivoting]);
  Log::printf("\tTerminate:              %8.4f (%4.1f%%)\n",
              times[kTimeFactoriseTerminate],
              times[kTimeFactoriseTerminate] / times[kTimeFactorise] * 100);
#endif

  Log::printf("----------------------------------------------------\n");
  Log::printf("Solve time              \t%8.4f\n", times[kTimeSolve]);

#if HPM_TIMING_LEVEL >= 2
  Log::printf("\tPrepare:                %8.4f (%4.1f%%)\n",
              times[kTimeSolvePrepare],
              times[kTimeSolvePrepare] / times[kTimeSolve] * 100);
  Log::printf("\tSolve:                  %8.4f (%4.1f%%)\n",
              times[kTimeSolveSolve],
              times[kTimeSolveSolve] / times[kTimeSolve] * 100);
  Log::printf("\t\tdense:          %8.4f\n", times[kTimeSolveSolve_dense]);
  Log::printf("\t\tsparse:         %8.4f\n", times[kTimeSolveSolve_sparse]);
  Log::printf("\t\tswap:           %8.4f\n", times[kTimeSolveSolve_swap]);
  Log::printf("\tResidual:               %8.4f (%4.1f%%)\n",
              times[kTimeSolveResidual],
              times[kTimeSolveResidual] / times[kTimeSolve] * 100);
  Log::printf("\tOmega:                  %8.4f (%4.1f%%)\n",
              times[kTimeSolveOmega],
              times[kTimeSolveOmega] / times[kTimeSolve] * 100);
#endif
  Log::printf("----------------------------------------------------\n");

#if HPM_TIMING_LEVEL >= 3

  const std::vector<Int>& blas_calls = counter_data_.blas_calls;

  double total_blas_time =
      times[kTimeBlas_copy] + times[kTimeBlas_axpy] + times[kTimeBlas_scal] +
      times[kTimeBlas_swap] + times[kTimeBlas_gemv] + times[kTimeBlas_trsv] +
      times[kTimeBlas_tpsv] + times[kTimeBlas_ger] + times[kTimeBlas_trsm] +
      times[kTimeBlas_syrk] + times[kTimeBlas_gemm];

  Log::printf("BLAS time               \t%8.4f\n", total_blas_time);
  Log::printf("\tcopy:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_copy],
              times[kTimeBlas_copy] / total_blas_time * 100,
              blas_calls[kTimeBlas_copy - kTimeBlasStart]);
  Log::printf("\taxpy:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_axpy],
              times[kTimeBlas_axpy] / total_blas_time * 100,
              blas_calls[kTimeBlas_axpy - kTimeBlasStart]);
  Log::printf("\tscal:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_scal],
              times[kTimeBlas_scal] / total_blas_time * 100,
              blas_calls[kTimeBlas_scal - kTimeBlasStart]);
  Log::printf("\tswap:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_swap],
              times[kTimeBlas_swap] / total_blas_time * 100,
              blas_calls[kTimeBlas_swap - kTimeBlasStart]);
  Log::printf("\tgemv:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_gemv],
              times[kTimeBlas_gemv] / total_blas_time * 100,
              blas_calls[kTimeBlas_gemv - kTimeBlasStart]);
  Log::printf("\ttrsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_trsv],
              times[kTimeBlas_trsv] / total_blas_time * 100,
              blas_calls[kTimeBlas_trsv - kTimeBlasStart]);
  Log::printf("\ttpsv:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_tpsv],
              times[kTimeBlas_tpsv] / total_blas_time * 100,
              blas_calls[kTimeBlas_tpsv - kTimeBlasStart]);
  Log::printf("\tger:            \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_ger],
              times[kTimeBlas_ger] / total_blas_time * 100,
              blas_calls[kTimeBlas_ger - kTimeBlasStart]);
  Log::printf("\ttrsm:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_trsm],
              times[kTimeBlas_trsm] / total_blas_time * 100,
              blas_calls[kTimeBlas_trsm - kTimeBlasStart]);
  Log::printf("\tsyrk:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_syrk],
              times[kTimeBlas_syrk] / total_blas_time * 100,
              blas_calls[kTimeBlas_syrk - kTimeBlasStart]);
  Log::printf("\tgemm:           \t%8.4f (%4.1f%%) in %10d calls\n",
              times[kTimeBlas_gemm],
              times[kTimeBlas_gemm] / total_blas_time * 100,
              blas_calls[kTimeBlas_gemm - kTimeBlasStart]);
  Log::printf("----------------------------------------------------\n");
#endif
#endif
}

void DataCollector::printIter() const {
#ifdef HPM_COLLECT_EXPENSIVE_DATA
  Log::printf(
      "\niter |    min D     max D     min L     max L  |"
      "    reg   swap    2x2     ws | "
      "  max_reg  solv   omega     nw_be     cw_be     cw_large  max_ws |\n");
  for (Int i = 0; i < iter_data_record_.size(); ++i) {
    const IterData& iter = iter_data_record_[i];
    Log::printf(
        "%3d  |"
        " %9.1e %9.1e %9.1e %9.1e |"
        " %6d %6d %6d %6d |"
        " %9.1e %4d %9.1e %9.1e %9.1e %9d %9.1e |\n",
        i, iter.minD, iter.maxD, iter.minL, iter.maxL, iter.n_reg_piv,
        iter.n_swap, iter.n_2x2, iter.n_wrong_sign, iter.max_reg,
        iter.num_solves, iter.omega, iter.nw_back_err, iter.cw_back_err,
        iter.large_components_cw, iter.max_wrong_sign);
  }
#endif
}

}  // namespace highspm