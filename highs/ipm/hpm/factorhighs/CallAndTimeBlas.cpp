#include "CallAndTimeBlas.h"

#include "DataCollector.h"
#include "DenseFact.h"
#include "Timing.h"
#include "ipm/hpm/auxiliary/Auxiliary.h"
#include "ipm/hpm/auxiliary/mycblas.h"
#include "FactorHiGHSSettings.h"

namespace highspm {

// macros to interface with CBlas
#define TRANS(x) (x) == 'N' ? CblasNoTrans : CblasTrans
#define UPLO(x) (x) == 'U' ? CblasUpper : CblasLower
#define DIAG(x) (x) == 'N' ? CblasNonUnit : CblasUnit
#define SIDE(x) (x) == 'L' ? CblasLeft : CblasRight

// level 1

void callAndTime_daxpy(Int n, double da, const double* dx, Int incx, double* dy,
                       Int incy) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_daxpy(n, da, dx, incx, dy, incy);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_axpy, clock.stop());
#endif
}

void callAndTime_dcopy(Int n, const double* dx, Int incx, double* dy,
                       Int incy) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dcopy(n, dx, incx, dy, incy);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_copy, clock.stop());
#endif
}

void callAndTime_dscal(Int n, const double da, double* dx, Int incx) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dscal(n, da, dx, incx);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_scal, clock.stop());
#endif
}

void callAndTime_dswap(Int n, double* dx, Int incx, double* dy, Int incy) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dswap(n, dx, incx, dy, incy);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_swap, clock.stop());
#endif
}

// level 2

void callAndTime_dgemv(char trans, Int m, Int n, double alpha, const double* A,
                       Int lda, const double* x, Int incx, double beta,
                       double* y, Int incy) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dgemv(CblasColMajor, TRANS(trans), m, n, alpha, A, lda, x, incx, beta,
              y, incy);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_gemv, clock.stop());
#endif
}

void callAndTime_dtpsv(char uplo, char trans, char diag, Int n,
                       const double* ap, double* x, Int incx) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dtpsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, ap, x,
              incx);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_tpsv, clock.stop());
#endif
}

void callAndTime_dtrsv(char uplo, char trans, char diag, Int n, const double* A,
                       Int lda, double* x, Int incx) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dtrsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, A, lda, x,
              incx);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_trsv, clock.stop());
#endif
}

void callAndTime_dger(Int m, Int n, double alpha, const double* x, Int incx,
                      const double* y, Int incy, double* A, Int lda) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_ger, clock.stop());
#endif
}

// level 3

void callAndTime_dgemm(char transa, char transb, Int m, Int n, Int k,
                       double alpha, const double* A, Int lda, const double* B,
                       Int ldb, double beta, double* C, Int ldc) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dgemm(CblasColMajor, TRANS(transa), TRANS(transb), m, n, k, alpha, A,
              lda, B, ldb, beta, C, ldc);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_gemm, clock.stop());
#endif
}

void callAndTime_dsyrk(char uplo, char trans, Int n, Int k, double alpha,
                       const double* A, Int lda, double beta, double* C,
                       Int ldc) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dsyrk(CblasColMajor, UPLO(uplo), TRANS(trans), n, k, alpha, A, lda,
              beta, C, ldc);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_syrk, clock.stop());
#endif
}

void callAndTime_dtrsm(char side, char uplo, char trans, char diag, Int m,
                       Int n, double alpha, const double* A, Int lda, double* B,
                       Int ldb) {
#if HPM_TIMING_LEVEL >= 3
  Clock clock;
#endif
  cblas_dtrsm(CblasColMajor, SIDE(side), UPLO(uplo), TRANS(trans), DIAG(diag),
              m, n, alpha, A, lda, B, ldb);
#if HPM_TIMING_LEVEL >= 3
  DataCollector::get()->sumTime(kTimeBlas_trsm, clock.stop());
#endif
}

}  // namespace highspm
