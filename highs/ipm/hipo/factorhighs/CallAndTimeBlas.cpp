#include "CallAndTimeBlas.h"

#include "DataCollector.h"
#include "DenseFact.h"
#include "FactorHiGHSSettings.h"
#include "Timing.h"
#include "ipm/hipo/auxiliary/Auxiliary.h"
#include "ipm/hipo/auxiliary/mycblas.h"

#ifdef HIPO_EXTRAS
#include "DynamicDepsLoader.h"
#endif

namespace hipo {

// macros to interface with CBlas
#define TRANS(x) (x) == 'N' ? CblasNoTrans : CblasTrans
#define UPLO(x) (x) == 'U' ? CblasUpper : CblasLower
#define DIAG(x) (x) == 'N' ? CblasNonUnit : CblasUnit
#define SIDE(x) (x) == 'L' ? CblasLeft : CblasRight

// level 1

void callAndTime_daxpy(Int n, double da, const double* dx, Int incx, double* dy,
                       Int incy, DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_daxpy(n, da, dx, incx, dy, incy);
#elif __APPLE__
  cblas_daxpy(n, da, dx, incx, dy, incy);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_daxpy_(n, da, dx, incx, dy, incy);
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_axpy);
}

void callAndTime_dcopy(Int n, const double* dx, Int incx, double* dy, Int incy,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dcopy(n, dx, incx, dy, incy);
#elif __APPLE__
  cblas_dcopy(n, dx, incx, dy, incy);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dcopy_(n, dx, incx, dy, incy);
  }
#endif
  HIPO_CLOCK_STOP(3, data, kTimeBlas_copy);
}

void callAndTime_dscal(Int n, const double da, double* dx, Int incx,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dscal(n, da, dx, incx);
#elif __APPLE__
  cblas_dscal(n, da, dx, incx);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dscal_(n, da, dx, incx);
  }
#endif
  HIPO_CLOCK_STOP(3, data, kTimeBlas_scal);
}

void callAndTime_dswap(Int n, double* dx, Int incx, double* dy, Int incy,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dswap(n, dx, incx, dy, incy);
#elif __APPLE__
  cblas_dswap(n, dx, incx, dy, incy);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dswap_(n, dx, incx, dy, incy);
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_swap);
}

// level 2

void callAndTime_dgemv(char trans, Int m, Int n, double alpha, const double* A,
                       Int lda, const double* x, Int incx, double beta,
                       double* y, Int incy, DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dgemv(CblasColMajor, TRANS(trans), m, n, alpha, A, lda, x, incx, beta,
      y, incy);
#elif __APPLE__
  cblas_dgemv(CblasColMajor, TRANS(trans), m, n, alpha, A, lda, x, incx, beta,
      y, incy);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dgemv_
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_gemv);
}

void callAndTime_dtpsv(char uplo, char trans, char diag, Int n,
                       const double* ap, double* x, Int incx,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dtpsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, ap, x,
     incx);
#elif __APPLE__
  cblas_dtpsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, ap, x,
     incx);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dtpsv_(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, ap, x,
     incx);
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_tpsv);
}

void callAndTime_dtrsv(char uplo, char trans, char diag, Int n, const double* A,
                       Int lda, double* x, Int incx, DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dtrsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, A, lda, x,
              incx);
#elif __APPLE__
  cblas_dtrsv(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, A, lda, x,
              incx);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dtrsv_(CblasColMajor, UPLO(uplo), TRANS(trans), DIAG(diag), n, A, lda, x,
              incx);
  }
#endif


  HIPO_CLOCK_STOP(3, data, kTimeBlas_trsv);
}

void callAndTime_dger(Int m, Int n, double alpha, const double* x, Int incx,
                      const double* y, Int incy, double* A, Int lda,
                      DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);
#elif __APPLE__
  cblas_dger(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dger_(CblasColMajor, m, n, alpha, x, incx, y, incy, A, lda);
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_ger);
}

// level 3

void callAndTime_dgemm(char transa, char transb, Int m, Int n, Int k,
                       double alpha, const double* A, Int lda, const double* B,
                       Int ldb, double beta, double* C, Int ldc,
                       DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dgemm(CblasColMajor, TRANS(transa), TRANS(transb), m, n, k, alpha, A,
              lda, B, ldb, beta, C, ldc);
#elif __APPLE__
  cblas_dgemm(CblasColMajor, TRANS(transa), TRANS(transb), m, n, k, alpha, A,
              lda, B, ldb, beta, C, ldc);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dgemm_(CblasColMajor, TRANS(transa), TRANS(transb), m, n, k, alpha, A,
              lda, B, ldb, beta, C, ldc);
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_gemm);
}

void callAndTime_dsyrk(char uplo, char trans, Int n, Int k, double alpha,
                       const double* A, Int lda, double beta, double* C,
                       Int ldc, DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dsyrk(CblasColMajor, UPLO(uplo), TRANS(trans), n, k, alpha, A, lda,
              beta, C, ldc);
#elif __APPLE__
  cblas_dsyrk(CblasColMajor, UPLO(uplo), TRANS(trans), n, k, alpha, A, lda,
              beta, C, ldc);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dsyrk_(CblasColMajor, UPLO(uplo), TRANS(trans), n, k, alpha, A, lda,
              beta, C, ldc);
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_syrk);
}

void callAndTime_dtrsm(char side, char uplo, char trans, char diag, Int m,
                       Int n, double alpha, const double* A, Int lda, double* B,
                       Int ldb, DataCollector& data) {
  HIPO_CLOCK_CREATE;

#ifndef HIPO_EXTRAS
  cblas_dtrsm(CblasColMajor, SIDE(side), UPLO(uplo), TRANS(trans), DIAG(diag),
              m, n, alpha, A, lda, B, ldb);
#elif __APPLE__
  cblas_dtrsm(CblasColMajor, SIDE(side), UPLO(uplo), TRANS(trans), DIAG(diag),
              m, n, alpha, A, lda, B, ldb);
#else
  DynamicDepsLoader& hipo_loader = DynamicDepsLoader::instance();
  if (hipo_loader.isAvailable()) {
    hipo_loader.fn_hipo_extras_dtrsm_(CblasColMajor, SIDE(side), UPLO(uplo), TRANS(trans), DIAG(diag),
              m, n, alpha, A, lda, B, ldb);
  }
#endif

  HIPO_CLOCK_STOP(3, data, kTimeBlas_trsm);
}

}  // namespace hipo
