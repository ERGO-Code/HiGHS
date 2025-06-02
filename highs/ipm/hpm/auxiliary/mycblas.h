#ifndef FACTORHIGHS_MY_CBLAS_H
#define FACTORHIGHS_MY_CBLAS_H

#include "IntConfig.h"

// Provide definition for cblas functions
// Based on Netlib implementation

enum CBLAS_ORDER { CblasRowMajor = 101, CblasColMajor = 102 };
enum CBLAS_TRANSPOSE {
  CblasNoTrans = 111,
  CblasTrans = 112,
  CblasConjTrans = 113
};
enum CBLAS_UPLO { CblasUpper = 121, CblasLower = 122 };
enum CBLAS_DIAG { CblasNonUnit = 131, CblasUnit = 132 };
enum CBLAS_SIDE { CblasLeft = 141, CblasRight = 142 };

#ifdef __cplusplus
extern "C" {
#endif

using hpmint = highspm::Int;

// level 1

void cblas_daxpy(const hpmint n, const double alpha, const double* x,
                 const hpmint incx, double* y, const hpmint incy);
void cblas_dcopy(const hpmint n, const double* x, const hpmint incx, double* y,
                 const hpmint incy);
void cblas_dscal(const hpmint n, const double alpha, double* x,
                 const hpmint incx);
void cblas_dswap(const hpmint n, double* x, const hpmint incx, double* y,
                 const hpmint incy);

// level 2

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa, const hpmint M,
                 const hpmint n, const double alpha, const double* A,
                 const hpmint lda, const double* x, const hpmint incx,
                 const double beta, double* y, const hpmint incy);

void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const hpmint n, const double* ap, double* x,
                 const hpmint incx);

void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const hpmint n, const double* a, const hpmint lda, double* x,
                 const hpmint incx);

void cblas_dger(const enum CBLAS_ORDER order, const hpmint m, const hpmint n,
                const double alpha, const double* x, const hpmint incx,
                const double* y, const hpmint incy, double* A,
                const hpmint lda);

// level 3

void cblas_dgemm(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_TRANSPOSE transb, const hpmint m,
                 const hpmint n, const hpmint k, const double alpha,
                 const double* A, const hpmint lda, const double* B,
                 const hpmint ldb, const double beta, double* C,
                 const hpmint ldc);

void cblas_dsyrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const hpmint n,
                 const hpmint k, const double alpha, const double* a,
                 const hpmint lda, const double beta, double* C,
                 const hpmint ldc);

void cblas_dtrsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
                 const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_DIAG diag, const hpmint m, const hpmint n,
                 const double alpha, const double* a, const hpmint lda,
                 double* b, const hpmint ldb);

#ifdef __cplusplus
}
#endif

#endif