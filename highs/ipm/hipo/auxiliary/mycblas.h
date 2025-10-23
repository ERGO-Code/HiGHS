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

using hipoint = hipo::Int;

// level 1

void cblas_daxpy(const hipoint n, const double alpha, const double* x,
                 const hipoint incx, double* y, const hipoint incy);
void cblas_dcopy(const hipoint n, const double* x, const hipoint incx,
                 double* y, const hipoint incy);
void cblas_dscal(const hipoint n, const double alpha, double* x,
                 const hipoint incx);
void cblas_dswap(const hipoint n, double* x, const hipoint incx, double* y,
                 const hipoint incy);

// level 2

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa, const hipoint M,
                 const hipoint n, const double alpha, const double* A,
                 const hipoint lda, const double* x, const hipoint incx,
                 const double beta, double* y, const hipoint incy);

void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const hipoint n, const double* ap, double* x,
                 const hipoint incx);

void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                 const hipoint n, const double* a, const hipoint lda, double* x,
                 const hipoint incx);

void cblas_dger(const enum CBLAS_ORDER order, const hipoint m, const hipoint n,
                const double alpha, const double* x, const hipoint incx,
                const double* y, const hipoint incy, double* A,
                const hipoint lda);

// level 3

void cblas_dgemm(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_TRANSPOSE transb, const hipoint m,
                 const hipoint n, const hipoint k, const double alpha,
                 const double* A, const hipoint lda, const double* B,
                 const hipoint ldb, const double beta, double* C,
                 const hipoint ldc);

void cblas_dsyrk(const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo,
                 const enum CBLAS_TRANSPOSE trans, const hipoint n,
                 const hipoint k, const double alpha, const double* a,
                 const hipoint lda, const double beta, double* C,
                 const hipoint ldc);

void cblas_dtrsm(const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
                 const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
                 const enum CBLAS_DIAG diag, const hipoint m, const hipoint n,
                 const double alpha, const double* a, const hipoint lda,
                 double* b, const hipoint ldb);

#ifdef __cplusplus
}
#endif

#endif