/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsExtrasApi.cpp
 * @brief C-style API implementation for the HiGHS library dependencies
 */

#include "HighsExtrasApi.h"

extern "C" {

// core
HIGHS_EXTRAS_API const char* highs_extras_get_version(void) {
  return HIGHS_EXTRAS_VERSION;
}

HIGHS_EXTRAS_API const char* highs_extras_get_copyright(void) {
  return "External Dependencies: Copyright (c) 2026 under Apache 2.0 license "
         "terms";
}

// metis
HIGHS_EXTRAS_API int highs_extras_metis_set_default_options(idx_t* options) {
  return Highs_METIS_SetDefaultOptions(options);
}

HIGHS_EXTRAS_API int highs_extras_metis_nodend(idx_t* nvtxs, const idx_t* xadj,
                                               const idx_t* adjncy, idx_t* vwgt,
                                               idx_t* options, idx_t* perm,
                                               idx_t* iperm) {
  return Highs_METIS_NodeND(nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
}

// amd
HIGHS_EXTRAS_API void highs_extras_amd_defaults(double Control[]) {
  Highs_amd_defaults(Control);
}

HIGHS_EXTRAS_API int highs_extras_amd_order(amd_int n, const amd_int Ap[],
                                            const amd_int Ai[], amd_int P[],
                                            double Control[], double Info[]) {
  int status = Highs_amd_order(n, Ap, Ai, P, Control, Info);
  return status;
}

// rcm
// HIGHS_EXTRAS_API int highs_extras_genrcm(HighsInt node_num,
//                                                   HighsInt adj_num,
//                                                   const HighsInt adj_row[],
//                                                   const HighsInt adj[],
//                                                   HighsInt perm[]) {
//   HighsInt status = Highs_genrcm(node_num, adj_num, adj_row, adj, perm);
//   return status;
// }

// blas
HIGHS_EXTRAS_API void highs_extras_daxpy(const blasint n, const double alpha,
                                         const double* x, const blasint incx,
                                         double* y, const blasint incy) {
  cblas_daxpy(n, alpha, x, incx, y, incy);
}

HIGHS_EXTRAS_API void highs_extras_dcopy(const blasint n, const double* x,
                                         const blasint incx, double* y,
                                         const blasint incy) {
  cblas_dcopy(n, x, incx, y, incy);
}

HIGHS_EXTRAS_API void highs_extras_dscal(const blasint n, const double alpha,
                                         double* x, const blasint incx) {
  cblas_dscal(n, alpha, x, incx);
}

HIGHS_EXTRAS_API void highs_extras_dswap(const blasint n, double* x,
                                         const blasint incx, double* y,
                                         const blasint incy) {
  cblas_dswap(n, x, incx, y, incy);
}

HIGHS_EXTRAS_API void highs_extras_dgemv(const enum CBLAS_ORDER order,
                                         const enum CBLAS_TRANSPOSE transa,
                                         const blasint M, const blasint n,
                                         const double alpha, const double* A,
                                         const blasint lda, const double* x,
                                         const blasint incx, const double beta,
                                         double* y, const blasint incy) {
  cblas_dgemv(order, transa, M, n, alpha, A, lda, x, incx, beta, y, incy);
}

HIGHS_EXTRAS_API void highs_extras_dtpsv(const enum CBLAS_ORDER order,
                                         const enum CBLAS_UPLO uplo,
                                         const enum CBLAS_TRANSPOSE transa,
                                         const enum CBLAS_DIAG diag,
                                         const blasint n, const double* ap,
                                         double* x, const blasint incx) {
  cblas_dtpsv(order, uplo, transa, diag, n, ap, x, incx);
}

HIGHS_EXTRAS_API void highs_extras_dtrsv(const enum CBLAS_ORDER order,
                                         const enum CBLAS_UPLO uplo,
                                         const enum CBLAS_TRANSPOSE transa,
                                         const enum CBLAS_DIAG diag,
                                         const blasint n, const double* a,
                                         const blasint lda, double* x,
                                         const blasint incx) {
  cblas_dtrsv(order, uplo, transa, diag, n, a, lda, x, incx);
}

HIGHS_EXTRAS_API void highs_extras_dger(const enum CBLAS_ORDER order,
                                        const blasint m, const blasint n,
                                        const double alpha, const double* x,
                                        const blasint incx, const double* y,
                                        const blasint incy, double* A,
                                        const blasint lda) {
  cblas_dger(order, m, n, alpha, x, incx, y, incy, A, lda);
}

HIGHS_EXTRAS_API void highs_extras_dgemm(
    const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE transa,
    const enum CBLAS_TRANSPOSE transb, const blasint m, const blasint n,
    const blasint k, const double alpha, const double* A, const blasint lda,
    const double* B, const blasint ldb, const double beta, double* C,
    const blasint ldc) {
  cblas_dgemm(order, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C,
              ldc);
}

HIGHS_EXTRAS_API void highs_extras_dsyrk(const enum CBLAS_ORDER order,
                                         const enum CBLAS_UPLO uplo,
                                         const enum CBLAS_TRANSPOSE trans,
                                         const blasint n, const blasint k,
                                         const double alpha, const double* a,
                                         const blasint lda, const double beta,
                                         double* C, const blasint ldc) {
  cblas_dsyrk(order, uplo, trans, n, k, alpha, a, lda, beta, C, ldc);
}

HIGHS_EXTRAS_API void highs_extras_dtrsm(
    const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
    const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
    const enum CBLAS_DIAG diag, const blasint m, const blasint n,
    const double alpha, const double* a, const blasint lda, double* b,
    const blasint ldb) {
  cblas_dtrsm(order, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

HIGHS_EXTRAS_API void highs_extras_openblas_set_num_threads(int num_threads) {
#if defined(HIPO_USES_OPENBLAS)
  openblas_set_num_threads(num_threads);
#endif
}

HIGHS_EXTRAS_API const char* highs_extras_blas_library() {
#ifdef BLAS_LIBRARIES
  return BLAS_LIBRARIES;
#else
#ifdef HIPO_USES_OPENBLAS
  return "OpenBLAS";
#else
  return "unknown";
#endif
#endif
}

}  // extern "C"
