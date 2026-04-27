/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsExtrasCApi.h
 * @brief C-style API for the HiGHS external deps library for dynamic loading
 */
#ifndef HIGHS_EXTRAS_C_API_H_
#define HIGHS_EXTRAS_C_API_H_

#include "HConfig.h"

// Export macro for shared library
// Building vs importing vs static linking is determined by build macros:
#if defined(_WIN32) || defined(_WIN64)
#ifdef HIGHS_EXTRAS_LIBRARY_BUILD
#define HIGHS_EXTRAS_API __declspec(dllexport)
#elif defined(HIGHS_SHARED_EXTRAS_LIBRARY)
#define HIGHS_EXTRAS_API __declspec(dllimport)
#else 
#define HIGHS_EXTRAS_API
#endif
#else
#if defined(HIGHS_EXTRAS_LIBRARY_BUILD) || defined(HIGHS_SHARED_EXTRAS_LIBRARY)
#define HIGHS_EXTRAS_API __attribute__((visibility("default")))
#else
#define HIGHS_EXTRAS_API
#endif
#endif

// Helper macros for stringification
#define HIGHS_EXTRAS_X_STRINGIFY(x) #x
#define HIGHS_EXTRAS_X_TOSTRING(x) HIGHS_EXTRAS_X_STRINGIFY(x)

// Version string
#define HIGHS_EXTRAS_VERSION                                                    \
  HIGHS_EXTRAS_X_TOSTRING(HIGHS_VERSION_MAJOR)                                  \
  "." HIGHS_EXTRAS_X_TOSTRING(HIGHS_VERSION_MINOR) "." HIGHS_EXTRAS_X_TOSTRING( \
      HIGHS_VERSION_PATCH)

// C++ API with actual types (outside extern "C" block)
// These use C++ references and HiGHS types directly
#ifdef __cplusplus

#include "amd/amd.h"
#include "metis/metis.h"
#include "mycblas.h"

// Note: We use extern "C" here to disable C++ name mangling, allowing
// GetProcAddress/dlsym to find the function by its simple name.
// 
// While extern "C" is typically for C-compatible interfaces,
// using C++ types (references, classes) works here because:
// 1. Both highspy and highspy-extras are built with the same C++ compiler/ABI
// 2. The ABI version check ensures struct layouts match between builds
// 3. We only need C linkage for symbol lookup, not C type compatibility
extern "C" {

/**
* Get the version string of the HiGHS extras library, used to verify
* compatibility with the main HiGHS library.
*
* @return Version string (e.g., "1.12.0")
*/
HIGHS_EXTRAS_API const char* highs_extras_get_version(void);

// If the external dependencies have a variety copyright statements, we can
// return them all here, or the most restrictive one.
HIGHS_EXTRAS_API const char* highs_extras_get_copyright(void);


// metis
HIGHS_EXTRAS_API int highs_extras_metis_set_default_options(idx_t* options);

HIGHS_EXTRAS_API int highs_extras_metis_nodend(idx_t* nvtxs, const idx_t* xadj,
                                               const idx_t* adjncy, idx_t* vwgt,
                                               idx_t* options, idx_t* perm,
                                               idx_t* iperm);


// amd
HIGHS_EXTRAS_API void highs_extras_amd_defaults(double Control[]);

HIGHS_EXTRAS_API int highs_extras_amd_order(amd_int n, const amd_int Ap[],
                                            const amd_int Ai[], amd_int P[],
                                            double Control[], double Info[]);

// blas
HIGHS_EXTRAS_API void highs_extras_daxpy(const blasint n, const double alpha,
                                         const double* x, const blasint incx,
                                         double* y, const blasint incy);

HIGHS_EXTRAS_API void highs_extras_dcopy(const blasint n, const double* x,
                                         const blasint incx, double* y,
                                         const blasint incy);

HIGHS_EXTRAS_API void highs_extras_dscal(const blasint n, const double alpha,
                                         double* x, const blasint incx);

HIGHS_EXTRAS_API void highs_extras_dswap(const blasint n, double* x,
                                         const blasint incx, double* y,
                                         const blasint incy);

HIGHS_EXTRAS_API void highs_extras_dgemv(const enum CBLAS_ORDER order,
                                         const enum CBLAS_TRANSPOSE transa,
                                         const blasint M, const blasint n,
                                         const double alpha, const double* A,
                                         const blasint lda, const double* x,
                                         const blasint incx, const double beta,
                                         double* y, const blasint incy);

HIGHS_EXTRAS_API void highs_extras_dtpsv(const enum CBLAS_ORDER order,
                                         const enum CBLAS_UPLO uplo,
                                         const enum CBLAS_TRANSPOSE transa,
                                         const enum CBLAS_DIAG diag,
                                         const blasint n, const double* ap,
                                         double* x, const blasint incx);

HIGHS_EXTRAS_API void highs_extras_dtrsv(const enum CBLAS_ORDER order,
                                         const enum CBLAS_UPLO uplo,
                                         const enum CBLAS_TRANSPOSE transa,
                                         const enum CBLAS_DIAG diag,
                                         const blasint n, const double* a,
                                         const blasint lda, double* x,
                                         const blasint incx);

HIGHS_EXTRAS_API void highs_extras_dger(const enum CBLAS_ORDER order,
                                        const blasint m, const blasint n,
                                        const double alpha, const double* x,
                                        const blasint incx, const double* y,
                                        const blasint incy, double* A,
                                        const blasint lda);

HIGHS_EXTRAS_API void highs_extras_dgemm(
    const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE transa,
    const enum CBLAS_TRANSPOSE transb, const blasint m, const blasint n,
    const blasint k, const double alpha, const double* A, const blasint lda,
    const double* B, const blasint ldb, const double beta, double* C,
    const blasint ldc);

HIGHS_EXTRAS_API void highs_extras_dsyrk(const enum CBLAS_ORDER order,
                                         const enum CBLAS_UPLO uplo,
                                         const enum CBLAS_TRANSPOSE trans,
                                         const blasint n, const blasint k,
                                         const double alpha, const double* a,
                                         const blasint lda, const double beta,
                                         double* C, const blasint ldc);

HIGHS_EXTRAS_API void highs_extras_dtrsm(
    const enum CBLAS_ORDER order, const enum CBLAS_SIDE side,
    const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE transa,
    const enum CBLAS_DIAG diag, const blasint m, const blasint n,
    const double alpha, const double* a, const blasint lda, double* b,
    const blasint ldb);

HIGHS_EXTRAS_API void highs_extras_openblas_set_num_threads(int num_threads);

HIGHS_EXTRAS_API const char* highs_extras_blas_library();

}  // extern "C"

namespace highs_extras_api {

struct core {
  using get_version_t = decltype(&highs_extras_get_version);
  using get_copyright_t = decltype(&highs_extras_get_copyright);
};

struct metis {
    using set_default_options_t =
            decltype(&highs_extras_metis_set_default_options);
    using nodend_t = decltype(&highs_extras_metis_nodend);
};

struct amd {
    using defaults_t = decltype(&highs_extras_amd_defaults);
    using order_t = decltype(&highs_extras_amd_order);
};

struct blas {
    using daxpy_t = decltype(&highs_extras_daxpy);
    using dcopy_t = decltype(&highs_extras_dcopy);
    using dscal_t = decltype(&highs_extras_dscal);
    using dswap_t = decltype(&highs_extras_dswap);
    using dgemv_t = decltype(&highs_extras_dgemv);
    using dtpsv_t = decltype(&highs_extras_dtpsv);
    using dtrsv_t = decltype(&highs_extras_dtrsv);
    using dger_t = decltype(&highs_extras_dger);
    using dgemm_t = decltype(&highs_extras_dgemm);
    using dsyrk_t = decltype(&highs_extras_dsyrk);
    using dtrsm_t = decltype(&highs_extras_dtrsm);
    using openblas_set_num_threads_t =
            decltype(&highs_extras_openblas_set_num_threads);
    using library_t = decltype(&highs_extras_blas_library);
};

}  // namespace highs_extras_api

#endif // __cplusplus

#endif  // HIGHS_EXTRAS_C_API_H_
