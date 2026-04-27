/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsExternalDeps.h
 * @brief Manages access (dynamic or static) to optional external dependencies
 */
#ifndef HIGHS_EXTERNAL_DEPS_H_
#define HIGHS_EXTERNAL_DEPS_H_
#include <string>

#include "HighsExtrasApi.h"
#include "rcm/rcm.h"
#include "util/HighsInt.h"

// support dynamic or static function calls
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
#define HIGHS_EXTERN_CALL(fn_member, fn_direct) fn_member
#else
#define HIGHS_EXTERN_CALL(fn_member, fn_direct) fn_direct
#endif

/**
 * External dependencies for HiGHS that can either be dynamically loaded
 * or linked at compile time, depending on the build configuration.
 *
 * This class handles runtime loading or static linking of external dependencies
 * specifically those that have a different license from HiGHS.
 */
struct HighsExternalDeps {
  struct amd {
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
    highs_extras_api::amd::defaults_t defaults_ = nullptr;
    highs_extras_api::amd::order_t order_ = nullptr;
#endif

    static inline void set_defaults(double Control[]) {
      HIGHS_EXTERN_CALL(amd_.defaults_, Highs_amd_defaults)(Control);
    }

    static inline int order(amd_int n, const amd_int Ap[], const amd_int Ai[],
                            amd_int P[], double Control[], double Info[]) {
      return HIGHS_EXTERN_CALL(amd_.order_, Highs_amd_order)(n, Ap, Ai, P,
                                                             Control, Info);
    }
  };

  struct blas {
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
    highs_extras_api::blas::daxpy_t daxpy_ = nullptr;
    highs_extras_api::blas::dcopy_t dcopy_ = nullptr;
    highs_extras_api::blas::dscal_t dscal_ = nullptr;
    highs_extras_api::blas::dswap_t dswap_ = nullptr;
    highs_extras_api::blas::dgemv_t dgemv_ = nullptr;
    highs_extras_api::blas::dtpsv_t dtpsv_ = nullptr;
    highs_extras_api::blas::dtrsv_t dtrsv_ = nullptr;
    highs_extras_api::blas::dger_t dger_ = nullptr;
    highs_extras_api::blas::dgemm_t dgemm_ = nullptr;
    highs_extras_api::blas::dsyrk_t dsyrk_ = nullptr;
    highs_extras_api::blas::dtrsm_t dtrsm_ = nullptr;
    highs_extras_api::blas::openblas_set_num_threads_t set_num_threads_ =
        nullptr;
    highs_extras_api::blas::library_t library_ = nullptr;
#endif

    static void daxpy(blasint n, double alpha, const double* x, blasint incx,
                      double* y, blasint incy) {
      HIGHS_EXTERN_CALL(blas_.daxpy_, cblas_daxpy)(n, alpha, x, incx, y, incy);
    }

    static inline void dcopy(blasint n, const double* x, blasint incx,
                             double* y, blasint incy) {
      HIGHS_EXTERN_CALL(blas_.dcopy_, cblas_dcopy)(n, x, incx, y, incy);
    }

    static inline void dscal(blasint n, double alpha, double* x, blasint incx) {
      HIGHS_EXTERN_CALL(blas_.dscal_, cblas_dscal)(n, alpha, x, incx);
    }

    static inline void dswap(blasint n, double* x, blasint incx, double* y,
                             blasint incy) {
      HIGHS_EXTERN_CALL(blas_.dswap_, cblas_dswap)(n, x, incx, y, incy);
    }

    static inline void dgemv(CBLAS_ORDER order, CBLAS_TRANSPOSE trans,
                             blasint m, blasint n, double alpha,
                             const double* A, blasint lda, const double* x,
                             blasint incx, double beta, double* y,
                             blasint incy) {
      HIGHS_EXTERN_CALL(blas_.dgemv_, cblas_dgemv)(order, trans, m, n, alpha, A,
                                                   lda, x, incx, beta, y, incy);
    }

    static inline void dtpsv(CBLAS_ORDER order, CBLAS_UPLO uplo,
                             CBLAS_TRANSPOSE transa, CBLAS_DIAG diag, blasint n,
                             const double* ap, double* x, blasint incx) {
      HIGHS_EXTERN_CALL(blas_.dtpsv_, cblas_dtpsv)(order, uplo, transa, diag, n,
                                                   ap, x, incx);
    }

    static inline void dtrsv(CBLAS_ORDER order, CBLAS_UPLO uplo,
                             CBLAS_TRANSPOSE transa, const enum CBLAS_DIAG diag,
                             blasint n, const double* a, blasint lda, double* x,
                             blasint incx) {
      HIGHS_EXTERN_CALL(blas_.dtrsv_, cblas_dtrsv)(order, uplo, transa, diag, n,
                                                   a, lda, x, incx);
    }

    static inline void dger(CBLAS_ORDER order, blasint m, blasint n,
                            double alpha, const double* x, blasint incx,
                            const double* y, blasint incy, double* A,
                            blasint lda) {
      HIGHS_EXTERN_CALL(blas_.dger_, cblas_dger)(order, m, n, alpha, x, incx, y,
                                                 incy, A, lda);
    }

    static inline void dgemm(CBLAS_ORDER order, CBLAS_TRANSPOSE transa,
                             CBLAS_TRANSPOSE transb, blasint m, blasint n,
                             blasint k, double alpha, const double* A,
                             blasint lda, const double* B, blasint ldb,
                             double beta, double* C, blasint ldc) {
      HIGHS_EXTERN_CALL(blas_.dgemm_, cblas_dgemm)(
          order, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    static inline void dsyrk(CBLAS_ORDER order, CBLAS_UPLO uplo,
                             CBLAS_TRANSPOSE trans, blasint n, blasint k,
                             double alpha, const double* A, blasint lda,
                             double beta, double* C, blasint ldc) {
      HIGHS_EXTERN_CALL(blas_.dsyrk_, cblas_dsyrk)(order, uplo, trans, n, k,
                                                   alpha, A, lda, beta, C, ldc);
    }

    static inline void dtrsm(CBLAS_ORDER order, CBLAS_SIDE side,
                             CBLAS_UPLO uplo, CBLAS_TRANSPOSE transa,
                             CBLAS_DIAG diag, blasint m, blasint n,
                             double alpha, const double* A, blasint lda,
                             double* B, blasint ldb) {
      HIGHS_EXTERN_CALL(blas_.dtrsm_, cblas_dtrsm)(
          order, side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);
    }

    static inline void openblas_set_num_threads(int num_threads) {
      HIGHS_EXTERN_CALL(blas_.set_num_threads_,
                        highs_extras_openblas_set_num_threads)(num_threads);
    }

    static inline std::string blas_library() {
      return HIGHS_EXTERN_CALL(blas_.library_, highs_extras_blas_library)();
    }
  };

  struct metis {
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
    highs_extras_api::metis::set_default_options_t set_default_options_ =
        nullptr;
    highs_extras_api::metis::nodend_t nodend_ = nullptr;
#endif

    static inline int set_default_options(idx_t* options) {
      return HIGHS_EXTERN_CALL(metis_.set_default_options_,
                               Highs_METIS_SetDefaultOptions)(options);
    }

    static inline int nodeND(idx_t* nvtxs, const idx_t* xadj,
                             const idx_t* adjncy, idx_t* vwgt, idx_t* options,
                             idx_t* perm, idx_t* iperm) {
      return HIGHS_EXTERN_CALL(metis_.nodend_, Highs_METIS_NodeND)(
          nvtxs, xadj, adjncy, vwgt, options, perm, iperm);
    }
  };

  // rcm has MIT license and is always statically linked
  // accessed via HighsExternalDeps for consistency
  struct rcm {
    static inline int genrcm(HighsInt node_num, HighsInt adj_num,
                             const HighsInt adj_row[], const HighsInt adj[],
                             HighsInt perm[]) {
      return Highs_genrcm(node_num, adj_num, adj_row, adj, perm);
    }
  };

  // static access to dependencies
  static amd amd_;
  static blas blas_;
  static metis metis_;
  static rcm rcm_;

  static HighsExternalDeps& instance();

  HighsExternalDeps() = default;
  ~HighsExternalDeps() { unload(); }

  // Prevent copying
  HighsExternalDeps(const HighsExternalDeps&) = delete;
  HighsExternalDeps& operator=(const HighsExternalDeps&) = delete;

  static bool tryLoad();

// Exclude dependencies
#ifndef HIPO
  static inline bool isAvailable() { return false; }
  static constexpr bool isAvailableAtCompile() { return false; }

  static inline bool tryLoad(const std::string& path) { return false; }
  static inline void unload() {}
  static inline std::string getCopyrightInfo() { return ""; }
  static inline const std::string getLoadStatus() {
    return "Extras: Unavailable";
  }

// Shared library support
#elif defined(HIGHS_SHARED_EXTRAS_LIBRARY)
  static inline bool isAvailable() { return tryLoad(); }
  static constexpr bool isAvailableAtCompile() { return false; }

  static bool tryLoad(const std::string& path);
  static void unload();
  static std::string getCopyrightInfo() { return instance().get_copyright_(); }
  static const std::string getLoadStatus() { return instance().status_; }

// Static library support
#else
  static inline bool isAvailable() { return true; }
  static constexpr bool isAvailableAtCompile() { return true; }

  static inline bool tryLoad(const std::string& path) { return true; }
  static inline void unload() {}
  static std::string getCopyrightInfo() { return highs_extras_get_copyright(); }

  static inline const std::string getLoadStatus() {
    return "Extras: Available at compile time";
  }

#endif

 private:
#ifdef HIGHS_SHARED_EXTRAS_LIBRARY
  void* lib_handle_ = nullptr;
  highs_extras_api::core::get_copyright_t get_copyright_ = nullptr;
  bool available_ = false;
  std::string status_;

  void clear();
#endif
};

#endif
