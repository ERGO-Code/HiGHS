/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/DynamicDepsLoader.h
 * @brief Dynamic loader for the optional HiPO library dependencies
 */
#ifndef DYNAMIC_DEPS_LOADER_H_
#define DYNAMIC_DEPS_LOADER_H_

#include <string>
#include <vector>

// Forward declaration
// class HighsLpSolverObject;

// #if IDXTYPEWIDTH == 32
// typedef int32_t idx_t;

#include "amd/amd.h"
#include "metis/metis.h"
#include "rcm/rcm.h"
#include "ipm/hipo/auxiliary/mycblas.h"

// ABI version for compatibility checking
// Increment this when the C API signature changes
constexpr int kHipoExtrasAbiVersion = 1;

// C-style function pointer types for dynamic loading
extern "C" {

// Get the ABI version of the loaded HiPO library
typedef int (*hipo_extras_get_abi_version_t)();

// Get the version string of the loaded HiPO library
typedef const char* (*hipo_extras_get_version_t)();

// typedef HighsStatus (*hipo_solve_lp_t)(
typedef int (*hipo_extras_metis_set_default_options_t)(idx_t* options);

typedef int (*hipo_extras_metis_nodend_t)(idx_t* nvtxs, const idx_t* xadj,
                                          const idx_t* adjncy, idx_t* vwgt,
                                          idx_t* options, idx_t* perm,
                                          idx_t* iperm);

typedef void (*hipo_extras_amd_defaults_t)(double Control[]);

typedef int (*hipo_extras_amd_order_t)(amd_int n, const amd_int Ap[],
                                       const amd_int Ai[], amd_int P[],
                                       double Control[], double Info[]);

typedef int (*hipo_extras_genrcm_t)(HighsInt node_num, HighsInt adj_num,
                                    const HighsInt adj_row[],
                                    const HighsInt adj[], HighsInt perm[]);

// Solve LP using HiPO - uses actual HiGHS types for type safety
// Note: ABI version check ensures struct layouts match between builds
// typedef HighsStatus (*hipo_solve_lp_t)(
//     const HighsOptions& options,
//     HighsTimer& timer,
//     const HighsLp& lp,
//     HighsBasis& highs_basis,
//     HighsSolution& highs_solution,
//     HighsModelStatus& model_status,
//     HighsInfo& highs_info,
//     HighsCallback& callback
// );

}  // extern "C"

/**
 * Dynamic loader for the optional HiPO library.
 *
 * This class handles runtime loading of the HiPO shared library,
 * allowing HiGHS to optionally use HiPO when it's installed via
 * the highspy-hipo package, without requiring HiPO at compile time.
 */
class DynamicDepsLoader {
 public:
  /**
   * Get the singleton instance of the loader.
   */
  static DynamicDepsLoader& instance();

  /**
   * Check if the HiPO library is available and compatible.
   * This performs lazy initialization on first call.
   *
   * @return true if HiPO is loaded and ABI-compatible
   */
  bool isAvailable();

  /**
   * Get the version string of the loaded HiPO library.
   *
   * @return Version string, or empty string if not loaded
   */
  std::string getVersion() const;

  /**
   * Solve LP using the dynamically loaded HiPO solver.
   *
   * @param solver_object The LP solver object containing problem data
   * @return HighsStatus indicating success or failure
   */
  // HighsStatus solveLp(HighsLpSolverObject& solver_object);
  // int hipo_extras_metis_set_default_options( idx_t *options);

  /**
   * Get the last error message if loading failed.
   */
  const std::string& getLastError() const { return last_error_; }

 public:
  DynamicDepsLoader();
  ~DynamicDepsLoader();

  // Prevent copying
  DynamicDepsLoader(const DynamicDepsLoader&) = delete;
  DynamicDepsLoader& operator=(const DynamicDepsLoader&) = delete;

  /**
   * Attempt to load the HiPO library from various locations.
   */
  bool tryLoad(const std::string path);

  /**
   * Load a library by path.
   * @return true if successful
   */
  bool loadLibrary(const std::string& path);

  /**
   * Resolve a symbol from the loaded library.
   */
  void* resolveSymbol(const char* name);

  /**
   * Resolve all required function pointers.
   * @return true if all functions were resolved
   */
  bool resolveFunctions();

  /**
   * Unload the library if loaded.
   */
  void unloadLibrary();

  /**
   * Get platform-specific library filename.
   */
  std::string getLibraryFilename() const;

  // Library handle (platform-specific)
#if defined(_WIN32) || defined(_WIN64)
  void* lib_handle_ = nullptr;  // HMODULE
#else
  void* lib_handle_ = nullptr;
#endif

  // State
  bool initialized_ = false;
  bool available_ = false;
  std::string last_error_;
  std::string version_;

  // Function pointers
  hipo_extras_get_abi_version_t fn_get_abi_version_ = nullptr;
  hipo_extras_get_version_t fn_get_version_ = nullptr;

  // # todo: add function pointers here
  // hipo_solve_lp_t fn_solve_lp_ = nullptr;
  hipo_extras_metis_set_default_options_t
      fn_hipo_extras_metis_set_default_options_ = nullptr;
  hipo_extras_metis_nodend_t fn_hipo_extras_metis_nodend_ = nullptr;

  hipo_extras_amd_defaults_t fn_hipo_extras_amd_defaults_ = nullptr;
  hipo_extras_amd_order_t fn_hipo_extras_amd_order_ = nullptr;

  hipo_extras_genrcm_t fn_hipo_extras_genrcm_ = nullptr;

  // hipo_extras_daxpy_t fn_hipo_extras_daxpy_=   nullptr;
  // hipo_extras_docopy_t fn_hipo_extras_docopy_=  nullptr;
  // hipo_extras_dscal_t fn_hipo_extras_dscal_=  nullptr;
  // hipo_extras_dswap_t fn_hipo_extras_dswap_=  nullptr;
  // hipo_extras_dtpsv_t fn_hipo_extras_dtpsv_=  nullptr;
  // hipo_extras_dtrsv_t fn_hipo_extras_dtrsv_=  nullptr;
  // hipo_extras_dger_t fn_hipo_extras_dger_= nullptr;
  // hipo_extras_dgemm_t  fn_hipo_extras_dgemm_ = nullptr;
  // hipo_extras_dsyrk_t fn_hipo_extras_dsyrk_=  nullptr;
  // hipo_extras_dtrsm_t fn_hipo_extras_dtrsm_=  nullptr;
};

#endif  // LP_DATA_DYNAMIC_HIPO_LOADER_H_
