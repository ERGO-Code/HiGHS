/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/DynamicHipoLoader.h
 * @brief Dynamic loader for the optional HiPO library
 */
#ifndef LP_DATA_DYNAMIC_HIPO_LOADER_H_
#define LP_DATA_DYNAMIC_HIPO_LOADER_H_

#include <string>
#include <vector>

#include "lp_data/HighsCallback.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsTimer.h"

// Forward declaration
class HighsLpSolverObject;

// ABI version for compatibility checking
// Increment this when the C API signature changes
constexpr int kHipoAbiVersion = 1;

// C-style function pointer types for dynamic loading
extern "C" {

// Get the ABI version of the loaded HiPO library
typedef int (*hipo_get_abi_version_t)();

// Get the version string of the loaded HiPO library
typedef const char* (*hipo_get_version_t)();

// Solve LP using HiPO - uses actual HiGHS types for type safety
// Note: ABI version check ensures struct layouts match between builds
typedef HighsStatus (*hipo_solve_lp_t)(
    const HighsOptions& options,
    HighsTimer& timer,
    const HighsLp& lp,
    HighsBasis& highs_basis,
    HighsSolution& highs_solution,
    HighsModelStatus& model_status,
    HighsInfo& highs_info,
    HighsCallback& callback
);

}  // extern "C"

/**
 * Dynamic loader for the optional HiPO library.
 * 
 * This class handles runtime loading of the HiPO shared library,
 * allowing HiGHS to optionally use HiPO when it's installed via
 * the highspy-hipo package, without requiring HiPO at compile time.
 */
class DynamicHipoLoader {
 public:
  /**
   * Get the singleton instance of the loader.
   */
  static DynamicHipoLoader& instance();

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
  HighsStatus solveLp(HighsLpSolverObject& solver_object);

  /**
   * Get the last error message if loading failed.
   */
  const std::string& getLastError() const { return last_error_; }

 private:
  DynamicHipoLoader();
  ~DynamicHipoLoader();

  // Prevent copying
  DynamicHipoLoader(const DynamicHipoLoader&) = delete;
  DynamicHipoLoader& operator=(const DynamicHipoLoader&) = delete;

  /**
   * Attempt to load the HiPO library from various locations.
   */
  bool tryLoad();

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
   * Get platform-specific library search paths.
   */
  std::vector<std::string> getSearchPaths() const;

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
  hipo_get_abi_version_t fn_get_abi_version_ = nullptr;
  hipo_get_version_t fn_get_version_ = nullptr;
  hipo_solve_lp_t fn_solve_lp_ = nullptr;
};

#endif  // LP_DATA_DYNAMIC_HIPO_LOADER_H_
