/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/hipo/HipoCApi.h
 * @brief C-style API for the HiPO library for dynamic loading
 */
#ifndef IPM_HIPO_C_API_H_
#define IPM_HIPO_C_API_H_

#include "HConfig.h"

// Export macro for shared library
#if defined(_WIN32) || defined(_WIN64)
#ifdef HIPO_LIBRARY_BUILD
#define HIPO_API __declspec(dllexport)
#else
#define HIPO_API __declspec(dllimport)
#endif
#else
#define HIPO_API __attribute__((visibility("default")))
#endif

// ABI version - increment when the C API signature changes
#define HIPO_ABI_VERSION 1

// Helper macros for stringification
#define HIPO_STRINGIFY(x) #x
#define HIPO_TOSTRING(x) HIPO_STRINGIFY(x)

// Version string derived from HConfig.h
#define HIPO_VERSION \
  HIPO_TOSTRING(HIGHS_VERSION_MAJOR) "." \
  HIPO_TOSTRING(HIGHS_VERSION_MINOR) "." \
  HIPO_TOSTRING(HIGHS_VERSION_PATCH)

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Get the ABI version of the HiPO library.
 * Used for compatibility checking at runtime.
 * 
 * @return ABI version number
 */
HIPO_API int hipo_get_abi_version(void);

/**
 * Get the version string of the HiPO library.
 * 
 * @return Version string (e.g., "1.12.0")
 */
HIPO_API const char* hipo_get_version(void);

#ifdef __cplusplus
}  // extern "C"
#endif

// C++ API with actual types (outside extern "C" block)
// These use C++ references and HiGHS types directly
#ifdef __cplusplus

#include "lp_data/HighsCallback.h"
#include "lp_data/HighsInfo.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsTimer.h"

/**
 * Solve an LP using the HiPO solver.
 * 
 * This is the main entry point for the dynamically loaded HiPO library.
 * Uses actual HiGHS types for type safety.
 * 
 * @param options Reference to HighsOptions
 * @param timer Reference to HighsTimer
 * @param lp Reference to HighsLp
 * @param highs_basis Reference to HighsBasis
 * @param highs_solution Reference to HighsSolution
 * @param model_status Reference to HighsModelStatus
 * @param highs_info Reference to HighsInfo
 * @param callback Reference to HighsCallback
 * @return HighsStatus indicating success or failure
 */
extern "C" HIPO_API HighsStatus hipo_solve_lp(
    const HighsOptions& options,
    HighsTimer& timer,
    const HighsLp& lp,
    HighsBasis& highs_basis,
    HighsSolution& highs_solution,
    HighsModelStatus& model_status,
    HighsInfo& highs_info,
    HighsCallback& callback
);

#endif  // __cplusplus

#endif  // IPM_HIPO_C_API_H_
