/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ipm/hipo/HipoCApi.cpp
 * @brief C-style API implementation for the HiPO library
 */

#include "ipm/hipo/HipoCApi.h"

#include "ipm/IpxWrapper.h"
#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsOptions.h"
#include "lp_data/HighsSolution.h"
#include "lp_data/HighsStatus.h"
#include "util/HighsTimer.h"

extern "C" {

HIPO_API int hipo_get_abi_version(void) { return HIPO_ABI_VERSION; }

HIPO_API const char* hipo_get_version(void) { return HIPO_VERSION; }

}  // extern "C"

// Note: We use extern "C" here to disable C++ name mangling, allowing
// GetProcAddress/dlsym to find the function by its simple name "hipo_solve_lp".
// While extern "C" is typically for C-compatible interfaces, using C++ types
// (references, classes) works here because:
// 1. Both highspy and highspy-hipo are built with the same C++ compiler/ABI
// 2. The ABI version check ensures struct layouts match between builds
// 3. We only need C linkage for symbol lookup, not C type compatibility
extern "C" HIPO_API HighsStatus hipo_solve_lp(
    const HighsOptions& options,
    HighsTimer& timer,
    const HighsLp& lp,
    HighsBasis& highs_basis,
    HighsSolution& highs_solution,
    HighsModelStatus& model_status,
    HighsInfo& highs_info,
    HighsCallback& callback) {
  // Call the actual HiPO solver implementation
  // Note: solveLpHipo is defined in IpxWrapper.cpp when HIPO is defined
  return solveLpHipo(options, timer, lp, highs_basis, highs_solution,
                     model_status, highs_info, callback);
}
