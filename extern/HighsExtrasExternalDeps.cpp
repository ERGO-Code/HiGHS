/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsExtrasExternalDeps.cpp
 * @brief Defines the set of external features available
 */

#include "HighsExtrasExternalDeps.h"

#include <sstream>
#include <string>

namespace HighsExtras {

#ifdef HIPO
constexpr bool __hipo_enabled = true;
#else
constexpr bool __hipo_enabled = false;
#define HIGHS_BLAS_VENDOR "unknown"
#define HIGHS_BLAS_VERSION "unknown"
#define HIGHS_BLAS_LICENSE "unknown"
#endif

// CMake does not always know the version of a system OpenBLAS, and the one
// loaded at runtime may differ from the one HiGHS was configured with, so ask
// OpenBLAS itself
static const char* blasVersion() {
#if defined(HIPO_USES_OPENBLAS)
  static const std::string version = [] {
    // e.g. "OpenBLAS 0.3.26 DYNAMIC_ARCH NO_AFFINITY Haswell MAX_THREADS=64"
    std::istringstream config(openblas_get_config());
    std::string name, number;
    config >> name >> number;
    return number.empty() ? std::string(HIGHS_BLAS_VERSION) : number;
  }();
  return version.c_str();
#else
  return HIGHS_BLAS_VERSION;
#endif
}

const HighsExtrasFeatureInfo* extras_family::getInfo() {
  // feature details need to be compiled into library for runtime access
  static const HighsExtrasFeatureInfo extras_feature_info[] = {
      {"SuiteSparse AMD", "7.12.1+", "BSD-3-Clause", __hipo_enabled},
      {HIGHS_BLAS_VENDOR, blasVersion(), HIGHS_BLAS_LICENSE, __hipo_enabled},
      {"METIS-GKlib", "5.2.1+", "Apache-2.0", __hipo_enabled},
      {"SPARSEPAK", "unversioned+", "MIT", __hipo_enabled}};
  return extras_feature_info;
}

}  // namespace HighsExtras
