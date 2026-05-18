/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsExternalDeps.cpp
 * @brief Defines the set of external features available
 */

#include "HighsExternalDeps.h"

namespace HighsExtras {

#ifdef ZLIB_FOUND
constexpr bool __zlib_enabled = true;
#else
constexpr bool __zlib_enabled = false;
#define ZLIB_VERSION "unknown"
#endif

#ifdef CUPDLP_GPU
constexpr bool __cuda_enabled = true;
#else
constexpr bool __cuda_enabled = false;
#endif

const HighsExtrasFeatureInfo highs_family_info_[] = {
    {"pdqsort", "git:b1ef26a", "Zlib", true},
    {"zstr", "1.0.6", "MIT", __zlib_enabled},
    {"ZLIB", ZLIB_VERSION, "Zlib", __zlib_enabled},
    {"NVIDIA Driver API", "runtime", "N/A (not redistributed)",
     __cuda_enabled}};

}  // namespace HighsExtras