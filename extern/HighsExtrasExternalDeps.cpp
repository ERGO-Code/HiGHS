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

namespace HighsExtras {

#ifdef HIPO
constexpr bool __hipo_enabled = true;
#else
constexpr bool __hipo_enabled = false;
#define HIGHS_BLAS_VENDOR "unknown"
#define HIGHS_BLAS_VERSION "unknown"
#define HIGHS_BLAS_LICENSE "unknown"
#endif

// feature details need to be compiled into library for runtime access
const HighsExtrasFeatureInfo extras_feature_info[] = {
    {"SuiteSparse AMD", "7.12.1+", "BSD-3-Clause", __hipo_enabled},
    {HIGHS_BLAS_VENDOR, HIGHS_BLAS_VERSION, HIGHS_BLAS_LICENSE, __hipo_enabled},
    {"METIS-GKlib", "5.2.1+", "Apache-2.0", __hipo_enabled},
    {"SPARSEPAK", "unversioned+", "MIT", __hipo_enabled}};

#ifdef HIPO

template <>
void bind_api(feature_api<amd_methods>& api) {
  bind_from_tuple(api, std::make_tuple(&Highs_amd_defaults, &Highs_amd_order));
}

template <>
void bind_api(feature_api<blas_methods>& api) {
  bind_from_tuple(api, std::make_tuple(&cblas_daxpy, &cblas_dcopy, &cblas_dscal,
                                       &cblas_dswap, &cblas_dgemv, &cblas_dtpsv,
                                       &cblas_dtrsv, &cblas_dger, &cblas_dgemm,
                                       &cblas_dsyrk, &cblas_dtrsm,
                                       &highs_openblas_set_num_threads));
}

template <>
void bind_api(feature_api<metis_methods>& api) {
  bind_from_tuple(api, std::make_tuple(&Highs_METIS_SetDefaultOptions,
                                       &Highs_METIS_NodeND));
}

template <>
void bind_api(feature_api<rcm_methods>& api) {
  bind_from_tuple(api, std::make_tuple(&Highs_genrcm));
}

#endif

}  // namespace HighsExtras
