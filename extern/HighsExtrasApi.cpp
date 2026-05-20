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

#include "HighsExtrasExternalDeps.h"

using namespace HighsExtras;

template <class Methods, class... Fn>
void bind_api(HighsExtrasApi* api, const std::tuple<Fn...>& funcs) {
  static_assert(std::tuple_size<Methods>::value == sizeof...(Fn),
                "bind_from_tuple requires one function per API entry");
  bind_from_tuple_impl<Methods, std::tuple<Fn...>, 0>::apply(
      api->template as<Methods>(), funcs);
}

extern "C" {

HIGHS_EXTRAS_API const char* HighsExtras_getVersion(void) {
  return HIGHS_EXTRAS_VERSION;
}

HIGHS_EXTRAS_API const HighsExtrasFeatureInfo* HighsExtras_getFeatureInfo() {
  return extras_feature_info;
}

HIGHS_EXTRAS_API bool HighsExtras_getApi(HighsExtrasApi* api) {
  if (!api) return false;

  *api = HighsExtrasApi{};

  // set function pointers for each feature API, if available
  // i.e., use direct descriptor to set the storage value
#ifdef HIPO

  bind_api<amd_methods>(api,
                        std::make_tuple(&Highs_amd_defaults, &Highs_amd_order));

  bind_api<blas_methods>(
      api,
      std::make_tuple(&cblas_daxpy, &cblas_dcopy, &cblas_dscal, &cblas_dswap,
                      &cblas_dgemv, &cblas_dtpsv, &cblas_dtrsv, &cblas_dger,
                      &cblas_dgemm, &cblas_dsyrk, &cblas_dtrsm,
                      &highs_openblas_set_num_threads));

  bind_api<metis_methods>(api, std::make_tuple(&Highs_METIS_SetDefaultOptions,
                                               &Highs_METIS_NodeND));

  bind_api<rcm_methods>(api, std::make_tuple(&Highs_genrcm));

#endif

  return true;
}

}  // extern "C"
