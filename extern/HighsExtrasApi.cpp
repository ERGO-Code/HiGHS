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
  bind_api(api->template as<amd_methods>());
  bind_api(api->template as<metis_methods>());
  bind_api(api->template as<blas_methods>());
  bind_api(api->template as<rcm_methods>());
#endif

  return true;
}

}  // extern "C"
