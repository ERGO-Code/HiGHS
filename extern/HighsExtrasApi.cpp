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

// template recursion to set function pointer to internal method
template <class Methods, std::size_t Index, std::size_t Count>
struct bind_methods {
  static void apply(feature_api<Methods>& api) {
    using desc_type = typename std::tuple_element<Index, Methods>::type;
    api.template method<Index>() = desc_type::direct();
    bind_methods<Methods, Index + 1, Count>::apply(api);
  }
};

// specialization to terminate recursion
template <class Methods, std::size_t Count>
struct bind_methods<Methods, Count, Count> {
  static void apply(feature_api<Methods>&) {}
};

// recursively set function pointers to the direct methods
template <class Methods>
void bind_api(feature_api<Methods>& api) {
  bind_methods<Methods, 0, std::tuple_size<Methods>::value>::apply(api);
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
  bind_api(api->template as<amd_methods>());
  bind_api(api->template as<metis_methods>());
  bind_api(api->template as<blas_methods>());
  bind_api(api->template as<rcm_methods>());
#endif

  return true;
}

}  // extern "C"
