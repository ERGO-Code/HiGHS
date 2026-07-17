/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsAppExternalDeps.h
 * @brief Defines the set of external features available
 */

#ifndef HIGHS_APP_EXTERNAL_DEPS_H_
#define HIGHS_APP_EXTERNAL_DEPS_H_

#include "HighsExternalApi.h"
#include "cli11/CLI11.hpp"

namespace HighsExtras {
struct app_family {};

template <>
inline const HighsExtrasFeatureInfo* wrapper_storage<app_family>::getInfo() {
  static const HighsExtrasFeatureInfo info = {"CLI11", CLI11_VERSION,
                                              "BSD-3-Clause", true};
  return &info;
}

template <int Index>
struct app_feature : feature_base<app_family, Index> {
  static const char* name() {
    return std::get<Index>(std::make_tuple("cli11"));
  }
};

using cli11 = app_feature<0>;
using appAll = require<all, cli11>;

}  // namespace HighsExtras

#endif  // HIGHS_APP_EXTERNAL_DEPS_H_