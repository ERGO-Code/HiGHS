/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/PresolveAnalysis.h
 * @brief
 */
#ifndef PRESOLVE_PRESOLVE_ANALYSIS_H_
#define PRESOLVE_PRESOLVE_ANALYSIS_H_

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsTimer.h"

namespace presolve {

using std::min;

constexpr double inf = std::numeric_limits<double>::infinity();

class PresolveTimer {
 public:
  PresolveTimer(HighsTimer& timer) : timer_(timer) {}

  HighsTimer& timer_;

  double time_limit = 0.0;

};

}  // namespace presolve

#endif
