/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkk.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "simplex/HEkk.h"

#include <cassert>
#include <iostream>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
//#include "simplex/HSimplex.h"
//#include "simplex/HSimplexDebug.h"
//#include "simplex/HFactorDebug.h"
//#include "simplex/SimplexTimer.h"
//#include "util/HighsRandom.h"
//#include "util/HighsUtils.h"

using std::cout;
using std::endl;

HighsStatus HEkk::init() { return HighsStatus::OK; }
HighsStatus HEkk::solve() {
  HighsLogMessage(
      options.logfile, HighsMessageType::INFO,
      "HEkk::solve called for LP with %d columns, %d rows and %d entries",
      simplex_lp.numCol_, simplex_lp.numRow_, simplex_lp.Astart_[simplex_lp.numCol_]);

  return HighsStatus::Error;
}

// Private methods
