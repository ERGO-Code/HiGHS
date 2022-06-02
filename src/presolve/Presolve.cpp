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
/**@file presolve/Presolve.cpp
 * @brief
 */
#include "presolve/Presolve.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <queue>
#include <sstream>

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "presolve/HighsLpPropagator.h"
#include "presolve/PresolveUtils.h"

namespace presolve {

using std::cout;
using std::endl;
using std::flush;
using std::get;
using std::ios;
using std::list;
using std::make_pair;
using std::max;
using std::min;
using std::ofstream;
using std::setprecision;
using std::setw;
using std::stringstream;

dev_kkt_check::State Presolve::initState(const bool intermediate) {
  // update row value
  rowValue.assign(numRowOriginal, 0);
  for (HighsInt i = 0; i < numRowOriginal; ++i) {
    if (flagRow[i])
      for (HighsInt k = ARstart.at(i); k < ARstart.at(i + 1); ++k) {
        const HighsInt col = ARindex[k];
        if (flagCol[col]) rowValue.at(i) += valuePrimal.at(col) * ARvalue.at(k);
      }
  }

  if (!intermediate)
    return dev_kkt_check::State(
        numCol, numRow, Astart, Aend, Aindex, Avalue, ARstart, ARindex, ARvalue,
        colCost, colLower, colUpper, rowLower, rowUpper, flagCol, flagRow,
        colValue, colDual, rowValue, rowDual, col_status, row_status);

  // if intermediate step use checker's row and col bounds and cost
  return chk2.initState(numColOriginal, numRowOriginal, Astart, Aend, Aindex,
                        Avalue, ARstart, ARindex, ARvalue, flagCol, flagRow,
                        valuePrimal, valueColDual, rowValue, valueRowDual,
                        col_status, row_status);
}

}  // namespace presolve
