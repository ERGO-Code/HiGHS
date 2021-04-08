/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsModelObjectUtil.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSMODELOBJECTUTILS_H_
#define LP_DATA_HIGHSMODELOBJECTUTILS_H_

#include <cassert>
#include <iostream>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsModelObject.h"

void report_row_vec_sol(HighsInt nrow, vector<double>& XrowLower,
                        vector<double>& XrowUpper, vector<double>& XrowPrimal,
                        vector<double>& XrowDual,
                        vector<HighsInt>& XrowStatus) {
  // Report the LP row data and solution passed to the method, where
  // XrowStatus is the SCIP-like basis status
  if (nrow <= 0) return;
  printf("Row    St      Primal       Lower       Upper        Dual\n");
  for (HighsInt row = 0; row < nrow; row++) {
    if (XrowStatus[row] == (HighsInt)HighsBasisStatus::BASIC)
      printf("%6" HIGHSINT_FORMAT " BC", row);
    else if (XrowStatus[row] == (HighsInt)HighsBasisStatus::ZERO)
      printf("%6" HIGHSINT_FORMAT " FR", row);
    else if (XrowStatus[row] == (HighsInt)HighsBasisStatus::LOWER) {
      if (XrowLower[row] == XrowUpper[row])
        printf("%6" HIGHSINT_FORMAT " FX", row);
      else
        printf("%6" HIGHSINT_FORMAT " LB", row);
    } else if (XrowStatus[row] == (HighsInt)HighsBasisStatus::UPPER)
      printf("%6" HIGHSINT_FORMAT " UB", row);
    else
      printf("%6" HIGHSINT_FORMAT " ??", row);
    printf(" %11g %11g %11g %11g\n", XrowPrimal[row], XrowLower[row],
           XrowUpper[row], XrowDual[row]);
  }
}

void report_row_matrix(HighsInt nrow, vector<HighsInt>& XARstart,
                       vector<HighsInt>& XARindex, vector<double>& XARvalue) {
  // Report the row-wise matrix passed to the method
  if (nrow <= 0) return;
  printf("Row    Index       Value\n");
  for (HighsInt row = 0; row < nrow; row++) {
    printf("%6" HIGHSINT_FORMAT " Start %8" HIGHSINT_FORMAT "\n", row,
           XARstart[row]);
    for (HighsInt el = XARstart[row]; el < XARstart[row + 1]; el++) {
      printf("      %6" HIGHSINT_FORMAT " %11g\n", XARindex[el], XARvalue[el]);
    }
  }
  printf("       Start %8" HIGHSINT_FORMAT "\n", XARstart[nrow]);
}

void report_col_vec_sol(HighsInt ncol, vector<double>& XcolCost,
                        vector<double>& XcolLower, vector<double>& XcolUpper,
                        vector<double>& XcolPrimal, vector<double>& XcolDual,
                        vector<HighsInt>& XcolStatus) {
  // Report the LP column data and solution passed to the method,
  // where XcolStatus is the SCIP-like basis status
  if (ncol <= 0) return;
  printf(
      "Col    St      Primal       Lower       Upper        Dual        "
      "Cost\n");
  for (HighsInt col = 0; col < ncol; col++) {
    if (XcolStatus[col] == (HighsInt)HighsBasisStatus::BASIC)
      printf("%6" HIGHSINT_FORMAT " BC", col);
    else if (XcolStatus[col] == (HighsInt)HighsBasisStatus::ZERO)
      printf("%6" HIGHSINT_FORMAT " FR", col);
    else if (XcolStatus[col] == (HighsInt)HighsBasisStatus::LOWER) {
      if (colLower[col] == XcolUpper[col])
        printf("%6" HIGHSINT_FORMAT " FX", col);
      else
        printf("%6" HIGHSINT_FORMAT " LB", col);
    } else if (XcolStatus[col] == (HighsInt)HighsBasisStatus::UPPER)
      printf("%6" HIGHSINT_FORMAT " UB", col);
    else
      printf("%6" HIGHSINT_FORMAT " ??", col);
    printf(" %11g %11g %11g %11g %11g\n", XcolPrimal[col], colLower[col],
           XcolUpper[col], XcolDual[col], XcolCost[col]);
  }
}

#endif
