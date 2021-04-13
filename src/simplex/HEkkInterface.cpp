/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkInterface.cpp
 * @brief
 */

#include "lp_data/HighsLpUtils.h"
#include "simplex/HEkk.h"

void HEkk::appendColsToVectors(const HighsInt num_new_col,
                               const vector<double>& colCost,
                               const vector<double>& colLower,
                               const vector<double>& colUpper) {
  appendColsToLpVectors(simplex_lp_, num_new_col, colCost, colLower, colUpper);
}

void HEkk::appendRowsToVectors(const HighsInt num_new_row,
                               const vector<double>& rowLower,
                               const vector<double>& rowUpper) {
  appendRowsToLpVectors(simplex_lp_, num_new_row, rowLower, rowUpper);
}

void HEkk::appendColsToMatrix(const HighsInt num_new_col,
                              const HighsInt num_new_nz,
                              const HighsInt* XAstart, const HighsInt* XAindex,
                              const double* XAvalue) {
  appendColsToLpMatrix(simplex_lp_, num_new_col, num_new_nz, XAstart, XAindex,
                       XAvalue);
}

void HEkk::appendRowsToMatrix(const HighsInt num_new_row,
                              const HighsInt num_new_nz,
                              const HighsInt* XARstart,
                              const HighsInt* XARindex,
                              const double* XARvalue) {
  appendRowsToLpMatrix(simplex_lp_, num_new_row, num_new_nz, XARstart, XARindex,
                       XARvalue);
}
