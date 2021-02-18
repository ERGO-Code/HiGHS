/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/HEkkInterface.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "lp_data/HighsLpUtils.h"
#include "simplex/HEkk.h"

void HEkk::appendColsToVectors(const int num_new_col,
                               const vector<double>& colCost,
                               const vector<double>& colLower,
                               const vector<double>& colUpper) {
  appendColsToLpVectors(simplex_lp_, num_new_col, colCost, colLower, colUpper);
}

void HEkk::appendRowsToVectors(const int num_new_row,
                               const vector<double>& rowLower,
                               const vector<double>& rowUpper) {
  appendRowsToLpVectors(simplex_lp_, num_new_row, rowLower, rowUpper);
}

void HEkk::appendColsToMatrix(const int num_new_col, const int num_new_nz,
                              const int* XAstart, const int* XAindex,
                              const double* XAvalue) {
  appendColsToLpMatrix(simplex_lp_, num_new_col, num_new_nz, XAstart, XAindex,
                       XAvalue);
}

void HEkk::appendRowsToMatrix(const int num_new_row, const int num_new_nz,
                              const int* XARstart, const int* XARindex,
                              const double* XARvalue) {
  appendRowsToLpMatrix(simplex_lp_, num_new_row, num_new_nz, XARstart, XARindex,
                       XARvalue);
}
