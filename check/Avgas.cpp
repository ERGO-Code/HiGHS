/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file simplex/Avgas.cpp
 * @brief Utilities for tests with AVGAS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Avgas.h"
#include <cassert>

void Avgas::rows(int &num_row, std::vector<double> &rowLower, std::vector<double> &rowUpper) {
  rowLower.resize(10);
  rowUpper.resize(10);
  num_row = 0;
  rowLower[num_row] = -1;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = -1;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = -1;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = -1;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = -2;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = -2;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = 0;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = 0;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = 0;
  rowUpper[num_row] = 1e31;
  num_row++;
  rowLower[num_row] = 0;
  rowUpper[num_row] = 1e31;
  num_row++;
  assert(num_row == 10);
}

void Avgas::col(int col, int &num_col, int &num_nz,
		std::vector<double> &colCost, std::vector<double> &colLower, std::vector<double> &colUpper, 
		std::vector<int> &Astart, std::vector<int> &Aindex, std::vector<double> &Avalue) {

  colCost.resize(num_col+1);
  colLower.resize(num_col+1);
  colUpper.resize(num_col+1);
  Astart.resize(num_col+1);

  Astart[num_col] = num_nz;
  int num_new_nz = 4;
  if (col == 0) {
    colCost[num_col] = 0;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 0; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 4; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 6; Avalue[num_nz] =  2; num_nz++;
    Aindex[num_nz] = 7; Avalue[num_nz] =  5; num_nz++;
  } else if (col == 1) {
    colCost[num_col] = -2;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 0; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 5; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 8; Avalue[num_nz] =  1; num_nz++;
    Aindex[num_nz] = 9; Avalue[num_nz] =  1; num_nz++;
  } else if (col == 2) {
    colCost[num_col] = -1;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 1; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 4; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 6; Avalue[num_nz] =  1; num_nz++;
    Aindex[num_nz] = 7; Avalue[num_nz] =  3; num_nz++;
  } else if (col == 3) {
    num_new_nz = 3;
    colCost[num_col] = -3;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 1; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 5; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 8; Avalue[num_nz] = -1; num_nz++;
  } else if (col == 4) {
    num_new_nz = 3;
    colCost[num_col] = -2;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 2; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 4; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 7; Avalue[num_nz] = -3; num_nz++;
  } else if (col == 5) {
    colCost[num_col] = -4;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 2; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 5; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 8; Avalue[num_nz] = -3; num_nz++;
    Aindex[num_nz] = 9; Avalue[num_nz] = -3; num_nz++;
  } else if (col == 6) {
    colCost[num_col] = -3;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 3; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 4; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 6; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 7; Avalue[num_nz] = -1; num_nz++;
  } else if (col == 7) {
    colCost[num_col] = -5;
    colLower[num_col] = 0;
    colUpper[num_col] = 1;
    Aindex.resize(num_nz+num_new_nz);
    Avalue.resize(num_nz+num_new_nz);
    Aindex[num_nz] = 3; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 5; Avalue[num_nz] = -1; num_nz++;
    Aindex[num_nz] = 8; Avalue[num_nz] = -5; num_nz++;
    Aindex[num_nz] = 9; Avalue[num_nz] = -2; num_nz++;
  } 
  num_col++;
}
