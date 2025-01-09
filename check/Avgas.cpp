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
 */
#include "Avgas.h"

#include <cassert>
#include <cstdio>  // For printf

#include "lp_data/HConst.h"

const bool dev_run = false;

void Avgas::addRow(const HighsInt row, HighsInt& num_row, HighsInt& num_row_nz,
	      std::vector<double>& rowLower, std::vector<double>& rowUpper,
	      std::vector<HighsInt>& ARstart, std::vector<HighsInt>& ARindex,
	      std::vector<double>& ARvalue) {
  double lower;
  double upper;
  std::vector<HighsInt>index;
  std::vector<double>value;
  getRow(row, lower, upper, index, value);
  rowLower.push_back(lower);
  rowUpper.push_back(upper);
  HighsInt num_nz = index.size();
  ARstart.push_back(num_row_nz);
  assert(HighsInt(ARstart.size()) == num_row+1); 
  for (HighsInt iEl = 0; iEl < num_nz; iEl++) {
    ARvalue.push_back(value[iEl]);
    ARindex.push_back(index[iEl]);
  }
  num_row++;
  num_row_nz += num_nz;
}

void Avgas::getRow(const HighsInt row, 
		   double& lower, double& upper,
		   std::vector<HighsInt>& index,
		   std::vector<double>& value) {
  index.clear();
  value.clear();  
  upper = kHighsInf;
  if (row == 0) {
    lower = -1;
    index.push_back(0);
    value.push_back(-1);
    index.push_back(1);
    value.push_back(-1);
  } else if (row == 1) {
    lower = -1;
    index.push_back(2);
    value.push_back(-1);
    index.push_back(3);
    value.push_back(-1);
  } else if (row == 2) {
    lower = -1;
    index.push_back(4);
    value.push_back(-1);
    index.push_back(5);
    value.push_back(-1);
  } else if (row == 3) {
    lower = -1;
    index.push_back(6);
    value.push_back(-1);
    index.push_back(7);
    value.push_back(-1);
  } else if (row == 4) {
    lower = -2;
    index.push_back(0);
    value.push_back(-1);
    index.push_back(2);
    value.push_back(-1);
    index.push_back(4);
    value.push_back(-1);
    index.push_back(6);
    value.push_back(-1);
  } else if (row == 5) {
    lower = -2;
    index.push_back(1);
    value.push_back(-1);
    index.push_back(3);
    value.push_back(-1);
    index.push_back(5);
    value.push_back(-1);
    index.push_back(7);
    value.push_back(-1);
  } else if (row == 6) {
    lower = 0;
    index.push_back(0);
    value.push_back(2);
    index.push_back(2);
    value.push_back(1);
    index.push_back(6);
    value.push_back(-1);
  } else if (row == 7) {
    lower = 0;
    index.push_back(0);
    value.push_back(5);
    index.push_back(2);
    value.push_back(3);
    index.push_back(4);
    value.push_back(-3);
    index.push_back(6);
    value.push_back(-1);
  } else if (row == 8) {
    lower = 0;
    index.push_back(1);
    value.push_back(1);
    index.push_back(3);
    value.push_back(-1);
    index.push_back(5);
    value.push_back(-3);
    index.push_back(7);
    value.push_back(-5);
  } else if (row == 9) {
    lower = 0;
    index.push_back(1);
    value.push_back(1);
    index.push_back(5);
    value.push_back(-3);
    index.push_back(7);
    value.push_back(-2);
  } else {
    if (dev_run) printf("Avgas: row %d out of range\n", HighsInt(row));
  }
}

void Avgas::addCol(HighsInt col, HighsInt& num_col, HighsInt& num_col_nz,
                std::vector<double>& colCost, std::vector<double>& colLower,
                std::vector<double>& colUpper, std::vector<HighsInt>& Astart,
                std::vector<HighsInt>& Aindex, std::vector<double>& Avalue) {

  double cost;
  double lower;
  double upper;
  std::vector<HighsInt>index;
  std::vector<double>value;
  getCol(col, cost, lower, upper, index, value);
  colCost.push_back(cost);
  colLower.push_back(lower);
  colUpper.push_back(upper);
  HighsInt num_nz = index.size();
  Astart.push_back(num_col_nz);
  assert(HighsInt(Astart.size()) == num_col+1); 
  for (HighsInt iEl = 0; iEl < num_nz; iEl++) {
    Avalue.push_back(value[iEl]);
    Aindex.push_back(index[iEl]);
  }
  num_col++;
  num_col_nz += num_nz;
}

void Avgas::getCol(const HighsInt col, 
		   double& cost, double& lower, double& upper,
		   std::vector<HighsInt>& index, std::vector<double>& value) {

  lower = 0;
  upper = 1;
  if (col == 0) {
    cost = 0;
    index.push_back(0);
    value.push_back(-1);
    index.push_back(4);
    value.push_back(-1);
    index.push_back(6);
    value.push_back(2);
    index.push_back(7);
    value.push_back(5);
  } else if (col == 1) {
    cost = -2;
    index.push_back(0);
    value.push_back(-1);
    index.push_back(5);
    value.push_back(-1);
    index.push_back(8);
    value.push_back(1);
    index.push_back(9);
    value.push_back(1);
  } else if (col == 2) {
    cost = -1;
    index.push_back(1);
    value.push_back(-1);
    index.push_back(4);
    value.push_back(-1);
    index.push_back(6);
    value.push_back(1);
    index.push_back(7);
    value.push_back(3);
  } else if (col == 3) {
    cost = -3;
    index.push_back(1);
    value.push_back(-1);
    index.push_back(5);
    value.push_back(-1);
    index.push_back(8);
    value.push_back(-1);
  } else if (col == 4) {
    cost = -2;
    index.push_back(2);
    value.push_back(-1);
    index.push_back(4);
    value.push_back(-1);
    index.push_back(7);
    value.push_back(-3);
  } else if (col == 5) {
    cost = -4;
    index.push_back(2);
    value.push_back(-1);
    index.push_back(5);
    value.push_back(-1);
    index.push_back(8);
    value.push_back(-3);
    index.push_back(9);
    value.push_back(-3);
  } else if (col == 6) {
    cost = -3;
    index.push_back(3);
    value.push_back(-1);
    index.push_back(4);
    value.push_back(-1);
    index.push_back(6);
    value.push_back(-1);
    index.push_back(7);
    value.push_back(-1);
  } else if (col == 7) {
    cost = -5;
    index.push_back(3);
    value.push_back(-1);
    index.push_back(5);
    value.push_back(-1);
    index.push_back(8);
    value.push_back(-5);
    index.push_back(9);
    value.push_back(-2);
  } else {
    if (dev_run) printf("Avgas: col %d out of range\n", HighsInt(col));
  }
}

