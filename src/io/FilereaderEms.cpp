/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file io/FilereaderEms.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include <fstream>
#include <iomanip>

#include "FilereaderEms.h"
#include "HConst.h"

FilereaderRetcode FilereaderEms::readModelFromFile(const char* filename,
                                                   HighsLp& model) {
  std::ifstream f;
  int i;

  f.open(filename, std::ios::in);
  int numCol, numRow, AcountX;

  // counts
  f >> numCol;
  f >> numRow;
  f >> AcountX;

  model.numCol_ = numCol;
  model.numRow_ = numRow;
  model.nnz_ = AcountX;

  // matrix
  model.Astart_.resize(numCol + 1);
  model.Aindex_.resize(AcountX);
  model.Avalue_.resize(AcountX);

  for (i = 0; i < numCol + 1; i++) f >> model.Astart_[i];

  for (i = 0; i < AcountX; i++) f >> model.Aindex_[i];

  for (i = 0; i < AcountX; i++) f >> model.Avalue_[i];

  // cost and bounds
  model.colCost_.reserve(numCol);
  model.colLower_.reserve(numCol);
  model.colUpper_.reserve(numCol);

  model.colCost_.assign(numCol, 0);
  model.colLower_.assign(numCol, -HIGHS_CONST_INF);
  model.colUpper_.assign(numCol, HIGHS_CONST_INF);

  for (i = 0; i < numCol; i++) {
    f >> model.colCost_[i];
  }

  for (i = 0; i < numCol; i++) {
    f >> model.colLower_[i];
  }

  for (i = 0; i < numCol; i++) {
    f >> model.colUpper_[i];
  }

  model.rowLower_.reserve(numRow);
  model.rowUpper_.reserve(numRow);
  model.rowLower_.assign(numRow, -HIGHS_CONST_INF);
  model.rowUpper_.assign(numRow, HIGHS_CONST_INF);

  for (i = 0; i < numRow; i++) {
    f >> model.rowLower_[i];
  }

  for (i = 0; i < numRow; i++) {
    f >> model.rowUpper_[i];
  }

  f.close();
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderEms::writeModelToFile(const char* filename,
                                                  HighsLp& model) {
  std::ofstream f;
  f.open(filename, std::ios::out);

  // counts
  f << model.numCol_ << std::endl;
  f << model.numRow_ << std::endl;
  f << model.nnz_ << std::endl;

  // matrix
  for (int i = 0; i < model.numCol_ + 1; i++) f << model.Astart_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.nnz_; i++) f << model.Aindex_[i] << " ";
  f << std::endl;

  f << std::setprecision(9);
  for (int i = 0; i < model.nnz_; i++) f << model.Avalue_[i] << " ";
  f << std::endl;

  // cost and bounds
  f << std::setprecision(9);
  for (int i = 0; i < model.numCol_; i++) f << model.colCost_[i] << " ";

  f << std::endl;

  for (int i = 0; i < model.numCol_; i++) f << model.colLower_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.numCol_; i++) f << model.colUpper_[i] << " ";
  f << std::endl;

  f << std::setprecision(9);
  for (int i = 0; i < model.numRow_; i++) f << model.rowLower_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.numRow_; i++) f << model.rowUpper_[i] << " ";
  f << std::endl;

  f.close();
  return FilereaderRetcode::OKAY;
}
