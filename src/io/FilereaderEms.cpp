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

std::string& ltrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}
 
std::string& rtrim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}
 
std::string& trim(std::string& str, const std::string& chars = "\t\n\v\f\r ")
{
    return ltrim(rtrim(str, chars), chars);
}

FilereaderRetcode FilereaderEms::readModelFromFile(const char* filename,
                                                   HighsLp& model) {
  std::ifstream f;
  int i;

  f.open(filename, std::ios::in);
  std::string line;
  int numCol, numRow, AcountX;

  // counts
  while (line != "n_rows") std::getline(f, line);
  f >> numRow;

  std::getline(f, line);
  while (trim(line) == "") std::getline(f, line);
  if (line != "n_columns") return FilereaderRetcode::PARSERERROR;
  f >> numCol;

  std::getline(f, line);
  while (trim(line) == "") std::getline(f, line);
  if (line != "n_matrix_elements") return FilereaderRetcode::PARSERERROR;
  f >> AcountX;

  model.numCol_ = numCol;
  model.numRow_ = numRow;
  model.nnz_ = AcountX;

  // matrix
  std::getline(f, line);
  while (trim(line) == "") std::getline(f, line);
  if (line != "matrix") return FilereaderRetcode::PARSERERROR;

  model.Astart_.resize(numCol + 1);
  model.Aindex_.resize(AcountX);
  model.Avalue_.resize(AcountX);

  for (i = 0; i < numCol + 1; i++) f >> model.Astart_[i];

  for (i = 0; i < AcountX; i++) f >> model.Aindex_[i];

  for (i = 0; i < AcountX; i++) f >> model.Avalue_[i];

  // cost and bounds
  std::getline(f, line);
  while (trim(line) == "") std::getline(f, line);
  if (line != "column_bounds") return FilereaderRetcode::PARSERERROR;
  model.colLower_.reserve(numCol);
  model.colUpper_.reserve(numCol);

  model.colLower_.assign(numCol, -HIGHS_CONST_INF);
  model.colUpper_.assign(numCol, HIGHS_CONST_INF);

  for (i = 0; i < numCol; i++) {
    f >> model.colLower_[i];
  }

  for (i = 0; i < numCol; i++) {
    f >> model.colUpper_[i];
  }

  std::getline(f, line);
  while (trim(line) == "") std::getline(f, line);
  if (line != "row_bounds") return FilereaderRetcode::PARSERERROR;
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

  std::getline(f, line);
  while (trim(line) == "") std::getline(f, line);
  if (line != "column_costs") return FilereaderRetcode::PARSERERROR;
  model.colCost_.reserve(numCol);
  model.colCost_.assign(numCol, 0);
  for (i = 0; i < numCol; i++) {
    f >> model.colCost_[i];
  }

  // todo:
  // while (line != "integer_variables" && line != "names") std::getline(f,
  // line);
  // ...

  f.close();
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderEms::writeModelToFile(const char* filename,
                                                  HighsLp& model) {
  std::ofstream f;
  f.open(filename, std::ios::out);

  // counts
  f << "n_rows" << std::endl;
  f << model.numCol_ << std::endl;
  f << "n_columns" << std::endl;
  f << model.numRow_ << std::endl;
  f << "n_matrix_elements" << std::endl;
  f << model.nnz_ << std::endl;

  // matrix
  f << "matrix" << std::endl;
  for (int i = 0; i < model.numCol_ + 1; i++) f << model.Astart_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.nnz_; i++) f << model.Aindex_[i] << " ";
  f << std::endl;

  f << std::setprecision(9);
  for (int i = 0; i < model.nnz_; i++) f << model.Avalue_[i] << " ";
  f << std::endl;

  // cost and bounds
  f << std::setprecision(9);

  f << "column_bounds" << std::endl;
  for (int i = 0; i < model.numCol_; i++) f << model.colLower_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.numCol_; i++) f << model.colUpper_[i] << " ";
  f << std::endl;

  f << "row_bounds" << std::endl;
  f << std::setprecision(9);
  for (int i = 0; i < model.numRow_; i++) f << model.rowLower_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.numRow_; i++) f << model.rowUpper_[i] << " ";
  f << std::endl;

  f << "column_costs" << std::endl;
  for (int i = 0; i < model.numCol_; i++) f << model.colCost_[i] << " ";
  f << std::endl;

  // todo: names & integer variables.
  f.close();
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderEms::readModelFromFile(const char* filename,
                                                   HighsModel& model) {
  return FilereaderRetcode::NOT_IMPLEMENTED;
}