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

#include "io/FilereaderEms.h"
#include "lp_data/HConst.h"
#include "util/stringutil.h"

FilereaderRetcode FilereaderEms::readModelFromFile(const HighsOptions& options,
                                                   HighsLp &model) {
  std::ifstream f;
  int i;
  
  const char* filename = options.filename.c_str();
  f.open(filename, std::ios::in);
  if (f.is_open()) {
 
    std::string line;
    int numCol, numRow, AcountX;
    bool indices_from_one = false;

    // counts
    std::getline(f, line);
    if (trim(line) != "n_rows") {
      while (trim(line) != "n_rows")
        std::getline(f, line);
      indices_from_one = true;
    }
    f >> numRow;

    std::getline(f, line);
    while (trim(line) == "")
      std::getline(f, line);
    if (trim(line) != "n_columns")
      return FilereaderRetcode::PARSERERROR;
    f >> numCol;

    std::getline(f, line);
    while (trim(line) == "")
      std::getline(f, line);
    if (trim(line) != "n_matrix_elements")
      return FilereaderRetcode::PARSERERROR;
    f >> AcountX;

    model.numCol_ = numCol;
    model.numRow_ = numRow;
    model.nnz_ = AcountX;

    // matrix
    std::getline(f, line);
    while (trim(line) == "")
      std::getline(f, line);
    if (trim(line) != "matrix")
      return FilereaderRetcode::PARSERERROR;

    model.Astart_.resize(numCol + 1);
    model.Aindex_.resize(AcountX);
    model.Avalue_.resize(AcountX);

    for (i = 0; i < numCol + 1; i++) {
      f >> model.Astart_[i];
      if (indices_from_one)
        model.Astart_[i]--;
    }

    for (i = 0; i < AcountX; i++) {
      f >> model.Aindex_[i];
      if (indices_from_one)
        model.Aindex_[i]--;
    }

    for (i = 0; i < AcountX; i++)
      f >> model.Avalue_[i];

    // cost and bounds
    std::getline(f, line);
    while (trim(line) == "")
      std::getline(f, line);
    if (trim(line) != "column_bounds")
      return FilereaderRetcode::PARSERERROR;
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
    while (trim(line) == "")
      std::getline(f, line);
    if (trim(line) != "row_bounds")
      return FilereaderRetcode::PARSERERROR;
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
    while (trim(line) == "")
      std::getline(f, line);
    if (trim(line) != "column_costs")
      return FilereaderRetcode::PARSERERROR;
    model.colCost_.reserve(numCol);
    model.colCost_.assign(numCol, 0);
    for (i = 0; i < numCol; i++) {
      f >> model.colCost_[i];
    }

    std::getline(f, line);
    while (trim(line) == "" && f)
      std::getline(f, line);
    if (trim(line) == "integer_columns") {
      f >> model.numInt_;
      if (model.numInt_) {
	model.integrality_.resize(model.numCol_, 0);
	int iCol;
	for (i = 0; i < model.numInt_; i++) {
	  f >> iCol;
	  model.integrality_[iCol] = 1;
	}
      }
    }
    // Get the next keyword - should be "names", if anything
    std::getline(f, line);
    while (trim(line) == "" && f)
      std::getline(f, line);

    if (f && trim(line) != "names")
      return FilereaderRetcode::PARSERERROR;
    if (line == "names") {
      // Ignore length since we support any length.
      std::getline(f, line);
      if (trim(line) != "columns")
        std::getline(f, line);
      if (trim(line) != "columns")
        return FilereaderRetcode::PARSERERROR;

      model.row_names_.resize(numRow);
      model.col_names_.resize(numCol);

      for (i = 0; i < numCol; i++) {
        std::getline(f, line);
        model.col_names_[i] = trim(line);
      }

      std::getline(f, line);
      if (trim(line) != "rows")
        return FilereaderRetcode::PARSERERROR;

      for (i = 0; i < numRow; i++) {
        std::getline(f, line);
        model.row_names_[i] = trim(line);
      }
    }

    // todo:
    // while (trim(line) != "integer_variables" && trim(line) != "names") std::getline(f,
    // line);
    // ...

    f.close();
  } else {
       return FilereaderRetcode::FILENOTFOUND;
  }
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderEms::writeModelToFile(const char *filename,
                                                  HighsLp &model) {
  std::ofstream f;
  f.open(filename, std::ios::out);

  // counts
  f << "n_rows" << std::endl;
  f << model.numRow_ << std::endl;
  f << "n_columns" << std::endl;
  f << model.numCol_ << std::endl;
  f << "n_matrix_elements" << std::endl;
  f << model.nnz_ << std::endl;

  // matrix
  f << "matrix" << std::endl;
  for (int i = 0; i < model.numCol_ + 1; i++)
    f << model.Astart_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.nnz_; i++)
    f << model.Aindex_[i] << " ";
  f << std::endl;

  f << std::setprecision(9);
  for (int i = 0; i < model.nnz_; i++)
    f << model.Avalue_[i] << " ";
  f << std::endl;

  // cost and bounds
  f << std::setprecision(9);

  f << "column_bounds" << std::endl;
  for (int i = 0; i < model.numCol_; i++)
    f << model.colLower_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.numCol_; i++)
    f << model.colUpper_[i] << " ";
  f << std::endl;

  f << "row_bounds" << std::endl;
  f << std::setprecision(9);
  for (int i = 0; i < model.numRow_; i++)
    f << model.rowLower_[i] << " ";
  f << std::endl;

  for (int i = 0; i < model.numRow_; i++)
    f << model.rowUpper_[i] << " ";
  f << std::endl;

  f << "column_costs" << std::endl;
  for (int i = 0; i < model.numCol_; i++)
    f << model.colCost_[i] << " ";
  f << std::endl;

  if (model.row_names_.size() > 0 && model.col_names_.size() > 0) {
    f << "names" << std::endl;

    f << "columns" << std::endl;
    for (int i = 0; i < model.col_names_.size(); i++)
      f << model.col_names_[i] << std::endl;

    f << "rows" << std::endl;
    for (int i = 0; i < model.row_names_.size(); i++)
      f << model.row_names_[i] << std::endl;
  }

  // todo: integer variables.

  if (model.offset_ != 0)
    f << "shift" << std::endl << model.offset_ << std::endl;

  f << std::endl;
  f.close();
  return FilereaderRetcode::OKAY;
}

FilereaderRetcode FilereaderEms::readModelFromFile(const char *filename,
                                                   HighsModelBuilder &model) {
  return FilereaderRetcode::NOT_IMPLEMENTED;
}
