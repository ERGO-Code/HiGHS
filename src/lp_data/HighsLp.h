/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_LP_H_
#define LP_DATA_HIGHS_LP_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "HConfig.h"
#include "lp_data/HConst.h"         // For HiGHS strategy options
#include "lp_data/HStruct.h"        // For HighsBasis and HighsSolution
#include "simplex/SimplexConst.h"   // For simplex strategy options
#include "simplex/SimplexStruct.h"  // For SimplexBasis

enum class LpAction {
  DUALISE = 0,
  PERMUTE,
  SCALE,
  NEW_COSTS,
  NEW_BOUNDS,
  NEW_BASIS,
  NEW_COLS,
  NEW_ROWS,
  DEL_COLS,
  DEL_ROWS,
  DEL_ROWS_BASIS_OK,
  SCALED_COL,
  SCALED_ROW,
  BACKTRACKING
};

class HighsLp;

class HighsLp {
 public:
  // Model data
  int numCol_ = 0;
  int numRow_ = 0;

  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;

  ObjSense sense_ = ObjSense::MINIMIZE;
  double offset_ = 0;

  std::string model_name_ = "";
  std::string lp_name_ = "";

  std::vector<std::string> row_names_;
  std::vector<std::string> col_names_;

  std::vector<int> integrality_;

  bool equalButForNames(const HighsLp& lp);
  bool operator==(const HighsLp& lp);
};

// Make sure the sizes of solution and basis vectors are consistent
// with numRow_ and numCol_
bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis);
bool isSolutionRightSize(const HighsLp& lp, const HighsSolution& solution);
bool isBasisRightSize(const HighsLp& lp, const HighsBasis& basis);

void clearSolutionUtil(HighsSolution& solution);
void clearBasisUtil(HighsBasis& solution);
void clearLp(HighsLp& lp);

#endif
