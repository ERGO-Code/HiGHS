#ifndef LP_DATA_H
#define LP_DATA_H

#include <vector>

#include "HConst.h"

enum class LpError {
  none,
  matrix_dimensions,
  matrix_indices,
  matrix_start,
  matrix_value,
  col_bounds,
  row_bounds,
  objective
};

class LpData {
 public:
  // Model data
  int numCol;
  int numRow;

  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;
  std::vector<double> colCost;
  std::vector<double> colLower;
  std::vector<double> colUpper;
  std::vector<double> rowLower;
  std::vector<double> rowUpper;

  LpError checkLp();
};

LpError LpData::checkLp() {

  // Check dimensions.
  if (numCol <= 0 || numRow <= 0) return LpError::matrix_dimensions;

  // Check vectors.
  if (colCost.size() != numCol) return LpError::objective;

  if (colLower.size() != numCol || colUpper.size() != numCol)
    return LpError::col_bounds;
  if (rowLower.size() != numRow || rowUpper.size() != numRow)
    return LpError::row_bounds;

  for (int i = 0; i < numRow; i++)
    if (rowLower[i] < -HSOL_CONST_INF || rowUpper[i] > HSOL_CONST_INF)
      return LpError::row_bounds;

  for (int j = 0; j < numCol; j++) {
    if (colCost[j] < -HSOL_CONST_INF || colCost[j] > HSOL_CONST_INF)
      return LpError::objective;

    if (colLower[j] < -HSOL_CONST_INF || colUpper[j] > HSOL_CONST_INF)
      return LpError::col_bounds;
    if (colLower[j] > colUpper[j] + kBoundTolerance) return LpError::col_bounds;
  }

  // Check matrix.
  const int nnz = Avalue.size();
  if (nnz <= 0) return LpError::matrix_value;
  if (Aindex.size() != nnz) return LpError::matrix_indices;

  if (Astart.size() != numCol + 1) return LpError::matrix_start;
  for (int i = 0; i < numCol; i++) {
    if (Astart[i] > Astart[i + 1] || Astart[i] >= nnz || Astart[i] < 0)
      return LpError::matrix_start;
  }

  for (int k = 0; k < nnz; k++) {
    if (Aindex[k] < 0 || Aindex[k] >= numRow) return LpError::matrix_indices;
    if (Avalue[k] < -HSOL_CONST_INF || Avalue[k] > HSOL_CONST_INF)
      return LpError::matrix_value;
  }

  return LpError::none;
}

#endif