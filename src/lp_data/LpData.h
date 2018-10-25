#ifndef LP_DATA_H
#define LP_DATA_H

#include <string>
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
};

LpError checkLp(const LpData& lp) const {
  // Check dimensions.
  if (lp.numCol <= 0 || lp.numRow <= 0) return LpError::matrix_dimensions;

  // Check vectors.
  if (lp.colCost.size() != lp.numCol) return LpError::objective;

  if (lp.colLower.size() != lp.numCol || lp.colUpper.size() != lp.numCol)
    return LpError::col_bounds;
  if (lp.rowLower.size() != lp.numRow || lp.rowUpper.size() != lp.numRow)
    return LpError::row_bounds;

  for (int i = 0; i < numRow; i++)
    if (lp.rowLower[i] < -HSOL_CONST_INF || lp.rowUpper[i] > HSOL_CONST_INF)
      return LpError::row_bounds;

  for (int j = 0; j < numCol; j++) {
    if (lp.colCost[j] < -HSOL_CONST_INF || lp.colCost[j] > HSOL_CONST_INF)
      return LpError::objective;

    if (lp.colLower[j] < -HSOL_CONST_INF || lp.colUpper[j] > HSOL_CONST_INF)
      return LpError::col_bounds;
    if (lp.colLower[j] > lp.colUpper[j] + kBoundTolerance)
      return LpError::col_bounds;
  }

  // Check matrix.
  const int nnz = lp.Avalue.size();
  if (nnz <= 0) return LpError::matrix_value;
  if (lp.Aindex.size() != nnz) return LpError::matrix_indices;

  if (lp.Astart.size() != numCol + 1) return LpError::matrix_start;
  for (int i = 0; i < numCol; i++) {
    if (lp.Astart[i] > lp.Astart[i + 1] || lp.Astart[i] >= nnz ||
        lp.Astart[i] < 0)
      return LpError::matrix_start;
  }

  for (int k = 0; k < nnz; k++) {
    if (lp.Aindex[k] < 0 || lp.Aindex[k] >= lp.numRow)
      return LpError::matrix_indices;
    if (lp.Avalue[k] < -HSOL_CONST_INF || lp.Avalue[k] > HSOL_CONST_INF)
      return LpError::matrix_value;
  }

  return LpError::none;
}

#endif