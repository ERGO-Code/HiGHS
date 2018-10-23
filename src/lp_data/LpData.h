#ifndef LP_DATA_H
#define LP_DATA_H

#include <vector>

enum class LpError {
  none,
  matrix_dimensions,
  matrix_indices,
  matrix_start,
  matrix_value,
  col_bounds,
  row_bounds,
  objective
}

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
    if (numCol <=0 || numRow <= 0) return LpError::matrix_dimensions;
     
    // Check vectors.
    if (colCost.size() != numCol) return LpError::objective;
    if (colLower.size() != numCol || colUpper.size() != numCol) return LpError::col_bounds;
    if (rowLower.size() != numRow || rowUpper.size() != numRow) return LpError::row_bounds;

    // Check matrix.
    if (Astart.size() != numCol + 1) return LpError::matrix_start;
    for (int i = 0; i < numCol; i++) 
      if (Astart[i] > Astart[i+1]) 
        return LpError::matrix_start;


    

  }


#endif