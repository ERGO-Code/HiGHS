#ifndef LP_DATA_H
#define LP_DATA_H

#include <vector>

class LpData {
  // Model data
  int numCol;
  int numRow;
  int numRowOriginal;
  int numColOriginal;
  int numTot;

  std::vector<int> Astart;
  std::vector<int> Aindex;
  std::vector<double> Avalue;
  std::vector<double> colCost;
  std::vector<double> colLower;
  std::vector<double> colUpper;
  std::vector<double> rowLower;
  std::vector<double> rowUpper;
};
#endif