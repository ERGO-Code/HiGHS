#ifndef HIGHS_MIP_H_
#define HIGHS_MIP_H_

#include <cstdint>
#include <string>
#include <vector>

#include "lp_data/HighsLp.h"

class HighsMip {
 public:
  // Model data
  int numCol_ = 0;
  int numRow_ = 0;

  std::vector<int> ARstart_;
  std::vector<int> ARindex_;
  std::vector<double> ARvalue_;
  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;
  mutable std::vector<double> maxAbsRowCoef_;
  std::vector<uint8_t> integrality_;

  std::vector<double> debugSolution_;

  ObjSense sense_ = ObjSense::MINIMIZE;
  double offset_ = 0;

  std::string model_name_ = "";
  std::string lp_name_ = "";

  std::vector<std::string> row_names_;
  std::vector<std::string> col_names_;

  void getRow(int row, int& rowlen, const int*& rowinds,
              const double*& rowvals) const {
    int start = ARstart_[row];
    rowlen = ARstart_[row + 1] - start;
    rowinds = ARindex_.data() + start;
    rowvals = ARvalue_.data() + start;
  }

  void computeMaxAbsCoefs() const {
    maxAbsRowCoef_.resize(numRow_);
    for (int i = 0; i != numRow_; ++i) {
      double maxabsval = 0.0;

      int start = ARstart_[i];
      int end = ARstart_[i + 1];

      for (int j = start; j != end; ++j)
        maxabsval = std::max(maxabsval, std::abs(ARvalue_[j]));

      maxAbsRowCoef_[i] = maxabsval;
    }
  }
};

void solveMip(const HighsMip& mip);

#endif