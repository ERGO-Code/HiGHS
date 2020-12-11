#ifndef HIGHS_SPARSE_VECTOR_SUM_H_
#define HIGHS_SPARSE_VECTOR_SUM_H_

#include <algorithm>
#include <cstdint>
#include <vector>

#include "util/HighsCDouble.h"

class HighsSparseVectorSum {
 public:
  std::vector<uint8_t> nonzeroflag;
  std::vector<HighsCDouble> values;
  std::vector<int> nonzeroinds;
  HighsSparseVectorSum() = default;

  HighsSparseVectorSum(int dimension) { setDimension(dimension); }

  void setDimension(int dimension) {
    values.resize(dimension);
    nonzeroflag.resize(dimension);
    nonzeroinds.reserve(dimension);
  }

  void add(int index, double value) {
    if (nonzeroflag[index]) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroflag[index] = true;
      nonzeroinds.push_back(index);
    }
  }

  void add(int index, HighsCDouble value) {
    if (nonzeroflag[index]) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroflag[index] = true;
      nonzeroinds.push_back(index);
    }
  }

  void set(int index, double value) {
    values[index] = value;
    nonzeroinds.push_back(index);
  }

  void set(int index, HighsCDouble value) {
    values[index] = value;
    nonzeroinds.push_back(index);
  }

  void chgValue(int index, double val) { values[index] = val; }
  void chgValue(int index, HighsCDouble val) { values[index] = val; }

  const std::vector<int>& getNonzeros() const { return nonzeroinds; }

  double getValue(int index) const { return double(values[index]); }

  void clear() {
    for (int i : nonzeroinds) nonzeroflag[i] = false;

    nonzeroinds.clear();
  }

  void sort() { std::sort(nonzeroinds.begin(), nonzeroinds.end()); }
};

#endif