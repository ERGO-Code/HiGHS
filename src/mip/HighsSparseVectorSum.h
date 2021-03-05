/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SPARSE_VECTOR_SUM_H_
#define HIGHS_SPARSE_VECTOR_SUM_H_

#include <algorithm>
#include <cassert>
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
    assert(index >= 0 && index < (int)nonzeroflag.size());
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

  template <typename Pred>
  int partition(Pred&& pred) {
    return std::partition(nonzeroinds.begin(), nonzeroinds.end(), pred) -
           nonzeroinds.begin();
  }

  template <typename IsZero>
  void cleanup(IsZero&& isZero) {
    int numNz = nonzeroinds.size();

    for (int i = numNz - 1; i >= 0; --i) {
      int pos = nonzeroinds[i];
      double val = double(values[pos]);

      if (isZero(pos, val)) {
        values[pos] = 0.0;
        nonzeroflag[pos] = 0;
        --numNz;
        std::swap(nonzeroinds[numNz], nonzeroinds[i]);
      }
    }

    nonzeroinds.resize(numNz);
  }
};

#endif