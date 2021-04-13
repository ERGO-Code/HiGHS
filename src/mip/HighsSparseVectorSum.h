/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_SPARSE_VECTOR_SUM_H_
#define HIGHS_SPARSE_VECTOR_SUM_H_

#include <algorithm>
#include <cassert>
#include <vector>

#include "util/HighsCDouble.h"
#include "util/HighsInt.h"

class HighsSparseVectorSum {
 public:
  std::vector<uint8_t> nonzeroflag;
  std::vector<HighsCDouble> values;
  std::vector<HighsInt> nonzeroinds;
  HighsSparseVectorSum() = default;

  HighsSparseVectorSum(HighsInt dimension) { setDimension(dimension); }

  void setDimension(HighsInt dimension) {
    values.resize(dimension);
    nonzeroflag.resize(dimension);
    nonzeroinds.reserve(dimension);
  }

  void add(HighsInt index, double value) {
    assert(index >= 0 && index < (HighsInt)nonzeroflag.size());
    if (nonzeroflag[index]) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroflag[index] = true;
      nonzeroinds.push_back(index);
    }
  }

  void add(HighsInt index, HighsCDouble value) {
    if (nonzeroflag[index]) {
      values[index] += value;
    } else {
      values[index] = value;
      nonzeroflag[index] = true;
      nonzeroinds.push_back(index);
    }
  }

  void set(HighsInt index, double value) {
    values[index] = value;
    nonzeroinds.push_back(index);
  }

  void set(HighsInt index, HighsCDouble value) {
    values[index] = value;
    nonzeroinds.push_back(index);
  }

  void chgValue(HighsInt index, double val) { values[index] = val; }
  void chgValue(HighsInt index, HighsCDouble val) { values[index] = val; }

  const std::vector<HighsInt>& getNonzeros() const { return nonzeroinds; }

  double getValue(HighsInt index) const { return double(values[index]); }

  void clear() {
    for (HighsInt i : nonzeroinds) nonzeroflag[i] = false;

    nonzeroinds.clear();
  }

  void sort() { std::sort(nonzeroinds.begin(), nonzeroinds.end()); }

  template <typename Pred>
  HighsInt partition(Pred&& pred) {
    return std::partition(nonzeroinds.begin(), nonzeroinds.end(), pred) -
           nonzeroinds.begin();
  }

  template <typename IsZero>
  void cleanup(IsZero&& isZero) {
    HighsInt numNz = nonzeroinds.size();

    for (HighsInt i = numNz - 1; i >= 0; --i) {
      HighsInt pos = nonzeroinds[i];
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
