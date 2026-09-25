#ifndef HIGHS_FOLD_AUXILIARY_H
#define HIGHS_FOLD_AUXILIARY_H

#include <unordered_map>

#include "FoldConstants.h"
#include "util/HighsType.h"

namespace highs {

namespace folding {

struct RowInfo {
  double bl;
  double bu;
  HighsInt colour;
};

struct ColInfo {
  double c;
  double l;
  double u;
  HighsInt colour;
};

struct Hasher {
  const double tolerance;

  Hasher(double tol = kHashingTolerance) : tolerance{tol} {}

  size_t hash(double d) const;
  size_t hash(HighsInt i) const;
  void combine(size_t& x, size_t y) const;
  size_t operator()(double x) const;
  size_t operator()(const RowInfo& ri) const;
  size_t operator()(const ColInfo& ci) const;
};

struct DoubleEqual {
  const double tolerance;

  DoubleEqual(double tol = kEqualityTolerance) : tolerance{tol} {}

  bool areClose(double x, double y) const;
  bool operator()(double x, double y) const;
  bool operator()(const RowInfo& x, const RowInfo& y) const;
  bool operator()(const ColInfo& x, const ColInfo& y) const;
};

template <typename T>
using Map = std::unordered_map<T, double, Hasher, DoubleEqual>;

}  // namespace folding
}  // namespace highs

#endif