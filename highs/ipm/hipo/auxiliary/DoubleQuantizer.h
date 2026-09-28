#ifndef HIGHS_DOUBLE_QUANTIZER_H
#define HIGHS_DOUBLE_QUANTIZER_H

#include <unordered_map>

#include "FoldConstants.h"
#include "lp_data/HighsLp.h"
#include "util/HighsType.h"

namespace highs {

namespace folding {

struct RowData {
  double bl;
  double bu;
  HighsInt colour;

  RowData(const HighsLp& lp, HighsInt row, HighsInt c)
      : bl{lp.row_lower_[row]}, bu{lp.row_upper_[row]}, colour{c} {}
};

struct ColData {
  double c;
  double l;
  double u;
  HighsInt colour;

  ColData(const HighsLp& lp, HighsInt col, HighsInt c)
      : c{lp.col_cost_[col]},
        l{lp.col_lower_[col]},
        u{lp.col_upper_[col]},
        colour{c} {}
};

/*
Define hash and equality operators, so that variables of type double can be
grouped together if they are "close enough". This is based on an absolute
tolerance.
Must guarantee that `equal(a, b) => hash(a) == hash(b)` for
std::unordered_map to work correctly.
Two doubles are considered equivalent if they fall in the same cell of a
fixed-width grid:
    bucket(x) = round(x / tolerance)
    a ~ b     <=>  bucket(a) == bucket(b)
The same `bucket()` function drives both the hash and the equality, so
`equal(a, b) => hash(a) == hash(b)` holds exactly. This is possible because
"same grid cell" is an equivalence relation (unlike "|a - b| < tol", which is
not transitive).
It's possible that numbers that are close together end up in different buckets
and are therefore not considered equivalent. In order to reduce this effect on
real data, the cells are centred on integer multiples of the tolerance, so
values typed with at most as many decimal places as the tolerance
(e.g. 2.0, 1.5, 3.14159 for a tolerance of 1e-6) sit at the centre of their
cell.
`std::round` leaves infinities unchanged, so all +inf values share one bucket
and all -inf values share another.
Since the tolerance is absolute, values larger than around 2^52 * tolerance are
never merged, because the double format can no longer represent values a
tolerance apart.
*/

struct DoubleQuantizer {
  const double tolerance;

  DoubleQuantizer(double tol = kQuantizationTolerance) : tolerance{tol} {}

  double bucket(double d) const;

  size_t hash(double d) const;
  size_t hash(HighsInt i) const;
  void combine(size_t& x, size_t y) const;

  bool equal(double x, double y) const;

  // hashing operators
  size_t operator()(double x) const;
  size_t operator()(const RowData& ri) const;
  size_t operator()(const ColData& ci) const;

  // equality operators
  bool operator()(double x, double y) const;
  bool operator()(const RowData& x, const RowData& y) const;
  bool operator()(const ColData& x, const ColData& y) const;
};

template <typename P, typename Q>
using QuantizedMap = std::unordered_map<P, Q, DoubleQuantizer, DoubleQuantizer>;

}  // namespace folding
}  // namespace highs

#endif