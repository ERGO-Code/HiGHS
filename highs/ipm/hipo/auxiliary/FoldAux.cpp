#include "FoldAux.h"

namespace highs {

namespace folding {

// Define hash and equality operators, so that variables of type double can be
// grouped together if they are close enough. This is achieved by rounding the
// number to an integer grid point (using multiples of tolerance).
// Notice that it is possible for two numbers to be close together (areClose
// returns true), but have different hash values. This may create weird
// behaviour in the unordered_map.
// By keeping kEqualityTolerance << kHashingTolerance, this risk should be
// acceptable.

size_t Hasher::hash(double d) const {
  // associate an integer "grid point" to d, and hash it
  if (std::abs(d) < tolerance) return 0;
  double offset = d < 0 ? -tolerance : tolerance;
  int64_t grid_point = std::round((d + offset) / tolerance);
  return std::hash<int64_t>()(grid_point);
}

size_t Hasher::hash(HighsInt i) const { return std::hash<HighsInt>()(i); }

void Hasher::combine(size_t& x, size_t y) const {
  // hash_combine from boost
  x ^= y + 0x9e3779b9 + (x << 6) + (x >> 2);
}

size_t Hasher::operator()(double d) const { return hash(d); }

size_t Hasher::operator()(const RowInfo& ri) const {
  size_t h = hash(ri.bl);
  combine(h, hash(ri.bu));
  combine(h, hash(ri.colour));
  return h;
}

size_t Hasher::operator()(const ColInfo& ci) const {
  size_t h = hash(ci.c);
  combine(h, hash(ci.l));
  combine(h, hash(ci.u));
  combine(h, hash(ci.colour));
  return h;
}

bool DoubleEqual::areClose(double x, double y) const {
  return std::abs(x - y) < tolerance;
}

bool DoubleEqual::operator()(double x, double y) const {
  return areClose(x, y);
}

bool DoubleEqual::operator()(const ColInfo& x, const ColInfo& y) const {
  return x.colour == y.colour && areClose(x.c, y.c) && areClose(x.l, y.l) &&
         areClose(x.u, y.u);
}

bool DoubleEqual::operator()(const RowInfo& x, const RowInfo& y) const {
  return x.colour == y.colour && areClose(x.bl, y.bl) && areClose(x.bu, y.bu);
}

Map<ColInfo> map;

}  // namespace folding
}  // namespace highs