#include "DoubleQuantizer.h"

namespace highs {

namespace folding {

double DoubleQuantizer::bucket(double d) const {
  // Do not use llround and return int64_t, because that is undefined
  // behaviour if d is infinite or too large.
  return std::round(d * inv_tol);
}

size_t DoubleQuantizer::hash(double d) const {
  return std::hash<double>()(bucket(d));
}

size_t DoubleQuantizer::hash(HighsInt i) const {
  return std::hash<HighsInt>()(i);
}

void DoubleQuantizer::combine(size_t& x, size_t y) const {
  // hash_combine from boost
  x ^= y + 0x9e3779b97f4a7c15 + (x << 6) + (x >> 2);
}

bool DoubleQuantizer::equal(double x, double y) const {
  return bucket(x) == bucket(y);
}

size_t DoubleQuantizer::operator()(double d) const { return hash(d); }

size_t DoubleQuantizer::operator()(const RowData& ri) const {
  size_t h = hash(ri.bl);
  combine(h, hash(ri.bu));
  combine(h, hash(ri.colour));
  return h;
}

size_t DoubleQuantizer::operator()(const ColData& ci) const {
  size_t h = hash(ci.c);
  combine(h, hash(ci.l));
  combine(h, hash(ci.u));
  combine(h, hash(ci.colour));
  return h;
}

bool DoubleQuantizer::operator()(double x, double y) const {
  return equal(x, y);
}

bool DoubleQuantizer::operator()(const ColData& x, const ColData& y) const {
  return x.colour == y.colour && equal(x.c, y.c) && equal(x.l, y.l) &&
         equal(x.u, y.u);
}

bool DoubleQuantizer::operator()(const RowData& x, const RowData& y) const {
  return x.colour == y.colour && equal(x.bl, y.bl) && equal(x.bu, y.bu);
}

}  // namespace folding
}  // namespace highs