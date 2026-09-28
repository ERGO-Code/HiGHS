#include "DoubleQuantizer.h"

namespace highs {

namespace folding {

double DoubleQuantizer::bucket(double d) {
  // Do not use llround and return int64_t, because it is undefined
  // behaviour if d is infinite or too large.
  return std::round(d * inv_tol);
}

void DoubleQuantizer::hash_combine(size_t& x, size_t y) const {
  // hash_combine from boost
  x ^= y + 0x9e3779b97f4a7c15 + (x << 6) + (x >> 2);
}

size_t DoubleQuantizer::operator()(double d) const {
  return std::hash<double>()(bucket(d));
}

size_t DoubleQuantizer::operator()(const RowData& ri) const {
  size_t h = std::hash<double>()(ri.quantized_bl);
  hash_combine(h, std::hash<double>()(ri.quantized_bu));
  hash_combine(h, std::hash<HighsInt>()(ri.colour));
  return h;
}

size_t DoubleQuantizer::operator()(const ColData& ci) const {
  size_t h = std::hash<double>()(ci.quantized_c);
  hash_combine(h, std::hash<double>()(ci.quantized_l));
  hash_combine(h, std::hash<double>()(ci.quantized_u));
  hash_combine(h, std::hash<HighsInt>()(ci.colour));
  return h;
}

bool DoubleQuantizer::operator()(double x, double y) const {
  return bucket(x) == bucket(y);
}

bool DoubleQuantizer::operator()(const ColData& x, const ColData& y) const {
  return x.colour == y.colour && x.quantized_c == y.quantized_c &&
         x.quantized_l == y.quantized_l && x.quantized_u == y.quantized_u;
}

bool DoubleQuantizer::operator()(const RowData& x, const RowData& y) const {
  return x.colour == y.colour && x.quantized_bl == y.quantized_bl &&
         x.quantized_bu == y.quantized_bu;
}

}  // namespace folding
}  // namespace highs