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
  x ^= y + kFibonacciMultiplier + (x << 6) + (x >> 2);
}

size_t DoubleQuantizer::operator()(const RowQuantizedData& ri) const {
  size_t h = std::hash<double>()(ri.quantized_bl);
  hash_combine(h, std::hash<double>()(ri.quantized_bu));
  return h;
}

size_t DoubleQuantizer::operator()(const ColQuantizedData& ci) const {
  size_t h = std::hash<double>()(ci.quantized_c);
  hash_combine(h, std::hash<double>()(ci.quantized_l));
  hash_combine(h, std::hash<double>()(ci.quantized_u));
  return h;
}

bool DoubleQuantizer::operator()(const ColQuantizedData& x,
                                 const ColQuantizedData& y) const {
  return x.quantized_c == y.quantized_c && x.quantized_l == y.quantized_l &&
         x.quantized_u == y.quantized_u;
}

bool DoubleQuantizer::operator()(const RowQuantizedData& x,
                                 const RowQuantizedData& y) const {
  return x.quantized_bl == y.quantized_bl && x.quantized_bu == y.quantized_bu;
}

}  // namespace folding
}  // namespace highs