#ifndef HIPO_CENTRING_H
#define HIPO_CENTRING_H

#include <algorithm>

#include "IterationHistory.h"
#include "Parameters.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// All products xj * zj / mu within [kCentredProductRatio,
// 1/kCentredProductRatio]
inline bool isCentred(const IterationRecord& record) {
  return record.min_prod >= kCentredProductRatio &&
         record.max_prod <= 1.0 / kCentredProductRatio;
}

// Centring mode: while active, sigma is at least sigma_min. It ends after a
// number of iterations or, if until_centred, as soon as the point is centred.
struct Centring {
  bool active = false;
  double sigma_min = 0.0;
  Int iters_left = 0;
  bool until_centred = false;

  void start(double sigma, Int iters, bool stop_when_centred) {
    active = true;
    sigma_min = sigma;
    iters_left = iters;
    until_centred = stop_when_centred;
  }

  // sigma to use, given the one chosen by the ipm
  double apply(double sigma) const {
    return active ? std::max(sigma, sigma_min) : sigma;
  }

  // Count an iteration done in centring mode, whose record is given, and end
  // the mode if it is finished. Return true if it ended.
  bool step(const IterationRecord& record) {
    if (!active) return false;
    --iters_left;
    if (iters_left > 0 && !(until_centred && isCentred(record))) return false;
    active = false;
    return true;
  }
};

}  // namespace hipo

#endif
