#ifndef FACTOR_HIGHS_OPTIONS_H
#define FACTOR_HIGHS_OPTIONS_H

#include "FactorHighsSettings.h"
#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

struct FHoptions {
  double reg_p = 0.0;
  double reg_d = 0.0;
  bool pivoting = true;
  Int nb = kBlockSize;
  bool one_indexing = false;
};

}  // namespace hipo

#endif