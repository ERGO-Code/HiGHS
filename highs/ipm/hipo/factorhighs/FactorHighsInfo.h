#ifndef FACTORHIGHS_INFO_H
#define FACTORHIGHS_INFO_H

#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// Information about the last factorisation
struct FHinfo {
  // 1x1 pivots before regularisation, not sign-adjusted (expected signs are
  // the signs passed to analyse). Same order as getRegularisation. 2x2: +inf.
  std::vector<double> unreg_pivots;

  Int num_2x2 = 0;

  // inertia, see FHsolver::inertia
  Int num_pos = 0;
  Int num_neg = 0;
  Int num_zero = 0;
};

}  // namespace hipo

#endif
