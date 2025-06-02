#ifndef FACTORHIGHS_SYM_SCALING_H
#define FACTORHIGHS_SYM_SCALING_H

#include <cmath>
#include <vector>

#include "ipm/hpm/auxiliary/IntConfig.h"

namespace highspm {

// Scalings for symmetric matrices, provided as lower triangular.

void CurtisReidScalingSym(const std::vector<Int>& ptr,
                          const std::vector<Int>& rows,
                          const std::vector<double>& val,
                          std::vector<double>& colscale);

void RuizScalingSym(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                    const std::vector<double>& val,
                    std::vector<double>& colscale);

void JacekScalingSym(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                     const std::vector<double>& val,
                     std::vector<double>& colscale);

}  // namespace highspm

#endif