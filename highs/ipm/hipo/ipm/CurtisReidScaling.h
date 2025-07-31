#ifndef HIPO_CURTIS_REID_SCALING_H
#define HIPO_CURTIS_REID_SCALING_H

#include <cmath>
#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"
#include "ipm/hipo/auxiliary/VectorOperations.h"

namespace hipo {

Int CurtisReidScaling(const std::vector<Int>& ptr, const std::vector<Int>& rows,
                      const std::vector<double>& val, std::vector<Int>& rowexp,
                      std::vector<Int>& colexp);

}

#endif