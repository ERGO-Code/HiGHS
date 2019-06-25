#ifndef PRESOVLE_FINDFEASIBILITY_H_
#define PRESOVLE_FINDFEASIBILITY_H_

#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"

enum class MinimizationType {
  kComponentWise,
  kExact,
};

HighsStatus runFeasibility(const HighsLp& lp, HighsSolution& solution,
                           const MinimizationType type);

#endif