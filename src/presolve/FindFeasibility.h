#ifndef PRESOVLE_FINDFEASIBILITY_H_
#define PRESOVLE_FINDFEASIBILITY_H_

#include "lp_data/HighsStatus.h"
#include "lp_data/HighsLp.h"

HighsStatus runFeasibility(const HighsLp& lp,
                           HighsSolution& solution);

#endif