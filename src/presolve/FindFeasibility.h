#ifndef PRESOVLE_FINDFEASIBILITY_H_
#define PRESOVLE_FINDFEASIBILITY_H_

#include "lp_data/HighsStatus.h"
#include "lp_data/HighsLp.h"

HighsStatus runIdiot(const HighsLp& lp, const HighsSolution& solution);

#endif