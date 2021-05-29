#ifndef PRESOLVE_ICRASHX_H_
#define PRESOLVE_ICRASHX_H_

#include <iostream>

#include "HConfig.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsSolution.h"

bool callCrossover(const HighsLp& lp, const HighsOptions& options,
                   HighsSolution& solution, HighsBasis& highs_basis);

#endif