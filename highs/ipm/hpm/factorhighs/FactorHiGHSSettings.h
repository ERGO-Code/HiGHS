#ifndef FACTOR_HIGHS_SETTINGS_H
#define FACTOR_HIGHS_SETTINGS_H

#include <cmath>

#include "ipm/hpm/auxiliary/IntConfig.h"

// ===========================================================================
// SWITCHES
// ===========================================================================

#define PIVOTING
// #define HPM_COLLECT_EXPENSIVE_DATA

// choose level of timing:
// - 0: no timing
// - 1: basic timing
// - 2: advanced timing
// - 3: extreme timing (timing of each BLAS call, considerably slower)
#define HPM_TIMING_LEVEL 0

// ===========================================================================
// PARAMETERS
// ===========================================================================

namespace highspm {

// supernode amalgamation
const Int kStartThreshRelax = 256;
const double kUpperRatioRelax = 0.02;
const double kLowerRatioRelax = 0.01;
const Int kMaxIterRelax = 10;
const Int kSnSizeRelax = 16;

// dense factorisation
const Int kBlockSize = 128;
const double kAlphaBK = 0.01;  //(sqrt(17.0) + 1.0) / 8.0;
const Int kBlockGrainSize = 1;
const Int kBlockParallelThreshold = 5;

// regularisation
const double kPrimalStaticRegularisation = 1e-12;
const double kDualStaticRegularisation = 1e-10;
const double kDynamicDiagCoeff = 1e-24;

// refinement
const Int kMaxRefinementIter = 3;
const double kRefinementTolerance = 1e-12;

}  // namespace highspm

#endif