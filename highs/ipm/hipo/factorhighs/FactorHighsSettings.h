#ifndef FACTOR_HIGHS_SETTINGS_H
#define FACTOR_HIGHS_SETTINGS_H

#include <cmath>

#include "ipm/hipo/auxiliary/IntConfig.h"

// ===========================================================================
// SWITCHES
// ===========================================================================

// Choose level of timing:
// - 0: no timing
// - 1: basic timing
// - 2: advanced timing
// - 3: extreme timing (timing of each BLAS call, considerably slower)
#define HIPO_TIMING_LEVEL 0

// ===========================================================================
// PARAMETERS
// ===========================================================================

namespace hipo {

const Int kBlockSize = 128;

// supernode amalgamation
const Int kStartThreshRelax = 256;
const double kUpperRatioRelax = 0.02;
const double kLowerRatioRelax = 0.01;
const Int kMaxIterRelax = 20;
const Int kSnSizeRelax = 16;
const double kSpopsWeightSn = 50.0;

// dense factorisation
const double kAlphaBK = 0.01;  //(sqrt(17.0) + 1.0) / 8.0;
const Int kBlockGrainSize = 1;
const Int kBlockParallelThreshold = 5;

const Int kMinConsecutiveSums = 1;

// regularisation
const double kDynamicDiagCoeff = 1e-24;
const double kDroppedPivot = 1e128;

// parallel solve
const double kLargeTaskRelativeThresh = 0.01;
const double kLargeTaskAbsoluteThres = 1e3;
const double kParallelDiagTargetNumTasks = 32;

struct Regul {
  double primal{};
  double dual{};
};

}  // namespace hipo

#endif