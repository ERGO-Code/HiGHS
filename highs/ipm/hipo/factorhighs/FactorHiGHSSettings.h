#ifndef FACTOR_HIGHS_SETTINGS_H
#define FACTOR_HIGHS_SETTINGS_H

#include <cmath>

#include "ipm/hipo/auxiliary/IntConfig.h"

// ===========================================================================
// SWITCHES
// ===========================================================================

// Switch on/off pivoting. It uses a static variation of Bunch-Kaufman pivoting,
// with potential dynamic regularisation. If pivoting is switched off, only
// static regularisation is applied.
#define HIPO_PIVOTING

// Collect data during regularisation, e.g. number of regularised pivots, 2x2
// pivots, pivot swaps, pivots with wrong sign, min and max entry of L and D.
// This can be quite expensive and should only be used for debugging.
// #define HIPO_COLLECT_EXPENSIVE_DATA

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
const double kDynamicDiagCoeff = 1e-24;

// refinement
const Int kMaxRefinementIter = 3;
const double kRefinementTolerance = 1e-12;

struct Regul {
  double primal{};
  double dual{};
};

}  // namespace hipo

#endif