#ifndef HIPO_CONST_H
#define HIPO_CONST_H

#include "ipm/hpm/auxiliary/IntConfig.h"

namespace hipo {

// parameters for termination
const Int kMaxIterDefault = 200;
const double kIpmTolDefault = 1e-8;
const Int kMaxBadIter = 5;
const double kDivergeTol = 1e3;

// parameters for correctors
const double kGammaCorrector = 0.1;
const double kSigmaAffine = 0.01;
const Int kMaxCorrectors = 10;
const double kMccIncreaseAlpha = 0.1;
const double kMccIncreaseMin = 0.1;
const double kSmallProduct = 1e-3;
const double kLargeProduct = 1e3;

// parameters for choice of AS or NE
const double kSpopsWeight = 30.0;
const double kRatioOpsThresh = 10.0;
const double kRatioSnThresh = 1.5;
const double kSymbNzMult = 5.0;

// parameters for choice of parallelism
const double kLargeFlopsThresh = 1e7;
const double kLargeSpeedupThresh = 1.5;
const double kLargeSnThresh = 20.0;
const double kSmallSnThresh = 5.0;
const Int kMinNumberSn = 10;

// parameters for dense columns
const double kDenseColThresh = 0.5;
const Int kMinRowsForDensity = 2000;

}  // namespace hipo

#endif