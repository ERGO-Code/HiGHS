#ifndef HIPO_PARAMETERS_H
#define HIPO_PARAMETERS_H

#include <vector>

#include "ipm/hipo/auxiliary/IntConfig.h"

namespace hipo {

// parameters for termination
const Int kMaxIterDefault = 200;
const double kIpmTolDefault = 1e-8;
const double kDivergeTol = 1e3;
const double kCertificateTol = 1e-6;
const double kCertificateMinRadius = 1e3;
const double kKktStagnationRelativeViolation = 10.0;

// parameters for residual reduction test
const double kResReductionFraction = 0.1;
const double kResReductionMinAlpha = 1e-3;
const double kResReductionLargeIncrease = 1e3;
const double kResReductionNoiseFactor = 1e2;
const Int kBadResWindow = 5;
const Int kBadResMaxIter = 5;

// parameters for progress test
const double kProgressFactor = 0.95;
const Int kProgressMaxIter = 5;
const Int kProgressBestMaxIter = 20;
const double kFrozenTol = 1e-3;

// parameters for correctors
const double kGammaCorrector = 0.1;
const double kSigmaAffine = 0.01;
const double kSigmaAffineSmallStep = 0.1;
const Int kMaxCorrectors = 10;
const double kMccIncreaseAlpha = 0.1;
const double kMccIncreaseMin = 0.1;
const double kSmallProduct = 1e-3;
const double kLargeProduct = 1e3;

// parameters for signals collected at each iteration
const double kPivotRegDominatedRatio = 1.0;
const double kNearFreeRatio = 1e3;

// parameters for numerical detectors
const double kInaccurateSolveTol = 1e-6;
const double kVeryInaccurateSolveTol = 1e-2;
const Int kWrongInertiaMinPivots = 20;
const double kWrongInertiaFraction = 1e-3;

// parameters for centrality detector
const double kPoorCentralityFraction = 0.01;
const double kPoorCentralityAlpha = 0.1;

// parameters for recovery policy
const Int kActionPatience = 3;
const Int kActionNumericalWindow = 5;
const Int kMaxRegChanges = 2;
const Int kCleanIterBeforeRegDecrease = 5;
const double kShortStepAlpha = 0.1;
const double kRegScaleFactor = 100.0;
const double kMinRegFactor = 1e-3;

// parameters for nearly free columns
const double kGuardAlpha = 0.01;
const double kGuardRatio = 0.1;
const double kGuardRegFactor = 30.0;
const Int kGuardPasses = 2;

// parameters for centring
const double kRecentreSigmaMin = 0.5;
const Int kRecentreIter = 2;
const double kCentredProductRatio = 0.1;

// parameters for choice of AS or NE
const double kSystemSpopsWeight = 30.0;
const double kSystemRatioOpsThresh = 10.0;
const double kSystemRatioSnThresh = 1.5;
const double kSystemSymbNzMult = 5.0;

// parameters for choice of parallelism
const double kParallelLargeFlopsThresh = 1e7;
const double kParallelLargeSpeedupThresh = 1;
const double kParallelLargeSnThresh = 20.0;
const double kParallelSmallSnThresh = 5.0;
const Int kParallelMinNumberSn = 10;
const double kParallelLargeStorageGB = 20.0;
const double kParallelMaxTreeDepth = 1000;

// parameters for choice of ordering
const double kOrderingFlopsThresh = 1.2;

// parameters for choice of factorisation
const double kUplookFlopsThresh = 1e6;
const double kUplookNzPerColLower = 10;
const double kUplookNzPerColUpper = 25;
const double kUplookSpopsRatioLower = 20;
const double kUplookSpopsRatioUpper = 100;

// parameters for uplooking factorisation
const Int64 kUplookTimeCheckWork = 10000000;

// parameters for skipping AS or NE
const double kSkipSystemNzBoundsRatio = 50.0;

// margins over the serial storage of the factorisation when checking memory.
// Measured peaks were up to 1.15 times the storage in serial, and 2.3 times
// with a parallel tree. These are conservative guesses, to be refined
const double kMemoryMarginSerial = 1.2;
const double kMemoryMarginParallelTree = 2.5;

// Below this many bytes required, memory is not checked, since reading the
// available memory costs about 30us of file reads on Linux
const double kMemoryCheckMinBytes = 256.0 * 1024 * 1024;

// parameters for iterative refinement
const Int kRefineMaxIter = 3;
const double kRefineTol = 1e-12;
const double kRefineMinDecrease = 0.5;

// parameters for scaling
const double kScalingSmallCoeff = 1e-4;
const double kScalingLargeCoeff = 1e4;
const double kScalingSmallBoundDiff = 1e-3;

// parameters for free variables
const double kFreeVarsInitialBound = 1e4;
const double kFreeVarsCloseRatio = 0.5;

// parameters for parallel NE
const Int kParallelNEStructTasks = 50;   // 32 < . <= 64
const Int kParallelNEValuesTasks = 100;  // 64 < . <= 128
const Int kParallelNEnzPerColThresh = 10;
const Int kParallelNEnzPerRowThresh = 30;
const Int kParallelNEsizeThresh = 1e4;

// parameters for parallel solve
const double kParallelSolveMinSize = 1e4;
const double kParallelForwardMinSpeedup = 2;
const double kParallelBackwardMinSpeedup = 1.2;

// default static regularisation
const double kStaticRegPrimal = 1e-12;
const double kStaticRegDual = 1e-10;

struct SolverControls {
  // static regularisation: the defaults times reg_x_scale, times a factor for
  // each row/column
  double reg_p_scale = 1.0;
  double reg_d_scale = 1.0;
  std::vector<double> reg_p_factor;
  std::vector<double> reg_d_factor;
  double baseRegP() const { return reg_p_scale * kStaticRegPrimal; }
  double baseRegD() const { return reg_d_scale * kStaticRegDual; }
  double regP(Int j) const { return baseRegP() * reg_p_factor[j]; }
  double regD(Int i) const { return baseRegD() * reg_d_factor[i]; }

  Int refine_max_iter = kRefineMaxIter;
  double refine_tol = kRefineTol;
  double step_fraction = 0.9;
  double max_step = 1.0 - 1e-4;
  Int max_correctors = kMaxCorrectors;
};

}  // namespace hipo

#endif