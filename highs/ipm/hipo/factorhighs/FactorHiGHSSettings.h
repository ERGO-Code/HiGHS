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

// Tuning parameters, historically compile-time constants. They are now held in
// a single process-wide instance so that they can be overridden at runtime via
// the advanced HiGHS options hipo_* (see HighsOptions.h). The instance is
// written once, from Solver::load() during single-threaded setup, and is
// read-only while a solve is in progress.
struct HipoTuning {
  // supernode amalgamation
  Int sn_start_thresh = 256;         // was kStartThreshRelax
  double sn_upper_ratio = 0.02;      // was kUpperRatioRelax
  double sn_lower_ratio = 0.01;      // was kLowerRatioRelax
  Int sn_max_iter_relax = 20;        // was kMaxIterRelax
  Int sn_size_relax = 16;            // was kSnSizeRelax
  double sn_spops_weight = 50.0;     // was kSpopsWeightSn

  // dense factorisation
  double alpha_bk = 0.01;            // was kAlphaBK; (sqrt(17.0) + 1.0) / 8.0
  Int block_grain_size = 1;          // was kBlockGrainSize
  Int block_parallel_threshold = 5;  // was kBlockParallelThreshold

  Int min_consecutive_sums = 1;      // was kMinConsecutiveSums

  // regularisation
  double dynamic_reg_coeff = 1e-24;  // was kDynamicDiagCoeff

  // metis
  Int metis_seed = 42;               // was kMetisSeed
};

HipoTuning& hipoTuning();

struct Regul {
  double primal{};
  double dual{};
};

}  // namespace hipo

#endif