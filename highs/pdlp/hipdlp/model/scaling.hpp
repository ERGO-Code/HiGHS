#ifndef PDLP_HIPDLP_MODEL_SCALING_HPP_
#define PDLP_HIPDLP_MODEL_SCALING_HPP_

#include <vector>

#include "lp_data/HConst.h"
#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace model {

const HighsInt kDefaultCurtisReidSweeps = 20;

// HPR-LP-C's auto guard: check each raw CR log-scale range before exponent
// clamping or applying CR. A rejected CR candidate leaves Ruiz/PC a pristine
// matrix; no second full scaling pipeline is needed. This bounds spread, not
// absolute scales or the change in costs/RHS/bounds, and is not a speedup test.
const double kCurtisReidMaximumScaleSpread = 1e8;

// Retained for the historical scaled-data-inflation diagnostic tools. Auto
// selection now uses the raw CR scale-spread guard above.
const double kCurtisReidInflationLimit = 1e3;

enum class ScalingMode { kAuto, kOff, kCurtisReid, kLegacy };

struct ScalingOptions {
  ScalingMode mode;
  HighsInt curtis_reid_sweeps;
  HighsInt legacy_scaling_mode;
  HighsInt ruiz_iterations;
  HighsInt features_off;

  ScalingOptions();
};

struct ScalingResult {
  api::Error error;
  ScalingMetadata scaling;

  bool ok() const;
};

struct CurtisReidResult {
  api::Error error;
  std::vector<double> row_scale;
  std::vector<double> col_scale;
  double row_log_range = 0.0;
  double col_log_range = 0.0;
  bool rejected_by_spread = false;

  bool ok() const;
  bool scaleSpreadIsSafe() const;
};

CurtisReidResult computeCurtisReidScaling(const SparseMatrix& matrix,
                                          HighsInt sweeps,
                                          bool guard_scale_spread = false);
api::Error applyCurtisReidScaling(Problem& problem, HighsInt sweeps);

ScalingResult computeRuizScaling(const SparseMatrix& matrix,
                                 HighsInt iterations,
                                 PreprocessingDiagnostics* diagnostics = nullptr);
ScalingResult computePockChambolleScaling(const SparseMatrix& matrix);
ScalingResult computeL2Scaling(const SparseMatrix& matrix);
ScalingResult computeScalingPipeline(const SparseMatrix& matrix,
                                     const ScalingOptions& options,
                                     PreprocessingDiagnostics* diagnostics = nullptr);
// Largest factor, over the cost vector, the finite row lower bounds and the
// finite column bounds of the unscaled `problem`, by which the 2-norm after
// applying `candidate` exceeds the 2-norm after applying `reference`.
double scaledNormInflation(const Problem& problem,
                           const ScalingMetadata& candidate,
                           const ScalingMetadata& reference);
api::Error applyScalingPipeline(Problem& problem,
                                const ScalingOptions& options,
                                PreprocessingDiagnostics* diagnostics = nullptr);

api::Error scalePrimal(const Problem& problem,
                       const std::vector<double>& unscaled,
                       std::vector<double>& scaled);
api::Error unscalePrimal(const Problem& problem,
                         const std::vector<double>& scaled,
                         std::vector<double>& unscaled);
api::Error scaleDual(const Problem& problem,
                     const std::vector<double>& unscaled,
                     std::vector<double>& scaled);
api::Error unscaleDual(const Problem& problem,
                       const std::vector<double>& scaled,
                       std::vector<double>& unscaled);
api::Error scaleReducedCosts(const Problem& problem,
                             const std::vector<double>& unscaled,
                             std::vector<double>& scaled);
api::Error unscaleReducedCosts(const Problem& problem,
                               const std::vector<double>& scaled,
                               std::vector<double>& unscaled);

}  // namespace model
}  // namespace hipdlp

#endif
