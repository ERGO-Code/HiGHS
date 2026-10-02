#ifndef PDLP_HIPDLP_MODEL_SCALING_HPP_
#define PDLP_HIPDLP_MODEL_SCALING_HPP_

#include <vector>

#include "lp_data/HConst.h"
#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace model {

const HighsInt kDefaultCurtisReidSweeps = 20;

// Auto scaling discards Curtis-Reid when, compared with the same pipeline
// without it, the scaled costs, right-hand sides or finite bounds would have a
// 2-norm more than this factor larger.
//
// Curtis-Reid equalises the geometric mean of every row and column. A row
// holding a few unit entries among thousands of tiny ones gets a scale that
// is set by the tiny entries, and the columns it meets are scaled the other
// way. The solver is then close to converged in scaled space while the
// unscaled residual, which is what the stopping test measures, is amplified by
// the same factors. Ruiz and Pock-Chambolle weight the large entries instead,
// and do not do this.
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

  bool ok() const;
};

CurtisReidResult computeCurtisReidScaling(const SparseMatrix& matrix,
                                          HighsInt sweeps);
api::Error applyCurtisReidScaling(Problem& problem, HighsInt sweeps);

ScalingResult computeRuizScaling(const SparseMatrix& matrix,
                                 HighsInt iterations);
ScalingResult computePockChambolleScaling(const SparseMatrix& matrix);
ScalingResult computeL2Scaling(const SparseMatrix& matrix);
ScalingResult computeScalingPipeline(const SparseMatrix& matrix,
                                     const ScalingOptions& options);
// Largest factor, over the cost vector, the finite row lower bounds and the
// finite column bounds of the unscaled `problem`, by which the 2-norm after
// applying `candidate` exceeds the 2-norm after applying `reference`.
double scaledNormInflation(const Problem& problem,
                           const ScalingMetadata& candidate,
                           const ScalingMetadata& reference);
api::Error applyScalingPipeline(Problem& problem,
                                const ScalingOptions& options);

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
