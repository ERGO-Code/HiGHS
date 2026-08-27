#ifndef PDLP_HIPDLP_MODEL_SCALING_HPP_
#define PDLP_HIPDLP_MODEL_SCALING_HPP_

#include <vector>

#include "lp_data/HConst.h"
#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace model {

const HighsInt kDefaultCurtisReidSweeps = 20;

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
