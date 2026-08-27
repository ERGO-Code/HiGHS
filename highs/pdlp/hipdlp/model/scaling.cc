#include "pdlp/hipdlp/model/scaling.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace hipdlp {
namespace model {

namespace {

const double kCoefficientFloor = 1e-300;
const double kMinimumScale = 1e-30;
const double kMaximumScale = 1e30;

api::Error invalidScaling(const std::string& message) {
  return api::Error(api::ErrorCode::kInvalidModel, message);
}

api::Error validateSparseMatrix(const SparseMatrix& matrix) {
  if (matrix.num_rows < 0 || matrix.num_cols < 0 ||
      matrix.starts.size() != static_cast<std::size_t>(matrix.num_rows) + 1 ||
      matrix.starts.empty() || matrix.starts[0] != 0)
    return invalidScaling("HiPDLP scaling sparse dimensions are malformed");
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    if (matrix.starts[row] < 0 || matrix.starts[row] > matrix.starts[row + 1])
      return invalidScaling("HiPDLP scaling sparse starts are not monotone");
  }
  if (matrix.starts[matrix.num_rows] < 0 ||
      matrix.indices.size() !=
          static_cast<std::size_t>(matrix.starts[matrix.num_rows]) ||
      matrix.values.size() != matrix.indices.size())
    return invalidScaling(
        "HiPDLP scaling sparse arrays have inconsistent "
        "lengths");
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    HighsInt previous_col = -1;
    for (HighsInt entry = matrix.starts[row]; entry < matrix.starts[row + 1];
         ++entry) {
      const HighsInt col = matrix.indices[entry];
      if (col < 0 || col >= matrix.num_cols || col <= previous_col)
        return invalidScaling(
            "HiPDLP scaling sparse indices are not canonical");
      if (!std::isfinite(matrix.values[entry]) || matrix.values[entry] == 0.0)
        return invalidScaling(
            "HiPDLP scaling requires finite stored nonzero values");
      previous_col = col;
    }
  }
  return api::Error();
}

double clampedExponential(const double log_scale) {
  const double minimum_log = std::log(kMinimumScale);
  const double maximum_log = std::log(kMaximumScale);
  return std::exp(std::max(minimum_log, std::min(log_scale, maximum_log)));
}

bool validScalingMode(const ScalingMode mode) {
  return mode == ScalingMode::kAuto || mode == ScalingMode::kOff ||
         mode == ScalingMode::kCurtisReid || mode == ScalingMode::kLegacy;
}

api::Error validateScalingOptions(const ScalingOptions& options) {
  if (!validScalingMode(options.mode))
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP scaling mode is invalid");
  if (options.curtis_reid_sweeps <= 0)
    return api::Error(api::ErrorCode::kInvalidOption,
                      "Curtis-Reid sweep count must be positive");
  if (options.legacy_scaling_mode < kPdlpScalingMin ||
      options.legacy_scaling_mode > kPdlpScalingMax)
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP legacy scaling bit mask is invalid");
  if (options.ruiz_iterations < 0)
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP Ruiz iteration count must be nonnegative");
  if (options.features_off < kPdlpAllFeaturesOn ||
      options.features_off > kPdlpAllFeaturesOff)
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP feature-disable bit mask is invalid");
  return api::Error();
}

bool positiveFiniteProduct(const double left, const double right,
                           double& result) {
  const long double product =
      static_cast<long double>(left) * static_cast<long double>(right);
  if (!std::isfinite(product) || product <= 0.0 ||
      product > std::numeric_limits<double>::max())
    return false;
  result = static_cast<double>(product);
  return std::isfinite(result) && result > 0.0;
}

bool finiteScaledValue(const double value, const double row_scale,
                       const double col_scale, double& result) {
  const long double scaled =
      static_cast<long double>(value) * row_scale * col_scale;
  if (!std::isfinite(scaled) ||
      std::abs(scaled) > std::numeric_limits<double>::max())
    return false;
  result = static_cast<double>(scaled);
  return std::isfinite(result) && result != 0.0;
}

api::Error composeScaling(SparseMatrix& matrix,
                          const std::vector<double>& row_scale,
                          const std::vector<double>& col_scale,
                          std::vector<double>& cumulative_row_scale,
                          std::vector<double>& cumulative_col_scale) {
  if (row_scale.size() != static_cast<std::size_t>(matrix.num_rows) ||
      col_scale.size() != static_cast<std::size_t>(matrix.num_cols) ||
      cumulative_row_scale.size() != row_scale.size() ||
      cumulative_col_scale.size() != col_scale.size())
    return invalidScaling("HiPDLP scale composition dimensions are invalid");

  std::vector<double> values(matrix.values.size());
  std::vector<double> composed_rows(row_scale.size());
  std::vector<double> composed_cols(col_scale.size());
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    if (!positiveFiniteProduct(cumulative_row_scale[row], row_scale[row],
                               composed_rows[row]))
      return invalidScaling("HiPDLP cumulative row scale is invalid");
    for (HighsInt entry = matrix.starts[row]; entry < matrix.starts[row + 1];
         ++entry) {
      const HighsInt col = matrix.indices[entry];
      if (!finiteScaledValue(matrix.values[entry], row_scale[row],
                             col_scale[col], values[entry]))
        return invalidScaling(
            "HiPDLP matrix value is unsafe during scale composition");
    }
  }
  for (HighsInt col = 0; col < matrix.num_cols; ++col) {
    if (!positiveFiniteProduct(cumulative_col_scale[col], col_scale[col],
                               composed_cols[col]))
      return invalidScaling("HiPDLP cumulative column scale is invalid");
  }
  matrix.values.swap(values);
  cumulative_row_scale.swap(composed_rows);
  cumulative_col_scale.swap(composed_cols);
  return api::Error();
}

ScalingResult makeScalingResult(const SparseMatrix& matrix,
                                const ScalingStage stage) {
  ScalingResult result;
  result.scaling.row_scale.assign(matrix.num_rows, 1.0);
  result.scaling.col_scale.assign(matrix.num_cols, 1.0);
  result.scaling.stages.push_back(stage);
  return result;
}

api::Error inverseSquareRoot(const double value, double& result) {
  if (!std::isfinite(value) || value < 0.0)
    return invalidScaling("HiPDLP scaling norm is invalid");
  result = value == 0.0 ? 1.0 : 1.0 / std::sqrt(value);
  if (!std::isfinite(result) || result <= 0.0)
    return invalidScaling("HiPDLP scaling factor is invalid");
  return api::Error();
}

void updateScaledSumOfSquares(const double value, double& scale,
                              double& sum_of_squares) {
  const double magnitude = std::abs(value);
  if (magnitude == 0.0) return;
  if (scale < magnitude) {
    const double ratio = scale / magnitude;
    sum_of_squares = 1.0 + sum_of_squares * ratio * ratio;
    scale = magnitude;
  } else {
    const double ratio = magnitude / scale;
    sum_of_squares += ratio * ratio;
  }
}

api::Error appendScaling(const ScalingResult& stage, SparseMatrix& matrix,
                         ScalingMetadata& cumulative) {
  if (!stage.ok()) return stage.error;
  api::Error error =
      composeScaling(matrix, stage.scaling.row_scale, stage.scaling.col_scale,
                     cumulative.row_scale, cumulative.col_scale);
  if (!error.ok()) return error;
  cumulative.stages.insert(cumulative.stages.end(),
                           stage.scaling.stages.begin(),
                           stage.scaling.stages.end());
  return api::Error();
}

enum class RoundTripKind { kPrimal, kDual, kReducedCost };

api::Error transformVector(const Problem& problem,
                           const std::vector<double>& input,
                           std::vector<double>& output,
                           const RoundTripKind kind, const bool scale) {
  if (!problem.initialized())
    return invalidScaling("HiPDLP problem is not initialized");
  const std::vector<double>& factors = kind == RoundTripKind::kDual
                                           ? problem.scaling().row_scale
                                           : problem.scaling().col_scale;
  if (input.size() != factors.size())
    return invalidScaling("HiPDLP scale round-trip dimension mismatch");
  std::vector<double> transformed(input.size());
  for (std::size_t index = 0; index < input.size(); ++index) {
    if (!std::isfinite(input[index]) || !std::isfinite(factors[index]) ||
        factors[index] <= 0.0)
      return invalidScaling("HiPDLP scale round-trip input is invalid");
    const bool multiply = (kind == RoundTripKind::kReducedCost) == scale;
    transformed[index] = multiply ? input[index] * factors[index]
                                  : input[index] / factors[index];
    if (!std::isfinite(transformed[index]))
      return invalidScaling("HiPDLP scale round-trip result is non-finite");
  }
  output.swap(transformed);
  return api::Error();
}

}  // namespace

ScalingOptions::ScalingOptions()
    : mode(ScalingMode::kAuto),
      curtis_reid_sweeps(kDefaultCurtisReidSweeps),
      legacy_scaling_mode(kPdlpScalingRuiz + kPdlpScalingPC),
      ruiz_iterations(10),
      features_off(kPdlpAllFeaturesOn) {}

bool ScalingResult::ok() const { return error.ok(); }

bool CurtisReidResult::ok() const { return error.ok(); }

CurtisReidResult computeCurtisReidScaling(const SparseMatrix& matrix,
                                          const HighsInt sweeps) {
  CurtisReidResult result;
  if (sweeps <= 0) {
    result.error = api::Error(api::ErrorCode::kInvalidOption,
                              "Curtis-Reid sweep count must be positive");
    return result;
  }
  result.error = validateSparseMatrix(matrix);
  if (!result.error.ok()) return result;

  std::vector<double> negative_logs(matrix.values.size());
  for (std::size_t entry = 0; entry < matrix.values.size(); ++entry) {
    negative_logs[entry] =
        -std::log(std::max(std::abs(matrix.values[entry]), kCoefficientFloor));
    if (!std::isfinite(negative_logs[entry])) {
      result.error =
          invalidScaling("Curtis-Reid coefficient logarithm is non-finite");
      return result;
    }
  }

  std::vector<double> row_log(matrix.num_rows, 0.0);
  std::vector<double> col_log(matrix.num_cols, 0.0);
  std::vector<double> col_sum(matrix.num_cols, 0.0);
  std::vector<HighsInt> col_count(matrix.num_cols, 0);
  for (HighsInt sweep = 0; sweep < sweeps; ++sweep) {
    for (HighsInt row = 0; row < matrix.num_rows; ++row) {
      const HighsInt count = matrix.starts[row + 1] - matrix.starts[row];
      if (count == 0) {
        row_log[row] = 0.0;
        continue;
      }
      double sum = 0.0;
      for (HighsInt entry = matrix.starts[row]; entry < matrix.starts[row + 1];
           ++entry) {
        sum += negative_logs[entry] - col_log[matrix.indices[entry]];
        if (!std::isfinite(sum)) {
          result.error =
              invalidScaling("Curtis-Reid row accumulation is non-finite");
          return result;
        }
      }
      row_log[row] = sum / static_cast<double>(count);
      if (!std::isfinite(row_log[row])) {
        result.error = invalidScaling("Curtis-Reid row scale is non-finite");
        return result;
      }
    }

    std::fill(col_sum.begin(), col_sum.end(), 0.0);
    std::fill(col_count.begin(), col_count.end(), 0);
    for (HighsInt row = 0; row < matrix.num_rows; ++row) {
      for (HighsInt entry = matrix.starts[row]; entry < matrix.starts[row + 1];
           ++entry) {
        const HighsInt col = matrix.indices[entry];
        col_sum[col] += negative_logs[entry] - row_log[row];
        ++col_count[col];
        if (!std::isfinite(col_sum[col])) {
          result.error =
              invalidScaling("Curtis-Reid column accumulation is non-finite");
          return result;
        }
      }
    }
    for (HighsInt col = 0; col < matrix.num_cols; ++col) {
      col_log[col] = col_count[col] == 0
                         ? 0.0
                         : col_sum[col] / static_cast<double>(col_count[col]);
      if (!std::isfinite(col_log[col])) {
        result.error = invalidScaling("Curtis-Reid column scale is non-finite");
        return result;
      }
    }
  }

  result.row_scale.resize(matrix.num_rows);
  result.col_scale.resize(matrix.num_cols);
  for (HighsInt row = 0; row < matrix.num_rows; ++row)
    result.row_scale[row] = clampedExponential(row_log[row]);
  for (HighsInt col = 0; col < matrix.num_cols; ++col)
    result.col_scale[col] = clampedExponential(col_log[col]);
  return result;
}

api::Error applyCurtisReidScaling(Problem& problem, const HighsInt sweeps) {
  if (!problem.initialized())
    return invalidScaling("HiPDLP problem is not initialized");
  CurtisReidResult scaling = computeCurtisReidScaling(problem.matrix(), sweeps);
  if (!scaling.ok()) return scaling.error;
  return problem.applyDiagonalScaling(scaling.row_scale, scaling.col_scale,
                                      sweeps);
}

ScalingResult computeRuizScaling(const SparseMatrix& matrix,
                                 const HighsInt iterations) {
  ScalingResult result = makeScalingResult(matrix, ScalingStage::kRuiz);
  result.scaling.ruiz_iterations = iterations;
  if (iterations < 0) {
    result.error = api::Error(api::ErrorCode::kInvalidOption,
                              "HiPDLP Ruiz iteration count is negative");
    return result;
  }
  result.error = validateSparseMatrix(matrix);
  if (!result.error.ok()) return result;

  SparseMatrix working = matrix;
  for (HighsInt iteration = 0; iteration < iterations; ++iteration) {
    std::vector<double> row_max(matrix.num_rows, 0.0);
    std::vector<double> col_max(matrix.num_cols, 0.0);
    for (HighsInt row = 0; row < working.num_rows; ++row) {
      for (HighsInt entry = working.starts[row];
           entry < working.starts[row + 1]; ++entry) {
        const HighsInt col = working.indices[entry];
        const double magnitude = std::abs(working.values[entry]);
        row_max[row] = std::max(row_max[row], magnitude);
        col_max[col] = std::max(col_max[col], magnitude);
      }
    }
    std::vector<double> row_scale(matrix.num_rows, 1.0);
    std::vector<double> col_scale(matrix.num_cols, 1.0);
    for (HighsInt row = 0; row < matrix.num_rows; ++row) {
      result.error = inverseSquareRoot(row_max[row], row_scale[row]);
      if (!result.error.ok()) return result;
    }
    for (HighsInt col = 0; col < matrix.num_cols; ++col) {
      result.error = inverseSquareRoot(col_max[col], col_scale[col]);
      if (!result.error.ok()) return result;
    }
    result.error =
        composeScaling(working, row_scale, col_scale, result.scaling.row_scale,
                       result.scaling.col_scale);
    if (!result.error.ok()) return result;
  }
  return result;
}

ScalingResult computePockChambolleScaling(const SparseMatrix& matrix) {
  ScalingResult result =
      makeScalingResult(matrix, ScalingStage::kPockChambolle);
  result.error = validateSparseMatrix(matrix);
  if (!result.error.ok()) return result;

  std::vector<double> row_sum(matrix.num_rows, 0.0);
  std::vector<double> col_sum(matrix.num_cols, 0.0);
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    for (HighsInt entry = matrix.starts[row]; entry < matrix.starts[row + 1];
         ++entry) {
      const HighsInt col = matrix.indices[entry];
      const double magnitude = std::abs(matrix.values[entry]);
      row_sum[row] += magnitude;
      col_sum[col] += magnitude;
      if (!std::isfinite(row_sum[row]) || !std::isfinite(col_sum[col])) {
        result.error =
            invalidScaling("Pock-Chambolle accumulation is non-finite");
        return result;
      }
    }
  }
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    result.error =
        inverseSquareRoot(row_sum[row], result.scaling.row_scale[row]);
    if (!result.error.ok()) return result;
  }
  for (HighsInt col = 0; col < matrix.num_cols; ++col) {
    result.error =
        inverseSquareRoot(col_sum[col], result.scaling.col_scale[col]);
    if (!result.error.ok()) return result;
  }
  return result;
}

ScalingResult computeL2Scaling(const SparseMatrix& matrix) {
  ScalingResult result = makeScalingResult(matrix, ScalingStage::kL2);
  result.error = validateSparseMatrix(matrix);
  if (!result.error.ok()) return result;

  std::vector<double> row_norm_scale(matrix.num_rows, 0.0);
  std::vector<double> row_sum_squares(matrix.num_rows, 1.0);
  std::vector<double> col_norm_scale(matrix.num_cols, 0.0);
  std::vector<double> col_sum_squares(matrix.num_cols, 1.0);
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    for (HighsInt entry = matrix.starts[row]; entry < matrix.starts[row + 1];
         ++entry) {
      const HighsInt col = matrix.indices[entry];
      updateScaledSumOfSquares(matrix.values[entry], row_norm_scale[row],
                               row_sum_squares[row]);
      updateScaledSumOfSquares(matrix.values[entry], col_norm_scale[col],
                               col_sum_squares[col]);
    }
  }
  for (HighsInt row = 0; row < matrix.num_rows; ++row) {
    const double norm =
        row_norm_scale[row] == 0.0
            ? 0.0
            : row_norm_scale[row] * std::sqrt(row_sum_squares[row]);
    result.error = inverseSquareRoot(norm, result.scaling.row_scale[row]);
    if (!result.error.ok()) return result;
  }
  for (HighsInt col = 0; col < matrix.num_cols; ++col) {
    const double norm =
        col_norm_scale[col] == 0.0
            ? 0.0
            : col_norm_scale[col] * std::sqrt(col_sum_squares[col]);
    result.error = inverseSquareRoot(norm, result.scaling.col_scale[col]);
    if (!result.error.ok()) return result;
  }
  return result;
}

ScalingResult computeScalingPipeline(const SparseMatrix& matrix,
                                     const ScalingOptions& options) {
  ScalingResult result;
  result.error = validateScalingOptions(options);
  if (!result.error.ok()) return result;
  result.error = validateSparseMatrix(matrix);
  if (!result.error.ok()) return result;
  result.scaling.row_scale.assign(matrix.num_rows, 1.0);
  result.scaling.col_scale.assign(matrix.num_cols, 1.0);

  if (options.mode == ScalingMode::kOff ||
      (options.features_off & kPdlpScalingOff) != 0)
    return result;

  SparseMatrix working = matrix;
  if (options.mode == ScalingMode::kAuto ||
      options.mode == ScalingMode::kCurtisReid) {
    const CurtisReidResult curtis_reid =
        computeCurtisReidScaling(working, options.curtis_reid_sweeps);
    if (!curtis_reid.ok()) {
      result.error = curtis_reid.error;
      return result;
    }
    result.error =
        composeScaling(working, curtis_reid.row_scale, curtis_reid.col_scale,
                       result.scaling.row_scale, result.scaling.col_scale);
    if (!result.error.ok()) return result;
    result.scaling.curtis_reid_sweeps = options.curtis_reid_sweeps;
    result.scaling.stages.push_back(ScalingStage::kCurtisReid);
  }
  if (options.mode == ScalingMode::kCurtisReid) return result;

  if ((options.legacy_scaling_mode & kPdlpScalingRuiz) != 0) {
    const ScalingResult ruiz =
        computeRuizScaling(working, options.ruiz_iterations);
    result.error = appendScaling(ruiz, working, result.scaling);
    if (!result.error.ok()) return result;
    result.scaling.ruiz_iterations = options.ruiz_iterations;
  }
  if ((options.legacy_scaling_mode & kPdlpScalingPC) != 0) {
    result.error = appendScaling(computePockChambolleScaling(working), working,
                                 result.scaling);
    if (!result.error.ok()) return result;
  }
  if ((options.legacy_scaling_mode & kPdlpScalingL2) != 0) {
    result.error =
        appendScaling(computeL2Scaling(working), working, result.scaling);
    if (!result.error.ok()) return result;
  }
  return result;
}

api::Error applyScalingPipeline(Problem& problem,
                                const ScalingOptions& options) {
  if (!problem.initialized())
    return invalidScaling("HiPDLP problem is not initialized");
  const ScalingResult result =
      computeScalingPipeline(problem.matrix(), options);
  if (!result.ok()) return result.error;
  if (result.scaling.stages.empty()) return api::Error();
  return problem.applyDiagonalScaling(result.scaling.row_scale,
                                      result.scaling.col_scale, result.scaling);
}

api::Error scalePrimal(const Problem& problem,
                       const std::vector<double>& unscaled,
                       std::vector<double>& scaled) {
  return transformVector(problem, unscaled, scaled, RoundTripKind::kPrimal,
                         true);
}

api::Error unscalePrimal(const Problem& problem,
                         const std::vector<double>& scaled,
                         std::vector<double>& unscaled) {
  return transformVector(problem, scaled, unscaled, RoundTripKind::kPrimal,
                         false);
}

api::Error scaleDual(const Problem& problem,
                     const std::vector<double>& unscaled,
                     std::vector<double>& scaled) {
  return transformVector(problem, unscaled, scaled, RoundTripKind::kDual, true);
}

api::Error unscaleDual(const Problem& problem,
                       const std::vector<double>& scaled,
                       std::vector<double>& unscaled) {
  return transformVector(problem, scaled, unscaled, RoundTripKind::kDual,
                         false);
}

api::Error scaleReducedCosts(const Problem& problem,
                             const std::vector<double>& unscaled,
                             std::vector<double>& scaled) {
  return transformVector(problem, unscaled, scaled, RoundTripKind::kReducedCost,
                         true);
}

api::Error unscaleReducedCosts(const Problem& problem,
                               const std::vector<double>& scaled,
                               std::vector<double>& unscaled) {
  return transformVector(problem, scaled, unscaled, RoundTripKind::kReducedCost,
                         false);
}

}  // namespace model
}  // namespace hipdlp
