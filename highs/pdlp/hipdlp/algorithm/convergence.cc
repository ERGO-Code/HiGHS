#include "pdlp/hipdlp/algorithm/convergence.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace algorithm {
namespace {

bool finiteNonnegative(const double value) {
  return std::isfinite(value) && value >= 0.0;
}

bool finitePositive(const double value) {
  return std::isfinite(value) && value > 0.0;
}

StoppingResult failure(const api::ErrorCode code, const char* message) {
  StoppingResult result;
  result.error = api::Error(code, message);
  return result;
}

void updateScaledSumOfSquares(const double value, double& scale,
                              double& sum_of_squares) {
  const double magnitude = std::fabs(value);
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

double scaledNorm(const double scale, const double sum_of_squares) {
  return scale == 0.0 ? 0.0 : scale * std::sqrt(sum_of_squares);
}

}  // namespace

RawConvergenceMetrics::RawConvergenceMetrics()
    : primal_residual(0.0),
      dual_residual(0.0),
      objective_gap(0.0),
      primal_objective(0.0),
      dual_objective(0.0),
      unscaled_rhs_norm(0.0),
      unscaled_cost_norm(0.0) {}

ConvergenceTolerances::ConvergenceTolerances()
    : primal(std::numeric_limits<double>::quiet_NaN()),
      dual(std::numeric_limits<double>::quiet_NaN()),
      gap(std::numeric_limits<double>::quiet_NaN()) {}

StoppingLimits::StoppingLimits() : iteration(kHighsIInf), time(kHighsInf) {}

StoppingInput::StoppingInput() : iteration(0), elapsed_time(0.0) {}

NormalizedConvergenceMetrics::NormalizedConvergenceMetrics()
    : relative_primal(0.0), relative_dual(0.0), relative_gap(0.0) {}

StoppingResult::StoppingResult()
    : decision(StoppingDecision::kNumericalError), metrics_valid(false) {}

bool StoppingResult::ok() const {
  return error.ok() && metrics_valid &&
         decision != StoppingDecision::kNumericalError;
}

api::Error computeRawConvergenceMetrics(
    const model::Problem& problem, const std::vector<double>& primal,
    const std::vector<double>& dual, const std::vector<double>& ax,
    const std::vector<double>& aty, RawConvergenceMetrics& raw) {
  const model::ScalingMetadata& scaling = problem.scaling();
  const std::vector<double>& row_lower = problem.rowLower();
  const std::vector<double>& col_lower = problem.colLower();
  const std::vector<double>& col_upper = problem.colUpper();
  const std::vector<double>& costs = problem.costs();
  const std::vector<HighsBool>& equality_rows = problem.equalityRows();

  double primal_objective = problem.objectiveOffset();
  double dual_objective = problem.objectiveOffset();
  double primal_scale = 0.0;
  double primal_sum_of_squares = 1.0;
  for (HighsInt row = 0; row < problem.numRows(); ++row) {
    double violation = ax[row] - row_lower[row];
    if (!equality_rows[row]) violation = std::min(0.0, violation);
    updateScaledSumOfSquares(violation / scaling.row_scale[row], primal_scale,
                             primal_sum_of_squares);
    dual_objective += row_lower[row] * dual[row];
  }

  double dual_scale = 0.0;
  double dual_sum_of_squares = 1.0;
  for (HighsInt col = 0; col < problem.numCols(); ++col) {
    const double reduced_cost = costs[col] - aty[col];
    const bool has_lower = col_lower[col] > -kHighsInf;
    const bool has_upper = col_upper[col] < kHighsInf;
    const double lower_slack = has_lower ? std::max(0.0, reduced_cost) : 0.0;
    const double upper_slack = has_upper ? std::max(0.0, -reduced_cost) : 0.0;
    updateScaledSumOfSquares(
        (reduced_cost - lower_slack + upper_slack) / scaling.col_scale[col],
        dual_scale, dual_sum_of_squares);
    primal_objective += costs[col] * primal[col];
    if (has_lower) dual_objective += col_lower[col] * lower_slack;
    if (has_upper) dual_objective -= col_upper[col] * upper_slack;
  }

  raw.primal_residual = scaledNorm(primal_scale, primal_sum_of_squares);
  raw.dual_residual = scaledNorm(dual_scale, dual_sum_of_squares);
  raw.primal_objective = primal_objective;
  raw.dual_objective = dual_objective;
  raw.objective_gap = std::fabs(primal_objective - dual_objective);
  raw.unscaled_rhs_norm = problem.unscaledNorms().row_bounds_l2;
  raw.unscaled_cost_norm = problem.unscaledNorms().objective_l2;
  if (!std::isfinite(raw.primal_residual) ||
      !std::isfinite(raw.dual_residual) ||
      !std::isfinite(raw.primal_objective) ||
      !std::isfinite(raw.dual_objective) || !std::isfinite(raw.objective_gap)) {
    return api::Error(api::ErrorCode::kInvalidArgument,
                      "fixed trace convergence metrics are non-finite");
  }
  return api::Error();
}

StoppingResult evaluateStopping(const StoppingInput& input) {
  if (!finitePositive(input.tolerances.primal) ||
      !finitePositive(input.tolerances.dual) ||
      !finitePositive(input.tolerances.gap)) {
    return failure(api::ErrorCode::kInvalidOption,
                   "HiPDLP convergence tolerances must be finite and positive");
  }

  if (!finiteNonnegative(input.raw.primal_residual) ||
      !finiteNonnegative(input.raw.dual_residual) ||
      !finiteNonnegative(input.raw.objective_gap) ||
      !finiteNonnegative(input.raw.unscaled_rhs_norm) ||
      !finiteNonnegative(input.raw.unscaled_cost_norm) ||
      !std::isfinite(input.raw.primal_objective) ||
      !std::isfinite(input.raw.dual_objective)) {
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP convergence metrics must be finite and nonnegative "
                   "where required");
  }

  if (input.iteration < 0 || !finiteNonnegative(input.elapsed_time)) {
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP stopping progress is invalid");
  }
  if (input.limits.iteration < 0 || (input.limits.time != kHighsInf &&
                                     !finiteNonnegative(input.limits.time))) {
    return failure(api::ErrorCode::kInvalidOption,
                   "HiPDLP stopping limits are invalid");
  }

  const double primal_denominator = 1.0 + input.raw.unscaled_rhs_norm;
  const double dual_denominator = 1.0 + input.raw.unscaled_cost_norm;
  const double gap_denominator = 1.0 + std::fabs(input.raw.primal_objective) +
                                 std::fabs(input.raw.dual_objective);
  if (!std::isfinite(primal_denominator) || !std::isfinite(dual_denominator) ||
      !std::isfinite(gap_denominator)) {
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP convergence normalization overflowed");
  }

  StoppingResult result;
  result.metrics.relative_primal =
      input.raw.primal_residual / primal_denominator;
  result.metrics.relative_dual = input.raw.dual_residual / dual_denominator;
  result.metrics.relative_gap = input.raw.objective_gap / gap_denominator;
  if (!std::isfinite(result.metrics.relative_primal) ||
      !std::isfinite(result.metrics.relative_dual) ||
      !std::isfinite(result.metrics.relative_gap)) {
    return failure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP convergence normalization produced a non-finite metric");
  }
  result.metrics_valid = true;

  if (result.metrics.relative_primal <= input.tolerances.primal &&
      result.metrics.relative_dual <= input.tolerances.dual &&
      result.metrics.relative_gap <= input.tolerances.gap) {
    result.decision = StoppingDecision::kOptimal;
  } else if (input.limits.iteration != kHighsIInf &&
             input.iteration >= input.limits.iteration) {
    result.decision = StoppingDecision::kIterationLimit;
  } else if (input.limits.time != kHighsInf &&
             input.elapsed_time >= input.limits.time) {
    result.decision = StoppingDecision::kTimeLimit;
  } else {
    result.decision = StoppingDecision::kContinue;
  }
  return result;
}

}  // namespace algorithm
}  // namespace hipdlp