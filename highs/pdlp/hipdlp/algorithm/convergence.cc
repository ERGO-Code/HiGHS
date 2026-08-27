#include "pdlp/hipdlp/algorithm/convergence.hpp"

#include <cmath>
#include <limits>

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