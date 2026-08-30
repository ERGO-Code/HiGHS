#include "pdlp/hipdlp/algorithm/step_size_policy.hpp"

#include <cmath>
#include <limits>

namespace hipdlp {
namespace algorithm {
namespace {

// Strict legacy guard windows from PDLPSolver::updatePrimalWeightAtRestart().
const double kMinimumDistance = 1e-16;
const double kMaximumDistance = 1e12;
const double kMinimumResidualRatio = 1e-8;
const double kMaximumResidualRatio = 1e8;
const double kUndefinedResidualRatio = 1e300;

bool finitePositive(const double value) {
  return std::isfinite(value) && value > 0.0;
}

bool finiteNonnegative(const double value) {
  return std::isfinite(value) && value >= 0.0;
}

StepSizeInitializationResult failure(const api::ErrorCode code,
                                     const char* message) {
  StepSizeInitializationResult result;
  result.error = api::Error(code, message);
  return result;
}

bool validPrimalWeightParameters(const PrimalWeightParameters& parameters) {
  return finiteNonnegative(parameters.integral_smoothing) &&
         parameters.integral_smoothing <= 1.0 &&
         finiteNonnegative(parameters.k_p) &&
         finiteNonnegative(parameters.k_i) && finiteNonnegative(parameters.k_d);
}

bool validPrimalWeightState(const PrimalWeightState& state,
                            const bool require_legacy_gap_flag) {
  if (!state.history_valid) return false;
  if (!finitePositive(state.primal_weight) ||
      !finitePositive(state.best_primal_weight))
    return false;
  if (!std::isfinite(state.error_sum) || !std::isfinite(state.previous_error) ||
      !finiteNonnegative(state.best_residual_gap))
    return false;
  if (require_legacy_gap_flag && !state.best_gap_valid)
    return state.best_residual_gap == 0.0 &&
           !std::signbit(state.best_residual_gap);
  return true;
}

PrimalWeightUpdateResult primalWeightFailure(const api::ErrorCode code,
                                             const char* message) {
  PrimalWeightUpdateResult result;
  result.error = api::Error(code, message);
  return result;
}

bool multiplicationInRange(const double left, const double right,
                           const double product) {
  return std::isfinite(product) &&
         (product != 0.0 || left == 0.0 || right == 0.0);
}

}  // namespace

StepSizeInitializationParameters::StepSizeInitializationParameters()
    : safety_factor(0.998) {}

StepSizeInitializationInput::StepSizeInitializationInput()
    : operator_norm_squared(std::numeric_limits<double>::quiet_NaN()),
      unscaled_rhs_norm(std::numeric_limits<double>::quiet_NaN()),
      unscaled_cost_norm(std::numeric_limits<double>::quiet_NaN()) {}

InitialStepSizes::InitialStepSizes()
    : primal_weight(0.0), base_step(0.0), primal_step(0.0), dual_step(0.0) {}

StepSizeInitializationResult::StepSizeInitializationResult() {}

bool StepSizeInitializationResult::ok() const {
  return error.ok() && finitePositive(step_sizes.primal_weight) &&
         finitePositive(step_sizes.base_step) &&
         finitePositive(step_sizes.primal_step) &&
         finitePositive(step_sizes.dual_step);
}

StepSizeInitializationResult initializeStepSizes(
    const StepSizeInitializationInput& input) {
  if (!finitePositive(input.parameters.safety_factor) ||
      input.parameters.safety_factor > 1.0)
    return failure(
        api::ErrorCode::kInvalidOption,
        "HiPDLP step-size safety factor must lie in (0, 1] and be finite");
  if (!finitePositive(input.operator_norm_squared))
    return failure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP step-size operator norm squared must be finite and positive");
  if (!finiteNonnegative(input.unscaled_rhs_norm) ||
      !finiteNonnegative(input.unscaled_cost_norm))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP step-size unscaled norms must be finite and "
                   "nonnegative");

  const double primal_denominator = 1.0 + input.unscaled_rhs_norm;
  const double dual_denominator = 1.0 + input.unscaled_cost_norm;
  if (!finitePositive(primal_denominator) || !finitePositive(dual_denominator))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP step-size norm offsets are outside FP64 range");

  const double primal_weight = dual_denominator / primal_denominator;
  if (!finitePositive(primal_weight))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP step-size primal weight is outside FP64 range");

  const double operator_norm = std::sqrt(input.operator_norm_squared);
  if (!finitePositive(operator_norm))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP step-size operator norm is outside FP64 range");

  const double base_step = input.parameters.safety_factor / operator_norm;
  if (!finitePositive(base_step))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP step-size base step is outside FP64 range");

  const double primal_step = base_step / primal_weight;
  const double dual_step = base_step * primal_weight;
  if (!finitePositive(primal_step) || !finitePositive(dual_step))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP initial primal/dual steps are outside FP64 range");

  StepSizeInitializationResult result;
  result.step_sizes.primal_weight = primal_weight;
  result.step_sizes.base_step = base_step;
  result.step_sizes.primal_step = primal_step;
  result.step_sizes.dual_step = dual_step;
  return result;
}

PrimalWeightParameters::PrimalWeightParameters()
    : integral_smoothing(0.3), k_p(0.99), k_i(0.01), k_d(0.0) {}

PrimalWeightState::PrimalWeightState()
    : primal_weight(0.0),
      error_sum(0.0),
      previous_error(0.0),
      best_primal_weight(0.0),
      best_residual_gap(0.0),
      best_gap_valid(false),
      history_valid(false) {}

PrimalWeightInput::PrimalWeightInput()
    : primal_step(std::numeric_limits<double>::quiet_NaN()),
      dual_step(std::numeric_limits<double>::quiet_NaN()),
      primal_distance(std::numeric_limits<double>::quiet_NaN()),
      dual_distance(std::numeric_limits<double>::quiet_NaN()),
      relative_primal_residual(std::numeric_limits<double>::quiet_NaN()),
      relative_dual_residual(std::numeric_limits<double>::quiet_NaN()) {}

PrimalWeightUpdateInput::PrimalWeightUpdateInput()
    : primal_anchor_distance(std::numeric_limits<double>::quiet_NaN()),
      dual_anchor_distance(std::numeric_limits<double>::quiet_NaN()),
      relative_primal_residual(std::numeric_limits<double>::quiet_NaN()),
      relative_dual_residual(std::numeric_limits<double>::quiet_NaN()),
      primal_step_size(std::numeric_limits<double>::quiet_NaN()),
      dual_step_size(std::numeric_limits<double>::quiet_NaN()) {}

PrimalWeightUpdateResult::PrimalWeightUpdateResult()
    : kind(PrimalWeightUpdateKind::kNone),
      outcome(PrimalWeightOutcome::kNone),
      base_step(0.0),
      primal_step(0.0),
      dual_step(0.0) {}

bool PrimalWeightUpdateResult::ok() const {
  return error.ok() && kind != PrimalWeightUpdateKind::kNone &&
         outcome != PrimalWeightOutcome::kNone &&
         validPrimalWeightState(next, false) &&
         finitePositive(step_sizes.primal_weight) &&
         finitePositive(step_sizes.base_step) &&
         finitePositive(step_sizes.primal_step) &&
         finitePositive(step_sizes.dual_step) &&
         step_sizes.primal_weight == next.primal_weight &&
         base_step == step_sizes.base_step &&
         primal_step == step_sizes.primal_step &&
         dual_step == step_sizes.dual_step;
}

namespace {

PrimalWeightUpdateResult updatePrimalWeightImpl(
    const PrimalWeightUpdateInput& input,
    const bool preserve_legacy_missing_best_gap) {
  if (!validPrimalWeightParameters(input.parameters))
    return primalWeightFailure(
        api::ErrorCode::kInvalidOption,
        "HiPDLP primal-weight PID parameters are invalid");
  if (!validPrimalWeightState(input.previous, preserve_legacy_missing_best_gap))
    return primalWeightFailure(api::ErrorCode::kInvalidArgument,
                               "HiPDLP primal-weight state is invalid");
  if (!finiteNonnegative(input.primal_anchor_distance) ||
      !finiteNonnegative(input.dual_anchor_distance))
    return primalWeightFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP primal-weight anchor distances must be finite and "
        "nonnegative");
  if (!finiteNonnegative(input.relative_primal_residual) ||
      !finiteNonnegative(input.relative_dual_residual))
    return primalWeightFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP primal-weight relative residuals must be finite and "
        "nonnegative");
  if (!finitePositive(input.primal_step_size) ||
      !finitePositive(input.dual_step_size))
    return primalWeightFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP primal-weight step sizes must be finite and positive");

  const bool residual_ratio_defined = input.relative_primal_residual > 0.0;
  const double residual_ratio =
      residual_ratio_defined
          ? input.relative_dual_residual / input.relative_primal_residual
          : kUndefinedResidualRatio;
  if (!std::isfinite(residual_ratio) ||
      (input.relative_dual_residual > 0.0 && residual_ratio == 0.0))
    return primalWeightFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP primal-weight residual ratio is outside FP64 range");

  PrimalWeightState next = input.previous;
  PrimalWeightUpdateKind kind = PrimalWeightUpdateKind::kRevertedToBest;
  if (input.primal_anchor_distance > kMinimumDistance &&
      input.dual_anchor_distance > kMinimumDistance &&
      input.primal_anchor_distance < kMaximumDistance &&
      input.dual_anchor_distance < kMaximumDistance &&
      residual_ratio > kMinimumResidualRatio &&
      residual_ratio < kMaximumResidualRatio) {
    const double dual_log = std::log(input.dual_anchor_distance);
    const double primal_log = std::log(input.primal_anchor_distance);
    const double weight_log = std::log(input.previous.primal_weight);
    if (!std::isfinite(dual_log) || !std::isfinite(primal_log) ||
        !std::isfinite(weight_log))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight logarithms are outside FP64 range");

    const double distance_log_ratio = dual_log - primal_log;
    const double error = distance_log_ratio - weight_log;
    if (!std::isfinite(distance_log_ratio) || !std::isfinite(error))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight PID error is outside FP64 range");

    const double smoothed_error_sum =
        input.parameters.integral_smoothing * input.previous.error_sum;
    if (!multiplicationInRange(input.parameters.integral_smoothing,
                               input.previous.error_sum, smoothed_error_sum))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight smoothed integral is outside FP64 range");
    const double error_sum = smoothed_error_sum + error;
    const double delta_error = error - input.previous.previous_error;
    if (!std::isfinite(error_sum) || !std::isfinite(delta_error))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight PID error terms are outside FP64 range");

    const double proportional_term = input.parameters.k_p * error;
    if (!multiplicationInRange(input.parameters.k_p, error, proportional_term))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight PID proportional term is outside FP64 range");
    const double integral_term = input.parameters.k_i * error_sum;
    if (!multiplicationInRange(input.parameters.k_i, error_sum, integral_term))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight PID integral term is outside FP64 range");
    const double derivative_term = input.parameters.k_d * delta_error;
    if (!multiplicationInRange(input.parameters.k_d, delta_error,
                               derivative_term))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight PID derivative term is outside FP64 range");

    const double proportional_integral = proportional_term + integral_term;
    const double exponent = proportional_integral + derivative_term;
    if (!std::isfinite(proportional_integral) || !std::isfinite(exponent))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight PID exponent is outside FP64 range");

    const double factor = std::exp(exponent);
    if (!finitePositive(factor))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight PID factor overflowed or underflowed");

    const double primal_weight = input.previous.primal_weight * factor;
    if (!finitePositive(primal_weight))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP updated primal weight is outside FP64 range");

    next.primal_weight = primal_weight;
    next.error_sum = error_sum;
    next.previous_error = error;
    kind = PrimalWeightUpdateKind::kControlled;
  } else {
    next.primal_weight = input.previous.best_primal_weight;
    next.error_sum = 0.0;
    next.previous_error = 0.0;
  }

  bool residual_gap_computed = false;
  double residual_gap = next.best_residual_gap;
  if (residual_ratio_defined && input.relative_dual_residual > 0.0 &&
      residual_ratio > 0.0) {
    residual_gap = std::fabs(std::log10(residual_ratio));
    if (!finiteNonnegative(residual_gap))
      return primalWeightFailure(
          api::ErrorCode::kInvalidArgument,
          "HiPDLP primal-weight residual balance is outside FP64 range");
    residual_gap_computed = true;
  }
  const bool strict_improvement = residual_gap < next.best_residual_gap;
  if (residual_gap_computed &&
      (strict_improvement ||
       (preserve_legacy_missing_best_gap && !next.best_gap_valid))) {
    next.best_residual_gap = residual_gap;
    next.best_gap_valid = true;
    next.best_primal_weight = next.primal_weight;
  }

  const double step_product = input.primal_step_size * input.dual_step_size;
  if (!finitePositive(step_product))
    return primalWeightFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP primal-weight step-size product is outside FP64 range");
  const double base_step = std::sqrt(step_product);
  if (!finitePositive(base_step))
    return primalWeightFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP primal-weight geometric mean step is outside FP64 range");

  const double primal_step = base_step / next.primal_weight;
  const double dual_step = base_step * next.primal_weight;
  if (!finitePositive(primal_step) || !finitePositive(dual_step))
    return primalWeightFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP primal-weight derived steps are outside FP64 range");

  PrimalWeightUpdateResult result;
  result.next = next;
  result.kind = kind;
  result.step_sizes.primal_weight = next.primal_weight;
  result.step_sizes.base_step = base_step;
  result.step_sizes.primal_step = primal_step;
  result.step_sizes.dual_step = dual_step;
  result.outcome = kind == PrimalWeightUpdateKind::kControlled
                       ? PrimalWeightOutcome::kPidUpdate
                       : PrimalWeightOutcome::kFallback;
  result.base_step = base_step;
  result.primal_step = primal_step;
  result.dual_step = dual_step;
  return result;
}

}  // namespace

PrimalWeightUpdateResult updatePrimalWeightAtRestart(
    const PrimalWeightUpdateInput& input) {
  return updatePrimalWeightImpl(input, false);
}

PrimalWeightUpdateResult updatePrimalWeight(const PrimalWeightInput& input) {
  PrimalWeightUpdateInput adapted;
  adapted.previous = input.previous;
  adapted.primal_anchor_distance = input.primal_distance;
  adapted.dual_anchor_distance = input.dual_distance;
  adapted.relative_primal_residual = input.relative_primal_residual;
  adapted.relative_dual_residual = input.relative_dual_residual;
  adapted.primal_step_size = input.primal_step;
  adapted.dual_step_size = input.dual_step;
  adapted.parameters = input.parameters;
  return updatePrimalWeightImpl(adapted, true);
}

PrimalWeightState resetPrimalWeight(const double primal_weight) {
  PrimalWeightState state;
  if (!finitePositive(primal_weight)) return state;
  state.primal_weight = primal_weight;
  state.best_primal_weight = primal_weight;
  state.history_valid = true;
  return state;
}

}  // namespace algorithm
}  // namespace hipdlp
