#include "pdlp/hipdlp/algorithm/fixed_point_error.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "pdlp/hipdlp/backend/backend.hpp"

namespace hipdlp {
namespace algorithm {
namespace {

const std::size_t kScalarReadbacksPerEvaluation = 3;

bool finitePositive(const double value) {
  return std::isfinite(value) && value > 0.0;
}

bool finiteNonnegative(const double value) {
  return std::isfinite(value) && value >= 0.0;
}

FixedPointErrorResult failure(const api::Error& error) {
  FixedPointErrorResult result;
  result.error = error;
  return result;
}

api::Error invalidArgument(const char* message) {
  return api::Error(api::ErrorCode::kInvalidArgument, message);
}

bool validDistinctHandles(const FixedPointErrorInput& input) {
  const backend::VectorHandle* handles[] = {
      input.next_primal,    input.next_dual,    input.reflected_primal,
      input.reflected_dual, input.delta_primal, input.delta_dual,
      input.aty_delta,
  };
  const std::size_t handle_count = sizeof(handles) / sizeof(handles[0]);
  for (std::size_t index = 0; index < handle_count; ++index) {
    if (handles[index] == nullptr || !handles[index]->valid()) return false;
    for (std::size_t previous = 0; previous < index; ++previous) {
      if (handles[index] == handles[previous]) return false;
    }
  }
  return true;
}

}  // namespace

FixedPointErrorInput::FixedPointErrorInput()
    : next_primal(nullptr),
      next_dual(nullptr),
      reflected_primal(nullptr),
      reflected_dual(nullptr),
      omega(std::numeric_limits<double>::quiet_NaN()),
      eta(std::numeric_limits<double>::quiet_NaN()),
      delta_primal(nullptr),
      delta_dual(nullptr),
      aty_delta(nullptr) {}

FixedPointErrorDiagnostics::FixedPointErrorDiagnostics()
    : device_to_host_scalar_readbacks(0) {}

FixedPointErrorResult::FixedPointErrorResult()
    : computed(false), fixed_point_error(0.0) {}

bool FixedPointErrorResult::ok() const {
  return error.ok() && computed && finiteNonnegative(fixed_point_error);
}

api::Error prepareFixedPointError(
    backend::Backend& backend, const backend::ProblemHandle& problem,
    const FixedPointErrorInput& input) {
  if (!problem.valid() || !validDistinctHandles(input)) {
    return invalidArgument(
        "HiPDLP fixed-point error handles must be valid and distinct");
  }
  if (!finitePositive(input.omega) || !finitePositive(input.eta)) {
    return invalidArgument(
        "HiPDLP fixed-point error omega and eta must be finite and positive");
  }

  api::Error error = backend.copy(*input.next_primal, *input.delta_primal);
  if (!error.ok()) return error;
  error = backend.axpy(-1.0, *input.reflected_primal, *input.delta_primal);
  if (!error.ok()) return error;
  error = backend.copy(*input.next_dual, *input.delta_dual);
  if (!error.ok()) return error;
  error = backend.axpy(-1.0, *input.reflected_dual, *input.delta_dual);
  if (!error.ok()) return error;
  return backend.multiplyTranspose(problem, *input.delta_dual, *input.aty_delta);
}

FixedPointErrorResult evaluateFixedPointError(
    backend::Backend& backend, const backend::ProblemHandle& problem,
    const FixedPointErrorInput& input) {
  api::Error error = prepareFixedPointError(backend, problem, input);
  if (!error.ok()) return failure(error);

  double primal_norm = 0.0;
  error = backend.norm2(*input.delta_primal, primal_norm);
  if (!error.ok()) return failure(error);
  if (!finiteNonnegative(primal_norm)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point primal displacement norm is invalid"));
  }

  double dual_norm = 0.0;
  error = backend.norm2(*input.delta_dual, dual_norm);
  if (!error.ok()) return failure(error);
  if (!finiteNonnegative(dual_norm)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point dual displacement norm is invalid"));
  }

  double cross_term = 0.0;
  error = backend.dot(*input.delta_primal, *input.aty_delta, cross_term);
  if (!error.ok()) return failure(error);
  FixedPointErrorResult result = computeFixedPointError(
      primal_norm, dual_norm, cross_term, input.omega, input.eta);
  if (!result.ok()) return result;
  result.diagnostics.device_to_host_scalar_readbacks =
      kScalarReadbacksPerEvaluation;
  return result;
}

FixedPointErrorResult computeFixedPointError(
    const double primal_norm, const double dual_norm, const double cross_term,
    const double omega, const double eta) {
  if (!finitePositive(omega) || !finitePositive(eta)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point error omega and eta must be finite and positive"));
  }
  if (!finiteNonnegative(primal_norm)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point primal displacement norm is invalid"));
  }
  if (!finiteNonnegative(dual_norm)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point dual displacement norm is invalid"));
  }
  if (!std::isfinite(cross_term)) {
    return failure(
        invalidArgument("HiPDLP fixed-point interaction dot is invalid"));
  }

  const double primal_norm_squared = primal_norm * primal_norm;
  const double dual_norm_squared = dual_norm * dual_norm;
  if (!finiteNonnegative(primal_norm_squared) ||
      !finiteNonnegative(dual_norm_squared)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point displacement norm square is outside FP64 range"));
  }

  const double weighted_primal = primal_norm_squared * omega;
  const double weighted_dual = dual_norm_squared / omega;
  if (!finiteNonnegative(weighted_primal) ||
      !finiteNonnegative(weighted_dual)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point weighted movement is outside FP64 range"));
  }
  const double movement = weighted_primal + weighted_dual;
  if (!finiteNonnegative(movement)) {
    return failure(
        invalidArgument("HiPDLP fixed-point movement is outside FP64 range"));
  }

  const double twice_eta = 2.0 * eta;
  if (!finitePositive(twice_eta)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point interaction factor is outside FP64 range"));
  }
  const double interaction = twice_eta * cross_term;
  if (!std::isfinite(interaction)) {
    return failure(invalidArgument(
        "HiPDLP fixed-point interaction is outside FP64 range"));
  }

  const double combined = movement + interaction;
  if (!std::isfinite(combined)) {
    return failure(
        invalidArgument("HiPDLP fixed-point error sum is outside FP64 range"));
  }
  const double clamped = std::max(0.0, combined);
  if (!finiteNonnegative(clamped)) {
    return failure(
        invalidArgument("HiPDLP fixed-point error clamp is invalid"));
  }
  const double fixed_point_error = std::sqrt(clamped);
  if (!finiteNonnegative(fixed_point_error)) {
    return failure(
        invalidArgument("HiPDLP fixed-point square root is invalid"));
  }

  FixedPointErrorResult result;
  result.computed = true;
  result.fixed_point_error = fixed_point_error;
  return result;
}

}  // namespace algorithm
}  // namespace hipdlp