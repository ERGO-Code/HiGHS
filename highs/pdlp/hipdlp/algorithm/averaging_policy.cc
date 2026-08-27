#include "pdlp/hipdlp/algorithm/averaging_policy.hpp"

#include <cmath>
#include <limits>

namespace hipdlp {
namespace algorithm {
namespace {

bool validState(const AveragingState& state) {
  if (!state.history_valid || state.sample_count < 0 ||
      !std::isfinite(state.total_weight))
    return false;
  if (state.sample_count == 0)
    return state.total_weight == 0.0 && !std::signbit(state.total_weight);
  return state.total_weight > 0.0;
}

AveragingResult failure(const char* message) {
  AveragingResult result;
  result.error = api::Error(api::ErrorCode::kInvalidArgument, message);
  return result;
}

}  // namespace

AveragingState::AveragingState()
    : sample_count(0), total_weight(0.0), history_valid(true) {}

bool AveragingState::hasAverage() const {
  return history_valid && sample_count > 0 && std::isfinite(total_weight) &&
         total_weight > 0.0;
}

AveragingInput::AveragingInput()
    : primal_step_size(std::numeric_limits<double>::quiet_NaN()),
      dual_step_size(std::numeric_limits<double>::quiet_NaN()),
      completed_halpern_iterate(false) {}

AveragingResult::AveragingResult() : sample_weight(0.0), normalization(0.0) {
  next.history_valid = false;
}

bool AveragingResult::ok() const {
  return error.ok() && next.hasAverage() && std::isfinite(sample_weight) &&
         sample_weight > 0.0 && std::isfinite(normalization) &&
         normalization > 0.0;
}

AveragingResult updateAveraging(const AveragingInput& input) {
  if (!input.completed_halpern_iterate)
    return failure("HiPDLP averaging requires a completed Halpern iterate");
  if (!validState(input.previous))
    return failure("HiPDLP averaging history is invalid");
  if (!std::isfinite(input.primal_step_size) || input.primal_step_size <= 0.0 ||
      !std::isfinite(input.dual_step_size) || input.dual_step_size <= 0.0)
    return failure("HiPDLP averaging step sizes must be finite and positive");
  if (input.previous.sample_count == std::numeric_limits<HighsInt>::max())
    return failure("HiPDLP averaging sample count overflowed");

  const double product = input.primal_step_size * input.dual_step_size;
  if (!std::isfinite(product) || product <= 0.0)
    return failure("HiPDLP averaging step-size product is outside FP64 range");
  const double sample_weight = std::sqrt(product);
  if (!std::isfinite(sample_weight) || sample_weight <= 0.0)
    return failure("HiPDLP averaging sample weight is invalid");
  if (sample_weight >
      std::numeric_limits<double>::max() - input.previous.total_weight)
    return failure("HiPDLP averaging total weight overflowed");

  const double next_total_weight = input.previous.total_weight + sample_weight;
  const double normalization = 1.0 / next_total_weight;
  if (!std::isfinite(next_total_weight) || next_total_weight <= 0.0 ||
      !std::isfinite(normalization) || normalization <= 0.0)
    return failure("HiPDLP averaging normalization is invalid");

  AveragingResult result;
  result.next.sample_count = input.previous.sample_count + 1;
  result.next.total_weight = next_total_weight;
  result.next.history_valid = true;
  result.sample_weight = sample_weight;
  result.normalization = normalization;
  return result;
}

AveragingState resetAveraging() { return AveragingState(); }

}  // namespace algorithm
}  // namespace hipdlp