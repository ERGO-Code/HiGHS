#include "pdlp/hipdlp/algorithm/restart_policy.hpp"

#include <cmath>
#include <limits>

namespace hipdlp {
namespace algorithm {
namespace {

bool finitePositive(const double value) {
  return std::isfinite(value) && value > 0.0;
}

bool validMetrics(const NormalizedConvergenceMetrics& metrics) {
  return std::isfinite(metrics.relative_primal) &&
         metrics.relative_primal >= 0.0 &&
         std::isfinite(metrics.relative_dual) && metrics.relative_dual >= 0.0 &&
         std::isfinite(metrics.relative_gap) && metrics.relative_gap >= 0.0;
}

bool zeroMetrics(const NormalizedConvergenceMetrics& metrics) {
  return metrics.relative_primal == 0.0 && metrics.relative_dual == 0.0 &&
         metrics.relative_gap == 0.0;
}

bool validAveragingState(const AveragingState& averaging) {
  if (!averaging.history_valid || averaging.sample_count < 0 ||
      !std::isfinite(averaging.total_weight))
    return false;
  if (averaging.sample_count == 0)
    return averaging.total_weight == 0.0 &&
           !std::signbit(averaging.total_weight);
  return averaging.total_weight > 0.0;
}

bool validRestartState(const RestartState& state) {
  if (!state.history_valid || state.restart_count < 0) return false;
  if (!state.initialized)
    return state.restart_count == 0 && zeroMetrics(state.restart_metrics) &&
           zeroMetrics(state.last_candidate_metrics);
  return validMetrics(state.restart_metrics) &&
         validMetrics(state.last_candidate_metrics);
}

bool validParameters(const RestartParameters& parameters) {
  return finitePositive(parameters.primal_weight) &&
         std::isfinite(parameters.sufficient_decay) &&
         parameters.sufficient_decay > 0.0 &&
         std::isfinite(parameters.necessary_decay) &&
         parameters.necessary_decay > parameters.sufficient_decay &&
         parameters.necessary_decay < 1.0 &&
         std::isfinite(parameters.artificial_restart_threshold) &&
         parameters.artificial_restart_threshold > 0.0 &&
         parameters.artificial_restart_threshold <= 1.0;
}

bool computeScore(const NormalizedConvergenceMetrics& metrics,
                  const double primal_weight, double& score) {
  const double weighted_primal = primal_weight * metrics.relative_primal;
  const double weighted_dual = metrics.relative_dual / primal_weight;
  if (!std::isfinite(weighted_primal) || !std::isfinite(weighted_dual))
    return false;
  score = std::hypot(std::hypot(weighted_primal, weighted_dual),
                     metrics.relative_gap);
  return std::isfinite(score);
}

RestartResult failure(const api::ErrorCode code, const char* message) {
  RestartResult result;
  result.error = api::Error(code, message);
  return result;
}

bool validOptionalError(const double value, const bool valid) {
  if (!valid) return value == 0.0 && !std::signbit(value);
  return std::isfinite(value) && value >= 0.0;
}

bool validFixedPointRestartState(const FixedPointRestartState& state) {
  return state.history_valid &&
         validOptionalError(state.reference_error, state.reference_valid) &&
         validOptionalError(state.previous_trial_error,
                            state.previous_trial_valid);
}

bool validFixedPointRestartParameters(
    const FixedPointRestartParameters& parameters) {
  return std::isfinite(parameters.sufficient_decay) &&
         parameters.sufficient_decay > 0.0 &&
         std::isfinite(parameters.necessary_decay) &&
         parameters.necessary_decay > parameters.sufficient_decay &&
         parameters.necessary_decay < 1.0 &&
         std::isfinite(parameters.artificial_threshold) &&
         parameters.artificial_threshold > 0.0 && parameters.check_interval > 0;
}

FixedPointRestartResult fixedPointFailure(const api::ErrorCode code,
                                          const char* message) {
  FixedPointRestartResult result;
  result.error = api::Error(code, message);
  return result;
}

}  // namespace

RestartParameters::RestartParameters()
    : primal_weight(1.0),
      sufficient_decay(0.2),
      necessary_decay(0.8),
      artificial_restart_threshold(0.36) {}

RestartState::RestartState()
    : restart_count(0), initialized(false), history_valid(true) {}

RestartInput::RestartInput() : iteration(0), iterations_since_restart(0) {}

RestartCandidate::RestartCandidate()
    : kind(RestartCandidateKind::kNone),
      score(0.0),
      averaged_available(false) {}

RestartResult::RestartResult() : decision(RestartDecision::kNone) {
  next.history_valid = false;
}

bool RestartResult::ok() const {
  return error.ok() && validRestartState(next) &&
         decision != RestartDecision::kNone &&
         candidate.kind != RestartCandidateKind::kNone &&
         validMetrics(candidate.metrics) && std::isfinite(candidate.score) &&
         candidate.score >= 0.0;
}

RestartResult evaluateRestart(const RestartInput& input) {
  if (!validParameters(input.parameters))
    return failure(api::ErrorCode::kInvalidOption,
                   "HiPDLP restart parameters are invalid");
  if (!validRestartState(input.previous))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP restart history is invalid");
  if (!validAveragingState(input.averaging))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP restart averaging state is invalid");
  if (input.iteration < 0 || input.iterations_since_restart < 0 ||
      input.iterations_since_restart > input.iteration ||
      (!input.previous.initialized && input.iterations_since_restart != 0))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP restart iteration counters are invalid");
  if (!validMetrics(input.current))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP current restart metrics are invalid");

  const bool averaged_available = input.averaging.hasAverage();
  if (averaged_available && !validMetrics(input.averaged))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP averaged restart metrics are invalid");

  double current_score = 0.0;
  if (!computeScore(input.current, input.parameters.primal_weight,
                    current_score))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP current restart score is invalid");

  if (!input.previous.initialized || input.iterations_since_restart == 0) {
    RestartResult result;
    result.next = input.previous;
    result.next.restart_metrics = input.current;
    result.next.last_candidate_metrics = input.current;
    result.next.initialized = true;
    result.next.history_valid = true;
    result.decision = RestartDecision::kContinue;
    result.candidate.kind = RestartCandidateKind::kCurrent;
    result.candidate.metrics = input.current;
    result.candidate.score = current_score;
    result.candidate.averaged_available = averaged_available;
    return result;
  }

  RestartCandidateKind candidate_kind = RestartCandidateKind::kCurrent;
  NormalizedConvergenceMetrics candidate_metrics = input.current;
  double candidate_score = current_score;
  if (averaged_available) {
    double averaged_score = 0.0;
    if (!computeScore(input.averaged, input.parameters.primal_weight,
                      averaged_score))
      return failure(api::ErrorCode::kInvalidArgument,
                     "HiPDLP averaged restart score is invalid");
    if (averaged_score <= current_score) {
      candidate_kind = RestartCandidateKind::kAveraged;
      candidate_metrics = input.averaged;
      candidate_score = averaged_score;
    }
  }

  double restart_score = 0.0;
  double last_candidate_score = 0.0;
  if (!computeScore(input.previous.restart_metrics,
                    input.parameters.primal_weight, restart_score) ||
      !computeScore(input.previous.last_candidate_metrics,
                    input.parameters.primal_weight, last_candidate_score))
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP restart history score is invalid");

  const double artificial_boundary =
      input.parameters.artificial_restart_threshold * input.iteration;
  bool should_restart = input.iterations_since_restart >= artificial_boundary;
  if (!should_restart &&
      candidate_score < input.parameters.sufficient_decay * restart_score)
    should_restart = true;
  if (!should_restart &&
      candidate_score < input.parameters.necessary_decay * restart_score &&
      candidate_score > last_candidate_score)
    should_restart = true;

  if (should_restart &&
      input.previous.restart_count == std::numeric_limits<HighsInt>::max())
    return failure(api::ErrorCode::kInvalidArgument,
                   "HiPDLP restart count overflowed");

  RestartResult result;
  result.next = input.previous;
  result.next.last_candidate_metrics = candidate_metrics;
  if (should_restart) {
    result.next.restart_metrics = candidate_metrics;
    ++result.next.restart_count;
  }
  result.decision =
      should_restart ? RestartDecision::kRestart : RestartDecision::kContinue;
  result.candidate.kind = candidate_kind;
  result.candidate.metrics = candidate_metrics;
  result.candidate.score = candidate_score;
  result.candidate.averaged_available = averaged_available;
  return result;
}

RestartState resetRestart() { return RestartState(); }

FixedPointRestartParameters::FixedPointRestartParameters()
    : sufficient_decay(0.2),
      necessary_decay(0.8),
      artificial_threshold(0.36),
      check_interval(40) {}

FixedPointRestartState::FixedPointRestartState()
    : reference_error(0.0),
      previous_trial_error(0.0),
      reference_valid(false),
      previous_trial_valid(false),
      history_valid(true) {}

FixedPointRestartInput::FixedPointRestartInput()
    : current_error(std::numeric_limits<double>::quiet_NaN()),
      halpern_iterations(0),
      total_iterations(0) {}

FixedPointRestartResult::FixedPointRestartResult()
    : decision(RestartDecision::kNone),
      condition(FixedPointRestartCondition::kNone) {
  next.history_valid = false;
}

bool FixedPointRestartResult::ok() const {
  return error.ok() && decision != RestartDecision::kNone &&
         condition != FixedPointRestartCondition::kNone &&
         validFixedPointRestartState(next);
}

FixedPointRestartResult evaluateFixedPointRestart(
    const FixedPointRestartInput& input) {
  if (!validFixedPointRestartParameters(input.parameters))
    return fixedPointFailure(
        api::ErrorCode::kInvalidOption,
        "HiPDLP fixed-point restart parameters are invalid");
  if (!validFixedPointRestartState(input.previous) ||
      !input.previous.reference_valid)
    return fixedPointFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP fixed-point restart history is invalid or has no reference "
        "error");
  if (!std::isfinite(input.current_error) || input.current_error < 0.0)
    return fixedPointFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP current fixed-point error must be finite and nonnegative");
  if (input.total_iterations < 0 || input.halpern_iterations < 0 ||
      input.halpern_iterations > input.total_iterations)
    return fixedPointFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP fixed-point restart iteration counters are invalid");

  const double sufficient_boundary =
      input.parameters.sufficient_decay * input.previous.reference_error;
  const double necessary_boundary =
      input.parameters.necessary_decay * input.previous.reference_error;
  const double artificial_boundary =
      input.parameters.artificial_threshold *
      static_cast<double>(input.total_iterations);
  if (!std::isfinite(sufficient_boundary) ||
      !std::isfinite(necessary_boundary) || !std::isfinite(artificial_boundary))
    return fixedPointFailure(
        api::ErrorCode::kInvalidArgument,
        "HiPDLP fixed-point restart thresholds are outside FP64 range");

  FixedPointRestartCondition condition =
      FixedPointRestartCondition::kNotTriggered;
  if (input.total_iterations == input.parameters.check_interval) {
    condition = FixedPointRestartCondition::kFirstCheckInterval;
  } else if (input.total_iterations > input.parameters.check_interval) {
    if (input.current_error <= sufficient_boundary) {
      condition = FixedPointRestartCondition::kSufficientDecay;
    } else if (input.current_error <= necessary_boundary &&
               input.previous.previous_trial_valid &&
               input.current_error > input.previous.previous_trial_error) {
      condition = FixedPointRestartCondition::kNecessaryDecay;
    } else if (static_cast<double>(input.halpern_iterations) >=
               artificial_boundary) {
      condition = FixedPointRestartCondition::kArtificialThreshold;
    }
  }

  const bool restart = condition != FixedPointRestartCondition::kNotTriggered;

  FixedPointRestartResult result;
  result.decision =
      restart ? RestartDecision::kRestart : RestartDecision::kContinue;
  result.condition = condition;
  result.next.history_valid = true;
  if (!restart) {
    result.next.reference_error = input.previous.reference_error;
    result.next.reference_valid = true;
    result.next.previous_trial_error = input.current_error;
    result.next.previous_trial_valid = true;
  }
  return result;
}

FixedPointRestartState resetFixedPointRestart(const double reference_error) {
  FixedPointRestartState state;
  state.reference_error = reference_error;
  state.reference_valid = true;
  return state;
}

}  // namespace algorithm
}  // namespace hipdlp