#include <cmath>
#include <cstring>
#include <limits>

#include "catch.hpp"
#include "pdlp/hipdlp/algorithm/restart_policy.hpp"

namespace {

hipdlp::algorithm::NormalizedConvergenceMetrics makeMetrics(const double primal,
                                                            const double dual,
                                                            const double gap) {
  hipdlp::algorithm::NormalizedConvergenceMetrics metrics;
  metrics.relative_primal = primal;
  metrics.relative_dual = dual;
  metrics.relative_gap = gap;
  return metrics;
}

hipdlp::algorithm::AveragingState makeAvailableAverage() {
  hipdlp::algorithm::AveragingState averaging;
  averaging.sample_count = 1;
  averaging.total_weight = 1.0;
  return averaging;
}

bool bitwiseEqual(const double left, const double right) {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

void checkZeroMetrics(
    const hipdlp::algorithm::NormalizedConvergenceMetrics& metrics) {
  CHECK(metrics.relative_primal == 0.0);
  CHECK(metrics.relative_dual == 0.0);
  CHECK(metrics.relative_gap == 0.0);
}

void checkInvalidResult(const hipdlp::algorithm::RestartResult& result,
                        const hipdlp::api::ErrorCode expected_code) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.error.code == expected_code);
  CHECK(result.decision == hipdlp::algorithm::RestartDecision::kNone);
  CHECK(result.candidate.kind ==
        hipdlp::algorithm::RestartCandidateKind::kNone);
  CHECK_FALSE(result.candidate.averaged_available);
  CHECK(result.candidate.score == 0.0);
  checkZeroMetrics(result.candidate.metrics);
  CHECK_FALSE(result.next.history_valid);
  CHECK_FALSE(result.next.initialized);
  CHECK(result.next.restart_count == 0);
  checkZeroMetrics(result.next.restart_metrics);
  checkZeroMetrics(result.next.last_candidate_metrics);
}

hipdlp::algorithm::FixedPointRestartInput makeFixedPointInput(
    const double reference_error, const double current_error,
    const HighsInt halpern_iterations, const HighsInt total_iterations) {
  hipdlp::algorithm::FixedPointRestartInput input;
  input.previous = hipdlp::algorithm::resetFixedPointRestart(reference_error);
  input.current_error = current_error;
  input.halpern_iterations = halpern_iterations;
  input.total_iterations = total_iterations;
  return input;
}

void checkInvalidFixedPointResult(
    const hipdlp::algorithm::FixedPointRestartResult& result,
    const hipdlp::api::ErrorCode expected_code) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.error.code == expected_code);
  CHECK(result.decision == hipdlp::algorithm::RestartDecision::kNone);
  CHECK(result.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNone);
  CHECK_FALSE(result.next.history_valid);
  CHECK_FALSE(result.next.reference_valid);
  CHECK_FALSE(result.next.previous_trial_valid);
  CHECK(result.next.reference_error == 0.0);
  CHECK(result.next.previous_trial_error == 0.0);
}

}  // namespace

TEST_CASE("hi-pdlp-restart-policy-triggers-on-an-exact-scalar-sequence",
          "[pdlp][hipdlp][phase7c][restart][policy][determinism]") {
  hipdlp::algorithm::RestartInput input;
  input.current = makeMetrics(10.0, 0.0, 0.0);
  input.iteration = 10;
  input.parameters.artificial_restart_threshold = 0.9;

  const hipdlp::algorithm::RestartResult initial =
      hipdlp::algorithm::evaluateRestart(input);
  INFO(initial.error.message);
  REQUIRE(initial.ok());
  CHECK(initial.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(initial.candidate.kind ==
        hipdlp::algorithm::RestartCandidateKind::kCurrent);
  CHECK(initial.candidate.score == 10.0);
  CHECK(initial.next.restart_count == 0);

  input.previous = initial.next;
  input.current = makeMetrics(5.0, 0.0, 0.0);
  input.averaged = makeMetrics(4.0, 0.0, 0.0);
  input.averaging = makeAvailableAverage();
  input.iteration = 11;
  input.iterations_since_restart = 1;
  const hipdlp::algorithm::RestartResult falling =
      hipdlp::algorithm::evaluateRestart(input);
  INFO(falling.error.message);
  REQUIRE(falling.ok());
  CHECK(falling.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(falling.candidate.kind ==
        hipdlp::algorithm::RestartCandidateKind::kAveraged);
  CHECK(falling.candidate.score == 4.0);

  input.previous = falling.next;
  input.current = makeMetrics(6.0, 0.0, 0.0);
  input.averaged = makeMetrics(7.0, 0.0, 0.0);
  input.iteration = 12;
  input.iterations_since_restart = 2;
  const hipdlp::algorithm::RestartResult triggered =
      hipdlp::algorithm::evaluateRestart(input);
  const hipdlp::algorithm::RestartResult repeated =
      hipdlp::algorithm::evaluateRestart(input);
  INFO(triggered.error.message);
  REQUIRE(triggered.ok());
  REQUIRE(repeated.ok());
  CHECK(triggered.decision == hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(triggered.candidate.kind ==
        hipdlp::algorithm::RestartCandidateKind::kCurrent);
  CHECK(triggered.candidate.score == 6.0);
  CHECK(triggered.next.restart_count == 1);
  CHECK(triggered.next.restart_metrics.relative_primal == 6.0);
  CHECK(repeated.decision == triggered.decision);
  CHECK(repeated.candidate.kind == triggered.candidate.kind);
  CHECK(bitwiseEqual(repeated.candidate.score, triggered.candidate.score));
  CHECK(repeated.next.restart_count == triggered.next.restart_count);
}

TEST_CASE("hi-pdlp-restart-policy-reports-available-candidates-explicitly",
          "[pdlp][hipdlp][phase7c][restart][policy][candidate]") {
  hipdlp::algorithm::RestartInput seed_input;
  seed_input.current = makeMetrics(10.0, 0.0, 0.0);
  seed_input.iteration = 20;
  const hipdlp::algorithm::RestartResult seed =
      hipdlp::algorithm::evaluateRestart(seed_input);
  REQUIRE(seed.ok());

  hipdlp::algorithm::RestartInput input;
  input.previous = seed.next;
  input.current = makeMetrics(3.0, 0.0, 0.0);
  input.averaged = makeMetrics(1.0, 0.0, 0.0);
  input.iteration = 21;
  input.iterations_since_restart = 1;
  const hipdlp::algorithm::RestartResult current_only =
      hipdlp::algorithm::evaluateRestart(input);
  INFO(current_only.error.message);
  REQUIRE(current_only.ok());
  CHECK_FALSE(current_only.candidate.averaged_available);
  CHECK(current_only.candidate.kind ==
        hipdlp::algorithm::RestartCandidateKind::kCurrent);
  CHECK(current_only.candidate.metrics.relative_primal == 3.0);
  CHECK(current_only.candidate.score == 3.0);

  input.averaging = makeAvailableAverage();
  input.averaged = makeMetrics(2.0, 0.0, 0.0);
  const hipdlp::algorithm::RestartResult with_average =
      hipdlp::algorithm::evaluateRestart(input);
  INFO(with_average.error.message);
  REQUIRE(with_average.ok());
  CHECK(with_average.candidate.averaged_available);
  CHECK(with_average.candidate.kind ==
        hipdlp::algorithm::RestartCandidateKind::kAveraged);
  CHECK(with_average.candidate.metrics.relative_primal == 2.0);
  CHECK(with_average.candidate.score == 2.0);
}

TEST_CASE("hi-pdlp-restart-policy-reset-and-counter-are-exact",
          "[pdlp][hipdlp][phase7c][restart][policy][reset][counter]") {
  const hipdlp::algorithm::RestartState first_reset =
      hipdlp::algorithm::resetRestart();
  const hipdlp::algorithm::RestartState second_reset =
      hipdlp::algorithm::resetRestart();
  CHECK(first_reset.history_valid);
  CHECK_FALSE(first_reset.initialized);
  CHECK(first_reset.restart_count == 0);
  checkZeroMetrics(first_reset.restart_metrics);
  checkZeroMetrics(first_reset.last_candidate_metrics);
  CHECK(second_reset.history_valid == first_reset.history_valid);
  CHECK(second_reset.initialized == first_reset.initialized);
  CHECK(second_reset.restart_count == first_reset.restart_count);
  CHECK(bitwiseEqual(second_reset.restart_metrics.relative_primal,
                     first_reset.restart_metrics.relative_primal));
  CHECK(bitwiseEqual(second_reset.restart_metrics.relative_dual,
                     first_reset.restart_metrics.relative_dual));
  CHECK(bitwiseEqual(second_reset.restart_metrics.relative_gap,
                     first_reset.restart_metrics.relative_gap));

  hipdlp::algorithm::RestartInput input;
  input.current = makeMetrics(10.0, 0.0, 0.0);
  input.iteration = 5;
  input.parameters.artificial_restart_threshold = 0.5;
  const hipdlp::algorithm::RestartResult seed =
      hipdlp::algorithm::evaluateRestart(input);
  REQUIRE(seed.ok());

  input.previous = seed.next;
  input.iteration = 10;
  input.iterations_since_restart = 5;
  const hipdlp::algorithm::RestartResult triggered =
      hipdlp::algorithm::evaluateRestart(input);
  REQUIRE(triggered.ok());
  CHECK(triggered.decision == hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(triggered.next.restart_count == 1);

  input.previous = triggered.next;
  input.current = makeMetrics(8.0, 0.0, 0.0);
  input.iterations_since_restart = 0;
  const hipdlp::algorithm::RestartResult after_restart =
      hipdlp::algorithm::evaluateRestart(input);
  REQUIRE(after_restart.ok());
  CHECK(after_restart.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(after_restart.next.restart_count == 1);
  CHECK(after_restart.next.restart_metrics.relative_primal == 8.0);
  CHECK(after_restart.next.last_candidate_metrics.relative_primal == 8.0);

  const hipdlp::algorithm::RestartState cleared =
      hipdlp::algorithm::resetRestart();
  CHECK(cleared.history_valid);
  CHECK_FALSE(cleared.initialized);
  CHECK(cleared.restart_count == 0);
  checkZeroMetrics(cleared.restart_metrics);
  checkZeroMetrics(cleared.last_candidate_metrics);
}

TEST_CASE("hi-pdlp-restart-policy-invalid-input-is-transactional",
          "[pdlp][hipdlp][phase7c][restart][policy][errors]"
          "[transactional]") {
  hipdlp::algorithm::RestartInput input;
  input.parameters.primal_weight = 0.0;
  checkInvalidResult(hipdlp::algorithm::evaluateRestart(input),
                     hipdlp::api::ErrorCode::kInvalidOption);

  input = hipdlp::algorithm::RestartInput();
  input.iterations_since_restart = -1;
  checkInvalidResult(hipdlp::algorithm::evaluateRestart(input),
                     hipdlp::api::ErrorCode::kInvalidArgument);

  input = hipdlp::algorithm::RestartInput();
  input.previous.restart_count = -1;
  checkInvalidResult(hipdlp::algorithm::evaluateRestart(input),
                     hipdlp::api::ErrorCode::kInvalidArgument);

  input = hipdlp::algorithm::RestartInput();
  input.averaging.sample_count = 1;
  input.averaging.total_weight = 0.0;
  checkInvalidResult(hipdlp::algorithm::evaluateRestart(input),
                     hipdlp::api::ErrorCode::kInvalidArgument);

  input = hipdlp::algorithm::RestartInput();
  input.current.relative_gap = std::numeric_limits<double>::quiet_NaN();
  checkInvalidResult(hipdlp::algorithm::evaluateRestart(input),
                     hipdlp::api::ErrorCode::kInvalidArgument);

  input = hipdlp::algorithm::RestartInput();
  input.averaging = makeAvailableAverage();
  input.averaged.relative_dual = std::numeric_limits<double>::infinity();
  checkInvalidResult(hipdlp::algorithm::evaluateRestart(input),
                     hipdlp::api::ErrorCode::kInvalidArgument);
}
TEST_CASE("hi-pdlp-fixed-point-restart-first-block-is-unconditional",
          "[pdlp][hipdlp][phase7g][restart][fixed-point][first-block]") {
  hipdlp::algorithm::FixedPointRestartInput input =
      makeFixedPointInput(1.0, 1.0e6, 0, 40);
  REQUIRE(input.parameters.check_interval == 40);
  REQUIRE(input.parameters.sufficient_decay == 0.2);
  REQUIRE(input.parameters.necessary_decay == 0.8);
  REQUIRE(input.parameters.artificial_threshold == 0.36);

  const hipdlp::algorithm::FixedPointRestartResult at_interval =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  INFO(at_interval.error.message);
  REQUIRE(at_interval.ok());
  CHECK(at_interval.decision == hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(at_interval.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kFirstCheckInterval);
  CHECK(at_interval.next.history_valid);
  CHECK_FALSE(at_interval.next.reference_valid);
  CHECK_FALSE(at_interval.next.previous_trial_valid);
  CHECK(at_interval.next.reference_error == 0.0);
  CHECK(at_interval.next.previous_trial_error == 0.0);

  // Below the interval nothing fires, even when every post-interval condition
  // would have been satisfied.
  input = makeFixedPointInput(10.0, 0.0, 39, 39);
  const hipdlp::algorithm::FixedPointRestartResult below_zero_error =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  INFO(below_zero_error.error.message);
  REQUIRE(below_zero_error.ok());
  CHECK(below_zero_error.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(below_zero_error.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
  CHECK(below_zero_error.next.reference_valid);
  CHECK(below_zero_error.next.reference_error == 10.0);
  CHECK(below_zero_error.next.previous_trial_valid);
  CHECK(below_zero_error.next.previous_trial_error == 0.0);

  input = makeFixedPointInput(10.0, 1.0e12, 39, 39);
  input.previous.previous_trial_error = 1.0;
  input.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartResult below_huge_error =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(below_huge_error.ok());
  CHECK(below_huge_error.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(below_huge_error.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);

  input = makeFixedPointInput(10.0, 5.0, 3, 4);
  input.parameters.check_interval = 4;
  const hipdlp::algorithm::FixedPointRestartResult custom_interval =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(custom_interval.ok());
  CHECK(custom_interval.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(custom_interval.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kFirstCheckInterval);

  input.total_iterations = 3;
  input.halpern_iterations = 3;
  const hipdlp::algorithm::FixedPointRestartResult custom_below =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(custom_below.ok());
  CHECK(custom_below.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(custom_below.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
}

TEST_CASE("hi-pdlp-fixed-point-restart-conditions-fire-in-isolation",
          "[pdlp][hipdlp][phase7g][restart][fixed-point][condition]") {
  hipdlp::algorithm::FixedPointRestartInput input =
      makeFixedPointInput(10.0, 1.0, 1, 1000);
  const hipdlp::algorithm::FixedPointRestartResult sufficient =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  INFO(sufficient.error.message);
  REQUIRE(sufficient.ok());
  CHECK(sufficient.decision == hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(sufficient.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kSufficientDecay);

  input = makeFixedPointInput(10.0, 7.5, 1, 1000);
  input.previous.previous_trial_error = 7.0;
  input.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartResult necessary =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  INFO(necessary.error.message);
  REQUIRE(necessary.ok());
  CHECK(necessary.decision == hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(necessary.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNecessaryDecay);

  input = makeFixedPointInput(10.0, 9.5, 400, 1000);
  input.previous.previous_trial_error = 1.0;
  input.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartResult artificial =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  INFO(artificial.error.message);
  REQUIRE(artificial.ok());
  CHECK(artificial.decision == hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(artificial.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kArtificialThreshold);

  // The necessary-decay window is entered but non-improvement fails, so the
  // legacy fall-through still reaches the artificial check.
  input = makeFixedPointInput(10.0, 7.5, 400, 1000);
  input.previous.previous_trial_error = 8.0;
  input.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartResult fall_through =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(fall_through.ok());
  CHECK(fall_through.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kArtificialThreshold);

  input = makeFixedPointInput(10.0, 9.5, 1, 1000);
  input.previous.previous_trial_error = 1.0;
  input.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartResult none =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(none.ok());
  CHECK(none.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(none.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
}

TEST_CASE("hi-pdlp-fixed-point-restart-boundaries-match-the-legacy-comparisons",
          "[pdlp][hipdlp][phase7g][restart][fixed-point][boundary]") {
  REQUIRE(0.2 * 10.0 == 2.0);
  REQUIRE(0.8 * 10.0 == 8.0);
  REQUIRE(0.36 * 1000.0 == 360.0);

  // Sufficient decay is inclusive.
  hipdlp::algorithm::FixedPointRestartInput input =
      makeFixedPointInput(10.0, 2.0, 1, 1000);
  const hipdlp::algorithm::FixedPointRestartResult at_sufficient =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(at_sufficient.ok());
  CHECK(at_sufficient.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kSufficientDecay);

  input.current_error = std::nextafter(2.0, 10.0);
  const hipdlp::algorithm::FixedPointRestartResult above_sufficient =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(above_sufficient.ok());
  CHECK(above_sufficient.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(above_sufficient.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);

  // Necessary decay is inclusive, but non-improvement is strict.
  input = makeFixedPointInput(10.0, 8.0, 1, 1000);
  input.previous.previous_trial_error = std::nextafter(8.0, 0.0);
  input.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartResult at_necessary =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(at_necessary.ok());
  CHECK(at_necessary.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNecessaryDecay);

  input.current_error = std::nextafter(8.0, 10.0);
  const hipdlp::algorithm::FixedPointRestartResult above_necessary =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(above_necessary.ok());
  CHECK(above_necessary.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);

  input = makeFixedPointInput(10.0, 8.0, 1, 1000);
  input.previous.previous_trial_error = 8.0;
  input.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartResult equal_trial =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(equal_trial.ok());
  CHECK(equal_trial.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);

  // An invalid previous trial can never satisfy the necessary-decay branch.
  input = makeFixedPointInput(10.0, 8.0, 1, 1000);
  REQUIRE_FALSE(input.previous.previous_trial_valid);
  const hipdlp::algorithm::FixedPointRestartResult no_trial =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(no_trial.ok());
  CHECK(no_trial.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(no_trial.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);

  // The artificial threshold product is inclusive.
  input = makeFixedPointInput(10.0, 9.5, 360, 1000);
  const hipdlp::algorithm::FixedPointRestartResult at_artificial =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(at_artificial.ok());
  CHECK(at_artificial.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kArtificialThreshold);

  input.halpern_iterations = 359;
  const hipdlp::algorithm::FixedPointRestartResult below_artificial =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(below_artificial.ok());
  CHECK(below_artificial.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(below_artificial.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
}

TEST_CASE("hi-pdlp-fixed-point-restart-chains-reference-and-trial-scalars",
          "[pdlp][hipdlp][phase7g][restart][fixed-point][chain]") {
  hipdlp::algorithm::FixedPointRestartState state =
      hipdlp::algorithm::resetFixedPointRestart(10.0);
  CHECK(state.history_valid);
  CHECK(state.reference_valid);
  CHECK(state.reference_error == 10.0);
  CHECK_FALSE(state.previous_trial_valid);
  CHECK(state.previous_trial_error == 0.0);

  hipdlp::algorithm::FixedPointRestartInput input;
  input.previous = state;
  input.current_error = 9.0;
  input.halpern_iterations = 10;
  input.total_iterations = 1000;
  const hipdlp::algorithm::FixedPointRestartResult first =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  INFO(first.error.message);
  REQUIRE(first.ok());
  CHECK(first.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(first.next.reference_valid);
  CHECK(first.next.reference_error == 10.0);
  CHECK(first.next.previous_trial_valid);
  CHECK(first.next.previous_trial_error == 9.0);

  input.previous = first.next;
  input.current_error = 7.0;
  input.halpern_iterations = 50;
  input.total_iterations = 1040;
  const hipdlp::algorithm::FixedPointRestartResult second =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(second.ok());
  CHECK(second.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(second.next.reference_error == 10.0);
  CHECK(second.next.previous_trial_error == 7.0);

  input.previous = second.next;
  input.current_error = 7.5;
  input.halpern_iterations = 90;
  input.total_iterations = 1080;
  const hipdlp::algorithm::FixedPointRestartResult third =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(third.ok());
  CHECK(third.decision == hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(third.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNecessaryDecay);
  CHECK(third.next.history_valid);
  CHECK_FALSE(third.next.reference_valid);
  CHECK_FALSE(third.next.previous_trial_valid);
  CHECK(third.next.reference_error == 0.0);
  CHECK(third.next.previous_trial_error == 0.0);

  // A restarted cycle cannot be chained until the caller supplies a new
  // reference fixed-point error.
  input.previous = third.next;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input.previous = hipdlp::algorithm::resetFixedPointRestart(7.5);
  input.current_error = 5.0;
  input.halpern_iterations = 40;
  input.total_iterations = 1120;
  const hipdlp::algorithm::FixedPointRestartResult fourth =
      hipdlp::algorithm::evaluateFixedPointRestart(input);
  REQUIRE(fourth.ok());
  CHECK(fourth.decision == hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(fourth.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
  CHECK(fourth.next.reference_error == 7.5);
  CHECK(fourth.next.previous_trial_valid);
  CHECK(fourth.next.previous_trial_error == 5.0);
}

TEST_CASE("hi-pdlp-fixed-point-restart-invalid-input-is-transactional",
          "[pdlp][hipdlp][phase7g][restart][fixed-point][errors]"
          "[transactional][determinism]") {
  hipdlp::algorithm::FixedPointRestartInput input =
      makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.parameters.sufficient_decay = 0.0;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidOption);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.parameters.necessary_decay = input.parameters.sufficient_decay;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidOption);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.parameters.necessary_decay = 1.0;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidOption);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.parameters.artificial_threshold = 0.0;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidOption);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.parameters.artificial_threshold =
      std::numeric_limits<double>::quiet_NaN();
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidOption);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.parameters.check_interval = 0;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidOption);

  input = hipdlp::algorithm::FixedPointRestartInput();
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.previous.history_valid = false;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(-1.0, 5.0, 1, 1000);
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(std::numeric_limits<double>::infinity(), 5.0, 1,
                              1000);
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.previous.previous_trial_error = 5.0;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, 5.0, 1, 1000);
  input.previous.previous_trial_valid = true;
  input.previous.previous_trial_error =
      -std::numeric_limits<double>::quiet_NaN();
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, -1.0, 1, 1000);
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, std::numeric_limits<double>::infinity(), 1,
                              1000);
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, 5.0, 1, -1);
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, 5.0, -1, 1000);
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  input = makeFixedPointInput(10.0, 5.0, 1001, 1000);
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  // The artificial threshold product overflows FP64.
  input = makeFixedPointInput(10.0, 9.5, 1, 100);
  input.parameters.artificial_threshold = 1.0e308;
  checkInvalidFixedPointResult(
      hipdlp::algorithm::evaluateFixedPointRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument);

  const hipdlp::algorithm::FixedPointRestartInput stable =
      makeFixedPointInput(10.0, 7.5, 90, 1080);
  hipdlp::algorithm::FixedPointRestartInput mutated = stable;
  mutated.previous.previous_trial_error = 7.0;
  mutated.previous.previous_trial_valid = true;
  const hipdlp::algorithm::FixedPointRestartInput expected = mutated;
  const hipdlp::algorithm::FixedPointRestartResult once =
      hipdlp::algorithm::evaluateFixedPointRestart(mutated);
  const hipdlp::algorithm::FixedPointRestartResult twice =
      hipdlp::algorithm::evaluateFixedPointRestart(mutated);
  REQUIRE(once.ok());
  REQUIRE(twice.ok());
  CHECK(once.decision == twice.decision);
  CHECK(once.condition == twice.condition);
  CHECK(once.next.reference_valid == twice.next.reference_valid);
  CHECK(once.next.previous_trial_valid == twice.next.previous_trial_valid);
  CHECK(bitwiseEqual(once.next.reference_error, twice.next.reference_error));
  CHECK(bitwiseEqual(once.next.previous_trial_error,
                     twice.next.previous_trial_error));

  CHECK(bitwiseEqual(mutated.current_error, expected.current_error));
  CHECK(mutated.halpern_iterations == expected.halpern_iterations);
  CHECK(mutated.total_iterations == expected.total_iterations);
  CHECK(mutated.previous.history_valid == expected.previous.history_valid);
  CHECK(mutated.previous.reference_valid == expected.previous.reference_valid);
  CHECK(mutated.previous.previous_trial_valid ==
        expected.previous.previous_trial_valid);
  CHECK(bitwiseEqual(mutated.previous.reference_error,
                     expected.previous.reference_error));
  CHECK(bitwiseEqual(mutated.previous.previous_trial_error,
                     expected.previous.previous_trial_error));
  CHECK(bitwiseEqual(mutated.parameters.sufficient_decay,
                     expected.parameters.sufficient_decay));
  CHECK(bitwiseEqual(mutated.parameters.necessary_decay,
                     expected.parameters.necessary_decay));
  CHECK(bitwiseEqual(mutated.parameters.artificial_threshold,
                     expected.parameters.artificial_threshold));
  CHECK(mutated.parameters.check_interval ==
        expected.parameters.check_interval);
}
