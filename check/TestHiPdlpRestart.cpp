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