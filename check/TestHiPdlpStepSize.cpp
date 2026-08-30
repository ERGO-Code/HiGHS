#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

#include "catch.hpp"
#include "pdlp/hipdlp/algorithm/step_size_policy.hpp"

namespace {

hipdlp::algorithm::StepSizeInitializationInput makeInput(
    const double operator_norm_squared, const double unscaled_rhs_norm,
    const double unscaled_cost_norm) {
  hipdlp::algorithm::StepSizeInitializationInput input;
  input.operator_norm_squared = operator_norm_squared;
  input.unscaled_rhs_norm = unscaled_rhs_norm;
  input.unscaled_cost_norm = unscaled_cost_norm;
  return input;
}

bool bitwiseEqual(const double left, const double right) {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

void checkInvalidResult(
    const hipdlp::algorithm::StepSizeInitializationResult& result,
    const hipdlp::api::ErrorCode expected_code) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.error.code == expected_code);
  CHECK(result.step_sizes.primal_weight == 0.0);
  CHECK(result.step_sizes.base_step == 0.0);
  CHECK(result.step_sizes.primal_step == 0.0);
  CHECK(result.step_sizes.dual_step == 0.0);
}

hipdlp::algorithm::PrimalWeightInput makeUpdateInput(
    const hipdlp::algorithm::PrimalWeightState& previous,
    const double primal_step, const double dual_step,
    const double primal_distance, const double dual_distance,
    const double relative_primal, const double relative_dual) {
  hipdlp::algorithm::PrimalWeightInput input;
  input.previous = previous;
  input.primal_step = primal_step;
  input.dual_step = dual_step;
  input.primal_distance = primal_distance;
  input.dual_distance = dual_distance;
  input.relative_primal_residual = relative_primal;
  input.relative_dual_residual = relative_dual;
  return input;
}

double pidError(const double dual_distance, const double primal_distance,
                const double primal_weight) {
  return std::log(dual_distance) - std::log(primal_distance) -
         std::log(primal_weight);
}

// Mirrors the legacy PID recurrence so the tests never reuse solver code.
hipdlp::algorithm::PrimalWeightState referenceUpdate(
    const hipdlp::algorithm::PrimalWeightState& previous,
    const double primal_distance, const double dual_distance,
    const hipdlp::algorithm::PrimalWeightParameters& parameters) {
  const double error =
      pidError(dual_distance, primal_distance, previous.primal_weight);
  const double error_sum =
      parameters.integral_smoothing * previous.error_sum + error;
  const double delta_error = error - previous.previous_error;
  hipdlp::algorithm::PrimalWeightState next = previous;
  next.primal_weight =
      previous.primal_weight *
      std::exp(parameters.k_p * error + parameters.k_i * error_sum +
               parameters.k_d * delta_error);
  next.error_sum = error_sum;
  next.previous_error = error;
  return next;
}

void checkInvalidUpdate(
    const hipdlp::algorithm::PrimalWeightUpdateResult& result,
    const hipdlp::api::ErrorCode expected_code) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.error.code == expected_code);
  CHECK(result.outcome == hipdlp::algorithm::PrimalWeightOutcome::kNone);
  CHECK(result.base_step == 0.0);
  CHECK(result.primal_step == 0.0);
  CHECK(result.dual_step == 0.0);
  CHECK(result.next.primal_weight == 0.0);
  CHECK(result.next.error_sum == 0.0);
  CHECK(result.next.previous_error == 0.0);
  CHECK(result.next.best_primal_weight == 0.0);
  CHECK(result.next.best_residual_gap == 0.0);
  CHECK_FALSE(result.next.best_gap_valid);
  CHECK_FALSE(result.next.history_valid);
}

void checkRejected(const hipdlp::algorithm::PrimalWeightInput& input,
                   const hipdlp::api::ErrorCode expected_code) {
  checkInvalidUpdate(hipdlp::algorithm::updatePrimalWeight(input),
                     expected_code);
}

void checkSameState(const hipdlp::algorithm::PrimalWeightState& actual,
                    const hipdlp::algorithm::PrimalWeightState& expected) {
  CHECK(bitwiseEqual(actual.primal_weight, expected.primal_weight));
  CHECK(bitwiseEqual(actual.error_sum, expected.error_sum));
  CHECK(bitwiseEqual(actual.previous_error, expected.previous_error));
  CHECK(bitwiseEqual(actual.best_primal_weight, expected.best_primal_weight));
  CHECK(bitwiseEqual(actual.best_residual_gap, expected.best_residual_gap));
  CHECK(actual.best_gap_valid == expected.best_gap_valid);
  CHECK(actual.history_valid == expected.history_valid);
}

void checkSameInput(const hipdlp::algorithm::PrimalWeightInput& actual,
                    const hipdlp::algorithm::PrimalWeightInput& expected) {
  checkSameState(actual.previous, expected.previous);
  CHECK(bitwiseEqual(actual.primal_step, expected.primal_step));
  CHECK(bitwiseEqual(actual.dual_step, expected.dual_step));
  CHECK(bitwiseEqual(actual.primal_distance, expected.primal_distance));
  CHECK(bitwiseEqual(actual.dual_distance, expected.dual_distance));
  CHECK(bitwiseEqual(actual.relative_primal_residual,
                     expected.relative_primal_residual));
  CHECK(bitwiseEqual(actual.relative_dual_residual,
                     expected.relative_dual_residual));
  CHECK(bitwiseEqual(actual.parameters.k_p, expected.parameters.k_p));
  CHECK(bitwiseEqual(actual.parameters.k_i, expected.parameters.k_i));
  CHECK(bitwiseEqual(actual.parameters.k_d, expected.parameters.k_d));
  CHECK(bitwiseEqual(actual.parameters.integral_smoothing,
                     expected.parameters.integral_smoothing));
}

void checkSameResult(
    const hipdlp::algorithm::PrimalWeightUpdateResult& actual,
    const hipdlp::algorithm::PrimalWeightUpdateResult& expected) {
  CHECK(actual.error.code == expected.error.code);
  CHECK(actual.outcome == expected.outcome);
  CHECK(bitwiseEqual(actual.base_step, expected.base_step));
  CHECK(bitwiseEqual(actual.primal_step, expected.primal_step));
  CHECK(bitwiseEqual(actual.dual_step, expected.dual_step));
  checkSameState(actual.next, expected.next);
}

}  // namespace

TEST_CASE("hi-pdlp-step-size-initialization-is-exact-on-a-scalar-fixture",
          "[pdlp][hipdlp][phase7d][stepsize][initialization]") {
  const hipdlp::algorithm::StepSizeInitializationParameters defaults;
  CHECK(defaults.safety_factor == 0.998);

  const hipdlp::algorithm::StepSizeInitializationResult result =
      hipdlp::algorithm::initializeStepSizes(makeInput(4.0, 1.0, 3.0));
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.error.code == hipdlp::api::ErrorCode::kNone);
  CHECK(result.step_sizes.primal_weight == 2.0);
  CHECK(result.step_sizes.base_step == 0.499);
  CHECK(result.step_sizes.primal_step == 0.2495);
  CHECK(result.step_sizes.dual_step == 0.998);

  hipdlp::algorithm::StepSizeInitializationInput unit_factor =
      makeInput(4.0, 1.0, 3.0);
  unit_factor.parameters.safety_factor = 1.0;
  const hipdlp::algorithm::StepSizeInitializationResult unit_result =
      hipdlp::algorithm::initializeStepSizes(unit_factor);
  INFO(unit_result.error.message);
  REQUIRE(unit_result.ok());
  CHECK(unit_result.step_sizes.base_step == 0.5);
  CHECK(unit_result.step_sizes.primal_step == 0.25);
  CHECK(unit_result.step_sizes.dual_step == 1.0);

  hipdlp::algorithm::StepSizeInitializationInput zero_norms =
      makeInput(0.25, 0.0, 0.0);
  zero_norms.parameters.safety_factor = 0.5;
  const hipdlp::algorithm::StepSizeInitializationResult zero_result =
      hipdlp::algorithm::initializeStepSizes(zero_norms);
  INFO(zero_result.error.message);
  REQUIRE(zero_result.ok());
  CHECK(zero_result.step_sizes.primal_weight == 1.0);
  CHECK(zero_result.step_sizes.base_step == 1.0);
  CHECK(zero_result.step_sizes.primal_step == 1.0);
  CHECK(zero_result.step_sizes.dual_step == 1.0);
}

TEST_CASE("hi-pdlp-step-size-primal-dual-decomposition-is-exact",
          "[pdlp][hipdlp][phase7d][stepsize][decomposition]") {
  hipdlp::algorithm::StepSizeInitializationInput input =
      makeInput(16.0, 3.0, 7.0);
  input.parameters.safety_factor = 0.5;
  const hipdlp::algorithm::StepSizeInitializationResult result =
      hipdlp::algorithm::initializeStepSizes(input);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.step_sizes.primal_weight == 2.0);
  CHECK(result.step_sizes.base_step == 0.125);
  CHECK(result.step_sizes.primal_step == 0.0625);
  CHECK(result.step_sizes.dual_step == 0.25);
  CHECK(bitwiseEqual(
      result.step_sizes.primal_step * result.step_sizes.primal_weight,
      result.step_sizes.base_step));
  CHECK(bitwiseEqual(
      result.step_sizes.dual_step / result.step_sizes.primal_weight,
      result.step_sizes.base_step));
  CHECK(bitwiseEqual(
      std::sqrt(result.step_sizes.primal_step * result.step_sizes.dual_step),
      result.step_sizes.base_step));

  hipdlp::algorithm::StepSizeInitializationInput dual_heavy =
      makeInput(4.0, 3.0, 0.0);
  dual_heavy.parameters.safety_factor = 0.5;
  const hipdlp::algorithm::StepSizeInitializationResult dual_heavy_result =
      hipdlp::algorithm::initializeStepSizes(dual_heavy);
  INFO(dual_heavy_result.error.message);
  REQUIRE(dual_heavy_result.ok());
  CHECK(dual_heavy_result.step_sizes.primal_weight == 0.25);
  CHECK(dual_heavy_result.step_sizes.base_step == 0.25);
  CHECK(dual_heavy_result.step_sizes.primal_step == 1.0);
  CHECK(dual_heavy_result.step_sizes.dual_step == 0.0625);
  CHECK(dual_heavy_result.step_sizes.primal_step >
        dual_heavy_result.step_sizes.dual_step);
  CHECK(bitwiseEqual(std::sqrt(dual_heavy_result.step_sizes.primal_step *
                               dual_heavy_result.step_sizes.dual_step),
                     dual_heavy_result.step_sizes.base_step));

  const hipdlp::algorithm::StepSizeInitializationResult default_result =
      hipdlp::algorithm::initializeStepSizes(makeInput(4.0, 1.0, 3.0));
  REQUIRE(default_result.ok());
  CHECK(bitwiseEqual(std::sqrt(default_result.step_sizes.primal_step *
                               default_result.step_sizes.dual_step),
                     default_result.step_sizes.base_step));
}

TEST_CASE("hi-pdlp-step-size-initialization-is-bitwise-deterministic",
          "[pdlp][hipdlp][phase7d][stepsize][determinism]") {
  const hipdlp::algorithm::StepSizeInitializationInput input =
      makeInput(3.7, 2.3, 5.9);
  const hipdlp::algorithm::StepSizeInitializationResult first =
      hipdlp::algorithm::initializeStepSizes(input);
  INFO(first.error.message);
  REQUIRE(first.ok());

  for (int repeat = 0; repeat < 4; ++repeat) {
    const hipdlp::algorithm::StepSizeInitializationResult repeated =
        hipdlp::algorithm::initializeStepSizes(makeInput(3.7, 2.3, 5.9));
    REQUIRE(repeated.ok());
    CHECK(bitwiseEqual(repeated.step_sizes.primal_weight,
                       first.step_sizes.primal_weight));
    CHECK(bitwiseEqual(repeated.step_sizes.base_step,
                       first.step_sizes.base_step));
    CHECK(bitwiseEqual(repeated.step_sizes.primal_step,
                       first.step_sizes.primal_step));
    CHECK(bitwiseEqual(repeated.step_sizes.dual_step,
                       first.step_sizes.dual_step));
  }

  const hipdlp::algorithm::StepSizeInitializationResult interleaved =
      hipdlp::algorithm::initializeStepSizes(makeInput(1.0, 0.0, 0.0));
  REQUIRE(interleaved.ok());
  const hipdlp::algorithm::StepSizeInitializationResult after =
      hipdlp::algorithm::initializeStepSizes(input);
  REQUIRE(after.ok());
  CHECK(bitwiseEqual(after.step_sizes.primal_weight,
                     first.step_sizes.primal_weight));
  CHECK(bitwiseEqual(after.step_sizes.base_step, first.step_sizes.base_step));
  CHECK(
      bitwiseEqual(after.step_sizes.primal_step, first.step_sizes.primal_step));
  CHECK(bitwiseEqual(after.step_sizes.dual_step, first.step_sizes.dual_step));
}

TEST_CASE("hi-pdlp-step-size-initialization-invalid-input-is-transactional",
          "[pdlp][hipdlp][phase7d][stepsize][errors][transactional]") {
  const double nan_value = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  const double huge = std::numeric_limits<double>::max();

  checkInvalidResult(hipdlp::algorithm::initializeStepSizes(
                         hipdlp::algorithm::StepSizeInitializationInput()),
                     hipdlp::api::ErrorCode::kInvalidArgument);

  hipdlp::algorithm::StepSizeInitializationInput input =
      makeInput(4.0, 1.0, 3.0);
  input.parameters.safety_factor = 0.0;
  checkInvalidResult(hipdlp::algorithm::initializeStepSizes(input),
                     hipdlp::api::ErrorCode::kInvalidOption);

  input = makeInput(4.0, 1.0, 3.0);
  input.parameters.safety_factor = -0.5;
  checkInvalidResult(hipdlp::algorithm::initializeStepSizes(input),
                     hipdlp::api::ErrorCode::kInvalidOption);

  input = makeInput(4.0, 1.0, 3.0);
  input.parameters.safety_factor = std::nextafter(1.0, 2.0);
  checkInvalidResult(hipdlp::algorithm::initializeStepSizes(input),
                     hipdlp::api::ErrorCode::kInvalidOption);

  input = makeInput(4.0, 1.0, 3.0);
  input.parameters.safety_factor = nan_value;
  checkInvalidResult(hipdlp::algorithm::initializeStepSizes(input),
                     hipdlp::api::ErrorCode::kInvalidOption);

  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(0.0, 1.0, 3.0)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(-4.0, 1.0, 3.0)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(infinity, 1.0, 3.0)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(nan_value, 1.0, 3.0)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(4.0, -1.0, 3.0)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(4.0, infinity, 3.0)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(4.0, 1.0, -3.0)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(4.0, 1.0, nan_value)),
      hipdlp::api::ErrorCode::kInvalidArgument);

  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(1e300, 0.0, huge)),
      hipdlp::api::ErrorCode::kInvalidArgument);
  checkInvalidResult(
      hipdlp::algorithm::initializeStepSizes(makeInput(1e-300, 0.0, huge)),
      hipdlp::api::ErrorCode::kInvalidArgument);
}

TEST_CASE("hi-pdlp-primal-weight-first-update-is-exact-on-a-scalar-fixture",
          "[pdlp][hipdlp][phase7f][stepsize][primalweight]") {
  const hipdlp::algorithm::PrimalWeightParameters defaults;
  CHECK(defaults.k_p == 0.99);
  CHECK(defaults.k_i == 0.01);
  CHECK(defaults.k_d == 0.0);
  CHECK(defaults.integral_smoothing == 0.3);

  const hipdlp::algorithm::PrimalWeightState initial =
      hipdlp::algorithm::resetPrimalWeight(2.0);
  CHECK(initial.primal_weight == 2.0);
  CHECK(initial.best_primal_weight == 2.0);
  CHECK(initial.error_sum == 0.0);
  CHECK(initial.previous_error == 0.0);
  CHECK(initial.best_residual_gap == 0.0);
  CHECK_FALSE(initial.best_gap_valid);
  CHECK(initial.history_valid);

  // eta = sqrt(0.25 * 1.0) = 0.5 for a primal weight of 2.
  const hipdlp::algorithm::PrimalWeightUpdateResult result =
      hipdlp::algorithm::updatePrimalWeight(
          makeUpdateInput(initial, 0.25, 1.0, 1.0, 8.0, 1e-3, 1e-3));
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.error.code == hipdlp::api::ErrorCode::kNone);
  CHECK(result.outcome == hipdlp::algorithm::PrimalWeightOutcome::kPidUpdate);
  CHECK(result.base_step == 0.5);
  CHECK(result.next.error_sum == Approx(std::log(4.0)).epsilon(1e-15));
  CHECK(bitwiseEqual(result.next.previous_error, result.next.error_sum));
  CHECK(result.next.primal_weight == Approx(8.0).epsilon(1e-14));
  CHECK(result.primal_step == Approx(0.0625).epsilon(1e-14));
  CHECK(result.dual_step == Approx(4.0).epsilon(1e-14));
  CHECK(result.next.history_valid);
  CHECK(result.next.best_gap_valid);
  CHECK(result.next.best_residual_gap == 0.0);
  CHECK(
      bitwiseEqual(result.next.best_primal_weight, result.next.primal_weight));

  const hipdlp::algorithm::PrimalWeightState reference =
      referenceUpdate(initial, 1.0, 8.0, defaults);
  CHECK(result.next.primal_weight ==
        Approx(reference.primal_weight).epsilon(1e-15));
  CHECK(result.next.error_sum == Approx(reference.error_sum).epsilon(1e-15));
  CHECK(result.primal_step ==
        Approx(result.base_step / reference.primal_weight).epsilon(1e-14));
  CHECK(result.dual_step ==
        Approx(result.base_step * reference.primal_weight).epsilon(1e-14));
}

TEST_CASE("hi-pdlp-primal-weight-updates-evolve-pid-state-and-track-best",
          "[pdlp][hipdlp][phase7f][stepsize][primalweight][evolution]") {
  const hipdlp::algorithm::PrimalWeightParameters defaults;
  hipdlp::algorithm::PrimalWeightState state =
      hipdlp::algorithm::resetPrimalWeight(2.0);
  hipdlp::algorithm::PrimalWeightState reference = state;
  double primal_step = 0.25;
  double dual_step = 1.0;

  const double primal_distances[3] = {1.0, 2.0, 0.5};
  const double dual_distances[3] = {8.0, 1.0, 3.0};
  const double relative_primals[3] = {1e-3, 1e-3, 1e-3};
  const double relative_duals[3] = {1e-2, 1.5e-3, 1e-1};

  double best_gap = 0.0;
  bool best_gap_valid = false;
  double best_weight = state.best_primal_weight;

  for (int step = 0; step < 3; ++step) {
    INFO("step " << step);
    const hipdlp::algorithm::PrimalWeightUpdateResult result =
        hipdlp::algorithm::updatePrimalWeight(
            makeUpdateInput(state, primal_step, dual_step,
                            primal_distances[step], dual_distances[step],
                            relative_primals[step], relative_duals[step]));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.outcome == hipdlp::algorithm::PrimalWeightOutcome::kPidUpdate);
    CHECK(result.base_step == Approx(0.5).epsilon(1e-12));

    reference = referenceUpdate(reference, primal_distances[step],
                                dual_distances[step], defaults);
    CHECK(result.next.primal_weight ==
          Approx(reference.primal_weight).epsilon(1e-13));
    CHECK(result.next.error_sum == Approx(reference.error_sum).epsilon(1e-13));
    CHECK(result.next.previous_error ==
          Approx(reference.previous_error).epsilon(1e-13));

    const double gap =
        std::fabs(std::log10(relative_duals[step] / relative_primals[step]));
    if (!best_gap_valid || gap < best_gap) {
      best_gap = gap;
      best_gap_valid = true;
      best_weight = result.next.primal_weight;
    }
    CHECK(result.next.best_gap_valid == best_gap_valid);
    CHECK(result.next.best_residual_gap == Approx(best_gap).epsilon(1e-13));
    CHECK(bitwiseEqual(result.next.best_primal_weight, best_weight));

    CHECK(result.primal_step ==
          Approx(result.base_step / result.next.primal_weight).epsilon(1e-14));
    CHECK(result.dual_step ==
          Approx(result.base_step * result.next.primal_weight).epsilon(1e-14));

    state = result.next;
    primal_step = result.primal_step;
    dual_step = result.dual_step;
  }
  CHECK(state.best_residual_gap ==
        Approx(std::fabs(std::log10(1.5))).epsilon(1e-13));
  CHECK_FALSE(bitwiseEqual(state.best_primal_weight, state.primal_weight));

  // The derivative term only contributes once the PID history is non-zero.
  const hipdlp::algorithm::PrimalWeightUpdateResult first =
      hipdlp::algorithm::updatePrimalWeight(
          makeUpdateInput(hipdlp::algorithm::resetPrimalWeight(2.0), 0.25, 1.0,
                          1.0, 8.0, 1e-3, 1e-3));
  REQUIRE(first.ok());
  hipdlp::algorithm::PrimalWeightInput derivative_input = makeUpdateInput(
      first.next, first.primal_step, first.dual_step, 2.0, 1.0, 1e-3, 1e-3);
  const hipdlp::algorithm::PrimalWeightUpdateResult without_derivative =
      hipdlp::algorithm::updatePrimalWeight(derivative_input);
  REQUIRE(without_derivative.ok());
  derivative_input.parameters.k_d = 0.5;
  const hipdlp::algorithm::PrimalWeightUpdateResult with_derivative =
      hipdlp::algorithm::updatePrimalWeight(derivative_input);
  REQUIRE(with_derivative.ok());
  CHECK(with_derivative.next.primal_weight <
        without_derivative.next.primal_weight);
  CHECK(
      with_derivative.next.primal_weight ==
      Approx(referenceUpdate(first.next, 2.0, 1.0, derivative_input.parameters)
                 .primal_weight)
          .epsilon(1e-13));
}

TEST_CASE("hi-pdlp-primal-weight-guard-windows-are-strict-and-fall-back",
          "[pdlp][hipdlp][phase7f][stepsize][primalweight][fallback]") {
  hipdlp::algorithm::PrimalWeightState history =
      hipdlp::algorithm::resetPrimalWeight(5.0);
  history.error_sum = 1.25;
  history.previous_error = -0.75;
  history.best_primal_weight = 3.0;
  history.best_residual_gap = 0.5;
  history.best_gap_valid = true;

  struct FallbackCase {
    double primal_distance;
    double dual_distance;
    double relative_primal;
    double relative_dual;
    double expected_best_gap;
  };
  const FallbackCase cases[] = {
      {1e-16, 8.0, 1.0, 1.0, 0.0}, {1e12, 8.0, 1.0, 1.0, 0.0},
      {1.0, 1e-16, 1.0, 1.0, 0.0}, {1.0, 1e12, 1.0, 1.0, 0.0},
      {0.0, 8.0, 1.0, 1.0, 0.0},   {1.0, 0.0, 1.0, 1.0, 0.0},
      {1.0, 8.0, 1.0, 1e8, 0.5},   {1.0, 8.0, 1.0, 1e-8, 0.5},
      {1.0, 8.0, 0.0, 1e-3, 0.5},  {1.0, 8.0, 1e-3, 0.0, 0.5}};

  for (std::size_t index = 0; index < sizeof(cases) / sizeof(cases[0]);
       ++index) {
    INFO("fallback case " << index);
    const hipdlp::algorithm::PrimalWeightUpdateResult result =
        hipdlp::algorithm::updatePrimalWeight(makeUpdateInput(
            history, 0.25, 1.0, cases[index].primal_distance,
            cases[index].dual_distance, cases[index].relative_primal,
            cases[index].relative_dual));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.outcome == hipdlp::algorithm::PrimalWeightOutcome::kFallback);
    CHECK(result.next.primal_weight == 3.0);
    CHECK(result.next.error_sum == 0.0);
    CHECK(result.next.previous_error == 0.0);
    CHECK(result.next.history_valid);
    CHECK(result.next.best_gap_valid);
    CHECK(result.next.best_primal_weight == 3.0);
    CHECK(result.next.best_residual_gap == cases[index].expected_best_gap);
    CHECK(result.base_step == 0.5);
    CHECK(result.primal_step == 0.5 / 3.0);
    CHECK(result.dual_step == 1.5);
  }

  // One ulp inside every window still performs the PID update.
  const double inside_cases[4][4] = {
      {std::nextafter(1e-16, 1.0), 8.0, 1.0, 1.0},
      {std::nextafter(1e12, 0.0), 8.0, 1.0, 1.0},
      {1.0, 8.0, 1.0, std::nextafter(1e8, 0.0)},
      {1.0, 8.0, 1.0, std::nextafter(1e-8, 1.0)}};
  for (std::size_t index = 0; index < 4; ++index) {
    INFO("inside case " << index);
    const hipdlp::algorithm::PrimalWeightUpdateResult result =
        hipdlp::algorithm::updatePrimalWeight(makeUpdateInput(
            history, 0.25, 1.0, inside_cases[index][0], inside_cases[index][1],
            inside_cases[index][2], inside_cases[index][3]));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.outcome == hipdlp::algorithm::PrimalWeightOutcome::kPidUpdate);
    CHECK(result.next.primal_weight != 3.0);
    CHECK(result.base_step == 0.5);
  }
}

TEST_CASE("hi-pdlp-primal-weight-invalid-input-is-transactional",
          "[pdlp][hipdlp][phase7f][stepsize][primalweight][errors]") {
  const double nan_value = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();

  checkInvalidUpdate(hipdlp::algorithm::updatePrimalWeight(
                         hipdlp::algorithm::PrimalWeightInput()),
                     hipdlp::api::ErrorCode::kInvalidArgument);

  const hipdlp::algorithm::PrimalWeightState valid =
      hipdlp::algorithm::resetPrimalWeight(2.0);
  const hipdlp::algorithm::PrimalWeightInput base =
      makeUpdateInput(valid, 0.25, 1.0, 1.0, 8.0, 1e-3, 1e-3);
  REQUIRE(hipdlp::algorithm::updatePrimalWeight(base).ok());

  hipdlp::algorithm::PrimalWeightInput input = base;
  input.parameters.k_p = nan_value;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidOption);
  input = base;
  input.parameters.k_p = -0.1;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidOption);
  input = base;
  input.parameters.k_i = infinity;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidOption);
  input = base;
  input.parameters.k_d = -1.0;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidOption);
  input = base;
  input.parameters.integral_smoothing = nan_value;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidOption);
  input = base;
  input.parameters.integral_smoothing = -0.1;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidOption);
  input = base;
  input.parameters.integral_smoothing = std::nextafter(1.0, 2.0);
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidOption);
  input = base;
  input.parameters.integral_smoothing = 0.0;
  CHECK(hipdlp::algorithm::updatePrimalWeight(input).ok());
  input.parameters.integral_smoothing = 1.0;
  CHECK(hipdlp::algorithm::updatePrimalWeight(input).ok());

  input = base;
  input.previous.history_valid = false;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  const double invalid_weights[4] = {0.0, -1.0, nan_value, infinity};
  for (std::size_t index = 0; index < 4; ++index) {
    input = base;
    input.previous.primal_weight = invalid_weights[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
    input = base;
    input.previous.best_primal_weight = invalid_weights[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  }
  input = base;
  input.previous.error_sum = nan_value;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  input = base;
  input.previous.previous_error = infinity;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  input = base;
  input.previous.best_residual_gap = 1.0;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  input = base;
  input.previous.best_gap_valid = true;
  input.previous.best_residual_gap = -1.0;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  input = base;
  input.previous.best_gap_valid = true;
  input.previous.best_residual_gap = nan_value;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);

  const double invalid_steps[4] = {0.0, -1.0, nan_value, infinity};
  for (std::size_t index = 0; index < 4; ++index) {
    input = base;
    input.primal_step = invalid_steps[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
    input = base;
    input.dual_step = invalid_steps[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  }
  const double invalid_magnitudes[3] = {-1.0, nan_value, infinity};
  for (std::size_t index = 0; index < 3; ++index) {
    input = base;
    input.primal_distance = invalid_magnitudes[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
    input = base;
    input.dual_distance = invalid_magnitudes[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
    input = base;
    input.relative_primal_residual = invalid_magnitudes[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
    input = base;
    input.relative_dual_residual = invalid_magnitudes[index];
    checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  }

  // The residual ratio itself must stay inside FP64 range.
  input = base;
  input.relative_primal_residual = 1e-300;
  input.relative_dual_residual = 1e300;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);

  // Exponential overflow and underflow.
  checkRejected(makeUpdateInput(hipdlp::algorithm::resetPrimalWeight(1e-308),
                                0.25, 1.0, 1e-15, 1e11, 1.0, 1.0),
                hipdlp::api::ErrorCode::kInvalidArgument);
  checkRejected(makeUpdateInput(hipdlp::algorithm::resetPrimalWeight(1e300),
                                0.25, 1.0, 1e11, 1e-15, 1.0, 1.0),
                hipdlp::api::ErrorCode::kInvalidArgument);

  // Updated-weight overflow and underflow.
  hipdlp::algorithm::PrimalWeightState integral_state =
      hipdlp::algorithm::resetPrimalWeight(1e300);
  integral_state.error_sum = 700.0;
  input = makeUpdateInput(integral_state, 0.25, 1.0, 1e-15, 1e11, 1.0, 1.0);
  input.parameters.k_p = 0.0;
  input.parameters.k_i = 1.0;
  input.parameters.integral_smoothing = 1.0;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);
  integral_state = hipdlp::algorithm::resetPrimalWeight(1e-300);
  integral_state.error_sum = -730.909;
  input = makeUpdateInput(integral_state, 0.25, 1.0, 1e11, 1e-15, 1.0, 1.0);
  input.parameters.k_p = 0.0;
  input.parameters.k_i = 1.0;
  input.parameters.integral_smoothing = 1.0;
  checkRejected(input, hipdlp::api::ErrorCode::kInvalidArgument);

  // The preserved geometric mean must stay inside FP64 range.
  checkRejected(makeUpdateInput(valid, 1e300, 1e300, 1.0, 8.0, 1e-3, 1e-3),
                hipdlp::api::ErrorCode::kInvalidArgument);
  checkRejected(makeUpdateInput(valid, 1e-320, 1e-320, 1.0, 8.0, 1e-3, 1e-3),
                hipdlp::api::ErrorCode::kInvalidArgument);

  // Derived step overflow and underflow, reached through the fallback weight.
  hipdlp::algorithm::PrimalWeightState tiny_best =
      hipdlp::algorithm::resetPrimalWeight(1.0);
  tiny_best.best_primal_weight = 1e-200;
  checkRejected(makeUpdateInput(tiny_best, 1e154, 1e154, 1e-16, 8.0, 1.0, 1.0),
                hipdlp::api::ErrorCode::kInvalidArgument);
  hipdlp::algorithm::PrimalWeightState huge_best =
      hipdlp::algorithm::resetPrimalWeight(1.0);
  huge_best.best_primal_weight = 1e200;
  checkRejected(makeUpdateInput(huge_best, 1e154, 1e154, 1e-16, 8.0, 1.0, 1.0),
                hipdlp::api::ErrorCode::kInvalidArgument);
  checkRejected(
      makeUpdateInput(huge_best, 1e-154, 1e-154, 1e-16, 8.0, 1.0, 1.0),
      hipdlp::api::ErrorCode::kInvalidArgument);
}

TEST_CASE("hi-pdlp-primal-weight-update-is-deterministic-and-input-immutable",
          "[pdlp][hipdlp][phase7f][stepsize][primalweight][determinism]") {
  hipdlp::algorithm::PrimalWeightState seed =
      hipdlp::algorithm::resetPrimalWeight(1.7);
  seed.error_sum = 0.31;
  seed.previous_error = -0.11;
  seed.best_primal_weight = 1.3;
  seed.best_residual_gap = 0.42;
  seed.best_gap_valid = true;

  const hipdlp::algorithm::PrimalWeightInput input =
      makeUpdateInput(seed, 0.37, 2.9, 1.9, 7.3, 3.1e-4, 8.7e-4);
  hipdlp::algorithm::PrimalWeightInput reused = input;
  const hipdlp::algorithm::PrimalWeightUpdateResult first =
      hipdlp::algorithm::updatePrimalWeight(reused);
  INFO(first.error.message);
  REQUIRE(first.ok());
  checkSameInput(reused, input);

  for (int repeat = 0; repeat < 4; ++repeat) {
    const hipdlp::algorithm::PrimalWeightUpdateResult repeated =
        hipdlp::algorithm::updatePrimalWeight(reused);
    REQUIRE(repeated.ok());
    checkSameResult(repeated, first);
    checkSameInput(reused, input);
  }

  const hipdlp::algorithm::PrimalWeightUpdateResult unrelated =
      hipdlp::algorithm::updatePrimalWeight(
          makeUpdateInput(hipdlp::algorithm::resetPrimalWeight(9.0), 1.0, 1.0,
                          4.0, 0.5, 1e-2, 1e-5));
  REQUIRE(unrelated.ok());
  const hipdlp::algorithm::PrimalWeightUpdateResult invalid =
      hipdlp::algorithm::updatePrimalWeight(
          hipdlp::algorithm::PrimalWeightInput());
  CHECK_FALSE(invalid.ok());

  const hipdlp::algorithm::PrimalWeightUpdateResult after =
      hipdlp::algorithm::updatePrimalWeight(input);
  REQUIRE(after.ok());
  checkSameResult(after, first);
}

namespace {

std::uint64_t doubleBits(const double value) {
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

hipdlp::algorithm::PrimalWeightUpdateInput makeRestartWeightInput(
    const hipdlp::algorithm::PrimalWeightState& previous,
    const double primal_anchor_distance, const double dual_anchor_distance,
    const double relative_primal_residual, const double relative_dual_residual,
    const double primal_step_size, const double dual_step_size) {
  hipdlp::algorithm::PrimalWeightUpdateInput input;
  input.previous = previous;
  input.primal_anchor_distance = primal_anchor_distance;
  input.dual_anchor_distance = dual_anchor_distance;
  input.relative_primal_residual = relative_primal_residual;
  input.relative_dual_residual = relative_dual_residual;
  input.primal_step_size = primal_step_size;
  input.dual_step_size = dual_step_size;
  return input;
}

hipdlp::algorithm::PrimalWeightState controlledReference(
    const hipdlp::algorithm::PrimalWeightState& previous,
    const double primal_anchor_distance, const double dual_anchor_distance,
    const hipdlp::algorithm::PrimalWeightParameters& parameters) {
  const double error = std::log(dual_anchor_distance) -
                       std::log(primal_anchor_distance) -
                       std::log(previous.primal_weight);
  const double error_sum =
      parameters.integral_smoothing * previous.error_sum + error;
  const double delta = error - previous.previous_error;
  const double exponent = parameters.k_p * error + parameters.k_i * error_sum +
                          parameters.k_d * delta;

  hipdlp::algorithm::PrimalWeightState next = previous;
  next.primal_weight *= std::exp(exponent);
  next.error_sum = error_sum;
  next.previous_error = error;
  return next;
}

void checkA12ZeroState(const hipdlp::algorithm::PrimalWeightState& state) {
  CHECK(doubleBits(state.primal_weight) == UINT64_C(0x0000000000000000));
  CHECK(doubleBits(state.error_sum) == UINT64_C(0x0000000000000000));
  CHECK(doubleBits(state.previous_error) == UINT64_C(0x0000000000000000));
  CHECK(doubleBits(state.best_primal_weight) == UINT64_C(0x0000000000000000));
  CHECK(doubleBits(state.best_residual_gap) == UINT64_C(0x0000000000000000));
  CHECK_FALSE(state.history_valid);
}

void checkA12Failure(const hipdlp::algorithm::PrimalWeightUpdateResult& result,
                     const hipdlp::api::ErrorCode expected_code,
                     const char* expected_message) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.error.code == expected_code);
  CHECK(result.error.message == expected_message);
  CHECK(result.kind == hipdlp::algorithm::PrimalWeightUpdateKind::kNone);
  checkA12ZeroState(result.next);
  CHECK(result.step_sizes.primal_weight == 0.0);
  CHECK(result.step_sizes.base_step == 0.0);
  CHECK(result.step_sizes.primal_step == 0.0);
  CHECK(result.step_sizes.dual_step == 0.0);
}

void checkA12StateBits(const hipdlp::algorithm::PrimalWeightState& actual,
                       const hipdlp::algorithm::PrimalWeightState& expected) {
  CHECK(doubleBits(actual.primal_weight) == doubleBits(expected.primal_weight));
  CHECK(doubleBits(actual.error_sum) == doubleBits(expected.error_sum));
  CHECK(doubleBits(actual.previous_error) ==
        doubleBits(expected.previous_error));
  CHECK(doubleBits(actual.best_primal_weight) ==
        doubleBits(expected.best_primal_weight));
  CHECK(doubleBits(actual.best_residual_gap) ==
        doubleBits(expected.best_residual_gap));
  CHECK(actual.history_valid == expected.history_valid);
}

}  // namespace

TEST_CASE("hi-pdlp-restart-primal-weight-pid-fixture-is-bitwise-exact",
          "[pdlp][hipdlp][phase7l][stepsize][primalweight][pid]") {
  const hipdlp::algorithm::PrimalWeightParameters defaults;
  CHECK(defaults.integral_smoothing == 0.3);
  CHECK(defaults.k_p == 0.99);
  CHECK(defaults.k_i == 0.01);
  CHECK(defaults.k_d == 0.0);

  hipdlp::algorithm::PrimalWeightState state =
      hipdlp::algorithm::resetPrimalWeight(1.7);
  state.error_sum = 0.31;
  state.previous_error = -0.11;
  state.best_primal_weight = 1.3;
  state.best_residual_gap = 0.6;

  const hipdlp::algorithm::PrimalWeightUpdateResult result =
      hipdlp::algorithm::updatePrimalWeightAtRestart(
          makeRestartWeightInput(state, 1.9, 7.3, 3.1e-4, 8.7e-4, 0.37, 2.9));
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.kind == hipdlp::algorithm::PrimalWeightUpdateKind::kControlled);

  // Hand-derived in the legacy evaluation order, then pinned as binary64.
  CHECK(doubleBits(result.next.previous_error) == UINT64_C(0x3fea17b167ea08ef));
  CHECK(doubleBits(result.next.error_sum) == UINT64_C(0x3fed118c8abaee4f));
  CHECK(doubleBits(result.next.primal_weight) == UINT64_C(0x400ec3f3eb5ddf2f));
  CHECK(doubleBits(result.next.best_residual_gap) ==
        UINT64_C(0x3fdcae9d0a9bc60c));
  CHECK(doubleBits(result.next.best_primal_weight) ==
        UINT64_C(0x400ec3f3eb5ddf2f));
  CHECK(doubleBits(result.step_sizes.primal_weight) ==
        UINT64_C(0x400ec3f3eb5ddf2f));
  CHECK(doubleBits(result.step_sizes.base_step) ==
        UINT64_C(0x3ff092deedd971ee));
  CHECK(doubleBits(result.step_sizes.primal_step) ==
        UINT64_C(0x3fd13d2116a2ef39));
  CHECK(doubleBits(result.step_sizes.dual_step) ==
        UINT64_C(0x400fde5ca552a336));
  CHECK(doubleBits(std::sqrt(0.37 * 2.9)) ==
        doubleBits(result.step_sizes.base_step));
}

TEST_CASE("hi-pdlp-restart-primal-weight-guard-boundaries-are-strict",
          "[pdlp][hipdlp][phase7l][stepsize][primalweight][guards]") {
  hipdlp::algorithm::PrimalWeightState state =
      hipdlp::algorithm::resetPrimalWeight(5.0);
  state.error_sum = 1.25;
  state.previous_error = -0.75;
  state.best_primal_weight = 3.0;
  state.best_residual_gap = 0.0;

  struct GuardCase {
    double primal_distance;
    double dual_distance;
    double relative_primal;
    double relative_dual;
    hipdlp::algorithm::PrimalWeightUpdateKind expected;
  };
  const GuardCase cases[] = {
      {1e-16, 8.0, 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest},
      {std::nextafter(1e-16, 1.0), 8.0, 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kControlled},
      {1.0, 1e-16, 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest},
      {1.0, std::nextafter(1e-16, 1.0), 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kControlled},
      {1e12, 8.0, 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest},
      {std::nextafter(1e12, 0.0), 8.0, 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kControlled},
      {1.0, 1e12, 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest},
      {1.0, std::nextafter(1e12, 0.0), 1.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kControlled},
      {1.0, 8.0, 1.0, 1e-8,
       hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest},
      {1.0, 8.0, 1.0, std::nextafter(1e-8, 1.0),
       hipdlp::algorithm::PrimalWeightUpdateKind::kControlled},
      {1.0, 8.0, 1.0, 1e8,
       hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest},
      {1.0, 8.0, 1.0, std::nextafter(1e8, 0.0),
       hipdlp::algorithm::PrimalWeightUpdateKind::kControlled},
      {1.0, 8.0, 0.0, 1.0,
       hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest}};

  for (std::size_t index = 0; index < sizeof(cases) / sizeof(cases[0]);
       ++index) {
    INFO("guard case " << index);
    const hipdlp::algorithm::PrimalWeightUpdateResult result =
        hipdlp::algorithm::updatePrimalWeightAtRestart(makeRestartWeightInput(
            state, cases[index].primal_distance, cases[index].dual_distance,
            cases[index].relative_primal, cases[index].relative_dual, 0.25,
            1.0));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.kind == cases[index].expected);
    if (cases[index].expected ==
        hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest) {
      CHECK(doubleBits(result.next.primal_weight) == doubleBits(3.0));
      CHECK(result.next.error_sum == 0.0);
      CHECK(result.next.previous_error == 0.0);
    } else {
      CHECK_FALSE(bitwiseEqual(result.next.primal_weight, 3.0));
    }
  }
}

TEST_CASE(
    "hi-pdlp-restart-primal-weight-reversion-restores-best-and-resets-pid",
    "[pdlp][hipdlp][phase7l][stepsize][primalweight][reversion]") {
  hipdlp::algorithm::PrimalWeightState state =
      hipdlp::algorithm::resetPrimalWeight(5.0);
  state.error_sum = 1.25;
  state.previous_error = -0.75;
  state.best_primal_weight = 3.0;
  state.best_residual_gap = 0.2;

  const hipdlp::algorithm::PrimalWeightUpdateResult result =
      hipdlp::algorithm::updatePrimalWeightAtRestart(
          makeRestartWeightInput(state, 1e-16, 8.0, 1.0, 10.0, 0.25, 1.0));
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.kind ==
        hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest);
  CHECK(doubleBits(result.next.primal_weight) == doubleBits(3.0));
  CHECK(doubleBits(result.next.error_sum) == UINT64_C(0x0000000000000000));
  CHECK(doubleBits(result.next.previous_error) == UINT64_C(0x0000000000000000));
  CHECK(doubleBits(result.next.best_primal_weight) == doubleBits(3.0));
  CHECK(doubleBits(result.next.best_residual_gap) == doubleBits(0.2));
  CHECK(doubleBits(result.step_sizes.base_step) == doubleBits(0.5));
  CHECK(doubleBits(result.step_sizes.primal_weight) == doubleBits(3.0));
  CHECK(doubleBits(result.step_sizes.primal_step) == doubleBits(0.5 / 3.0));
  CHECK(doubleBits(result.step_sizes.dual_step) == doubleBits(1.5));
}

TEST_CASE("hi-pdlp-restart-primal-weight-best-tracking-is-strict",
          "[pdlp][hipdlp][phase7l][stepsize][primalweight][best]") {
  hipdlp::algorithm::PrimalWeightState state =
      hipdlp::algorithm::resetPrimalWeight(2.0);
  state.best_primal_weight = 1.25;
  state.best_residual_gap = 0.5;

  const hipdlp::algorithm::PrimalWeightUpdateResult improved =
      hipdlp::algorithm::updatePrimalWeightAtRestart(
          makeRestartWeightInput(state, 1.0, 4.0, 1.0, 2.0, 0.25, 1.0));
  INFO(improved.error.message);
  REQUIRE(improved.ok());
  REQUIRE(improved.kind ==
          hipdlp::algorithm::PrimalWeightUpdateKind::kControlled);
  CHECK(improved.next.best_residual_gap < 0.5);
  CHECK(doubleBits(improved.next.best_primal_weight) ==
        doubleBits(improved.next.primal_weight));

  const double improved_best_weight = improved.next.best_primal_weight;
  const double improved_best_gap = improved.next.best_residual_gap;
  const hipdlp::algorithm::PrimalWeightUpdateResult equal =
      hipdlp::algorithm::updatePrimalWeightAtRestart(makeRestartWeightInput(
          improved.next, 1.0, improved.next.primal_weight, 1.0, 2.0,
          improved.step_sizes.primal_step, improved.step_sizes.dual_step));
  INFO(equal.error.message);
  REQUIRE(equal.ok());
  CHECK(equal.kind == hipdlp::algorithm::PrimalWeightUpdateKind::kControlled);
  CHECK(doubleBits(equal.next.best_residual_gap) ==
        doubleBits(improved_best_gap));
  CHECK(doubleBits(equal.next.best_primal_weight) ==
        doubleBits(improved_best_weight));

  const hipdlp::algorithm::PrimalWeightUpdateResult zero_residual =
      hipdlp::algorithm::updatePrimalWeightAtRestart(makeRestartWeightInput(
          equal.next, 1.0, 4.0, 0.0, 2.0, equal.step_sizes.primal_step,
          equal.step_sizes.dual_step));
  INFO(zero_residual.error.message);
  REQUIRE(zero_residual.ok());
  CHECK(zero_residual.kind ==
        hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest);
  CHECK(doubleBits(zero_residual.next.best_residual_gap) ==
        doubleBits(improved_best_gap));
  CHECK(doubleBits(zero_residual.next.best_primal_weight) ==
        doubleBits(improved_best_weight));
}

TEST_CASE("hi-pdlp-restart-primal-weight-chaining-and-reversion-reset-history",
          "[pdlp][hipdlp][phase7l][stepsize][primalweight][chaining]") {
  hipdlp::algorithm::PrimalWeightState state =
      hipdlp::algorithm::resetPrimalWeight(2.0);
  state.error_sum = 0.4;
  state.previous_error = -0.2;
  state.best_primal_weight = 1.5;
  state.best_residual_gap = 0.0;
  hipdlp::algorithm::PrimalWeightState reference = state;
  double primal_step = 0.25;
  double dual_step = 1.0;
  const double primal_distances[] = {1.0, 2.0, 0.5};
  const double dual_distances[] = {8.0, 1.0, 3.0};

  hipdlp::algorithm::PrimalWeightParameters parameters;
  parameters.integral_smoothing = 0.5;
  parameters.k_p = 0.4;
  parameters.k_i = 0.2;
  parameters.k_d = 0.3;
  for (std::size_t index = 0; index < 3; ++index) {
    INFO("controlled update " << index);
    hipdlp::algorithm::PrimalWeightUpdateInput input = makeRestartWeightInput(
        state, primal_distances[index], dual_distances[index], 1.0, 1.0,
        primal_step, dual_step);
    input.parameters = parameters;
    const hipdlp::algorithm::PrimalWeightUpdateResult result =
        hipdlp::algorithm::updatePrimalWeightAtRestart(input);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.kind ==
            hipdlp::algorithm::PrimalWeightUpdateKind::kControlled);
    reference = controlledReference(reference, primal_distances[index],
                                    dual_distances[index], parameters);
    CHECK(doubleBits(result.next.primal_weight) ==
          doubleBits(reference.primal_weight));
    CHECK(doubleBits(result.next.error_sum) == doubleBits(reference.error_sum));
    CHECK(doubleBits(result.next.previous_error) ==
          doubleBits(reference.previous_error));
    state = result.next;
    primal_step = result.step_sizes.primal_step;
    dual_step = result.step_sizes.dual_step;
  }

  hipdlp::algorithm::PrimalWeightState reset_chain =
      hipdlp::algorithm::resetPrimalWeight(2.0);
  reset_chain.best_primal_weight = 1.25;
  reset_chain.best_residual_gap = 0.0;
  hipdlp::algorithm::PrimalWeightUpdateInput first_input =
      makeRestartWeightInput(reset_chain, 1.0, 8.0, 1.0, 1.0, 0.25, 1.0);
  first_input.parameters = parameters;
  const hipdlp::algorithm::PrimalWeightUpdateResult first =
      hipdlp::algorithm::updatePrimalWeightAtRestart(first_input);
  REQUIRE(first.ok());

  hipdlp::algorithm::PrimalWeightUpdateInput revert_input =
      makeRestartWeightInput(first.next, 1e-16, 8.0, 1.0, 1.0,
                             first.step_sizes.primal_step,
                             first.step_sizes.dual_step);
  revert_input.parameters = parameters;
  const hipdlp::algorithm::PrimalWeightUpdateResult reverted =
      hipdlp::algorithm::updatePrimalWeightAtRestart(revert_input);
  REQUIRE(reverted.ok());
  REQUIRE(reverted.kind ==
          hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest);
  CHECK(reverted.next.error_sum == 0.0);
  CHECK(reverted.next.previous_error == 0.0);
  CHECK(doubleBits(reverted.next.primal_weight) == doubleBits(1.25));

  hipdlp::algorithm::PrimalWeightUpdateInput after_input =
      makeRestartWeightInput(reverted.next, 1.0, 5.0, 1.0, 1.0,
                             reverted.step_sizes.primal_step,
                             reverted.step_sizes.dual_step);
  after_input.parameters = parameters;
  const hipdlp::algorithm::PrimalWeightUpdateResult after =
      hipdlp::algorithm::updatePrimalWeightAtRestart(after_input);
  REQUIRE(after.ok());

  hipdlp::algorithm::PrimalWeightState fresh_state =
      hipdlp::algorithm::resetPrimalWeight(1.25);
  fresh_state.best_residual_gap = 0.0;
  hipdlp::algorithm::PrimalWeightUpdateInput fresh_input =
      makeRestartWeightInput(fresh_state, 1.0, 5.0, 1.0, 1.0,
                             reverted.step_sizes.primal_step,
                             reverted.step_sizes.dual_step);
  fresh_input.parameters = parameters;
  const hipdlp::algorithm::PrimalWeightUpdateResult fresh =
      hipdlp::algorithm::updatePrimalWeightAtRestart(fresh_input);
  REQUIRE(fresh.ok());
  checkA12StateBits(after.next, fresh.next);
  CHECK(doubleBits(after.step_sizes.primal_step) ==
        doubleBits(fresh.step_sizes.primal_step));
  CHECK(doubleBits(after.step_sizes.dual_step) ==
        doubleBits(fresh.step_sizes.dual_step));
}

TEST_CASE("hi-pdlp-restart-primal-weight-fp64-errors-are-transactional",
          "[pdlp][hipdlp][phase7l][stepsize][primalweight][errors]") {
  const double nan_value = std::numeric_limits<double>::quiet_NaN();
  const double maximum = std::numeric_limits<double>::max();
  hipdlp::algorithm::PrimalWeightUpdateInput input = makeRestartWeightInput(
      hipdlp::algorithm::resetPrimalWeight(2.0), 1.0, 8.0, 1.0, 1.0, 0.25, 1.0);

  input.parameters.k_p = nan_value;
  input.previous.primal_weight = 0.0;
  checkA12Failure(hipdlp::algorithm::updatePrimalWeightAtRestart(input),
                  hipdlp::api::ErrorCode::kInvalidOption,
                  "HiPDLP primal-weight PID parameters are invalid");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(2.0), 1.0,
                                 8.0, 1.0, 1.0, 0.25, 1.0);
  input.previous.primal_weight = 0.0;
  checkA12Failure(hipdlp::algorithm::updatePrimalWeightAtRestart(input),
                  hipdlp::api::ErrorCode::kInvalidArgument,
                  "HiPDLP primal-weight state is invalid");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(2.0), 1.0,
                                 8.0, 1.0, 1.0, 0.25, 1.0);
  input.previous.primal_weight = -1.0;
  checkA12Failure(hipdlp::algorithm::updatePrimalWeightAtRestart(input),
                  hipdlp::api::ErrorCode::kInvalidArgument,
                  "HiPDLP primal-weight state is invalid");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(2.0),
                                 nan_value, 8.0, 1.0, 1.0, 0.25, 1.0);
  checkA12Failure(
      hipdlp::algorithm::updatePrimalWeightAtRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument,
      "HiPDLP primal-weight anchor distances must be finite and nonnegative");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(1.0),
                                 1e-15, 1e11, 1.0, 1.0, 0.25, 1.0);
  input.parameters.k_p = maximum;
  input.parameters.k_i = 0.0;
  checkA12Failure(
      hipdlp::algorithm::updatePrimalWeightAtRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument,
      "HiPDLP primal-weight PID proportional term is outside FP64 range");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(1e-300),
                                 1e-15, 1e11, 1.0, 1.0, 0.25, 1.0);
  checkA12Failure(hipdlp::algorithm::updatePrimalWeightAtRestart(input),
                  hipdlp::api::ErrorCode::kInvalidArgument,
                  "HiPDLP primal-weight PID factor overflowed or underflowed");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(1e300),
                                 1e11, 1e-15, 1.0, 1.0, 0.25, 1.0);
  checkA12Failure(hipdlp::algorithm::updatePrimalWeightAtRestart(input),
                  hipdlp::api::ErrorCode::kInvalidArgument,
                  "HiPDLP primal-weight PID factor overflowed or underflowed");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(2.0), 1.0,
                                 8.0, 1.0, 1.0, maximum, 2.0);
  checkA12Failure(
      hipdlp::algorithm::updatePrimalWeightAtRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument,
      "HiPDLP primal-weight step-size product is outside FP64 range");

  input = makeRestartWeightInput(hipdlp::algorithm::resetPrimalWeight(2.0), 1.0,
                                 8.0, 1.0, 1.0,
                                 std::numeric_limits<double>::denorm_min(),
                                 std::numeric_limits<double>::denorm_min());
  checkA12Failure(
      hipdlp::algorithm::updatePrimalWeightAtRestart(input),
      hipdlp::api::ErrorCode::kInvalidArgument,
      "HiPDLP primal-weight step-size product is outside FP64 range");
}

TEST_CASE("hi-pdlp-restart-primal-weight-is-bitwise-deterministic",
          "[pdlp][hipdlp][phase7l][stepsize][primalweight][determinism]") {
  hipdlp::algorithm::PrimalWeightState state =
      hipdlp::algorithm::resetPrimalWeight(2.0);
  state.error_sum = 0.25;
  state.previous_error = -0.5;
  state.best_primal_weight = 1.5;
  state.best_residual_gap = 0.5;
  const hipdlp::algorithm::PrimalWeightUpdateInput input =
      makeRestartWeightInput(state, 2.0, 12.0, 1e-3, 1e-2, 0.125, 2.0);

  const hipdlp::algorithm::PrimalWeightUpdateResult first =
      hipdlp::algorithm::updatePrimalWeightAtRestart(input);
  INFO(first.error.message);
  REQUIRE(first.ok());
  CHECK(doubleBits(first.next.previous_error) == UINT64_C(0x3ff193ea7aad030a));
  CHECK(doubleBits(first.next.error_sum) == UINT64_C(0x3ff2c71dade0363d));
  CHECK(doubleBits(first.next.primal_weight) == UINT64_C(0x4018049c17299b39));
  CHECK(doubleBits(first.step_sizes.base_step) == UINT64_C(0x3fe0000000000000));
  CHECK(doubleBits(first.step_sizes.primal_step) ==
        UINT64_C(0x3fb5513d2683f5bb));
  CHECK(doubleBits(first.step_sizes.dual_step) == UINT64_C(0x4008049c17299b39));

  for (int repeat = 0; repeat < 8; ++repeat) {
    const hipdlp::algorithm::PrimalWeightUpdateResult repeated =
        hipdlp::algorithm::updatePrimalWeightAtRestart(input);
    REQUIRE(repeated.ok());
    CHECK(repeated.kind == first.kind);
    checkA12StateBits(repeated.next, first.next);
    CHECK(doubleBits(repeated.step_sizes.primal_weight) ==
          doubleBits(first.step_sizes.primal_weight));
    CHECK(doubleBits(repeated.step_sizes.base_step) ==
          doubleBits(first.step_sizes.base_step));
    CHECK(doubleBits(repeated.step_sizes.primal_step) ==
          doubleBits(first.step_sizes.primal_step));
    CHECK(doubleBits(repeated.step_sizes.dual_step) ==
          doubleBits(first.step_sizes.dual_step));
  }
}
