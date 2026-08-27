#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "pdlp/hipdlp/algorithm/convergence.hpp"

namespace {

hipdlp::algorithm::StoppingInput makeValidInput() {
  hipdlp::algorithm::StoppingInput input;
  input.raw.primal_residual = 0.0;
  input.raw.dual_residual = 0.0;
  input.raw.objective_gap = 0.0;
  input.raw.primal_objective = 0.0;
  input.raw.dual_objective = 0.0;
  input.raw.unscaled_rhs_norm = 0.0;
  input.raw.unscaled_cost_norm = 0.0;
  input.tolerances.primal = 0.25;
  input.tolerances.dual = 0.25;
  input.tolerances.gap = 0.25;
  input.iteration = 0;
  input.elapsed_time = 0.0;
  input.limits.iteration = kHighsIInf;
  input.limits.time = kHighsInf;
  return input;
}

hipdlp::algorithm::StoppingInput makeNonOptimalInput() {
  hipdlp::algorithm::StoppingInput input = makeValidInput();
  input.raw.primal_residual = 1.0;
  return input;
}

bool bitwiseEqual(const double left, const double right) {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

void checkNumericalError(const hipdlp::algorithm::StoppingInput& input) {
  const hipdlp::algorithm::StoppingResult result =
      hipdlp::algorithm::evaluateStopping(input);
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK_FALSE(result.metrics_valid);
  CHECK(result.decision ==
        hipdlp::algorithm::StoppingDecision::kNumericalError);
  CHECK(result.metrics.relative_primal == 0.0);
  CHECK(result.metrics.relative_dual == 0.0);
  CHECK(result.metrics.relative_gap == 0.0);
}

std::string readSourceFile(const char* relative_path) {
  const std::string path = std::string(HIGHS_DIR) + "/" + relative_path;
  std::ifstream stream(path.c_str(), std::ios::in | std::ios::binary);
  INFO(path);
  REQUIRE(stream.good());
  std::ostringstream contents;
  contents << stream.rdbuf();
  return contents.str();
}

}  // namespace

TEST_CASE("hi-pdlp-convergence-policy-normalizes-exact-scalars",
          "[pdlp][hipdlp][phase7a][convergence][policy][normalization]") {
  hipdlp::algorithm::StoppingInput input = makeNonOptimalInput();
  input.raw.primal_residual = 6.0;
  input.raw.dual_residual = 16.0;
  input.raw.objective_gap = 20.0;
  input.raw.primal_objective = 3.0;
  input.raw.dual_objective = -4.0;
  input.raw.unscaled_rhs_norm = 3.0;
  input.raw.unscaled_cost_norm = 7.0;

  const hipdlp::algorithm::StoppingResult result =
      hipdlp::algorithm::evaluateStopping(input);
  REQUIRE(result.ok());
  REQUIRE(result.metrics_valid);
  CHECK(result.metrics.relative_primal == 1.5);
  CHECK(result.metrics.relative_dual == 2.0);
  CHECK(result.metrics.relative_gap == 2.5);
  CHECK(result.decision == hipdlp::algorithm::StoppingDecision::kContinue);
}

TEST_CASE("hi-pdlp-convergence-policy-requires-every-threshold",
          "[pdlp][hipdlp][phase7a][convergence][policy][thresholds]") {
  hipdlp::algorithm::StoppingInput input = makeValidInput();
  input.raw.primal_residual = 1.0;
  input.raw.dual_residual = 2.0;
  input.raw.objective_gap = 2.25;
  input.raw.primal_objective = 4.0;
  input.raw.dual_objective = 4.0;
  input.raw.unscaled_rhs_norm = 3.0;
  input.raw.unscaled_cost_norm = 7.0;

  hipdlp::algorithm::StoppingResult result =
      hipdlp::algorithm::evaluateStopping(input);
  REQUIRE(result.ok());
  CHECK(result.metrics.relative_primal == input.tolerances.primal);
  CHECK(result.metrics.relative_dual == input.tolerances.dual);
  CHECK(result.metrics.relative_gap == input.tolerances.gap);
  CHECK(result.decision == hipdlp::algorithm::StoppingDecision::kOptimal);

  input.raw.primal_residual = 2.0;
  result = hipdlp::algorithm::evaluateStopping(input);
  REQUIRE(result.ok());
  CHECK(result.metrics.relative_primal > input.tolerances.primal);
  CHECK(result.metrics.relative_dual <= input.tolerances.dual);
  CHECK(result.metrics.relative_gap <= input.tolerances.gap);
  CHECK(result.decision == hipdlp::algorithm::StoppingDecision::kContinue);

  input.raw.primal_residual = 1.0;
  input.raw.dual_residual = 4.0;
  result = hipdlp::algorithm::evaluateStopping(input);
  REQUIRE(result.ok());
  CHECK(result.metrics.relative_primal <= input.tolerances.primal);
  CHECK(result.metrics.relative_dual > input.tolerances.dual);
  CHECK(result.metrics.relative_gap <= input.tolerances.gap);
  CHECK(result.decision == hipdlp::algorithm::StoppingDecision::kContinue);

  input.raw.dual_residual = 2.0;
  input.raw.objective_gap = 4.5;
  result = hipdlp::algorithm::evaluateStopping(input);
  REQUIRE(result.ok());
  CHECK(result.metrics.relative_primal <= input.tolerances.primal);
  CHECK(result.metrics.relative_dual <= input.tolerances.dual);
  CHECK(result.metrics.relative_gap > input.tolerances.gap);
  CHECK(result.decision == hipdlp::algorithm::StoppingDecision::kContinue);
}

TEST_CASE("hi-pdlp-convergence-policy-handles-zero-scale-objectives",
          "[pdlp][hipdlp][phase7a][convergence][policy][normalization]") {
  hipdlp::algorithm::StoppingInput zero = makeNonOptimalInput();
  zero.raw.objective_gap = 0.5;
  const hipdlp::algorithm::StoppingResult zero_result =
      hipdlp::algorithm::evaluateStopping(zero);
  REQUIRE(zero_result.ok());
  CHECK(zero_result.metrics.relative_gap == 0.5);

  hipdlp::algorithm::StoppingInput near_zero = zero;
  near_zero.raw.primal_objective = std::numeric_limits<double>::denorm_min();
  near_zero.raw.dual_objective = -std::numeric_limits<double>::denorm_min();
  const hipdlp::algorithm::StoppingResult near_zero_result =
      hipdlp::algorithm::evaluateStopping(near_zero);
  REQUIRE(near_zero_result.ok());
  CHECK(near_zero_result.metrics.relative_gap == 0.5);
}

TEST_CASE("hi-pdlp-convergence-policy-rejects-invalid-scalars",
          "[pdlp][hipdlp][phase7a][convergence][policy][errors]") {
  typedef double hipdlp::algorithm::RawConvergenceMetrics::* RawField;
  const RawField nonnegative_fields[] = {
      &hipdlp::algorithm::RawConvergenceMetrics::primal_residual,
      &hipdlp::algorithm::RawConvergenceMetrics::dual_residual,
      &hipdlp::algorithm::RawConvergenceMetrics::objective_gap,
      &hipdlp::algorithm::RawConvergenceMetrics::unscaled_rhs_norm,
      &hipdlp::algorithm::RawConvergenceMetrics::unscaled_cost_norm};
  for (const RawField field : nonnegative_fields) {
    hipdlp::algorithm::StoppingInput input = makeValidInput();
    input.raw.*field = -1.0;
    checkNumericalError(input);
  }

  const RawField finite_fields[] = {
      &hipdlp::algorithm::RawConvergenceMetrics::primal_residual,
      &hipdlp::algorithm::RawConvergenceMetrics::dual_residual,
      &hipdlp::algorithm::RawConvergenceMetrics::objective_gap,
      &hipdlp::algorithm::RawConvergenceMetrics::primal_objective,
      &hipdlp::algorithm::RawConvergenceMetrics::dual_objective,
      &hipdlp::algorithm::RawConvergenceMetrics::unscaled_rhs_norm,
      &hipdlp::algorithm::RawConvergenceMetrics::unscaled_cost_norm};
  const double nonfinite_values[] = {std::numeric_limits<double>::quiet_NaN(),
                                     std::numeric_limits<double>::infinity()};
  for (const RawField field : finite_fields) {
    for (const double value : nonfinite_values) {
      hipdlp::algorithm::StoppingInput input = makeValidInput();
      input.raw.*field = value;
      checkNumericalError(input);
    }
  }

  typedef double hipdlp::algorithm::ConvergenceTolerances::* ToleranceField;
  const ToleranceField tolerance_fields[] = {
      &hipdlp::algorithm::ConvergenceTolerances::primal,
      &hipdlp::algorithm::ConvergenceTolerances::dual,
      &hipdlp::algorithm::ConvergenceTolerances::gap};
  const double invalid_tolerances[] = {0.0, -1.0,
                                       std::numeric_limits<double>::quiet_NaN(),
                                       std::numeric_limits<double>::infinity()};
  for (const ToleranceField field : tolerance_fields) {
    for (const double value : invalid_tolerances) {
      hipdlp::algorithm::StoppingInput input = makeValidInput();
      input.tolerances.*field = value;
      checkNumericalError(input);
    }
  }

  hipdlp::algorithm::StoppingInput input = makeValidInput();
  input.iteration = -1;
  checkNumericalError(input);
  input = makeValidInput();
  input.elapsed_time = -1.0;
  checkNumericalError(input);
  input = makeValidInput();
  input.elapsed_time = std::numeric_limits<double>::quiet_NaN();
  checkNumericalError(input);
  input = makeValidInput();
  input.elapsed_time = std::numeric_limits<double>::infinity();
  checkNumericalError(input);
  input = makeValidInput();
  input.limits.iteration = -1;
  checkNumericalError(input);
  input = makeValidInput();
  input.limits.time = -1.0;
  checkNumericalError(input);
  input = makeValidInput();
  input.limits.time = std::numeric_limits<double>::quiet_NaN();
  checkNumericalError(input);
  input = makeValidInput();
  input.limits.time = -std::numeric_limits<double>::infinity();
  checkNumericalError(input);

  input = makeValidInput();
  input.raw.primal_objective = std::numeric_limits<double>::max();
  input.raw.dual_objective = std::numeric_limits<double>::max();
  checkNumericalError(input);
}

TEST_CASE("hi-pdlp-convergence-policy-uses-inclusive-limits",
          "[pdlp][hipdlp][phase7a][convergence][policy][limits]") {
  hipdlp::algorithm::StoppingInput input = makeNonOptimalInput();
  input.limits.iteration = 5;
  input.iteration = 4;
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kContinue);
  input.iteration = 5;
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kIterationLimit);

  input = makeNonOptimalInput();
  input.limits.time = 2.0;
  input.elapsed_time = 1.5;
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kContinue);
  input.elapsed_time = 2.0;
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kTimeLimit);

  input = makeNonOptimalInput();
  input.iteration = kHighsIInf;
  input.elapsed_time = std::numeric_limits<double>::max();
  REQUIRE(input.limits.iteration == kHighsIInf);
  REQUIRE(input.limits.time == kHighsInf);
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kContinue);
}

TEST_CASE("hi-pdlp-convergence-policy-applies-decision-precedence",
          "[pdlp][hipdlp][phase7a][convergence][policy][precedence]") {
  hipdlp::algorithm::StoppingInput input = makeValidInput();
  input.iteration = 3;
  input.elapsed_time = 4.0;
  input.limits.iteration = 3;
  input.limits.time = 4.0;
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kOptimal);

  input = makeNonOptimalInput();
  input.iteration = 3;
  input.elapsed_time = 4.0;
  input.limits.iteration = 3;
  input.limits.time = 4.0;
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kIterationLimit);

  input.limits.iteration = 4;
  CHECK(hipdlp::algorithm::evaluateStopping(input).decision ==
        hipdlp::algorithm::StoppingDecision::kTimeLimit);

  input.raw.primal_residual = std::numeric_limits<double>::quiet_NaN();
  checkNumericalError(input);
}

TEST_CASE("hi-pdlp-convergence-policy-is-bitwise-deterministic",
          "[pdlp][hipdlp][phase7a][convergence][policy][determinism]") {
  std::vector<hipdlp::algorithm::StoppingInput> sequence(5,
                                                         makeNonOptimalInput());
  sequence[0].raw.primal_residual = 1.0;
  sequence[1].raw.primal_residual = 0.5;
  sequence[2].raw.primal_residual = 0.25;
  sequence[3].iteration = 3;
  sequence[3].limits.iteration = 3;
  sequence[4].iteration = 2;
  sequence[4].elapsed_time = 4.0;
  sequence[4].limits.iteration = 3;
  sequence[4].limits.time = 4.0;
  const hipdlp::algorithm::StoppingDecision expected[] = {
      hipdlp::algorithm::StoppingDecision::kContinue,
      hipdlp::algorithm::StoppingDecision::kContinue,
      hipdlp::algorithm::StoppingDecision::kOptimal,
      hipdlp::algorithm::StoppingDecision::kIterationLimit,
      hipdlp::algorithm::StoppingDecision::kTimeLimit};

  for (std::size_t point = 0; point < sequence.size(); ++point) {
    CAPTURE(point);
    const hipdlp::algorithm::StoppingResult first =
        hipdlp::algorithm::evaluateStopping(sequence[point]);
    const hipdlp::algorithm::StoppingResult second =
        hipdlp::algorithm::evaluateStopping(sequence[point]);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    CHECK(first.decision == expected[point]);
    CHECK(second.decision == first.decision);
    CHECK(first.metrics_valid == second.metrics_valid);
    CHECK(first.error.code == second.error.code);
    CHECK(first.error.message == second.error.message);
    CHECK(bitwiseEqual(first.metrics.relative_primal,
                       second.metrics.relative_primal));
    CHECK(bitwiseEqual(first.metrics.relative_dual,
                       second.metrics.relative_dual));
    CHECK(
        bitwiseEqual(first.metrics.relative_gap, second.metrics.relative_gap));
  }
}

TEST_CASE("hi-pdlp-convergence-policy-has-no-resource-surface",
          "[pdlp][hipdlp][phase7a][convergence][policy][structure]") {
  const char* policy_files[] = {"highs/pdlp/hipdlp/algorithm/convergence.hpp",
                                "highs/pdlp/hipdlp/algorithm/convergence.cc"};
  const char* forbidden_terms[] = {"backend/",
                                   "backend::",
                                   "ProblemHandle",
                                   "VectorHandle",
                                   "WorkspaceHandle",
                                   "CommandBufferHandle",
                                   "DeviceInfo",
                                   "Resource",
                                   "cuda",
                                   "CUDA",
                                   "rocm",
                                   "ROCm",
                                   "#include <hip/",
                                   "hip_runtime",
                                   "hipblas",
                                   "hipsparse",
                                   "synchronize",
                                   "malloc(",
                                   "calloc(",
                                   "realloc(",
                                   "new ",
                                   "verifyKkt",
                                   "kkt_verifier"};
  for (const char* policy_file : policy_files) {
    CAPTURE(policy_file);
    const std::string contents = readSourceFile(policy_file);
    for (const char* forbidden_term : forbidden_terms) {
      CAPTURE(forbidden_term);
      CHECK(contents.find(forbidden_term) == std::string::npos);
    }
  }

  const hipdlp::algorithm::StoppingResult result =
      hipdlp::algorithm::evaluateStopping(makeValidInput());
  REQUIRE(result.ok());
  CHECK(result.decision == hipdlp::algorithm::StoppingDecision::kOptimal);
}