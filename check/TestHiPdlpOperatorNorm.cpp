#include <cstring>
#include <limits>
#include <memory>
#include <vector>

#include "HConfig.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "pdlp/hipdlp/algorithm/operator_norm.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/model/problem.hpp"

namespace {

HighsLp makeEqualityLp(const HighsInt rows, const HighsInt cols,
                       const std::vector<HighsInt>& starts,
                       const std::vector<HighsInt>& indices,
                       const std::vector<double>& values) {
  HighsLp lp;
  lp.num_row_ = rows;
  lp.num_col_ = cols;
  lp.col_cost_.assign(cols, 0.0);
  lp.col_lower_.assign(cols, -kHighsInf);
  lp.col_upper_.assign(cols, kHighsInf);
  lp.row_lower_.assign(rows, 0.0);
  lp.row_upper_.assign(rows, 0.0);
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_row_ = rows;
  lp.a_matrix_.num_col_ = cols;
  lp.a_matrix_.start_ = starts;
  lp.a_matrix_.index_ = indices;
  lp.a_matrix_.value_ = values;
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0.0;
  return lp;
}

HighsLp makeExactLp() { return makeEqualityLp(1, 1, {0, 1}, {0}, {3.0}); }

HighsLp makeConvergingLp() {
  return makeEqualityLp(2, 2, {0, 1, 2}, {0, 1}, {2.0, 1.0});
}

HighsLp makeOrthogonalDominantDirectionLp() {
    return makeEqualityLp(2, 2, {0, 2, 4}, {0, 1, 0, 1},
                                                {2.0, -2.0, 1.0, 1.0});
}

hipdlp::algorithm::OperatorNormResult runEstimate(
    const hipdlp::backend::Kind kind, const HighsLp& lp,
    const hipdlp::algorithm::OperatorNormParameters& parameters =
        hipdlp::algorithm::OperatorNormParameters()) {
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  INFO(built.error.message);
  REQUIRE(built.ok());
  hipdlp::backend::FactoryResult factory = hipdlp::backend::createBackend(kind);
  INFO(factory.error.message);
  REQUIRE(factory.ok());
  hipdlp::backend::ProblemHandle problem;
  REQUIRE(factory.backend->bindProblem(built.problem, problem).ok());
  return hipdlp::algorithm::estimateOperatorNormSquared(
      *factory.backend, problem, built.problem.numRows(),
      built.problem.numCols(), parameters);
}

bool bitwiseEqual(const double left, const double right) {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

void checkZeroEstimate(const hipdlp::algorithm::OperatorNormResult& result,
                       const hipdlp::api::ErrorCode expected_code) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.error.code == expected_code);
  CHECK(result.operator_norm_squared == 0.0);
  CHECK(result.performed_iterations == 0);
  CHECK(result.final_relative_change == 0.0);
  CHECK_FALSE(result.tolerance_met);
  CHECK_FALSE(result.empty_matrix);
}

void checkEmptyEstimate(const hipdlp::algorithm::OperatorNormResult& result) {
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.error.code == hipdlp::api::ErrorCode::kNone);
  CHECK(result.operator_norm_squared == 0.0);
  CHECK(result.performed_iterations == 0);
  CHECK(result.final_relative_change == 0.0);
  CHECK_FALSE(result.tolerance_met);
  CHECK(result.empty_matrix);
  CHECK(result.diagnostics.resource_creation_calls_before_iterations == 0);
  CHECK(result.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks_per_iteration == 0);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks_during_iterations ==
        0);
}

}  // namespace

TEST_CASE("hi-pdlp-operator-norm-is-exact-on-a-closed-form-fixture",
          "[pdlp][hipdlp][phase7e][operator-norm][exact]") {
  const hipdlp::algorithm::OperatorNormParameters defaults;
  CHECK(defaults.iteration_budget == 20);
  CHECK(defaults.relative_tolerance == 1e-6);

  const hipdlp::algorithm::OperatorNormResult result =
      runEstimate(hipdlp::backend::Kind::kCpu, makeExactLp());
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.error.code == hipdlp::api::ErrorCode::kNone);
  CHECK(result.operator_norm_squared == 9.0);
  CHECK(result.performed_iterations == 2);
  CHECK(result.final_relative_change == 0.0);
  CHECK(result.tolerance_met);
  CHECK_FALSE(result.empty_matrix);
  CHECK(result.diagnostics.resource_creation_calls_before_iterations == 3);
  CHECK(result.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks_per_iteration == 2);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks_during_iterations ==
        4);
}

TEST_CASE("hi-pdlp-operator-norm-upper-bound-covers-an-orthogonal-start",
          "[pdlp][hipdlp][phase7e][operator-norm][upper-bound]") {
  hipdlp::model::ProblemBuildResult built =
      hipdlp::model::buildProblem(makeOrthogonalDominantDirectionLp());
  REQUIRE(built.ok());
  double upper_bound_squared = 0.0;
  const hipdlp::api::Error bound_error =
      hipdlp::algorithm::computeOperatorNormSquaredUpperBound(
          built.problem.matrix(), upper_bound_squared);
  INFO(bound_error.message);
  REQUIRE(bound_error.ok());
  CHECK(upper_bound_squared >= 8.0);

  const hipdlp::algorithm::OperatorNormResult estimate = runEstimate(
      hipdlp::backend::Kind::kCpu, makeOrthogonalDominantDirectionLp());
  REQUIRE(estimate.ok());
  CHECK(estimate.operator_norm_squared == Approx(2.0));
  CHECK(upper_bound_squared > estimate.operator_norm_squared);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-operator-norm-cpu-cuda-agree-on-the-same-fixture",
          "[pdlp][hipdlp][phase7e][operator-norm][cuda][parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }

  const double cpu_cuda_relative_tolerance = 1e-13;
  const hipdlp::algorithm::OperatorNormResult cpu =
      runEstimate(hipdlp::backend::Kind::kCpu, makeExactLp());
  const hipdlp::algorithm::OperatorNormResult cuda =
      runEstimate(hipdlp::backend::Kind::kCuda, makeExactLp());
  INFO(cpu.error.message);
  INFO(cuda.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(cuda.ok());
  CHECK(cuda.operator_norm_squared ==
        Approx(cpu.operator_norm_squared).epsilon(cpu_cuda_relative_tolerance));
  CHECK(cuda.performed_iterations == cpu.performed_iterations);
  CHECK(cuda.tolerance_met == cpu.tolerance_met);
  CHECK(cuda.empty_matrix == cpu.empty_matrix);
}
#endif

TEST_CASE("hi-pdlp-operator-norm-repeated-evaluation-is-deterministic",
          "[pdlp][hipdlp][phase7e][operator-norm][determinism]") {
  hipdlp::algorithm::OperatorNormParameters parameters;
  parameters.iteration_budget = 6;
  parameters.relative_tolerance = std::numeric_limits<double>::denorm_min();
  const hipdlp::algorithm::OperatorNormResult first =
      runEstimate(hipdlp::backend::Kind::kCpu, makeConvergingLp(), parameters);
  INFO(first.error.message);
  REQUIRE(first.ok());

  for (int repetition = 0; repetition < 4; ++repetition) {
    const hipdlp::algorithm::OperatorNormResult repeated = runEstimate(
        hipdlp::backend::Kind::kCpu, makeConvergingLp(), parameters);
    INFO(repeated.error.message);
    REQUIRE(repeated.ok());
    CHECK(repeated.error.code == first.error.code);
    CHECK(repeated.error.message == first.error.message);
    CHECK(bitwiseEqual(repeated.operator_norm_squared,
                       first.operator_norm_squared));
    CHECK(repeated.performed_iterations == first.performed_iterations);
    CHECK(bitwiseEqual(repeated.final_relative_change,
                       first.final_relative_change));
    CHECK(repeated.tolerance_met == first.tolerance_met);
    CHECK(repeated.empty_matrix == first.empty_matrix);
  }
}

TEST_CASE("hi-pdlp-operator-norm-outcomes-are-transactional",
          "[pdlp][hipdlp][phase7e][operator-norm][errors]"
          "[transactional]") {
  hipdlp::model::ProblemBuildResult built =
      hipdlp::model::buildProblem(makeExactLp());
  REQUIRE(built.ok());
  hipdlp::backend::FactoryResult factory =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
  REQUIRE(factory.ok());
  hipdlp::backend::ProblemHandle problem;
  REQUIRE(factory.backend->bindProblem(built.problem, problem).ok());

  hipdlp::algorithm::OperatorNormParameters parameters;
  parameters.iteration_budget = 0;
  checkZeroEstimate(hipdlp::algorithm::estimateOperatorNormSquared(
                        *factory.backend, problem, 1, 1, parameters),
                    hipdlp::api::ErrorCode::kInvalidOption);
  parameters.iteration_budget = -1;
  checkZeroEstimate(hipdlp::algorithm::estimateOperatorNormSquared(
                        *factory.backend, problem, 1, 1, parameters),
                    hipdlp::api::ErrorCode::kInvalidOption);
  parameters.iteration_budget = kHighsIInf;
  checkZeroEstimate(hipdlp::algorithm::estimateOperatorNormSquared(
                        *factory.backend, problem, 1, 1, parameters),
                    hipdlp::api::ErrorCode::kInvalidOption);

  const double invalid_tolerances[] = {0.0, -1.0,
                                       std::numeric_limits<double>::quiet_NaN(),
                                       std::numeric_limits<double>::infinity()};
  for (const double tolerance : invalid_tolerances) {
    parameters = hipdlp::algorithm::OperatorNormParameters();
    parameters.relative_tolerance = tolerance;
    checkZeroEstimate(hipdlp::algorithm::estimateOperatorNormSquared(
                          *factory.backend, problem, 1, 1, parameters),
                      hipdlp::api::ErrorCode::kInvalidOption);
  }

  parameters = hipdlp::algorithm::OperatorNormParameters();
  checkZeroEstimate(hipdlp::algorithm::estimateOperatorNormSquared(
                        *factory.backend, problem, -1, 1, parameters),
                    hipdlp::api::ErrorCode::kInvalidArgument);
  checkZeroEstimate(hipdlp::algorithm::estimateOperatorNormSquared(
                        *factory.backend, problem, 1, -1, parameters),
                    hipdlp::api::ErrorCode::kInvalidArgument);
  hipdlp::backend::ProblemHandle invalid_problem;
  checkZeroEstimate(hipdlp::algorithm::estimateOperatorNormSquared(
                        *factory.backend, invalid_problem, 1, 1, parameters),
                    hipdlp::api::ErrorCode::kInvalidArgument);

  checkEmptyEstimate(runEstimate(hipdlp::backend::Kind::kCpu,
                                 makeEqualityLp(0, 1, {0, 0}, {}, {}),
                                 parameters));
  checkEmptyEstimate(runEstimate(hipdlp::backend::Kind::kCpu,
                                 makeEqualityLp(1, 0, {0}, {}, {}),
                                 parameters));

  parameters.iteration_budget = 1;
  parameters.relative_tolerance = std::numeric_limits<double>::denorm_min();
  const hipdlp::algorithm::OperatorNormResult exhausted =
      runEstimate(hipdlp::backend::Kind::kCpu, makeConvergingLp(), parameters);
  INFO(exhausted.error.message);
  REQUIRE(exhausted.ok());
  CHECK(exhausted.performed_iterations == 1);
  CHECK_FALSE(exhausted.tolerance_met);
  CHECK(exhausted.final_relative_change == 1.0);
  CHECK(exhausted.operator_norm_squared == Approx(65.0 / 17.0));

  const hipdlp::algorithm::OperatorNormResult zero_matrix =
      runEstimate(hipdlp::backend::Kind::kCpu,
                  makeEqualityLp(1, 1, {0, 0}, {}, {}), parameters);
  checkZeroEstimate(zero_matrix, hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(zero_matrix.diagnostics.resource_creation_calls_before_iterations == 3);
  CHECK(zero_matrix.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(zero_matrix.diagnostics.device_to_host_scalar_readbacks_per_iteration ==
        2);
  CHECK(zero_matrix.diagnostics
            .device_to_host_scalar_readbacks_during_iterations == 1);
}

TEST_CASE("hi-pdlp-operator-norm-allocates-nothing-inside-iteration",
          "[pdlp][hipdlp][phase7e][operator-norm][resources][readbacks]") {
  hipdlp::algorithm::OperatorNormParameters parameters;
  parameters.iteration_budget = 3;
  parameters.relative_tolerance = std::numeric_limits<double>::denorm_min();
  const hipdlp::algorithm::OperatorNormResult result =
      runEstimate(hipdlp::backend::Kind::kCpu, makeConvergingLp(), parameters);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.performed_iterations == parameters.iteration_budget);
  CHECK_FALSE(result.tolerance_met);
  CHECK(result.diagnostics.resource_creation_calls_before_iterations == 3);
  CHECK(result.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks_per_iteration == 2);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks_during_iterations ==
        2 * static_cast<std::size_t>(result.performed_iterations));
}