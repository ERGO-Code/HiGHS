#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "HConfig.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "pdlp/hipdlp/algorithm/fixed_point_error.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/backend/cpu/cpu_backend.hpp"
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
#include "pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp"
#endif
#include "pdlp/hipdlp/model/problem.hpp"

namespace {

HighsLp makeRowwiseEqualityLp(const HighsInt rows, const HighsInt cols,
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
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.num_row_ = rows;
  lp.a_matrix_.num_col_ = cols;
  lp.a_matrix_.start_ = starts;
  lp.a_matrix_.index_ = indices;
  lp.a_matrix_.value_ = values;
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0.0;
  return lp;
}

HighsLp makeScalarLp(const double coefficient) {
  return makeRowwiseEqualityLp(1, 1, {0, 1}, {0}, {coefficient});
}

struct EvaluationFixture {
  HighsLp lp;
  std::unique_ptr<hipdlp::backend::Backend> backend;
  hipdlp::backend::ProblemHandle problem;
  hipdlp::backend::VectorHandle next_primal;
  hipdlp::backend::VectorHandle next_dual;
  hipdlp::backend::VectorHandle reflected_primal;
  hipdlp::backend::VectorHandle reflected_dual;
  hipdlp::backend::VectorHandle delta_primal;
  hipdlp::backend::VectorHandle delta_dual;
  hipdlp::backend::VectorHandle aty_delta;
};

std::unique_ptr<EvaluationFixture> makeFixtureWithBackend(
    std::unique_ptr<hipdlp::backend::Backend> backend, const HighsLp& lp,
    const std::vector<double>& next_primal,
    const std::vector<double>& next_dual,
    const std::vector<double>& reflected_primal,
    const std::vector<double>& reflected_dual) {
  hipdlp::model::ProblemBuildResult built = hipdlp::model::buildProblem(lp);
  INFO(built.error.message);
  REQUIRE(built.ok());
  REQUIRE(backend.get() != nullptr);

  std::unique_ptr<EvaluationFixture> fixture(new EvaluationFixture());
  fixture->lp = lp;
  fixture->backend = std::move(backend);
  REQUIRE(fixture->backend->bindProblem(built.problem, fixture->problem).ok());
  REQUIRE(
      fixture->backend->createVector(lp.num_col_, fixture->next_primal).ok());
  REQUIRE(fixture->backend->createVector(lp.num_row_, fixture->next_dual).ok());
  REQUIRE(fixture->backend->createVector(lp.num_col_, fixture->reflected_primal)
              .ok());
  REQUIRE(fixture->backend->createVector(lp.num_row_, fixture->reflected_dual)
              .ok());
  REQUIRE(
      fixture->backend->createVector(lp.num_col_, fixture->delta_primal).ok());
  REQUIRE(
      fixture->backend->createVector(lp.num_row_, fixture->delta_dual).ok());
  REQUIRE(fixture->backend->createVector(lp.num_col_, fixture->aty_delta).ok());
  REQUIRE(
      fixture->backend->uploadVector(next_primal, fixture->next_primal).ok());
  REQUIRE(fixture->backend->uploadVector(next_dual, fixture->next_dual).ok());
  REQUIRE(fixture->backend
              ->uploadVector(reflected_primal, fixture->reflected_primal)
              .ok());
  REQUIRE(
      fixture->backend->uploadVector(reflected_dual, fixture->reflected_dual)
          .ok());
  REQUIRE(fixture->backend->fill(fixture->delta_primal, 71.0).ok());
  REQUIRE(fixture->backend->fill(fixture->delta_dual, 73.0).ok());
  REQUIRE(fixture->backend->fill(fixture->aty_delta, 79.0).ok());
  return fixture;
}

std::unique_ptr<EvaluationFixture> makeFixture(
    const hipdlp::backend::Kind kind, const HighsLp& lp,
    const std::vector<double>& next_primal,
    const std::vector<double>& next_dual,
    const std::vector<double>& reflected_primal,
    const std::vector<double>& reflected_dual) {
  hipdlp::backend::FactoryResult factory = hipdlp::backend::createBackend(kind);
  INFO(factory.error.message);
  REQUIRE(factory.ok());
  return makeFixtureWithBackend(std::move(factory.backend), lp, next_primal,
                                next_dual, reflected_primal, reflected_dual);
}

std::vector<double> download(hipdlp::backend::Backend& backend,
                             const hipdlp::backend::VectorHandle& vector) {
  std::vector<double> values;
  REQUIRE(backend.downloadVector(vector, values).ok());
  return values;
}

std::array<hipdlp::backend::VectorHandle*, 7> fixtureHandles(
    EvaluationFixture& fixture) {
  return {{&fixture.next_primal, &fixture.next_dual, &fixture.reflected_primal,
           &fixture.reflected_dual, &fixture.delta_primal, &fixture.delta_dual,
           &fixture.aty_delta}};
}

hipdlp::algorithm::FixedPointErrorInput makeInput(
    const std::array<hipdlp::backend::VectorHandle*, 7>& handles,
    const double omega, const double eta) {
  hipdlp::algorithm::FixedPointErrorInput input;
  input.next_primal = handles[0];
  input.next_dual = handles[1];
  input.reflected_primal = handles[2];
  input.reflected_dual = handles[3];
  input.omega = omega;
  input.eta = eta;
  input.delta_primal = handles[4];
  input.delta_dual = handles[5];
  input.aty_delta = handles[6];
  return input;
}

hipdlp::algorithm::FixedPointErrorResult evaluate(EvaluationFixture& fixture,
                                                  const double omega,
                                                  const double eta) {
  return hipdlp::algorithm::evaluateFixedPointError(
      *fixture.backend, fixture.problem,
      makeInput(fixtureHandles(fixture), omega, eta));
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-cuda-fixed-point-reuses-prepared-products",
      "[hipdlp][cuda][fixed-point-error][check-optimization]") {
  auto fixture = makeFixture(hipdlp::backend::Kind::kCuda, makeScalarLp(1.0),
               {1.0}, {2.0}, {0.0}, {0.0});
  const auto input = makeInput(fixtureHandles(*fixture), 1.0, 1.0);
  hipdlp::backend::cuda::testing::resetResourceAudit();
  REQUIRE(hipdlp::algorithm::prepareFixedPointError(
    *fixture->backend, fixture->problem, input).ok());
  const auto preparation = hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(preparation.device_to_device_transfers == 0);
  CHECK(preparation.direct_primitive_enqueue_calls == 4);
  hipdlp::backend::cuda::testing::resetResourceAudit();
  const auto result = hipdlp::algorithm::evaluateFixedPointError(
    *fixture->backend, fixture->problem, input, true);
  REQUIRE(result.ok());
  CHECK(result.fixed_point_error == 3.0);
  const auto audit = hipdlp::backend::cuda::testing::resourceAudit();
  CHECK(audit.direct_primitive_enqueue_calls == 4);
  CHECK(audit.device_to_device_transfers == 0);
  CHECK(audit.device_to_host_transfers == 2);
  CHECK(audit.synchronization_calls == 1);
  CHECK(audit.allocation_calls == 0);
}
#endif

TEST_CASE("hi-pdlp-fixed-point-validity-is-independent-of-readback-count",
          "[hipdlp][phase8][fixed-point-error]") {
  hipdlp::algorithm::FixedPointErrorResult uncomputed;
  uncomputed.diagnostics.device_to_host_scalar_readbacks = 3;
  CHECK_FALSE(uncomputed.ok());

  auto fixture = makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(1.0),
                             {1.0}, {2.0}, {0.0}, {0.0});
  auto result = evaluate(*fixture, 1.0, 1.0);
  REQUIRE(result.ok());
  CHECK(result.fixed_point_error == 3.0);
  for (const std::size_t readbacks : {0u, 1u, 7u}) {
    result.diagnostics.device_to_host_scalar_readbacks = readbacks;
    CHECK(result.ok());
  }
  result.error = hipdlp::api::Error(hipdlp::api::ErrorCode::kBackendFailure,
                                   "failed result must remain invalid");
  CHECK_FALSE(result.ok());
}

double hostFixedPointError(const HighsLp& lp,
                           const std::vector<double>& next_primal,
                           const std::vector<double>& next_dual,
                           const std::vector<double>& reflected_primal,
                           const std::vector<double>& reflected_dual,
                           const double omega, const double eta) {
  REQUIRE(lp.a_matrix_.format_ == MatrixFormat::kRowwise);
  REQUIRE(next_primal.size() == reflected_primal.size());
  REQUIRE(next_dual.size() == reflected_dual.size());

  std::vector<double> delta_primal(next_primal.size(), 0.0);
  std::vector<double> delta_dual(next_dual.size(), 0.0);
  double primal_norm_squared = 0.0;
  double dual_norm_squared = 0.0;
  for (std::size_t col = 0; col < next_primal.size(); ++col) {
    delta_primal[col] = next_primal[col] - reflected_primal[col];
    primal_norm_squared += delta_primal[col] * delta_primal[col];
  }
  for (std::size_t row = 0; row < next_dual.size(); ++row) {
    delta_dual[row] = next_dual[row] - reflected_dual[row];
    dual_norm_squared += delta_dual[row] * delta_dual[row];
  }

  std::vector<double> aty_delta(next_primal.size(), 0.0);
  for (HighsInt row = 0; row < lp.num_row_; ++row) {
    for (HighsInt entry = lp.a_matrix_.start_[row];
         entry < lp.a_matrix_.start_[row + 1]; ++entry) {
      aty_delta[lp.a_matrix_.index_[entry]] +=
          lp.a_matrix_.value_[entry] * delta_dual[row];
    }
  }
  double cross_term = 0.0;
  for (std::size_t col = 0; col < delta_primal.size(); ++col)
    cross_term += delta_primal[col] * aty_delta[col];

  const double movement =
      primal_norm_squared * omega + dual_norm_squared / omega;
  const double interaction = 2.0 * eta * cross_term;
  return std::sqrt(std::max(0.0, movement + interaction));
}

bool bitwiseEqual(const double left, const double right) {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

void checkZeroFailure(const hipdlp::algorithm::FixedPointErrorResult& result,
                      const hipdlp::api::ErrorCode expected_code) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.error.code == expected_code);
  CHECK(result.fixed_point_error == 0.0);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks == 0);
}

class ScriptedReductionBackend : public hipdlp::backend::cpu::CpuBackend {
 public:
  ScriptedReductionBackend()
      : norm_values{{1.0, 1.0}},
        dot_value(0.0),
        norm_calls(0),
        dot_calls(0),
        fail_first_norm(false) {}

  hipdlp::api::Error norm2(const hipdlp::backend::VectorHandle& vector,
                           double& result) override {
    (void)vector;
    ++norm_calls;
    if (fail_first_norm && norm_calls == 1) {
      return hipdlp::api::Error(hipdlp::api::ErrorCode::kBackendFailure,
                                "injected first reduction failure");
    }
    result = norm_values[norm_calls - 1];
    return hipdlp::api::Error();
  }

  hipdlp::api::Error dot(const hipdlp::backend::VectorHandle& left,
                         const hipdlp::backend::VectorHandle& right,
                         double& result) override {
    (void)left;
    (void)right;
    ++dot_calls;
    result = dot_value;
    return hipdlp::api::Error();
  }

  std::array<double, 2> norm_values;
  double dot_value;
  std::size_t norm_calls;
  std::size_t dot_calls;
  bool fail_first_norm;
};

}  // namespace

TEST_CASE("hi-pdlp-fixed-point-scalar-finalizer-without-readbacks",
      "[pdlp][hipdlp][phase8][fixed-point-error]") {
  const hipdlp::algorithm::FixedPointErrorResult result =
    hipdlp::algorithm::computeFixedPointError(1.0, 2.0, 2.0, 1.0, 1.0);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.computed);
  CHECK(result.fixed_point_error == 3.0);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks == 0);
  const hipdlp::algorithm::FixedPointErrorResult clamped =
    hipdlp::algorithm::computeFixedPointError(2.0, 1.0, -6.0, 1.0, 1.0);
  REQUIRE(clamped.ok());
  CHECK(clamped.fixed_point_error == 0.0);
  const hipdlp::algorithm::FixedPointErrorResult zero =
    hipdlp::algorithm::computeFixedPointError(0.0, 0.0, 0.0, 7.0, 0.25);
  REQUIRE(zero.ok());
  CHECK(zero.fixed_point_error == 0.0);
}

TEST_CASE("hi-pdlp-fixed-point-scalar-finalizer-failures-are-uncomputed",
      "[pdlp][hipdlp][phase8][fixed-point-error][errors]") {
  const double infinity = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double maximum = std::numeric_limits<double>::max();
  const std::array<std::array<double, 5>, 13> invalid = {{
    {{-1.0, 1.0, 0.0, 1.0, 1.0}},
    {{nan, 1.0, 0.0, 1.0, 1.0}},
    {{1.0, -1.0, 0.0, 1.0, 1.0}},
    {{1.0, infinity, 0.0, 1.0, 1.0}},
    {{1.0, 1.0, infinity, 1.0, 1.0}},
    {{1.0, 1.0, 0.0, 0.0, 1.0}},
    {{1.0, 1.0, 0.0, 1.0, nan}},
    {{1e200, 0.0, 0.0, 1.0, 1.0}},
    {{1e154, 0.0, 0.0, 2.0, 1.0}},
    {{1e154, 1e154, 0.0, 1.0, 1.0}},
    {{0.0, 0.0, 0.0, 1.0, maximum}},
    {{1.0, 1.0, 1e308, 1.0, 1.0}},
    {{1.0, 1.0, 5e307, 1e308, 1.0}},
  }};
  for (const std::array<double, 5>& values : invalid) {
  const hipdlp::algorithm::FixedPointErrorResult result =
    hipdlp::algorithm::computeFixedPointError(
      values[0], values[1], values[2], values[3], values[4]);
  checkZeroFailure(result, hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK_FALSE(result.computed);
  }
}

TEST_CASE("hi-pdlp-fixed-point-error-matches-the-legacy-formula",
          "[pdlp][hipdlp][phase7h][fixed-point-error][exact]") {
  std::unique_ptr<EvaluationFixture> fixture =
      makeFixture(hipdlp::backend::Kind::kCpu,
                  makeRowwiseEqualityLp(2, 3, {0, 3, 6}, {0, 1, 2, 0, 1, 2},
                                        {-1.0, 1.0, 2.0, -1.0, -1.0, 3.0}),
                  {4.0, 6.0, -1.0}, {5.0, 3.0}, {1.0, 2.0, -1.0}, {2.0, -1.0});
  hipdlp::backend::Backend& backend = *fixture->backend;
  const std::vector<double> next_primal =
      download(backend, fixture->next_primal);
  const std::vector<double> next_dual = download(backend, fixture->next_dual);
  const std::vector<double> reflected_primal =
      download(backend, fixture->reflected_primal);
  const std::vector<double> reflected_dual =
      download(backend, fixture->reflected_dual);
  const double expected =
      hostFixedPointError(fixture->lp, next_primal, next_dual, reflected_primal,
                          reflected_dual, 1.0, 0.5);

  const hipdlp::algorithm::FixedPointErrorResult result =
      evaluate(*fixture, 1.0, 0.5);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.error.code == hipdlp::api::ErrorCode::kNone);
  CHECK(expected == 5.0);
  CHECK(result.fixed_point_error == expected);
  CHECK(result.diagnostics.device_to_host_scalar_readbacks == 3);

  CHECK(download(backend, fixture->next_primal) == next_primal);
  CHECK(download(backend, fixture->next_dual) == next_dual);
  CHECK(download(backend, fixture->reflected_primal) == reflected_primal);
  CHECK(download(backend, fixture->reflected_dual) == reflected_dual);
  CHECK(download(backend, fixture->delta_primal) ==
        std::vector<double>{3.0, 4.0, 0.0});
  CHECK(download(backend, fixture->delta_dual) ==
        std::vector<double>{3.0, 4.0});
  CHECK(download(backend, fixture->aty_delta) ==
        std::vector<double>{-7.0, -1.0, 18.0});
}

TEST_CASE("hi-pdlp-fixed-point-error-clamps-negative-and-zero-displacement",
          "[pdlp][hipdlp][phase7h][fixed-point-error][clamp]") {
  SECTION("negative interaction clamps to zero") {
    std::unique_ptr<EvaluationFixture> fixture =
        makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(3.0), {3.0},
                    {0.0}, {1.0}, {1.0});
    const hipdlp::algorithm::FixedPointErrorResult result =
        evaluate(*fixture, 1.0, 1.0);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.fixed_point_error == 0.0);
    CHECK(result.diagnostics.device_to_host_scalar_readbacks == 3);
  }

  SECTION("zero displacement is exactly zero") {
    std::unique_ptr<EvaluationFixture> fixture =
        makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(-4.0), {2.0},
                    {-3.0}, {2.0}, {-3.0});
    const hipdlp::algorithm::FixedPointErrorResult result =
        evaluate(*fixture, 7.0, 0.25);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.fixed_point_error == 0.0);
    CHECK(result.diagnostics.device_to_host_scalar_readbacks == 3);
  }
}

TEST_CASE("hi-pdlp-fixed-point-error-consumes-current-omega-and-stored-eta",
          "[pdlp][hipdlp][phase7h][fixed-point-error][parameters]") {
  std::unique_ptr<EvaluationFixture> fixture =
      makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(3.0), {3.0}, {2.0},
                  {1.0}, {1.0});
  hipdlp::backend::Backend& backend = *fixture->backend;
  const std::vector<double> next_primal =
      download(backend, fixture->next_primal);
  const std::vector<double> next_dual = download(backend, fixture->next_dual);
  const std::vector<double> reflected_primal =
      download(backend, fixture->reflected_primal);
  const std::vector<double> reflected_dual =
      download(backend, fixture->reflected_dual);
  const std::array<std::array<double, 2>, 4> parameters = {
      {{1.0, 0.25}, {2.0, 0.5}, {0.5, 0.125}, {4.0, 1.375}}};

  for (const std::array<double, 2>& parameter : parameters) {
    const double expected = hostFixedPointError(
        fixture->lp, next_primal, next_dual, reflected_primal, reflected_dual,
        parameter[0], parameter[1]);
    const hipdlp::algorithm::FixedPointErrorResult result =
        evaluate(*fixture, parameter[0], parameter[1]);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.fixed_point_error == expected);
    CHECK(result.diagnostics.device_to_host_scalar_readbacks == 3);
  }
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-fixed-point-error-cpu-cuda-agree-without-new-resources",
          "[pdlp][hipdlp][phase7h][fixed-point-error][cuda][resources]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  const HighsLp lp = makeRowwiseEqualityLp(2, 3, {0, 3, 6}, {0, 1, 2, 0, 1, 2},
                                           {-1.0, 1.0, 2.0, -1.0, -1.0, 3.0});
  const std::vector<double> next_primal = {4.0, 6.0, -1.0};
  const std::vector<double> next_dual = {5.0, 3.0};
  const std::vector<double> reflected_primal = {1.0, 2.0, -1.0};
  const std::vector<double> reflected_dual = {2.0, -1.0};
  const double expected = hostFixedPointError(
      lp, next_primal, next_dual, reflected_primal, reflected_dual, 1.0, 0.5);
  REQUIRE(expected == 5.0);
  std::unique_ptr<EvaluationFixture> cpu =
      makeFixture(hipdlp::backend::Kind::kCpu, lp, next_primal, next_dual,
                  reflected_primal, reflected_dual);
  std::unique_ptr<EvaluationFixture> cuda =
      makeFixture(hipdlp::backend::Kind::kCuda, lp, next_primal, next_dual,
                  reflected_primal, reflected_dual);

  const hipdlp::algorithm::FixedPointErrorResult cpu_result =
      evaluate(*cpu, 1.0, 0.5);
  hipdlp::backend::cuda::testing::resetResourceAudit();
  const hipdlp::algorithm::FixedPointErrorResult cuda_result =
      evaluate(*cuda, 1.0, 0.5);
  const hipdlp::algorithm::FixedPointErrorResult repeated_cuda_result =
      evaluate(*cuda, 1.0, 0.5);
  const hipdlp::backend::cuda::testing::ResourceAudit audit =
      hipdlp::backend::cuda::testing::resourceAudit();
  INFO(cpu_result.error.message);
  INFO(cuda_result.error.message);
  INFO(repeated_cuda_result.error.message);
  REQUIRE(cpu_result.ok());
  REQUIRE(cuda_result.ok());
  REQUIRE(repeated_cuda_result.ok());
  CHECK(cpu_result.fixed_point_error == expected);
  CHECK(cuda_result.fixed_point_error == Approx(expected).epsilon(5e-13));
  CHECK(cuda_result.fixed_point_error ==
        Approx(cpu_result.fixed_point_error).epsilon(5e-13));
  CHECK(bitwiseEqual(repeated_cuda_result.fixed_point_error,
                     cuda_result.fixed_point_error));
  CHECK(cpu_result.diagnostics.device_to_host_scalar_readbacks == 3);
  CHECK(cuda_result.diagnostics.device_to_host_scalar_readbacks == 3);
  CHECK(repeated_cuda_result.diagnostics.device_to_host_scalar_readbacks == 3);
  CHECK(audit.allocation_calls == 0);
  CHECK(audit.free_calls == 0);
  CHECK(audit.descriptor_creations == 0);
  CHECK(audit.descriptor_destructions == 0);
  CHECK(audit.host_to_device_transfers == 0);
  CHECK(audit.device_to_host_transfers == 4);
  CHECK(audit.device_to_host_bytes == 2 * (3 * sizeof(double) + sizeof(int)));
  CHECK(audit.synchronization_calls == 2);
  CHECK(download(*cuda->backend, cuda->next_primal) == next_primal);
  CHECK(download(*cuda->backend, cuda->next_dual) == next_dual);
  CHECK(download(*cuda->backend, cuda->reflected_primal) == reflected_primal);
  CHECK(download(*cuda->backend, cuda->reflected_dual) == reflected_dual);
  CHECK(download(*cuda->backend, cuda->delta_primal) ==
        std::vector<double>{3.0, 4.0, 0.0});
  CHECK(download(*cuda->backend, cuda->delta_dual) ==
        std::vector<double>{3.0, 4.0});
  CHECK(download(*cuda->backend, cuda->aty_delta) ==
        std::vector<double>{-7.0, -1.0, 18.0});
}
#endif

TEST_CASE("hi-pdlp-fixed-point-error-validates-handles-scalars-and-dimensions",
          "[pdlp][hipdlp][phase7h][fixed-point-error][errors]"
          "[transactional]") {
  std::unique_ptr<EvaluationFixture> fixture =
      makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(3.0), {3.0}, {2.0},
                  {1.0}, {1.0});
  const std::array<hipdlp::backend::VectorHandle*, 7> valid_handles =
      fixtureHandles(*fixture);

  hipdlp::backend::ProblemHandle invalid_problem;
  checkZeroFailure(hipdlp::algorithm::evaluateFixedPointError(
                       *fixture->backend, invalid_problem,
                       makeInput(valid_handles, 1.0, 1.0)),
                   hipdlp::api::ErrorCode::kInvalidArgument);

  for (std::size_t missing = 0; missing < valid_handles.size(); ++missing) {
    std::array<hipdlp::backend::VectorHandle*, 7> handles = valid_handles;
    handles[missing] = nullptr;
    INFO(missing);
    checkZeroFailure(
        hipdlp::algorithm::evaluateFixedPointError(
            *fixture->backend, fixture->problem, makeInput(handles, 1.0, 1.0)),
        hipdlp::api::ErrorCode::kInvalidArgument);
  }

  hipdlp::backend::VectorHandle invalid_vector;
  for (std::size_t invalid = 0; invalid < valid_handles.size(); ++invalid) {
    std::array<hipdlp::backend::VectorHandle*, 7> handles = valid_handles;
    handles[invalid] = &invalid_vector;
    INFO(invalid);
    checkZeroFailure(
        hipdlp::algorithm::evaluateFixedPointError(
            *fixture->backend, fixture->problem, makeInput(handles, 1.0, 1.0)),
        hipdlp::api::ErrorCode::kInvalidArgument);
  }

  for (std::size_t left = 0; left < valid_handles.size(); ++left) {
    for (std::size_t right = left + 1; right < valid_handles.size(); ++right) {
      std::array<hipdlp::backend::VectorHandle*, 7> handles = valid_handles;
      handles[right] = handles[left];
      INFO(left);
      INFO(right);
      checkZeroFailure(hipdlp::algorithm::evaluateFixedPointError(
                           *fixture->backend, fixture->problem,
                           makeInput(handles, 1.0, 1.0)),
                       hipdlp::api::ErrorCode::kInvalidArgument);
    }
  }

  hipdlp::backend::VectorHandle wrong_size;
  REQUIRE(fixture->backend->createVector(2, wrong_size).ok());
  REQUIRE(fixture->backend->fill(wrong_size, 0.0).ok());
  for (std::size_t wrong = 0; wrong < valid_handles.size(); ++wrong) {
    std::array<hipdlp::backend::VectorHandle*, 7> handles = valid_handles;
    handles[wrong] = &wrong_size;
    INFO(wrong);
    checkZeroFailure(
        hipdlp::algorithm::evaluateFixedPointError(
            *fixture->backend, fixture->problem, makeInput(handles, 1.0, 1.0)),
        hipdlp::api::ErrorCode::kInvalidArgument);
  }

  const double invalid_scalars[] = {0.0, -1.0,
                                    std::numeric_limits<double>::quiet_NaN(),
                                    std::numeric_limits<double>::infinity()};
  for (const double invalid : invalid_scalars) {
    INFO(invalid);
    checkZeroFailure(hipdlp::algorithm::evaluateFixedPointError(
                         *fixture->backend, fixture->problem,
                         makeInput(valid_handles, invalid, 1.0)),
                     hipdlp::api::ErrorCode::kInvalidArgument);
    checkZeroFailure(hipdlp::algorithm::evaluateFixedPointError(
                         *fixture->backend, fixture->problem,
                         makeInput(valid_handles, 1.0, invalid)),
                     hipdlp::api::ErrorCode::kInvalidArgument);
  }
}

TEST_CASE("hi-pdlp-fixed-point-error-rejects-invalid-reductions-and-overflow",
          "[pdlp][hipdlp][phase7h][fixed-point-error][errors][overflow]") {
  SECTION("reduction scalars are independently validated") {
    const std::array<std::array<double, 3>, 4> reductions = {{
        {{-1.0, 1.0, 0.0}},
        {{std::numeric_limits<double>::quiet_NaN(), 1.0, 0.0}},
        {{1.0, -1.0, 0.0}},
        {{1.0, 1.0, std::numeric_limits<double>::infinity()}},
    }};
    for (const std::array<double, 3>& reduction : reductions) {
      std::unique_ptr<ScriptedReductionBackend> scripted(
          new ScriptedReductionBackend());
      scripted->norm_values = {{reduction[0], reduction[1]}};
      scripted->dot_value = reduction[2];
      std::unique_ptr<EvaluationFixture> fixture = makeFixtureWithBackend(
          std::unique_ptr<hipdlp::backend::Backend>(scripted.release()),
          makeScalarLp(1.0), {1.0}, {1.0}, {0.0}, {0.0});
      checkZeroFailure(evaluate(*fixture, 1.0, 1.0),
                       hipdlp::api::ErrorCode::kInvalidArgument);
    }
  }

  SECTION("the first backend error is preserved") {
    std::unique_ptr<ScriptedReductionBackend> scripted(
        new ScriptedReductionBackend());
    ScriptedReductionBackend* scripted_pointer = scripted.get();
    scripted->fail_first_norm = true;
    std::unique_ptr<EvaluationFixture> fixture = makeFixtureWithBackend(
        std::unique_ptr<hipdlp::backend::Backend>(scripted.release()),
        makeScalarLp(1.0), {1.0}, {1.0}, {0.0}, {0.0});
    const hipdlp::algorithm::FixedPointErrorResult result =
        evaluate(*fixture, 1.0, 1.0);
    checkZeroFailure(result, hipdlp::api::ErrorCode::kBackendFailure);
    CHECK(result.error.message == "injected first reduction failure");
    CHECK(scripted_pointer->norm_calls == 1);
    CHECK(scripted_pointer->dot_calls == 0);
  }

  SECTION("norm square overflow is transactional") {
    std::unique_ptr<EvaluationFixture> fixture =
        makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(1.0), {1e200},
                    {0.0}, {0.0}, {0.0});
    checkZeroFailure(evaluate(*fixture, 1.0, 1.0),
                     hipdlp::api::ErrorCode::kInvalidArgument);
  }

  SECTION("weighted movement overflow is transactional") {
    std::unique_ptr<EvaluationFixture> fixture =
        makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(1.0), {1e154},
                    {0.0}, {0.0}, {0.0});
    checkZeroFailure(evaluate(*fixture, 2.0, 1.0),
                     hipdlp::api::ErrorCode::kInvalidArgument);
  }

  SECTION("interaction and final sum overflow are transactional") {
    std::unique_ptr<EvaluationFixture> interaction =
        makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(1e308), {1.0},
                    {1.0}, {0.0}, {0.0});
    checkZeroFailure(evaluate(*interaction, 1.0, 1.0),
                     hipdlp::api::ErrorCode::kInvalidArgument);

    std::unique_ptr<EvaluationFixture> sum =
        makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(5e307), {1.0},
                    {1.0}, {0.0}, {0.0});
    checkZeroFailure(evaluate(*sum, 1e308, 1.0),
                     hipdlp::api::ErrorCode::kInvalidArgument);
  }

  SECTION("each multiplicative factor is validated") {
    std::unique_ptr<EvaluationFixture> fixture =
        makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(0.0), {0.0},
                    {0.0}, {0.0}, {0.0});
    checkZeroFailure(
        evaluate(*fixture, 1.0, std::numeric_limits<double>::max()),
        hipdlp::api::ErrorCode::kInvalidArgument);
  }
}

TEST_CASE("hi-pdlp-fixed-point-error-is-deterministic-and-mutation-is-local",
          "[pdlp][hipdlp][phase7h][fixed-point-error][determinism]") {
  std::unique_ptr<EvaluationFixture> fixture =
      makeFixture(hipdlp::backend::Kind::kCpu, makeScalarLp(3.0), {3.0}, {2.0},
                  {1.0}, {1.0});
  hipdlp::backend::Backend& backend = *fixture->backend;
  const std::vector<double> next_primal =
      download(backend, fixture->next_primal);
  const std::vector<double> next_dual = download(backend, fixture->next_dual);
  const std::vector<double> reflected_primal =
      download(backend, fixture->reflected_primal);
  const std::vector<double> reflected_dual =
      download(backend, fixture->reflected_dual);
  const hipdlp::algorithm::FixedPointErrorResult first =
      evaluate(*fixture, 2.0, 1.375);
  INFO(first.error.message);
  REQUIRE(first.ok());

  for (int repetition = 0; repetition < 4; ++repetition) {
    const hipdlp::algorithm::FixedPointErrorResult repeated =
        evaluate(*fixture, 2.0, 1.375);
    INFO(repeated.error.message);
    REQUIRE(repeated.ok());
    CHECK(repeated.error.code == first.error.code);
    CHECK(repeated.error.message == first.error.message);
    CHECK(bitwiseEqual(repeated.fixed_point_error, first.fixed_point_error));
    CHECK(repeated.diagnostics.device_to_host_scalar_readbacks ==
          first.diagnostics.device_to_host_scalar_readbacks);
  }

  CHECK(download(backend, fixture->next_primal) == next_primal);
  CHECK(download(backend, fixture->next_dual) == next_dual);
  CHECK(download(backend, fixture->reflected_primal) == reflected_primal);
  CHECK(download(backend, fixture->reflected_dual) == reflected_dual);
  CHECK(download(backend, fixture->delta_primal) == std::vector<double>{2.0});
  CHECK(download(backend, fixture->delta_dual) == std::vector<double>{1.0});
  CHECK(download(backend, fixture->aty_delta) == std::vector<double>{3.0});
}