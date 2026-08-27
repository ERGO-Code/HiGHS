#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#include "HConfig.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "lp_data/HighsOptions.h"
#include "pdlp/hipdlp/algorithm/fixed_trace.hpp"
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
#include "pdlp/hipdlp/backend/cuda/cuda_backend_testing.hpp"
#endif
#ifndef CUPDLP_GPU
#include "pdlp/hipdlp/pdhg.hpp"
#endif
#include "util/HighsTimer.h"

namespace {

HighsLp makeOneStepLp() {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {-1.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {2.0};
  lp.row_lower_ = {2.0};
  lp.row_upper_ = {2.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1.0};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0.0;
  return lp;
}

HighsLp makeProjectionLp() {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 2;
  lp.col_cost_ = {0.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {0.0};
  lp.row_lower_ = {1.0, 1.0};
  lp.row_upper_ = {1.0, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 2;
  lp.a_matrix_.start_ = {0, 0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

HighsLp makeColumnBoundsLp() {
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 0;
  lp.col_cost_ = {-2.0, 1.0, -1.0};
  lp.col_lower_ = {1.0, -kHighsInf, 0.0};
  lp.col_upper_ = {1.0, kHighsInf, 0.5};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 0;
  lp.a_matrix_.start_ = {0, 0, 0, 0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

HighsLp makeFreeColumnLp() {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 0;
  lp.col_cost_ = {-1.0};
  lp.col_lower_ = {-kHighsInf};
  lp.col_upper_ = {kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 0;
  lp.a_matrix_.start_ = {0, 0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

HighsLp makeScalingLp() {
  HighsLp lp = makeOneStepLp();
  lp.col_cost_ = {-2.0};
  lp.col_upper_ = {10.0};
  lp.row_lower_ = {8.0};
  lp.row_upper_ = {8.0};
  lp.a_matrix_.value_ = {4.0};
  return lp;
}

HighsLp makeRectangularLp(const bool skewed) {
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {-1.0, 2.0, -0.5};
  lp.col_lower_ = {-kHighsInf, 0.0, 1.0};
  lp.col_upper_ = {kHighsInf, 3.0, 1.0};
  lp.row_lower_ = {1.0, -2.0};
  lp.row_upper_ = {1.0, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 2;
  lp.a_matrix_.start_ = {0, 2, 3, 5};
  lp.a_matrix_.index_ = {0, 1, 0, 0, 1};
  lp.a_matrix_.value_ = skewed
                            ? std::vector<double>({1e-9, -2.0, 3.0, -1.0, 1e9})
                            : std::vector<double>({1.0, -2.0, 3.0, -1.0, 4.0});
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

HighsLp makeEmptyLp() {
  HighsLp lp;
  lp.num_col_ = 0;
  lp.num_row_ = 0;
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 0;
  lp.a_matrix_.num_row_ = 0;
  lp.a_matrix_.start_ = {0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

hipdlp::algorithm::FixedTraceOptions makeOptions(
    const HighsInt iteration_count) {
  hipdlp::algorithm::FixedTraceOptions options;
  options.backend = hipdlp::backend::Kind::kCpu;
  options.scaling.mode = hipdlp::model::ScalingMode::kOff;
  options.primal_step_size = 0.5;
  options.dual_step_size = 0.25;
  options.iteration_count = iteration_count;
  options.scaled_primal_start = {0.0};
  options.scaled_dual_start = {0.0};
  return options;
}

void checkExactPoint(const hipdlp::algorithm::FixedTracePoint& point,
                     const std::vector<double>& primal,
                     const std::vector<double>& dual,
                     const std::vector<double>& ax,
                     const std::vector<double>& aty) {
  CHECK(point.primal == primal);
  CHECK(point.dual == dual);
  CHECK(point.ax == ax);
  CHECK(point.aty == aty);
}

void checkNoAverage(const hipdlp::algorithm::FixedTraceResult& result) {
  CHECK_FALSE(result.averaging.hasAverage());
  CHECK(result.averaged_primal.empty());
  CHECK(result.averaged_dual.empty());
}

void checkNoObservation(
    const hipdlp::algorithm::FixedTraceConvergenceObservation& observation) {
  CHECK_FALSE(observation.observed);
  CHECK_FALSE(observation.stopping.metrics_valid);
  CHECK(observation.stopping.decision ==
        hipdlp::algorithm::StoppingDecision::kNumericalError);
  CHECK(observation.raw.primal_residual == 0.0);
  CHECK(observation.raw.dual_residual == 0.0);
  CHECK(observation.raw.objective_gap == 0.0);
  CHECK(observation.stopping.metrics.relative_primal == 0.0);
  CHECK(observation.stopping.metrics.relative_dual == 0.0);
  CHECK(observation.stopping.metrics.relative_gap == 0.0);
}

void checkNoConvergence(const hipdlp::algorithm::FixedTraceResult& result) {
  checkNoObservation(result.current_convergence);
  checkNoObservation(result.averaged_convergence);
}

bool bitwiseEqual(const double left, const double right) {
  return std::memcmp(&left, &right, sizeof(double)) == 0;
}

bool bitwiseEqual(const std::vector<double>& left,
                  const std::vector<double>& right) {
  if (left.size() != right.size()) return false;
  if (left.empty()) return true;
  return std::memcmp(left.data(), right.data(), left.size() * sizeof(double)) ==
         0;
}

void checkBitwiseEqualObservation(
    const hipdlp::algorithm::FixedTraceConvergenceObservation& left,
    const hipdlp::algorithm::FixedTraceConvergenceObservation& right) {
  REQUIRE(left.observed == right.observed);
  CHECK(left.stopping.decision == right.stopping.decision);
  CHECK(left.stopping.metrics_valid == right.stopping.metrics_valid);
  CHECK(bitwiseEqual(left.raw.primal_residual, right.raw.primal_residual));
  CHECK(bitwiseEqual(left.raw.dual_residual, right.raw.dual_residual));
  CHECK(bitwiseEqual(left.raw.objective_gap, right.raw.objective_gap));
  CHECK(bitwiseEqual(left.raw.primal_objective, right.raw.primal_objective));
  CHECK(bitwiseEqual(left.raw.dual_objective, right.raw.dual_objective));
  CHECK(bitwiseEqual(left.stopping.metrics.relative_primal,
                     right.stopping.metrics.relative_primal));
  CHECK(bitwiseEqual(left.stopping.metrics.relative_dual,
                     right.stopping.metrics.relative_dual));
  CHECK(bitwiseEqual(left.stopping.metrics.relative_gap,
                     right.stopping.metrics.relative_gap));
}

void checkBitwiseEqual(const hipdlp::algorithm::FixedTraceResult& left,
                       const hipdlp::algorithm::FixedTraceResult& right) {
  REQUIRE(left.points.size() == right.points.size());
  for (std::size_t point = 0; point < left.points.size(); ++point) {
    CHECK(bitwiseEqual(left.points[point].primal, right.points[point].primal));
    CHECK(bitwiseEqual(left.points[point].dual, right.points[point].dual));
    CHECK(bitwiseEqual(left.points[point].ax, right.points[point].ax));
    CHECK(bitwiseEqual(left.points[point].aty, right.points[point].aty));
  }
  CHECK(bitwiseEqual(left.averaged_primal, right.averaged_primal));
  CHECK(bitwiseEqual(left.averaged_dual, right.averaged_dual));
  CHECK(left.averaging.history_valid == right.averaging.history_valid);
  CHECK(left.averaging.sample_count == right.averaging.sample_count);
  CHECK(std::memcmp(&left.averaging.total_weight, &right.averaging.total_weight,
                    sizeof(double)) == 0);
  checkBitwiseEqualObservation(left.current_convergence,
                               right.current_convergence);
  checkBitwiseEqualObservation(left.averaged_convergence,
                               right.averaged_convergence);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
void checkApproxObservation(
    const hipdlp::algorithm::FixedTraceConvergenceObservation& actual,
    const hipdlp::algorithm::FixedTraceConvergenceObservation& expected,
    const double margin) {
  REQUIRE(actual.observed == expected.observed);
  if (!expected.observed) return;
  CHECK(actual.stopping.metrics_valid == expected.stopping.metrics_valid);
  CHECK(actual.stopping.decision == expected.stopping.decision);
  const double actual_values[] = {actual.raw.primal_residual,
                                  actual.raw.dual_residual,
                                  actual.raw.objective_gap,
                                  actual.raw.primal_objective,
                                  actual.raw.dual_objective,
                                  actual.stopping.metrics.relative_primal,
                                  actual.stopping.metrics.relative_dual,
                                  actual.stopping.metrics.relative_gap};
  const double expected_values[] = {expected.raw.primal_residual,
                                    expected.raw.dual_residual,
                                    expected.raw.objective_gap,
                                    expected.raw.primal_objective,
                                    expected.raw.dual_objective,
                                    expected.stopping.metrics.relative_primal,
                                    expected.stopping.metrics.relative_dual,
                                    expected.stopping.metrics.relative_gap};
  for (std::size_t value = 0; value < 8; ++value) {
    CAPTURE(value);
    CHECK(actual_values[value] ==
          Approx(expected_values[value]).margin(margin));
  }
  CHECK(actual.raw.unscaled_rhs_norm == expected.raw.unscaled_rhs_norm);
  CHECK(actual.raw.unscaled_cost_norm == expected.raw.unscaled_cost_norm);
}

void checkApproxTrace(const hipdlp::algorithm::FixedTraceResult& actual,
                      const hipdlp::algorithm::FixedTraceResult& expected,
                      const double margin) {
  REQUIRE(actual.points.size() == expected.points.size());
  for (std::size_t point = 0; point < actual.points.size(); ++point) {
    CAPTURE(point);
    const std::vector<double>* actual_vectors[] = {
        &actual.points[point].primal, &actual.points[point].dual,
        &actual.points[point].ax, &actual.points[point].aty};
    const std::vector<double>* expected_vectors[] = {
        &expected.points[point].primal, &expected.points[point].dual,
        &expected.points[point].ax, &expected.points[point].aty};
    for (std::size_t vector = 0; vector < 4; ++vector) {
      REQUIRE(actual_vectors[vector]->size() ==
              expected_vectors[vector]->size());
      for (std::size_t index = 0; index < actual_vectors[vector]->size();
           ++index) {
        CAPTURE(vector, index);
        CHECK((*actual_vectors[vector])[index] ==
              Approx((*expected_vectors[vector])[index]).margin(margin));
      }
    }
  }
  CHECK(actual.averaging.history_valid == expected.averaging.history_valid);
  CHECK(actual.averaging.sample_count == expected.averaging.sample_count);
  CHECK(actual.averaging.total_weight == expected.averaging.total_weight);
  CHECK(actual.averaging.hasAverage() == expected.averaging.hasAverage());
  const std::vector<double>* actual_averages[] = {&actual.averaged_primal,
                                                  &actual.averaged_dual};
  const std::vector<double>* expected_averages[] = {&expected.averaged_primal,
                                                    &expected.averaged_dual};
  for (std::size_t vector = 0; vector < 2; ++vector) {
    REQUIRE(actual_averages[vector]->size() ==
            expected_averages[vector]->size());
    for (std::size_t index = 0; index < actual_averages[vector]->size();
         ++index) {
      CAPTURE(vector, index);
      CHECK((*actual_averages[vector])[index] ==
            Approx((*expected_averages[vector])[index]).margin(margin));
    }
  }
  checkApproxObservation(actual.current_convergence,
                         expected.current_convergence, margin);
  checkApproxObservation(actual.averaged_convergence,
                         expected.averaged_convergence, margin);
}
#endif

hipdlp::algorithm::FixedTraceResult runTrace(
    const HighsLp& lp, const hipdlp::algorithm::FixedTraceOptions& options) {
  return hipdlp::algorithm::runFixedTrace(hipdlp::model::OriginalModel(lp),
                                          options);
}

}  // namespace

TEST_CASE("hi-pdlp-fixed-trace-one-step",
          "[pdlp][hipdlp][phase3][phase3c][trace][one-step]") {
  const HighsLp lp = makeOneStepLp();
  hipdlp::algorithm::FixedTraceOptions options;
  options.backend = hipdlp::backend::Kind::kCpu;
  options.scaling.mode = hipdlp::model::ScalingMode::kOff;
  options.primal_step_size = 0.5;
  options.dual_step_size = 0.25;
  options.iteration_count = 1;
  options.scaled_primal_start = {0.0};
  options.scaled_dual_start = {0.0};

  const hipdlp::algorithm::FixedTraceResult result =
      hipdlp::algorithm::runFixedTrace(hipdlp::model::OriginalModel(lp),
                                       options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 2);
  CHECK(result.points[0].primal == std::vector<double>({0.0}));
  CHECK(result.points[0].dual == std::vector<double>({0.0}));
  CHECK(result.points[0].ax == std::vector<double>({0.0}));
  CHECK(result.points[0].aty == std::vector<double>({0.0}));
  CHECK(result.points[1].primal == std::vector<double>({0.5}));
  CHECK(result.points[1].dual == std::vector<double>({0.25}));
  CHECK(result.points[1].ax == std::vector<double>({0.5}));
  CHECK(result.points[1].aty == std::vector<double>({0.25}));
}

TEST_CASE("hi-pdlp-fixed-trace-one-step-weighted-average",
          "[pdlp][hipdlp][phase7b][averaging][integration][weights]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
  options.primal_step_size = 1.0;
  options.dual_step_size = 1.0;

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.averaging.history_valid);
  CHECK(result.averaging.sample_count == 1);
  CHECK(result.averaging.total_weight == 1.0);
  CHECK(result.averaging.hasAverage());
  CHECK(result.averaged_primal == std::vector<double>({1.0}));
  CHECK(result.averaged_dual == std::vector<double>({0.0}));
}

TEST_CASE("hi-pdlp-fixed-trace-variable-step-weighted-average",
          "[pdlp][hipdlp][phase7b][averaging][integration][weights]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
  options.scaled_primal_start = {0.0};
  options.scaled_dual_start.clear();
  options.primal_step_sizes = {0.25, 2.25};
  options.dual_step_sizes = {4.0, 4.0};

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeFreeColumnLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.averaging.history_valid);
  CHECK(result.averaging.sample_count == 2);
  CHECK(result.averaging.total_weight == 4.0);
  CHECK(result.averaging.hasAverage());
  REQUIRE(result.averaged_primal.size() == 1);
  CHECK(result.averaged_primal[0] == Approx(39.0 / 16.0).margin(1e-15));
  CHECK(result.averaged_dual.empty());
}

TEST_CASE("hi-pdlp-fixed-trace-averaging-empty-and-transactional",
          "[pdlp][hipdlp][phase7b][averaging][integration][failure]") {
  SECTION("zero iterations produce no valid average") {
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), makeOptions(0));
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.points.size() == 1);
    CHECK(result.averaging.history_valid);
    CHECK(result.averaging.sample_count == 0);
    CHECK(result.averaging.total_weight == 0.0);
    checkNoAverage(result);
  }

  SECTION("scalar failure discards an earlier weighted sample") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    options.primal_step_sizes = {1.0,
                                 std::numeric_limits<double>::denorm_min()};
    options.dual_step_sizes = {1.0, 0.5};
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    CHECK_FALSE(result.ok());
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.points.empty());
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
    checkNoAverage(result);
  }

  SECTION("backend failure discards an earlier weighted sample") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    options.scaled_primal_start = {1e200};
    options.scaled_dual_start.clear();
    options.primal_step_sizes = {1.0, 1e154};
    options.dual_step_sizes = {1.0, 1e154};
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeFreeColumnLp(), options);
    CHECK_FALSE(result.ok());
    CHECK_FALSE(result.error.ok());
    CHECK(result.points.empty());
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
    checkNoAverage(result);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-final-current-convergence",
          "[pdlp][hipdlp][phase7a][convergence][integration][current]") {
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), makeOptions(1));
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 2);
  checkExactPoint(result.points[1], {0.5}, {0.25}, {0.5}, {0.25});

  const hipdlp::algorithm::FixedTraceConvergenceObservation& current =
      result.current_convergence;
  REQUIRE(current.observed);
  REQUIRE(current.stopping.ok());
  CHECK(current.raw.primal_residual == 1.5);
  CHECK(current.raw.dual_residual == 0.0);
  CHECK(current.raw.primal_objective == -0.5);
  CHECK(current.raw.dual_objective == -2.0);
  CHECK(current.raw.objective_gap == 1.5);
  CHECK(current.raw.unscaled_rhs_norm == 2.0);
  CHECK(current.raw.unscaled_cost_norm == 1.0);
  CHECK(current.stopping.metrics_valid);
  CHECK(current.stopping.metrics.relative_primal == 0.5);
  CHECK(current.stopping.metrics.relative_dual == 0.0);
  CHECK(current.stopping.metrics.relative_gap == 1.5 / 3.5);
  CHECK(current.stopping.decision ==
        hipdlp::algorithm::StoppingDecision::kContinue);
}

TEST_CASE("hi-pdlp-fixed-trace-averaged-convergence",
          "[pdlp][hipdlp][phase7a][convergence][integration][averaging]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
  options.primal_step_size = 1.0;
  options.dual_step_size = 1.0;
  options.scaled_dual_start.clear();
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeFreeColumnLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 3);
  CHECK(result.points[2].primal == std::vector<double>({2.0}));
  REQUIRE(result.averaging.hasAverage());
  CHECK(result.averaging.sample_count == 2);
  CHECK(result.averaged_primal == std::vector<double>({1.5}));

  const hipdlp::algorithm::FixedTraceConvergenceObservation& current =
      result.current_convergence;
  REQUIRE(current.observed);
  REQUIRE(current.stopping.ok());
  CHECK(current.raw.primal_residual == 0.0);
  CHECK(current.raw.dual_residual == 1.0);
  CHECK(current.raw.primal_objective == -2.0);
  CHECK(current.raw.dual_objective == 0.0);
  CHECK(current.raw.objective_gap == 2.0);
  CHECK(current.stopping.metrics.relative_primal == 0.0);
  CHECK(current.stopping.metrics.relative_dual == 0.5);
  CHECK(current.stopping.metrics.relative_gap == 2.0 / 3.0);

  const hipdlp::algorithm::FixedTraceConvergenceObservation& averaged =
      result.averaged_convergence;
  REQUIRE(averaged.observed);
  REQUIRE(averaged.stopping.ok());
  CHECK(averaged.raw.primal_residual == 0.0);
  CHECK(averaged.raw.dual_residual == 1.0);
  CHECK(averaged.raw.primal_objective == -1.5);
  CHECK(averaged.raw.dual_objective == 0.0);
  CHECK(averaged.raw.objective_gap == 1.5);
  CHECK(averaged.raw.unscaled_rhs_norm == 0.0);
  CHECK(averaged.raw.unscaled_cost_norm == 1.0);
  CHECK(averaged.stopping.metrics.relative_primal == 0.0);
  CHECK(averaged.stopping.metrics.relative_dual == 0.5);
  CHECK(averaged.stopping.metrics.relative_gap == 0.6);
  CHECK(averaged.stopping.decision ==
        hipdlp::algorithm::StoppingDecision::kContinue);
}

TEST_CASE("hi-pdlp-fixed-trace-convergence-empty-and-transactional",
          "[pdlp][hipdlp][phase7a][convergence][integration][failure]") {
  SECTION("zero iterations observe only the current iterate") {
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), makeOptions(0));
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.points.size() == 1);
    checkNoAverage(result);
    checkNoObservation(result.averaged_convergence);
    REQUIRE(result.current_convergence.observed);
    REQUIRE(result.current_convergence.stopping.ok());
    CHECK(result.current_convergence.raw.primal_residual == 2.0);
    CHECK(result.current_convergence.raw.dual_residual == 0.0);
    CHECK(result.current_convergence.raw.primal_objective == 0.0);
    CHECK(result.current_convergence.raw.dual_objective == -2.0);
    CHECK(result.current_convergence.stopping.metrics.relative_primal ==
          2.0 / 3.0);
    CHECK(result.current_convergence.stopping.metrics.relative_gap ==
          2.0 / 3.0);
  }

  SECTION("an optimal iterate is reported without selecting between iterates") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    options.scaled_primal_start = {2.0};
    options.scaled_dual_start = {-1.0};
    options.primal_step_size = 1e-12;
    options.dual_step_size = 1e-12;
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.current_convergence.observed);
    REQUIRE(result.averaged_convergence.observed);
    CHECK(result.current_convergence.stopping.decision ==
          hipdlp::algorithm::StoppingDecision::kOptimal);
    CHECK(result.averaged_convergence.stopping.decision ==
          hipdlp::algorithm::StoppingDecision::kOptimal);
  }

  SECTION("a scalar failure exposes no convergence observation") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    options.primal_step_sizes = {1.0,
                                 std::numeric_limits<double>::denorm_min()};
    options.dual_step_sizes = {1.0, 0.5};
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    CHECK_FALSE(result.ok());
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.points.empty());
    checkNoConvergence(result);
  }

  SECTION("an invalid tolerance exposes no convergence observation") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    options.convergence_tolerances.dual = 0.0;
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    CHECK_FALSE(result.ok());
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.points.empty());
    checkNoAverage(result);
    checkNoConvergence(result);
  }

  SECTION("a backend failure exposes no convergence observation") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    options.scaled_primal_start = {1e200};
    options.scaled_dual_start.clear();
    options.primal_step_sizes = {1.0, 1e154};
    options.dual_step_sizes = {1.0, 1e154};
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeFreeColumnLp(), options);
    CHECK_FALSE(result.ok());
    CHECK_FALSE(result.error.ok());
    CHECK(result.points.empty());
    checkNoConvergence(result);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-exact-counts",
          "[pdlp][hipdlp][phase3][phase3c][trace]") {
  const HighsLp lp = makeOneStepLp();
  const std::vector<std::vector<double>> primal = {
      {0.0}, {0.5}, {7.0 / 6.0}, {1.8125000000000002}, {1.75}};
  const std::vector<std::vector<double>> dual = {
      {0.0}, {0.25}, {0.25}, {0.03124999999999989}, {-0.05000000000000009}};

  for (const HighsInt iteration_count : {0, 1, 2, 4}) {
    CAPTURE(iteration_count);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(lp, makeOptions(iteration_count));
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.points.size() ==
            static_cast<std::size_t>(iteration_count) + 1);
    for (HighsInt point = 0; point <= iteration_count; ++point) {
      checkExactPoint(result.points[point], primal[point], dual[point],
                      primal[point], dual[point]);
    }
  }
}

TEST_CASE("hi-pdlp-fixed-trace-projections",
          "[pdlp][hipdlp][phase3][phase3c][trace]") {
  SECTION("equality duals are free and inequality duals are nonnegative") {
    hipdlp::algorithm::FixedTraceOptions options;
    options.backend = hipdlp::backend::Kind::kCpu;
    options.scaling.mode = hipdlp::model::ScalingMode::kOff;
    options.primal_step_size = 1.0;
    options.dual_step_size = 1.0;
    options.iteration_count = 1;
    options.scaled_primal_start = {0.0};
    options.scaled_dual_start = {-2.0, -2.0};
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeProjectionLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    checkExactPoint(result.points[0], {0.0}, {-2.0, -2.0}, {0.0, 0.0}, {0.0});
    checkExactPoint(result.points[1], {0.0}, {-1.0, 0.0}, {0.0, 0.0}, {0.0});
  }

  SECTION("fixed free and boxed columns use their distinct projections") {
    hipdlp::algorithm::FixedTraceOptions options;
    options.backend = hipdlp::backend::Kind::kCpu;
    options.scaling.mode = hipdlp::model::ScalingMode::kOff;
    options.primal_step_size = 1.0;
    options.dual_step_size = 1.0;
    options.iteration_count = 1;
    options.scaled_primal_start = {0.0, 0.0, 0.0};
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeColumnBoundsLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    checkExactPoint(result.points[0], {1.0, 0.0, 0.0}, {}, {}, {0.0, 0.0, 0.0});
    checkExactPoint(result.points[1], {1.0, -1.0, 0.5}, {}, {},
                    {0.0, 0.0, 0.0});
  }
}

TEST_CASE("hi-pdlp-fixed-trace-halpern-anchor-weights",
          "[pdlp][hipdlp][phase3][phase3c][trace]") {
  hipdlp::algorithm::FixedTraceOptions options;
  options.backend = hipdlp::backend::Kind::kCpu;
  options.scaling.mode = hipdlp::model::ScalingMode::kOff;
  options.primal_step_size = 0.5;
  options.dual_step_size = 1.0;
  options.iteration_count = 2;
  options.scaled_primal_start = {2.0};
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeFreeColumnLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  checkExactPoint(result.points[0], {2.0}, {}, {}, {0.0});
  checkExactPoint(result.points[1], {2.5}, {}, {}, {0.0});
  checkExactPoint(result.points[2], {3.0}, {}, {}, {0.0});
}

TEST_CASE("hi-pdlp-fixed-trace-determinism-and-loop-audit",
          "[pdlp][hipdlp][phase3][phase3c][trace]") {
  const HighsInt iteration_count = 8;
  const hipdlp::algorithm::FixedTraceResult first =
      runTrace(makeOneStepLp(), makeOptions(iteration_count));
  const hipdlp::algorithm::FixedTraceResult second =
      runTrace(makeOneStepLp(), makeOptions(iteration_count));
  REQUIRE(first.ok());
  REQUIRE(second.ok());
  checkBitwiseEqual(first, second);
  CHECK(first.diagnostics.resource_creation_calls_before_iterations > 0);
  CHECK(first.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(first.diagnostics.transfer_calls_before_iterations == 2);
  CHECK(first.diagnostics.transfer_calls_during_iterations == 0);
  CHECK(first.diagnostics.transfer_calls_after_iterations ==
        4 * (static_cast<std::size_t>(iteration_count) + 1) + 4);
}

TEST_CASE("hi-pdlp-fixed-trace-scaling-pipeline",
          "[pdlp][hipdlp][phase3][phase3b][phase3c][trace][scaling]") {
  struct Selection {
    hipdlp::model::ScalingMode mode;
    HighsInt features_off;
    std::vector<hipdlp::model::ScalingStage> stages;
  };
  const std::vector<Selection> selections = {
      {hipdlp::model::ScalingMode::kAuto,
       kPdlpAllFeaturesOn,
       {hipdlp::model::ScalingStage::kCurtisReid,
        hipdlp::model::ScalingStage::kRuiz,
        hipdlp::model::ScalingStage::kPockChambolle,
        hipdlp::model::ScalingStage::kL2}},
      {hipdlp::model::ScalingMode::kCurtisReid,
       kPdlpAllFeaturesOn,
       {hipdlp::model::ScalingStage::kCurtisReid}},
      {hipdlp::model::ScalingMode::kLegacy,
       kPdlpAllFeaturesOn,
       {hipdlp::model::ScalingStage::kRuiz,
        hipdlp::model::ScalingStage::kPockChambolle,
        hipdlp::model::ScalingStage::kL2}},
      {hipdlp::model::ScalingMode::kOff, kPdlpAllFeaturesOn, {}},
      {hipdlp::model::ScalingMode::kAuto, kPdlpScalingOff, {}}};

  for (const Selection& selection : selections) {
    CAPTURE(static_cast<int>(selection.mode), selection.features_off);
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    options.scaling.mode = selection.mode;
    options.scaling.features_off = selection.features_off;
    options.scaling.legacy_scaling_mode = kPdlpScalingMax;
    options.scaling.ruiz_iterations = 2;
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeScalingLp(), options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeScalingLp(), options);
    INFO(first.error.message);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    CHECK(first.scaling.stages == selection.stages);
    CHECK(first.scaling.applied == !selection.stages.empty());
    checkBitwiseEqual(first, second);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-transactional-errors",
          "[pdlp][hipdlp][phase3][phase3c][trace]") {
  const HighsLp lp = makeOneStepLp();

  SECTION("invalid counts steps dimensions and values") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    options.iteration_count = -1;
    hipdlp::algorithm::FixedTraceResult result = runTrace(lp, options);
    CHECK_FALSE(result.ok());
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.points.empty());

    options = makeOptions(1);
    options.primal_step_size = std::numeric_limits<double>::quiet_NaN();
    result = runTrace(lp, options);
    CHECK_FALSE(result.ok());
    CHECK(result.points.empty());

    options = makeOptions(1);
    options.dual_step_size = std::numeric_limits<double>::infinity();
    result = runTrace(lp, options);
    CHECK_FALSE(result.ok());
    CHECK(result.points.empty());

    options = makeOptions(1);
    options.primal_step_size = -1.0;
    result = runTrace(lp, options);
    CHECK_FALSE(result.ok());
    CHECK(result.points.empty());

    options = makeOptions(1);
    options.scaled_primal_start.clear();
    result = runTrace(lp, options);
    CHECK_FALSE(result.ok());
    CHECK(result.points.empty());

    options = makeOptions(1);
    options.scaled_dual_start[0] = std::numeric_limits<double>::quiet_NaN();
    result = runTrace(lp, options);
    CHECK_FALSE(result.ok());
    CHECK(result.points.empty());
  }

  SECTION("invalid models and unsupported backends") {
    HighsLp invalid = lp;
    invalid.a_matrix_.start_ = {0};
    hipdlp::algorithm::FixedTraceResult result =
        runTrace(invalid, makeOptions(1));
    CHECK_FALSE(result.ok());
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidModel);
    CHECK(result.points.empty());

    invalid = lp;
    invalid.col_cost_[0] = std::numeric_limits<double>::infinity();
    result = runTrace(invalid, makeOptions(1));
    CHECK_FALSE(result.ok());
    CHECK(result.points.empty());

    std::vector<hipdlp::backend::Kind> unavailable = {
        hipdlp::backend::Kind::kOpenmp, hipdlp::backend::Kind::kRocm};
#ifndef HIGHS_HIPDLP_CUDA_DETECTED
    unavailable.push_back(hipdlp::backend::Kind::kCuda);
#endif
    for (const hipdlp::backend::Kind kind : unavailable) {
      hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
      options.backend = kind;
      result = runTrace(lp, options);
      CHECK_FALSE(result.ok());
      CHECK(result.points.empty());
      CHECK(result.error.code != hipdlp::api::ErrorCode::kNone);
    }
  }

  SECTION("backend failure exposes no partial trace") {
    HighsLp overflow = makeOneStepLp();
    overflow.col_lower_[0] = -kHighsInf;
    overflow.col_upper_[0] = kHighsInf;
    overflow.row_lower_[0] = 0.0;
    overflow.row_upper_[0] = 0.0;
    overflow.a_matrix_.value_[0] = std::numeric_limits<double>::max();
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(overflow, options);
    CHECK_FALSE(result.ok());
    CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendFailure);
    CHECK(result.points.empty());
    CHECK(result.diagnostics.resource_creation_calls_during_iterations == 0);
    CHECK(result.diagnostics.transfer_calls_during_iterations == 0);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-openmp-remains-transactional",
          "[pdlp][hipdlp][phase4][trace][openmp]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
  options.backend = hipdlp::backend::Kind::kOpenmp;

  const hipdlp::algorithm::FixedTraceResult first =
      runTrace(makeOneStepLp(), options);
  CHECK_FALSE(first.ok());
  CHECK(first.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
  CHECK(first.error.message == "requested HiPDLP backend is unavailable");
  CHECK(first.points.empty());
  CHECK(first.diagnostics.resource_creation_calls_before_iterations == 0);
  CHECK(first.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(first.diagnostics.transfer_calls_before_iterations == 0);
  CHECK(first.diagnostics.transfer_calls_during_iterations == 0);
  CHECK(first.diagnostics.transfer_calls_after_iterations == 0);

  const hipdlp::algorithm::FixedTraceResult second =
      runTrace(makeOneStepLp(), options);
  CHECK_FALSE(second.ok());
  CHECK(second.error.code == first.error.code);
  CHECK(second.error.message == first.error.message);
  CHECK(second.points.empty());
}

TEST_CASE("hi-pdlp-fixed-trace-rocm-remains-transactional",
          "[pdlp][hipdlp][phase6][rocm][trace]") {
  const hipdlp::algorithm::FixedTraceExecutionMode execution_modes[] = {
      hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded};
  for (const hipdlp::algorithm::FixedTraceExecutionMode execution_mode :
       execution_modes) {
    CAPTURE(static_cast<int>(execution_mode));
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    options.backend = hipdlp::backend::Kind::kRocm;
    options.execution_mode = execution_mode;

    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), options);
    CHECK_FALSE(first.ok());
    CHECK(first.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK(first.error.message == "requested HiPDLP backend is unavailable");
    CHECK(first.points.empty());
    CHECK(first.diagnostics.resource_creation_calls_before_iterations == 0);
    CHECK(first.diagnostics.resource_creation_calls_during_iterations == 0);
    CHECK(first.diagnostics.transfer_calls_before_iterations == 0);
    CHECK(first.diagnostics.transfer_calls_during_iterations == 0);
    CHECK(first.diagnostics.transfer_calls_after_iterations == 0);
    CHECK(first.diagnostics.command_buffer.replay_calls == 0);
    CHECK(first.diagnostics.command_buffer.direct_operation_enqueue_calls == 0);
    CHECK(first.diagnostics.command_buffer.allocation_calls == 0);
    CHECK(first.diagnostics.command_buffer.descriptor_creation_calls == 0);
    CHECK(first.diagnostics.command_buffer.host_to_device_transfer_calls == 0);
    CHECK(first.diagnostics.command_buffer.device_to_host_transfer_calls == 0);
    CHECK(first.diagnostics.command_buffer.synchronization_calls == 0);

    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), options);
    CHECK_FALSE(second.ok());
    CHECK(second.error.code == first.error.code);
    CHECK(second.error.message == first.error.message);
    CHECK(second.points.empty());
    CHECK(second.diagnostics.resource_creation_calls_before_iterations == 0);
    CHECK(second.diagnostics.resource_creation_calls_during_iterations == 0);
    CHECK(second.diagnostics.transfer_calls_before_iterations == 0);
    CHECK(second.diagnostics.transfer_calls_during_iterations == 0);
    CHECK(second.diagnostics.transfer_calls_after_iterations == 0);
    CHECK(second.diagnostics.command_buffer.replay_calls == 0);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-recording-requires-capability",
          "[pdlp][hipdlp][phase5g][cpu][recording]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
  options.backend = hipdlp::backend::Kind::kCpu;
  options.execution_mode =
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  CHECK_FALSE(result.ok());
  CHECK(result.error.code == hipdlp::api::ErrorCode::kUnsupportedCapability);
  CHECK(result.points.empty());
  CHECK(result.diagnostics.resource_creation_calls_before_iterations == 0);
  CHECK(result.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(result.diagnostics.transfer_calls_before_iterations == 0);
  CHECK(result.diagnostics.transfer_calls_during_iterations == 0);
  CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-shared-step-five-iteration-parity",
          "[pdlp][hipdlp][phase5h][algorithm][shared-step]") {
  hipdlp::algorithm::FixedTraceOptions cpu_options = makeOptions(5);
  hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
  direct_options.backend = hipdlp::backend::Kind::kCuda;
  hipdlp::algorithm::FixedTraceOptions recorded_options = direct_options;
  recorded_options.execution_mode =
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;

  const hipdlp::algorithm::FixedTraceResult cpu =
      runTrace(makeOneStepLp(), cpu_options);
  const hipdlp::algorithm::FixedTraceResult direct =
      runTrace(makeOneStepLp(), direct_options);
  hipdlp::backend::cuda::testing::resetResourceAudit();
  const hipdlp::algorithm::FixedTraceResult recorded =
      runTrace(makeOneStepLp(), recorded_options);
  const hipdlp::backend::cuda::testing::ResourceAudit recorded_audit =
      hipdlp::backend::cuda::testing::resourceAudit();
  INFO(cpu.error.message);
  INFO(direct.error.message);
  INFO(recorded.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(direct.ok());
  REQUIRE(recorded.ok());
  checkApproxTrace(direct, cpu, 5e-13);
  checkApproxTrace(recorded, direct, 5e-13);
  CHECK(recorded.diagnostics.command_buffer.replay_calls == 3);
  CHECK(recorded_audit.graph_launch_calls == 3);
}

TEST_CASE("hi-pdlp-cuda-graph-fixed-trace-basic-parity",
          "[pdlp][hipdlp][phase5g][phase7a][phase7b][averaging][convergence]"
          "[cuda][graph][trace][parity]") {
  const std::size_t expected_launches[] = {0, 1, 1, 2, 3};
  const HighsInt iteration_counts[] = {0, 1, 2, 4, 5};
  const std::vector<double> primal_steps = {0.25, 2.25, 1.0, 4.0, 0.5};
  const std::vector<double> dual_steps = {4.0, 4.0, 4.0, 1.0, 8.0};
  for (std::size_t selection = 0; selection < 5; ++selection) {
    const HighsInt iteration_count = iteration_counts[selection];
    CAPTURE(iteration_count);
    hipdlp::algorithm::FixedTraceOptions cpu_options =
        makeOptions(iteration_count);
    if (iteration_count > 0) {
      cpu_options.primal_step_sizes.assign(
          primal_steps.begin(), primal_steps.begin() + iteration_count);
      cpu_options.dual_step_sizes.assign(dual_steps.begin(),
                                         dual_steps.begin() + iteration_count);
    }
    hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
    direct_options.backend = hipdlp::backend::Kind::kCuda;
    hipdlp::algorithm::FixedTraceOptions graph_options = direct_options;
    graph_options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    const hipdlp::algorithm::FixedTraceResult cpu =
        runTrace(makeOneStepLp(), cpu_options);
    const hipdlp::algorithm::FixedTraceResult direct =
        runTrace(makeOneStepLp(), direct_options);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult graph =
        runTrace(makeOneStepLp(), graph_options);
    const hipdlp::backend::cuda::testing::ResourceAudit graph_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    INFO(cpu.error.message);
    INFO(direct.error.message);
    INFO(graph.error.message);
    REQUIRE(cpu.ok());
    REQUIRE(direct.ok());
    REQUIRE(graph.ok());
    checkApproxTrace(graph, direct, 5e-13);
    checkApproxTrace(graph, cpu, 5e-13);
    CHECK(graph.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(graph_audit.graph_launch_calls == expected_launches[selection]);
    CHECK(graph.diagnostics.command_buffer.direct_operation_enqueue_calls == 0);
    CHECK(graph.diagnostics.command_buffer.allocation_calls == 0);
    CHECK(graph.diagnostics.command_buffer.descriptor_creation_calls == 0);
    CHECK(graph.diagnostics.command_buffer.host_to_device_transfer_calls == 0);
    CHECK(graph.diagnostics.command_buffer.device_to_host_transfer_calls == 0);
    CHECK(graph.diagnostics.command_buffer.synchronization_calls == 0);
  }
}

TEST_CASE("hi-pdlp-cuda-graph-fixed-trace-scaling-parity",
          "[pdlp][hipdlp][phase5g][cuda][graph][trace][parity][scaling]") {
  struct Selection {
    hipdlp::model::ScalingMode mode;
    HighsInt features_off;
  };
  const Selection selections[] = {
      {hipdlp::model::ScalingMode::kAuto, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kCurtisReid, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kLegacy, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kOff, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kAuto, kPdlpScalingOff}};
  for (const Selection& selection : selections) {
    for (const HighsInt iteration_count : {0, 1, 2, 4}) {
      CAPTURE(static_cast<int>(selection.mode), selection.features_off,
              iteration_count);
      hipdlp::algorithm::FixedTraceOptions cpu_options =
          makeOptions(iteration_count);
      cpu_options.scaling.mode = selection.mode;
      cpu_options.scaling.features_off = selection.features_off;
      cpu_options.scaling.legacy_scaling_mode = kPdlpScalingMax;
      cpu_options.scaling.ruiz_iterations = 2;
      hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
      direct_options.backend = hipdlp::backend::Kind::kCuda;
      hipdlp::algorithm::FixedTraceOptions graph_options = direct_options;
      graph_options.execution_mode =
          hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
      const hipdlp::algorithm::FixedTraceResult cpu =
          runTrace(makeScalingLp(), cpu_options);
      const hipdlp::algorithm::FixedTraceResult direct =
          runTrace(makeScalingLp(), direct_options);
      const hipdlp::algorithm::FixedTraceResult graph =
          runTrace(makeScalingLp(), graph_options);
      INFO(cpu.error.message);
      INFO(direct.error.message);
      INFO(graph.error.message);
      REQUIRE(cpu.ok());
      REQUIRE(direct.ok());
      REQUIRE(graph.ok());
      CHECK(graph.scaling.stages == direct.scaling.stages);
      CHECK(graph.scaling.stages == cpu.scaling.stages);
      checkApproxTrace(graph, direct, 5e-13);
      checkApproxTrace(graph, cpu, 5e-13);
      CHECK(graph.diagnostics.command_buffer.replay_calls ==
            static_cast<std::size_t>((iteration_count + 1) / 2));
    }
  }
}

TEST_CASE("hi-pdlp-cuda-graph-fixed-trace-structural-parity",
          "[pdlp][hipdlp][phase5g][cuda][graph][trace][parity][structure]") {
  struct Fixture {
    HighsLp lp;
    std::vector<double> primal_start;
    std::vector<double> dual_start;
    double primal_step;
    double dual_step;
  };
  const Fixture fixtures[] = {
      {makeProjectionLp(), {0.0}, {-2.0, -2.0}, 1.0, 1.0},
      {makeColumnBoundsLp(), {0.0, 0.0, 0.0}, {}, 1.0, 1.0},
      {makeRectangularLp(false), {0.0, 0.0, 1.0}, {0.0, 0.0}, 0.05, 0.05},
      {makeRectangularLp(true), {0.0, 0.0, 1.0}, {0.0, 0.0}, 1e-12, 1e-12},
      {makeEmptyLp(), {}, {}, 1.0, 1.0}};
  for (std::size_t fixture = 0;
       fixture < sizeof(fixtures) / sizeof(fixtures[0]); ++fixture) {
    for (const HighsInt iteration_count : {0, 1, 2, 4}) {
      CAPTURE(fixture, iteration_count);
      hipdlp::algorithm::FixedTraceOptions cpu_options;
      cpu_options.backend = hipdlp::backend::Kind::kCpu;
      cpu_options.scaling.mode = hipdlp::model::ScalingMode::kOff;
      cpu_options.primal_step_size = fixtures[fixture].primal_step;
      cpu_options.dual_step_size = fixtures[fixture].dual_step;
      cpu_options.iteration_count = iteration_count;
      cpu_options.scaled_primal_start = fixtures[fixture].primal_start;
      cpu_options.scaled_dual_start = fixtures[fixture].dual_start;
      hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
      direct_options.backend = hipdlp::backend::Kind::kCuda;
      hipdlp::algorithm::FixedTraceOptions graph_options = direct_options;
      graph_options.execution_mode =
          hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
      const hipdlp::algorithm::FixedTraceResult cpu =
          runTrace(fixtures[fixture].lp, cpu_options);
      const hipdlp::algorithm::FixedTraceResult direct =
          runTrace(fixtures[fixture].lp, direct_options);
      const hipdlp::algorithm::FixedTraceResult graph =
          runTrace(fixtures[fixture].lp, graph_options);
      INFO(cpu.error.message);
      INFO(direct.error.message);
      INFO(graph.error.message);
      REQUIRE(cpu.ok());
      REQUIRE(direct.ok());
      REQUIRE(graph.ok());
      checkApproxTrace(graph, direct, 5e-13);
      checkApproxTrace(graph, cpu, 5e-13);
      CHECK(graph.diagnostics.command_buffer.replay_calls ==
            static_cast<std::size_t>((iteration_count + 1) / 2));
    }
  }
}

TEST_CASE("hi-pdlp-cuda-graph-fixed-trace-one-step-remainder",
          "[pdlp][hipdlp][phase5g][cuda][graph][trace][remainder]") {
  hipdlp::algorithm::FixedTraceOptions cpu_options = makeOptions(5);
  hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
  direct_options.backend = hipdlp::backend::Kind::kCuda;
  hipdlp::algorithm::FixedTraceOptions graph_options = direct_options;
  graph_options.execution_mode =
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
  const hipdlp::algorithm::FixedTraceResult cpu =
      runTrace(makeOneStepLp(), cpu_options);
  const hipdlp::algorithm::FixedTraceResult direct =
      runTrace(makeOneStepLp(), direct_options);
  const hipdlp::algorithm::FixedTraceResult graph =
      runTrace(makeOneStepLp(), graph_options);
  INFO(cpu.error.message);
  INFO(direct.error.message);
  INFO(graph.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(direct.ok());
  REQUIRE(graph.ok());
  checkApproxTrace(graph, direct, 5e-13);
  checkApproxTrace(graph, cpu, 5e-13);
  CHECK(graph.diagnostics.command_buffer.replay_calls == 3);
}

TEST_CASE("hi-pdlp-cuda-graph-fixed-trace-failures-are-transactional",
          "[pdlp][hipdlp][phase5g][cuda][graph][trace][failure]") {
  const hipdlp::backend::cuda::testing::FailurePoint failure_points[] = {
      hipdlp::backend::cuda::testing::FailurePoint::kGraphCapture,
      hipdlp::backend::cuda::testing::FailurePoint::kGraphInstantiation,
      hipdlp::backend::cuda::testing::FailurePoint::kGraphLaunch,
      hipdlp::backend::cuda::testing::FailurePoint::
          kGraphAsynchronousExecution};
  for (const hipdlp::backend::cuda::testing::FailurePoint failure_point :
       failure_points) {
    CAPTURE(static_cast<int>(failure_point));
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(4);
    options.backend = hipdlp::backend::Kind::kCuda;
    options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    hipdlp::backend::cuda::testing::setFailurePoint(failure_point);
    const hipdlp::algorithm::FixedTraceResult failed =
        runTrace(makeOneStepLp(), options);
    hipdlp::backend::cuda::testing::setFailurePoint(
        hipdlp::backend::cuda::testing::FailurePoint::kNone);
    CHECK_FALSE(failed.ok());
    CHECK(failed.error.code == hipdlp::api::ErrorCode::kBackendFailure);
    CHECK_FALSE(failed.error.message.empty());
    CHECK(failed.points.empty());
    const hipdlp::backend::cuda::testing::ResourceAudit audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    CHECK(audit.allocation_calls == audit.free_calls);
    CHECK(audit.descriptor_creations == audit.descriptor_destructions);
    CHECK(audit.graph_creations == audit.graph_destructions);
    CHECK(audit.graph_exec_creations == audit.graph_exec_destructions);

    const hipdlp::algorithm::FixedTraceResult recovered =
        runTrace(makeOneStepLp(), options);
    INFO(recovered.error.message);
    REQUIRE(recovered.ok());
  }
}

TEST_CASE("hi-pdlp-cuda-fixed-trace-parity-across-scaling",
          "[pdlp][hipdlp][phase3c][phase5][cuda][trace][parity]") {
  struct Selection {
    hipdlp::model::ScalingMode mode;
    HighsInt features_off;
  };
  const std::vector<Selection> selections = {
      {hipdlp::model::ScalingMode::kAuto, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kCurtisReid, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kLegacy, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kOff, kPdlpAllFeaturesOn},
      {hipdlp::model::ScalingMode::kAuto, kPdlpScalingOff}};

  for (const Selection& selection : selections) {
    for (const HighsInt iteration_count : {0, 1, 2, 4}) {
      CAPTURE(static_cast<int>(selection.mode), selection.features_off,
              iteration_count);
      hipdlp::algorithm::FixedTraceOptions cpu_options =
          makeOptions(iteration_count);
      cpu_options.scaling.mode = selection.mode;
      cpu_options.scaling.features_off = selection.features_off;
      cpu_options.scaling.legacy_scaling_mode = kPdlpScalingMax;
      cpu_options.scaling.ruiz_iterations = 2;
      hipdlp::algorithm::FixedTraceOptions cuda_options = cpu_options;
      cuda_options.backend = hipdlp::backend::Kind::kCuda;
      const hipdlp::algorithm::FixedTraceResult cpu =
          runTrace(makeScalingLp(), cpu_options);
      const hipdlp::algorithm::FixedTraceResult cuda =
          runTrace(makeScalingLp(), cuda_options);
      INFO(cpu.error.message);
      INFO(cuda.error.message);
      REQUIRE(cpu.ok());
      REQUIRE(cuda.ok());
      CHECK(cuda.scaling.stages == cpu.scaling.stages);
      checkApproxTrace(cuda, cpu, 5e-13);
      CHECK(cuda.diagnostics.resource_creation_calls_during_iterations == 0);
      CHECK(cuda.diagnostics.transfer_calls_during_iterations == 0);
    }
  }
}

TEST_CASE("hi-pdlp-cuda-fixed-trace-determinism-and-loop-audit",
          "[pdlp][hipdlp][phase3c][phase5][cuda][trace][determinism]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(8);
  options.backend = hipdlp::backend::Kind::kCuda;
  const hipdlp::algorithm::FixedTraceResult first =
      runTrace(makeOneStepLp(), options);
  const hipdlp::algorithm::FixedTraceResult second =
      runTrace(makeOneStepLp(), options);
  INFO(first.error.message);
  INFO(second.error.message);
  REQUIRE(first.ok());
  REQUIRE(second.ok());
  checkBitwiseEqual(first, second);
  CHECK(first.diagnostics.resource_creation_calls_before_iterations > 0);
  CHECK(first.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(first.diagnostics.transfer_calls_before_iterations == 2);
  CHECK(first.diagnostics.transfer_calls_during_iterations == 0);
  CHECK(first.diagnostics.transfer_calls_after_iterations == 40);
}

TEST_CASE("hi-pdlp-cuda-fixed-trace-failure-is-transactional",
          "[pdlp][hipdlp][phase3c][phase5][cuda][trace][failure]") {
  HighsLp overflow = makeOneStepLp();
  overflow.col_lower_[0] = -kHighsInf;
  overflow.col_upper_[0] = kHighsInf;
  overflow.row_lower_[0] = 0.0;
  overflow.row_upper_[0] = 0.0;
  overflow.a_matrix_.value_[0] = std::numeric_limits<double>::max();
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
  options.backend = hipdlp::backend::Kind::kCuda;
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(overflow, options);
  CHECK_FALSE(result.ok());
  CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendFailure);
  CHECK(result.points.empty());
  CHECK(result.diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(result.diagnostics.transfer_calls_during_iterations == 0);
}
#endif

#ifndef CUPDLP_GPU
TEST_CASE("hi-pdlp-fixed-trace-legacy-differential",
          "[pdlp][hipdlp][phase3][phase3c][trace][legacy]") {
  const HighsLp lp = makeOneStepLp();
  const HighsInt iteration_count = 4;
  const hipdlp::algorithm::FixedTraceResult rewrite =
      runTrace(lp, makeOptions(iteration_count));
  REQUIRE(rewrite.ok());

  HighsOptions options;
  options.output_flag = false;
  options.pdlp_features_off = kPdlpScalingOff;
  HighsTimer timer;
  PDLPSolver legacy;
  legacy.setup(options, timer);
  legacy.passLp(&lp);
  legacy.preprocessLp();
  legacy.scaleProblem();
  std::vector<std::vector<double>> primal;
  std::vector<std::vector<double>> dual;
  std::vector<std::vector<double>> ax;
  std::vector<std::vector<double>> aty;
  legacy.runFixedTraceForTesting(0.5, 0.25, iteration_count, primal, dual, ax,
                                 aty);

  REQUIRE(rewrite.points.size() == primal.size());
  REQUIRE(rewrite.points.size() == dual.size());
  REQUIRE(rewrite.points.size() == ax.size());
  REQUIRE(rewrite.points.size() == aty.size());
  for (std::size_t point = 0; point < rewrite.points.size(); ++point) {
    CAPTURE(point);
    REQUIRE(rewrite.points[point].primal.size() == primal[point].size());
    REQUIRE(rewrite.points[point].dual.size() == dual[point].size());
    REQUIRE(rewrite.points[point].ax.size() == ax[point].size());
    REQUIRE(rewrite.points[point].aty.size() == aty[point].size());
    for (std::size_t col = 0; col < primal[point].size(); ++col) {
      CHECK(rewrite.points[point].primal[col] ==
            Approx(primal[point][col]).margin(1e-15));
      CHECK(rewrite.points[point].aty[col] ==
            Approx(aty[point][col]).margin(1e-15));
    }
    for (std::size_t row = 0; row < dual[point].size(); ++row) {
      CHECK(rewrite.points[point].dual[row] ==
            Approx(dual[point][row]).margin(1e-15));
      CHECK(rewrite.points[point].ax[row] ==
            Approx(ax[point][row]).margin(1e-15));
    }
  }
}
#endif