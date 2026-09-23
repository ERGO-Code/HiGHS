#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#include "HConfig.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "lp_data/HighsOptions.h"
#include "pdlp/hipdlp/algorithm/fixed_trace.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
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

// The same iterates as makeOneStepLp while the two free bounds expose a
// nonzero normalized dual residual.
HighsLp makeFreeStepLp() {
  HighsLp lp = makeOneStepLp();
  lp.col_lower_ = {-kHighsInf};
  lp.col_upper_ = {kHighsInf};
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

void checkNoFixedPointErrorObservation(
    const hipdlp::algorithm::FixedTraceFixedPointErrorObservation&
        observation) {
  CHECK_FALSE(observation.observed);
  CHECK(observation.evaluation.error.ok());
  CHECK_FALSE(observation.evaluation.ok());
  CHECK(observation.evaluation.fixed_point_error == 0.0);
  CHECK(observation.evaluation.diagnostics.device_to_host_scalar_readbacks ==
        0);
}

void checkNoFixedPointError(const hipdlp::algorithm::FixedTraceResult& result) {
  checkNoFixedPointErrorObservation(result.fixed_point_error);
}

void checkNoFixedPointRestartObservation(
    const hipdlp::algorithm::FixedTraceFixedPointRestartObservation&
        observation) {
  CHECK_FALSE(observation.observed);
  CHECK(observation.evaluation.error.ok());
  CHECK_FALSE(observation.evaluation.ok());
  CHECK(observation.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kNone);
  CHECK(observation.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNone);
  CHECK_FALSE(observation.evaluation.next.history_valid);
  CHECK_FALSE(observation.evaluation.next.reference_valid);
  CHECK_FALSE(observation.evaluation.next.previous_trial_valid);
  CHECK(observation.evaluation.next.reference_error == 0.0);
  CHECK(observation.evaluation.next.previous_trial_error == 0.0);
}

void checkNoFixedPointRestart(
    const hipdlp::algorithm::FixedTraceResult& result) {
  checkNoFixedPointRestartObservation(result.fixed_point_restart);
  CHECK(result.fixed_point_restart_checks.empty());
}

void checkNoPrimalWeight(const hipdlp::algorithm::FixedTraceResult& result) {
  CHECK(result.primal_weight_updates.empty());
}

void checkTransactionalFailure(
    const hipdlp::algorithm::FixedTraceResult& result) {
  CHECK_FALSE(result.ok());
  CHECK_FALSE(result.error.ok());
  CHECK(result.points.empty());
  CHECK(result.convergence_checks.empty());
  checkNoAverage(result);
  checkNoConvergence(result);
  checkNoFixedPointError(result);
  checkNoFixedPointRestart(result);
  checkNoPrimalWeight(result);
}

void enableFixedPointError(hipdlp::algorithm::FixedTraceOptions& options,
                           const double omega, const double eta) {
  options.fixed_point_error.enabled = true;
  options.fixed_point_error.omega = omega;
  options.fixed_point_error.eta = eta;
}

hipdlp::algorithm::FixedPointRestartState makeRestartReference(
    const double reference_error) {
  hipdlp::algorithm::FixedPointRestartState state;
  state.reference_error = reference_error;
  state.reference_valid = true;
  return state;
}

hipdlp::algorithm::FixedPointRestartState makeRestartHistory(
    const double reference_error, const double previous_trial_error) {
  hipdlp::algorithm::FixedPointRestartState state =
      makeRestartReference(reference_error);
  state.previous_trial_error = previous_trial_error;
  state.previous_trial_valid = true;
  return state;
}

void enableFixedPointRestart(
    hipdlp::algorithm::FixedTraceOptions& options,
    const hipdlp::algorithm::FixedPointRestartState& previous,
    const HighsInt check_interval) {
  options.fixed_point_restart.enabled = true;
  options.fixed_point_restart.previous = previous;
  options.fixed_point_restart.parameters.check_interval = check_interval;
}

void enableFixedPointRestartCadence(
    hipdlp::algorithm::FixedTraceOptions& options,
    const hipdlp::algorithm::FixedPointRestartState& previous,
    const HighsInt check_interval, const HighsInt check_cadence) {
  enableFixedPointRestart(options, previous, check_interval);
  options.fixed_point_restart.check_cadence = check_cadence;
}

void enablePeriodicConvergence(hipdlp::algorithm::FixedTraceOptions& options,
                               const HighsInt check_cadence) {
  options.convergence.periodic_enabled = true;
  options.convergence.check_cadence = check_cadence;
}

void enablePrimalWeight(hipdlp::algorithm::FixedTraceOptions& options,
                        const double primal_weight) {
  options.primal_weight.enabled = true;
  options.primal_weight.previous =
      hipdlp::algorithm::resetPrimalWeight(primal_weight);
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

void checkBitwiseEqualFixedPointError(
    const hipdlp::algorithm::FixedTraceFixedPointErrorObservation& left,
    const hipdlp::algorithm::FixedTraceFixedPointErrorObservation& right) {
  REQUIRE(left.observed == right.observed);
  CHECK(left.evaluation.error.code == right.evaluation.error.code);
  CHECK(left.evaluation.error.message == right.evaluation.error.message);
  CHECK(bitwiseEqual(left.evaluation.fixed_point_error,
                     right.evaluation.fixed_point_error));
  CHECK(left.evaluation.diagnostics.device_to_host_scalar_readbacks ==
        right.evaluation.diagnostics.device_to_host_scalar_readbacks);
}

void checkBitwiseEqualTraceOutputs(
    const hipdlp::algorithm::FixedTraceResult& left,
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

void checkBitwiseEqualFixedPointRestart(
    const hipdlp::algorithm::FixedTraceFixedPointRestartObservation& left,
    const hipdlp::algorithm::FixedTraceFixedPointRestartObservation& right) {
  REQUIRE(left.observed == right.observed);
  CHECK(left.evaluation.error.code == right.evaluation.error.code);
  CHECK(left.evaluation.error.message == right.evaluation.error.message);
  CHECK(left.evaluation.decision == right.evaluation.decision);
  CHECK(left.evaluation.condition == right.evaluation.condition);
  CHECK(left.evaluation.next.history_valid ==
        right.evaluation.next.history_valid);
  CHECK(left.evaluation.next.reference_valid ==
        right.evaluation.next.reference_valid);
  CHECK(left.evaluation.next.previous_trial_valid ==
        right.evaluation.next.previous_trial_valid);
  CHECK(bitwiseEqual(left.evaluation.next.reference_error,
                     right.evaluation.next.reference_error));
  CHECK(bitwiseEqual(left.evaluation.next.previous_trial_error,
                     right.evaluation.next.previous_trial_error));
}

void checkBitwiseEqualFixedPointRestartChecks(
    const std::vector<hipdlp::algorithm::FixedTraceFixedPointRestartCheck>&
        left,
    const std::vector<hipdlp::algorithm::FixedTraceFixedPointRestartCheck>&
        right) {
  REQUIRE(left.size() == right.size());
  for (std::size_t check = 0; check < left.size(); ++check) {
    CAPTURE(check);
    CHECK(left[check].iteration == right[check].iteration);
    CHECK(bitwiseEqual(left[check].error.fixed_point_error,
                       right[check].error.fixed_point_error));
    CHECK(left[check].error.diagnostics.device_to_host_scalar_readbacks ==
          right[check].error.diagnostics.device_to_host_scalar_readbacks);
    CHECK(left[check].evaluation.decision == right[check].evaluation.decision);
    CHECK(left[check].evaluation.condition ==
          right[check].evaluation.condition);
    CHECK(left[check].evaluation.next.reference_valid ==
          right[check].evaluation.next.reference_valid);
    CHECK(left[check].evaluation.next.previous_trial_valid ==
          right[check].evaluation.next.previous_trial_valid);
    CHECK(bitwiseEqual(left[check].evaluation.next.reference_error,
                       right[check].evaluation.next.reference_error));
    CHECK(bitwiseEqual(left[check].evaluation.next.previous_trial_error,
                       right[check].evaluation.next.previous_trial_error));
  }
}

void checkBitwiseEqualConvergenceChecks(
    const std::vector<hipdlp::algorithm::FixedTraceConvergenceCheck>& left,
    const std::vector<hipdlp::algorithm::FixedTraceConvergenceCheck>& right) {
  REQUIRE(left.size() == right.size());
  for (std::size_t check = 0; check < left.size(); ++check) {
    CAPTURE(check);
    CHECK(left[check].iteration == right[check].iteration);
    CHECK(left[check].stopping.error.code == right[check].stopping.error.code);
    CHECK(left[check].stopping.error.message ==
          right[check].stopping.error.message);
    CHECK(left[check].stopping.decision == right[check].stopping.decision);
    CHECK(left[check].stopping.metrics_valid ==
          right[check].stopping.metrics_valid);
    const double left_values[] = {left[check].raw.primal_residual,
                                  left[check].raw.dual_residual,
                                  left[check].raw.objective_gap,
                                  left[check].raw.primal_objective,
                                  left[check].raw.dual_objective,
                                  left[check].raw.unscaled_rhs_norm,
                                  left[check].raw.unscaled_cost_norm,
                                  left[check].stopping.metrics.relative_primal,
                                  left[check].stopping.metrics.relative_dual,
                                  left[check].stopping.metrics.relative_gap};
    const double right_values[] = {
        right[check].raw.primal_residual,
        right[check].raw.dual_residual,
        right[check].raw.objective_gap,
        right[check].raw.primal_objective,
        right[check].raw.dual_objective,
        right[check].raw.unscaled_rhs_norm,
        right[check].raw.unscaled_cost_norm,
        right[check].stopping.metrics.relative_primal,
        right[check].stopping.metrics.relative_dual,
        right[check].stopping.metrics.relative_gap};
    for (std::size_t value = 0; value < 10; ++value) {
      CAPTURE(value);
      CHECK(bitwiseEqual(left_values[value], right_values[value]));
    }
  }
}

void checkBitwiseEqualPrimalWeightResult(
    const hipdlp::algorithm::PrimalWeightUpdateResult& left,
    const hipdlp::algorithm::PrimalWeightUpdateResult& right) {
  CHECK(left.error.code == right.error.code);
  CHECK(left.error.message == right.error.message);
  CHECK(left.kind == right.kind);
  CHECK(left.next.best_gap_valid == right.next.best_gap_valid);
  CHECK(left.next.history_valid == right.next.history_valid);
  const double left_values[] = {
      left.next.primal_weight,     left.next.error_sum,
      left.next.previous_error,    left.next.best_primal_weight,
      left.next.best_residual_gap, left.step_sizes.primal_weight,
      left.step_sizes.base_step,   left.step_sizes.primal_step,
      left.step_sizes.dual_step};
  const double right_values[] = {
      right.next.primal_weight,     right.next.error_sum,
      right.next.previous_error,    right.next.best_primal_weight,
      right.next.best_residual_gap, right.step_sizes.primal_weight,
      right.step_sizes.base_step,   right.step_sizes.primal_step,
      right.step_sizes.dual_step};
  for (std::size_t value = 0; value < 9; ++value) {
    CAPTURE(value);
    CHECK(bitwiseEqual(left_values[value], right_values[value]));
  }
}

void checkBitwiseEqualPrimalWeightUpdates(
    const std::vector<hipdlp::algorithm::FixedTracePrimalWeightUpdate>& left,
    const std::vector<hipdlp::algorithm::FixedTracePrimalWeightUpdate>& right) {
  REQUIRE(left.size() == right.size());
  for (std::size_t update = 0; update < left.size(); ++update) {
    CAPTURE(update);
    CHECK(left[update].iteration == right[update].iteration);
    CHECK(bitwiseEqual(left[update].primal_anchor_distance,
                       right[update].primal_anchor_distance));
    CHECK(bitwiseEqual(left[update].dual_anchor_distance,
                       right[update].dual_anchor_distance));
    checkBitwiseEqualPrimalWeightResult(left[update].result,
                                        right[update].result);
  }
}

void checkBitwiseEqual(const hipdlp::algorithm::FixedTraceResult& left,
                       const hipdlp::algorithm::FixedTraceResult& right) {
  checkBitwiseEqualTraceOutputs(left, right);
  checkBitwiseEqualFixedPointError(left.fixed_point_error,
                                   right.fixed_point_error);
  checkBitwiseEqualFixedPointRestart(left.fixed_point_restart,
                                     right.fixed_point_restart);
  checkBitwiseEqualFixedPointRestartChecks(left.fixed_point_restart_checks,
                                           right.fixed_point_restart_checks);
  checkBitwiseEqualPrimalWeightUpdates(left.primal_weight_updates,
                                       right.primal_weight_updates);
}

void checkSameRestartDecisions(
    const hipdlp::algorithm::FixedTraceResult& left,
    const hipdlp::algorithm::FixedTraceResult& right) {
  REQUIRE(left.fixed_point_restart_checks.size() ==
          right.fixed_point_restart_checks.size());
  for (std::size_t check = 0; check < left.fixed_point_restart_checks.size();
       ++check) {
    CAPTURE(check);
    CHECK(left.fixed_point_restart_checks[check].iteration ==
          right.fixed_point_restart_checks[check].iteration);
    CHECK(left.fixed_point_restart_checks[check].evaluation.decision ==
          right.fixed_point_restart_checks[check].evaluation.decision);
    CHECK(left.fixed_point_restart_checks[check].evaluation.condition ==
          right.fixed_point_restart_checks[check].evaluation.condition);
  }
}

void checkSameDiagnostics(const hipdlp::algorithm::FixedTraceResult& left,
                          const hipdlp::algorithm::FixedTraceResult& right) {
  CHECK(left.diagnostics.resource_creation_calls_before_iterations ==
        right.diagnostics.resource_creation_calls_before_iterations);
  CHECK(left.diagnostics.resource_creation_calls_during_iterations ==
        right.diagnostics.resource_creation_calls_during_iterations);
  CHECK(left.diagnostics.transfer_calls_before_iterations ==
        right.diagnostics.transfer_calls_before_iterations);
  CHECK(left.diagnostics.transfer_calls_during_iterations ==
        right.diagnostics.transfer_calls_during_iterations);
  CHECK(left.diagnostics.transfer_calls_after_iterations ==
        right.diagnostics.transfer_calls_after_iterations);
  CHECK(left.diagnostics.command_buffer.replay_calls ==
        right.diagnostics.command_buffer.replay_calls);
  CHECK(left.diagnostics.command_buffer.direct_operation_enqueue_calls ==
        right.diagnostics.command_buffer.direct_operation_enqueue_calls);
  CHECK(left.diagnostics.command_buffer.allocation_calls ==
        right.diagnostics.command_buffer.allocation_calls);
  CHECK(left.diagnostics.command_buffer.descriptor_creation_calls ==
        right.diagnostics.command_buffer.descriptor_creation_calls);
  CHECK(left.diagnostics.command_buffer.host_to_device_transfer_calls ==
        right.diagnostics.command_buffer.host_to_device_transfer_calls);
  CHECK(left.diagnostics.command_buffer.device_to_host_transfer_calls ==
        right.diagnostics.command_buffer.device_to_host_transfer_calls);
  CHECK(left.diagnostics.command_buffer.synchronization_calls ==
        right.diagnostics.command_buffer.synchronization_calls);
}

void checkPeriodicConvergenceDiagnosticsDelta(
    const hipdlp::algorithm::FixedTraceResult& disabled,
    const hipdlp::algorithm::FixedTraceResult& enabled,
    const std::size_t evaluations) {
  CHECK(enabled.diagnostics.resource_creation_calls_before_iterations ==
        disabled.diagnostics.resource_creation_calls_before_iterations + 2);
  CHECK(enabled.diagnostics.resource_creation_calls_during_iterations ==
        disabled.diagnostics.resource_creation_calls_during_iterations);
  CHECK(enabled.diagnostics.transfer_calls_before_iterations ==
        disabled.diagnostics.transfer_calls_before_iterations);
  CHECK(enabled.diagnostics.transfer_calls_during_iterations ==
        disabled.diagnostics.transfer_calls_during_iterations +
            4 * evaluations);
  CHECK(enabled.diagnostics.transfer_calls_after_iterations ==
        disabled.diagnostics.transfer_calls_after_iterations);
  CHECK(enabled.diagnostics.command_buffer.replay_calls ==
        disabled.diagnostics.command_buffer.replay_calls);
  CHECK(enabled.diagnostics.command_buffer.direct_operation_enqueue_calls ==
        disabled.diagnostics.command_buffer.direct_operation_enqueue_calls);
  CHECK(enabled.diagnostics.command_buffer.allocation_calls ==
        disabled.diagnostics.command_buffer.allocation_calls);
  CHECK(enabled.diagnostics.command_buffer.descriptor_creation_calls ==
        disabled.diagnostics.command_buffer.descriptor_creation_calls);
  CHECK(enabled.diagnostics.command_buffer.host_to_device_transfer_calls ==
        disabled.diagnostics.command_buffer.host_to_device_transfer_calls);
  CHECK(enabled.diagnostics.command_buffer.device_to_host_transfer_calls ==
        disabled.diagnostics.command_buffer.device_to_host_transfer_calls);
  CHECK(enabled.diagnostics.command_buffer.synchronization_calls ==
        disabled.diagnostics.command_buffer.synchronization_calls);
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
void checkApproxConvergenceChecks(
    const std::vector<hipdlp::algorithm::FixedTraceConvergenceCheck>& actual,
    const std::vector<hipdlp::algorithm::FixedTraceConvergenceCheck>& expected,
    const double margin) {
  REQUIRE(actual.size() == expected.size());
  for (std::size_t check = 0; check < actual.size(); ++check) {
    CAPTURE(check);
    CHECK(actual[check].iteration == expected[check].iteration);
    CHECK(actual[check].stopping.error.code ==
          expected[check].stopping.error.code);
    CHECK(actual[check].stopping.decision == expected[check].stopping.decision);
    CHECK(actual[check].stopping.metrics_valid ==
          expected[check].stopping.metrics_valid);
    const double actual_values[] = {
        actual[check].raw.primal_residual,
        actual[check].raw.dual_residual,
        actual[check].raw.objective_gap,
        actual[check].raw.primal_objective,
        actual[check].raw.dual_objective,
        actual[check].stopping.metrics.relative_primal,
        actual[check].stopping.metrics.relative_dual,
        actual[check].stopping.metrics.relative_gap};
    const double expected_values[] = {
        expected[check].raw.primal_residual,
        expected[check].raw.dual_residual,
        expected[check].raw.objective_gap,
        expected[check].raw.primal_objective,
        expected[check].raw.dual_objective,
        expected[check].stopping.metrics.relative_primal,
        expected[check].stopping.metrics.relative_dual,
        expected[check].stopping.metrics.relative_gap};
    for (std::size_t value = 0; value < 8; ++value) {
      CAPTURE(value);
      CHECK(actual_values[value] ==
            Approx(expected_values[value]).margin(margin));
    }
    CHECK(actual[check].raw.unscaled_rhs_norm ==
          expected[check].raw.unscaled_rhs_norm);
    CHECK(actual[check].raw.unscaled_cost_norm ==
          expected[check].raw.unscaled_cost_norm);
  }
}

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
  REQUIRE(actual.fixed_point_error.observed ==
          expected.fixed_point_error.observed);
  CHECK(actual.fixed_point_error.evaluation.error.code ==
        expected.fixed_point_error.evaluation.error.code);
  CHECK(actual.fixed_point_error.evaluation.diagnostics
            .device_to_host_scalar_readbacks ==
        expected.fixed_point_error.evaluation.diagnostics
            .device_to_host_scalar_readbacks);
  if (expected.fixed_point_error.observed) {
    REQUIRE(actual.fixed_point_error.evaluation.ok());
    REQUIRE(expected.fixed_point_error.evaluation.ok());
    CHECK(actual.fixed_point_error.evaluation.fixed_point_error ==
          Approx(expected.fixed_point_error.evaluation.fixed_point_error)
              .margin(margin));
  } else {
    checkNoFixedPointError(actual);
    checkNoFixedPointError(expected);
  }
}

void checkFixedPointErrorAuditDelta(
    const hipdlp::backend::cuda::testing::ResourceAudit& disabled,
    const hipdlp::backend::cuda::testing::ResourceAudit& enabled) {
  CHECK(enabled.stream_creations == disabled.stream_creations);
  CHECK(enabled.stream_destructions == disabled.stream_destructions);
  CHECK(enabled.cublas_creations == disabled.cublas_creations);
  CHECK(enabled.cublas_destructions == disabled.cublas_destructions);
  CHECK(enabled.cusparse_creations == disabled.cusparse_creations);
  CHECK(enabled.cusparse_destructions == disabled.cusparse_destructions);
  CHECK(enabled.allocation_calls == disabled.allocation_calls);
  CHECK(enabled.free_calls == disabled.free_calls);
  CHECK(enabled.descriptor_creations == disabled.descriptor_creations);
  CHECK(enabled.descriptor_destructions == disabled.descriptor_destructions);
  CHECK(enabled.host_to_device_transfers == disabled.host_to_device_transfers);
  CHECK(enabled.device_to_host_transfers ==
      disabled.device_to_host_transfers + 2);
  CHECK(enabled.device_to_device_transfers ==
        disabled.device_to_device_transfers);
  CHECK(enabled.graph_capture_calls == disabled.graph_capture_calls);
  CHECK(enabled.graph_creations == disabled.graph_creations);
  CHECK(enabled.graph_destructions == disabled.graph_destructions);
  CHECK(enabled.graph_instantiation_calls ==
        disabled.graph_instantiation_calls);
  CHECK(enabled.graph_exec_creations == disabled.graph_exec_creations);
  CHECK(enabled.graph_exec_destructions == disabled.graph_exec_destructions);
  CHECK(enabled.graph_launch_calls == disabled.graph_launch_calls);
  CHECK(enabled.graph_update_calls == disabled.graph_update_calls);
  CHECK(enabled.graph_capture_primitive_enqueue_calls ==
        disabled.graph_capture_primitive_enqueue_calls);
  CHECK(enabled.direct_primitive_enqueue_calls >
        disabled.direct_primitive_enqueue_calls);
  CHECK(enabled.synchronization_calls == disabled.synchronization_calls + 1);
}

void checkSameResourceAudit(
    const hipdlp::backend::cuda::testing::ResourceAudit& left,
    const hipdlp::backend::cuda::testing::ResourceAudit& right) {
  CHECK(left.stream_creations == right.stream_creations);
  CHECK(left.stream_destructions == right.stream_destructions);
  CHECK(left.cublas_creations == right.cublas_creations);
  CHECK(left.cublas_destructions == right.cublas_destructions);
  CHECK(left.cusparse_creations == right.cusparse_creations);
  CHECK(left.cusparse_destructions == right.cusparse_destructions);
  CHECK(left.allocation_calls == right.allocation_calls);
  CHECK(left.free_calls == right.free_calls);
  CHECK(left.descriptor_creations == right.descriptor_creations);
  CHECK(left.descriptor_destructions == right.descriptor_destructions);
  CHECK(left.host_to_device_transfers == right.host_to_device_transfers);
  CHECK(left.device_to_host_transfers == right.device_to_host_transfers);
  CHECK(left.device_to_device_transfers == right.device_to_device_transfers);
  CHECK(left.graph_capture_calls == right.graph_capture_calls);
  CHECK(left.graph_creations == right.graph_creations);
  CHECK(left.graph_destructions == right.graph_destructions);
  CHECK(left.graph_instantiation_calls == right.graph_instantiation_calls);
  CHECK(left.graph_exec_creations == right.graph_exec_creations);
  CHECK(left.graph_exec_destructions == right.graph_exec_destructions);
  CHECK(left.graph_launch_calls == right.graph_launch_calls);
  CHECK(left.graph_update_calls == right.graph_update_calls);
  CHECK(left.graph_capture_primitive_enqueue_calls ==
        right.graph_capture_primitive_enqueue_calls);
  CHECK(left.direct_primitive_enqueue_calls ==
        right.direct_primitive_enqueue_calls);
  CHECK(left.synchronization_calls == right.synchronization_calls);
}

void checkPeriodicConvergenceAuditDelta(
    const hipdlp::backend::cuda::testing::ResourceAudit& disabled,
    const hipdlp::backend::cuda::testing::ResourceAudit& enabled,
    const std::size_t evaluations) {
  CHECK(enabled.stream_creations == disabled.stream_creations);
  CHECK(enabled.stream_destructions == disabled.stream_destructions);
  CHECK(enabled.cublas_creations == disabled.cublas_creations);
  CHECK(enabled.cublas_destructions == disabled.cublas_destructions);
  CHECK(enabled.cusparse_creations == disabled.cusparse_creations);
  CHECK(enabled.cusparse_destructions == disabled.cusparse_destructions);
  CHECK(enabled.allocation_calls == disabled.allocation_calls + 2);
  CHECK(enabled.free_calls == disabled.free_calls + 2);
  CHECK(enabled.descriptor_creations == disabled.descriptor_creations + 2);
  CHECK(enabled.descriptor_destructions ==
        disabled.descriptor_destructions + 2);
  CHECK(enabled.host_to_device_transfers == disabled.host_to_device_transfers);
  CHECK(enabled.device_to_host_transfers ==
        disabled.device_to_host_transfers + 8 * evaluations);
  CHECK(enabled.device_to_device_transfers ==
        disabled.device_to_device_transfers);
  CHECK(enabled.graph_capture_calls == disabled.graph_capture_calls);
  CHECK(enabled.graph_creations == disabled.graph_creations);
  CHECK(enabled.graph_destructions == disabled.graph_destructions);
  CHECK(enabled.graph_instantiation_calls ==
        disabled.graph_instantiation_calls);
  CHECK(enabled.graph_exec_creations == disabled.graph_exec_creations);
  CHECK(enabled.graph_exec_destructions == disabled.graph_exec_destructions);
  CHECK(enabled.graph_launch_calls == disabled.graph_launch_calls);
  CHECK(enabled.graph_update_calls == disabled.graph_update_calls);
  CHECK(enabled.graph_capture_primitive_enqueue_calls ==
        disabled.graph_capture_primitive_enqueue_calls);
  CHECK(enabled.direct_primitive_enqueue_calls ==
        disabled.direct_primitive_enqueue_calls + 4 * evaluations);
  CHECK(enabled.synchronization_calls ==
        disabled.synchronization_calls + 4 * evaluations);
}

void checkFixedPointCheckAuditDelta(
    const hipdlp::backend::cuda::testing::ResourceAudit& base,
    const hipdlp::backend::cuda::testing::ResourceAudit& extended,
    const std::size_t extra_checks) {
  CHECK(extended.stream_creations == base.stream_creations);
  CHECK(extended.stream_destructions == base.stream_destructions);
  CHECK(extended.cublas_creations == base.cublas_creations);
  CHECK(extended.cublas_destructions == base.cublas_destructions);
  CHECK(extended.cusparse_creations == base.cusparse_creations);
  CHECK(extended.cusparse_destructions == base.cusparse_destructions);
  CHECK(extended.allocation_calls == base.allocation_calls);
  CHECK(extended.free_calls == base.free_calls);
  CHECK(extended.descriptor_creations == base.descriptor_creations);
  CHECK(extended.descriptor_destructions == base.descriptor_destructions);
  CHECK(extended.host_to_device_transfers == base.host_to_device_transfers);
  CHECK(extended.device_to_host_transfers ==
      base.device_to_host_transfers + 2 * extra_checks);
  CHECK(extended.graph_capture_calls == base.graph_capture_calls);
  CHECK(extended.graph_creations == base.graph_creations);
  CHECK(extended.graph_destructions == base.graph_destructions);
  CHECK(extended.graph_instantiation_calls == base.graph_instantiation_calls);
  CHECK(extended.graph_exec_creations == base.graph_exec_creations);
  CHECK(extended.graph_exec_destructions == base.graph_exec_destructions);
  CHECK(extended.graph_launch_calls == base.graph_launch_calls);
  CHECK(extended.graph_update_calls == base.graph_update_calls);
  CHECK(extended.graph_capture_primitive_enqueue_calls ==
        base.graph_capture_primitive_enqueue_calls);
  CHECK(extended.synchronization_calls ==
      base.synchronization_calls + extra_checks);
  if (extra_checks == 0) {
    CHECK(extended.device_to_device_transfers ==
          base.device_to_device_transfers);
    CHECK(extended.direct_primitive_enqueue_calls ==
          base.direct_primitive_enqueue_calls);
  } else {
    CHECK(extended.device_to_device_transfers ==
          base.device_to_device_transfers);
    CHECK(extended.direct_primitive_enqueue_calls >
          base.direct_primitive_enqueue_calls);
  }
}

void checkApproxFixedPointRestart(
    const hipdlp::algorithm::FixedTraceFixedPointRestartObservation& actual,
    const hipdlp::algorithm::FixedTraceFixedPointRestartObservation& expected,
    const double margin) {
  REQUIRE(actual.observed == expected.observed);
  CHECK(actual.evaluation.error.code == expected.evaluation.error.code);
  if (!expected.observed) {
    checkNoFixedPointRestartObservation(actual);
    return;
  }
  REQUIRE(actual.evaluation.ok());
  CHECK(actual.evaluation.decision == expected.evaluation.decision);
  CHECK(actual.evaluation.condition == expected.evaluation.condition);
  CHECK(actual.evaluation.next.history_valid ==
        expected.evaluation.next.history_valid);
  CHECK(actual.evaluation.next.reference_valid ==
        expected.evaluation.next.reference_valid);
  CHECK(actual.evaluation.next.previous_trial_valid ==
        expected.evaluation.next.previous_trial_valid);
  CHECK(actual.evaluation.next.reference_error ==
        expected.evaluation.next.reference_error);
  CHECK(actual.evaluation.next.previous_trial_error ==
        Approx(expected.evaluation.next.previous_trial_error).margin(margin));
}

// An executed restart adds only device-to-device copies and direct kernel
// enqueues. A restart that clears the average also skips its four downloads,
// each of which costs two device-to-host events and one synchronization.
void checkRestartTransitionAuditDelta(
    const hipdlp::backend::cuda::testing::ResourceAudit& base,
    const hipdlp::backend::cuda::testing::ResourceAudit& restarted,
    const std::size_t skipped_downloads) {
  CHECK(restarted.stream_creations == base.stream_creations);
  CHECK(restarted.stream_destructions == base.stream_destructions);
  CHECK(restarted.cublas_creations == base.cublas_creations);
  CHECK(restarted.cublas_destructions == base.cublas_destructions);
  CHECK(restarted.cusparse_creations == base.cusparse_creations);
  CHECK(restarted.cusparse_destructions == base.cusparse_destructions);
  CHECK(restarted.allocation_calls == base.allocation_calls);
  CHECK(restarted.free_calls == base.free_calls);
  CHECK(restarted.descriptor_creations == base.descriptor_creations);
  CHECK(restarted.descriptor_destructions == base.descriptor_destructions);
  CHECK(restarted.host_to_device_transfers == base.host_to_device_transfers);
  CHECK(restarted.device_to_host_transfers + 2 * skipped_downloads ==
        base.device_to_host_transfers);
  CHECK(restarted.synchronization_calls + skipped_downloads ==
        base.synchronization_calls);
  CHECK(restarted.graph_capture_calls == base.graph_capture_calls);
  CHECK(restarted.graph_creations == base.graph_creations);
  CHECK(restarted.graph_destructions == base.graph_destructions);
  CHECK(restarted.graph_instantiation_calls == base.graph_instantiation_calls);
  CHECK(restarted.graph_exec_creations == base.graph_exec_creations);
  CHECK(restarted.graph_exec_destructions == base.graph_exec_destructions);
  CHECK(restarted.graph_launch_calls == base.graph_launch_calls);
  CHECK(restarted.graph_update_calls == base.graph_update_calls);
  CHECK(restarted.graph_capture_primitive_enqueue_calls ==
        base.graph_capture_primitive_enqueue_calls);
  CHECK(restarted.device_to_device_transfers > base.device_to_device_transfers);
  CHECK(restarted.direct_primitive_enqueue_calls >
        base.direct_primitive_enqueue_calls);
}

void checkApproxPrimalWeightUpdates(
    const std::vector<hipdlp::algorithm::FixedTracePrimalWeightUpdate>& actual,
    const std::vector<hipdlp::algorithm::FixedTracePrimalWeightUpdate>&
        expected,
    const double margin) {
  REQUIRE(actual.size() == expected.size());
  for (std::size_t update = 0; update < actual.size(); ++update) {
    CAPTURE(update);
    CHECK(actual[update].iteration == expected[update].iteration);
    CHECK(actual[update].result.error.code ==
          expected[update].result.error.code);
    CHECK(actual[update].result.kind == expected[update].result.kind);
    const double actual_values[] = {
        actual[update].primal_anchor_distance,
        actual[update].dual_anchor_distance,
        actual[update].result.next.primal_weight,
        actual[update].result.next.error_sum,
        actual[update].result.step_sizes.base_step,
        actual[update].result.step_sizes.primal_step,
        actual[update].result.step_sizes.dual_step};
    const double expected_values[] = {
        expected[update].primal_anchor_distance,
        expected[update].dual_anchor_distance,
        expected[update].result.next.primal_weight,
        expected[update].result.next.error_sum,
        expected[update].result.step_sizes.base_step,
        expected[update].result.step_sizes.primal_step,
        expected[update].result.step_sizes.dual_step};
    for (std::size_t value = 0; value < 7; ++value) {
      CAPTURE(value);
      CHECK(actual_values[value] ==
            Approx(expected_values[value]).margin(margin));
    }
  }
}

void checkPrimalWeightAuditDelta(
    const hipdlp::backend::cuda::testing::ResourceAudit& disabled,
    const hipdlp::backend::cuda::testing::ResourceAudit& enabled,
    const std::size_t updates) {
  CHECK(enabled.stream_creations == disabled.stream_creations);
  CHECK(enabled.stream_destructions == disabled.stream_destructions);
  CHECK(enabled.cublas_creations == disabled.cublas_creations);
  CHECK(enabled.cublas_destructions == disabled.cublas_destructions);
  CHECK(enabled.cusparse_creations == disabled.cusparse_creations);
  CHECK(enabled.cusparse_destructions == disabled.cusparse_destructions);
  CHECK(enabled.allocation_calls == disabled.allocation_calls);
  CHECK(enabled.free_calls == disabled.free_calls);
  CHECK(enabled.descriptor_creations == disabled.descriptor_creations);
  CHECK(enabled.descriptor_destructions == disabled.descriptor_destructions);
  CHECK(enabled.host_to_device_transfers == disabled.host_to_device_transfers);
  CHECK(enabled.device_to_host_transfers ==
        disabled.device_to_host_transfers + 4 * updates);
  CHECK(enabled.device_to_device_transfers ==
        disabled.device_to_device_transfers);
  CHECK(enabled.graph_capture_calls == disabled.graph_capture_calls);
  CHECK(enabled.graph_creations == disabled.graph_creations);
  CHECK(enabled.graph_destructions == disabled.graph_destructions);
  CHECK(enabled.graph_instantiation_calls ==
        disabled.graph_instantiation_calls);
  CHECK(enabled.graph_exec_creations == disabled.graph_exec_creations);
  CHECK(enabled.graph_exec_destructions == disabled.graph_exec_destructions);
  CHECK(enabled.graph_launch_calls == disabled.graph_launch_calls);
  CHECK(enabled.graph_update_calls == disabled.graph_update_calls);
  CHECK(enabled.graph_capture_primitive_enqueue_calls ==
        disabled.graph_capture_primitive_enqueue_calls);
  CHECK(enabled.direct_primitive_enqueue_calls ==
      disabled.direct_primitive_enqueue_calls + 6 * updates);
  CHECK(enabled.synchronization_calls ==
        disabled.synchronization_calls + 2 * updates);
}
#endif

// Exact restart-at-two fixture: makeOneStepLp with tau = 1/2, sigma = 1/4,
// omega = eta = 1, a zero start, check interval 2, and cadence 2.
const double kCandidatePrimalAtTwo = 9.0 / 8.0;
const double kCandidateDualAtTwo = 5.0 / 16.0;
const double kErrorAtTwo = 11.0 / 16.0;
const double kRestartedPrimalAtThree = 57.0 / 32.0;
const double kRestartedDualAtThree = 13.0 / 64.0;
const double kRestartedErrorAtThree = 35.0 / 64.0;
const double kRestartedErrorAtFour = 21.0 / 128.0;

hipdlp::algorithm::FixedTraceOptions makeRestartAtTwoOptions(
    const HighsInt iteration_count) {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(iteration_count);
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestartCadence(options, makeRestartReference(10.0), 2, 2);
  return options;
}

struct RestartExecutionSelection {
  hipdlp::backend::Kind backend;
  hipdlp::algorithm::FixedTraceExecutionMode mode;
};

std::vector<RestartExecutionSelection> restartExecutionSelections() {
  std::vector<RestartExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif
  return selections;
}

// Both anchors still hold the zero start when the restart at two is decided,
// so the candidate-to-anchor distances are the candidate magnitudes.
const double kAnchorPrimalDistanceAtTwo = 9.0 / 8.0;
const double kAnchorDualDistanceAtTwo = 5.0 / 16.0;
const double kRelativePrimalAtTwo = 7.0 / 24.0;
const double kBoundedRelativeDualAtTwo = 0.0;
const double kFreeRelativeDualAtTwo = 21.0 / 32.0;

hipdlp::algorithm::FixedTraceOptions makePrimalWeightOptions(
    const HighsInt iteration_count) {
  hipdlp::algorithm::FixedTraceOptions options =
      makeRestartAtTwoOptions(iteration_count);
  enablePeriodicConvergence(options, 2);
  enablePrimalWeight(options, 1.0);
  return options;
}

// Independent oracle: the Tier A12 policy called directly on the scalars the
// shared boundary is required to forward.
hipdlp::algorithm::PrimalWeightUpdateResult expectedRestartPrimalWeightAtTwo(
    const double initial_weight, const double relative_dual) {
  hipdlp::algorithm::PrimalWeightUpdateInput input;
  input.previous = hipdlp::algorithm::resetPrimalWeight(initial_weight);
  input.primal_anchor_distance = kAnchorPrimalDistanceAtTwo;
  input.dual_anchor_distance = kAnchorDualDistanceAtTwo;
  input.relative_primal_residual = kRelativePrimalAtTwo;
  input.relative_dual_residual = relative_dual;
  input.primal_step_size = 0.5;
  input.dual_step_size = 0.25;
  return hipdlp::algorithm::updatePrimalWeightAtRestart(input);
}

void checkPrimalWeightUpdatesFollowRestarts(
    const hipdlp::algorithm::FixedTraceResult& result) {
  std::size_t update = 0;
  for (std::size_t check = 0; check < result.fixed_point_restart_checks.size();
       ++check) {
    if (result.fixed_point_restart_checks[check].evaluation.decision !=
        hipdlp::algorithm::RestartDecision::kRestart)
      continue;
    CAPTURE(check);
    REQUIRE(update < result.primal_weight_updates.size());
    CHECK(result.primal_weight_updates[update].iteration ==
          result.fixed_point_restart_checks[check].iteration);
    CHECK(result.primal_weight_updates[update].result.ok());
    ++update;
  }
  CHECK(update == result.primal_weight_updates.size());
}

hipdlp::algorithm::FixedTraceResult runTrace(
    const HighsLp& lp, const hipdlp::algorithm::FixedTraceOptions& options) {
  return hipdlp::algorithm::runFixedTrace(hipdlp::model::OriginalModel(lp),
                                          options);
}

// Independent oracle: the fixed-point error a separate final-only trace of the
// same prefix length reports.
double prefixFixedPointError(const HighsInt iteration_count) {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(iteration_count);
  enableFixedPointError(options, 1.0, 1.0);
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.fixed_point_error.observed);
  REQUIRE(result.fixed_point_error.evaluation.ok());
  return result.fixed_point_error.evaluation.fixed_point_error;
}

}  // namespace

TEST_CASE("hi-pdlp-fixed-trace-observes-the-exact-final-fixed-point-error",
          "[pdlp][hipdlp][phase7i][fixed-point-error][integration][exact]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
  enableFixedPointError(options, 1.0, 1.0);

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 2);
  checkExactPoint(result.points[1], {0.5}, {0.25}, {0.5}, {0.25});

  const double delta_primal = 0.5 - 1.0;
  const double delta_dual = 0.25 - 0.5;
  REQUIRE(delta_primal == -0.5);
  REQUIRE(delta_dual == -0.25);
  const double expected =
      std::sqrt(delta_primal * delta_primal + delta_dual * delta_dual +
                2.0 * delta_primal * delta_dual);
  REQUIRE(expected == 0.75);

  const hipdlp::algorithm::FixedTraceFixedPointErrorObservation& observation =
      result.fixed_point_error;
  INFO(observation.evaluation.error.message);
  REQUIRE(observation.observed);
  REQUIRE(observation.evaluation.ok());
  CHECK(observation.evaluation.fixed_point_error == expected);
  CHECK(observation.evaluation.diagnostics.device_to_host_scalar_readbacks ==
        3);
}

TEST_CASE("hi-pdlp-fixed-trace-forwards-fixed-point-error-scalars",
          "[pdlp][hipdlp][phase7i][fixed-point-error][integration]"
          "[parameters]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
  enableFixedPointError(options, 4.0, 1.0);

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());

  const double delta_primal = 0.5 - 1.0;
  const double delta_dual = 0.25 - 0.5;
  const double expected = std::sqrt(4.0 * delta_primal * delta_primal +
                                    delta_dual * delta_dual / 4.0 +
                                    2.0 * delta_primal * delta_dual);
  REQUIRE(expected == 1.125);
  REQUIRE(result.fixed_point_error.observed);
  REQUIRE(result.fixed_point_error.evaluation.ok());
  CHECK(result.fixed_point_error.evaluation.fixed_point_error == expected);
  CHECK(result.fixed_point_error.evaluation.diagnostics
            .device_to_host_scalar_readbacks == 3);
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-error-is-opt-in-and-zero-safe",
          "[pdlp][hipdlp][phase7i][fixed-point-error][integration][disabled]"
          "[zero]") {
  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
    bool supported;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
                        true});
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded,
                        false});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
                        true});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded,
                        true});
#endif

  const hipdlp::algorithm::FixedTraceOptions defaults;
  CHECK_FALSE(defaults.fixed_point_error.enabled);
  CHECK(std::isnan(defaults.fixed_point_error.omega));
  CHECK(std::isnan(defaults.fixed_point_error.eta));

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(1);
    disabled_options.backend = selection.backend;
    disabled_options.execution_mode = selection.mode;
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    checkNoFixedPointError(disabled);
    if (!selection.supported) {
      CHECK_FALSE(disabled.ok());
      CHECK(disabled.error.code ==
            hipdlp::api::ErrorCode::kUnsupportedCapability);
      continue;
    }
    INFO(disabled.error.message);
    REQUIRE(disabled.ok());

    hipdlp::algorithm::FixedTraceOptions zero_options = makeOptions(0);
    zero_options.backend = selection.backend;
    zero_options.execution_mode = selection.mode;
    enableFixedPointError(zero_options, 4.0, 1.0);
    const hipdlp::algorithm::FixedTraceResult zero =
        runTrace(makeOneStepLp(), zero_options);
    INFO(zero.error.message);
    REQUIRE(zero.ok());
    REQUIRE(zero.points.size() == 1);
    checkNoFixedPointError(zero);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-error-failures-are-transactional",
          "[pdlp][hipdlp][phase7i][fixed-point-error][integration][errors]"
          "[transactional]") {
  SECTION("enabled scalars are validated for a zero-iteration trace") {
    const double invalid_scalars[] = {0.0, -1.0,
                                      std::numeric_limits<double>::quiet_NaN(),
                                      std::numeric_limits<double>::infinity()};
    for (const double invalid : invalid_scalars) {
      CAPTURE(invalid);
      hipdlp::algorithm::FixedTraceOptions omega_options = makeOptions(0);
      enableFixedPointError(omega_options, invalid, 1.0);
      const hipdlp::algorithm::FixedTraceResult omega_result =
          runTrace(makeOneStepLp(), omega_options);
      checkTransactionalFailure(omega_result);
      CHECK(omega_result.error.code ==
            hipdlp::api::ErrorCode::kInvalidArgument);
      CHECK(
          omega_result.diagnostics.resource_creation_calls_before_iterations ==
          0);

      hipdlp::algorithm::FixedTraceOptions eta_options = makeOptions(0);
      enableFixedPointError(eta_options, 1.0, invalid);
      const hipdlp::algorithm::FixedTraceResult eta_result =
          runTrace(makeOneStepLp(), eta_options);
      checkTransactionalFailure(eta_result);
      CHECK(eta_result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
      CHECK(eta_result.diagnostics.resource_creation_calls_before_iterations ==
            0);
    }
  }

  SECTION("weighted movement overflow discards the completed trace") {
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(1);
    disabled_options.primal_step_size = 4.0;
    disabled_options.dual_step_size = 0.25;
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    INFO(disabled.error.message);
    REQUIRE(disabled.ok());

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enableFixedPointError(enabled_options, std::numeric_limits<double>::max(),
                          1.0);
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    checkTransactionalFailure(enabled);
    CHECK(enabled.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(enabled.error.message ==
          "HiPDLP fixed-point weighted movement is outside FP64 range");
    CHECK(enabled.diagnostics.transfer_calls_after_iterations == 0);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-error-is-bitwise-deterministic",
          "[pdlp][hipdlp][phase7i][fixed-point-error][integration]"
          "[determinism]") {
  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(4);
    disabled_options.backend = selection.backend;
    disabled_options.execution_mode = selection.mode;
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enableFixedPointError(enabled_options, 4.0, 1.0);
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), enabled_options);
    INFO(disabled.error.message);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    checkBitwiseEqual(first, second);
    checkBitwiseEqualTraceOutputs(disabled, first);
    checkSameDiagnostics(disabled, first);
  }
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-fixed-trace-fixed-point-error-cpu-cuda-resource-parity",
          "[pdlp][hipdlp][phase7i][fixed-point-error][integration][cuda]"
          "[resources][parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  hipdlp::algorithm::FixedTraceOptions cpu_options = makeOptions(1);
  enableFixedPointError(cpu_options, 4.0, 1.0);
  const hipdlp::algorithm::FixedTraceResult cpu =
      runTrace(makeOneStepLp(), cpu_options);
  INFO(cpu.error.message);
  REQUIRE(cpu.ok());

  const hipdlp::algorithm::FixedTraceExecutionMode modes[] = {
      hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded};
  for (const hipdlp::algorithm::FixedTraceExecutionMode mode : modes) {
    CAPTURE(static_cast<int>(mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(1);
    disabled_options.backend = hipdlp::backend::Kind::kCuda;
    disabled_options.execution_mode = mode;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enableFixedPointError(enabled_options, 4.0, 1.0);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    checkApproxTrace(enabled, cpu, 5e-13);
    checkBitwiseEqualTraceOutputs(disabled, enabled);
    checkSameDiagnostics(disabled, enabled);
    checkFixedPointErrorAuditDelta(disabled_audit, enabled_audit);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-error-preserves-graph-launches",
          "[pdlp][hipdlp][phase7i][fixed-point-error][integration][cuda]"
          "[graph][replay]") {
  const std::size_t expected_launches[] = {0, 1, 1, 2, 3};
  const HighsInt iteration_counts[] = {0, 1, 2, 4, 5};
  for (std::size_t selection = 0; selection < 5; ++selection) {
    const HighsInt iteration_count = iteration_counts[selection];
    CAPTURE(iteration_count);
    hipdlp::algorithm::FixedTraceOptions disabled_options =
        makeOptions(iteration_count);
    disabled_options.backend = hipdlp::backend::Kind::kCuda;
    disabled_options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enableFixedPointError(enabled_options, 1.0, 1.0);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    CHECK(disabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(enabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(disabled_audit.graph_launch_calls == expected_launches[selection]);
    CHECK(enabled_audit.graph_launch_calls == expected_launches[selection]);
    if (iteration_count == 0) {
      checkNoFixedPointError(enabled);
      CHECK(enabled_audit.device_to_host_transfers ==
            disabled_audit.device_to_host_transfers);
      CHECK(enabled_audit.synchronization_calls ==
            disabled_audit.synchronization_calls);
    } else {
      REQUIRE(enabled.fixed_point_error.observed);
      REQUIRE(enabled.fixed_point_error.evaluation.ok());
      CHECK(enabled.fixed_point_error.evaluation.diagnostics
                .device_to_host_scalar_readbacks == 3);
      checkFixedPointErrorAuditDelta(disabled_audit, enabled_audit);
    }
  }
}
#endif

TEST_CASE("hi-pdlp-fixed-trace-observes-the-final-fixed-point-restart-trigger",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration][exact]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestart(options, makeRestartReference(10.0), 1);

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 2);
  checkExactPoint(result.points[1], {0.5}, {0.25}, {0.5}, {0.25});
  // The executed restart leaves no post-restart sample to average.
  checkNoAverage(result);
  checkNoObservation(result.averaged_convergence);
  REQUIRE(result.current_convergence.observed);
  REQUIRE(result.fixed_point_error.observed);
  REQUIRE(result.fixed_point_error.evaluation.ok());
  CHECK(result.fixed_point_error.evaluation.fixed_point_error == 0.75);

  const hipdlp::algorithm::FixedTraceFixedPointRestartObservation& observation =
      result.fixed_point_restart;
  INFO(observation.evaluation.error.message);
  REQUIRE(observation.observed);
  REQUIRE(observation.evaluation.ok());
  CHECK(observation.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(observation.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kFirstCheckInterval);
  // The reported next state is data; the caller-supplied history is not
  // refreshed in place.
  CHECK(observation.evaluation.next.history_valid);
  CHECK_FALSE(observation.evaluation.next.reference_valid);
  CHECK_FALSE(observation.evaluation.next.previous_trial_valid);
  CHECK(observation.evaluation.next.reference_error == 0.0);
  CHECK(observation.evaluation.next.previous_trial_error == 0.0);
  CHECK(options.fixed_point_restart.previous.reference_valid);
  CHECK(options.fixed_point_restart.previous.reference_error == 10.0);
  CHECK_FALSE(options.fixed_point_restart.previous.previous_trial_valid);
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-and-optimality-coexist",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration]"
          "[convergence]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
  options.scaled_primal_start = {2.0};
  options.scaled_dual_start = {-1.0};
  options.primal_step_size = 1e-12;
  options.dual_step_size = 1e-12;
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestart(options, makeRestartReference(10.0), 1);

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.current_convergence.observed);
  CHECK(result.current_convergence.stopping.decision ==
        hipdlp::algorithm::StoppingDecision::kOptimal);
  // The executed restart discards the only average sample.
  checkNoAverage(result);
  checkNoObservation(result.averaged_convergence);
  REQUIRE(result.fixed_point_restart.observed);
  REQUIRE(result.fixed_point_restart.evaluation.ok());
  CHECK(result.fixed_point_restart.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(result.fixed_point_restart.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kFirstCheckInterval);
}

TEST_CASE("hi-pdlp-fixed-trace-forwards-fixed-point-restart-scalars",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration]"
          "[parameters]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestart(options, makeRestartReference(10.0), 2);

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.fixed_point_error.observed);
  REQUIRE(result.fixed_point_error.evaluation.fixed_point_error == 0.75);

  const hipdlp::algorithm::FixedTraceFixedPointRestartObservation& observation =
      result.fixed_point_restart;
  INFO(observation.evaluation.error.message);
  REQUIRE(observation.observed);
  REQUIRE(observation.evaluation.ok());
  CHECK(observation.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(observation.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
  CHECK(observation.evaluation.next.history_valid);
  CHECK(observation.evaluation.next.reference_valid);
  CHECK(observation.evaluation.next.reference_error == 10.0);
  CHECK(observation.evaluation.next.previous_trial_valid);
  CHECK(observation.evaluation.next.previous_trial_error == 0.75);
  CHECK(options.fixed_point_restart.previous.reference_error == 10.0);
  CHECK_FALSE(options.fixed_point_restart.previous.previous_trial_valid);
  CHECK(options.fixed_point_restart.previous.previous_trial_error == 0.0);
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-conditions-fire-in-order",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration]"
          "[condition]") {
  // The second projected candidate is 9/8 primal and 5/16 dual, while the
  // reflected point is 7/4 primal and 3/8 dual. The legacy fixed-point error
  // is therefore sqrt(121/256) = 11/16.
  const double anchor_weight = 1.0 / 3.0;
  const double candidate_primal = 9.0 / 8.0;
  const double candidate_dual = 5.0 / 16.0;
  const double reflected_primal = 1.75;
  const double reflected_dual = 0.375;
  const double delta_primal = candidate_primal - reflected_primal;
  const double delta_dual = candidate_dual - reflected_dual;
  REQUIRE(delta_primal == -5.0 / 8.0);
  REQUIRE(delta_dual == -1.0 / 16.0);
  const double expected_error =
      std::sqrt(delta_primal * delta_primal + delta_dual * delta_dual +
                2.0 * delta_primal * delta_dual);
  REQUIRE(expected_error == 11.0 / 16.0);

  struct Expectation {
    hipdlp::algorithm::FixedPointRestartState previous;
    double artificial_threshold;
    hipdlp::algorithm::RestartDecision decision;
    hipdlp::algorithm::FixedPointRestartCondition condition;
  };
  std::vector<Expectation> expectations;
  // 11/16 <= 0.2 * 10.
  expectations.push_back(
      {makeRestartHistory(10.0, 0.5), 0.36,
       hipdlp::algorithm::RestartDecision::kRestart,
       hipdlp::algorithm::FixedPointRestartCondition::kSufficientDecay});
  // 0.2 < 11/16 <= 0.8 and the trial error worsens.
  expectations.push_back(
      {makeRestartHistory(1.0, 0.5), 0.36,
       hipdlp::algorithm::RestartDecision::kRestart,
       hipdlp::algorithm::FixedPointRestartCondition::kNecessaryDecay});
  // The necessary-decay window is entered but non-improvement fails, so the
  // fall-through still reaches the artificial check with 2 >= 1.0 * 2.
  expectations.push_back(
      {makeRestartHistory(1.0, 1.0), 1.0,
       hipdlp::algorithm::RestartDecision::kRestart,
       hipdlp::algorithm::FixedPointRestartCondition::kArtificialThreshold});
  // No decay window is entered and 2 >= 2.0 * 2 is false.
  expectations.push_back(
      {makeRestartHistory(0.8, 0.5), 2.0,
       hipdlp::algorithm::RestartDecision::kContinue,
       hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered});

  for (std::size_t selection = 0; selection < expectations.size();
       ++selection) {
    CAPTURE(selection);
    const Expectation& expectation = expectations[selection];
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    enableFixedPointError(options, 1.0, 1.0);
    enableFixedPointRestart(options, expectation.previous, 1);
    options.fixed_point_restart.parameters.artificial_threshold =
        expectation.artificial_threshold;

    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.points.size() == 3);
    if (expectation.decision == hipdlp::algorithm::RestartDecision::kRestart) {
      CHECK(result.points[2].primal == std::vector<double>({candidate_primal}));
      CHECK(result.points[2].dual == std::vector<double>({candidate_dual}));
    } else {
      CHECK(result.points[2].primal ==
            std::vector<double>({(1.0 - anchor_weight) * reflected_primal}));
      CHECK(result.points[2].dual == std::vector<double>({0.25}));
    }
    REQUIRE(result.fixed_point_error.observed);
    REQUIRE(result.fixed_point_error.evaluation.ok());
    CHECK(result.fixed_point_error.evaluation.fixed_point_error ==
          expected_error);

    const hipdlp::algorithm::FixedTraceFixedPointRestartObservation&
        observation = result.fixed_point_restart;
    INFO(observation.evaluation.error.message);
    REQUIRE(observation.observed);
    REQUIRE(observation.evaluation.ok());
    CHECK(observation.evaluation.decision == expectation.decision);
    CHECK(observation.evaluation.condition == expectation.condition);
    if (expectation.decision == hipdlp::algorithm::RestartDecision::kRestart) {
      CHECK_FALSE(observation.evaluation.next.reference_valid);
      CHECK_FALSE(observation.evaluation.next.previous_trial_valid);
    } else {
      CHECK(observation.evaluation.next.reference_error ==
            expectation.previous.reference_error);
      CHECK(observation.evaluation.next.previous_trial_error == expected_error);
    }
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-is-opt-in-and-zero-safe",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration]"
          "[disabled][zero]") {
  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif

  const hipdlp::algorithm::FixedTraceOptions defaults;
  CHECK_FALSE(defaults.fixed_point_restart.enabled);
  CHECK(defaults.fixed_point_restart.parameters.sufficient_decay == 0.2);
  CHECK(defaults.fixed_point_restart.parameters.necessary_decay == 0.8);
  CHECK(defaults.fixed_point_restart.parameters.artificial_threshold == 0.36);
  CHECK(defaults.fixed_point_restart.parameters.check_interval == 40);
  CHECK(defaults.fixed_point_restart.previous.history_valid);
  CHECK_FALSE(defaults.fixed_point_restart.previous.reference_valid);
  CHECK_FALSE(defaults.fixed_point_restart.previous.previous_trial_valid);

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(1);
    disabled_options.backend = selection.backend;
    disabled_options.execution_mode = selection.mode;
    enableFixedPointError(disabled_options, 1.0, 1.0);
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    INFO(disabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(disabled.fixed_point_error.observed);
    checkNoFixedPointRestart(disabled);

    hipdlp::algorithm::FixedTraceOptions inconsistent_options = makeOptions(1);
    inconsistent_options.backend = selection.backend;
    inconsistent_options.execution_mode = selection.mode;
    enableFixedPointRestart(inconsistent_options, makeRestartReference(10.0),
                            1);
    const hipdlp::algorithm::FixedTraceResult inconsistent =
        runTrace(makeOneStepLp(), inconsistent_options);
    checkTransactionalFailure(inconsistent);
    CHECK(inconsistent.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(inconsistent.diagnostics.resource_creation_calls_before_iterations ==
          0);

    hipdlp::algorithm::FixedTraceOptions zero_options = makeOptions(0);
    zero_options.backend = selection.backend;
    zero_options.execution_mode = selection.mode;
    enableFixedPointError(zero_options, 1.0, 1.0);
    enableFixedPointRestart(zero_options, makeRestartReference(10.0), 1);
    const hipdlp::algorithm::FixedTraceResult zero =
        runTrace(makeOneStepLp(), zero_options);
    INFO(zero.error.message);
    REQUIRE(zero.ok());
    REQUIRE(zero.points.size() == 1);
    checkNoFixedPointError(zero);
    checkNoFixedPointRestart(zero);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-failures-are-transactional",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration][errors]"
          "[transactional]") {
  SECTION("an invalid previous state preserves the trigger error") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    enableFixedPointError(options, 1.0, 1.0);
    options.fixed_point_restart.enabled = true;
    REQUIRE_FALSE(options.fixed_point_restart.previous.reference_valid);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.error.message ==
          "HiPDLP fixed-point restart history is invalid or has no reference "
          "error");
  }

  SECTION("an invalid trial state preserves the trigger error") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    enableFixedPointError(options, 1.0, 1.0);
    hipdlp::algorithm::FixedPointRestartState previous =
        makeRestartReference(10.0);
    previous.previous_trial_error = 1.0;
    enableFixedPointRestart(options, previous, 1);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
  }

  SECTION("invalid parameters preserve the trigger error") {
    const HighsInt invalid_intervals[] = {0, -1};
    for (const HighsInt interval : invalid_intervals) {
      CAPTURE(interval);
      hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
      enableFixedPointError(options, 1.0, 1.0);
      enableFixedPointRestart(options, makeRestartReference(10.0), interval);
      const hipdlp::algorithm::FixedTraceResult result =
          runTrace(makeOneStepLp(), options);
      checkTransactionalFailure(result);
      CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
      CHECK(result.error.message ==
            "HiPDLP fixed-point restart parameters are invalid");
    }

    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    enableFixedPointError(options, 1.0, 1.0);
    enableFixedPointRestart(options, makeRestartReference(10.0), 1);
    options.fixed_point_restart.parameters.artificial_threshold =
        std::numeric_limits<double>::quiet_NaN();
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
  }

  SECTION("an overflowing artificial product preserves the trigger error") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    enableFixedPointError(options, 1.0, 1.0);
    enableFixedPointRestart(options, makeRestartReference(10.0), 1);
    options.fixed_point_restart.parameters.artificial_threshold =
        std::numeric_limits<double>::max();
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.error.message ==
          "HiPDLP fixed-point restart thresholds are outside FP64 range");
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-is-bitwise-deterministic",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration]"
          "[determinism]") {
  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(4);
    disabled_options.backend = selection.backend;
    disabled_options.execution_mode = selection.mode;
    enableFixedPointError(disabled_options, 4.0, 1.0);
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enableFixedPointRestart(enabled_options, makeRestartReference(10.0), 40);
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), enabled_options);
    INFO(disabled.error.message);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    checkBitwiseEqual(first, second);
    checkBitwiseEqualTraceOutputs(disabled, first);
    checkBitwiseEqualFixedPointError(disabled.fixed_point_error,
                                     first.fixed_point_error);
    checkSameDiagnostics(disabled, first);
    checkNoFixedPointRestart(disabled);
    REQUIRE(first.fixed_point_restart.observed);
    CHECK(first.fixed_point_restart.evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kContinue);
    CHECK(first.fixed_point_restart.evaluation.next.previous_trial_error ==
          first.fixed_point_error.evaluation.fixed_point_error);
  }
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-cpu-cuda-resource-parity",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration][cuda]"
          "[resources][parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  hipdlp::algorithm::FixedTraceOptions cpu_options = makeOptions(2);
  enableFixedPointError(cpu_options, 1.0, 1.0);
  enableFixedPointRestart(cpu_options, makeRestartHistory(1.0, 0.5), 1);
  const hipdlp::algorithm::FixedTraceResult cpu =
      runTrace(makeOneStepLp(), cpu_options);
  INFO(cpu.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(cpu.fixed_point_restart.observed);
  CHECK(cpu.fixed_point_restart.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNecessaryDecay);

  const hipdlp::algorithm::FixedTraceExecutionMode modes[] = {
      hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded};
  for (const hipdlp::algorithm::FixedTraceExecutionMode mode : modes) {
    CAPTURE(static_cast<int>(mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options = cpu_options;
    disabled_options.backend = hipdlp::backend::Kind::kCuda;
    disabled_options.execution_mode = mode;
    disabled_options.fixed_point_restart =
        hipdlp::algorithm::FixedTraceFixedPointRestartOptions();
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enableFixedPointRestart(enabled_options, makeRestartHistory(1.0, 0.5), 1);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    checkApproxTrace(enabled, cpu, 5e-13);
    checkApproxFixedPointRestart(enabled.fixed_point_restart,
                                 cpu.fixed_point_restart, 5e-13);
    checkBitwiseEqualFixedPointError(disabled.fixed_point_error,
                                     enabled.fixed_point_error);
    checkNoFixedPointRestart(disabled);
    // Only the executed restart separates the runs: the boundary point becomes
    // the candidate and the cleared average drops its four downloads.
    CHECK(disabled.points[2].primal != enabled.points[2].primal);
    REQUIRE(disabled.averaging.hasAverage());
    checkNoAverage(enabled);
    CHECK(enabled.diagnostics.resource_creation_calls_before_iterations ==
          disabled.diagnostics.resource_creation_calls_before_iterations);
    CHECK(enabled.diagnostics.resource_creation_calls_during_iterations == 0);
    CHECK(enabled.diagnostics.transfer_calls_before_iterations ==
          disabled.diagnostics.transfer_calls_before_iterations);
    CHECK(enabled.diagnostics.transfer_calls_during_iterations == 0);
    CHECK(enabled.diagnostics.transfer_calls_after_iterations + 4 ==
          disabled.diagnostics.transfer_calls_after_iterations);
    CHECK(enabled.diagnostics.command_buffer.replay_calls ==
          disabled.diagnostics.command_buffer.replay_calls);
    checkRestartTransitionAuditDelta(disabled_audit, enabled_audit, 4);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-preserves-graph-launches",
          "[pdlp][hipdlp][phase7j][restart][fixed-point][integration][cuda]"
          "[graph][replay]") {
  const std::size_t expected_launches[] = {0, 1, 1, 2, 3};
  const HighsInt iteration_counts[] = {0, 1, 2, 4, 5};
  for (std::size_t selection = 0; selection < 5; ++selection) {
    const HighsInt iteration_count = iteration_counts[selection];
    CAPTURE(iteration_count);
    hipdlp::algorithm::FixedTraceOptions disabled_options =
        makeOptions(iteration_count);
    disabled_options.backend = hipdlp::backend::Kind::kCuda;
    disabled_options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    enableFixedPointError(disabled_options, 1.0, 1.0);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enableFixedPointRestart(enabled_options, makeRestartReference(10.0), 40);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    CHECK(disabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(enabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(disabled_audit.graph_launch_calls == expected_launches[selection]);
    CHECK(enabled_audit.graph_launch_calls == expected_launches[selection]);
    checkSameResourceAudit(disabled_audit, enabled_audit);
    checkBitwiseEqualTraceOutputs(disabled, enabled);
    checkBitwiseEqualFixedPointError(disabled.fixed_point_error,
                                     enabled.fixed_point_error);
    checkSameDiagnostics(disabled, enabled);
    if (iteration_count == 0) {
      checkNoFixedPointRestart(enabled);
    } else {
      REQUIRE(enabled.fixed_point_restart.observed);
      REQUIRE(enabled.fixed_point_restart.evaluation.ok());
      CHECK(enabled.fixed_point_restart.evaluation.decision ==
            hipdlp::algorithm::RestartDecision::kContinue);
    }
  }
}
#endif

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-cadence-zero-is-final-only",
          "[pdlp][hipdlp][phase7k][restart][fixed-point][cadence][disabled]") {
  const hipdlp::algorithm::FixedTraceOptions defaults;
  CHECK(defaults.fixed_point_restart.check_cadence == 0);

  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions error_only_options = makeOptions(4);
    error_only_options.backend = selection.backend;
    error_only_options.execution_mode = selection.mode;
    enableFixedPointError(error_only_options, 1.0, 1.0);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult error_only =
        runTrace(makeOneStepLp(), error_only_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit error_only_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif

    hipdlp::algorithm::FixedTraceOptions default_options = error_only_options;
    enableFixedPointRestart(default_options, makeRestartReference(10.0), 40);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult defaulted =
        runTrace(makeOneStepLp(), default_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit defaulted_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif

    hipdlp::algorithm::FixedTraceOptions explicit_options = error_only_options;
    enableFixedPointRestartCadence(explicit_options, makeRestartReference(10.0),
                                   40, 0);
    const hipdlp::algorithm::FixedTraceResult explicitly_zero =
        runTrace(makeOneStepLp(), explicit_options);

    INFO(error_only.error.message);
    INFO(defaulted.error.message);
    INFO(explicitly_zero.error.message);
    REQUIRE(error_only.ok());
    REQUIRE(defaulted.ok());
    REQUIRE(explicitly_zero.ok());
    checkNoFixedPointRestart(error_only);

    REQUIRE(defaulted.fixed_point_restart_checks.size() == 1);
    CHECK(defaulted.fixed_point_restart_checks[0].iteration == 4);
    CHECK(bitwiseEqual(
        defaulted.fixed_point_restart_checks[0].error.fixed_point_error,
        defaulted.fixed_point_error.evaluation.fixed_point_error));
    CHECK(defaulted.fixed_point_restart_checks[0].evaluation.decision ==
          defaulted.fixed_point_restart.evaluation.decision);
    CHECK(defaulted.fixed_point_restart_checks[0].evaluation.condition ==
          defaulted.fixed_point_restart.evaluation.condition);

    checkBitwiseEqual(defaulted, explicitly_zero);
    checkSameDiagnostics(defaulted, explicitly_zero);
    checkBitwiseEqualTraceOutputs(error_only, defaulted);
    checkBitwiseEqualFixedPointError(error_only.fixed_point_error,
                                     defaulted.fixed_point_error);
    checkSameDiagnostics(error_only, defaulted);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    if (selection.backend == hipdlp::backend::Kind::kCuda)
      checkSameResourceAudit(error_only_audit, defaulted_audit);
#endif
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-cadence-one-records-five",
          "[pdlp][hipdlp][phase7k][restart][fixed-point][cadence][exact]") {
  // The first restart lands at the third check, so a restart-free prefix trace
  // is still an independent oracle for the first three errors.
  const double errors[] = {prefixFixedPointError(1), prefixFixedPointError(2),
                           prefixFixedPointError(3)};
  const double expected_errors[] = {3.0 / 4.0, 11.0 / 16.0, 25.0 / 48.0};
  for (std::size_t index = 0; index < 3; ++index) {
    CAPTURE(index);
    REQUIRE(errors[index] == Approx(expected_errors[index]).margin(1e-15));
  }
  // The fourth error is measured from the restarted iterate and stays outside
  // every decay window; the fifth collapses onto the bound and enters the
  // sufficient-decay window of the reset reference.
  const double restarted_errors[] = {5.0 / 32.0, 0.0};
  REQUIRE(restarted_errors[0] > 0.2 * expected_errors[2]);
  REQUIRE(restarted_errors[1] <= 0.2 * expected_errors[2]);

  hipdlp::algorithm::FixedTraceOptions options = makeOptions(5);
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestartCadence(options, makeRestartReference(10.0), 3, 1);

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 6);
  REQUIRE(result.fixed_point_restart_checks.size() == 5);

  const double expected_check_errors[] = {
      expected_errors[0], expected_errors[1], expected_errors[2],
      restarted_errors[0], restarted_errors[1]};
  const hipdlp::algorithm::RestartDecision expected_decisions[] = {
      hipdlp::algorithm::RestartDecision::kContinue,
      hipdlp::algorithm::RestartDecision::kContinue,
      hipdlp::algorithm::RestartDecision::kRestart,
      hipdlp::algorithm::RestartDecision::kContinue,
      hipdlp::algorithm::RestartDecision::kRestart};
  const hipdlp::algorithm::FixedPointRestartCondition expected_conditions[] = {
      hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered,
      hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered,
      hipdlp::algorithm::FixedPointRestartCondition::kFirstCheckInterval,
      hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered,
      hipdlp::algorithm::FixedPointRestartCondition::kSufficientDecay};

  for (std::size_t check = 0; check < 5; ++check) {
    CAPTURE(check);
    const hipdlp::algorithm::FixedTraceFixedPointRestartCheck& record =
        result.fixed_point_restart_checks[check];
    INFO(record.error.error.message);
    INFO(record.evaluation.error.message);
    CHECK(record.iteration == static_cast<HighsInt>(check) + 1);
    REQUIRE(record.error.ok());
    CHECK(record.error.fixed_point_error ==
          Approx(expected_check_errors[check]).margin(1e-15));
    CHECK(record.error.diagnostics.device_to_host_scalar_readbacks == 3);
    REQUIRE(record.evaluation.ok());
    CHECK(record.evaluation.decision == expected_decisions[check]);
    CHECK(record.evaluation.condition == expected_conditions[check]);
  }

  CHECK(result.fixed_point_restart_checks[0].evaluation.next.reference_error ==
        10.0);
  CHECK(result.fixed_point_restart_checks[1]
            .evaluation.next.previous_trial_error == errors[1]);
  CHECK(result.fixed_point_restart_checks[3].evaluation.next.reference_error ==
        errors[2]);

  REQUIRE(result.fixed_point_error.observed);
  CHECK(result.fixed_point_error.evaluation.fixed_point_error ==
        Approx(restarted_errors[1]).margin(1e-15));
  REQUIRE(result.fixed_point_restart.observed);
  CHECK(result.fixed_point_restart.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(result.fixed_point_restart.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kSufficientDecay);
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-chains-the-reset-reference",
          "[pdlp][hipdlp][phase7k][restart][fixed-point][cadence][chaining]") {
  const double first_error = prefixFixedPointError(1);
  const double second_error = prefixFixedPointError(2);
  REQUIRE(first_error == 0.75);
  // The original reference would take the sufficient-decay branch outright.
  REQUIRE(second_error <= 0.2 * 10.0);
  // The reset reference leaves every decay window closed.
  REQUIRE(second_error > 0.8 * first_error);

  hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestartCadence(options, makeRestartReference(10.0), 1, 1);
  options.fixed_point_restart.parameters.artificial_threshold = 0.6;

  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.fixed_point_restart_checks.size() == 2);

  const hipdlp::algorithm::FixedTraceFixedPointRestartCheck& first =
      result.fixed_point_restart_checks[0];
  CHECK(first.iteration == 1);
  CHECK(first.error.fixed_point_error == first_error);
  CHECK(first.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(first.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kFirstCheckInterval);
  CHECK_FALSE(first.evaluation.next.reference_valid);

  const hipdlp::algorithm::FixedTraceFixedPointRestartCheck& second =
      result.fixed_point_restart_checks[1];
  CHECK(second.iteration == 2);
  CHECK(second.error.fixed_point_error == second_error);
  CHECK(second.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(second.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
  CHECK(second.evaluation.next.reference_valid);
  CHECK(second.evaluation.next.reference_error == first_error);
  CHECK(second.evaluation.next.reference_error != 10.0);
  CHECK(second.evaluation.next.previous_trial_error == second_error);
  CHECK(options.fixed_point_restart.previous.reference_error == 10.0);
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-halpern-and-total-diverge",
          "[pdlp][hipdlp][phase7k][restart][fixed-point][cadence][counters]") {
  // A restart at iteration 1 leaves halpern_iterations at 1 while
  // total_iterations reaches 2, so the artificial boundary 0.36 * total
  // straddles the two counters between thresholds 0.4 and 0.6.
  struct Expectation {
    double artificial_threshold;
    hipdlp::algorithm::RestartDecision decision;
    hipdlp::algorithm::FixedPointRestartCondition condition;
  };
  const Expectation expectations[] = {
      {0.4, hipdlp::algorithm::RestartDecision::kRestart,
       hipdlp::algorithm::FixedPointRestartCondition::kArtificialThreshold},
      {0.6, hipdlp::algorithm::RestartDecision::kContinue,
       hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered}};

  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(selection);
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    enableFixedPointError(options, 1.0, 1.0);
    enableFixedPointRestartCadence(options, makeRestartReference(10.0), 1, 1);
    options.fixed_point_restart.parameters.artificial_threshold =
        expectations[selection].artificial_threshold;

    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.fixed_point_restart_checks.size() == 2);
    CHECK(result.fixed_point_restart_checks[0].evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kRestart);
    CHECK(result.fixed_point_restart_checks[1].iteration == 2);
    CHECK(result.fixed_point_restart_checks[1].evaluation.decision ==
          expectations[selection].decision);
    CHECK(result.fixed_point_restart_checks[1].evaluation.condition ==
          expectations[selection].condition);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-cadence-is-transactional",
          "[pdlp][hipdlp][phase7k][restart][fixed-point][cadence][errors]"
          "[transactional]") {
  hipdlp::algorithm::FixedTraceOptions survivor = makeOptions(1);
  enableFixedPointError(survivor, 1.0, 1.0);
  enableFixedPointRestartCadence(survivor, makeRestartReference(10.0), 1, 1);
  survivor.fixed_point_restart.parameters.artificial_threshold =
      std::numeric_limits<double>::max();
  const hipdlp::algorithm::FixedTraceResult first_check_only =
      runTrace(makeOneStepLp(), survivor);
  INFO(first_check_only.error.message);
  REQUIRE(first_check_only.ok());
  REQUIRE(first_check_only.fixed_point_restart_checks.size() == 1);

  // The same parameters overflow at the second of three checks because the
  // artificial boundary scales with total_iterations.
  hipdlp::algorithm::FixedTraceOptions options = survivor;
  options.iteration_count = 3;
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  checkTransactionalFailure(result);
  CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
  CHECK(result.error.message ==
        "HiPDLP fixed-point restart thresholds are outside FP64 range");
  CHECK(result.diagnostics.transfer_calls_after_iterations == 0);

  SECTION("a negative cadence is rejected before any backend resource") {
    hipdlp::algorithm::FixedTraceOptions negative = makeOptions(2);
    enableFixedPointError(negative, 1.0, 1.0);
    enableFixedPointRestartCadence(negative, makeRestartReference(10.0), 1, -1);
    const hipdlp::algorithm::FixedTraceResult rejected =
        runTrace(makeOneStepLp(), negative);
    checkTransactionalFailure(rejected);
    CHECK(rejected.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(rejected.diagnostics.resource_creation_calls_before_iterations == 0);
  }

  SECTION("a mid-block cadence is rejected in recorded mode") {
    hipdlp::algorithm::FixedTraceOptions recorded = makeOptions(5);
    recorded.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    enableFixedPointError(recorded, 1.0, 1.0);
    enableFixedPointRestartCadence(recorded, makeRestartReference(10.0), 40, 1);
    const hipdlp::algorithm::FixedTraceResult rejected =
        runTrace(makeOneStepLp(), recorded);
    checkTransactionalFailure(rejected);
    CHECK(rejected.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(rejected.diagnostics.resource_creation_calls_before_iterations == 0);

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::algorithm::FixedTraceOptions aligned = recorded;
    aligned.backend = hipdlp::backend::Kind::kCuda;
    aligned.fixed_point_restart.check_cadence = 2;
    const hipdlp::algorithm::FixedTraceResult accepted =
        runTrace(makeOneStepLp(), aligned);
    INFO(accepted.error.message);
    REQUIRE(accepted.ok());
    REQUIRE(accepted.fixed_point_restart_checks.size() == 3);
#endif
  }
}

TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-cadence-is-deterministic",
          "[pdlp][hipdlp][phase7k][restart][fixed-point][cadence]"
          "[determinism]") {
  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions final_only_options = makeOptions(4);
    final_only_options.backend = selection.backend;
    final_only_options.execution_mode = selection.mode;
    enableFixedPointError(final_only_options, 1.0, 1.0);
    enableFixedPointRestart(final_only_options, makeRestartReference(10.0), 2);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult final_only =
        runTrace(makeOneStepLp(), final_only_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit final_only_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif

    hipdlp::algorithm::FixedTraceOptions cadence_options = final_only_options;
    cadence_options.fixed_point_restart.check_cadence = 2;
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), cadence_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit cadence_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), cadence_options);

    INFO(final_only.error.message);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(final_only.ok());
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    REQUIRE(final_only.fixed_point_restart_checks.size() == 1);
    REQUIRE(first.fixed_point_restart_checks.size() == 2);
    checkBitwiseEqual(first, second);
    checkSameDiagnostics(first, second);
    // The extra periodic restart moves the trace off the final-only path.
    CHECK(final_only.points[2].primal != first.points[2].primal);
    checkSameDiagnostics(final_only, first);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    if (selection.backend == hipdlp::backend::Kind::kCuda)
      checkFixedPointCheckAuditDelta(final_only_audit, cadence_audit, 1);
#endif
  }
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-fixed-trace-fixed-point-restart-cadence-backend-agreement",
          "[pdlp][hipdlp][phase7k][restart][fixed-point][cadence][cuda]"
          "[graph][replay][parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  const std::size_t expected_launches[] = {0, 1, 1, 2, 3};
  const HighsInt iteration_counts[] = {0, 1, 2, 4, 5};
  const std::size_t expected_checks[] = {0, 1, 1, 2, 3};
  for (std::size_t selection = 0; selection < 5; ++selection) {
    const HighsInt iteration_count = iteration_counts[selection];
    CAPTURE(iteration_count);
    hipdlp::algorithm::FixedTraceOptions cpu_options =
        makeOptions(iteration_count);
    enableFixedPointError(cpu_options, 1.0, 1.0);
    enableFixedPointRestartCadence(cpu_options, makeRestartReference(10.0), 2,
                                   2);
    const hipdlp::algorithm::FixedTraceResult cpu =
        runTrace(makeOneStepLp(), cpu_options);
    INFO(cpu.error.message);
    REQUIRE(cpu.ok());
    REQUIRE(cpu.fixed_point_restart_checks.size() ==
            expected_checks[selection]);

    hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
    direct_options.backend = hipdlp::backend::Kind::kCuda;
    const hipdlp::algorithm::FixedTraceResult direct =
        runTrace(makeOneStepLp(), direct_options);

    hipdlp::algorithm::FixedTraceOptions recorded_options = direct_options;
    recorded_options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult recorded =
        runTrace(makeOneStepLp(), recorded_options);
    const hipdlp::backend::cuda::testing::ResourceAudit recorded_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions baseline_options = recorded_options;
    baseline_options.fixed_point_restart.check_cadence = 0;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult baseline =
        runTrace(makeOneStepLp(), baseline_options);
    const hipdlp::backend::cuda::testing::ResourceAudit baseline_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    INFO(direct.error.message);
    INFO(recorded.error.message);
    INFO(baseline.error.message);
    REQUIRE(direct.ok());
    REQUIRE(recorded.ok());
    REQUIRE(baseline.ok());
    checkSameRestartDecisions(cpu, direct);
    checkSameRestartDecisions(cpu, recorded);
    checkApproxTrace(recorded, cpu, 5e-13);
    for (std::size_t check = 0; check < cpu.fixed_point_restart_checks.size();
         ++check) {
      CAPTURE(check);
      CHECK(
          recorded.fixed_point_restart_checks[check].error.fixed_point_error ==
          Approx(cpu.fixed_point_restart_checks[check].error.fixed_point_error)
              .margin(5e-13));
    }

    CHECK(recorded.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(recorded_audit.graph_launch_calls == expected_launches[selection]);
    CHECK(baseline.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(baseline_audit.graph_launch_calls == expected_launches[selection]);
    checkSameDiagnostics(baseline, recorded);
    const std::size_t baseline_checks = expected_checks[selection] == 0 ? 0 : 1;
    REQUIRE(baseline.fixed_point_restart_checks.size() == baseline_checks);
    if (expected_checks[selection] == baseline_checks) {
      checkBitwiseEqualTraceOutputs(baseline, recorded);
    } else {
      // The extra periodic restarts move the cadence trace off the final-only
      // path from their boundary onward.
      CHECK(baseline.points[2].primal != recorded.points[2].primal);
    }
    checkFixedPointCheckAuditDelta(
        baseline_audit, recorded_audit,
        expected_checks[selection] - baseline_checks);
  }
}
#endif
TEST_CASE("hi-pdlp-fixed-trace-uses-projected-candidates-for-legacy-error",
          "[pdlp][hipdlp][phase7m][fixed-point-error][candidate][exact]") {
  SECTION("one step") {
    const double candidate_primal = 0.5;
    const double candidate_dual = 0.25;
    const double reflected_primal = 1.0;
    const double reflected_dual = 0.5;
    const double delta_primal = candidate_primal - reflected_primal;
    const double delta_dual = candidate_dual - reflected_dual;
    const double expected =
        std::sqrt(delta_primal * delta_primal + delta_dual * delta_dual +
                  2.0 * delta_primal * delta_dual);
    REQUIRE(expected == 3.0 / 4.0);

    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    enableFixedPointError(options, 1.0, 1.0);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.fixed_point_error.observed);
    REQUIRE(result.fixed_point_error.evaluation.ok());
    CHECK(result.fixed_point_error.evaluation.fixed_point_error == expected);
    CHECK(result.fixed_point_error.evaluation.diagnostics
              .device_to_host_scalar_readbacks == 3);
  }

  SECTION("two steps separate the candidate from the Halpern current") {
    const double candidate_primal = 9.0 / 8.0;
    const double candidate_dual = 5.0 / 16.0;
    const double reflected_primal = 7.0 / 4.0;
    const double reflected_dual = 3.0 / 8.0;
    const double candidate_delta_primal = candidate_primal - reflected_primal;
    const double candidate_delta_dual = candidate_dual - reflected_dual;
    const double expected =
        std::sqrt(candidate_delta_primal * candidate_delta_primal +
                  candidate_delta_dual * candidate_delta_dual +
                  2.0 * candidate_delta_primal * candidate_delta_dual);
    REQUIRE(expected == 11.0 / 16.0);

    const double halpern_primal = (2.0 / 3.0) * reflected_primal;
    const double halpern_dual = (2.0 / 3.0) * reflected_dual;
    const double halpern_delta_primal = halpern_primal - reflected_primal;
    const double halpern_delta_dual = halpern_dual - reflected_dual;
    const double current_based_error =
        std::sqrt(halpern_delta_primal * halpern_delta_primal +
                  halpern_delta_dual * halpern_delta_dual +
                  2.0 * halpern_delta_primal * halpern_delta_dual);
    REQUIRE(current_based_error != expected);

    hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
    enableFixedPointError(options, 1.0, 1.0);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.fixed_point_error.observed);
    REQUIRE(result.fixed_point_error.evaluation.ok());
    CHECK(result.fixed_point_error.evaluation.fixed_point_error == expected);
    CHECK(result.fixed_point_error.evaluation.fixed_point_error !=
          current_based_error);
    CHECK(result.fixed_point_error.evaluation.diagnostics
              .device_to_host_scalar_readbacks == 3);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-cadence-uses-each-boundary-candidate",
          "[pdlp][hipdlp][phase7m][fixed-point-error][candidate][cadence]"
          "[exact]") {
  hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(4);
  const hipdlp::algorithm::FixedTraceResult disabled =
      runTrace(makeOneStepLp(), disabled_options);

  hipdlp::algorithm::FixedTraceOptions options = disabled_options;
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestartCadence(options, makeRestartReference(10.0), 100, 2);
  options.fixed_point_restart.parameters.artificial_threshold = 100.0;
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);

  INFO(disabled.error.message);
  INFO(result.error.message);
  REQUIRE(disabled.ok());
  REQUIRE(result.ok());
  checkBitwiseEqualTraceOutputs(disabled, result);
  checkSameDiagnostics(disabled, result);
  REQUIRE(result.fixed_point_restart_checks.size() == 2);

  const HighsInt expected_iterations[] = {2, 4};
  const double expected_errors[] = {11.0 / 16.0, 9.0 / 64.0};
  for (std::size_t check = 0; check < 2; ++check) {
    CAPTURE(check);
    const hipdlp::algorithm::FixedTraceFixedPointRestartCheck& observation =
        result.fixed_point_restart_checks[check];
    CHECK(observation.iteration == expected_iterations[check]);
    REQUIRE(observation.error.ok());
    CHECK(observation.error.fixed_point_error == expected_errors[check]);
    CHECK(observation.error.diagnostics.device_to_host_scalar_readbacks == 3);
  }
  REQUIRE(result.fixed_point_error.observed);
  CHECK(result.fixed_point_error.evaluation.fixed_point_error ==
        expected_errors[1]);
}

TEST_CASE("hi-pdlp-fixed-trace-candidate-error-is-bitwise-deterministic",
          "[pdlp][hipdlp][phase7m][fixed-point-error][candidate]"
          "[determinism]") {
  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(4);
    options.backend = selection.backend;
    options.execution_mode = selection.mode;
    enableFixedPointError(options, 1.0, 1.0);
    enableFixedPointRestartCadence(options, makeRestartReference(10.0), 100, 2);
    options.fixed_point_restart.parameters.artificial_threshold = 100.0;
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), options);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    checkBitwiseEqual(first, second);
    checkSameDiagnostics(first, second);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-candidate-error-failure-is-transactional",
          "[pdlp][hipdlp][phase7m][fixed-point-error][candidate][errors]"
          "[transactional]") {
  struct ExecutionSelection {
    hipdlp::backend::Kind backend;
    hipdlp::algorithm::FixedTraceExecutionMode mode;
  };
  std::vector<ExecutionSelection> selections;
  selections.push_back({hipdlp::backend::Kind::kCpu,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kDirect});
  selections.push_back({hipdlp::backend::Kind::kCuda,
                        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded});
#endif

  for (const ExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    options.backend = selection.backend;
    options.execution_mode = selection.mode;
    options.primal_step_size = 4.0;
    enableFixedPointError(options, std::numeric_limits<double>::max(), 1.0);
    enableFixedPointRestart(options, makeRestartReference(10.0), 1);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.error.message ==
          "HiPDLP fixed-point weighted movement is outside FP64 range");
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
  }
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-fixed-trace-candidate-error-backend-audit-parity",
          "[pdlp][hipdlp][phase7m][fixed-point-error][candidate][cuda]"
          "[resources][graph][replay][parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  hipdlp::algorithm::FixedTraceOptions cpu_options = makeOptions(4);
  enableFixedPointError(cpu_options, 1.0, 1.0);
  enableFixedPointRestartCadence(cpu_options, makeRestartReference(10.0), 100,
                                 2);
  cpu_options.fixed_point_restart.parameters.artificial_threshold = 100.0;
  const hipdlp::algorithm::FixedTraceResult cpu =
      runTrace(makeOneStepLp(), cpu_options);
  INFO(cpu.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(cpu.fixed_point_restart_checks.size() == 2);

  const hipdlp::algorithm::FixedTraceExecutionMode modes[] = {
      hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded};
  const std::size_t expected_launches[] = {0, 2};
  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(selection);
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(4);
    disabled_options.backend = hipdlp::backend::Kind::kCuda;
    disabled_options.execution_mode = modes[selection];
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions enabled_options = cpu_options;
    enabled_options.backend = hipdlp::backend::Kind::kCuda;
    enabled_options.execution_mode = modes[selection];
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    checkNoFixedPointError(disabled);
    checkNoFixedPointRestart(disabled);
    checkApproxTrace(enabled, cpu, 5e-13);
    checkSameRestartDecisions(enabled, cpu);
    REQUIRE(enabled.fixed_point_restart_checks.size() == 2);
    for (std::size_t check = 0; check < 2; ++check) {
      CAPTURE(check);
      CHECK(
          enabled.fixed_point_restart_checks[check].error.fixed_point_error ==
          Approx(cpu.fixed_point_restart_checks[check].error.fixed_point_error)
              .margin(5e-13));
      CHECK(enabled.fixed_point_restart_checks[check]
                .error.diagnostics.device_to_host_scalar_readbacks == 3);
    }
    checkBitwiseEqualTraceOutputs(disabled, enabled);
    checkSameDiagnostics(disabled, enabled);
    checkFixedPointCheckAuditDelta(disabled_audit, enabled_audit, 2);
    CHECK(disabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(enabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(disabled_audit.graph_launch_calls == expected_launches[selection]);
    CHECK(enabled_audit.graph_launch_calls == expected_launches[selection]);
  }
}
#endif

TEST_CASE("hi-pdlp-fixed-trace-applies-the-restart-state-transition",
          "[pdlp][hipdlp][phase7n][restart][transition][exact]") {
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), makeRestartAtTwoOptions(3));
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 4);
  REQUIRE(result.fixed_point_restart_checks.size() == 2);

  const hipdlp::algorithm::FixedTraceFixedPointRestartCheck& at_two =
      result.fixed_point_restart_checks[0];
  CHECK(at_two.iteration == 2);
  REQUIRE(at_two.error.ok());
  // The record keeps the pre-reset candidate error and decision.
  CHECK(at_two.error.fixed_point_error == kErrorAtTwo);
  CHECK(at_two.error.diagnostics.device_to_host_scalar_readbacks == 3);
  CHECK(at_two.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(at_two.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kFirstCheckInterval);

  // The boundary point is the post-decision state, not the Halpern blend
  // 2/3 * 7/4 that would otherwise continue.
  checkExactPoint(result.points[2], {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo}, {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo});
  CHECK(result.points[2].primal[0] != (1.0 - 1.0 / 3.0) * (7.0 / 4.0));
  checkExactPoint(result.points[3], {kRestartedPrimalAtThree},
                  {kRestartedDualAtThree}, {kRestartedPrimalAtThree},
                  {kRestartedDualAtThree});

  const hipdlp::algorithm::FixedTraceFixedPointRestartCheck& at_three =
      result.fixed_point_restart_checks[1];
  CHECK(at_three.iteration == 3);
  REQUIRE(at_three.error.ok());
  CHECK(at_three.error.fixed_point_error == kRestartedErrorAtThree);
  CHECK(at_three.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);
  CHECK(at_three.evaluation.condition ==
        hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered);
  CHECK(at_three.evaluation.next.reference_error == kErrorAtTwo);
  REQUIRE(result.fixed_point_error.observed);
  CHECK(result.fixed_point_error.evaluation.fixed_point_error ==
        kRestartedErrorAtThree);

  // Observation without application keeps running from the Halpern blend.
  hipdlp::algorithm::FixedTraceOptions observed_options = makeOptions(3);
  enableFixedPointError(observed_options, 1.0, 1.0);
  const hipdlp::algorithm::FixedTraceResult observed =
      runTrace(makeOneStepLp(), observed_options);
  INFO(observed.error.message);
  REQUIRE(observed.ok());
  CHECK(observed.points[1].primal == result.points[1].primal);
  CHECK(observed.points[2].primal != result.points[2].primal);
  CHECK(observed.points[3].primal != result.points[3].primal);
  CHECK(observed.fixed_point_error.evaluation.fixed_point_error !=
        kRestartedErrorAtThree);
}

TEST_CASE("hi-pdlp-fixed-trace-restart-adopts-candidates-as-current-and-anchor",
          "[pdlp][hipdlp][phase7n][restart][transition][anchor]") {
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), makeRestartAtTwoOptions(3));
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 4);

  // Current, Ax, and A^T y all follow the candidate at the boundary.
  CHECK(result.points[2].primal ==
        std::vector<double>({kCandidatePrimalAtTwo}));
  CHECK(result.points[2].dual == std::vector<double>({kCandidateDualAtTwo}));
  CHECK(result.points[2].ax == std::vector<double>({kCandidatePrimalAtTwo}));
  CHECK(result.points[2].aty == std::vector<double>({kCandidateDualAtTwo}));

  // The next step reflects about the adopted current and blends against the
  // adopted anchor with weight 1/2, so a stale anchor stays observable.
  const double reflected_primal =
      2.0 * kRestartedPrimalAtThree - kCandidatePrimalAtTwo;
  const double reflected_dual =
      2.0 * kRestartedDualAtThree - kCandidateDualAtTwo;
  REQUIRE(reflected_primal == 39.0 / 16.0);
  REQUIRE(reflected_dual == 3.0 / 32.0);
  CHECK(result.points[3].primal ==
        std::vector<double>(
            {0.5 * kCandidatePrimalAtTwo + 0.5 * reflected_primal}));
  CHECK(
      result.points[3].dual ==
      std::vector<double>({0.5 * kCandidateDualAtTwo + 0.5 * reflected_dual}));
  CHECK(result.points[3].primal[0] != 0.5 * reflected_primal);
  CHECK(result.points[3].dual[0] != 0.5 * reflected_dual);
}

TEST_CASE("hi-pdlp-fixed-trace-restart-separates-total-and-halpern-counters",
          "[pdlp][hipdlp][phase7n][restart][transition][counters]") {
  // A restart at iteration 2 leaves the Halpern age at 2 when the total count
  // reaches 4. The artificial boundary 0.4 * 4 = 1.6 is below that age while
  // 0.6 * 4 = 2.4 is above it; an unreset age of 4 or a reset total of 2 would
  // restart in both cases.
  struct Expectation {
    double artificial_threshold;
    hipdlp::algorithm::RestartDecision decision;
    hipdlp::algorithm::FixedPointRestartCondition condition;
  };
  const Expectation expectations[] = {
      {0.4, hipdlp::algorithm::RestartDecision::kRestart,
       hipdlp::algorithm::FixedPointRestartCondition::kArtificialThreshold},
      {0.6, hipdlp::algorithm::RestartDecision::kContinue,
       hipdlp::algorithm::FixedPointRestartCondition::kNotTriggered}};

  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(selection);
    hipdlp::algorithm::FixedTraceOptions options = makeRestartAtTwoOptions(4);
    options.fixed_point_restart.parameters.artificial_threshold =
        expectations[selection].artificial_threshold;
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    // Trace indices stay tied to the monotone total.
    REQUIRE(result.points.size() == 5);
    REQUIRE(result.fixed_point_restart_checks.size() == 2);
    CHECK(result.fixed_point_restart_checks[0].iteration == 2);
    CHECK(result.fixed_point_restart_checks[0].evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kRestart);
    CHECK(result.fixed_point_restart_checks[1].iteration == 4);
    REQUIRE(result.fixed_point_restart_checks[1].error.ok());
    CHECK(result.fixed_point_restart_checks[1].error.fixed_point_error ==
          kRestartedErrorAtFour);
    CHECK(result.fixed_point_restart_checks[1].evaluation.decision ==
          expectations[selection].decision);
    CHECK(result.fixed_point_restart_checks[1].evaluation.condition ==
          expectations[selection].condition);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-resets-weighted-averaging",
          "[pdlp][hipdlp][phase7n][restart][transition][averaging]") {
  hipdlp::algorithm::FixedTraceOptions options = makeRestartAtTwoOptions(3);
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  const double sample_weight =
      std::sqrt(options.primal_step_size * options.dual_step_size);
  REQUIRE(result.averaging.hasAverage());
  // Only the single post-restart sample survives.
  CHECK(result.averaging.sample_count == 1);
  CHECK(result.averaging.total_weight == sample_weight);
  REQUIRE(result.averaged_primal.size() == 1);
  REQUIRE(result.averaged_dual.size() == 1);
  CHECK(result.averaged_primal[0] ==
        Approx(kRestartedPrimalAtThree).margin(1e-15));
  CHECK(result.averaged_dual[0] == Approx(kRestartedDualAtThree).margin(1e-15));
  REQUIRE(result.averaged_convergence.observed);

  hipdlp::algorithm::FixedTraceOptions observed_options = makeOptions(3);
  enableFixedPointError(observed_options, 1.0, 1.0);
  const hipdlp::algorithm::FixedTraceResult observed =
      runTrace(makeOneStepLp(), observed_options);
  INFO(observed.error.message);
  REQUIRE(observed.ok());
  CHECK(observed.averaging.sample_count == 3);
  CHECK(observed.averaging.total_weight == 3.0 * sample_weight);
  CHECK(observed.averaged_primal[0] != result.averaged_primal[0]);
}

TEST_CASE("hi-pdlp-fixed-trace-final-restart-applies-and-clears-the-average",
          "[pdlp][hipdlp][phase7n][restart][transition][final]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
  enableFixedPointError(options, 1.0, 1.0);
  enableFixedPointRestart(options, makeRestartReference(10.0), 2);
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 3);
  REQUIRE(result.fixed_point_restart_checks.size() == 1);
  CHECK(result.fixed_point_restart_checks[0].iteration == 2);
  CHECK(result.fixed_point_restart_checks[0].error.fixed_point_error ==
        kErrorAtTwo);
  REQUIRE(result.fixed_point_restart.observed);
  CHECK(result.fixed_point_restart.evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  checkExactPoint(result.points[2], {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo}, {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo});
  // No post-restart sample exists yet, so no average is valid.
  checkNoAverage(result);
  checkNoObservation(result.averaged_convergence);
  REQUIRE(result.current_convergence.observed);
  REQUIRE(result.fixed_point_error.observed);
  CHECK(result.fixed_point_error.evaluation.fixed_point_error == kErrorAtTwo);
}

TEST_CASE("hi-pdlp-fixed-trace-restart-mixes-periodic-continue-and-restart",
          "[pdlp][hipdlp][phase7n][restart][transition][cadence]") {
  hipdlp::algorithm::FixedTraceOptions options = makeRestartAtTwoOptions(4);
  options.fixed_point_restart.parameters.artificial_threshold = 0.6;
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 5);
  REQUIRE(result.fixed_point_restart_checks.size() == 2);
  CHECK(result.fixed_point_restart_checks[0].error.fixed_point_error ==
        kErrorAtTwo);
  CHECK(result.fixed_point_restart_checks[0].evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(result.fixed_point_restart_checks[1].error.fixed_point_error ==
        kRestartedErrorAtFour);
  CHECK(result.fixed_point_restart_checks[1].evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);

  checkExactPoint(result.points[2], {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo}, {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo});
  CHECK(result.points[3].primal ==
        std::vector<double>({kRestartedPrimalAtThree}));
  // The continued fourth step still blends against the anchor adopted at two.
  const double weight = 1.0 / 3.0;
  const double reflected_primal = 71.0 / 32.0;
  const double reflected_dual = 3.0 / 32.0;
  CHECK(result.points[4].primal[0] == Approx(weight * kCandidatePrimalAtTwo +
                                             (1.0 - weight) * reflected_primal)
                                          .margin(1e-15));
  CHECK(result.points[4].dual[0] ==
        Approx(weight * kCandidateDualAtTwo + (1.0 - weight) * reflected_dual)
            .margin(1e-15));
  // Two post-restart samples remain in the average.
  REQUIRE(result.averaging.hasAverage());
  CHECK(result.averaging.sample_count == 2);
}

TEST_CASE("hi-pdlp-fixed-trace-restart-transition-preserves-untriggered-traces",
          "[pdlp][hipdlp][phase7n][restart][transition][disabled][bitwise]") {
  const std::vector<RestartExecutionSelection> selections =
      restartExecutionSelections();
  for (const RestartExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(4);
    disabled_options.backend = selection.backend;
    disabled_options.execution_mode = selection.mode;
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);

    hipdlp::algorithm::FixedTraceOptions error_only_options = disabled_options;
    enableFixedPointError(error_only_options, 1.0, 1.0);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult error_only =
        runTrace(makeOneStepLp(), error_only_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit error_only_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif

    hipdlp::algorithm::FixedTraceOptions untriggered_options =
        error_only_options;
    enableFixedPointRestartCadence(untriggered_options,
                                   makeRestartReference(10.0), 100, 2);
    untriggered_options.fixed_point_restart.parameters.artificial_threshold =
        100.0;
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult untriggered =
        runTrace(makeOneStepLp(), untriggered_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit untriggered_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif

    INFO(disabled.error.message);
    INFO(error_only.error.message);
    INFO(untriggered.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(error_only.ok());
    REQUIRE(untriggered.ok());
    checkNoFixedPointError(disabled);
    checkNoFixedPointRestart(disabled);
    checkNoFixedPointRestart(error_only);
    REQUIRE(untriggered.fixed_point_restart_checks.size() == 2);
    for (std::size_t check = 0; check < 2; ++check) {
      CAPTURE(check);
      CHECK(untriggered.fixed_point_restart_checks[check].evaluation.decision ==
            hipdlp::algorithm::RestartDecision::kContinue);
      CHECK(untriggered.fixed_point_restart_checks[check]
                .error.diagnostics.device_to_host_scalar_readbacks == 3);
    }
    checkBitwiseEqualTraceOutputs(disabled, error_only);
    checkBitwiseEqualTraceOutputs(error_only, untriggered);
    checkBitwiseEqualFixedPointError(error_only.fixed_point_error,
                                     untriggered.fixed_point_error);
    checkSameDiagnostics(disabled, error_only);
    checkSameDiagnostics(error_only, untriggered);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    if (selection.backend == hipdlp::backend::Kind::kCuda)
      checkFixedPointCheckAuditDelta(error_only_audit, untriggered_audit, 1);
#endif
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-transition-is-deterministic",
          "[pdlp][hipdlp][phase7n][restart][transition][determinism]") {
  const std::vector<RestartExecutionSelection> selections =
      restartExecutionSelections();
  for (const RestartExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions options = makeRestartAtTwoOptions(4);
    options.backend = selection.backend;
    options.execution_mode = selection.mode;
    options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), options);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    REQUIRE(first.fixed_point_restart_checks.size() == 2);
    CHECK(first.fixed_point_restart_checks[0].evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kRestart);
    CHECK(first.fixed_point_restart_checks[1].evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kContinue);
    checkBitwiseEqual(first, second);
    checkSameDiagnostics(first, second);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-transition-failure-is-transactional",
          "[pdlp][hipdlp][phase7n][restart][transition][errors]"
          "[transactional]") {
  const std::vector<RestartExecutionSelection> selections =
      restartExecutionSelections();
  for (const RestartExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    // The artificial boundary scales with the total count, so the same
    // parameters survive the restarting check at two and overflow at four.
    hipdlp::algorithm::FixedTraceOptions survivor = makeRestartAtTwoOptions(2);
    survivor.backend = selection.backend;
    survivor.execution_mode = selection.mode;
    survivor.fixed_point_restart.parameters.artificial_threshold =
        std::numeric_limits<double>::max() / 3.0;
    const hipdlp::algorithm::FixedTraceResult restarted =
        runTrace(makeOneStepLp(), survivor);
    INFO(restarted.error.message);
    REQUIRE(restarted.ok());
    REQUIRE(restarted.fixed_point_restart_checks.size() == 1);
    CHECK(restarted.fixed_point_restart_checks[0].evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kRestart);

    hipdlp::algorithm::FixedTraceOptions options = survivor;
    options.iteration_count = 4;
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.error.message ==
          "HiPDLP fixed-point restart thresholds are outside FP64 range");
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
  }

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  const hipdlp::backend::cuda::testing::FailurePoint failure_points[] = {
      hipdlp::backend::cuda::testing::FailurePoint::kGraphCapture,
      hipdlp::backend::cuda::testing::FailurePoint::kGraphInstantiation,
      hipdlp::backend::cuda::testing::FailurePoint::kGraphLaunch,
      hipdlp::backend::cuda::testing::FailurePoint::
          kGraphAsynchronousExecution};
  for (const hipdlp::backend::cuda::testing::FailurePoint failure_point :
       failure_points) {
    CAPTURE(static_cast<int>(failure_point));
    hipdlp::algorithm::FixedTraceOptions options = makeRestartAtTwoOptions(4);
    options.backend = hipdlp::backend::Kind::kCuda;
    options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    hipdlp::backend::cuda::testing::setFailurePoint(failure_point);
    const hipdlp::algorithm::FixedTraceResult failed =
        runTrace(makeOneStepLp(), options);
    hipdlp::backend::cuda::testing::setFailurePoint(
        hipdlp::backend::cuda::testing::FailurePoint::kNone);
    checkTransactionalFailure(failed);
    CHECK(failed.error.code == hipdlp::api::ErrorCode::kBackendFailure);
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
    REQUIRE(recovered.fixed_point_restart_checks.size() == 2);
  }
#endif
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-fixed-trace-restart-transition-backend-agreement",
          "[pdlp][hipdlp][phase7n][restart][transition][cuda][graph][replay]"
          "[parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  const HighsInt iteration_counts[] = {0, 1, 2, 4, 5};
  const std::size_t expected_launches[] = {0, 1, 1, 2, 3};
  const std::size_t expected_checks[] = {0, 1, 1, 2, 3};
  for (std::size_t selection = 0; selection < 5; ++selection) {
    const HighsInt iteration_count = iteration_counts[selection];
    CAPTURE(iteration_count);
    hipdlp::algorithm::FixedTraceOptions cpu_options =
        makeRestartAtTwoOptions(iteration_count);
    cpu_options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    const hipdlp::algorithm::FixedTraceResult cpu =
        runTrace(makeOneStepLp(), cpu_options);
    INFO(cpu.error.message);
    REQUIRE(cpu.ok());
    REQUIRE(cpu.fixed_point_restart_checks.size() ==
            expected_checks[selection]);

    hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
    direct_options.backend = hipdlp::backend::Kind::kCuda;
    const hipdlp::algorithm::FixedTraceResult direct =
        runTrace(makeOneStepLp(), direct_options);

    hipdlp::algorithm::FixedTraceOptions recorded_options = direct_options;
    recorded_options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult recorded =
        runTrace(makeOneStepLp(), recorded_options);
    const hipdlp::backend::cuda::testing::ResourceAudit recorded_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    INFO(direct.error.message);
    INFO(recorded.error.message);
    REQUIRE(direct.ok());
    REQUIRE(recorded.ok());
    checkSameRestartDecisions(cpu, direct);
    checkSameRestartDecisions(cpu, recorded);
    checkApproxTrace(direct, cpu, 5e-13);
    checkApproxTrace(recorded, cpu, 5e-13);
    for (std::size_t check = 0; check < cpu.fixed_point_restart_checks.size();
         ++check) {
      CAPTURE(check);
      CHECK(
          recorded.fixed_point_restart_checks[check].error.fixed_point_error ==
          Approx(cpu.fixed_point_restart_checks[check].error.fixed_point_error)
              .margin(5e-13));
      CHECK(recorded.fixed_point_restart_checks[check]
                .error.diagnostics.device_to_host_scalar_readbacks == 3);
    }
    CHECK(recorded.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(recorded_audit.graph_launch_calls == expected_launches[selection]);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-transition-adds-only-restart-work",
          "[pdlp][hipdlp][phase7n][restart][transition][cuda][resources]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  const hipdlp::algorithm::FixedTraceExecutionMode modes[] = {
      hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded};
  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(selection);
    hipdlp::algorithm::FixedTraceOptions untriggered_options =
        makeRestartAtTwoOptions(4);
    untriggered_options.backend = hipdlp::backend::Kind::kCuda;
    untriggered_options.execution_mode = modes[selection];
    untriggered_options.fixed_point_restart.parameters.check_interval = 100;
    untriggered_options.fixed_point_restart.parameters.artificial_threshold =
        100.0;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult untriggered =
        runTrace(makeOneStepLp(), untriggered_options);
    const hipdlp::backend::cuda::testing::ResourceAudit untriggered_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions restarted_options =
        untriggered_options;
    restarted_options.fixed_point_restart.parameters.check_interval = 2;
    restarted_options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult restarted =
        runTrace(makeOneStepLp(), restarted_options);
    const hipdlp::backend::cuda::testing::ResourceAudit restarted_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    INFO(untriggered.error.message);
    INFO(restarted.error.message);
    REQUIRE(untriggered.ok());
    REQUIRE(restarted.ok());
    REQUIRE(untriggered.fixed_point_restart_checks.size() == 2);
    REQUIRE(restarted.fixed_point_restart_checks.size() == 2);
    CHECK(restarted.fixed_point_restart_checks[0].evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kRestart);
    CHECK(restarted.fixed_point_restart_checks[1].evaluation.decision ==
          hipdlp::algorithm::RestartDecision::kContinue);
    REQUIRE(untriggered.averaging.hasAverage());
    REQUIRE(restarted.averaging.hasAverage());
    checkSameDiagnostics(untriggered, restarted);
    checkRestartTransitionAuditDelta(untriggered_audit, restarted_audit, 0);
  }
}
#endif

TEST_CASE("hi-pdlp-fixed-trace-periodic-convergence-uses-candidates",
          "[pdlp][hipdlp][phase7o][convergence][cadence][candidate][exact]") {
  hipdlp::algorithm::FixedTraceOptions options = makeOptions(2);
  enablePeriodicConvergence(options, 1);
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.convergence_checks.size() == 2);
  REQUIRE(result.points.size() == 3);

  const hipdlp::algorithm::FixedTraceConvergenceCheck& at_one =
      result.convergence_checks[0];
  CHECK(at_one.iteration == 1);
  CHECK(at_one.raw.primal_residual == 3.0 / 2.0);
  CHECK(at_one.raw.dual_residual == 0.0);
  CHECK(at_one.raw.primal_objective == -1.0 / 2.0);
  CHECK(at_one.raw.dual_objective == -2.0);
  CHECK(at_one.raw.objective_gap == 3.0 / 2.0);
  CHECK(at_one.raw.unscaled_rhs_norm == 2.0);
  CHECK(at_one.raw.unscaled_cost_norm == 1.0);
  REQUIRE(at_one.stopping.ok());
  CHECK(at_one.stopping.metrics.relative_primal == 1.0 / 2.0);
  CHECK(at_one.stopping.metrics.relative_dual == 0.0);
  CHECK(at_one.stopping.metrics.relative_gap == 3.0 / 7.0);
  CHECK(at_one.stopping.decision ==
        hipdlp::algorithm::StoppingDecision::kContinue);
  // At the first step the Halpern blend is algebraically the candidate.
  checkExactPoint(result.points[1], {1.0 / 2.0}, {1.0 / 4.0}, {1.0 / 2.0},
                  {1.0 / 4.0});

  const hipdlp::algorithm::FixedTraceConvergenceCheck& at_two =
      result.convergence_checks[1];
  CHECK(at_two.iteration == 2);
  CHECK(at_two.raw.primal_residual == 7.0 / 8.0);
  CHECK(at_two.raw.dual_residual == 0.0);
  CHECK(at_two.raw.primal_objective == -9.0 / 8.0);
  CHECK(at_two.raw.dual_objective == -2.0);
  CHECK(at_two.raw.objective_gap == 7.0 / 8.0);
  CHECK(at_two.raw.unscaled_rhs_norm == 2.0);
  CHECK(at_two.raw.unscaled_cost_norm == 1.0);
  REQUIRE(at_two.stopping.ok());
  CHECK(at_two.stopping.metrics.relative_primal == 7.0 / 24.0);
  CHECK(at_two.stopping.metrics.relative_dual == 0.0);
  CHECK(at_two.stopping.metrics.relative_gap == 7.0 / 33.0);
  CHECK(at_two.stopping.decision ==
        hipdlp::algorithm::StoppingDecision::kContinue);

  // The same boundary exposes the Halpern current when no restart fires.
  checkExactPoint(result.points[2], {7.0 / 6.0}, {1.0 / 4.0}, {7.0 / 6.0},
                  {1.0 / 4.0});
  CHECK(at_two.raw.primal_residual != 5.0 / 6.0);
  CHECK(at_two.raw.primal_objective != -7.0 / 6.0);
  CHECK(at_two.raw.objective_gap != 5.0 / 6.0);
  CHECK(at_two.stopping.metrics.relative_primal != 5.0 / 18.0);
  CHECK(at_two.stopping.metrics.relative_gap != 1.0 / 5.0);
}

TEST_CASE("hi-pdlp-fixed-trace-periodic-convergence-cadence-is-observing",
          "[pdlp][hipdlp][phase7o][convergence][cadence][stopping]") {
  SECTION("positive cadence includes a nonmultiple final boundary") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(5);
    enablePeriodicConvergence(options, 2);
    options.convergence_tolerances.primal = 1.0;
    options.convergence_tolerances.dual = 1.0;
    options.convergence_tolerances.gap = 1.0;
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.points.size() == 6);
    REQUIRE(result.convergence_checks.size() == 3);
    CHECK(result.convergence_checks[0].iteration == 2);
    CHECK(result.convergence_checks[1].iteration == 4);
    CHECK(result.convergence_checks[2].iteration == 5);
    CHECK(result.convergence_checks[0].stopping.decision ==
          hipdlp::algorithm::StoppingDecision::kOptimal);
    CHECK(result.convergence_checks[2].stopping.ok());
    REQUIRE(result.current_convergence.observed);
    REQUIRE(result.averaged_convergence.observed);
  }

  SECTION("a final multiple is not duplicated") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(4);
    enablePeriodicConvergence(options, 2);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.convergence_checks.size() == 2);
    CHECK(result.convergence_checks[0].iteration == 2);
    CHECK(result.convergence_checks[1].iteration == 4);
  }

  SECTION("zero cadence observes only the final positive iteration") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(5);
    enablePeriodicConvergence(options, 0);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.convergence_checks.size() == 1);
    CHECK(result.convergence_checks[0].iteration == 5);
  }

  SECTION("zero iterations produce no periodic observation") {
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(0);
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enablePeriodicConvergence(enabled_options, 0);
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    CHECK(enabled.convergence_checks.empty());
    checkBitwiseEqual(disabled, enabled);
    checkPeriodicConvergenceDiagnosticsDelta(disabled, enabled, 0);
  }
}

TEST_CASE("hi-pdlp-periodic-convergence-precedes-the-restart-transition",
          "[pdlp][hipdlp][phase7o][convergence][restart][ordering]") {
  hipdlp::algorithm::FixedTraceOptions baseline_options =
      makeRestartAtTwoOptions(4);
  baseline_options.fixed_point_restart.parameters.artificial_threshold = 0.6;
  const hipdlp::algorithm::FixedTraceResult baseline =
      runTrace(makeOneStepLp(), baseline_options);

  hipdlp::algorithm::FixedTraceOptions enabled_options = baseline_options;
  enablePeriodicConvergence(enabled_options, 2);
  const hipdlp::algorithm::FixedTraceResult enabled =
      runTrace(makeOneStepLp(), enabled_options);

  INFO(baseline.error.message);
  INFO(enabled.error.message);
  REQUIRE(baseline.ok());
  REQUIRE(enabled.ok());
  checkBitwiseEqual(baseline, enabled);
  checkPeriodicConvergenceDiagnosticsDelta(baseline, enabled, 2);
  REQUIRE(enabled.convergence_checks.size() == 2);
  REQUIRE(enabled.fixed_point_restart_checks.size() == 2);
  CHECK(enabled.convergence_checks[0].iteration == 2);
  CHECK(enabled.convergence_checks[0].raw.primal_residual == 7.0 / 8.0);
  CHECK(enabled.convergence_checks[0].raw.primal_objective == -9.0 / 8.0);
  CHECK(enabled.convergence_checks[0].raw.objective_gap == 7.0 / 8.0);
  CHECK(enabled.fixed_point_restart_checks[0].iteration == 2);
  CHECK(enabled.fixed_point_restart_checks[0].evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  // Tier A14 exposes the adopted candidate at a restarting boundary.
  checkExactPoint(enabled.points[2], {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo}, {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo});

  hipdlp::algorithm::FixedTraceOptions observed_options = makeOptions(4);
  enablePeriodicConvergence(observed_options, 2);
  const hipdlp::algorithm::FixedTraceResult observed =
      runTrace(makeOneStepLp(), observed_options);
  INFO(observed.error.message);
  REQUIRE(observed.ok());
  REQUIRE(observed.convergence_checks.size() == 2);
  CHECK(bitwiseEqual(enabled.convergence_checks[0].raw.primal_residual,
                     observed.convergence_checks[0].raw.primal_residual));
  CHECK(bitwiseEqual(enabled.convergence_checks[0].raw.primal_objective,
                     observed.convergence_checks[0].raw.primal_objective));
  CHECK(bitwiseEqual(enabled.convergence_checks[0].raw.objective_gap,
                     observed.convergence_checks[0].raw.objective_gap));
  // The reset still changes the trajectory after the observed boundary.
  CHECK(enabled.points[3].primal != observed.points[3].primal);
}

TEST_CASE("hi-pdlp-periodic-convergence-default-and-validation",
          "[pdlp][hipdlp][phase7o][convergence][disabled][errors]") {
  SECTION("the default and explicitly disabled modes are bitwise identical") {
    hipdlp::algorithm::FixedTraceOptions default_options =
        makeRestartAtTwoOptions(4);
    default_options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    hipdlp::algorithm::FixedTraceOptions disabled_options = default_options;
    disabled_options.convergence.periodic_enabled = false;
    disabled_options.convergence.check_cadence = 1;
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), default_options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), disabled_options);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    CHECK(first.convergence_checks.empty());
    CHECK(second.convergence_checks.empty());
    checkBitwiseEqual(first, second);
    checkSameDiagnostics(first, second);
  }

  SECTION("negative cadence is rejected before backend resources") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(4);
    enablePeriodicConvergence(options, -1);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "fixed trace convergence check cadence must be nonnegative");
    CHECK(result.diagnostics.resource_creation_calls_before_iterations == 0);
    CHECK(result.diagnostics.transfer_calls_before_iterations == 0);
  }

  SECTION("recorded mid-block cadence is rejected before backend resources") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(4);
    options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    enablePeriodicConvergence(options, 1);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "fixed trace convergence check cadence must align with recorded "
          "iteration blocks");
    CHECK(result.diagnostics.resource_creation_calls_before_iterations == 0);
    CHECK(result.diagnostics.transfer_calls_before_iterations == 0);
    CHECK(result.diagnostics.command_buffer.replay_calls == 0);
  }
}

TEST_CASE("hi-pdlp-periodic-convergence-is-deterministic-and-transactional",
          "[pdlp][hipdlp][phase7o][convergence][determinism]"
          "[transactional]") {
  const std::vector<RestartExecutionSelection> selections =
      restartExecutionSelections();
  for (const RestartExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(5);
    options.backend = selection.backend;
    options.execution_mode = selection.mode;
    enablePeriodicConvergence(options, 2);
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeOneStepLp(), options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeOneStepLp(), options);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    checkBitwiseEqual(first, second);
    checkBitwiseEqualConvergenceChecks(first.convergence_checks,
                                       second.convergence_checks);
    checkSameDiagnostics(first, second);
  }

  SECTION("stopping failure wins over a shared failing restart check") {
    hipdlp::algorithm::FixedTraceOptions options = makeRestartAtTwoOptions(2);
    enablePeriodicConvergence(options, 2);
    options.convergence_tolerances.primal =
        std::numeric_limits<double>::quiet_NaN();
    options.fixed_point_restart.parameters.artificial_threshold =
        std::numeric_limits<double>::max();
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "HiPDLP convergence tolerances must be finite and positive");
    CHECK(result.diagnostics.transfer_calls_during_iterations == 4);
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
  }

  SECTION("non-finite candidate metrics expose no partial result") {
    HighsLp lp = makeOneStepLp();
    lp.row_lower_[0] = std::numeric_limits<double>::max();
    lp.row_upper_[0] = std::numeric_limits<double>::max();
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(1);
    enablePeriodicConvergence(options, 1);
    const hipdlp::algorithm::FixedTraceResult result = runTrace(lp, options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.error.message ==
          "fixed trace convergence metrics are non-finite");
    CHECK(result.diagnostics.transfer_calls_during_iterations == 4);
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
  }
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-periodic-convergence-backend-and-budget-agreement",
          "[pdlp][hipdlp][phase7o][convergence][cuda][graph][replay]"
          "[resources][parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  hipdlp::algorithm::FixedTraceOptions cpu_options = makeOptions(4);
  enablePeriodicConvergence(cpu_options, 2);
  const hipdlp::algorithm::FixedTraceResult cpu =
      runTrace(makeOneStepLp(), cpu_options);
  INFO(cpu.error.message);
  REQUIRE(cpu.ok());
  REQUIRE(cpu.convergence_checks.size() == 2);

  const hipdlp::algorithm::FixedTraceExecutionMode modes[] = {
      hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded};
  const std::size_t expected_launches[] = {0, 2};
  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(selection);
    hipdlp::algorithm::FixedTraceOptions disabled_options = makeOptions(4);
    disabled_options.backend = hipdlp::backend::Kind::kCuda;
    disabled_options.execution_mode = modes[selection];
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions explicit_disabled_options =
        disabled_options;
    explicit_disabled_options.convergence.periodic_enabled = false;
    explicit_disabled_options.convergence.check_cadence = 1;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult explicit_disabled =
        runTrace(makeOneStepLp(), explicit_disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit
        explicit_disabled_audit =
            hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enablePeriodicConvergence(enabled_options, 2);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    INFO(disabled.error.message);
    INFO(explicit_disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(explicit_disabled.ok());
    REQUIRE(enabled.ok());
    CHECK(disabled.convergence_checks.empty());
    CHECK(explicit_disabled.convergence_checks.empty());
    REQUIRE(enabled.convergence_checks.size() == 2);
    checkBitwiseEqual(disabled, explicit_disabled);
    checkSameDiagnostics(disabled, explicit_disabled);
    checkSameResourceAudit(disabled_audit, explicit_disabled_audit);
    checkBitwiseEqual(disabled, enabled);
    checkApproxTrace(enabled, cpu, 5e-13);
    checkApproxConvergenceChecks(enabled.convergence_checks,
                                 cpu.convergence_checks, 5e-13);
    checkPeriodicConvergenceDiagnosticsDelta(disabled, enabled, 2);
    checkPeriodicConvergenceAuditDelta(disabled_audit, enabled_audit, 2);
    CHECK(disabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(enabled.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(disabled_audit.graph_launch_calls == expected_launches[selection]);
    CHECK(enabled_audit.graph_launch_calls == expected_launches[selection]);
  }
}

TEST_CASE("hi-pdlp-recorded-periodic-convergence-failure-is-transactional",
          "[pdlp][hipdlp][phase7o][convergence][cuda][graph][replay]"
          "[errors][transactional]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  hipdlp::algorithm::FixedTraceOptions options = makeOptions(4);
  options.backend = hipdlp::backend::Kind::kCuda;
  options.execution_mode =
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
  enablePeriodicConvergence(options, 2);
  hipdlp::backend::cuda::testing::resetResourceAudit();
  hipdlp::backend::cuda::testing::setFailurePoint(
      hipdlp::backend::cuda::testing::FailurePoint::
          kGraphAsynchronousExecution);
  const hipdlp::algorithm::FixedTraceResult failed =
      runTrace(makeOneStepLp(), options);
  hipdlp::backend::cuda::testing::setFailurePoint(
      hipdlp::backend::cuda::testing::FailurePoint::kNone);
  checkTransactionalFailure(failed);
  CHECK(failed.error.code == hipdlp::api::ErrorCode::kBackendFailure);
  CHECK(failed.diagnostics.transfer_calls_after_iterations == 0);
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
  REQUIRE(recovered.convergence_checks.size() == 2);
}
#endif

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-uses-anchor-distances",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][exact]") {
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), makePrimalWeightOptions(3));
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 4);
  REQUIRE(result.convergence_checks.size() == 2);
  REQUIRE(result.fixed_point_restart_checks.size() == 2);
  REQUIRE(result.primal_weight_updates.size() == 1);
  checkPrimalWeightUpdatesFollowRestarts(result);

  const hipdlp::algorithm::FixedTracePrimalWeightUpdate& update =
      result.primal_weight_updates[0];
  CHECK(update.iteration == 2);
  CHECK(update.primal_anchor_distance == kAnchorPrimalDistanceAtTwo);
  CHECK(update.dual_anchor_distance == kAnchorDualDistanceAtTwo);
  // The transition then makes both anchors the candidate, so the same
  // distances measured afterwards would be zero.
  checkExactPoint(result.points[2], {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo}, {kCandidatePrimalAtTwo},
                  {kCandidateDualAtTwo});
  CHECK(update.primal_anchor_distance != 0.0);
  CHECK(update.dual_anchor_distance != 0.0);

  // The shared boundary records are the ones the policy consumed.
  CHECK(result.convergence_checks[0].iteration == 2);
  CHECK(result.convergence_checks[0].stopping.metrics.relative_primal ==
        kRelativePrimalAtTwo);
  CHECK(result.convergence_checks[0].stopping.metrics.relative_dual ==
        kBoundedRelativeDualAtTwo);
  CHECK(result.fixed_point_restart_checks[0].iteration == 2);
  CHECK(result.fixed_point_restart_checks[0].error.fixed_point_error ==
        kErrorAtTwo);
  CHECK(result.fixed_point_restart_checks[0]
            .error.diagnostics.device_to_host_scalar_readbacks == 3);
  CHECK(result.fixed_point_restart_checks[0].evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  CHECK(result.fixed_point_restart_checks[1].iteration == 3);
  CHECK(result.fixed_point_restart_checks[1].evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kContinue);

  const hipdlp::algorithm::PrimalWeightUpdateResult expected =
      expectedRestartPrimalWeightAtTwo(1.0, kBoundedRelativeDualAtTwo);
  REQUIRE(expected.ok());
  REQUIRE(update.result.ok());
  checkBitwiseEqualPrimalWeightResult(update.result, expected);
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-controls-and-falls-back",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][policy][exact]") {
  const hipdlp::algorithm::FixedTraceResult bounded =
      runTrace(makeOneStepLp(), makePrimalWeightOptions(3));
  const hipdlp::algorithm::FixedTraceResult free_column =
      runTrace(makeFreeStepLp(), makePrimalWeightOptions(3));
  INFO(bounded.error.message);
  INFO(free_column.error.message);
  REQUIRE(bounded.ok());
  REQUIRE(free_column.ok());
  REQUIRE(bounded.primal_weight_updates.size() == 1);
  REQUIRE(free_column.primal_weight_updates.size() == 1);
  REQUIRE(bounded.convergence_checks.size() == 2);
  REQUIRE(free_column.convergence_checks.size() == 2);

  // Both fixtures share the boundary candidate, so both share the distances.
  CHECK(free_column.primal_weight_updates[0].primal_anchor_distance ==
        kAnchorPrimalDistanceAtTwo);
  CHECK(free_column.primal_weight_updates[0].dual_anchor_distance ==
        kAnchorDualDistanceAtTwo);
  CHECK(free_column.fixed_point_restart_checks[0].error.fixed_point_error ==
        kErrorAtTwo);

  // A zero relative dual residual makes the ratio inadmissible and reverts to
  // the best weight; the free column keeps the ratio inside its window.
  CHECK(bounded.convergence_checks[0].stopping.metrics.relative_dual ==
        kBoundedRelativeDualAtTwo);
  CHECK(bounded.primal_weight_updates[0].result.kind ==
        hipdlp::algorithm::PrimalWeightUpdateKind::kRevertedToBest);
  CHECK(bounded.primal_weight_updates[0].result.next.primal_weight == 1.0);
  CHECK(bounded.primal_weight_updates[0].result.next.error_sum == 0.0);

  CHECK(free_column.convergence_checks[0].raw.dual_residual == 21.0 / 16.0);
  CHECK(free_column.convergence_checks[0].stopping.metrics.relative_primal ==
        kRelativePrimalAtTwo);
  CHECK(free_column.convergence_checks[0].stopping.metrics.relative_dual ==
        kFreeRelativeDualAtTwo);
  CHECK(free_column.primal_weight_updates[0].result.kind ==
        hipdlp::algorithm::PrimalWeightUpdateKind::kControlled);
  CHECK(free_column.primal_weight_updates[0].result.next.primal_weight != 1.0);

  checkBitwiseEqualPrimalWeightResult(
      bounded.primal_weight_updates[0].result,
      expectedRestartPrimalWeightAtTwo(1.0, kBoundedRelativeDualAtTwo));
  checkBitwiseEqualPrimalWeightResult(
      free_column.primal_weight_updates[0].result,
      expectedRestartPrimalWeightAtTwo(1.0, kFreeRelativeDualAtTwo));
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-preserves-eta",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][eta]") {
  const HighsLp fixtures[] = {makeOneStepLp(), makeFreeStepLp()};
  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(selection);
    hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(3);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(fixtures[selection], options);
    INFO(result.error.message);
    REQUIRE(result.ok());
    REQUIRE(result.primal_weight_updates.size() == 1);
    const hipdlp::algorithm::InitialStepSizes& steps =
        result.primal_weight_updates[0].result.step_sizes;
    const double eta =
        std::sqrt(options.primal_step_size * options.dual_step_size);
    CHECK(steps.base_step == eta);
    CHECK(std::sqrt(steps.primal_step * steps.dual_step) ==
          Approx(eta).margin(1e-15));
    CHECK(steps.primal_step == Approx(eta / steps.primal_weight).margin(1e-15));
    CHECK(steps.dual_step == Approx(eta * steps.primal_weight).margin(1e-15));
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-installs-later-steps",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][steps]") {
  hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(4);
  options.fixed_point_restart.parameters.artificial_threshold = 0.6;
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.points.size() == 5);
  REQUIRE(result.primal_weight_updates.size() == 1);
  checkPrimalWeightUpdatesFollowRestarts(result);
  CHECK(result.primal_weight_updates[0].iteration == 2);
  const hipdlp::algorithm::InitialStepSizes installed =
      result.primal_weight_updates[0].result.step_sizes;
  CHECK(installed.primal_step != options.primal_step_size);
  CHECK(installed.dual_step != options.dual_step_size);

  // Independent oracle: the same run with the installed steps supplied as a
  // caller sequence that starts only after the boundary.
  hipdlp::algorithm::FixedTraceOptions oracle_options = options;
  oracle_options.primal_weight.enabled = false;
  oracle_options.primal_step_sizes.assign(2, options.primal_step_size);
  oracle_options.primal_step_sizes.resize(4, installed.primal_step);
  oracle_options.dual_step_sizes.assign(2, options.dual_step_size);
  oracle_options.dual_step_sizes.resize(4, installed.dual_step);
  const hipdlp::algorithm::FixedTraceResult oracle =
      runTrace(makeOneStepLp(), oracle_options);
  INFO(oracle.error.message);
  REQUIRE(oracle.ok());
  CHECK(oracle.primal_weight_updates.empty());
  checkBitwiseEqualTraceOutputs(result, oracle);
  checkBitwiseEqualFixedPointError(result.fixed_point_error,
                                   oracle.fixed_point_error);
  checkBitwiseEqualFixedPointRestartChecks(result.fixed_point_restart_checks,
                                           oracle.fixed_point_restart_checks);
  checkBitwiseEqualConvergenceChecks(result.convergence_checks,
                                     oracle.convergence_checks);

  // Only iterations after the boundary move.
  hipdlp::algorithm::FixedTraceOptions disabled_options = options;
  disabled_options.primal_weight.enabled = false;
  const hipdlp::algorithm::FixedTraceResult disabled =
      runTrace(makeOneStepLp(), disabled_options);
  INFO(disabled.error.message);
  REQUIRE(disabled.ok());
  REQUIRE(disabled.convergence_checks.size() == 2);
  REQUIRE(disabled.fixed_point_restart_checks.size() == 2);
  for (std::size_t point = 0; point <= 2; ++point) {
    CAPTURE(point);
    CHECK(bitwiseEqual(disabled.points[point].primal,
                       result.points[point].primal));
    CHECK(bitwiseEqual(disabled.points[point].dual, result.points[point].dual));
  }
  CHECK_FALSE(bitwiseEqual(disabled.points[3].primal, result.points[3].primal));
  CHECK_FALSE(bitwiseEqual(disabled.points[4].primal, result.points[4].primal));
  CHECK(bitwiseEqual(disabled.convergence_checks[0].raw.primal_residual,
                     result.convergence_checks[0].raw.primal_residual));
  CHECK(bitwiseEqual(
      disabled.fixed_point_restart_checks[0].error.fixed_point_error,
      result.fixed_point_restart_checks[0].error.fixed_point_error));
  CHECK(result.convergence_checks[1].iteration == 4);
  CHECK(result.fixed_point_restart_checks[1].iteration == 4);
  checkSameDiagnostics(disabled, result);
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-records-final-boundary",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][final]") {
  hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(2);
  const hipdlp::algorithm::FixedTraceResult result =
      runTrace(makeOneStepLp(), options);
  hipdlp::algorithm::FixedTraceOptions disabled_options = options;
  disabled_options.primal_weight.enabled = false;
  const hipdlp::algorithm::FixedTraceResult disabled =
      runTrace(makeOneStepLp(), disabled_options);
  INFO(result.error.message);
  INFO(disabled.error.message);
  REQUIRE(result.ok());
  REQUIRE(disabled.ok());
  REQUIRE(result.fixed_point_restart_checks.size() == 1);
  CHECK(result.fixed_point_restart_checks[0].iteration == 2);
  CHECK(result.fixed_point_restart_checks[0].evaluation.decision ==
        hipdlp::algorithm::RestartDecision::kRestart);
  REQUIRE(result.primal_weight_updates.size() == 1);
  CHECK(result.primal_weight_updates[0].iteration == 2);
  CHECK(result.primal_weight_updates[0].primal_anchor_distance ==
        kAnchorPrimalDistanceAtTwo);
  CHECK(result.primal_weight_updates[0].dual_anchor_distance ==
        kAnchorDualDistanceAtTwo);
  REQUIRE(result.primal_weight_updates[0].result.ok());

  // No later iteration consumes the returned steps, so nothing else moves.
  CHECK(disabled.primal_weight_updates.empty());
  checkBitwiseEqualTraceOutputs(result, disabled);
  checkBitwiseEqualFixedPointError(result.fixed_point_error,
                                   disabled.fixed_point_error);
  checkBitwiseEqualFixedPointRestart(result.fixed_point_restart,
                                     disabled.fixed_point_restart);
  checkBitwiseEqualFixedPointRestartChecks(result.fixed_point_restart_checks,
                                           disabled.fixed_point_restart_checks);
  checkBitwiseEqualConvergenceChecks(result.convergence_checks,
                                     disabled.convergence_checks);
  checkSameDiagnostics(result, disabled);
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-skips-continue",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][disabled]"
          "[bitwise]") {
  const std::vector<RestartExecutionSelection> selections =
      restartExecutionSelections();
  for (const RestartExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions disabled_options =
        makePrimalWeightOptions(4);
    disabled_options.backend = selection.backend;
    disabled_options.execution_mode = selection.mode;
    disabled_options.primal_weight.enabled = false;
    disabled_options.fixed_point_restart.parameters.check_interval = 100;
    disabled_options.fixed_point_restart.parameters.artificial_threshold =
        100.0;
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeOneStepLp(), disabled_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enablePrimalWeight(enabled_options, 1.0);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    hipdlp::backend::cuda::testing::resetResourceAudit();
#endif
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeOneStepLp(), enabled_options);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();
#endif

    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    REQUIRE(enabled.fixed_point_restart_checks.size() == 2);
    for (std::size_t check = 0; check < 2; ++check) {
      CAPTURE(check);
      CHECK(enabled.fixed_point_restart_checks[check].evaluation.decision ==
            hipdlp::algorithm::RestartDecision::kContinue);
    }
    CHECK(enabled.primal_weight_updates.empty());
    checkBitwiseEqual(disabled, enabled);
    checkBitwiseEqualConvergenceChecks(disabled.convergence_checks,
                                       enabled.convergence_checks);
    checkSameDiagnostics(disabled, enabled);
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    if (selection.backend == hipdlp::backend::Kind::kCuda)
      checkPrimalWeightAuditDelta(disabled_audit, enabled_audit, 0);
#endif
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-is-deterministic",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][determinism]") {
  const std::vector<RestartExecutionSelection> selections =
      restartExecutionSelections();
  for (const RestartExecutionSelection& selection : selections) {
    CAPTURE(static_cast<int>(selection.backend),
            static_cast<int>(selection.mode));
    hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(4);
    options.backend = selection.backend;
    options.execution_mode = selection.mode;
    options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    const hipdlp::algorithm::FixedTraceResult first =
        runTrace(makeFreeStepLp(), options);
    const hipdlp::algorithm::FixedTraceResult second =
        runTrace(makeFreeStepLp(), options);
    INFO(first.error.message);
    INFO(second.error.message);
    REQUIRE(first.ok());
    REQUIRE(second.ok());
    REQUIRE(first.primal_weight_updates.size() == 1);
    checkPrimalWeightUpdatesFollowRestarts(first);
    checkBitwiseEqual(first, second);
    checkBitwiseEqualConvergenceChecks(first.convergence_checks,
                                       second.convergence_checks);
    checkSameDiagnostics(first, second);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-validation",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][errors]") {
  SECTION("updates without the fixed-point restart trigger are rejected") {
    hipdlp::algorithm::FixedTraceOptions options = makeOptions(4);
    enablePeriodicConvergence(options, 2);
    enablePrimalWeight(options, 1.0);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "fixed trace primal-weight updates require the fixed-point restart "
          "trigger");
    CHECK(result.diagnostics.resource_creation_calls_before_iterations == 0);
    CHECK(result.diagnostics.transfer_calls_before_iterations == 0);
  }

  SECTION("an unobserved restart boundary is rejected") {
    const HighsInt convergence_cadences[] = {-1, 0, 4};
    const HighsInt iteration_counts[] = {4, 4, 8};
    const bool periodic_enabled[] = {false, true, true};
    for (std::size_t selection = 0; selection < 3; ++selection) {
      CAPTURE(selection);
      hipdlp::algorithm::FixedTraceOptions options =
          makeRestartAtTwoOptions(iteration_counts[selection]);
      if (periodic_enabled[selection])
        enablePeriodicConvergence(options, convergence_cadences[selection]);
      enablePrimalWeight(options, 1.0);
      const hipdlp::algorithm::FixedTraceResult result =
          runTrace(makeOneStepLp(), options);
      checkTransactionalFailure(result);
      CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
      CHECK(result.error.message ==
            "fixed trace primal-weight updates require a convergence "
            "observation at every restart boundary");
      CHECK(result.diagnostics.resource_creation_calls_before_iterations == 0);
      CHECK(result.diagnostics.transfer_calls_before_iterations == 0);
    }
  }

  SECTION("caller-provided step-size sequences are rejected") {
    hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(4);
    options.primal_step_sizes.assign(4, options.primal_step_size);
    options.dual_step_sizes.assign(4, options.dual_step_size);
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "fixed trace primal-weight updates cannot combine with "
          "caller-provided step-size sequences");
    CHECK(result.diagnostics.resource_creation_calls_before_iterations == 0);
    CHECK(result.diagnostics.transfer_calls_before_iterations == 0);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-is-transactional",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][errors]"
          "[transactional]") {
  SECTION("an invalid initial state fails at the boundary") {
    hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(3);
    options.primal_weight.previous = hipdlp::algorithm::PrimalWeightState();
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidArgument);
    CHECK(result.error.message == "HiPDLP primal-weight state is invalid");
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
  }

  SECTION("invalid PID parameters fail at the boundary") {
    hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(3);
    options.primal_weight.parameters.k_p = -1.0;
    const hipdlp::algorithm::FixedTraceResult result =
        runTrace(makeOneStepLp(), options);
    checkTransactionalFailure(result);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "HiPDLP primal-weight PID parameters are invalid");
    CHECK(result.diagnostics.transfer_calls_after_iterations == 0);
  }

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
  SECTION("injected backend failures expose no partial result") {
    hipdlp::backend::FactoryResult cuda_probe =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
    if (!cuda_probe.ok()) {
      CHECK(cuda_probe.error.code ==
            hipdlp::api::ErrorCode::kBackendUnavailable);
      return;
    }
    cuda_probe.backend.reset();

    hipdlp::algorithm::FixedTraceOptions options = makePrimalWeightOptions(4);
    options.backend = hipdlp::backend::Kind::kCuda;
    options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    hipdlp::backend::cuda::testing::setFailurePoint(
        hipdlp::backend::cuda::testing::FailurePoint::
            kGraphAsynchronousExecution);
    const hipdlp::algorithm::FixedTraceResult failed =
        runTrace(makeOneStepLp(), options);
    hipdlp::backend::cuda::testing::setFailurePoint(
        hipdlp::backend::cuda::testing::FailurePoint::kNone);
    checkTransactionalFailure(failed);
    CHECK(failed.error.code == hipdlp::api::ErrorCode::kBackendFailure);
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
    REQUIRE(recovered.primal_weight_updates.size() == 1);
  }
#endif
}

#ifdef HIGHS_HIPDLP_CUDA_DETECTED
TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-backend-agreement",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][cuda][graph]"
          "[replay][parity]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  const HighsInt iteration_counts[] = {2, 4};
  const std::size_t expected_launches[] = {1, 2};
  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(iteration_counts[selection]);
    hipdlp::algorithm::FixedTraceOptions cpu_options =
        makePrimalWeightOptions(iteration_counts[selection]);
    cpu_options.fixed_point_restart.parameters.artificial_threshold = 0.6;
    const hipdlp::algorithm::FixedTraceResult cpu =
        runTrace(makeFreeStepLp(), cpu_options);
    INFO(cpu.error.message);
    REQUIRE(cpu.ok());
    REQUIRE(cpu.primal_weight_updates.size() == 1);

    hipdlp::algorithm::FixedTraceOptions direct_options = cpu_options;
    direct_options.backend = hipdlp::backend::Kind::kCuda;
    const hipdlp::algorithm::FixedTraceResult direct =
        runTrace(makeFreeStepLp(), direct_options);

    hipdlp::algorithm::FixedTraceOptions recorded_options = direct_options;
    recorded_options.execution_mode =
        hipdlp::algorithm::FixedTraceExecutionMode::kRecorded;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult recorded =
        runTrace(makeFreeStepLp(), recorded_options);
    const hipdlp::backend::cuda::testing::ResourceAudit recorded_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    INFO(direct.error.message);
    INFO(recorded.error.message);
    REQUIRE(direct.ok());
    REQUIRE(recorded.ok());
    checkSameRestartDecisions(cpu, direct);
    checkSameRestartDecisions(cpu, recorded);
    checkPrimalWeightUpdatesFollowRestarts(direct);
    checkPrimalWeightUpdatesFollowRestarts(recorded);
    checkApproxPrimalWeightUpdates(direct.primal_weight_updates,
                                   cpu.primal_weight_updates, 5e-13);
    checkApproxPrimalWeightUpdates(recorded.primal_weight_updates,
                                   cpu.primal_weight_updates, 5e-13);
    checkApproxTrace(direct, cpu, 5e-13);
    checkApproxTrace(recorded, cpu, 5e-13);
    checkApproxConvergenceChecks(recorded.convergence_checks,
                                 cpu.convergence_checks, 5e-13);
    CHECK(recorded.diagnostics.command_buffer.replay_calls ==
          expected_launches[selection]);
    CHECK(recorded_audit.graph_launch_calls == expected_launches[selection]);
  }
}

TEST_CASE("hi-pdlp-fixed-trace-restart-primal-weight-adds-only-two-reductions",
          "[pdlp][hipdlp][phase7p][restart][primal-weight][cuda][resources]") {
  hipdlp::backend::FactoryResult cuda_probe =
      hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
  if (!cuda_probe.ok()) {
    CHECK(cuda_probe.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK_FALSE(cuda_probe.error.message.empty());
    return;
  }
  cuda_probe.backend.reset();

  const hipdlp::algorithm::FixedTraceExecutionMode modes[] = {
      hipdlp::algorithm::FixedTraceExecutionMode::kDirect,
      hipdlp::algorithm::FixedTraceExecutionMode::kRecorded};
  for (std::size_t selection = 0; selection < 2; ++selection) {
    CAPTURE(selection);
    // A final-boundary restart has no later iteration, so the enabled and
    // disabled runs differ only by the two reductions.
    hipdlp::algorithm::FixedTraceOptions disabled_options =
        makePrimalWeightOptions(2);
    disabled_options.backend = hipdlp::backend::Kind::kCuda;
    disabled_options.execution_mode = modes[selection];
    disabled_options.primal_weight.enabled = false;
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult disabled =
        runTrace(makeFreeStepLp(), disabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit disabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    hipdlp::algorithm::FixedTraceOptions enabled_options = disabled_options;
    enablePrimalWeight(enabled_options, 1.0);
    hipdlp::backend::cuda::testing::resetResourceAudit();
    const hipdlp::algorithm::FixedTraceResult enabled =
        runTrace(makeFreeStepLp(), enabled_options);
    const hipdlp::backend::cuda::testing::ResourceAudit enabled_audit =
        hipdlp::backend::cuda::testing::resourceAudit();

    INFO(disabled.error.message);
    INFO(enabled.error.message);
    REQUIRE(disabled.ok());
    REQUIRE(enabled.ok());
    CHECK(disabled.primal_weight_updates.empty());
    REQUIRE(enabled.primal_weight_updates.size() == 1);
    CHECK(enabled.primal_weight_updates[0].result.kind ==
          hipdlp::algorithm::PrimalWeightUpdateKind::kControlled);
    checkBitwiseEqualTraceOutputs(disabled, enabled);
    checkSameDiagnostics(disabled, enabled);
    checkPrimalWeightAuditDelta(disabled_audit, enabled_audit, 1);
  }
}
#endif

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