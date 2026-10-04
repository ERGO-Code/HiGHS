// Test-only iteration replay and oracle. Never linked into the solver library.
#ifndef CHECK_HIPDLP_FIXED_TRACE_HPP_
#define CHECK_HIPDLP_FIXED_TRACE_HPP_

#include <cstddef>
#include <vector>

#include "pdlp/hipdlp/algorithm/averaging_policy.hpp"
#include "pdlp/hipdlp/algorithm/convergence.hpp"
#include "pdlp/hipdlp/algorithm/fixed_point_error.hpp"
#include "pdlp/hipdlp/algorithm/restart_policy.hpp"
#include "pdlp/hipdlp/algorithm/step_size_policy.hpp"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/command_buffer.hpp"
#include "pdlp/hipdlp/backend/types.hpp"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

namespace hipdlp {
namespace algorithm {

enum class FixedTraceExecutionMode { kDirect, kRecorded };

struct FixedTraceConvergenceOptions {
  bool periodic_enabled;
  HighsInt check_cadence;

  FixedTraceConvergenceOptions();
};

struct FixedTraceFixedPointErrorOptions {
  bool enabled;
  double omega;
  double eta;

  FixedTraceFixedPointErrorOptions();
};

struct FixedTraceFixedPointRestartOptions {
  bool enabled;
  // 0 observes only the final iteration; a positive cadence also observes every
  // multiple of it.
  HighsInt check_cadence;
  FixedPointRestartState previous;
  FixedPointRestartParameters parameters;

  FixedTraceFixedPointRestartOptions();
};

struct FixedTracePrimalWeightOptions {
  bool enabled;
  PrimalWeightState previous;
  PrimalWeightParameters parameters;

  FixedTracePrimalWeightOptions();
};

struct FixedTraceOptions {
  backend::Kind backend;
  FixedTraceExecutionMode execution_mode;
  model::ScalingOptions scaling;
  double primal_step_size;
  double dual_step_size;
  std::vector<double> primal_step_sizes;
  std::vector<double> dual_step_sizes;
  HighsInt iteration_count;
  std::vector<double> scaled_primal_start;
  std::vector<double> scaled_dual_start;
  FixedTraceConvergenceOptions convergence;
  ConvergenceTolerances convergence_tolerances;
  FixedTraceFixedPointErrorOptions fixed_point_error;
  FixedTraceFixedPointRestartOptions fixed_point_restart;
  FixedTracePrimalWeightOptions primal_weight;

  FixedTraceOptions();
};

struct FixedTraceConvergenceObservation {
  bool observed;
  RawConvergenceMetrics raw;
  StoppingResult stopping;

  FixedTraceConvergenceObservation();
};

struct FixedTraceConvergenceCheck {
  HighsInt iteration;
  RawConvergenceMetrics raw;
  StoppingResult stopping;

  FixedTraceConvergenceCheck();
};

struct FixedTraceFixedPointErrorObservation {
  bool observed;
  FixedPointErrorResult evaluation;

  FixedTraceFixedPointErrorObservation();
};

struct FixedTraceFixedPointRestartObservation {
  bool observed;
  FixedPointRestartResult evaluation;

  FixedTraceFixedPointRestartObservation();
};

// A kRestart decision carries resetFixedPointRestart(error.fixed_point_error)
// into the following check, so the reset is recoverable from this record.
struct FixedTraceFixedPointRestartCheck {
  HighsInt iteration;
  FixedPointErrorResult error;
  FixedPointRestartResult evaluation;

  FixedTraceFixedPointRestartCheck();
};

struct FixedTracePoint {
  std::vector<double> primal;
  std::vector<double> dual;
  std::vector<double> ax;
  std::vector<double> aty;
};

// Recorded only where the shared boundary decided kRestart; the distances are
// measured before either anchor is overwritten.
struct FixedTracePrimalWeightUpdate {
  HighsInt iteration;
  double primal_anchor_distance;
  double dual_anchor_distance;
  PrimalWeightUpdateResult result;

  FixedTracePrimalWeightUpdate();
};
struct FixedTraceDiagnostics {
  std::size_t resource_creation_calls_before_iterations;
  std::size_t resource_creation_calls_during_iterations;
  std::size_t transfer_calls_before_iterations;
  std::size_t transfer_calls_during_iterations;
  std::size_t transfer_calls_after_iterations;
  backend::CommandBufferStatistics command_buffer;

  FixedTraceDiagnostics();
};

struct FixedTraceResult {
  api::Error error;
  std::vector<FixedTracePoint> points;
  std::vector<double> averaged_primal;
  std::vector<double> averaged_dual;
  AveragingState averaging;
  FixedTraceConvergenceObservation current_convergence;
  FixedTraceConvergenceObservation averaged_convergence;
  std::vector<FixedTraceConvergenceCheck> convergence_checks;
  FixedTraceFixedPointErrorObservation fixed_point_error;
  FixedTraceFixedPointRestartObservation fixed_point_restart;
  std::vector<FixedTraceFixedPointRestartCheck> fixed_point_restart_checks;
  std::vector<FixedTracePrimalWeightUpdate> primal_weight_updates;
  FixedTraceDiagnostics diagnostics;
  model::ScalingMetadata scaling;

  bool ok() const;
};

FixedTraceResult runFixedTrace(const model::OriginalModel& original_model,
                               const FixedTraceOptions& options);

}  // namespace algorithm
}  // namespace hipdlp

#endif
