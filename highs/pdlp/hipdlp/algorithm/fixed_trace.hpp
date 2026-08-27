#ifndef PDLP_HIPDLP_ALGORITHM_FIXED_TRACE_HPP_
#define PDLP_HIPDLP_ALGORITHM_FIXED_TRACE_HPP_

#include <cstddef>
#include <vector>

#include "pdlp/hipdlp/algorithm/averaging_policy.hpp"
#include "pdlp/hipdlp/algorithm/convergence.hpp"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/command_buffer.hpp"
#include "pdlp/hipdlp/backend/types.hpp"
#include "pdlp/hipdlp/model/problem.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

namespace hipdlp {
namespace algorithm {

enum class FixedTraceExecutionMode { kDirect, kRecorded };

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
  ConvergenceTolerances convergence_tolerances;

  FixedTraceOptions();
};

struct FixedTraceConvergenceObservation {
  bool observed;
  RawConvergenceMetrics raw;
  StoppingResult stopping;

  FixedTraceConvergenceObservation();
};

struct FixedTracePoint {
  std::vector<double> primal;
  std::vector<double> dual;
  std::vector<double> ax;
  std::vector<double> aty;
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
  FixedTraceDiagnostics diagnostics;
  model::ScalingMetadata scaling;

  bool ok() const;
};

FixedTraceResult runFixedTrace(const model::OriginalModel& original_model,
                               const FixedTraceOptions& options);

}  // namespace algorithm
}  // namespace hipdlp

#endif