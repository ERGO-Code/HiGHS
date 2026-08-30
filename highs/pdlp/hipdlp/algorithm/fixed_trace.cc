#include "pdlp/hipdlp/algorithm/fixed_trace.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <new>
#include <utility>

#include "pdlp/hipdlp/algorithm/halpern_pdhg_step.hpp"
#include "pdlp/hipdlp/backend/backend.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"

namespace hipdlp {
namespace algorithm {

namespace {

const double kDefaultFixedTraceTolerance = 1e-4;

// Recorded traces enqueue this many iterations per command buffer.
const HighsInt kRecordedBlockIterations = 2;

struct TraceHandles {
  backend::VectorHandle primal;
  backend::VectorHandle dual;
  backend::VectorHandle ax;
  backend::VectorHandle aty;
};

// The retired next-iterate slots, which alternate as the Halpern step swaps.
struct IterateSlots {
  backend::VectorHandle* current_primal;
  backend::VectorHandle* current_dual;
  backend::VectorHandle* next_primal;
  backend::VectorHandle* next_dual;
};

bool finiteVector(const std::vector<double>& values) {
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (!std::isfinite(values[index])) return false;
  }
  return true;
}

FixedTraceResult failure(const api::Error& error,
                         const FixedTraceDiagnostics& diagnostics) {
  FixedTraceResult result;
  result.error = error;
  result.diagnostics = diagnostics;
  return result;
}

api::Error createVector(backend::Backend& backend, const HighsInt size,
                        backend::VectorHandle& vector,
                        FixedTraceDiagnostics& diagnostics) {
  ++diagnostics.resource_creation_calls_before_iterations;
  return backend.createVector(size, vector);
}

api::Error recordPoint(backend::Backend& backend,
                       const backend::VectorHandle& primal,
                       const backend::VectorHandle& dual,
                       const backend::VectorHandle& ax,
                       const backend::VectorHandle& aty, TraceHandles& point) {
  api::Error error = backend.copy(primal, point.primal);
  if (!error.ok()) return error;
  error = backend.copy(dual, point.dual);
  if (!error.ok()) return error;
  error = backend.copy(ax, point.ax);
  if (!error.ok()) return error;
  return backend.copy(aty, point.aty);
}

api::Error downloadPoint(backend::Backend& backend, const TraceHandles& handles,
                         FixedTracePoint& point,
                         FixedTraceDiagnostics& diagnostics) {
  api::Error error = backend.downloadVector(handles.primal, point.primal);
  ++diagnostics.transfer_calls_after_iterations;
  if (!error.ok()) return error;
  error = backend.downloadVector(handles.dual, point.dual);
  ++diagnostics.transfer_calls_after_iterations;
  if (!error.ok()) return error;
  error = backend.downloadVector(handles.ax, point.ax);
  ++diagnostics.transfer_calls_after_iterations;
  if (!error.ok()) return error;
  error = backend.downloadVector(handles.aty, point.aty);
  ++diagnostics.transfer_calls_after_iterations;
  return error;
}

api::Error requireTraceCapabilities(backend::Backend& backend) {
  const backend::Capability required[] = {
      backend::Capability::kSparseMatrixVector,
      backend::Capability::kFusedPrimalDualUpdate,
      backend::Capability::kFusedAverageUpdate,
      backend::Capability::kVectorOperations,
      backend::Capability::kHostDownload};
  for (std::size_t index = 0; index < sizeof(required) / sizeof(required[0]);
       ++index) {
    const api::Error error = backend.requireCapability(required[index]);
    if (!error.ok()) return error;
  }
  return api::Error();
}

void updateScaledSumOfSquares(const double value, double& scale,
                              double& sum_of_squares) {
  const double magnitude = std::fabs(value);
  if (magnitude == 0.0) return;
  if (scale < magnitude) {
    const double ratio = scale / magnitude;
    sum_of_squares = 1.0 + sum_of_squares * ratio * ratio;
    scale = magnitude;
  } else {
    const double ratio = magnitude / scale;
    sum_of_squares += ratio * ratio;
  }
}

double scaledNorm(const double scale, const double sum_of_squares) {
  return scale == 0.0 ? 0.0 : scale * std::sqrt(sum_of_squares);
}

// Raw metrics stay in the transformed problem's unscaled units, which is the
// space measured by Problem::unscaledNorms().
api::Error computeRawMetrics(const model::Problem& problem,
                             const std::vector<double>& primal,
                             const std::vector<double>& dual,
                             const std::vector<double>& ax,
                             const std::vector<double>& aty,
                             RawConvergenceMetrics& raw) {
  const model::ScalingMetadata& scaling = problem.scaling();
  const std::vector<double>& row_lower = problem.rowLower();
  const std::vector<double>& col_lower = problem.colLower();
  const std::vector<double>& col_upper = problem.colUpper();
  const std::vector<double>& costs = problem.costs();
  const std::vector<HighsBool>& equality_rows = problem.equalityRows();

  double primal_objective = problem.objectiveOffset();
  double dual_objective = problem.objectiveOffset();
  double primal_scale = 0.0;
  double primal_sum_of_squares = 1.0;
  for (HighsInt row = 0; row < problem.numRows(); ++row) {
    double violation = ax[row] - row_lower[row];
    if (!equality_rows[row]) violation = std::min(0.0, violation);
    updateScaledSumOfSquares(violation / scaling.row_scale[row], primal_scale,
                             primal_sum_of_squares);
    dual_objective += row_lower[row] * dual[row];
  }

  double dual_scale = 0.0;
  double dual_sum_of_squares = 1.0;
  for (HighsInt col = 0; col < problem.numCols(); ++col) {
    const double reduced_cost = costs[col] - aty[col];
    const bool has_lower = col_lower[col] > -kHighsInf;
    const bool has_upper = col_upper[col] < kHighsInf;
    const double lower_slack = has_lower ? std::max(0.0, reduced_cost) : 0.0;
    const double upper_slack = has_upper ? std::max(0.0, -reduced_cost) : 0.0;
    updateScaledSumOfSquares(
        (reduced_cost - lower_slack + upper_slack) / scaling.col_scale[col],
        dual_scale, dual_sum_of_squares);
    primal_objective += costs[col] * primal[col];
    if (has_lower) dual_objective += col_lower[col] * lower_slack;
    if (has_upper) dual_objective -= col_upper[col] * upper_slack;
  }

  raw.primal_residual = scaledNorm(primal_scale, primal_sum_of_squares);
  raw.dual_residual = scaledNorm(dual_scale, dual_sum_of_squares);
  raw.primal_objective = primal_objective;
  raw.dual_objective = dual_objective;
  raw.objective_gap = std::fabs(primal_objective - dual_objective);
  raw.unscaled_rhs_norm = problem.unscaledNorms().row_bounds_l2;
  raw.unscaled_cost_norm = problem.unscaledNorms().objective_l2;
  if (!std::isfinite(raw.primal_residual) ||
      !std::isfinite(raw.dual_residual) ||
      !std::isfinite(raw.primal_objective) ||
      !std::isfinite(raw.dual_objective) || !std::isfinite(raw.objective_gap)) {
    return api::Error(api::ErrorCode::kInvalidArgument,
                      "fixed trace convergence metrics are non-finite");
  }
  return api::Error();
}

api::Error observeConvergence(
    const model::Problem& problem, const FixedTraceOptions& options,
    const std::vector<double>& primal, const std::vector<double>& dual,
    const std::vector<double>& ax, const std::vector<double>& aty,
    const HighsInt iteration, FixedTraceConvergenceObservation& observation) {
  RawConvergenceMetrics raw;
  const api::Error error =
      computeRawMetrics(problem, primal, dual, ax, aty, raw);
  if (!error.ok()) return error;
  StoppingInput input;
  input.raw = raw;
  input.tolerances = options.convergence_tolerances;
  input.iteration = iteration;
  const StoppingResult stopping = evaluateStopping(input);
  if (!stopping.ok()) return stopping.error;
  observation.observed = true;
  observation.raw = raw;
  observation.stopping = stopping;
  return api::Error();
}

bool periodicConvergenceCheck(const FixedTraceOptions& options,
                              const HighsInt iteration) {
  if (!options.convergence.periodic_enabled || iteration <= 0) return false;
  if (iteration == options.iteration_count) return true;
  const HighsInt cadence = options.convergence.check_cadence;
  return cadence > 0 && iteration % cadence == 0;
}

api::Error downloadConvergenceVector(backend::Backend& backend,
                                     const backend::VectorHandle& vector,
                                     std::vector<double>& values,
                                     FixedTraceDiagnostics& diagnostics) {
  const api::Error error = backend.downloadVector(vector, values);
  ++diagnostics.transfer_calls_during_iterations;
  return error;
}

api::Error appendConvergenceCheck(
    backend::Backend& backend, const model::Problem& problem,
    const FixedTraceOptions& options, const HalpernPdhgState& state,
    backend::VectorHandle& candidate_ax, backend::VectorHandle& candidate_aty,
    const HighsInt iteration, FixedTraceDiagnostics& diagnostics,
    std::vector<FixedTraceConvergenceCheck>& checks) {
  api::Error error =
      backend.multiply(*state.problem, *state.primal_candidate, candidate_ax);
  if (!error.ok()) return error;
  error = backend.multiplyTranspose(*state.problem, *state.dual_candidate,
                                    candidate_aty);
  if (!error.ok()) return error;

  std::vector<double> primal;
  std::vector<double> dual;
  std::vector<double> ax;
  std::vector<double> aty;
  error = downloadConvergenceVector(backend, *state.primal_candidate, primal,
                                    diagnostics);
  if (!error.ok()) return error;
  error = downloadConvergenceVector(backend, *state.dual_candidate, dual,
                                    diagnostics);
  if (!error.ok()) return error;
  error = downloadConvergenceVector(backend, candidate_ax, ax, diagnostics);
  if (!error.ok()) return error;
  error = downloadConvergenceVector(backend, candidate_aty, aty, diagnostics);
  if (!error.ok()) return error;

  FixedTraceConvergenceObservation observation;
  error = observeConvergence(problem, options, primal, dual, ax, aty, iteration,
                             observation);
  if (!error.ok()) return error;
  FixedTraceConvergenceCheck check;
  check.iteration = iteration;
  check.raw = observation.raw;
  check.stopping = observation.stopping;
  checks.push_back(check);
  return api::Error();
}

bool intermediateFixedPointCheck(const FixedTraceOptions& options,
                                 const HighsInt iteration) {
  const HighsInt cadence = options.fixed_point_restart.check_cadence;
  return options.fixed_point_restart.enabled && cadence > 0 && iteration > 0 &&
         iteration < options.iteration_count && iteration % cadence == 0;
}

FixedPointRestartState chainedFixedPointRestartState(
    const std::vector<FixedTraceFixedPointRestartCheck>& checks,
    const FixedPointRestartState& initial) {
  if (checks.empty()) return initial;
  const FixedTraceFixedPointRestartCheck& last = checks.back();
  if (last.evaluation.decision == RestartDecision::kRestart)
    return resetFixedPointRestart(last.error.fixed_point_error);
  return last.evaluation.next;
}

api::Error appendFixedPointRestartCheck(
    backend::Backend& backend, const FixedTraceOptions& options,
    const HalpernPdhgState& state, const IterateSlots& slots,
    backend::VectorHandle& aty_delta_scratch, const HighsInt iteration,
    std::vector<FixedTraceFixedPointRestartCheck>& checks) {
  FixedPointErrorInput error_input;
  error_input.next_primal = state.primal_candidate;
  error_input.next_dual = state.dual_candidate;
  error_input.reflected_primal = state.reflected_primal;
  error_input.reflected_dual = state.reflected_dual;
  error_input.omega = options.fixed_point_error.omega;
  error_input.eta = options.fixed_point_error.eta;
  error_input.delta_primal = slots.next_primal;
  error_input.delta_dual = slots.next_dual;
  error_input.aty_delta = &aty_delta_scratch;
  const FixedPointErrorResult error =
      evaluateFixedPointError(backend, *state.problem, error_input);
  if (!error.ok()) return error.error;

  FixedPointRestartInput restart_input;
  restart_input.previous = chainedFixedPointRestartState(
      checks, options.fixed_point_restart.previous);
  restart_input.current_error = error.fixed_point_error;
  // HalpernPdhgState::iteration is the resettable age; the trace owns the
  // monotone total separately.
  restart_input.halpern_iterations = state.iteration;
  restart_input.total_iterations = iteration;
  restart_input.parameters = options.fixed_point_restart.parameters;
  const FixedPointRestartResult evaluation =
      evaluateFixedPointRestart(restart_input);
  if (!evaluation.ok()) return evaluation.error;

  FixedTraceFixedPointRestartCheck check;
  check.iteration = iteration;
  check.error = error;
  check.evaluation = evaluation;
  checks.push_back(check);
  return api::Error();
}

// Models the legacy PDLPSolver::solve() restart body: both anchors and the
// current iterate take the projected candidates and the cached products are
// rebuilt, without touching either step size or the primal weight.
api::Error applyFixedPointRestart(backend::Backend& backend,
                                  HalpernPdhgState& state,
                                  backend::VectorHandle& primal_sum,
                                  backend::VectorHandle& dual_sum,
                                  AveragingState& averaging,
                                  double& average_normalization,
                                  TraceHandles& point) {
  api::Error error =
      backend.copy(*state.primal_candidate, *state.current_primal);
  if (!error.ok()) return error;
  error = backend.copy(*state.dual_candidate, *state.current_dual);
  if (!error.ok()) return error;
  error = backend.copy(*state.primal_candidate, *state.primal_anchor);
  if (!error.ok()) return error;
  error = backend.copy(*state.dual_candidate, *state.dual_anchor);
  if (!error.ok()) return error;
  error = backend.multiply(*state.problem, *state.current_primal,
                           *state.current_ax);
  if (!error.ok()) return error;
  error = backend.multiplyTranspose(*state.problem, *state.current_dual,
                                    *state.current_aty);
  if (!error.ok()) return error;
  error = backend.fill(primal_sum, 0.0);
  if (!error.ok()) return error;
  error = backend.fill(dual_sum, 0.0);
  if (!error.ok()) return error;
  state.iteration = 0;
  averaging = resetAveraging();
  average_normalization = 0.0;
  // The boundary point becomes the state that continues into the next step.
  return recordPoint(backend, *state.current_primal, *state.current_dual,
                     *state.current_ax, *state.current_aty, point);
}

// The check is evaluated and stored before any iterate, anchor, age, average,
// or sum is mutated.
api::Error observeAndApplyFixedPointRestart(
    backend::Backend& backend, const FixedTraceOptions& options,
    HalpernPdhgState& state, backend::VectorHandle& aty_delta_scratch,
    backend::VectorHandle& primal_sum, backend::VectorHandle& dual_sum,
    AveragingState& averaging, double& average_normalization,
    const HighsInt iteration, TraceHandles& point,
    std::vector<FixedTraceFixedPointRestartCheck>& checks) {
  const IterateSlots slots = {state.current_primal, state.current_dual,
                              state.next_primal, state.next_dual};
  const api::Error error = appendFixedPointRestartCheck(
      backend, options, state, slots, aty_delta_scratch, iteration, checks);
  if (!error.ok()) return error;
  if (checks.back().evaluation.decision != RestartDecision::kRestart)
    return api::Error();
  return applyFixedPointRestart(backend, state, primal_sum, dual_sum, averaging,
                                average_normalization, point);
}

}  // namespace

FixedTraceConvergenceOptions::FixedTraceConvergenceOptions()
    : periodic_enabled(false), check_cadence(0) {}

FixedTraceFixedPointErrorOptions::FixedTraceFixedPointErrorOptions()
    : enabled(false),
      omega(std::numeric_limits<double>::quiet_NaN()),
      eta(std::numeric_limits<double>::quiet_NaN()) {}

FixedTraceFixedPointRestartOptions::FixedTraceFixedPointRestartOptions()
    : enabled(false), check_cadence(0) {}

FixedTraceOptions::FixedTraceOptions()
    : backend(backend::Kind::kCpu),
      execution_mode(FixedTraceExecutionMode::kDirect),
      primal_step_size(std::numeric_limits<double>::quiet_NaN()),
      dual_step_size(std::numeric_limits<double>::quiet_NaN()),
      iteration_count(-1),
      convergence() {
  convergence_tolerances.primal = kDefaultFixedTraceTolerance;
  convergence_tolerances.dual = kDefaultFixedTraceTolerance;
  convergence_tolerances.gap = kDefaultFixedTraceTolerance;
}

FixedTraceConvergenceObservation::FixedTraceConvergenceObservation()
    : observed(false) {}

FixedTraceConvergenceCheck::FixedTraceConvergenceCheck() : iteration(0) {}

FixedTraceFixedPointErrorObservation::FixedTraceFixedPointErrorObservation()
    : observed(false) {}

FixedTraceFixedPointRestartObservation::FixedTraceFixedPointRestartObservation()
    : observed(false) {}

FixedTraceFixedPointRestartCheck::FixedTraceFixedPointRestartCheck()
    : iteration(0) {}

FixedTraceDiagnostics::FixedTraceDiagnostics()
    : resource_creation_calls_before_iterations(0),
      resource_creation_calls_during_iterations(0),
      transfer_calls_before_iterations(0),
      transfer_calls_during_iterations(0),
      transfer_calls_after_iterations(0) {}

bool FixedTraceResult::ok() const { return error.ok(); }

FixedTraceResult runFixedTrace(const model::OriginalModel& original_model,
                               const FixedTraceOptions& options) {
  FixedTraceDiagnostics diagnostics;
  if (options.execution_mode != FixedTraceExecutionMode::kDirect &&
      options.execution_mode != FixedTraceExecutionMode::kRecorded) {
    return failure(api::Error(api::ErrorCode::kInvalidArgument,
                              "fixed trace execution mode is invalid"),
                   diagnostics);
  }
  if (options.convergence.periodic_enabled &&
      options.convergence.check_cadence < 0) {
    return failure(api::Error(api::ErrorCode::kInvalidOption,
                              "fixed trace convergence check cadence must be "
                              "nonnegative"),
                   diagnostics);
  }
  if (options.convergence.periodic_enabled &&
      options.execution_mode == FixedTraceExecutionMode::kRecorded &&
      options.convergence.check_cadence > 0 &&
      options.convergence.check_cadence < options.iteration_count &&
      options.convergence.check_cadence % kRecordedBlockIterations != 0) {
    return failure(
        api::Error(api::ErrorCode::kInvalidOption,
                   "fixed trace convergence check cadence must align with "
                   "recorded iteration blocks"),
        diagnostics);
  }
  if (options.fixed_point_error.enabled &&
      (!std::isfinite(options.fixed_point_error.omega) ||
       options.fixed_point_error.omega <= 0.0 ||
       !std::isfinite(options.fixed_point_error.eta) ||
       options.fixed_point_error.eta <= 0.0)) {
    return failure(api::Error(api::ErrorCode::kInvalidArgument,
                              "fixed trace fixed-point-error omega and eta "
                              "must be finite and positive"),
                   diagnostics);
  }
  if (options.fixed_point_restart.enabled &&
      !options.fixed_point_error.enabled) {
    return failure(api::Error(api::ErrorCode::kInvalidOption,
                              "fixed trace fixed-point restart trigger "
                              "requires the fixed-point-error observation"),
                   diagnostics);
  }
  if (options.fixed_point_restart.enabled &&
      options.fixed_point_restart.check_cadence < 0) {
    return failure(api::Error(api::ErrorCode::kInvalidOption,
                              "fixed trace fixed-point restart check cadence "
                              "must be nonnegative"),
                   diagnostics);
  }
  // Recorded blocks are never split, so a cadence that is not a multiple of the
  // block size would demand a check inside a block.
  if (options.fixed_point_restart.enabled &&
      options.execution_mode == FixedTraceExecutionMode::kRecorded &&
      options.fixed_point_restart.check_cadence > 0 &&
      options.fixed_point_restart.check_cadence < options.iteration_count &&
      options.fixed_point_restart.check_cadence % kRecordedBlockIterations !=
          0) {
    return failure(api::Error(api::ErrorCode::kInvalidOption,
                              "fixed trace fixed-point restart check cadence "
                              "must align with recorded iteration blocks"),
                   diagnostics);
  }
  const bool uses_step_sequences =
      !options.primal_step_sizes.empty() || !options.dual_step_sizes.empty();
  if (!uses_step_sequences && (!std::isfinite(options.primal_step_size) ||
                               options.primal_step_size <= 0.0 ||
                               !std::isfinite(options.dual_step_size) ||
                               options.dual_step_size <= 0.0)) {
    return failure(api::Error(api::ErrorCode::kInvalidArgument,
                              "fixed trace step sizes must be finite and "
                              "positive"),
                   diagnostics);
  }
  if (options.iteration_count < 0) {
    return failure(api::Error(api::ErrorCode::kInvalidArgument,
                              "fixed trace iteration count must be "
                              "nonnegative"),
                   diagnostics);
  }
  if (options.iteration_count == std::numeric_limits<HighsInt>::max()) {
    return failure(api::Error(api::ErrorCode::kInvalidArgument,
                              "fixed trace iteration count is too large"),
                   diagnostics);
  }
  if (uses_step_sequences) {
    const std::size_t iteration_count =
        static_cast<std::size_t>(options.iteration_count);
    if (options.primal_step_sizes.size() != iteration_count ||
        options.dual_step_sizes.size() != iteration_count) {
      return failure(
          api::Error(api::ErrorCode::kInvalidArgument,
                     "fixed trace step-size sequences must match the "
                     "iteration count"),
          diagnostics);
    }
    for (std::size_t iteration = 0; iteration < iteration_count; ++iteration) {
      if (!std::isfinite(options.primal_step_sizes[iteration]) ||
          options.primal_step_sizes[iteration] <= 0.0 ||
          !std::isfinite(options.dual_step_sizes[iteration]) ||
          options.dual_step_sizes[iteration] <= 0.0) {
        return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                  "fixed trace step sizes must be finite and "
                                  "positive"),
                       diagnostics);
      }
    }
  }

  try {
    model::ProblemBuildResult built = model::buildProblem(original_model.lp());
    if (!built.ok()) return failure(built.error, diagnostics);
    api::Error error =
        model::applyScalingPipeline(built.problem, options.scaling);
    if (!error.ok()) return failure(error, diagnostics);

    const HighsInt num_cols = built.problem.numCols();
    const HighsInt num_rows = built.problem.numRows();
    if (options.scaled_primal_start.size() !=
            static_cast<std::size_t>(num_cols) ||
        options.scaled_dual_start.size() !=
            static_cast<std::size_t>(num_rows)) {
      return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                "fixed trace initial vector dimensions do "
                                "not match the transformed problem"),
                     diagnostics);
    }
    if (!finiteVector(options.scaled_primal_start) ||
        !finiteVector(options.scaled_dual_start)) {
      return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                "fixed trace initial vectors must be finite"),
                     diagnostics);
    }

    const std::size_t point_count =
        static_cast<std::size_t>(options.iteration_count) + 1;
    std::vector<TraceHandles> trace_handles(point_count);
    std::vector<FixedTracePoint> host_trace(point_count);
    for (std::size_t point = 0; point < point_count; ++point) {
      host_trace[point].primal.resize(num_cols);
      host_trace[point].dual.resize(num_rows);
      host_trace[point].ax.resize(num_rows);
      host_trace[point].aty.resize(num_cols);
    }

    backend::FactoryResult factory = backend::createBackend(options.backend);
    if (!factory.ok()) return failure(factory.error, diagnostics);
    std::unique_ptr<backend::Backend> backend = std::move(factory.backend);
    const backend::Descriptor descriptor = backend->descriptor();
    if (!descriptor.precision.isFp64()) {
      return failure(
          api::Error(api::ErrorCode::kUnsupportedCapability,
                     "fixed traces require deterministic FP64 storage, "
                     "compute, and reduction precision"),
          diagnostics);
    }
    error = requireTraceCapabilities(*backend);
    if (!error.ok()) return failure(error, diagnostics);
    if (options.fixed_point_error.enabled && options.iteration_count > 0) {
      error = backend->requireCapability(backend::Capability::kReductions);
      if (!error.ok()) return failure(error, diagnostics);
    }
    if (options.execution_mode == FixedTraceExecutionMode::kRecorded) {
      error = backend->requireCapability(backend::Capability::kCommandBuffer);
      if (!error.ok()) return failure(error, diagnostics);
    }

    backend::ProblemHandle problem;
    ++diagnostics.resource_creation_calls_before_iterations;
    error = backend->bindProblem(built.problem, problem);
    if (!error.ok()) return failure(error, diagnostics);
    backend::WorkspaceHandle workspace;
    ++diagnostics.resource_creation_calls_before_iterations;
    error = backend->createWorkspace(problem, workspace);
    if (!error.ok()) return failure(error, diagnostics);

    backend::VectorHandle current_primal;
    backend::VectorHandle current_dual;
    backend::VectorHandle next_primal;
    backend::VectorHandle next_dual;
    backend::VectorHandle primal_anchor;
    backend::VectorHandle dual_anchor;
    backend::VectorHandle primal_candidate;
    backend::VectorHandle dual_candidate;
    backend::VectorHandle reflected_primal;
    backend::VectorHandle reflected_dual;
    backend::VectorHandle current_ax;
    backend::VectorHandle current_aty;
    backend::VectorHandle reflected_ax;
    backend::VectorHandle primal_sum;
    backend::VectorHandle dual_sum;
    backend::VectorHandle averaged_ax;
    backend::VectorHandle averaged_aty;
    backend::VectorHandle candidate_ax;
    backend::VectorHandle candidate_aty;
    backend::VectorHandle* col_vectors[] = {
        &current_primal,   &next_primal, &primal_anchor, &primal_candidate,
        &reflected_primal, &current_aty, &primal_sum,    &averaged_aty};
    backend::VectorHandle* row_vectors[] = {
        &current_dual,   &next_dual,      &dual_anchor,
        &dual_candidate, &reflected_dual, &current_ax,
        &reflected_ax,   &dual_sum,       &averaged_ax};
    for (std::size_t index = 0;
         index < sizeof(col_vectors) / sizeof(col_vectors[0]); ++index) {
      error =
          createVector(*backend, num_cols, *col_vectors[index], diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
    }
    for (std::size_t index = 0;
         index < sizeof(row_vectors) / sizeof(row_vectors[0]); ++index) {
      error =
          createVector(*backend, num_rows, *row_vectors[index], diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
    }
    for (std::size_t point = 0; point < point_count; ++point) {
      error = createVector(*backend, num_cols, trace_handles[point].primal,
                           diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
      error = createVector(*backend, num_rows, trace_handles[point].dual,
                           diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
      error = createVector(*backend, num_rows, trace_handles[point].ax,
                           diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
      error = createVector(*backend, num_cols, trace_handles[point].aty,
                           diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
    }
    if (options.convergence.periodic_enabled) {
      error = createVector(*backend, num_rows, candidate_ax, diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
      error = createVector(*backend, num_cols, candidate_aty, diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
    }

    error = backend->uploadVector(options.scaled_primal_start, current_primal);
    ++diagnostics.transfer_calls_before_iterations;
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->uploadVector(options.scaled_dual_start, current_dual);
    ++diagnostics.transfer_calls_before_iterations;
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->fill(primal_sum, 0.0);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->fill(dual_sum, 0.0);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->fill(current_aty, 0.0);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->fusedProjectedPrimalUpdate(
        problem, workspace, current_primal, current_aty, 0.0, next_primal);
    if (!error.ok()) return failure(error, diagnostics);
    std::swap(current_primal, next_primal);
    error = backend->copy(current_primal, primal_anchor);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->copy(current_dual, dual_anchor);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->multiply(problem, current_primal, current_ax);
    if (!error.ok()) return failure(error, diagnostics);
    error = backend->multiplyTranspose(problem, current_dual, current_aty);
    if (!error.ok()) return failure(error, diagnostics);
    error = recordPoint(*backend, current_primal, current_dual, current_ax,
                        current_aty, trace_handles[0]);
    if (!error.ok()) return failure(error, diagnostics);

    HalpernPdhgState state;
    state.problem = &problem;
    state.workspace = &workspace;
    state.current_primal = &current_primal;
    state.current_dual = &current_dual;
    state.next_primal = &next_primal;
    state.next_dual = &next_dual;
    state.primal_anchor = &primal_anchor;
    state.dual_anchor = &dual_anchor;
    state.primal_candidate = &primal_candidate;
    state.dual_candidate = &dual_candidate;
    state.reflected_primal = &reflected_primal;
    state.reflected_dual = &reflected_dual;
    state.current_ax = &current_ax;
    state.current_aty = &current_aty;
    state.reflected_ax = &reflected_ax;
    state.primal_step_size = options.primal_step_size;
    state.dual_step_size = options.dual_step_size;

    AveragingState averaging;
    double average_normalization = 0.0;
    // The monotone total is independent of the resettable Halpern age in
    // HalpernPdhgState::iteration.
    HighsInt total_iterations = 0;
    std::vector<FixedTraceConvergenceCheck> convergence_checks;
    if (options.convergence.periodic_enabled && options.iteration_count > 0) {
      const HighsInt cadence = options.convergence.check_cadence;
      const HighsInt check_count =
          cadence == 0 ? 1
                       : options.iteration_count / cadence +
                             (options.iteration_count % cadence == 0 ? 0 : 1);
      convergence_checks.reserve(static_cast<std::size_t>(check_count));
    }
    std::vector<FixedTraceFixedPointRestartCheck> restart_checks;
    if (options.execution_mode == FixedTraceExecutionMode::kDirect) {
      while (total_iterations < options.iteration_count) {
        if (uses_step_sequences) {
          const std::size_t iteration =
              static_cast<std::size_t>(total_iterations);
          state.primal_step_size = options.primal_step_sizes[iteration];
          state.dual_step_size = options.dual_step_sizes[iteration];
        }
        error = enqueueHalpernPdhgStep(*backend, state);
        if (!error.ok()) return failure(error, diagnostics);
        ++total_iterations;
        AveragingInput averaging_input;
        averaging_input.previous = averaging;
        averaging_input.primal_step_size = state.primal_step_size;
        averaging_input.dual_step_size = state.dual_step_size;
        averaging_input.completed_halpern_iterate = true;
        const AveragingResult averaging_result =
            updateAveraging(averaging_input);
        if (!averaging_result.ok())
          return failure(averaging_result.error, diagnostics);
        error = backend->fusedAverageUpdate(
            workspace, *state.current_primal, *state.current_dual,
            averaging_result.sample_weight, primal_sum, dual_sum);
        if (!error.ok()) return failure(error, diagnostics);
        averaging = averaging_result.next;
        average_normalization = averaging_result.normalization;
        TraceHandles& point =
            trace_handles[static_cast<std::size_t>(total_iterations)];
        error =
            recordPoint(*backend, *state.current_primal, *state.current_dual,
                        *state.current_ax, *state.current_aty, point);
        if (!error.ok()) return failure(error, diagnostics);
        if (periodicConvergenceCheck(options, total_iterations)) {
          error = appendConvergenceCheck(
              *backend, built.problem, options, state, candidate_ax,
              candidate_aty, total_iterations, diagnostics, convergence_checks);
          if (!error.ok()) return failure(error, diagnostics);
        }
        if (intermediateFixedPointCheck(options, total_iterations)) {
          error = observeAndApplyFixedPointRestart(
              *backend, options, state, averaged_aty, primal_sum, dual_sum,
              averaging, average_normalization, total_iterations, point,
              restart_checks);
          if (!error.ok()) return failure(error, diagnostics);
        }
      }
    } else {
      backend::CommandBufferRecorder* recorder =
          backend::commandBufferRecorder(*backend);
      if (recorder == nullptr) {
        return failure(
            api::Error(api::ErrorCode::kUnsupportedCapability,
                       "selected backend has no command-buffer recorder"),
            diagnostics);
      }
      std::vector<backend::CommandBufferHandle> command_buffers;
      command_buffers.reserve(
          static_cast<std::size_t>((options.iteration_count + 1) / 2));
      // Each block is replayed, and its restart decision applied, before the
      // next block is captured, because a captured anchor weight bakes in the
      // Halpern age this block's decision may reset.
      while (total_iterations < options.iteration_count) {
        error = recorder->beginRecording();
        if (!error.ok()) return failure(error, diagnostics);
        const HighsInt block_end =
            std::min(options.iteration_count,
                     total_iterations + kRecordedBlockIterations);
        while (total_iterations < block_end) {
          if (uses_step_sequences) {
            const std::size_t iteration =
                static_cast<std::size_t>(total_iterations);
            state.primal_step_size = options.primal_step_sizes[iteration];
            state.dual_step_size = options.dual_step_sizes[iteration];
          }
          error = enqueueHalpernPdhgStep(*backend, state);
          if (!error.ok()) {
            recorder->cancelRecording();
            return failure(error, diagnostics);
          }
          ++total_iterations;
          AveragingInput averaging_input;
          averaging_input.previous = averaging;
          averaging_input.primal_step_size = state.primal_step_size;
          averaging_input.dual_step_size = state.dual_step_size;
          averaging_input.completed_halpern_iterate = true;
          const AveragingResult averaging_result =
              updateAveraging(averaging_input);
          if (!averaging_result.ok()) {
            recorder->cancelRecording();
            return failure(averaging_result.error, diagnostics);
          }
          error = backend->fusedAverageUpdate(
              workspace, *state.current_primal, *state.current_dual,
              averaging_result.sample_weight, primal_sum, dual_sum);
          if (!error.ok()) {
            recorder->cancelRecording();
            return failure(error, diagnostics);
          }
          averaging = averaging_result.next;
          average_normalization = averaging_result.normalization;
          error = recordPoint(
              *backend, *state.current_primal, *state.current_dual,
              *state.current_ax, *state.current_aty,
              trace_handles[static_cast<std::size_t>(total_iterations)]);
          if (!error.ok()) {
            recorder->cancelRecording();
            return failure(error, diagnostics);
          }
        }
        backend::CommandBufferHandle command_buffer;
        error = recorder->endRecording(command_buffer);
        if (!error.ok()) return failure(error, diagnostics);
        ++diagnostics.resource_creation_calls_before_iterations;
        command_buffers.push_back(std::move(command_buffer));
        backend::CommandBufferStatistics block_statistics;
        error = recorder->replay(command_buffers.back(), block_statistics);
        backend::addCommandBufferStatistics(block_statistics,
                                            diagnostics.command_buffer);
        if (!error.ok()) return failure(error, diagnostics);
        if (periodicConvergenceCheck(options, total_iterations)) {
          error = appendConvergenceCheck(
              *backend, built.problem, options, state, candidate_ax,
              candidate_aty, total_iterations, diagnostics, convergence_checks);
          if (!error.ok()) return failure(error, diagnostics);
        }
        if (intermediateFixedPointCheck(options, total_iterations)) {
          error = observeAndApplyFixedPointRestart(
              *backend, options, state, averaged_aty, primal_sum, dual_sum,
              averaging, average_normalization, total_iterations,
              trace_handles[static_cast<std::size_t>(total_iterations)],
              restart_checks);
          if (!error.ok()) return failure(error, diagnostics);
        }
      }
    }

    FixedTraceFixedPointErrorObservation fixed_point_error;
    FixedTraceFixedPointRestartObservation fixed_point_restart;
    if (options.fixed_point_error.enabled && options.iteration_count > 0) {
      if (options.fixed_point_restart.enabled) {
        error = observeAndApplyFixedPointRestart(
            *backend, options, state, averaged_aty, primal_sum, dual_sum,
            averaging, average_normalization, total_iterations,
            trace_handles[point_count - 1], restart_checks);
        if (!error.ok()) return failure(error, diagnostics);
        fixed_point_error.observed = true;
        fixed_point_error.evaluation = restart_checks.back().error;
        fixed_point_restart.observed = true;
        fixed_point_restart.evaluation = restart_checks.back().evaluation;
      } else {
        FixedPointErrorInput input;
        input.next_primal = state.primal_candidate;
        input.next_dual = state.dual_candidate;
        input.reflected_primal = state.reflected_primal;
        input.reflected_dual = state.reflected_dual;
        input.omega = options.fixed_point_error.omega;
        input.eta = options.fixed_point_error.eta;
        input.delta_primal = state.next_primal;
        input.delta_dual = state.next_dual;
        input.aty_delta = &averaged_aty;
        const FixedPointErrorResult evaluation =
            evaluateFixedPointError(*backend, problem, input);
        if (!evaluation.ok()) return failure(evaluation.error, diagnostics);
        fixed_point_error.observed = true;
        fixed_point_error.evaluation = evaluation;
      }
    }

    if (averaging.hasAverage()) {
      error = backend->scale(average_normalization, primal_sum);
      if (!error.ok()) return failure(error, diagnostics);
      error = backend->scale(average_normalization, dual_sum);
      if (!error.ok()) return failure(error, diagnostics);
      error = backend->multiply(problem, primal_sum, averaged_ax);
      if (!error.ok()) return failure(error, diagnostics);
      error = backend->multiplyTranspose(problem, dual_sum, averaged_aty);
      if (!error.ok()) return failure(error, diagnostics);
    }
    error = backend->synchronize();
    if (!error.ok()) return failure(error, diagnostics);
    for (std::size_t point = 0; point < point_count; ++point) {
      error = downloadPoint(*backend, trace_handles[point], host_trace[point],
                            diagnostics);
      if (!error.ok()) return failure(error, diagnostics);
    }

    FixedTraceResult result;
    if (averaging.hasAverage()) {
      error = backend->downloadVector(primal_sum, result.averaged_primal);
      ++diagnostics.transfer_calls_after_iterations;
      if (!error.ok()) return failure(error, diagnostics);
      error = backend->downloadVector(dual_sum, result.averaged_dual);
      ++diagnostics.transfer_calls_after_iterations;
      if (!error.ok()) return failure(error, diagnostics);
      std::vector<double> averaged_ax_values;
      std::vector<double> averaged_aty_values;
      error = backend->downloadVector(averaged_ax, averaged_ax_values);
      ++diagnostics.transfer_calls_after_iterations;
      if (!error.ok()) return failure(error, diagnostics);
      error = backend->downloadVector(averaged_aty, averaged_aty_values);
      ++diagnostics.transfer_calls_after_iterations;
      if (!error.ok()) return failure(error, diagnostics);
      error = observeConvergence(built.problem, options, result.averaged_primal,
                                 result.averaged_dual, averaged_ax_values,
                                 averaged_aty_values, options.iteration_count,
                                 result.averaged_convergence);
      if (!error.ok()) return failure(error, diagnostics);
    }
    const FixedTracePoint& final_point = host_trace[point_count - 1];
    error =
        observeConvergence(built.problem, options, final_point.primal,
                           final_point.dual, final_point.ax, final_point.aty,
                           options.iteration_count, result.current_convergence);
    if (!error.ok()) return failure(error, diagnostics);
    result.points.swap(host_trace);
    result.averaging = averaging;
    result.convergence_checks.swap(convergence_checks);
    result.fixed_point_error = fixed_point_error;
    result.fixed_point_restart = fixed_point_restart;
    result.fixed_point_restart_checks.swap(restart_checks);
    result.diagnostics = diagnostics;
    result.scaling = built.problem.scaling();
    return result;
  } catch (const std::bad_alloc&) {
    return failure(api::Error(api::ErrorCode::kAllocationFailure,
                              "fixed trace allocation failed"),
                   diagnostics);
  }
}

}  // namespace algorithm
}  // namespace hipdlp