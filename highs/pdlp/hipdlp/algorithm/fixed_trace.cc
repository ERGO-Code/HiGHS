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

struct TraceHandles {
  backend::VectorHandle primal;
  backend::VectorHandle dual;
  backend::VectorHandle ax;
  backend::VectorHandle aty;
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

api::Error observeConvergence(const model::Problem& problem,
                              const FixedTraceOptions& options,
                              const std::vector<double>& primal,
                              const std::vector<double>& dual,
                              const std::vector<double>& ax,
                              const std::vector<double>& aty,
                              FixedTraceConvergenceObservation& observation) {
  RawConvergenceMetrics raw;
  const api::Error error =
      computeRawMetrics(problem, primal, dual, ax, aty, raw);
  if (!error.ok()) return error;
  StoppingInput input;
  input.raw = raw;
  input.tolerances = options.convergence_tolerances;
  input.iteration = options.iteration_count;
  const StoppingResult stopping = evaluateStopping(input);
  if (!stopping.ok()) return stopping.error;
  observation.observed = true;
  observation.raw = raw;
  observation.stopping = stopping;
  return api::Error();
}

}  // namespace

FixedTraceOptions::FixedTraceOptions()
    : backend(backend::Kind::kCpu),
      execution_mode(FixedTraceExecutionMode::kDirect),
      primal_step_size(std::numeric_limits<double>::quiet_NaN()),
      dual_step_size(std::numeric_limits<double>::quiet_NaN()),
      iteration_count(-1) {
  convergence_tolerances.primal = kDefaultFixedTraceTolerance;
  convergence_tolerances.dual = kDefaultFixedTraceTolerance;
  convergence_tolerances.gap = kDefaultFixedTraceTolerance;
}

FixedTraceConvergenceObservation::FixedTraceConvergenceObservation()
    : observed(false) {}

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
    if (options.execution_mode == FixedTraceExecutionMode::kDirect) {
      while (state.iteration < options.iteration_count) {
        if (uses_step_sequences) {
          const std::size_t iteration =
              static_cast<std::size_t>(state.iteration);
          state.primal_step_size = options.primal_step_sizes[iteration];
          state.dual_step_size = options.dual_step_sizes[iteration];
        }
        error = enqueueHalpernPdhgStep(*backend, state);
        if (!error.ok()) return failure(error, diagnostics);
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
        error = recordPoint(
            *backend, *state.current_primal, *state.current_dual,
            *state.current_ax, *state.current_aty,
            trace_handles[static_cast<std::size_t>(state.iteration)]);
        if (!error.ok()) return failure(error, diagnostics);
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
      while (state.iteration < options.iteration_count) {
        error = recorder->beginRecording();
        if (!error.ok()) return failure(error, diagnostics);
        const HighsInt block_end =
            std::min(options.iteration_count, state.iteration + 2);
        while (state.iteration < block_end) {
          if (uses_step_sequences) {
            const std::size_t iteration =
                static_cast<std::size_t>(state.iteration);
            state.primal_step_size = options.primal_step_sizes[iteration];
            state.dual_step_size = options.dual_step_sizes[iteration];
          }
          error = enqueueHalpernPdhgStep(*backend, state);
          if (!error.ok()) {
            recorder->cancelRecording();
            return failure(error, diagnostics);
          }
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
              trace_handles[static_cast<std::size_t>(state.iteration)]);
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
      }
      for (std::size_t block = 0; block < command_buffers.size(); ++block) {
        backend::CommandBufferStatistics block_statistics;
        error = recorder->replay(command_buffers[block], block_statistics);
        backend::addCommandBufferStatistics(block_statistics,
                                            diagnostics.command_buffer);
        if (!error.ok()) return failure(error, diagnostics);
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
      error = observeConvergence(
          built.problem, options, result.averaged_primal, result.averaged_dual,
          averaged_ax_values, averaged_aty_values, result.averaged_convergence);
      if (!error.ok()) return failure(error, diagnostics);
    }
    const FixedTracePoint& final_point = host_trace[point_count - 1];
    error = observeConvergence(built.problem, options, final_point.primal,
                               final_point.dual, final_point.ax,
                               final_point.aty, result.current_convergence);
    if (!error.ok()) return failure(error, diagnostics);
    result.points.swap(host_trace);
    result.averaging = averaging;
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