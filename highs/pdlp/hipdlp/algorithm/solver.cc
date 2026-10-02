#include "pdlp/hipdlp/algorithm/solver.hpp"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <limits>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

#include "lp_data/HConst.h"
#include "pdlp/hipdlp/algorithm/convergence.hpp"
#include "pdlp/hipdlp/algorithm/fixed_point_error.hpp"
#include "pdlp/hipdlp/algorithm/halpern_pdhg_step.hpp"
#include "pdlp/hipdlp/algorithm/operator_norm.hpp"
#include "pdlp/hipdlp/algorithm/restart_policy.hpp"
#include "pdlp/hipdlp/algorithm/step_size_policy.hpp"
#include "pdlp/hipdlp/backend/backend.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/backend/check.hpp"
#include "pdlp/hipdlp/backend/command_buffer.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"
#include "pdlp/hipdlp/verification/certificate_verifier.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

namespace hipdlp {
namespace algorithm {
namespace {

const HighsInt kConvergenceCadence = 40;
const HighsInt kRecordedStepCount = kConvergenceCadence - 1;
// On a device backend a host-side original-KKT verification costs about as
// much as a boundary of iterations, so rejected verifications back off to at
// most this many boundaries, which bounds how late a candidate is confirmed.
const HighsInt kMaximumVerificationGap = 16;

bool iteratesOnDevice(const backend::Kind kind) {
  return kind == backend::Kind::kCuda;
}

class RecordedIterationBlock {
 public:
  bool ready() const { return commands_.valid(); }

  api::Error prepare(backend::Backend& selected_backend,
                      const HalpernPdhgState& initial) {
    recorder_ = backend::commandBufferRecorder(selected_backend);
    if (recorder_ == nullptr) return api::Error();
    values_.resize(3 * kRecordedStepCount);
    api::Error error = selected_backend.createVector(
        static_cast<HighsInt>(values_.size()), parameters_);
    if (!error.ok()) return error;
    HalpernPdhgState recorded = initial;
    std::swap(recorded.current_primal, recorded.next_primal);
    std::swap(recorded.current_dual, recorded.next_dual);
    ++recorded.iteration;
    updateValues(recorded);
    error = selected_backend.uploadVector(values_, parameters_);
    if (!error.ok()) return error;
    error = recorder_->beginRecording();
    if (!error.ok()) return error;
    try {
      for (HighsInt offset = 0; offset < kRecordedStepCount; ++offset) {
        error = enqueueHalpernPdhgStep(selected_backend, recorded, &parameters_,
                                       3 * offset);
        if (!error.ok()) {
          recorder_->cancelRecording();
          return error;
        }
      }
    } catch (...) {
      recorder_->cancelRecording();
      throw;
    }
    return recorder_->endRecording(commands_);
  }

  api::Error replay(backend::Backend& selected_backend, HalpernPdhgState& state) {
    updateValues(state);
    api::Error error = selected_backend.uploadVector(values_, parameters_);
    if (!error.ok()) return error;
    backend::CommandBufferStatistics statistics;
    error = recorder_->replay(commands_, statistics);
    if (!error.ok()) return error;
    std::swap(state.current_primal, state.next_primal);
    std::swap(state.current_dual, state.next_dual);
    state.iteration += kRecordedStepCount;
    return api::Error();
  }

 private:
  void updateValues(const HalpernPdhgState& state) {
    for (HighsInt offset = 0; offset < kRecordedStepCount; ++offset) {
      values_[3 * offset] = state.primal_step_size;
      values_[3 * offset + 1] = state.dual_step_size;
      values_[3 * offset + 2] =
          1.0 / (static_cast<double>(state.iteration + offset + 1) + 1.0);
    }
  }

  backend::CommandBufferRecorder* recorder_ = nullptr;
  std::vector<double> values_;
  backend::VectorHandle parameters_;
  backend::CommandBufferHandle commands_;
};

typedef std::chrono::steady_clock SolveClock;

double* phaseSeconds(SolveDiagnostics* diagnostics, const SolvePhase phase) {
  if (diagnostics == nullptr) return nullptr;
  switch (phase) {
    case SolvePhase::kSolve:
      return &diagnostics->total_wall_seconds;
    case SolvePhase::kPreprocessing:
      return &diagnostics->preprocessing_seconds;
    case SolvePhase::kBackendSetup:
      return &diagnostics->backend_setup_seconds;
    case SolvePhase::kOperatorNorm:
      return &diagnostics->operator_norm_seconds;
    case SolvePhase::kIterations:
      return &diagnostics->iterations_seconds;
    case SolvePhase::kIterationBlock:
      return &diagnostics->iteration_block_seconds;
    case SolvePhase::kStep:
      return &diagnostics->step_seconds;
    case SolvePhase::kBoundary:
      return &diagnostics->boundary_seconds;
    case SolvePhase::kDownload:
      return &diagnostics->download_seconds;
    case SolvePhase::kHostPreparation:
      return &diagnostics->host_preparation_seconds;
    case SolvePhase::kHostVerification:
      return &diagnostics->host_verification_seconds;
    case SolvePhase::kReferenceFixedPoint:
      return &diagnostics->reference_fixed_point_seconds;
    case SolvePhase::kFixedPointFallback:
      return &diagnostics->fixed_point_fallback_seconds;
    case SolvePhase::kCount:
      return nullptr;
  }
  return nullptr;
}

class SolvePhaseGuard {
 public:
  SolvePhaseGuard(const SolveContext& context, const SolvePhase phase,
                  const HighsInt& iteration, const bool begin_now = true)
      : observer_(context.observer),
        seconds_(phaseSeconds(context.diagnostics, phase)),
        phase_(phase), iteration_(iteration),
        active_(false) {
    if (begin_now) begin();
  }

  ~SolvePhaseGuard() { finish(); }

  void begin() {
    if (active_) return;
    active_ = true;
    if (observer_ != nullptr)
      observer_->onEvent({phase_, true, iteration_, false});
    if (seconds_ != nullptr) start_ = SolveClock::now();
  }

  void finish(const bool continued = false) {
    if (!active_) return;
    if (seconds_ != nullptr)
      *seconds_ +=
          std::chrono::duration<double>(SolveClock::now() - start_).count();
    active_ = false;
    if (observer_ != nullptr)
      observer_->onEvent({phase_, false, iteration_, continued});
  }

  SolvePhaseGuard(const SolvePhaseGuard&) = delete;
  SolvePhaseGuard& operator=(const SolvePhaseGuard&) = delete;

 private:
  SolveObserver* observer_;
  double* seconds_;
  SolveClock::time_point start_;
  SolvePhase phase_;
  const HighsInt& iteration_;
  bool active_;
};

struct IterationVectors {
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
  backend::VectorHandle candidate_ax;
  backend::VectorHandle candidate_aty;
  backend::VectorHandle aty_delta;
};

struct HostIterate {
  std::vector<double> primal;
  std::vector<double> dual;
  std::vector<double> ax;
  std::vector<double> aty;
  std::vector<double> reduced_costs;
};

// Reused across boundaries so certificate detection allocates nothing new
// once the loop is running.
struct CertificateBuffers {
  std::vector<double> scaled_column_multiplier;
  std::vector<double> primal_direction;
  std::vector<double> dual_direction;
  std::vector<double> column_multiplier;
};

struct CertificateOutcome {
  bool available;
  api::SolverStatus status;
  api::Certificate certificate;

  CertificateOutcome()
      : available(false), status(api::SolverStatus::kNotSet) {}
};

api::SolveResult failure(const api::Error& error) {
  return api::SolveResult::failure(error.code, error.message);
}

api::Error validateOptions(const api::SolverOptions& options) {
  if (options.implementation != api::SolverImplementation::kRewrite) {
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP rewrite solve requires the rewrite selector");
  }
  if (!backend::isValidKind(options.backend)) {
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP rewrite solve received an invalid backend");
  }
  if (options.policy.step_size != api::StepSizePolicy::kFixed &&
      options.policy.step_size != api::StepSizePolicy::kPid) {
    return api::Error(
        api::ErrorCode::kInvalidOption,
        "HiPDLP rewrite solve received an invalid step-size policy");
  }
  if (options.policy.restart != api::RestartPolicy::kHalpern) {
    return api::Error(
        api::ErrorCode::kInvalidOption,
        "HiPDLP rewrite solve received an invalid restart policy");
  }
  if ((options.time_limit != kHighsInf && !std::isfinite(options.time_limit)) ||
      options.time_limit < 0.0) {
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP rewrite time limit is invalid");
  }
  if (options.iteration_limit < 0) {
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP rewrite iteration limit is invalid");
  }
  if (!std::isfinite(options.optimality_tolerance) ||
      options.optimality_tolerance <= 0.0) {
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP rewrite optimality tolerance is invalid");
  }
  if (!std::isfinite(options.kkt_tolerance) || options.kkt_tolerance <= 0.0) {
    return api::Error(api::ErrorCode::kInvalidOption,
                      "HiPDLP rewrite KKT tolerance is invalid");
  }
  return api::Error();
}

api::Error requireSolveCapabilities(backend::Backend& selected_backend) {
  const backend::Capability required[] = {
      backend::Capability::kSparseMatrixVector,
      backend::Capability::kFusedPrimalDualUpdate,
      backend::Capability::kVectorOperations, backend::Capability::kReductions,
      backend::Capability::kHostDownload};
  for (std::size_t index = 0; index < sizeof(required) / sizeof(required[0]);
       ++index) {
    const api::Error error =
        selected_backend.requireCapability(required[index]);
    if (!error.ok()) return error;
  }
  return api::Error();
}

api::Error createIterationVectors(backend::Backend& selected_backend,
                                  const HighsInt num_rows,
                                  const HighsInt num_cols,
                                  IterationVectors& vectors,
                                  SolveDiagnostics* diagnostics) {
  backend::VectorHandle* col_vectors[] = {
      &vectors.current_primal,   &vectors.next_primal,
      &vectors.primal_anchor,    &vectors.primal_candidate,
      &vectors.reflected_primal, &vectors.current_aty,
      &vectors.candidate_aty,    &vectors.aty_delta};
  backend::VectorHandle* row_vectors[] = {
      &vectors.current_dual,   &vectors.next_dual,      &vectors.dual_anchor,
      &vectors.dual_candidate, &vectors.reflected_dual, &vectors.current_ax,
      &vectors.reflected_ax,   &vectors.candidate_ax};
  for (std::size_t index = 0;
       index < sizeof(col_vectors) / sizeof(col_vectors[0]); ++index) {
    const api::Error error =
        selected_backend.createVector(num_cols, *col_vectors[index]);
    if (!error.ok()) return error;
    if (diagnostics != nullptr)
      ++diagnostics->resource_creation_calls_before_iterations;
  }
  for (std::size_t index = 0;
       index < sizeof(row_vectors) / sizeof(row_vectors[0]); ++index) {
    const api::Error error =
        selected_backend.createVector(num_rows, *row_vectors[index]);
    if (!error.ok()) return error;
    if (diagnostics != nullptr)
      ++diagnostics->resource_creation_calls_before_iterations;
  }
  return api::Error();
}

double elapsedSeconds(const SolveClock::time_point& start) {
  return std::chrono::duration<double>(SolveClock::now() - start).count();
}

bool timeLimitReached(const api::SolverOptions& options,
                      const SolveClock::time_point& start) {
  return options.time_limit != kHighsInf &&
         elapsedSeconds(start) >= options.time_limit;
}

api::Error downloadIterate(backend::Backend& selected_backend,
                           const backend::VectorHandle& primal,
                           const backend::VectorHandle& dual,
                           const backend::VectorHandle& ax,
                           const backend::VectorHandle& aty,
                           HostIterate& host, const SolveContext& context,
                           const HighsInt iteration) {
  SolvePhaseGuard download(context, SolvePhase::kDownload, iteration);
  const backend::VectorDownload downloads[] = {{&primal, &host.primal},
                                               {&dual, &host.dual},
                                               {&ax, &host.ax},
                                               {&aty, &host.aty}};
  return selected_backend.downloadVectors(downloads, 4);
}

api::Error prepareCandidateProducts(
    backend::Backend& selected_backend,
    const backend::ProblemHandle& problem_handle,
    const HalpernPdhgState& state, IterationVectors& vectors) {
  api::Error error = selected_backend.multiply(
      problem_handle, *state.primal_candidate, vectors.candidate_ax);
  if (!error.ok()) return error;
  return selected_backend.multiplyTranspose(
      problem_handle, *state.dual_candidate, vectors.candidate_aty);
}

api::Error prepareCandidate(backend::Backend& selected_backend,
                            const backend::ProblemHandle& problem_handle,
                            const HalpernPdhgState& state,
                            IterationVectors& vectors, HostIterate& host,
                            const SolveContext& context,
                            const HighsInt iteration) {
  const api::Error error = prepareCandidateProducts(
      selected_backend, problem_handle, state, vectors);
  if (!error.ok()) return error;
  if (context.diagnostics != nullptr)
    ++context.diagnostics->candidate_downloads;
  return downloadIterate(selected_backend, *state.primal_candidate,
                         *state.dual_candidate, vectors.candidate_ax,
                         vectors.candidate_aty, host, context, iteration);
}

api::Error evaluateResidentCheck(
    backend::Backend& selected_backend, const model::Problem& problem,
    backend::WorkspaceHandle& workspace, const backend::CheckInput& input,
    const bool ordinary_boundary, SolveDiagnostics* diagnostics,
    backend::CheckResult& result) {
  if (diagnostics != nullptr) ++diagnostics->resident_check_calls;
  const api::Error error =
      selected_backend.evaluateCheck(problem, workspace, input, result);
  if (diagnostics != nullptr && ordinary_boundary) {
    diagnostics->ordinary_check_d2h_bytes += result.device_to_host_bytes;
    diagnostics->ordinary_check_transfers += result.device_to_host_transfers;
    diagnostics->ordinary_check_synchronizations += result.synchronizations;
  }
  return error;
}

api::Error postsolveIterate(const model::Problem& problem, HostIterate& host,
                            HighsSolution& solution,
                            const SolveContext& context,
                            const HighsInt iteration) {
  {
    SolvePhaseGuard preparation(context, SolvePhase::kHostPreparation,
                                iteration);
    host.reduced_costs.resize(static_cast<std::size_t>(problem.numCols()));
    for (HighsInt col = 0; col < problem.numCols(); ++col) {
      host.reduced_costs[col] = problem.costs()[col] - host.aty[col];
    }
  }
  SolvePhaseGuard verification(context, SolvePhase::kHostVerification,
                               iteration);
  const model::PostsolveResult postsolved =
      model::postsolve(problem, host.primal, host.dual, host.reduced_costs);
  if (!postsolved.ok()) return postsolved.error;
  solution = postsolved.solution;
  return api::Error();
}

FixedPointErrorInput candidateFixedPointErrorInput(
  const HalpernPdhgState& state,
    backend::VectorHandle& aty_delta, const InitialStepSizes& step_sizes) {
  FixedPointErrorInput input;
  input.next_primal = state.primal_candidate;
  input.next_dual = state.dual_candidate;
  input.reflected_primal = state.reflected_primal;
  input.reflected_dual = state.reflected_dual;
  input.omega = step_sizes.primal_weight;
  input.eta = step_sizes.base_step;
  input.delta_primal = state.next_primal;
  input.delta_dual = state.next_dual;
  input.aty_delta = &aty_delta;
  return input;
}

FixedPointErrorResult evaluateCandidateFixedPointError(
    backend::Backend& selected_backend,
    const backend::ProblemHandle& problem_handle, HalpernPdhgState& state,
    backend::VectorHandle& aty_delta, const InitialStepSizes& step_sizes,
    const SolveContext& context, const HighsInt iteration,
    const SolvePhase phase, const bool products_prepared = false) {
  SolvePhaseGuard evaluation(context, phase, iteration);
  return evaluateFixedPointError(
      selected_backend, problem_handle,
      candidateFixedPointErrorInput(state, aty_delta, step_sizes), products_prepared);
}

api::Error invalidKktResult(
    const verification::KktVerificationResult& verification) {
  const std::string detail = verification.error.empty()
                                 ? "independent KKT verification failed"
                                 : verification.error;
  return api::Error(api::ErrorCode::kInvalidModel, detail);
}

double euclideanNorm(const std::vector<double>& values) {
  double norm = 0.0;
  for (std::size_t index = 0; index < values.size(); ++index)
    norm = std::hypot(norm, values[index]);
  return norm;
}

bool normalizeInPlace(std::vector<double>& values, const double norm) {
  if (!std::isfinite(norm) || norm <= 0.0) return false;
  for (std::size_t index = 0; index < values.size(); ++index) {
    values[index] /= norm;
    if (!std::isfinite(values[index])) return false;
  }
  return true;
}

double largest(const double left, const double right) {
  return left > right ? left : right;
}

// Turns the host candidate into normalized original-space directions and asks
// the independent verifier what, if anything, they prove. A recovery or
// normalization failure yields no certificate at all rather than a partial
// one, so the loop simply keeps iterating.
CertificateOutcome detectCertificate(const SolveContext& context,
                                     const model::Problem& problem,
                                     const HighsLp& lp,
                                     const HostIterate& host,
                                     const double tolerance,
                                     CertificateBuffers& buffers,
                                     const HighsInt iteration) {
                          SolvePhaseGuard verification(context, SolvePhase::kHostVerification,
                                   iteration);
  if (context.diagnostics != nullptr)
    ++context.diagnostics->certificate_checks;
  CertificateOutcome outcome;
  if (lp.num_col_ == 0 && lp.num_row_ == 0) return outcome;

  bool ray_available = false;
  verification::PrimalRayVerificationResult ray;
  if (problem.recoverPrimalDirection(host.primal, buffers.primal_direction)
          .ok()) {
    const double norm = euclideanNorm(buffers.primal_direction);
    if (normalizeInPlace(buffers.primal_direction, norm)) {
      ray = verification::verifyPrimalRay(lp, buffers.primal_direction,
                                          tolerance);
      ray.metrics.direction_norm = norm;
      ray_available = ray.pass;
    }
  }

  bool farkas_available = false;
  verification::FarkasVerificationResult farkas;
  buffers.scaled_column_multiplier.resize(host.aty.size());
  for (std::size_t col = 0; col < host.aty.size(); ++col)
    buffers.scaled_column_multiplier[col] = -host.aty[col];
  if (problem.recoverDualDirection(host.dual, buffers.dual_direction).ok() &&
      problem
          .recoverReducedCostDirection(buffers.scaled_column_multiplier,
                                       buffers.column_multiplier)
          .ok()) {
    const double norm = std::hypot(euclideanNorm(buffers.dual_direction),
                                   euclideanNorm(buffers.column_multiplier));
    if (normalizeInPlace(buffers.dual_direction, norm) &&
        normalizeInPlace(buffers.column_multiplier, norm)) {
      farkas = verification::verifyFarkasCertificate(
          lp, buffers.dual_direction, buffers.column_multiplier, tolerance);
      farkas.metrics.direction_norm = norm;
      farkas_available = farkas.pass;
    }
  }

  if (!ray_available && !farkas_available) return outcome;
  outcome.available = true;
  outcome.certificate.valid = true;
  if (ray_available && farkas_available) {
    // Both conclusions hold, so the two cannot be told apart. The stored pair
    // is the Farkas certificate and the margin is the weaker of the two.
    outcome.status = api::SolverStatus::kUnboundedOrInfeasible;
    outcome.certificate.kind = api::CertificateKind::kAmbiguous;
    outcome.certificate.row_direction = buffers.dual_direction;
    outcome.certificate.column_direction = buffers.column_multiplier;
    outcome.certificate.metrics.direction_norm =
        farkas.metrics.direction_norm;
    outcome.certificate.metrics.max_normalized_violation = largest(
        largest(farkas.metrics.max_normalized_row_sign_violation,
                farkas.metrics.max_normalized_column_sign_violation),
        largest(ray.metrics.max_normalized_column_violation,
                ray.metrics.max_normalized_row_violation));
    outcome.certificate.metrics.max_normalized_residual =
        farkas.metrics.max_normalized_stationarity_residual;
    outcome.certificate.metrics.margin =
        farkas.metrics.certificate_margin < ray.metrics.objective_margin
            ? farkas.metrics.certificate_margin
            : ray.metrics.objective_margin;
    outcome.certificate.metrics.worst_row = farkas.metrics.worst_row;
    outcome.certificate.metrics.worst_col = farkas.metrics.worst_col;
  } else if (farkas_available) {
    outcome.status = api::SolverStatus::kInfeasible;
    outcome.certificate.kind = api::CertificateKind::kPrimalInfeasibility;
    outcome.certificate.row_direction = buffers.dual_direction;
    outcome.certificate.column_direction = buffers.column_multiplier;
    outcome.certificate.metrics.direction_norm =
        farkas.metrics.direction_norm;
    outcome.certificate.metrics.max_normalized_violation =
        largest(farkas.metrics.max_normalized_row_sign_violation,
                farkas.metrics.max_normalized_column_sign_violation);
    outcome.certificate.metrics.max_normalized_residual =
        farkas.metrics.max_normalized_stationarity_residual;
    outcome.certificate.metrics.margin = farkas.metrics.certificate_margin;
    outcome.certificate.metrics.worst_row = farkas.metrics.worst_row;
    outcome.certificate.metrics.worst_col = farkas.metrics.worst_col;
  } else {
    outcome.status = api::SolverStatus::kUnbounded;
    outcome.certificate.kind = api::CertificateKind::kPrimalUnboundedness;
    outcome.certificate.column_direction = buffers.primal_direction;
    outcome.certificate.row_direction.assign(
        static_cast<std::size_t>(lp.num_row_), 0.0);
    outcome.certificate.metrics.direction_norm = ray.metrics.direction_norm;
    outcome.certificate.metrics.max_normalized_violation =
        largest(ray.metrics.max_normalized_column_violation,
                ray.metrics.max_normalized_row_violation);
    outcome.certificate.metrics.margin = ray.metrics.objective_margin;
    outcome.certificate.metrics.worst_row = ray.metrics.worst_row;
    outcome.certificate.metrics.worst_col = ray.metrics.worst_col;
  }
  return outcome;
}

// The registry owns backend names; the control loop never branches on kinds.
const char* requestedBackendName(const backend::Kind kind) {
  const backend::Descriptor* registration = backend::findRegistration(kind);
  return registration != nullptr ? registration->name : "auto";
}

const char* scalingName(const model::ScalingMode mode) {
  switch (mode) {
    case model::ScalingMode::kAuto:
      return "auto";
    case model::ScalingMode::kOff:
      return "off";
    case model::ScalingMode::kCurtisReid:
      return "curtis_reid";
    case model::ScalingMode::kLegacy:
      return "legacy";
  }
  return "unknown";
}

const char* statusName(const api::SolverStatus status) {
  switch (status) {
    case api::SolverStatus::kOptimal:
      return "optimal";
    case api::SolverStatus::kInfeasible:
      return "infeasible";
    case api::SolverStatus::kUnbounded:
      return "unbounded";
    case api::SolverStatus::kUnboundedOrInfeasible:
      return "unbounded or infeasible";
    case api::SolverStatus::kIterationLimit:
      return "iteration limit";
    case api::SolverStatus::kTimeLimit:
      return "time limit";
    case api::SolverStatus::kError:
      return "error";
    case api::SolverStatus::kNotSet:
      return "not set";
  }
  return "unknown";
}

const char* certificateKindName(const api::CertificateKind kind) {
  switch (kind) {
    case api::CertificateKind::kNone:
      return "none";
    case api::CertificateKind::kPrimalInfeasibility:
      return "primal infeasibility (Farkas)";
    case api::CertificateKind::kPrimalUnboundedness:
      return "primal unboundedness (recession ray)";
    case api::CertificateKind::kAmbiguous:
      return "ambiguous";
  }
  return "unknown";
}

std::string formatted(const char* format, ...) {
  char buffer[512];
  va_list arguments;
  va_start(arguments, format);
  vsnprintf(buffer, sizeof(buffer), format, arguments);
  va_end(arguments);
  return std::string(buffer);
}

const char* restartConditionName(const FixedPointRestartCondition condition) {
  switch (condition) {
    case FixedPointRestartCondition::kNone:
      return "none";
    case FixedPointRestartCondition::kNotTriggered:
      return "not-triggered";
    case FixedPointRestartCondition::kFirstCheckInterval:
      return "first-check-interval";
    case FixedPointRestartCondition::kSufficientDecay:
      return "sufficient-decay";
    case FixedPointRestartCondition::kNecessaryDecay:
      return "necessary-decay";
    case FixedPointRestartCondition::kArtificialThreshold:
      return "artificial-threshold";
  }
  return "unknown";
}

const char* pidOutcomeName(const bool use_pid,
                           const PrimalWeightUpdateKind kind) {
  if (!use_pid) return "disabled";
  switch (kind) {
    case PrimalWeightUpdateKind::kNone:
      return "not-run";
    case PrimalWeightUpdateKind::kControlled:
      return "controlled";
    case PrimalWeightUpdateKind::kRevertedToBest:
      return "reverted-to-best";
  }
  return "unknown";
}

void logBoundary(const SolveContext& context, const HighsInt iteration,
                 const NormalizedConvergenceMetrics& metrics,
                 const InitialStepSizes& steps,
                 const FixedPointRestartCondition restart_condition,
                 const bool use_pid,
                 const PrimalWeightUpdateKind pid_outcome) {
  if (context.log == nullptr) return;
  context.log->info(formatted(
      "HiPDLP rewrite boundary: iteration %lld, relative primal %.17g, dual "
      "%.17g, gap %.17g, tau %.17g, sigma %.17g, omega %.17g, eta %.17g, "
      "restart condition %s, PID outcome %s",
      static_cast<long long>(iteration), metrics.relative_primal,
      metrics.relative_dual, metrics.relative_gap, steps.primal_step,
      steps.dual_step, steps.primal_weight, steps.base_step,
      restartConditionName(restart_condition),
      pidOutcomeName(use_pid, pid_outcome)));
}

// Reports the requested and effective policies exactly once per solve.
void logStartup(const SolveContext& context, const api::SolverOptions& options,
                backend::Backend& selected_backend) {
  if (context.log == nullptr) return;
  backend::DeviceInfo device;
  const api::Error device_error = selected_backend.queryDevice(device);
  context.log->info(formatted(
      "HiPDLP rewrite backend: requested %s, effective %s (%s)",
      requestedBackendName(options.backend), selected_backend.descriptor().name,
      device_error.ok() && device.available ? device.name.c_str()
                                            : "device unavailable"));
  if (options.policy.step_size_substituted) {
    context.log->warning(formatted(
        "HiPDLP rewrite step size: requested pdlp_step_size_strategy=%d "
        "(adaptive), effective fixed; the adaptive strategy is not implemented",
        static_cast<int>(options.policy.requested_step_size_strategy)));
  } else {
    context.log->info(formatted(
        "HiPDLP rewrite step size: requested pdlp_step_size_strategy=%d, "
        "effective %s",
        static_cast<int>(options.policy.requested_step_size_strategy),
        options.policy.step_size == api::StepSizePolicy::kPid ? "PID"
                                                              : "fixed"));
  }
  if (options.policy.restart_substituted) {
    context.log->warning(formatted(
        "HiPDLP rewrite restart: requested pdlp_restart_strategy=%d, effective "
        "Halpern fixed-point restart; no other restart strategy is implemented",
        static_cast<int>(options.policy.requested_restart_strategy)));
  } else {
    context.log->info(
        "HiPDLP rewrite restart: effective Halpern fixed-point restart");
  }
  context.log->info(formatted("HiPDLP rewrite scaling: effective %s",
                              scalingName(options.scaling.mode)));
}

void logCompletion(const SolveContext& context, const HighsLp& lp,
                   const api::SolverStatus status,
                   const HighsSolution& solution,
                   const HighsInt iteration_count, const double elapsed,
                   const double kkt_tolerance) {
  if (context.log == nullptr) return;
  const verification::KktVerificationResult verification =
      verification::verifyKkt(lp, solution, kkt_tolerance);
  context.log->info(formatted(
      "HiPDLP rewrite finished: status %s, %d iterations, %.4gs, relative KKT "
      "primal %.3e, dual %.3e, gap %.3e, complementarity %.3e",
      statusName(status), static_cast<int>(iteration_count), elapsed,
      verification.metrics.relative_primal_residual,
      verification.metrics.relative_dual_residual,
      verification.metrics.relative_objective_gap,
      verification.metrics.relative_complementarity));
}

// A certificate status never reports optimality metrics it did not establish.
void logCertificateCompletion(const SolveContext& context,
                              const api::SolverStatus status,
                              const api::Certificate& certificate,
                              const HighsInt iteration_count,
                              const double elapsed) {
  if (context.log == nullptr) return;
  context.log->info(formatted(
      "HiPDLP rewrite finished: status %s, %d iterations, %.4gs, certificate "
      "%s, normalized violation %.3e, residual %.3e, margin %.3e, worst row "
      "%d, worst column %d",
      statusName(status), static_cast<int>(iteration_count), elapsed,
      certificateKindName(certificate.kind),
      certificate.metrics.max_normalized_violation,
      certificate.metrics.max_normalized_residual, certificate.metrics.margin,
      static_cast<int>(certificate.metrics.worst_row),
      static_cast<int>(certificate.metrics.worst_col)));
}

api::SolveResult certifiedResult(const api::SolverStatus status,
                                 const HighsSolution& solution,
                                 const HighsInt iteration_count,
                                 const api::Certificate& certificate) {
  api::SolveResult result =
      api::SolveResult::success(status, solution, iteration_count);
  result.certificate = certificate;
  return result;
}

}  // namespace

SolveDiagnostics::SolveDiagnostics()
    : resource_creation_calls_before_iterations(0),
      resource_creation_calls_during_iterations(0),
      primal_weight_updates(0),
      restarts(0),
      check_boundaries(0),
      resident_check_calls(0),
      ordinary_check_boundaries(0),
      ordinary_check_d2h_bytes(0),
      ordinary_check_transfers(0),
      ordinary_check_synchronizations(0),
      candidate_downloads(0),
      host_oracle_checks(0),
      certificate_checks(0),
      original_kkt_verifications(0),
      original_kkt_rejections(0),
      restart_reference_checks(0),
      fixed_point_evaluations(0),
      command_buffer_replays(0),
      recorded_iterations(0),
      total_wall_seconds(0.0),
      preprocessing_seconds(0.0),
      backend_setup_seconds(0.0),
      operator_norm_seconds(0.0),
      iterations_seconds(0.0),
      iteration_block_seconds(0.0),
      step_seconds(0.0),
      boundary_seconds(0.0),
      download_seconds(0.0),
      host_preparation_seconds(0.0),
      host_verification_seconds(0.0),
      reference_fixed_point_seconds(0.0),
      fixed_point_fallback_seconds(0.0) {}

SolveLogSink::~SolveLogSink() {}

SolveContext::SolveContext()
    : log(nullptr), diagnostics(nullptr), observer(nullptr),
      use_command_buffers(true) {}

api::SolveResult solve(const api::SolverOptions& options,
                       const model::OriginalModel& original_model,
                       const SolveContext& context) {
  HighsInt total_iterations = 0;
  SolvePhaseGuard solve_phase(context, SolvePhase::kSolve, total_iterations);
  const SolveClock::time_point start = SolveClock::now();
  api::Error error = validateOptions(options);
  if (!error.ok()) return failure(error);

  try {
    SolvePhaseGuard backend_setup(context, SolvePhase::kBackendSetup,
                                  total_iterations);
    // Backend availability is resolved first and never falls back.
    backend::FactoryResult factory = backend::createBackend(options.backend);
    if (!factory.ok()) return failure(factory.error);
    std::unique_ptr<backend::Backend> selected_backend =
        std::move(factory.backend);
    error = requireSolveCapabilities(*selected_backend);
    if (!error.ok()) return failure(error);
    logStartup(context, options, *selected_backend);
    backend_setup.finish();

    SolvePhaseGuard preprocessing(context, SolvePhase::kPreprocessing,
                    total_iterations);
    model::ProblemBuildResult built = model::buildProblem(original_model.lp());
    if (!built.ok()) return failure(built.error);
    error = model::applyScalingPipeline(built.problem, options.scaling);
    if (!error.ok()) return failure(error);
    if (context.log != nullptr && built.problem.scaling().curtis_reid_rejected)
      context.log->info(
          "HiPDLP rewrite scaling: Curtis-Reid rejected, it would inflate the "
          "scaled costs, right-hand sides or bounds; using Ruiz and "
          "Pock-Chambolle");
    if (built.problem.matrix().numNonzeros() == 0) {
      return failure(
          api::Error(api::ErrorCode::kInvalidModel,
                     "HiPDLP rewrite does not yet support zero-operator LPs"));
    }
            preprocessing.finish();

            backend_setup.begin();
    backend::ProblemHandle problem_handle;
    error = selected_backend->bindProblem(built.problem, problem_handle);
    if (!error.ok()) return failure(error);
    if (context.diagnostics != nullptr)
      ++context.diagnostics->resource_creation_calls_before_iterations;
    backend::WorkspaceHandle workspace;
    error = selected_backend->createWorkspace(problem_handle, workspace);
    if (!error.ok()) return failure(error);
    if (context.diagnostics != nullptr)
      ++context.diagnostics->resource_creation_calls_before_iterations;

    IterationVectors vectors;
    error = createIterationVectors(*selected_backend, built.problem.numRows(),
                                   built.problem.numCols(), vectors,
                                   context.diagnostics);
    if (!error.ok()) return failure(error);
    error = selected_backend->initialize(vectors.current_primal,
                                         vectors.current_dual);
    if (!error.ok()) return failure(error);
    error = selected_backend->fill(vectors.current_aty, 0.0);
    if (!error.ok()) return failure(error);
    error = selected_backend->fusedProjectedPrimalUpdate(
        problem_handle, workspace, vectors.current_primal, vectors.current_aty,
        0.0, vectors.next_primal);
    if (!error.ok()) return failure(error);
    std::swap(vectors.current_primal, vectors.next_primal);
    error =
        selected_backend->copy(vectors.current_primal, vectors.primal_anchor);
    if (!error.ok()) return failure(error);
    error = selected_backend->copy(vectors.current_dual, vectors.dual_anchor);
    if (!error.ok()) return failure(error);
    error = selected_backend->multiply(problem_handle, vectors.current_primal,
                                       vectors.current_ax);
    if (!error.ok()) return failure(error);
    error = selected_backend->multiplyTranspose(
        problem_handle, vectors.current_dual, vectors.current_aty);
    if (!error.ok()) return failure(error);
    backend_setup.finish();

    ConvergenceTolerances tolerances;
    tolerances.primal = options.optimality_tolerance;
    tolerances.dual = options.optimality_tolerance;
    tolerances.gap = options.optimality_tolerance;

    HostIterate host;
    HighsSolution solution;
    CertificateBuffers certificate_buffers;
    NormalizedConvergenceMetrics boundary_metrics;
    bool boundary_metrics_valid = false;

    // Iteration 0: the starting point is checked like any later candidate.
    {
      SolvePhaseGuard boundary(context, SolvePhase::kBoundary,
                               total_iterations);
      error = downloadIterate(*selected_backend, vectors.current_primal,
                              vectors.current_dual, vectors.current_ax,
                              vectors.current_aty, host, context,
                              total_iterations);
      if (!error.ok()) return failure(error);
      RawConvergenceMetrics raw;
      if (context.diagnostics != nullptr) {
        ++context.diagnostics->check_boundaries;
        ++context.diagnostics->host_oracle_checks;
      }
      {
        SolvePhaseGuard verification(context, SolvePhase::kHostVerification,
                                     total_iterations);
        error = computeRawConvergenceMetrics(built.problem, host.primal,
                                             host.dual, host.ax, host.aty, raw);
      }
      if (!error.ok()) return failure(error);
      StoppingInput stopping_input;
      stopping_input.raw = raw;
      stopping_input.tolerances = tolerances;
      stopping_input.iteration = 0;
      stopping_input.elapsed_time = elapsedSeconds(start);
      stopping_input.limits.iteration = options.iteration_limit;
      stopping_input.limits.time = options.time_limit;
      const StoppingResult stopping = evaluateStopping(stopping_input);
      if (!stopping.ok()) return failure(stopping.error);

      bool postsolved = false;
      if (stopping.decision == StoppingDecision::kOptimal) {
        error = postsolveIterate(built.problem, host, solution, context,
                                total_iterations);
        if (!error.ok()) return failure(error);
        postsolved = true;
        SolvePhaseGuard host_verification(
          context, SolvePhase::kHostVerification, total_iterations);
        const verification::KktVerificationResult verification =
            verification::verifyKkt(original_model.lp(), solution,
                                    options.kkt_tolerance);
        if (!verification.input_valid || !verification.metrics_finite)
          return failure(invalidKktResult(verification));
        if (verification.pass) {
          logCompletion(context, original_model.lp(),
                        api::SolverStatus::kOptimal, solution, 0,
                        elapsedSeconds(start), options.kkt_tolerance);
          return api::SolveResult::success(api::SolverStatus::kOptimal,
                                           solution, 0);
        }
      }
        const bool initial_iteration_limit = options.iteration_limit == 0;
        if (initial_iteration_limit || timeLimitReached(options, start) ||
          stopping.decision == StoppingDecision::kIterationLimit ||
          stopping.decision == StoppingDecision::kTimeLimit) {
        if (!postsolved) {
          error = postsolveIterate(built.problem, host, solution, context,
                                  total_iterations);
          if (!error.ok()) return failure(error);
        }
        const CertificateOutcome outcome = detectCertificate(
          context, built.problem, original_model.lp(), host,
          options.kkt_tolerance,
            certificate_buffers, total_iterations);
        if (outcome.available) {
          logCertificateCompletion(context, outcome.status,
                                   outcome.certificate, 0,
                                   elapsedSeconds(start));
          return certifiedResult(outcome.status, solution, 0,
                                 outcome.certificate);
        }
        const api::SolverStatus status =
          initial_iteration_limit ||
              stopping.decision == StoppingDecision::kIterationLimit
                ? api::SolverStatus::kIterationLimit
                : api::SolverStatus::kTimeLimit;
        logCompletion(context, original_model.lp(), status, solution, 0,
                      elapsedSeconds(start), options.kkt_tolerance);
        return api::SolveResult::success(status, solution, 0);
      }

      // Iteration 0 is a convergence boundary too, so it gets the same
      // certificate treatment as every later candidate.
      const CertificateOutcome outcome = detectCertificate(
          context, built.problem, original_model.lp(), host,
          options.kkt_tolerance,
          certificate_buffers, total_iterations);
      if (outcome.available) {
        if (!postsolved) {
          error = postsolveIterate(built.problem, host, solution, context,
                                  total_iterations);
          if (!error.ok()) return failure(error);
        }
        logCertificateCompletion(context, outcome.status, outcome.certificate,
                                 0, elapsedSeconds(start));
        return certifiedResult(outcome.status, solution, 0,
                               outcome.certificate);
      }
      boundary.finish(true);
    }

    // The remaining budget is checked on both sides of the power iteration.
    if (timeLimitReached(options, start)) {
      error = postsolveIterate(built.problem, host, solution, context,
                              total_iterations);
      if (!error.ok()) return failure(error);
      logCompletion(context, original_model.lp(), api::SolverStatus::kTimeLimit,
                    solution, 0, elapsedSeconds(start), options.kkt_tolerance);
      return api::SolveResult::success(api::SolverStatus::kTimeLimit, solution,
                                       0);
    }
      SolvePhaseGuard norm_phase(context, SolvePhase::kOperatorNorm,
                     total_iterations);
    const OperatorNormResult operator_norm = estimateOperatorNormSquared(
        *selected_backend, problem_handle, built.problem.numRows(),
        built.problem.numCols());
    if (!operator_norm.ok()) return failure(operator_norm.error);
    if (operator_norm.empty_matrix ||
        !std::isfinite(operator_norm.operator_norm_squared) ||
        operator_norm.operator_norm_squared <= 0.0) {
      return failure(
          api::Error(api::ErrorCode::kInvalidModel,
                     "HiPDLP rewrite does not yet support zero-operator LPs"));
    }
            norm_phase.finish();
    if (timeLimitReached(options, start)) {
      error = postsolveIterate(built.problem, host, solution, context,
                              total_iterations);
      if (!error.ok()) return failure(error);
      logCompletion(context, original_model.lp(), api::SolverStatus::kTimeLimit,
                    solution, 0, elapsedSeconds(start), options.kkt_tolerance);
      return api::SolveResult::success(api::SolverStatus::kTimeLimit, solution,
                                       0);
    }

                      backend_setup.begin();
    StepSizeInitializationInput step_input;
    double operator_norm_upper_bound_squared = 0.0;
    const std::vector<model::ScalingStage>& scaling_stages =
        built.problem.scaling().stages;
    if (!scaling_stages.empty() &&
        scaling_stages.back() == model::ScalingStage::kPockChambolle) {
      // Pock-Chambolle row/column sums give the Schur bound ||A||_2 <= 1.
      operator_norm_upper_bound_squared = 1.0;
    } else {
      error = computeOperatorNormSquaredUpperBound(
          built.problem.matrix(), operator_norm_upper_bound_squared);
      if (!error.ok()) return failure(error);
    }
    step_input.operator_norm_squared =
        operator_norm.operator_norm_squared > operator_norm_upper_bound_squared
            ? operator_norm.operator_norm_squared
            : operator_norm_upper_bound_squared;
    step_input.unscaled_rhs_norm = built.problem.unscaledNorms().row_bounds_l2;
    step_input.unscaled_cost_norm = built.problem.unscaledNorms().objective_l2;
    const StepSizeInitializationResult initialized_steps =
        initializeStepSizes(step_input);
    if (!initialized_steps.ok()) return failure(initialized_steps.error);

    HalpernPdhgState state;
    state.problem = &problem_handle;
    state.workspace = &workspace;
    state.current_primal = &vectors.current_primal;
    state.current_dual = &vectors.current_dual;
    state.next_primal = &vectors.next_primal;
    state.next_dual = &vectors.next_dual;
    state.primal_anchor = &vectors.primal_anchor;
    state.dual_anchor = &vectors.dual_anchor;
    state.primal_candidate = &vectors.primal_candidate;
    state.dual_candidate = &vectors.dual_candidate;
    state.reflected_primal = &vectors.reflected_primal;
    state.reflected_dual = &vectors.reflected_dual;
    state.current_aty = &vectors.current_aty;
    state.reflected_ax = &vectors.reflected_ax;

    // omega follows the PID controller; eta stays the initialized base step.
    InitialStepSizes current_steps = initialized_steps.step_sizes;
    state.primal_step_size = current_steps.primal_step;
    state.dual_step_size = current_steps.dual_step;
    const bool use_pid = options.policy.step_size == api::StepSizePolicy::kPid;
    PrimalWeightState primal_weight_state =
        resetPrimalWeight(current_steps.primal_weight);
    const PrimalWeightParameters primal_weight_parameters;

    FixedPointRestartState restart_state = resetFixedPointRestart(0.0);
    const FixedPointRestartParameters restart_parameters;
    bool restart_reference_pending = false;
    const bool back_off_verifications =
        iteratesOnDevice(selected_backend->descriptor().kind);
    HighsInt rejected_verifications = 0;
    HighsInt next_verification_iteration = 0;
    error = selected_backend->prepareCheckWorkspace(
      built.problem, original_model.lp(), workspace);
    if (!error.ok()) return failure(error);
    const backend::VectorDownload candidate_downloads[] = {
      {state.primal_candidate, &host.primal},
      {state.dual_candidate, &host.dual},
      {&vectors.candidate_ax, &host.ax},
      {&vectors.candidate_aty, &host.aty}};
    error = selected_backend->prepareVectorDownloads(candidate_downloads, 4);
    if (!error.ok()) return failure(error);
    host.reduced_costs.resize(static_cast<std::size_t>(built.problem.numCols()));
    certificate_buffers.scaled_column_multiplier.reserve(
      static_cast<std::size_t>(built.problem.numCols()));
    certificate_buffers.primal_direction.reserve(
      static_cast<std::size_t>(original_model.lp().num_col_));
    certificate_buffers.dual_direction.reserve(
      static_cast<std::size_t>(original_model.lp().num_row_));
    certificate_buffers.column_multiplier.reserve(
      static_cast<std::size_t>(original_model.lp().num_col_));
    RecordedIterationBlock recorded_block;
    if (context.use_command_buffers &&
        (options.iteration_limit == kHighsIInf ||
         options.iteration_limit >= kConvergenceCadence)) {
      error = recorded_block.prepare(*selected_backend, state);
      if (!error.ok()) return failure(error);
      if (recorded_block.ready() && context.diagnostics != nullptr)
        context.diagnostics->resource_creation_calls_before_iterations += 2;
    }
    backend_setup.finish();
    SolvePhaseGuard iterations(context, SolvePhase::kIterations,
                               total_iterations);
    SolvePhaseGuard iteration_block(context, SolvePhase::kIterationBlock,
                                    total_iterations);
    while (true) {
      if (timeLimitReached(options, start)) {
        SolvePhaseGuard boundary(context, SolvePhase::kBoundary,
                                 total_iterations);
        if (total_iterations == 0) {
          error = downloadIterate(*selected_backend, *state.current_primal,
                                  *state.current_dual, vectors.current_ax,
                                  *state.current_aty, host, context,
                                  total_iterations);
        } else {
          error = prepareCandidate(*selected_backend, problem_handle, state,
                                   vectors, host, context, total_iterations);
        }
        if (!error.ok()) return failure(error);
        error = postsolveIterate(built.problem, host, solution, context,
              total_iterations);
        if (!error.ok()) return failure(error);
        const CertificateOutcome outcome = detectCertificate(
          context, built.problem, original_model.lp(), host,
          options.kkt_tolerance,
            certificate_buffers, total_iterations);
        if (outcome.available) {
          logCertificateCompletion(context, outcome.status,
                                   outcome.certificate, total_iterations,
                                   elapsedSeconds(start));
          return certifiedResult(outcome.status, solution, total_iterations,
                                 outcome.certificate);
        }
        logCompletion(context, original_model.lp(),
                      api::SolverStatus::kTimeLimit, solution, total_iterations,
                      elapsedSeconds(start), options.kkt_tolerance);
        return api::SolveResult::success(api::SolverStatus::kTimeLimit,
                                         solution, total_iterations);
      }

      {
        SolvePhaseGuard step(context, SolvePhase::kStep, total_iterations);
        const bool replay_block = recorded_block.ready() &&
            total_iterations % kConvergenceCadence == 1 &&
            !restart_reference_pending &&
            total_iterations <= std::numeric_limits<HighsInt>::max() - kRecordedStepCount &&
            state.iteration <= std::numeric_limits<HighsInt>::max() - kRecordedStepCount &&
            (options.iteration_limit == kHighsIInf ||
             options.iteration_limit - total_iterations >= kRecordedStepCount);
        const HighsInt completed_steps = replay_block ? kRecordedStepCount : 1;
        error = replay_block
            ? recorded_block.replay(*selected_backend, state)
            : enqueueHalpernPdhgStep(*selected_backend, state);
        if (!error.ok()) return failure(error);
        if (total_iterations > std::numeric_limits<HighsInt>::max() - completed_steps) {
          return failure(api::Error(api::ErrorCode::kInvalidArgument,
                                    "HiPDLP rewrite iteration count overflowed"));
        }
        total_iterations += completed_steps;
        if (replay_block && context.diagnostics != nullptr) {
          ++context.diagnostics->command_buffer_replays;
          context.diagnostics->recorded_iterations += completed_steps;
        }
      }

      const bool cadence_boundary = total_iterations % kConvergenceCadence == 0;
      const bool iteration_limit = options.iteration_limit != kHighsIInf &&
                                   total_iterations >= options.iteration_limit;
      const bool time_limit = timeLimitReached(options, start);
      boundary_metrics_valid = false;
      backend::CheckResult check;
      bool fixed_point_from_packet = false;
      SolvePhaseGuard boundary(context, SolvePhase::kBoundary,
                               total_iterations, false);
      if (cadence_boundary || iteration_limit || time_limit) {
        boundary.begin();
        if (context.diagnostics != nullptr) {
          ++context.diagnostics->check_boundaries;
          if (cadence_boundary)
            ++context.diagnostics->ordinary_check_boundaries;
        }
        error = prepareCandidateProducts(*selected_backend, problem_handle,
                                         state, vectors);
        if (!error.ok()) return failure(error);

        backend::CheckInput check_input;
        check_input.primal = state.primal_candidate;
        check_input.dual = state.dual_candidate;
        check_input.ax = &vectors.candidate_ax;
        check_input.aty = &vectors.candidate_aty;
        check_input.certificate_tolerance = options.kkt_tolerance;
        check_input.primal_tolerance = options.optimality_tolerance;
        check_input.dual_tolerance = options.optimality_tolerance;
        check_input.gap_tolerance = options.optimality_tolerance;
        if (cadence_boundary) {
          const FixedPointErrorInput fixed_point_input =
              candidateFixedPointErrorInput(state, vectors.aty_delta,
                                            current_steps);
            error = selected_backend->prepareCheckFixedPoint(
              problem_handle, workspace, *fixed_point_input.next_primal,
              *fixed_point_input.next_dual, *fixed_point_input.reflected_primal,
              *fixed_point_input.reflected_dual, *fixed_point_input.delta_primal,
              *fixed_point_input.delta_dual, *fixed_point_input.aty_delta);
          if (error.ok()) {
            check_input.fixed_point = true;
            check_input.delta_primal = fixed_point_input.delta_primal;
            check_input.delta_dual = fixed_point_input.delta_dual;
            check_input.aty_delta = fixed_point_input.aty_delta;
            check_input.anchor_distances = use_pid;
            check_input.primal_anchor = state.primal_anchor;
            check_input.dual_anchor = state.dual_anchor;
          }
        }
        error = evaluateResidentCheck(
            *selected_backend, built.problem, workspace, check_input,
            cadence_boundary, context.diagnostics, check);
        if (!error.ok() &&
            (check_input.fixed_point || check_input.anchor_distances)) {
          check_input.fixed_point = false;
          check_input.anchor_distances = false;
          error = evaluateResidentCheck(
              *selected_backend, built.problem, workspace, check_input,
              cadence_boundary, context.diagnostics, check);
        }
        if (!error.ok()) return failure(error);
        if (!check.valid && !check.requires_host_oracle) {
          return failure(api::Error(
              api::ErrorCode::kInternal,
              "HiPDLP rewrite check returned no valid metrics or host oracle"));
        }
        fixed_point_from_packet = check_input.fixed_point && check.fixed_point_valid;

        bool candidate_downloaded = false;
        const auto download_candidate = [&]() -> api::Error {
          if (candidate_downloaded) return api::Error();
          if (context.diagnostics != nullptr)
            ++context.diagnostics->candidate_downloads;
          const api::Error download_error = downloadIterate(
              *selected_backend, *state.primal_candidate,
              *state.dual_candidate, vectors.candidate_ax,
              vectors.candidate_aty, host, context, total_iterations);
          candidate_downloaded = download_error.ok();
          return download_error;
        };

        RawConvergenceMetrics raw = check.raw;
        if (check.requires_host_oracle) {
          error = download_candidate();
          if (!error.ok()) return failure(error);
          if (context.diagnostics != nullptr)
            ++context.diagnostics->host_oracle_checks;
          SolvePhaseGuard verification(context, SolvePhase::kHostVerification,
                                       total_iterations);
          error = computeRawConvergenceMetrics(
              built.problem, host.primal, host.dual, host.ax, host.aty, raw);
          if (!error.ok()) return failure(error);
        }
        StoppingInput stopping_input;
        stopping_input.raw = raw;
        stopping_input.tolerances = tolerances;
        stopping_input.iteration = total_iterations;
        stopping_input.elapsed_time = elapsedSeconds(start);
        stopping_input.limits.iteration = options.iteration_limit;
        stopping_input.limits.time = options.time_limit;
        const StoppingResult stopping = evaluateStopping(stopping_input);
        if (!stopping.ok()) return failure(stopping.error);
        boundary_metrics = stopping.metrics;
        boundary_metrics_valid = true;

        bool postsolved = false;
        const auto postsolve_candidate = [&]() -> api::Error {
          if (postsolved) return api::Error();
          const api::Error download_error = download_candidate();
          if (!download_error.ok()) return download_error;
          const api::Error postsolve_error =
              postsolveIterate(built.problem, host, solution, context,
                      total_iterations);
          postsolved = postsolve_error.ok();
          return postsolve_error;
        };
        const bool certificate_possible =
            check.primal_ray_possible || check.farkas_possible;
        // A terminal boundary always verifies, so a limit never hides an
        // acceptable candidate behind the backoff.
        const bool verification_due =
            iteration_limit || time_limit ||
            total_iterations >= next_verification_iteration;
        if (verification_due &&
            (stopping.decision == StoppingDecision::kOptimal ||
             (options.early_kkt_check &&
              shouldVerifyOriginalKkt(stopping, total_iterations, tolerances,
                                      options.kkt_tolerance)))) {
          error = postsolve_candidate();
          if (!error.ok()) return failure(error);
          SolvePhaseGuard host_verification(
              context, SolvePhase::kHostVerification, total_iterations);
          if (context.diagnostics != nullptr)
            ++context.diagnostics->original_kkt_verifications;
          const verification::KktVerificationResult verification =
              verification::verifyKkt(original_model.lp(), solution,
                                      options.kkt_tolerance);
          if (!verification.input_valid || !verification.metrics_finite)
            return failure(invalidKktResult(verification));
          if (verification.pass) {
            logCompletion(context, original_model.lp(),
                          api::SolverStatus::kOptimal, solution,
                          total_iterations, elapsedSeconds(start),
                          options.kkt_tolerance);
            return api::SolveResult::success(api::SolverStatus::kOptimal,
                                             solution, total_iterations);
          }
          if (context.diagnostics != nullptr)
            ++context.diagnostics->original_kkt_rejections;
          // The transformed metrics can pass long before the original ones;
          // wait 1, 2, 4, ... boundaries before paying for another check.
          if (back_off_verifications) {
            const HighsInt gap =
                rejected_verifications >= 4
                    ? kMaximumVerificationGap
                    : std::min(kMaximumVerificationGap,
                               HighsInt(1) << rejected_verifications);
            ++rejected_verifications;
            next_verification_iteration =
                total_iterations + gap * kConvergenceCadence;
          }
        }

        if (iteration_limit ||
            stopping.decision == StoppingDecision::kIterationLimit) {
          error = postsolve_candidate();
          if (!error.ok()) return failure(error);
          CertificateOutcome outcome;
          if (certificate_possible)
            outcome = detectCertificate(
                context, built.problem, original_model.lp(), host,
              options.kkt_tolerance, certificate_buffers, total_iterations);
          if (outcome.available) {
            logCertificateCompletion(context, outcome.status,
                                     outcome.certificate, total_iterations,
                                     elapsedSeconds(start));
            return certifiedResult(outcome.status, solution, total_iterations,
                                   outcome.certificate);
          }
          logCompletion(context, original_model.lp(),
                        api::SolverStatus::kIterationLimit, solution,
                        total_iterations, elapsedSeconds(start),
                        options.kkt_tolerance);
          return api::SolveResult::success(api::SolverStatus::kIterationLimit,
                                           solution, total_iterations);
        }
        if (time_limit || stopping.decision == StoppingDecision::kTimeLimit ||
            timeLimitReached(options, start)) {
          error = postsolve_candidate();
          if (!error.ok()) return failure(error);
          CertificateOutcome outcome;
          if (certificate_possible)
            outcome = detectCertificate(
                context, built.problem, original_model.lp(), host,
              options.kkt_tolerance, certificate_buffers, total_iterations);
          if (outcome.available) {
            logCertificateCompletion(context, outcome.status,
                                     outcome.certificate, total_iterations,
                                     elapsedSeconds(start));
            return certifiedResult(outcome.status, solution, total_iterations,
                                   outcome.certificate);
          }
          logCompletion(context, original_model.lp(),
                        api::SolverStatus::kTimeLimit, solution,
                        total_iterations, elapsedSeconds(start),
                        options.kkt_tolerance);
          return api::SolveResult::success(api::SolverStatus::kTimeLimit,
                                           solution, total_iterations);
        }

        // A verified certificate outranks continuing, but never outranks a
        // verified optimal iterate or a numerical error.
        if (certificate_possible) {
          error = download_candidate();
          if (!error.ok()) return failure(error);
          const CertificateOutcome outcome = detectCertificate(
              context, built.problem, original_model.lp(), host,
              options.kkt_tolerance, certificate_buffers, total_iterations);
          if (outcome.available) {
            error = postsolve_candidate();
            if (!error.ok()) return failure(error);
            logCertificateCompletion(context, outcome.status,
                                     outcome.certificate, total_iterations,
                                     elapsedSeconds(start));
            return certifiedResult(outcome.status, solution, total_iterations,
                                   outcome.certificate);
          }
        }
      }

      if (restart_reference_pending) {
        if (context.diagnostics != nullptr) {
          ++context.diagnostics->restart_reference_checks;
          ++context.diagnostics->fixed_point_evaluations;
        }
        const FixedPointErrorResult reference =
            evaluateCandidateFixedPointError(*selected_backend, problem_handle,
                                             state, vectors.aty_delta,
                                             current_steps, context,
                                             total_iterations,
                                             SolvePhase::kReferenceFixedPoint);
        if (!reference.ok()) return failure(reference.error);
        restart_state = resetFixedPointRestart(reference.fixed_point_error);
        restart_reference_pending = false;
      }

      if (cadence_boundary) {
        if (context.diagnostics != nullptr)
          ++context.diagnostics->fixed_point_evaluations;
        const FixedPointErrorResult fixed_point = fixed_point_from_packet
            ? computeFixedPointError(
                  check.primal_delta_norm, check.dual_delta_norm,
                  check.delta_dot, current_steps.primal_weight,
                  current_steps.base_step)
            : evaluateCandidateFixedPointError(
                  *selected_backend, problem_handle, state, vectors.aty_delta,
              current_steps, context, total_iterations,
              SolvePhase::kFixedPointFallback, check.fixed_point_prepared);
        if (!fixed_point.ok()) return failure(fixed_point.error);
        FixedPointRestartInput restart_input;
        restart_input.previous = restart_state;
        restart_input.current_error = fixed_point.fixed_point_error;
        restart_input.halpern_iterations = state.iteration;
        restart_input.total_iterations = total_iterations;
        restart_input.parameters = restart_parameters;
        const FixedPointRestartResult restart =
            evaluateFixedPointRestart(restart_input);
        if (!restart.ok()) return failure(restart.error);
        PrimalWeightUpdateKind pid_outcome = PrimalWeightUpdateKind::kNone;
        if (restart.decision == RestartDecision::kRestart) {
          PrimalWeightUpdateResult primal_weight_update;
          bool install_steps = false;
          if (use_pid) {
            if (!boundary_metrics_valid) {
              return failure(api::Error(
                  api::ErrorCode::kInternal,
                  "HiPDLP rewrite restart boundary has no convergence "
                  "observation"));
            }
            // Both anchor distances are measured before any anchor moves.
            double primal_distance = check.primal_anchor_distance;
            double dual_distance = check.dual_anchor_distance;
            if (!check.anchor_distances_valid) {
              SolvePhaseGuard fallback(context, SolvePhase::kFixedPointFallback,
                                       total_iterations);
              error = selected_backend->differenceNorm2(
                *state.primal_candidate, *state.primal_anchor,
                primal_distance);
              if (!error.ok()) return failure(error);
              error = selected_backend->differenceNorm2(
                *state.dual_candidate, *state.dual_anchor, dual_distance);
              if (!error.ok()) return failure(error);
            }

            PrimalWeightUpdateInput weight_input;
            weight_input.previous = primal_weight_state;
            weight_input.primal_anchor_distance = primal_distance;
            weight_input.dual_anchor_distance = dual_distance;
            weight_input.relative_primal_residual =
                boundary_metrics.relative_primal;
            weight_input.relative_dual_residual =
                boundary_metrics.relative_dual;
            weight_input.primal_step_size = current_steps.primal_step;
            weight_input.dual_step_size = current_steps.dual_step;
            weight_input.parameters = primal_weight_parameters;
            primal_weight_update = updatePrimalWeightAtRestart(weight_input);
            if (!primal_weight_update.ok())
              return failure(primal_weight_update.error);
            pid_outcome = primal_weight_update.kind;
            install_steps = true;
          }

          error =
              restartHalpernPdhgStateFromCandidate(*selected_backend, state);
          if (!error.ok()) return failure(error);
          if (install_steps) {
            // Installed steps only affect the iterations after this boundary.
            primal_weight_state = primal_weight_update.next;
            current_steps = primal_weight_update.step_sizes;
            state.primal_step_size = current_steps.primal_step;
            state.dual_step_size = current_steps.dual_step;
            if (context.diagnostics != nullptr)
              ++context.diagnostics->primal_weight_updates;
          }
          if (context.diagnostics != nullptr) ++context.diagnostics->restarts;
          restart_reference_pending = true;
        } else {
          restart_state = restart.next;
        }
        const HighsInt convergence_boundary =
            total_iterations / kConvergenceCadence;
        if (convergence_boundary % 40 == 0) {
          logBoundary(context, total_iterations, boundary_metrics,
                      current_steps, restart.condition, use_pid, pid_outcome);
        }
      }
      boundary.finish(true);
      if (cadence_boundary) {
        iteration_block.finish(true);
        iteration_block.begin();
      }
    }
  } catch (const std::bad_alloc&) {
    return failure(api::Error(api::ErrorCode::kAllocationFailure,
                              "HiPDLP rewrite solve allocation failed"));
  }
}

}  // namespace algorithm
}  // namespace hipdlp
