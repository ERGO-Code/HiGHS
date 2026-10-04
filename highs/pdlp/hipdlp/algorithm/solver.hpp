#ifndef PDLP_HIPDLP_ALGORITHM_SOLVER_HPP_
#define PDLP_HIPDLP_ALGORITHM_SOLVER_HPP_

#include <cstddef>
#include <string>

#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace algorithm {

// Self-audit of the resources the control loop itself asks the backend for.
struct SolveDiagnostics {
  std::size_t resource_creation_calls_before_iterations = 0;
  std::size_t resource_creation_calls_during_iterations = 0;
  std::size_t primal_weight_updates = 0;
  std::size_t restarts = 0;
  std::size_t check_boundaries = 0;
  std::size_t resident_check_calls = 0;
  std::size_t ordinary_check_boundaries = 0;
  std::size_t ordinary_check_d2h_bytes = 0;
  std::size_t ordinary_check_transfers = 0;
  std::size_t ordinary_check_synchronizations = 0;
  std::size_t candidate_downloads = 0;
  std::size_t host_oracle_checks = 0;
  std::size_t certificate_checks = 0;
  std::size_t original_kkt_verifications = 0;
  std::size_t original_kkt_rejections = 0;
  std::size_t kkt_screens = 0;
  std::size_t kkt_screen_rejections = 0;
  std::size_t restart_reference_checks = 0;
  std::size_t fixed_point_evaluations = 0;
  std::size_t command_buffer_replays = 0;
  std::size_t recorded_iterations = 0;

  double total_wall_seconds = 0.0;
  double preprocessing_seconds = 0.0;
  model::PreprocessingDiagnostics preprocessing;
  double backend_setup_seconds = 0.0;
  double operator_norm_seconds = 0.0;
  double iterations_seconds = 0.0;
  double iteration_block_seconds = 0.0;
  double step_seconds = 0.0;
  double boundary_seconds = 0.0;
  double download_seconds = 0.0;
  double host_preparation_seconds = 0.0;
  double host_verification_seconds = 0.0;
  double reference_fixed_point_seconds = 0.0;
  double fixed_point_fallback_seconds = 0.0;
};

// The solve loop reports through this sink instead of owning HiGHS log options.
class SolveLogSink {
 public:
  virtual ~SolveLogSink() = default;

  virtual void info(const std::string& message) = 0;
  virtual void warning(const std::string& message) = 0;
};

enum class SolvePhase {
  kSolve,
  kPreprocessing,
  kBackendSetup,
  kOperatorNorm,
  kIterations,
  kIterationBlock,
  kStep,
  kBoundary,
  kDownload,
  kHostPreparation,
  kHostVerification,
  kReferenceFixedPoint,
  kFixedPointFallback,
  kCount
};

struct SolveObservation {
  SolvePhase phase;
  bool begin;
  HighsInt iteration;
  bool continued;
};

class SolveObserver {
 public:
  virtual ~SolveObserver() = default;
  virtual void onEvent(const SolveObservation& observation) noexcept = 0;
};

// Internal experiment hook. Observes original-coordinate candidates at ordinary
// boundaries through iteration 200. Returning true requests an iteration-limit
// result; returning false leaves all iterate/restart state in place. Terminal
// convergence, certificates and limits retain precedence. No public option
// enables this hook; the default path performs no extra downloads or checks.
class PilotObserver {
 public:
  virtual ~PilotObserver() = default;
  virtual bool checkpoint(HighsInt iteration, const HighsSolution& solution,
                          double elapsed_seconds, double fixed_point_error,
                          std::size_t restarts) = 0;
};

struct SolveContext {
  SolveLogSink* log = nullptr;
  SolveDiagnostics* diagnostics = nullptr;
  SolveObserver* observer = nullptr;
  PilotObserver* pilot_observer = nullptr;
  bool use_command_buffers = true;
};

// options.time_limit is the remaining local budget for this call alone.
api::SolveResult solve(const api::SolverOptions& options,
                       const model::OriginalModel& original_model,
                       const SolveContext& context = SolveContext());

}  // namespace algorithm
}  // namespace hipdlp

#endif
