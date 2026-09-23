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
  std::size_t resource_creation_calls_before_iterations;
  std::size_t resource_creation_calls_during_iterations;
  std::size_t primal_weight_updates;
  std::size_t restarts;
  std::size_t check_boundaries;
  std::size_t resident_check_calls;
  std::size_t ordinary_check_boundaries;
  std::size_t ordinary_check_d2h_bytes;
  std::size_t ordinary_check_transfers;
  std::size_t ordinary_check_synchronizations;
  std::size_t candidate_downloads;
  std::size_t host_oracle_checks;
  std::size_t certificate_checks;
  std::size_t restart_reference_checks;
  std::size_t fixed_point_evaluations;
  std::size_t command_buffer_replays;
  std::size_t recorded_iterations;

  double total_wall_seconds;
  double preprocessing_seconds;
  double backend_setup_seconds;
  double operator_norm_seconds;
  double iterations_seconds;
  double iteration_block_seconds;
  double step_seconds;
  double boundary_seconds;
  double download_seconds;
  double host_preparation_seconds;
  double host_verification_seconds;
  double reference_fixed_point_seconds;
  double fixed_point_fallback_seconds;

  SolveDiagnostics();
};

// The solve loop reports through this sink instead of owning HiGHS log options.
class SolveLogSink {
 public:
  virtual ~SolveLogSink();

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
  virtual ~SolveObserver() {}
  virtual void onEvent(const SolveObservation& observation) noexcept = 0;
};

struct SolveContext {
  SolveLogSink* log;
  SolveDiagnostics* diagnostics;
  SolveObserver* observer;
  bool use_command_buffers;

  SolveContext();
};

// options.time_limit is the remaining local budget for this call alone.
api::SolveResult solve(const api::SolverOptions& options,
                       const model::OriginalModel& original_model,
                       const SolveContext& context = SolveContext());

}  // namespace algorithm
}  // namespace hipdlp

#endif