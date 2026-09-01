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

  SolveDiagnostics();
};

// The solve loop reports through this sink instead of owning HiGHS log options.
class SolveLogSink {
 public:
  virtual ~SolveLogSink();

  virtual void info(const std::string& message) = 0;
  virtual void warning(const std::string& message) = 0;
};

struct SolveContext {
  SolveLogSink* log;
  SolveDiagnostics* diagnostics;

  SolveContext();
};

// options.time_limit is the remaining local budget for this call alone.
api::SolveResult solve(const api::SolverOptions& options,
                       const model::OriginalModel& original_model,
                       const SolveContext& context = SolveContext());

}  // namespace algorithm
}  // namespace hipdlp

#endif