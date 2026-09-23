#ifndef PDLP_HIPDLP_ALGORITHM_CONVERGENCE_HPP_
#define PDLP_HIPDLP_ALGORITHM_CONVERGENCE_HPP_

#include <vector>

#include "lp_data/HConst.h"
#include "pdlp/hipdlp/api/status.hpp"

namespace hipdlp {
namespace model {
class Problem;
}
namespace algorithm {

enum class StoppingDecision {
  kContinue,
  kOptimal,
  kIterationLimit,
  kTimeLimit,
  kNumericalError
};

struct RawConvergenceMetrics {
  double primal_residual;
  double dual_residual;
  double objective_gap;
  double primal_objective;
  double dual_objective;
  double unscaled_rhs_norm;
  double unscaled_cost_norm;

  RawConvergenceMetrics();
};

struct ConvergenceTolerances {
  double primal;
  double dual;
  double gap;

  ConvergenceTolerances();
};

struct StoppingLimits {
  HighsInt iteration;
  double time;

  StoppingLimits();
};

struct StoppingInput {
  RawConvergenceMetrics raw;
  ConvergenceTolerances tolerances;
  HighsInt iteration;
  double elapsed_time;
  StoppingLimits limits;

  StoppingInput();
};

struct NormalizedConvergenceMetrics {
  double relative_primal;
  double relative_dual;
  double relative_gap;

  NormalizedConvergenceMetrics();
};

struct StoppingResult {
  api::Error error;
  StoppingDecision decision;
  NormalizedConvergenceMetrics metrics;
  bool metrics_valid;

  StoppingResult();

  bool ok() const;
};

api::Error computeRawConvergenceMetrics(
  const model::Problem& problem, const std::vector<double>& primal,
  const std::vector<double>& dual, const std::vector<double>& ax,
  const std::vector<double>& aty, RawConvergenceMetrics& raw);

StoppingResult evaluateStopping(const StoppingInput& input);

bool shouldVerifyOriginalKkt(const StoppingResult& stopping, HighsInt iteration,
                             const ConvergenceTolerances& tolerances,
                             double kkt_tolerance);

}  // namespace algorithm
}  // namespace hipdlp

#endif