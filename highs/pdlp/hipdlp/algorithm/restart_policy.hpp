#ifndef PDLP_HIPDLP_ALGORITHM_RESTART_POLICY_HPP_
#define PDLP_HIPDLP_ALGORITHM_RESTART_POLICY_HPP_

#include "pdlp/hipdlp/algorithm/averaging_policy.hpp"
#include "pdlp/hipdlp/algorithm/convergence.hpp"
#include "pdlp/hipdlp/api/status.hpp"

namespace hipdlp {
namespace algorithm {

enum class RestartDecision { kNone, kContinue, kRestart };

enum class RestartCandidateKind { kNone, kCurrent, kAveraged };

struct RestartParameters {
  double primal_weight;
  double sufficient_decay;
  double necessary_decay;
  double artificial_restart_threshold;

  RestartParameters();
};

struct RestartState {
  NormalizedConvergenceMetrics restart_metrics;
  NormalizedConvergenceMetrics last_candidate_metrics;
  HighsInt restart_count;
  bool initialized;
  bool history_valid;

  RestartState();
};

struct RestartInput {
  RestartState previous;
  NormalizedConvergenceMetrics current;
  NormalizedConvergenceMetrics averaged;
  AveragingState averaging;
  HighsInt iteration;
  HighsInt iterations_since_restart;
  RestartParameters parameters;

  RestartInput();
};

struct RestartCandidate {
  RestartCandidateKind kind;
  NormalizedConvergenceMetrics metrics;
  double score;
  bool averaged_available;

  RestartCandidate();
};

struct RestartResult {
  api::Error error;
  RestartState next;
  RestartDecision decision;
  RestartCandidate candidate;

  RestartResult();

  bool ok() const;
};

RestartResult evaluateRestart(const RestartInput& input);
RestartState resetRestart();

}  // namespace algorithm
}  // namespace hipdlp

#endif