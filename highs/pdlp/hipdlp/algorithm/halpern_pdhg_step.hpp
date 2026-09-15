#ifndef PDLP_HIPDLP_ALGORITHM_HALPERN_PDHG_STEP_HPP_
#define PDLP_HIPDLP_ALGORITHM_HALPERN_PDHG_STEP_HPP_

#include "pdlp/hipdlp/backend/backend.hpp"

namespace hipdlp {
namespace algorithm {

struct HalpernPdhgState {
  backend::ProblemHandle* problem;
  backend::WorkspaceHandle* workspace;
  backend::VectorHandle* current_primal;
  backend::VectorHandle* current_dual;
  backend::VectorHandle* next_primal;
  backend::VectorHandle* next_dual;
  backend::VectorHandle* primal_anchor;
  backend::VectorHandle* dual_anchor;
  backend::VectorHandle* primal_candidate;
  backend::VectorHandle* dual_candidate;
  backend::VectorHandle* reflected_primal;
  backend::VectorHandle* reflected_dual;
  backend::VectorHandle* current_aty;
  backend::VectorHandle* reflected_ax;
  double primal_step_size;
  double dual_step_size;
  HighsInt iteration;

  HalpernPdhgState();
};

api::Error enqueueHalpernPdhgStep(backend::Backend& backend,
                                  HalpernPdhgState& state);
api::Error restartHalpernPdhgStateFromCandidate(backend::Backend& backend,
                                                HalpernPdhgState& state);

}  // namespace algorithm
}  // namespace hipdlp

#endif