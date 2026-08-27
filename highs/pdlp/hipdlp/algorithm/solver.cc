#include "pdlp/hipdlp/algorithm/solver.hpp"

namespace hipdlp {
namespace algorithm {

api::SolveResult solve(const api::SolverOptions& options,
                       const model::OriginalModel& original_model) {
  (void)options;
  (void)original_model;
  return api::SolveResult::failure(
      api::ErrorCode::kNotReady,
      "the rewritten HiPDLP implementation is not ready");
}

}  // namespace algorithm
}  // namespace hipdlp