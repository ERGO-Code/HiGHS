#ifndef PDLP_HIPDLP_ALGORITHM_SOLVER_HPP_
#define PDLP_HIPDLP_ALGORITHM_SOLVER_HPP_

#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/model/problem.hpp"

namespace hipdlp {
namespace algorithm {

api::SolveResult solve(const api::SolverOptions& options,
                       const model::OriginalModel& original_model);

}  // namespace algorithm
}  // namespace hipdlp

#endif