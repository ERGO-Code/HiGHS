#ifndef PDLP_HIPDLP_API_SOLVER_HPP_
#define PDLP_HIPDLP_API_SOLVER_HPP_

#include "pdlp/hipdlp/api/options.hpp"

namespace hipdlp {
namespace api {

// rewrite_result, when given, receives the structured result of the rewrite
// path so internal callers and benchmarks can see the certificate. It is left
// untouched by the legacy path.
HighsStatus solve(const HighsOptions& options, HighsTimer& timer,
                  const HighsLp& lp, HighsBasis& highs_basis,
                  HighsSolution& highs_solution, HighsModelStatus& model_status,
                  HighsInfo& highs_info, HighsCallback& callback,
                  const MigrationOptions& migration_options,
                  SolveResult* rewrite_result = nullptr);

}  // namespace api
}  // namespace hipdlp

#endif