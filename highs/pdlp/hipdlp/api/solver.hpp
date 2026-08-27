#ifndef PDLP_HIPDLP_API_SOLVER_HPP_
#define PDLP_HIPDLP_API_SOLVER_HPP_

#include "pdlp/hipdlp/api/options.hpp"

namespace hipdlp {
namespace api {

HighsStatus solve(const HighsOptions& options, HighsTimer& timer,
                  const HighsLp& lp, HighsBasis& highs_basis,
                  HighsSolution& highs_solution, HighsModelStatus& model_status,
                  HighsInfo& highs_info, HighsCallback& callback,
                  const MigrationOptions& migration_options);

}  // namespace api
}  // namespace hipdlp

#endif