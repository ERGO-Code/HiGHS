#ifndef PDLP_HIPDLP_LEGACY_LEGACY_SOLVER_ADAPTER_HPP_
#define PDLP_HIPDLP_LEGACY_LEGACY_SOLVER_ADAPTER_HPP_

#include "lp_data/HighsSolution.h"

namespace hipdlp {
namespace legacy {

HighsStatus solve(const HighsOptions& options, HighsTimer& timer,
                  const HighsLp& lp, HighsBasis& highs_basis,
                  HighsSolution& highs_solution, HighsModelStatus& model_status,
                  HighsInfo& highs_info, HighsCallback& callback);

}  // namespace legacy
}  // namespace hipdlp

#endif