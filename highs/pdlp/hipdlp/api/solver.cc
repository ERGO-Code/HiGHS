#include "pdlp/hipdlp/api/solver.hpp"

#include "io/HighsIO.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/legacy/legacy_solver_adapter.hpp"

namespace hipdlp {
namespace api {

namespace {

HighsStatus applyRewriteResult(const HighsOptions& options,
                               const SolveResult& result,
                               HighsBasis& highs_basis,
                               HighsSolution& highs_solution,
                               HighsModelStatus& model_status,
                               HighsInfo& highs_info) {
  if (!result.error.ok()) {
    highsLogUser(options.log_options, HighsLogType::kError, "%s\n",
                 result.error.message.c_str());
  }
  return applySolveResult(result, highs_basis, highs_solution, model_status,
                          highs_info);
}

}  // namespace

HighsStatus solve(const HighsOptions& options, HighsTimer& timer,
                  const HighsLp& lp, HighsBasis& highs_basis,
                  HighsSolution& highs_solution, HighsModelStatus& model_status,
                  HighsInfo& highs_info, HighsCallback& callback,
                  const MigrationOptions& migration_options) {
  const OptionTranslationResult translation =
      translateOptions(options, migration_options);
  if (!translation.ok()) {
    return applyRewriteResult(
        options,
        SolveResult::failure(translation.error.code, translation.error.message),
        highs_basis, highs_solution, model_status, highs_info);
  }

  if (translation.options.implementation == SolverImplementation::kLegacy) {
    return legacy::solve(options, timer, lp, highs_basis, highs_solution,
                         model_status, highs_info, callback);
  }

  const model::OriginalModel original_model(lp);
  return applyRewriteResult(
      options, algorithm::solve(translation.options, original_model),
      highs_basis, highs_solution, model_status, highs_info);
}

}  // namespace api
}  // namespace hipdlp