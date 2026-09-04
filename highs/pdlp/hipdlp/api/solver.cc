#include "pdlp/hipdlp/api/solver.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "io/HighsIO.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/legacy/legacy_solver_adapter.hpp"

namespace hipdlp {
namespace api {

namespace {

class HighsSolveLogSink : public algorithm::SolveLogSink {
 public:
  explicit HighsSolveLogSink(const HighsLogOptions& log_options)
      : log_options_(log_options) {}

  void info(const std::string& message) override {
    highsLogUser(log_options_, HighsLogType::kInfo, "%s\n", message.c_str());
  }

  void warning(const std::string& message) override {
    highsLogUser(log_options_, HighsLogType::kWarning, "%s\n", message.c_str());
  }

 private:
  const HighsLogOptions& log_options_;
};

// The rewrite is given what is left of the global budget, never a second copy.
double remainingTimeLimit(const double time_limit, HighsTimer& timer) {
  if (time_limit == kHighsInf || !std::isfinite(time_limit)) return kHighsInf;
  return std::max(0.0, time_limit - timer.read());
}

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
                  const MigrationOptions& migration_options,
                  SolveResult* rewrite_result) {
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

  SolverOptions rewrite_options = translation.options;
  rewrite_options.time_limit =
      remainingTimeLimit(translation.options.time_limit, timer);

  HighsSolveLogSink log_sink(options.log_options);
  algorithm::SolveContext context;
  context.log = &log_sink;

  const model::OriginalModel original_model(lp);
  const SolveResult result =
      algorithm::solve(rewrite_options, original_model, context);
  if (rewrite_result != nullptr) *rewrite_result = result;
  return applyRewriteResult(options, result, highs_basis, highs_solution,
                            model_status, highs_info);
}

}  // namespace api
}  // namespace hipdlp
