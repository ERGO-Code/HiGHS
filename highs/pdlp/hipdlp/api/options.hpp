#ifndef PDLP_HIPDLP_API_OPTIONS_HPP_
#define PDLP_HIPDLP_API_OPTIONS_HPP_

#include <string>

#include "lp_data/HighsOptions.h"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/types.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

namespace hipdlp {
namespace api {

enum class SolverImplementation { kLegacy, kRewrite };

// The only step-size and restart behaviours the rewrite actually implements.
enum class StepSizePolicy { kFixed, kPid };
enum class RestartPolicy { kHalpern };

struct EffectivePolicy {
  StepSizePolicy step_size;
  RestartPolicy restart;
  HighsInt requested_step_size_strategy;
  HighsInt requested_restart_strategy;
  bool step_size_substituted;
  bool restart_substituted;

  EffectivePolicy();
};

struct MigrationOptions {
  SolverImplementation implementation;
  backend::Kind backend;

  MigrationOptions();
  MigrationOptions(SolverImplementation implementation, backend::Kind backend);
};

struct SolverOptions {
  SolverImplementation implementation;
  backend::Kind backend;
  double time_limit;
  double kkt_tolerance;
  HighsInt features_off;
  HighsInt iteration_limit;
  model::ScalingOptions scaling;
  EffectivePolicy policy;
  double optimality_tolerance;
  bool early_kkt_check;
  bool scaling_selector;

  SolverOptions();
};

struct OptionTranslationResult {
  SolverOptions options;
  Error error;

  bool ok() const;
};

MigrationOptions defaultMigrationOptions();
bool isValidImplementation(SolverImplementation implementation);
bool backendKindFromString(const std::string& value, backend::Kind& kind);
bool scalingModeFromString(const std::string& value, model::ScalingMode& mode);
OptionTranslationResult translateOptions(const HighsOptions& highs_options,
                                         const MigrationOptions& migration);

}  // namespace api
}  // namespace hipdlp

#endif
