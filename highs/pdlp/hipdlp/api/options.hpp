#ifndef PDLP_HIPDLP_API_OPTIONS_HPP_
#define PDLP_HIPDLP_API_OPTIONS_HPP_

#include "lp_data/HighsOptions.h"
#include "pdlp/hipdlp/api/status.hpp"
#include "pdlp/hipdlp/backend/types.hpp"
#include "pdlp/hipdlp/model/scaling.hpp"

namespace hipdlp {
namespace api {

enum class SolverImplementation { kLegacy, kRewrite };

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
  HighsInt restart_strategy;
  HighsInt step_size_strategy;
  double optimality_tolerance;

  SolverOptions();
};

struct OptionTranslationResult {
  SolverOptions options;
  Error error;

  bool ok() const;
};

MigrationOptions defaultMigrationOptions();
bool isValidImplementation(SolverImplementation implementation);
OptionTranslationResult translateOptions(const HighsOptions& highs_options,
                                         const MigrationOptions& migration);

}  // namespace api
}  // namespace hipdlp

#endif
