#include "pdlp/hipdlp/api/options.hpp"

#include <cmath>

namespace hipdlp {
namespace api {

MigrationOptions::MigrationOptions()
    : implementation(SolverImplementation::kLegacy),
      backend(backend::Kind::kAuto) {}

MigrationOptions::MigrationOptions(const SolverImplementation implementation,
                                   const backend::Kind backend)
    : implementation(implementation), backend(backend) {}

SolverOptions::SolverOptions()
    : implementation(SolverImplementation::kLegacy),
      backend(backend::Kind::kAuto),
      time_limit(0.0),
      kkt_tolerance(0.0),
      features_off(0),
      iteration_limit(0),
      restart_strategy(0),
      step_size_strategy(0),
      optimality_tolerance(0.0) {}

bool OptionTranslationResult::ok() const { return error.ok(); }

MigrationOptions defaultMigrationOptions() { return MigrationOptions(); }

bool isValidImplementation(const SolverImplementation implementation) {
  return implementation == SolverImplementation::kLegacy ||
         implementation == SolverImplementation::kRewrite;
}

OptionTranslationResult translateOptions(const HighsOptions& highs_options,
                                         const MigrationOptions& migration) {
  OptionTranslationResult result;
  if (!isValidImplementation(migration.implementation)) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid internal HiPDLP implementation selector");
    return result;
  }
  if (!backend::isValidKind(migration.backend)) {
    result.error =
        Error(ErrorCode::kInvalidOption, "invalid internal HiPDLP backend");
    return result;
  }
  if (migration.implementation == SolverImplementation::kLegacy &&
      migration.backend != backend::Kind::kAuto) {
    result.error = Error(
        ErrorCode::kInvalidOption,
        "the legacy HiPDLP implementation does not accept a rewrite backend");
    return result;
  }
  if (std::isnan(highs_options.time_limit) || highs_options.time_limit < 0.0) {
    result.error =
        Error(ErrorCode::kInvalidOption, "invalid HiPDLP rewrite time limit");
    return result;
  }
  if (!std::isfinite(highs_options.kkt_tolerance) ||
      highs_options.kkt_tolerance <= 0.0) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid HiPDLP rewrite KKT tolerance");
    return result;
  }
  if (highs_options.pdlp_iteration_limit < 0) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid HiPDLP rewrite iteration limit");
    return result;
  }
  if (highs_options.pdlp_scaling_mode < kPdlpScalingMin ||
      highs_options.pdlp_scaling_mode > kPdlpScalingMax) {
    result.error =
        Error(ErrorCode::kInvalidOption, "invalid HiPDLP rewrite scaling mode");
    return result;
  }
  if (highs_options.pdlp_ruiz_iterations < 0) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid HiPDLP rewrite Ruiz iteration count");
    return result;
  }
  if (highs_options.pdlp_restart_strategy < kPdlpRestartStrategyMin ||
      highs_options.pdlp_restart_strategy > kPdlpRestartStrategyMax) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid HiPDLP rewrite restart strategy");
    return result;
  }
  if (highs_options.pdlp_step_size_strategy < kPdlpStepSizeStrategyMin ||
      highs_options.pdlp_step_size_strategy > kPdlpStepSizeStrategyMax) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid HiPDLP rewrite step-size strategy");
    return result;
  }
  if (!std::isfinite(highs_options.pdlp_optimality_tolerance) ||
      highs_options.pdlp_optimality_tolerance <= 0.0) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid HiPDLP rewrite optimality tolerance");
    return result;
  }

  result.options.implementation = migration.implementation;
  result.options.backend = migration.backend;
  result.options.time_limit = highs_options.time_limit;
  result.options.kkt_tolerance = highs_options.kkt_tolerance;
  result.options.features_off = highs_options.pdlp_features_off;
  result.options.iteration_limit = highs_options.pdlp_iteration_limit;
  result.options.scaling.legacy_scaling_mode = highs_options.pdlp_scaling_mode;
  result.options.scaling.ruiz_iterations = highs_options.pdlp_ruiz_iterations;
  result.options.scaling.features_off = highs_options.pdlp_features_off;
  result.options.restart_strategy = highs_options.pdlp_restart_strategy;
  result.options.step_size_strategy = highs_options.pdlp_step_size_strategy;
  result.options.optimality_tolerance = highs_options.pdlp_optimality_tolerance;
  return result;
}

}  // namespace api
}  // namespace hipdlp
