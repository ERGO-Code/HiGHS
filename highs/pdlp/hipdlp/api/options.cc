#include "pdlp/hipdlp/api/options.hpp"

#include <cmath>

namespace hipdlp {
namespace api {

EffectivePolicy::EffectivePolicy()
    : step_size(StepSizePolicy::kFixed),
      restart(RestartPolicy::kHalpern),
      requested_step_size_strategy(kPdlpStepSizeStrategyFixed),
      requested_restart_strategy(kPdlpRestartStrategyHalpern),
      step_size_substituted(false),
      restart_substituted(false) {}

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
      optimality_tolerance(0.0),
      early_kkt_check(false),
      scaling_selector(false) {}

bool OptionTranslationResult::ok() const { return error.ok(); }

MigrationOptions defaultMigrationOptions() {
  return MigrationOptions(SolverImplementation::kRewrite, backend::Kind::kAuto);
}

bool isValidImplementation(const SolverImplementation implementation) {
  return implementation == SolverImplementation::kLegacy ||
         implementation == SolverImplementation::kRewrite;
}

bool backendKindFromString(const std::string& value, backend::Kind& kind) {
  if (value == kHiPdlpAutoString) {
    kind = backend::Kind::kAuto;
  } else if (value == kHiPdlpBackendCpuString) {
    kind = backend::Kind::kCpu;
  } else if (value == kHiPdlpBackendCudaString) {
    kind = backend::Kind::kCuda;
  } else if (value == kHiPdlpBackendOpenmpString) {
    kind = backend::Kind::kOpenmp;
  } else if (value == kHiPdlpBackendRocmString) {
    kind = backend::Kind::kRocm;
  } else {
    return false;
  }
  return true;
}

bool scalingModeFromString(const std::string& value, model::ScalingMode& mode) {
  if (value == kHiPdlpAutoString) {
    mode = model::ScalingMode::kAuto;
  } else if (value == kHighsOffString) {
    mode = model::ScalingMode::kOff;
  } else if (value == kHiPdlpScalingCurtisReidString) {
    mode = model::ScalingMode::kCurtisReid;
  } else if (value == kHiPdlpScalingLegacyString) {
    mode = model::ScalingMode::kLegacy;
  } else {
    return false;
  }
  return true;
}

namespace {

Error resolveStepSizePolicy(const HighsInt requested, EffectivePolicy& policy) {
  policy.requested_step_size_strategy = requested;
  policy.step_size_substituted = false;
  if (requested == kPdlpStepSizeStrategyFixed) {
    policy.step_size = StepSizePolicy::kFixed;
    return Error();
  }
  if (requested == kPdlpStepSizeStrategyPid) {
    policy.step_size = StepSizePolicy::kPid;
    return Error();
  }
  if (requested == kPdlpStepSizeStrategyAdaptive) {
    // Never silently promoted to PID the way the legacy solver does.
    policy.step_size = StepSizePolicy::kFixed;
    policy.step_size_substituted = true;
    return Error();
  }
  return Error(ErrorCode::kInvalidOption,
               "the HiPDLP rewrite does not implement the Malitsky-Pock "
               "step-size strategy; use pdlp_step_size_strategy in {0: fixed, "
               "3: PID}");
}

void resolveRestartPolicy(const HighsInt requested, const HighsInt features_off,
                          EffectivePolicy& policy) {
  policy.requested_restart_strategy = requested;
  policy.restart = RestartPolicy::kHalpern;
  policy.restart_substituted = requested != kPdlpRestartStrategyHalpern ||
                               (features_off & kPdlpRestartOff) != 0;
}

}  // namespace

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
  backend::Kind public_backend = backend::Kind::kAuto;
  if (!backendKindFromString(highs_options.hipdlp_backend, public_backend)) {
    result.error =
        Error(ErrorCode::kInvalidOption, "invalid HiPDLP backend option");
    return result;
  }
  model::ScalingMode public_scaling = model::ScalingMode::kAuto;
  if (!scalingModeFromString(highs_options.hipdlp_scaling, public_scaling)) {
    result.error =
        Error(ErrorCode::kInvalidOption, "invalid HiPDLP scaling option");
    return result;
  }
  if (highs_options.hipdlp_curtis_reid_iterations <= 0) {
    result.error = Error(ErrorCode::kInvalidOption,
                         "invalid HiPDLP Curtis-Reid iteration count");
    return result;
  }

  result.options.implementation = migration.implementation;
  // A non-auto migration backend is the test and migration override; otherwise
  // the public option decides.
  result.options.backend = migration.backend != backend::Kind::kAuto
                               ? migration.backend
                               : public_backend;
  result.options.time_limit = highs_options.time_limit;
  result.options.kkt_tolerance = highs_options.kkt_tolerance;
  result.options.features_off = highs_options.pdlp_features_off;
  result.options.iteration_limit = highs_options.pdlp_iteration_limit;
  result.options.scaling.mode = public_scaling;
  result.options.scaling.curtis_reid_sweeps =
      highs_options.hipdlp_curtis_reid_iterations;
  result.options.scaling.legacy_scaling_mode = highs_options.pdlp_scaling_mode;
  result.options.scaling.ruiz_iterations = highs_options.pdlp_ruiz_iterations;
  result.options.scaling.features_off = highs_options.pdlp_features_off;
  result.options.optimality_tolerance = highs_options.pdlp_optimality_tolerance;
  result.options.early_kkt_check = highs_options.hipdlp_early_kkt_check;
  result.options.scaling_selector = highs_options.hipdlp_scaling_selector;

  // Only the rewrite is constrained by the policies it actually implements.
  if (migration.implementation == SolverImplementation::kRewrite) {
    const Error step_size_error = resolveStepSizePolicy(
        highs_options.pdlp_step_size_strategy, result.options.policy);
    if (!step_size_error.ok()) {
      result.error = step_size_error;
      result.options = SolverOptions();
      return result;
    }
    resolveRestartPolicy(highs_options.pdlp_restart_strategy,
                         highs_options.pdlp_features_off,
                         result.options.policy);
  } else {
    result.options.policy.requested_step_size_strategy =
        highs_options.pdlp_step_size_strategy;
    result.options.policy.requested_restart_strategy =
        highs_options.pdlp_restart_strategy;
  }
  return result;
}

}  // namespace api
}  // namespace hipdlp
