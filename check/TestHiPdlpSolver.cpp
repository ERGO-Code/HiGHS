#include <cstdio>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/api/options.hpp"
#include "pdlp/hipdlp/api/solver.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

namespace {

const double kSolverTolerance = 1e-4;
const double kCudaTolerance = 1e-9;

bool cudaBackendUsable() {
  return hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda).ok();
}

HighsLp makeTinyFeasibleLp() {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {1.0};
  lp.col_lower_ = {0.0};
  lp.col_upper_ = {kHighsInf};
  lp.row_lower_ = {1.0};
  lp.row_upper_ = {1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1.0};
  lp.sense_ = ObjSense::kMinimize;
  return lp;
}

HighsLp makeScaledFeasibleLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_cost_[0] = 3.0;
  lp.row_lower_[0] = 8.0;
  lp.row_upper_[0] = 8.0;
  lp.a_matrix_.value_[0] = 4.0;
  return lp;
}

HighsLp makeMaximizationLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_cost_[0] = 3.0;
  lp.row_lower_[0] = -kHighsInf;
  lp.row_upper_[0] = 4.0;
  lp.a_matrix_.value_[0] = 2.0;
  lp.sense_ = ObjSense::kMaximize;
  lp.offset_ = 5.0;
  return lp;
}

HighsLp makeAlreadyOptimalLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_cost_ = {0.0};
  lp.row_lower_ = {0.0};
  lp.row_upper_ = {kHighsInf};
  return lp;
}

// Every column is boxed, so the relative dual residual is identically zero and
// the PID primal-weight ratio is always undefined.
HighsLp makeBoxedColumnLp() {
  HighsLp lp = makeTinyFeasibleLp();
  lp.col_upper_ = {5.0};
  return lp;
}

class RecordingLogSink : public hipdlp::algorithm::SolveLogSink {
 public:
  void info(const std::string& message) override {
    messages.push_back(message);
  }
  void warning(const std::string& message) override {
    messages.push_back(message);
    warnings.push_back(message);
  }

  std::size_t countContaining(const std::string& needle) const {
    std::size_t count = 0;
    for (std::size_t index = 0; index < messages.size(); ++index) {
      if (messages[index].find(needle) != std::string::npos) ++count;
    }
    return count;
  }

  std::vector<std::string> messages;
  std::vector<std::string> warnings;
};

hipdlp::api::SolverOptions makeCpuFixedOptions() {
  hipdlp::api::SolverOptions options;
  options.implementation = hipdlp::api::SolverImplementation::kRewrite;
  options.backend = hipdlp::backend::Kind::kCpu;
  options.time_limit = kHighsInf;
  options.kkt_tolerance = kSolverTolerance;
  options.features_off = kPdlpAllFeaturesOn;
  options.iteration_limit = 10000;
  options.scaling.mode = hipdlp::model::ScalingMode::kOff;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kFixed;
  options.policy.restart = hipdlp::api::RestartPolicy::kHalpern;
  options.optimality_tolerance = kSolverTolerance;
  return options;
}

}  // namespace

TEST_CASE("hi-pdlp-rewrite-solver-solves-tiny-feasible-lp",
          "[pdlp][hipdlp][phase7][solver][optimal]") {
  const HighsLp lp = makeTinyFeasibleLp();
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      makeCpuFixedOptions(), hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  REQUIRE(result.solution.col_value.size() == 1);
  REQUIRE(result.solution.row_value.size() == 1);
  REQUIRE(result.solution.col_dual.size() == 1);
  REQUIRE(result.solution.row_dual.size() == 1);
  CHECK(result.solution.col_value[0] == Approx(1.0).margin(kSolverTolerance));
  CHECK(result.solution.row_value[0] == Approx(1.0).margin(kSolverTolerance));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-returns-final-candidate-at-iteration-limit",
          "[pdlp][hipdlp][phase7][solver][iteration-limit]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 3;
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kIterationLimit);
  CHECK(result.iteration_count == 3);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  REQUIRE(result.solution.col_value.size() == 1);
  REQUIRE(result.solution.row_dual.size() == 1);
  CHECK(result.solution.col_value[0] == Approx(0.994008).margin(1e-12));
  CHECK(result.solution.row_value[0] == Approx(0.994008).margin(1e-12));
  CHECK(result.solution.row_dual[0] == Approx(1.009960032).margin(1e-12));
  CHECK(result.solution.col_dual[0] == Approx(-0.009960032).margin(1e-12));
}

TEST_CASE("hi-pdlp-rewrite-solver-honors-zero-time-limit-without-iteration",
          "[pdlp][hipdlp][phase7][solver][time-limit]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.time_limit = 0.0;
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kTimeLimit);
  CHECK(result.iteration_count == 0);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  REQUIRE(result.solution.col_value.size() == 1);
  REQUIRE(result.solution.row_dual.size() == 1);
  CHECK(result.solution.col_value[0] == 0.0);
  CHECK(result.solution.row_value[0] == 0.0);
  CHECK(result.solution.row_dual[0] == 0.0);
  CHECK(result.solution.col_dual[0] == 1.0);

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  CHECK_FALSE(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-postsolves-scaled-solution",
          "[pdlp][hipdlp][phase7][solver][scaling][postsolve]") {
  const HighsLp lp = makeScaledFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.scaling.mode = hipdlp::model::ScalingMode::kCurtisReid;
  const hipdlp::api::SolveResult result =
      hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  CHECK(result.solution.col_value[0] == Approx(2.0).margin(kSolverTolerance));
  CHECK(result.solution.row_value[0] ==
        Approx(8.0).margin(4.0 * kSolverTolerance));
  CHECK(result.solution.row_dual[0] == Approx(0.75).margin(kSolverTolerance));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-scaling-stages-solve-paired-inequalities",
          "[pdlp][hipdlp][phase7][solver][scaling][sctest]") {
  Highs highs;
  REQUIRE(highs.readModel(std::string(HIGHS_DIR) +
                          "/check/instances/sctest.mps") ==
          HighsStatus::kOk);
  const HighsLp lp = highs.getLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 100000;
  options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  options.scaling.mode = hipdlp::model::ScalingMode::kLegacy;

  SECTION("Ruiz") {
    options.scaling.legacy_scaling_mode = kPdlpScalingRuiz;
  }
  SECTION("Pock-Chambolle") {
    options.scaling.legacy_scaling_mode = kPdlpScalingPC;
  }

  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp));
  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.iteration_count < options.iteration_limit);
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution,
                                      options.kkt_tolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-restores-maximization-sense-and-offset",
          "[pdlp][hipdlp][phase7][solver][postsolve][maximization]") {
  const HighsLp lp = makeMaximizationLp();
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      makeCpuFixedOptions(), hipdlp::model::OriginalModel(lp));

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  CHECK(result.solution.col_value[0] == Approx(2.0).margin(kSolverTolerance));
  CHECK(result.solution.row_value[0] ==
        Approx(4.0).margin(2.0 * kSolverTolerance));
  CHECK(result.solution.row_dual[0] == Approx(1.5).margin(kSolverTolerance));
  const double objective =
      lp.offset_ + lp.col_cost_[0] * result.solution.col_value[0];
  CHECK(objective == Approx(11.0).margin(3.0 * kSolverTolerance));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-requires-independent-kkt-acceptance",
          "[pdlp][hipdlp][phase7][solver][kkt-gate]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions loose = makeCpuFixedOptions();
  loose.iteration_limit = 41;
  loose.optimality_tolerance = 0.4;
  loose.kkt_tolerance = 1.0;
  const hipdlp::api::SolveResult accepted =
      hipdlp::algorithm::solve(loose, hipdlp::model::OriginalModel(lp));
  INFO(accepted.error.message);
  REQUIRE(accepted.ok());
  REQUIRE(accepted.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(accepted.iteration_count == 40);

  hipdlp::api::SolverOptions strict = loose;
  strict.kkt_tolerance = 1e-12;
  const hipdlp::api::SolveResult rejected =
      hipdlp::algorithm::solve(strict, hipdlp::model::OriginalModel(lp));
  INFO(rejected.error.message);
  REQUIRE(rejected.ok());
  CHECK(rejected.status == hipdlp::api::SolverStatus::kIterationLimit);
  CHECK(rejected.iteration_count == 41);
  REQUIRE(rejected.has_solution);
  REQUIRE(rejected.solution.value_valid);
  REQUIRE(rejected.solution.dual_valid);
  CHECK(rejected.solution.col_value != accepted.solution.col_value);
  CHECK(rejected.solution.row_dual != accepted.solution.row_dual);

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, rejected.solution,
                                      strict.kkt_tolerance);
  REQUIRE(verification.input_valid);
  REQUIRE(verification.metrics_finite);
  CHECK_FALSE(verification.pass);
}

TEST_CASE("hi-pdlp-rewrite-solver-rejects-unsupported-execution-scope",
          "[pdlp][hipdlp][phase7][solver][options][errors]") {
  HighsLp malformed = makeTinyFeasibleLp();
  malformed.a_matrix_.start_.clear();

  SECTION("an unavailable backend is rejected before model processing") {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.backend = hipdlp::backend::Kind::kOpenmp;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
  }

  SECTION("an invalid backend value is rejected before model processing") {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.backend = static_cast<hipdlp::backend::Kind>(-1);
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "HiPDLP rewrite solve received an invalid backend");
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
  }

  SECTION("an invalid step-size policy is rejected before model processing") {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.policy.step_size = static_cast<hipdlp::api::StepSizePolicy>(7);
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(result.error.message ==
          "HiPDLP rewrite solve received an invalid step-size policy");
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
  }
}

TEST_CASE("hi-pdlp-rewrite-solver-returns-structured-model-errors",
          "[pdlp][hipdlp][phase7][solver][model][errors]") {
  SECTION("malformed model preserves the model error") {
    HighsLp malformed = makeTinyFeasibleLp();
    malformed.col_cost_.clear();
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        makeCpuFixedOptions(), hipdlp::model::OriginalModel(malformed));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidModel);
    CHECK_FALSE(result.error.message.empty());
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
    CHECK_FALSE(result.solution.value_valid);
    CHECK_FALSE(result.solution.dual_valid);
  }

  SECTION("zero operator reports the deferred trivial-LP path") {
    HighsLp zero_operator = makeTinyFeasibleLp();
    zero_operator.a_matrix_.start_ = {0, 0};
    zero_operator.a_matrix_.index_.clear();
    zero_operator.a_matrix_.value_.clear();
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        makeCpuFixedOptions(), hipdlp::model::OriginalModel(zero_operator));
    CHECK_FALSE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kError);
    CHECK(result.error.code == hipdlp::api::ErrorCode::kInvalidModel);
    CHECK(result.error.message ==
          "HiPDLP rewrite does not yet support zero-operator LPs");
    CHECK(result.iteration_count == 0);
    CHECK_FALSE(result.has_solution);
    CHECK_FALSE(result.solution.value_valid);
    CHECK_FALSE(result.solution.dual_valid);
  }
}

TEST_CASE("hi-pdlp-public-options-register-defaults-and-round-trip",
          "[pdlp][hipdlp][phase7][options][public]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);

  std::string backend;
  REQUIRE(highs.getOptionValue(kHiPdlpBackendString, backend) ==
          HighsStatus::kOk);
  CHECK(backend == kHiPdlpAutoString);
  std::string scaling;
  REQUIRE(highs.getOptionValue(kHiPdlpScalingString, scaling) ==
          HighsStatus::kOk);
  CHECK(scaling == kHiPdlpAutoString);
  HighsInt sweeps = 0;
  REQUIRE(highs.getOptionValue(kHiPdlpCurtisReidIterationsString, sweeps) ==
          HighsStatus::kOk);
  CHECK(sweeps == kHiPdlpDefaultCurtisReidIterations);

  SECTION("every documented backend string is accepted") {
    const std::string values[] = {
        kHiPdlpAutoString, kHiPdlpBackendCpuString, kHiPdlpBackendCudaString,
        kHiPdlpBackendOpenmpString, kHiPdlpBackendRocmString};
    const hipdlp::backend::Kind kinds[] = {
        hipdlp::backend::Kind::kAuto, hipdlp::backend::Kind::kCpu,
        hipdlp::backend::Kind::kCuda, hipdlp::backend::Kind::kOpenmp,
        hipdlp::backend::Kind::kRocm};
    for (std::size_t index = 0; index < 5; ++index) {
      REQUIRE(highs.setOptionValue(kHiPdlpBackendString, values[index]) ==
              HighsStatus::kOk);
      std::string read_back;
      REQUIRE(highs.getOptionValue(kHiPdlpBackendString, read_back) ==
              HighsStatus::kOk);
      CHECK(read_back == values[index]);
      hipdlp::backend::Kind kind = hipdlp::backend::Kind::kAuto;
      REQUIRE(hipdlp::api::backendKindFromString(read_back, kind));
      CHECK(kind == kinds[index]);
    }
  }

  SECTION("every documented scaling string is accepted") {
    const std::string values[] = {kHiPdlpAutoString, kHighsOffString,
                                  kHiPdlpScalingCurtisReidString,
                                  kHiPdlpScalingLegacyString};
    const hipdlp::model::ScalingMode modes[] = {
        hipdlp::model::ScalingMode::kAuto, hipdlp::model::ScalingMode::kOff,
        hipdlp::model::ScalingMode::kCurtisReid,
        hipdlp::model::ScalingMode::kLegacy};
    for (std::size_t index = 0; index < 4; ++index) {
      REQUIRE(highs.setOptionValue(kHiPdlpScalingString, values[index]) ==
              HighsStatus::kOk);
      std::string read_back;
      REQUIRE(highs.getOptionValue(kHiPdlpScalingString, read_back) ==
              HighsStatus::kOk);
      CHECK(read_back == values[index]);
      hipdlp::model::ScalingMode mode = hipdlp::model::ScalingMode::kAuto;
      REQUIRE(hipdlp::api::scalingModeFromString(read_back, mode));
      CHECK(mode == modes[index]);
    }
  }

  SECTION("illegal option values are rejected and leave the option alone") {
    CHECK(highs.setOptionValue(kHiPdlpBackendString, "gpu") ==
          HighsStatus::kError);
    CHECK(highs.setOptionValue(kHiPdlpScalingString, "ruiz") ==
          HighsStatus::kError);
    CHECK(highs.setOptionValue(kHiPdlpCurtisReidIterationsString,
                               HighsInt{0}) == HighsStatus::kError);
    CHECK(highs.setOptionValue(kHiPdlpCurtisReidIterationsString,
                               HighsInt{-1}) == HighsStatus::kError);
    CHECK(highs.getOptions().hipdlp_backend == kHiPdlpAutoString);
    CHECK(highs.getOptions().hipdlp_scaling == kHiPdlpAutoString);
    CHECK(highs.getOptions().hipdlp_curtis_reid_iterations ==
          kHiPdlpDefaultCurtisReidIterations);
    hipdlp::backend::Kind kind = hipdlp::backend::Kind::kCpu;
    CHECK_FALSE(hipdlp::api::backendKindFromString("gpu", kind));
    hipdlp::model::ScalingMode mode = hipdlp::model::ScalingMode::kOff;
    CHECK_FALSE(hipdlp::api::scalingModeFromString("ruiz", mode));
  }

  SECTION("a positive Curtis-Reid count reaches the rewrite scaling options") {
    REQUIRE(highs.setOptionValue(kHiPdlpCurtisReidIterationsString,
                                 HighsInt{7}) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue(kHiPdlpScalingString,
                                 kHiPdlpScalingCurtisReidString) ==
            HighsStatus::kOk);
    const hipdlp::api::MigrationOptions rewrite(
        hipdlp::api::SolverImplementation::kRewrite,
        hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    INFO(translated.error.message);
    REQUIRE(translated.ok());
    CHECK(translated.options.scaling.mode ==
          hipdlp::model::ScalingMode::kCurtisReid);
    CHECK(translated.options.scaling.curtis_reid_sweeps == 7);
  }
}

TEST_CASE("hi-pdlp-effective-policy-translation",
          "[pdlp][hipdlp][phase7][options][migration][policy]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite,
      hipdlp::backend::Kind::kAuto);

  SECTION("the public default migration selects rewrite") {
    const hipdlp::api::MigrationOptions defaults =
        hipdlp::api::defaultMigrationOptions();
    CHECK(defaults.implementation ==
          hipdlp::api::SolverImplementation::kRewrite);
    CHECK(defaults.backend == hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), defaults);
    INFO(translated.error.message);
    REQUIRE(translated.ok());
    CHECK(translated.options.implementation ==
          hipdlp::api::SolverImplementation::kRewrite);
  }

  SECTION("default HiGHS options translate for the rewrite selector") {
    REQUIRE(highs.getOptions().pdlp_step_size_strategy ==
            kPdlpStepSizeStrategyAdaptive);
    REQUIRE(highs.getOptions().pdlp_restart_strategy ==
            kPdlpRestartStrategyAdaptive);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    INFO(translated.error.message);
    REQUIRE(translated.ok());
    CHECK(translated.options.implementation ==
          hipdlp::api::SolverImplementation::kRewrite);
    CHECK(translated.options.backend == hipdlp::backend::Kind::kAuto);
    CHECK(translated.options.scaling.mode == hipdlp::model::ScalingMode::kAuto);
    CHECK(translated.options.scaling.curtis_reid_sweeps ==
          kHiPdlpDefaultCurtisReidIterations);
    CHECK(translated.options.policy.step_size ==
          hipdlp::api::StepSizePolicy::kFixed);
    CHECK_FALSE(translated.options.policy.step_size ==
                hipdlp::api::StepSizePolicy::kPid);
    CHECK(translated.options.policy.step_size_substituted);
    CHECK(translated.options.policy.requested_step_size_strategy ==
          kPdlpStepSizeStrategyAdaptive);
    CHECK(translated.options.policy.restart ==
          hipdlp::api::RestartPolicy::kHalpern);
    CHECK(translated.options.policy.restart_substituted);
    CHECK(translated.options.policy.requested_restart_strategy ==
          kPdlpRestartStrategyAdaptive);
  }

  SECTION("explicit fixed and PID step sizes are preserved") {
    REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                                 kPdlpStepSizeStrategyFixed) ==
            HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult fixed =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    REQUIRE(fixed.ok());
    CHECK(fixed.options.policy.step_size ==
          hipdlp::api::StepSizePolicy::kFixed);
    CHECK_FALSE(fixed.options.policy.step_size_substituted);

    REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                                 kPdlpStepSizeStrategyPid) == HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult pid =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    REQUIRE(pid.ok());
    CHECK(pid.options.policy.step_size == hipdlp::api::StepSizePolicy::kPid);
    CHECK_FALSE(pid.options.policy.step_size_substituted);
  }

  SECTION("Malitsky-Pock is an explicit unsupported-option error") {
    REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                                 kPdlpStepSizeStrategyMalitskyPock) ==
            HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult translated =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    CHECK_FALSE(translated.ok());
    CHECK(translated.error.code == hipdlp::api::ErrorCode::kInvalidOption);
    CHECK(translated.error.message.find("Malitsky-Pock") != std::string::npos);

    // The retained internal legacy selector keeps its own semantics.
    const hipdlp::api::MigrationOptions legacy_selector(
        hipdlp::api::SolverImplementation::kLegacy,
        hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult legacy =
        hipdlp::api::translateOptions(highs.getOptions(), legacy_selector);
    INFO(legacy.error.message);
    CHECK(legacy.ok());
  }

  SECTION("restart requests all resolve to Halpern") {
    const HighsInt requested[] = {
        kPdlpRestartStrategyOff, kPdlpRestartStrategyFixed,
        kPdlpRestartStrategyAdaptive, kPdlpRestartStrategyHalpern};
    for (std::size_t index = 0; index < 4; ++index) {
      REQUIRE(highs.setOptionValue("pdlp_restart_strategy", requested[index]) ==
              HighsStatus::kOk);
      const hipdlp::api::OptionTranslationResult translated =
          hipdlp::api::translateOptions(highs.getOptions(), rewrite);
      INFO(translated.error.message);
      REQUIRE(translated.ok());
      CHECK(translated.options.policy.restart ==
            hipdlp::api::RestartPolicy::kHalpern);
      CHECK(translated.options.policy.restart_substituted ==
            (requested[index] != kPdlpRestartStrategyHalpern));
    }
    REQUIRE(highs.setOptionValue("pdlp_restart_strategy",
                                 kPdlpRestartStrategyHalpern) ==
            HighsStatus::kOk);
    HighsOptions restart_off_options = highs.getOptions();
    restart_off_options.pdlp_features_off = kPdlpRestartOff;
    const hipdlp::api::OptionTranslationResult restart_off =
        hipdlp::api::translateOptions(restart_off_options, rewrite);
    REQUIRE(restart_off.ok());
    CHECK(restart_off.options.policy.restart ==
          hipdlp::api::RestartPolicy::kHalpern);
    CHECK(restart_off.options.policy.restart_substituted);
  }

  SECTION("a migration backend override outranks the public option") {
    REQUIRE(highs.setOptionValue(kHiPdlpBackendString,
                                 kHiPdlpBackendCudaString) == HighsStatus::kOk);
    const hipdlp::api::OptionTranslationResult from_public =
        hipdlp::api::translateOptions(highs.getOptions(), rewrite);
    REQUIRE(from_public.ok());
    CHECK(from_public.options.backend == hipdlp::backend::Kind::kCuda);

    const hipdlp::api::MigrationOptions override_cpu(
        hipdlp::api::SolverImplementation::kRewrite,
        hipdlp::backend::Kind::kCpu);
    const hipdlp::api::OptionTranslationResult from_override =
        hipdlp::api::translateOptions(highs.getOptions(), override_cpu);
    REQUIRE(from_override.ok());
    CHECK(from_override.options.backend == hipdlp::backend::Kind::kCpu);
  }

  SECTION("an invalid stored option string is rejected") {
    HighsOptions options = highs.getOptions();
    options.hipdlp_backend = "gpu";
    CHECK_FALSE(hipdlp::api::translateOptions(options, rewrite).ok());
    options = highs.getOptions();
    options.hipdlp_scaling = "ruiz";
    CHECK_FALSE(hipdlp::api::translateOptions(options, rewrite).ok());
    options = highs.getOptions();
    options.hipdlp_curtis_reid_iterations = 0;
    CHECK_FALSE(hipdlp::api::translateOptions(options, rewrite).ok());
  }
}

TEST_CASE("hi-pdlp-rewrite-solves-with-default-style-highs-options",
          "[pdlp][hipdlp][phase7][solver][options][defaults]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  // Nothing else is changed: adaptive step size and auto scaling are defaults.
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite, hipdlp::backend::Kind::kCpu);
  const hipdlp::api::OptionTranslationResult translated =
      hipdlp::api::translateOptions(highs.getOptions(), rewrite);
  INFO(translated.error.message);
  REQUIRE(translated.ok());
  REQUIRE(translated.options.policy.step_size ==
          hipdlp::api::StepSizePolicy::kFixed);

  const HighsLp lp = makeTinyFeasibleLp();
  RecordingLogSink sink;
  hipdlp::algorithm::SolveContext context;
  context.log = &sink;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      translated.options, hipdlp::model::OriginalModel(lp), context);
  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  REQUIRE(result.has_solution);
  CHECK(result.solution.col_value[0] == Approx(1.0).margin(1e-5));

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution,
                                      translated.options.kkt_tolerance);
  INFO(verification.error);
  CHECK(verification.pass);

  // The adaptive substitution is reported once and is never PID.
  CHECK(sink.countContaining("adaptive") == 1);
  CHECK(sink.countContaining("effective fixed") == 1);
  CHECK(sink.countContaining("PID") == 0);
  CHECK(sink.countContaining("HiPDLP rewrite backend") == 1);
  CHECK(sink.countContaining("HiPDLP rewrite restart") == 1);
  CHECK(sink.countContaining("HiPDLP rewrite scaling") == 1);
  CHECK(sink.countContaining("HiPDLP rewrite finished") == 1);
}

TEST_CASE("hi-pdlp-rewrite-logs-only-every-fortieth-convergence-boundary",
          "[pdlp][hipdlp][phase7][solver][diagnostics]") {
  Highs highs;
  REQUIRE(highs.readModel(std::string(HIGHS_DIR) +
                          "/check/instances/sctest.mps") ==
          HighsStatus::kOk);
  const HighsLp lp = highs.getLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 1601;
  options.optimality_tolerance = std::numeric_limits<double>::denorm_min();
  options.kkt_tolerance = std::numeric_limits<double>::denorm_min();
  options.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
  options.scaling.mode = hipdlp::model::ScalingMode::kAuto;
  RecordingLogSink sink;
  hipdlp::algorithm::SolveContext context;
  context.log = &sink;

  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);

  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.status == hipdlp::api::SolverStatus::kIterationLimit);
  REQUIRE(result.iteration_count == 1601);
  REQUIRE(sink.countContaining("HiPDLP rewrite boundary") == 1);
  std::string boundary;
  for (std::size_t index = 0; index < sink.messages.size(); ++index) {
    if (sink.messages[index].find("HiPDLP rewrite boundary") !=
        std::string::npos) {
      boundary = sink.messages[index];
    }
  }
  CHECK(boundary.find("iteration 1600") != std::string::npos);
  CHECK(boundary.find("relative primal ") != std::string::npos);
  CHECK(boundary.find("dual ") != std::string::npos);
  CHECK(boundary.find("gap ") != std::string::npos);
  CHECK(boundary.find("tau ") != std::string::npos);
  CHECK(boundary.find("sigma ") != std::string::npos);
  CHECK(boundary.find("omega ") != std::string::npos);
  CHECK(boundary.find("eta ") != std::string::npos);
  CHECK(boundary.find("restart condition ") != std::string::npos);
  CHECK(boundary.find("PID outcome ") != std::string::npos);
}

TEST_CASE("hi-pdlp-rewrite-accepts-an-already-optimal-starting-point",
          "[pdlp][hipdlp][phase7][solver][iteration-zero]") {
  const HighsLp lp = makeAlreadyOptimalLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);

  INFO(result.error.message);
  REQUIRE(result.ok());
  CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
  CHECK(result.iteration_count == 0);
  REQUIRE(result.has_solution);
  REQUIRE(result.solution.value_valid);
  REQUIRE(result.solution.dual_valid);
  CHECK(result.solution.col_value[0] == 0.0);

  // The initial point still had to survive postsolve and the KKT verifier.
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, result.solution,
                                      options.kkt_tolerance);
  INFO(verification.error);
  CHECK(verification.pass);
  CHECK(diagnostics.restarts == 0);
  CHECK(diagnostics.resource_creation_calls_during_iterations == 0);
}

TEST_CASE("hi-pdlp-rewrite-allocates-every-backend-vector-before-iterating",
          "[pdlp][hipdlp][phase7][solver][resources]") {
  const HighsLp lp = makeTinyFeasibleLp();
  hipdlp::api::SolverOptions options = makeCpuFixedOptions();
  options.iteration_limit = 200;
  options.optimality_tolerance = 1e-9;
  hipdlp::algorithm::SolveDiagnostics diagnostics;
  hipdlp::algorithm::SolveContext context;
  context.diagnostics = &diagnostics;
  const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
      options, hipdlp::model::OriginalModel(lp), context);

  INFO(result.error.message);
  REQUIRE(result.ok());
  REQUIRE(result.iteration_count > 40);
  // One problem binding, one workspace and sixteen iteration vectors.
  CHECK(diagnostics.resource_creation_calls_before_iterations == 18);
  CHECK(diagnostics.resource_creation_calls_during_iterations == 0);
  CHECK(diagnostics.restarts > 0);
  CHECK(diagnostics.primal_weight_updates == 0);
}

TEST_CASE("hi-pdlp-rewrite-consumes-the-pid-primal-weight-update",
          "[pdlp][hipdlp][phase7][solver][pid][restart]") {
  hipdlp::api::SolverOptions fixed = makeCpuFixedOptions();
  fixed.iteration_limit = 200;
  fixed.optimality_tolerance = 1e-9;
  hipdlp::api::SolverOptions pid = fixed;
  pid.policy.step_size = hipdlp::api::StepSizePolicy::kPid;

  SECTION("a defined residual ratio changes the later trajectory") {
    const HighsLp lp = makeTinyFeasibleLp();
    hipdlp::algorithm::SolveDiagnostics fixed_diagnostics;
    hipdlp::algorithm::SolveContext fixed_context;
    fixed_context.diagnostics = &fixed_diagnostics;
    const hipdlp::api::SolveResult fixed_result = hipdlp::algorithm::solve(
        fixed, hipdlp::model::OriginalModel(lp), fixed_context);
    hipdlp::algorithm::SolveDiagnostics pid_diagnostics;
    hipdlp::algorithm::SolveContext pid_context;
    pid_context.diagnostics = &pid_diagnostics;
    const hipdlp::api::SolveResult pid_result = hipdlp::algorithm::solve(
        pid, hipdlp::model::OriginalModel(lp), pid_context);

    INFO(fixed_result.error.message);
    REQUIRE(fixed_result.ok());
    INFO(pid_result.error.message);
    REQUIRE(pid_result.ok());
    CHECK(fixed_diagnostics.restarts == pid_diagnostics.restarts);
    CHECK(fixed_diagnostics.primal_weight_updates == 0);
    CHECK(pid_diagnostics.primal_weight_updates == pid_diagnostics.restarts);
    CHECK(pid_diagnostics.primal_weight_updates > 0);
    CHECK(pid_result.iteration_count == fixed_result.iteration_count);
    CHECK(pid_result.solution.col_value != fixed_result.solution.col_value);
  }

  SECTION("an undefined residual ratio falls back to the best weight") {
    const HighsLp lp = makeBoxedColumnLp();
    const hipdlp::api::SolveResult fixed_result =
        hipdlp::algorithm::solve(fixed, hipdlp::model::OriginalModel(lp));
    hipdlp::algorithm::SolveDiagnostics pid_diagnostics;
    hipdlp::algorithm::SolveContext pid_context;
    pid_context.diagnostics = &pid_diagnostics;
    const hipdlp::api::SolveResult pid_result = hipdlp::algorithm::solve(
        pid, hipdlp::model::OriginalModel(lp), pid_context);

    INFO(fixed_result.error.message);
    REQUIRE(fixed_result.ok());
    INFO(pid_result.error.message);
    REQUIRE(pid_result.ok());
    CHECK(pid_diagnostics.primal_weight_updates > 0);
    // Reverting to the best weight leaves the steps, and the run, unchanged.
    CHECK(pid_result.iteration_count == fixed_result.iteration_count);
    CHECK(pid_result.solution.col_value == fixed_result.solution.col_value);
    CHECK(pid_result.solution.row_dual == fixed_result.solution.row_dual);
  }
}

TEST_CASE("hi-pdlp-rewrite-postsolves-every-scaling-mode",
          "[pdlp][hipdlp][phase7][solver][scaling][postsolve]") {
  const HighsLp lp = makeScaledFeasibleLp();
  const hipdlp::model::ScalingMode modes[] = {
      hipdlp::model::ScalingMode::kAuto, hipdlp::model::ScalingMode::kOff,
      hipdlp::model::ScalingMode::kCurtisReid,
      hipdlp::model::ScalingMode::kLegacy};
  for (std::size_t index = 0; index < 4; ++index) {
    hipdlp::api::SolverOptions options = makeCpuFixedOptions();
    options.scaling.mode = modes[index];
    const hipdlp::api::SolveResult result =
        hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
    REQUIRE(result.has_solution);
    CHECK(result.solution.col_value[0] ==
          Approx(2.0).margin(4.0 * kSolverTolerance));
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                        options.kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }
}

TEST_CASE("hi-pdlp-rewrite-uses-only-the-remaining-global-time-budget",
          "[pdlp][hipdlp][phase7][solver][time-limit][api]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                               kPdlpStepSizeStrategyFixed) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_optimality_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("kkt_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  const HighsLp lp = makeTinyFeasibleLp();
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite, hipdlp::backend::Kind::kCpu);

  SECTION("an exhausted global budget leaves nothing for the rewrite") {
    REQUIRE(highs.setOptionValue("time_limit", 0.005) == HighsStatus::kOk);
    HighsTimer timer;
    timer.setPrintfFlag(false, false);
    timer.start();
    while (timer.read() < 0.01) {
    }
    HighsBasis basis;
    basis.valid = true;
    HighsSolution solution;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    HighsCallback callback(&highs);
    REQUIRE(hipdlp::api::solve(highs.getOptions(), timer, lp, basis, solution,
                               model_status, info, callback,
                               rewrite) == HighsStatus::kOk);
    timer.stop();
    CHECK(model_status == HighsModelStatus::kTimeLimit);
    CHECK(info.pdlp_iteration_count == 0);
    REQUIRE(solution.value_valid);
    REQUIRE(solution.dual_valid);
  }

  SECTION("the infinite sentinel is preserved") {
    REQUIRE(highs.getOptions().time_limit == kHighsInf);
    HighsTimer timer;
    timer.setPrintfFlag(false, false);
    timer.start();
    while (timer.read() < 0.001) {
    }
    HighsBasis basis;
    basis.valid = true;
    HighsSolution solution;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    HighsCallback callback(&highs);
    REQUIRE(hipdlp::api::solve(highs.getOptions(), timer, lp, basis, solution,
                               model_status, info, callback,
                               rewrite) == HighsStatus::kOk);
    timer.stop();
    CHECK(model_status == HighsModelStatus::kOptimal);
    CHECK(info.pdlp_iteration_count > 0);
  }
}

TEST_CASE("hi-pdlp-rewrite-emits-no-user-log-when-output-is-off",
          "[pdlp][hipdlp][phase7][solver][logging]") {
  const HighsLp lp = makeTinyFeasibleLp();
  const hipdlp::api::MigrationOptions rewrite(
      hipdlp::api::SolverImplementation::kRewrite, hipdlp::backend::Kind::kCpu);
  const bool output_flags[] = {false, true};
  const char* log_files[] = {"hipdlp-rewrite-log-off.log",
                             "hipdlp-rewrite-log-on.log"};
  for (std::size_t index = 0; index < 2; ++index) {
    std::remove(log_files[index]);
    Highs highs;
    REQUIRE(highs.setOptionValue("log_to_console", false) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue(kLogFileString, log_files[index]) ==
            HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("output_flag", output_flags[index]) ==
            HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("pdlp_optimality_tolerance",
                                 kSolverTolerance) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("kkt_tolerance", kSolverTolerance) ==
            HighsStatus::kOk);
    HighsTimer timer;
    timer.setPrintfFlag(false, false);
    HighsBasis basis;
    basis.valid = true;
    HighsSolution solution;
    HighsModelStatus model_status = HighsModelStatus::kNotset;
    HighsInfo info;
    HighsCallback callback(&highs);
    REQUIRE(hipdlp::api::solve(highs.getOptions(), timer, lp, basis, solution,
                               model_status, info, callback,
                               rewrite) == HighsStatus::kOk);
    CHECK(model_status == HighsModelStatus::kOptimal);
    REQUIRE(highs.setOptionValue(kLogFileString, "") == HighsStatus::kOk);

    std::ifstream stream(log_files[index]);
    const std::string contents((std::istreambuf_iterator<char>(stream)),
                               std::istreambuf_iterator<char>());
    stream.close();
    if (output_flags[index]) {
      CHECK(contents.find("HiPDLP rewrite") != std::string::npos);
    } else {
      CHECK(contents.empty());
    }
    std::remove(log_files[index]);
  }
}

TEST_CASE("hi-pdlp-public-path-never-enters-pdhg-for-a-zero-operator-lp",
          "[pdlp][hipdlp][phase7][solver][trivial]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", false) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("solver", kHiPdlpString) == HighsStatus::kOk);

  SECTION("a zero constraint matrix is solved directly") {
    HighsLp lp = makeTinyFeasibleLp();
    lp.row_lower_ = {-kHighsInf};
    lp.row_upper_ = {kHighsInf};
    lp.a_matrix_.start_ = {0, 0};
    lp.a_matrix_.index_.clear();
    lp.a_matrix_.value_.clear();
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    CHECK(highs.getModelStatus() == HighsModelStatus::kOptimal);
    CHECK(highs.getInfo().pdlp_iteration_count == 0);
    CHECK(highs.getSolution().col_value[0] == Approx(0.0).margin(1e-12));
  }

  SECTION("a row-free model is solved directly") {
    HighsLp lp;
    lp.num_col_ = 1;
    lp.num_row_ = 0;
    lp.col_cost_ = {1.0};
    lp.col_lower_ = {0.0};
    lp.col_upper_ = {kHighsInf};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 1;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0};
    lp.sense_ = ObjSense::kMinimize;
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    CHECK(highs.getModelStatus() == HighsModelStatus::kOptimal);
    CHECK(highs.getInfo().pdlp_iteration_count == 0);
  }
}

TEST_CASE("hi-pdlp-public-path-routes-explicit-cuda-pid-to-rewrite",
          "[pdlp][hipdlp][phase7][solver][public][cuda][pid]") {
  const char* log_file = "hipdlp-public-cuda-pid.log";
  std::remove(log_file);
  Highs highs;
  REQUIRE(highs.setOptionValue("log_to_console", false) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue(kLogFileString, log_file) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("output_flag", true) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("solver", kHiPdlpString) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("hipdlp_backend", "cuda") ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("hipdlp_scaling", "auto") ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                               kPdlpStepSizeStrategyPid) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_optimality_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("kkt_tolerance", kSolverTolerance) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("time_limit", 60.0) == HighsStatus::kOk);
  REQUIRE(highs.readModel(std::string(HIGHS_DIR) +
                          "/check/instances/sctest.mps") ==
          HighsStatus::kOk);

  const bool cuda_usable = cudaBackendUsable();
  const HighsStatus run_status = highs.run();
  REQUIRE(highs.setOptionValue(kLogFileString, "") == HighsStatus::kOk);
  std::ifstream stream(log_file);
  const std::string contents((std::istreambuf_iterator<char>(stream)),
                             std::istreambuf_iterator<char>());
  stream.close();
  std::remove(log_file);

  if (!cuda_usable) {
    CHECK(run_status == HighsStatus::kError);
    CHECK(highs.getModelStatus() == HighsModelStatus::kSolveError);
    CHECK(contents.find("requested HiPDLP backend is unavailable") !=
          std::string::npos);
    CHECK(contents.find("effective cpu") == std::string::npos);
    return;
  }

  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  const HighsLp& lp = highs.getLp();
  const HighsSolution& solution = highs.getSolution();
  REQUIRE(solution.value_valid);
  REQUIRE(solution.dual_valid);
  CHECK(solution.col_value.size() == static_cast<std::size_t>(lp.num_col_));
  CHECK(solution.row_value.size() == static_cast<std::size_t>(lp.num_row_));
  CHECK(solution.col_dual.size() == static_cast<std::size_t>(lp.num_col_));
  CHECK(solution.row_dual.size() == static_cast<std::size_t>(lp.num_row_));
  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(lp, solution, kSolverTolerance);
  INFO(verification.error);
  CHECK(verification.pass);
  CHECK(contents.find(
            "HiPDLP rewrite backend: requested cuda, effective cuda") !=
        std::string::npos);
  CHECK(contents.find("effective PID") != std::string::npos);
}

TEST_CASE("hi-pdlp-rewrite-runs-one-control-loop-on-every-backend",
          "[pdlp][hipdlp][phase7][solver][backend][cuda]") {
  const bool cuda_usable = cudaBackendUsable();
  const HighsLp lp = makeScaledFeasibleLp();

  hipdlp::api::SolverOptions cpu_options = makeCpuFixedOptions();
  cpu_options.scaling.mode = hipdlp::model::ScalingMode::kCurtisReid;
  const hipdlp::api::SolveResult cpu_result =
      hipdlp::algorithm::solve(cpu_options, hipdlp::model::OriginalModel(lp));
  INFO(cpu_result.error.message);
  REQUIRE(cpu_result.ok());
  REQUIRE(cpu_result.status == hipdlp::api::SolverStatus::kOptimal);

  SECTION("an unimplemented backend is never a silent fallback") {
    const hipdlp::backend::Kind unavailable[] = {hipdlp::backend::Kind::kOpenmp,
                                                 hipdlp::backend::Kind::kRocm};
    for (std::size_t index = 0; index < 2; ++index) {
      hipdlp::api::SolverOptions options = cpu_options;
      options.backend = unavailable[index];
      const hipdlp::api::SolveResult result =
          hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
      CHECK_FALSE(result.ok());
      CHECK(result.status == hipdlp::api::SolverStatus::kError);
      CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
      CHECK_FALSE(result.has_solution);
      CHECK(result.iteration_count == 0);
    }
  }

  SECTION("auto selects CUDA when it is usable and CPU otherwise") {
    hipdlp::api::SolverOptions options = cpu_options;
    options.backend = hipdlp::backend::Kind::kAuto;
    RecordingLogSink sink;
    hipdlp::algorithm::SolveContext context;
    context.log = &sink;
    const hipdlp::api::SolveResult result = hipdlp::algorithm::solve(
        options, hipdlp::model::OriginalModel(lp), context);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == cpu_result.status);
    CHECK(result.iteration_count == cpu_result.iteration_count);
    REQUIRE(result.has_solution);
    CHECK(result.solution.col_value[0] ==
          Approx(cpu_result.solution.col_value[0]).margin(kCudaTolerance));
    CHECK(result.solution.row_dual[0] ==
          Approx(cpu_result.solution.row_dual[0]).margin(kCudaTolerance));
    CHECK(sink.countContaining("requested auto, effective ") == 1);
    CHECK(sink.countContaining(cuda_usable ? "effective cuda"
                                           : "effective cpu") == 1);
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                        options.kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }

  SECTION("explicit CUDA either runs the same loop or reports unavailable") {
    hipdlp::api::SolverOptions options = cpu_options;
    options.backend = hipdlp::backend::Kind::kCuda;
    const hipdlp::api::SolveResult result =
        hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
    if (!cuda_usable) {
      CHECK_FALSE(result.ok());
      CHECK(result.status == hipdlp::api::SolverStatus::kError);
      CHECK(result.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
      CHECK_FALSE(result.has_solution);
      CHECK(result.iteration_count == 0);
      return;
    }
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == cpu_result.status);
    CHECK(result.iteration_count == cpu_result.iteration_count);
    REQUIRE(result.has_solution);
    CHECK(result.solution.col_value[0] ==
          Approx(cpu_result.solution.col_value[0]).margin(kCudaTolerance));
    CHECK(result.solution.row_value[0] ==
          Approx(cpu_result.solution.row_value[0]).margin(kCudaTolerance));
    CHECK(result.solution.row_dual[0] ==
          Approx(cpu_result.solution.row_dual[0]).margin(kCudaTolerance));
    CHECK(result.solution.col_dual[0] ==
          Approx(cpu_result.solution.col_dual[0]).margin(kCudaTolerance));
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                        options.kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }

  SECTION("explicit CUDA takes the same limit and restart decisions") {
    if (!cuda_usable) return;
    const HighsLp tiny = makeTinyFeasibleLp();

    hipdlp::api::SolverOptions cpu_limit = makeCpuFixedOptions();
    cpu_limit.iteration_limit = 137;
    cpu_limit.optimality_tolerance = 1e-16;
    hipdlp::algorithm::SolveDiagnostics cpu_diagnostics;
    hipdlp::algorithm::SolveContext cpu_context;
    cpu_context.diagnostics = &cpu_diagnostics;
    const hipdlp::api::SolveResult cpu_limited = hipdlp::algorithm::solve(
        cpu_limit, hipdlp::model::OriginalModel(tiny), cpu_context);

    hipdlp::api::SolverOptions cuda_limit = cpu_limit;
    cuda_limit.backend = hipdlp::backend::Kind::kCuda;
    hipdlp::algorithm::SolveDiagnostics cuda_diagnostics;
    hipdlp::algorithm::SolveContext cuda_context;
    cuda_context.diagnostics = &cuda_diagnostics;
    const hipdlp::api::SolveResult cuda_limited = hipdlp::algorithm::solve(
        cuda_limit, hipdlp::model::OriginalModel(tiny), cuda_context);

    INFO(cpu_limited.error.message);
    REQUIRE(cpu_limited.ok());
    INFO(cuda_limited.error.message);
    REQUIRE(cuda_limited.ok());
    // A non-multiple-of-40 limit is honoured exactly on both backends.
    CHECK(cpu_limited.status == hipdlp::api::SolverStatus::kIterationLimit);
    CHECK(cuda_limited.status == cpu_limited.status);
    CHECK(cpu_limited.iteration_count == 137);
    CHECK(cuda_limited.iteration_count == 137);
    CHECK(cuda_diagnostics.restarts == cpu_diagnostics.restarts);
    CHECK(cpu_diagnostics.restarts > 0);
    CHECK(cuda_diagnostics.resource_creation_calls_during_iterations == 0);
    CHECK(cuda_limited.solution.col_value[0] ==
          Approx(cpu_limited.solution.col_value[0]).margin(kCudaTolerance));

    hipdlp::api::SolverOptions cpu_pid = cpu_limit;
    cpu_pid.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
    hipdlp::api::SolverOptions cuda_pid = cuda_limit;
    cuda_pid.policy.step_size = hipdlp::api::StepSizePolicy::kPid;
    hipdlp::algorithm::SolveDiagnostics cpu_pid_diagnostics;
    hipdlp::algorithm::SolveContext cpu_pid_context;
    cpu_pid_context.diagnostics = &cpu_pid_diagnostics;
    const hipdlp::api::SolveResult cpu_pid_result = hipdlp::algorithm::solve(
        cpu_pid, hipdlp::model::OriginalModel(tiny), cpu_pid_context);
    hipdlp::algorithm::SolveDiagnostics cuda_pid_diagnostics;
    hipdlp::algorithm::SolveContext cuda_pid_context;
    cuda_pid_context.diagnostics = &cuda_pid_diagnostics;
    const hipdlp::api::SolveResult cuda_pid_result = hipdlp::algorithm::solve(
        cuda_pid, hipdlp::model::OriginalModel(tiny), cuda_pid_context);
    INFO(cpu_pid_result.error.message);
    REQUIRE(cpu_pid_result.ok());
    INFO(cuda_pid_result.error.message);
    REQUIRE(cuda_pid_result.ok());
    CHECK(cpu_pid_diagnostics.primal_weight_updates > 0);
    CHECK(cuda_pid_diagnostics.primal_weight_updates ==
          cpu_pid_diagnostics.primal_weight_updates);
    CHECK(cuda_pid_result.status == cpu_pid_result.status);
    CHECK(cuda_pid_result.iteration_count == cpu_pid_result.iteration_count);
    CHECK(cuda_pid_result.solution.col_value[0] ==
          Approx(cpu_pid_result.solution.col_value[0]).margin(kCudaTolerance));
  }

  SECTION("explicit CUDA honours a zero time budget without iterating") {
    if (!cuda_usable) return;
    hipdlp::api::SolverOptions options = cpu_options;
    options.backend = hipdlp::backend::Kind::kCuda;
    options.time_limit = 0.0;
    const hipdlp::api::SolveResult result =
        hipdlp::algorithm::solve(options, hipdlp::model::OriginalModel(lp));
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kTimeLimit);
    CHECK(result.iteration_count == 0);
    REQUIRE(result.has_solution);
    REQUIRE(result.solution.value_valid);
    REQUIRE(result.solution.dual_valid);
  }
}
