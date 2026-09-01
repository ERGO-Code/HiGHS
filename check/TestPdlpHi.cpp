#include <algorithm>
#include <chrono>
#include <limits>
#include <type_traits>
#include <vector>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HighsSolution.h"
#include "pdlp/hipdlp/algorithm/solver.hpp"
#include "pdlp/hipdlp/api/solver.hpp"
#include "pdlp/hipdlp/backend/backend_factory.hpp"
#include "pdlp/hipdlp/verification/kkt_verifier.hpp"

#ifdef CUPDLP_GPU
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#endif

const bool dev_run = false;
const double kkt_tolerance = 1e-4;

namespace {

HighsLp makeOneRowLp(const ObjSense sense, const double col_lower,
                     const double col_upper, const double row_lower,
                     const double row_upper) {
  HighsLp lp;
  lp.num_col_ = 1;
  lp.num_row_ = 1;
  lp.col_cost_ = {2.0};
  lp.col_lower_ = {col_lower};
  lp.col_upper_ = {col_upper};
  lp.row_lower_ = {row_lower};
  lp.row_upper_ = {row_upper};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 1;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1};
  lp.a_matrix_.index_ = {0};
  lp.a_matrix_.value_ = {1.0};
  lp.sense_ = sense;
  lp.offset_ = 3.0;
  return lp;
}

HighsSolution makeOneRowSolution(const double col_value, const double col_dual,
                                 const double row_dual) {
  HighsSolution solution;
  solution.value_valid = true;
  solution.dual_valid = true;
  solution.col_value = {col_value};
  solution.col_dual = {col_dual};
  solution.row_value = {col_value};
  solution.row_dual = {row_dual};
  return solution;
}

HighsInfo getCoreKktInfo(const Highs& highs) {
  HighsInfo info = highs.getInfo();
  HighsPrimalDualErrors errors;
  getLpKktFailures(highs.getOptions(), highs.getLp(), highs.getSolution(),
                   highs.getBasis(), info, errors, true);
  return info;
}

bool coreKktPass(const HighsInfo& info, const double tolerance) {
  return std::max({info.max_relative_primal_infeasibility,
                   info.max_relative_dual_infeasibility,
                   info.max_relative_primal_residual_error,
                   info.max_relative_dual_residual_error,
                   info.primal_dual_objective_error}) <= tolerance;
}

#ifndef CUPDLP_GPU
void configureLegacyHiPdlp(Highs& highs, const std::string& model,
                           const HighsInt iteration_limit) {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.setOptionValue("output_flag", dev_run) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("solver", kHiPdlpString) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("kkt_tolerance", kkt_tolerance) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("pdlp_iteration_limit", iteration_limit) ==
          HighsStatus::kOk);
  REQUIRE(highs.readModel(model_file) != HighsStatus::kError);
}
#endif

}  // namespace

TEST_CASE("hi-pdlp-migration-boundary", "[pdlp][hipdlp][migration][phase4]") {
  const hipdlp::api::MigrationOptions defaults =
      hipdlp::api::defaultMigrationOptions();
  CHECK(defaults.implementation == hipdlp::api::SolverImplementation::kLegacy);
  CHECK(defaults.backend == hipdlp::backend::Kind::kAuto);

  Highs highs;
  const hipdlp::api::OptionTranslationResult translated =
      hipdlp::api::translateOptions(highs.getOptions(), defaults);
  INFO(translated.error.message);
  REQUIRE(translated.ok());
  CHECK(translated.options.implementation ==
        hipdlp::api::SolverImplementation::kLegacy);
  CHECK(translated.options.iteration_limit ==
        highs.getOptions().pdlp_iteration_limit);
  CHECK(translated.options.scaling.mode == hipdlp::model::ScalingMode::kAuto);
  CHECK(translated.options.scaling.curtis_reid_sweeps == 20);
  CHECK(translated.options.scaling.legacy_scaling_mode ==
        highs.getOptions().pdlp_scaling_mode);
  CHECK(translated.options.scaling.ruiz_iterations ==
        highs.getOptions().pdlp_ruiz_iterations);
  CHECK(translated.options.scaling.features_off ==
        highs.getOptions().pdlp_features_off);

  SECTION("rewrite solves through the algorithm boundary") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    const hipdlp::model::OriginalModel original_model(lp);
    hipdlp::api::SolverOptions rewrite_options = translated.options;
    rewrite_options.implementation =
        hipdlp::api::SolverImplementation::kRewrite;
    rewrite_options.backend = hipdlp::backend::Kind::kCpu;
    rewrite_options.policy.step_size = hipdlp::api::StepSizePolicy::kFixed;
    const hipdlp::api::SolveResult result =
        hipdlp::algorithm::solve(rewrite_options, original_model);
    INFO(result.error.message);
    REQUIRE(result.ok());
    CHECK(result.status == hipdlp::api::SolverStatus::kOptimal);
    REQUIRE(result.has_solution);
    REQUIRE(result.solution.value_valid);
    REQUIRE(result.solution.dual_valid);
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, result.solution,
                                        rewrite_options.kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }

  SECTION("rewrite maps a solution through the internal API") {
    REQUIRE(highs.setOptionValue("pdlp_step_size_strategy",
                                 kPdlpStepSizeStrategyFixed) ==
            HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("pdlp_iteration_limit", 10000) ==
            HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("pdlp_optimality_tolerance", kkt_tolerance) ==
            HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("kkt_tolerance", kkt_tolerance) ==
            HighsStatus::kOk);
    HighsTimer timer;
    HighsBasis basis;
    basis.valid = true;
    HighsSolution solution = makeOneRowSolution(1.0, 0.0, 2.0);
    HighsModelStatus model_status = HighsModelStatus::kOptimal;
    HighsInfo info;
    info.pdlp_iteration_count = 99;
    HighsCallback callback(&highs);
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    const hipdlp::api::MigrationOptions rewrite(
        hipdlp::api::SolverImplementation::kRewrite,
        hipdlp::backend::Kind::kCpu);

    REQUIRE(hipdlp::api::solve(highs.getOptions(), timer, lp, basis, solution,
                               model_status, info, callback,
                               rewrite) == HighsStatus::kOk);
    CHECK(model_status == HighsModelStatus::kOptimal);
    CHECK(info.pdlp_iteration_count > 0);
    CHECK_FALSE(basis.valid);
    REQUIRE(solution.value_valid);
    REQUIRE(solution.dual_valid);
    const hipdlp::verification::KktVerificationResult verification =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    INFO(verification.error);
    CHECK(verification.pass);
  }

  SECTION("invalid internal selectors are rejected") {
    const hipdlp::api::MigrationOptions invalid_implementation(
        static_cast<hipdlp::api::SolverImplementation>(-1),
        hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult invalid_implementation_result =
        hipdlp::api::translateOptions(highs.getOptions(),
                                      invalid_implementation);
    CHECK_FALSE(invalid_implementation_result.ok());
    CHECK(invalid_implementation_result.error.code ==
          hipdlp::api::ErrorCode::kInvalidOption);

    const hipdlp::api::MigrationOptions invalid_backend(
        hipdlp::api::SolverImplementation::kRewrite,
        static_cast<hipdlp::backend::Kind>(-1));
    const hipdlp::api::OptionTranslationResult invalid_backend_result =
        hipdlp::api::translateOptions(highs.getOptions(), invalid_backend);
    CHECK_FALSE(invalid_backend_result.ok());
    CHECK(invalid_backend_result.error.code ==
          hipdlp::api::ErrorCode::kInvalidOption);

    const hipdlp::api::MigrationOptions legacy_backend(
        hipdlp::api::SolverImplementation::kLegacy,
        hipdlp::backend::Kind::kCpu);
    CHECK_FALSE(
        hipdlp::api::translateOptions(highs.getOptions(), legacy_backend).ok());

    HighsOptions invalid_options = highs.getOptions();
    invalid_options.pdlp_scaling_mode = kPdlpScalingMax + 1;
    const hipdlp::api::MigrationOptions rewrite(
        hipdlp::api::SolverImplementation::kRewrite,
        hipdlp::backend::Kind::kAuto);
    const hipdlp::api::OptionTranslationResult invalid_option_result =
        hipdlp::api::translateOptions(invalid_options, rewrite);
    CHECK_FALSE(invalid_option_result.ok());
    CHECK(invalid_option_result.error.code ==
          hipdlp::api::ErrorCode::kInvalidOption);
  }

  SECTION("result statuses map without inventing a solution") {
    struct StatusCase {
      hipdlp::api::SolverStatus solver_status;
      HighsModelStatus model_status;
    };
    const StatusCase cases[] = {
        {hipdlp::api::SolverStatus::kOptimal, HighsModelStatus::kOptimal},
        {hipdlp::api::SolverStatus::kInfeasible, HighsModelStatus::kInfeasible},
        {hipdlp::api::SolverStatus::kUnbounded, HighsModelStatus::kUnbounded},
        {hipdlp::api::SolverStatus::kIterationLimit,
         HighsModelStatus::kIterationLimit},
        {hipdlp::api::SolverStatus::kTimeLimit, HighsModelStatus::kTimeLimit}};
    const HighsSolution mapped_solution = makeOneRowSolution(1.0, 0.0, 2.0);
    for (const StatusCase& test_case : cases) {
      HighsBasis basis;
      HighsSolution solution;
      HighsModelStatus model_status = HighsModelStatus::kNotset;
      HighsInfo info;
      const hipdlp::api::SolveResult result = hipdlp::api::SolveResult::success(
          test_case.solver_status, mapped_solution, 17);
      REQUIRE(hipdlp::api::applySolveResult(result, basis, solution,
                                            model_status,
                                            info) == HighsStatus::kOk);
      CHECK(model_status == test_case.model_status);
      CHECK(solution.value_valid);
      CHECK(info.pdlp_iteration_count == 17);
    }

    HighsBasis basis;
    HighsSolution solution = mapped_solution;
    HighsModelStatus model_status = HighsModelStatus::kOptimal;
    HighsInfo info;
    const hipdlp::api::SolveResult not_ready =
        hipdlp::api::SolveResult::failure(hipdlp::api::ErrorCode::kNotReady,
                                          "not ready");
    CHECK(hipdlp::api::applySolveResult(not_ready, basis, solution,
                                        model_status,
                                        info) == HighsStatus::kError);
    CHECK(model_status == HighsModelStatus::kSolveError);
    CHECK_FALSE(solution.value_valid);
  }

  SECTION("backend contracts are move-only with runtime selection") {
    CHECK_FALSE(
        std::is_copy_constructible<hipdlp::backend::ProblemHandle>::value);
    CHECK(std::is_move_constructible<hipdlp::backend::ProblemHandle>::value);
    CHECK_FALSE(
        std::is_copy_constructible<hipdlp::backend::VectorHandle>::value);
    CHECK(std::is_move_constructible<hipdlp::backend::VectorHandle>::value);
    CHECK_FALSE(
        std::is_copy_constructible<hipdlp::backend::WorkspaceHandle>::value);
    CHECK(std::is_move_constructible<hipdlp::backend::WorkspaceHandle>::value);
    CHECK_FALSE(std::is_copy_constructible<hipdlp::model::Problem>::value);
    CHECK(std::is_move_constructible<hipdlp::model::Problem>::value);

    const hipdlp::backend::PrecisionDescriptor precision;
    CHECK(precision.valid());
    CHECK(precision.isFp64());

    std::size_t registration_count = 0;
    const hipdlp::backend::Descriptor* registrations =
        hipdlp::backend::registrations(registration_count);
    REQUIRE(registrations != nullptr);
    CHECK(registration_count == 4);
    for (std::size_t index = 0; index < registration_count; ++index) {
      if (registrations[index].kind == hipdlp::backend::Kind::kCpu
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
          || registrations[index].kind == hipdlp::backend::Kind::kCuda
#endif
      ) {
        CHECK(registrations[index].implemented);
        CHECK(registrations[index].capabilities.mask != 0);
      } else {
        CHECK_FALSE(registrations[index].implemented);
        CHECK(registrations[index].capabilities.mask == 0);
      }
    }

    const hipdlp::backend::FactoryResult cpu =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kCpu);
    CHECK(cpu.ok());
    REQUIRE(cpu.backend != nullptr);
    CHECK(cpu.backend->descriptor().kind == hipdlp::backend::Kind::kCpu);
    bool cuda_usable = false;
#ifdef HIGHS_HIPDLP_CUDA_DETECTED
    const hipdlp::backend::FactoryResult cuda =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kCuda);
    cuda_usable = cuda.ok();
#endif
    const hipdlp::backend::FactoryResult automatic =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kAuto);
    CHECK(automatic.ok());
    REQUIRE(automatic.backend != nullptr);
    CHECK(automatic.backend->descriptor().kind ==
          (cuda_usable ? hipdlp::backend::Kind::kCuda
                       : hipdlp::backend::Kind::kCpu));
    const hipdlp::backend::FactoryResult openmp =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kOpenmp);
    CHECK_FALSE(openmp.ok());
    CHECK(openmp.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
    const hipdlp::backend::FactoryResult rocm =
        hipdlp::backend::createBackend(hipdlp::backend::Kind::kRocm);
    CHECK_FALSE(rocm.ok());
    CHECK(rocm.error.code == hipdlp::api::ErrorCode::kBackendUnavailable);
  }
}

TEST_CASE("hi-pdlp-kkt-verifier", "[pdlp][hipdlp][verification]") {
  SECTION("minimization with offset and lower row") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    const HighsSolution solution = makeOneRowSolution(1.0, 0.0, 2.0);
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK(result.pass);
    CHECK(result.metrics.primal_objective == 5.0);
    CHECK(result.metrics.dual_objective == 5.0);
    CHECK(result.metrics.max_complementarity == 0.0);
  }

  SECTION("maximization with offset and boxed row") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMaximize, 0.0, kHighsInf, 0.0, 1.0);
    const HighsSolution solution = makeOneRowSolution(1.0, 0.0, 2.0);
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK(result.pass);
    CHECK(result.metrics.primal_objective == 5.0);
    CHECK(result.metrics.dual_objective == 5.0);
  }

  SECTION("free column and fixed row") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, -kHighsInf, kHighsInf, 2.0, 2.0);
    const HighsSolution solution = makeOneRowSolution(2.0, 0.0, 2.0);
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK(result.pass);
    CHECK(result.metrics.primal_objective == 7.0);
    CHECK(result.metrics.dual_objective == 7.0);
  }

  SECTION("fixed column and empty row set") {
    HighsLp lp;
    lp.num_col_ = 1;
    lp.num_row_ = 0;
    lp.col_cost_ = {2.0};
    lp.col_lower_ = {1.0};
    lp.col_upper_ = {1.0};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 1;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 3.0;

    HighsSolution solution;
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_value = {1.0};
    solution.col_dual = {2.0};

    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK(result.pass);
    CHECK(result.metrics.primal_objective == 5.0);
    CHECK(result.metrics.dual_objective == 5.0);
  }

  SECTION("boxed column has the correct dual sign at its upper bound") {
    HighsLp lp;
    lp.num_col_ = 1;
    lp.num_row_ = 0;
    lp.col_cost_ = {1.0};
    lp.col_lower_ = {0.0};
    lp.col_upper_ = {2.0};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 1;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0};
    lp.sense_ = ObjSense::kMinimize;

    HighsSolution solution;
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_value = {2.0};
    solution.col_dual = {1.0};

    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.pass);
    CHECK_FALSE(result.dual_feasible);
    CHECK(result.metrics.max_dual_violation == 1.0);
    CHECK(result.worst_dual_col == 0);
  }

  SECTION("narrow boxed column cannot hide complementarity in the gap") {
    HighsLp lp;
    lp.num_col_ = 1;
    lp.num_row_ = 0;
    lp.col_cost_ = {-1e9};
    lp.col_lower_ = {0.0};
    lp.col_upper_ = {0.005};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 1;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0};
    lp.sense_ = ObjSense::kMinimize;

    HighsSolution solution;
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_value = {0.0};
    solution.col_dual = {-1e9};

    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.pass);
    CHECK_FALSE(result.objective_gap_ok);
    CHECK(result.metrics.max_complementarity == 5e6);
    CHECK(result.metrics.dual_objective == -5e6);
  }

  SECTION("primal violation identifies the original row") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    const HighsSolution solution = makeOneRowSolution(0.5, 0.0, 2.0);
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.pass);
    CHECK_FALSE(result.primal_feasible);
    CHECK(result.worst_primal_row == 0);
    CHECK(result.metrics.max_primal_violation == 0.5);
  }

  SECTION("row activity residual identifies the original row") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    HighsSolution solution = makeOneRowSolution(1.0, 0.0, 2.0);
    solution.row_value[0] = 1.5;
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.pass);
    CHECK_FALSE(result.primal_feasible);
    CHECK(result.worst_primal_row == 0);
    CHECK(result.metrics.max_row_activity_residual == 0.5);
  }

  SECTION("inactive large bound does not dilute primal violation") {
    HighsLp lp;
    lp.num_col_ = 2;
    lp.num_row_ = 0;
    lp.col_cost_ = {0.0, 0.0};
    lp.col_lower_ = {0.0, 1.0};
    lp.col_upper_ = {1e12, kHighsInf};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 2;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0, 0};
    lp.sense_ = ObjSense::kMinimize;

    HighsSolution solution;
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_value = {0.0, 0.999};
    solution.col_dual = {0.0, 0.0};

    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.primal_feasible);
    CHECK(result.metrics.relative_primal_residual == Approx(0.0005));
    CHECK(result.worst_primal_col == 1);
  }

  SECTION("inactive large cost does not dilute dual violation") {
    HighsLp lp;
    lp.num_col_ = 2;
    lp.num_row_ = 0;
    lp.col_cost_ = {1e12, -0.001};
    lp.col_lower_ = {0.0, 0.0};
    lp.col_upper_ = {1.0, kHighsInf};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 2;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0, 0};
    lp.sense_ = ObjSense::kMinimize;

    HighsSolution solution;
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_value = {0.0, 0.0};
    solution.col_dual = {1e12, -0.001};

    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.dual_feasible);
    CHECK(result.metrics.relative_dual_residual > kkt_tolerance);
    CHECK(result.worst_dual_col == 1);
  }

  SECTION("stationarity violation identifies the original column") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    const HighsSolution solution = makeOneRowSolution(1.0, 1.0, 2.0);
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.pass);
    CHECK_FALSE(result.dual_feasible);
    CHECK(result.worst_dual_col == 0);
    CHECK(result.metrics.max_stationarity_residual == 1.0);
  }

  SECTION("worst dual column identifies the largest dual defect") {
    HighsLp lp;
    lp.num_col_ = 2;
    lp.num_row_ = 0;
    lp.col_cost_ = {10.0, 0.0};
    lp.col_lower_ = {-kHighsInf, 0.0};
    lp.col_upper_ = {0.0, 0.0};
    lp.a_matrix_.format_ = MatrixFormat::kColwise;
    lp.a_matrix_.num_col_ = 2;
    lp.a_matrix_.num_row_ = 0;
    lp.a_matrix_.start_ = {0, 0, 0};
    lp.sense_ = ObjSense::kMinimize;
    lp.offset_ = 0.0;

    HighsSolution solution;
    solution.value_valid = true;
    solution.dual_valid = true;
    solution.col_value = {0.0, 0.0};
    solution.col_dual = {10.0, 1.0};

    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK_FALSE(result.dual_feasible);
    CHECK(result.metrics.max_dual_violation == 10.0);
    CHECK(result.metrics.max_stationarity_residual == 1.0);
    CHECK(result.worst_dual_col == 0);
  }

  SECTION("worst complementarity index identifies its original kind") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    const HighsSolution solution = makeOneRowSolution(10.0, 0.1, 1.9);
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    REQUIRE(result.input_valid);
    REQUIRE(result.metrics_finite);
    CHECK(result.metrics.max_complementarity == Approx(17.1));
    CHECK(result.worst_complementarity_col == -1);
    CHECK(result.worst_complementarity_row == 0);
  }

  SECTION("non-finite solution is rejected before sparse products") {
    const HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    HighsSolution solution = makeOneRowSolution(1.0, 0.0, 2.0);
    solution.col_value[0] = std::numeric_limits<double>::quiet_NaN();
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    CHECK_FALSE(result.input_valid);
    CHECK_FALSE(result.metrics_finite);
    CHECK_FALSE(result.pass);
  }

  SECTION("wrong-direction infinite bound is rejected") {
    HighsLp lp =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    lp.col_lower_[0] = std::numeric_limits<double>::infinity();
    const HighsSolution solution = makeOneRowSolution(1.0, 0.0, 2.0);
    const auto result =
        hipdlp::verification::verifyKkt(lp, solution, kkt_tolerance);
    CHECK_FALSE(result.input_valid);
    CHECK_FALSE(result.metrics_finite);
    CHECK_FALSE(result.pass);
  }

  SECTION("malformed sparse storage is rejected before products") {
    HighsLp short_start =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    short_start.a_matrix_.start_ = {0};
    const HighsSolution solution = makeOneRowSolution(1.0, 0.0, 2.0);
    CHECK_FALSE(
        hipdlp::verification::verifyKkt(short_start, solution, kkt_tolerance)
            .input_valid);

    HighsLp invalid_index =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    invalid_index.a_matrix_.index_[0] = 1;
    CHECK_FALSE(
        hipdlp::verification::verifyKkt(invalid_index, solution, kkt_tolerance)
            .input_valid);

    HighsLp missing_partition =
        makeOneRowLp(ObjSense::kMinimize, 0.0, kHighsInf, 1.0, kHighsInf);
    missing_partition.a_matrix_.format_ = MatrixFormat::kRowwisePartitioned;
    CHECK_FALSE(hipdlp::verification::verifyKkt(missing_partition, solution,
                                                kkt_tolerance)
                    .input_valid);
  }
}

TEST_CASE("hi-pdlp", "[pdlp][hipdlp]") {
  std::string model = "afiro";
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(model_file) != HighsStatus::kError);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = h.run();
  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const hipdlp::verification::KktVerificationResult verification =
      hipdlp::verification::verifyKkt(h.getLp(), h.getSolution(),
                                      kkt_tolerance);
  INFO(verification.error);
  CAPTURE(verification.metrics.relative_primal_residual,
          verification.metrics.relative_dual_residual,
          verification.metrics.relative_objective_gap,
          verification.metrics.max_primal_violation,
          verification.metrics.max_dual_violation);
  CHECK(verification.input_valid);
  CHECK(verification.metrics_finite);
  CHECK(verification.pass);

  const HighsInfo core_kkt_info = getCoreKktInfo(h);
  CHECK(core_kkt_info.max_relative_primal_infeasibility <= kkt_tolerance);
  CHECK(core_kkt_info.max_relative_dual_infeasibility <= kkt_tolerance);
  CHECK(core_kkt_info.primal_dual_objective_error <= kkt_tolerance);
  CHECK(coreKktPass(core_kkt_info, kkt_tolerance) == verification.pass);
#ifndef CUPDLP_GPU
  CHECK(h.getInfo().pdlp_iteration_count == 320);
#endif

  const bool cupdlp_test = true;
  if (cupdlp_test) {
    h.clearSolver();
    h.setOptionValue("solver", kPdlpString);
    run_status = h.run();
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  }
  h.resetGlobalScheduler(true);
}

#ifndef CUPDLP_GPU
TEST_CASE("hi-pdlp-legacy-optimal-kkt-characterization",
          "[pdlp][hipdlp][legacy]") {
  struct LegacyKktCase {
    const char* model;
    bool strict_kkt_pass;
  };
  const std::vector<LegacyKktCase> cases = {
      {"25fv47", false},  {"adlittle", false}, {"afiro", true},
      {"avgas", false},   {"blending", true},  {"chip", true},
      {"e226", true},     {"shell", true},     {"standata", false},
      {"standgub", false}};

  for (const LegacyKktCase& test_case : cases) {
    CAPTURE(test_case.model, test_case.strict_kkt_pass);
    Highs h;
    configureLegacyHiPdlp(h, test_case.model, 500000);
    REQUIRE(h.run() != HighsStatus::kError);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

    const auto verification = hipdlp::verification::verifyKkt(
        h.getLp(), h.getSolution(), kkt_tolerance);
    CAPTURE(verification.error, verification.metrics.relative_primal_residual,
            verification.metrics.relative_dual_residual,
            verification.metrics.relative_objective_gap,
            h.getInfo().pdlp_iteration_count);
    REQUIRE(verification.input_valid);
    REQUIRE(verification.metrics_finite);
    CHECK(verification.pass == test_case.strict_kkt_pass);

    const HighsInfo core_kkt_info = getCoreKktInfo(h);
    CHECK(coreKktPass(core_kkt_info, kkt_tolerance) ==
          test_case.strict_kkt_pass);
    CHECK(coreKktPass(core_kkt_info, kkt_tolerance) == verification.pass);
    h.resetGlobalScheduler(true);
  }
}

TEST_CASE("hi-pdlp-scrs8-legacy-characterization", "[pdlp][hipdlp][legacy]") {
  const HighsInt characterization_iteration_limit = 40000;

  SECTION("default adaptive option is forced to divergent PID") {
    Highs h;
    configureLegacyHiPdlp(h, "scrs8", characterization_iteration_limit);
    HighsInt configured_strategy = -1;
    REQUIRE(h.getOptionValue("pdlp_step_size_strategy", configured_strategy) ==
            HighsStatus::kOk);
    REQUIRE(configured_strategy == kPdlpStepSizeStrategyAdaptive);

    REQUIRE(h.run() == HighsStatus::kWarning);
    CHECK(h.getModelStatus() == HighsModelStatus::kIterationLimit);
    CHECK(h.getInfo().pdlp_iteration_count == characterization_iteration_limit);

    const auto verification = hipdlp::verification::verifyKkt(
        h.getLp(), h.getSolution(), kkt_tolerance);
    CAPTURE(verification.error, verification.metrics.relative_primal_residual,
            verification.metrics.relative_dual_residual,
            verification.metrics.relative_objective_gap);
    REQUIRE(verification.input_valid);
    REQUIRE(verification.metrics_finite);
    CHECK_FALSE(verification.pass);

    const HighsInfo core_kkt_info = getCoreKktInfo(h);
    CHECK_FALSE(coreKktPass(core_kkt_info, kkt_tolerance));
    CHECK(coreKktPass(core_kkt_info, kkt_tolerance) == verification.pass);

    Highs explicit_pid;
    configureLegacyHiPdlp(explicit_pid, "scrs8",
                          characterization_iteration_limit);
    REQUIRE(explicit_pid.setOptionValue("pdlp_step_size_strategy",
                                        kPdlpStepSizeStrategyPid) ==
            HighsStatus::kOk);
    REQUIRE(explicit_pid.run() == HighsStatus::kWarning);
    CHECK(explicit_pid.getModelStatus() == HighsModelStatus::kIterationLimit);
    CHECK(explicit_pid.getInfo().pdlp_iteration_count ==
          characterization_iteration_limit);
    CHECK(explicit_pid.getSolution().col_value == h.getSolution().col_value);
    CHECK(explicit_pid.getSolution().row_dual == h.getSolution().row_dual);
    explicit_pid.resetGlobalScheduler(true);
    h.resetGlobalScheduler(true);
  }

  SECTION("fixed step reports the frozen false optimum") {
    Highs h;
    configureLegacyHiPdlp(h, "scrs8", characterization_iteration_limit);
    REQUIRE(h.setOptionValue("pdlp_step_size_strategy",
                             kPdlpStepSizeStrategyFixed) == HighsStatus::kOk);

    REQUIRE(h.run() == HighsStatus::kOk);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
    CHECK(h.getInfo().pdlp_iteration_count == 21040);
    CHECK(h.getInfo().objective_function_value ==
          Approx(907.71412156).epsilon(1e-8));

    const auto verification = hipdlp::verification::verifyKkt(
        h.getLp(), h.getSolution(), kkt_tolerance);
    CAPTURE(verification.error, verification.metrics.relative_primal_residual,
            verification.metrics.relative_dual_residual,
            verification.metrics.relative_objective_gap);
    REQUIRE(verification.input_valid);
    REQUIRE(verification.metrics_finite);
    CHECK_FALSE(verification.pass);
    CHECK_FALSE(verification.primal_feasible);
    CHECK(verification.metrics.relative_primal_residual > kkt_tolerance);

    const HighsInfo core_kkt_info = getCoreKktInfo(h);
    CHECK(core_kkt_info.max_relative_primal_infeasibility > kkt_tolerance);
    CHECK_FALSE(coreKktPass(core_kkt_info, kkt_tolerance));
    CHECK(coreKktPass(core_kkt_info, kkt_tolerance) == verification.pass);
    h.resetGlobalScheduler(true);
  }
}
#endif
