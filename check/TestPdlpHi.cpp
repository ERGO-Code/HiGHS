#include <iostream>
#include <cmath>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "ExperimentalPslpHiPdlp.h"
#include "Highs.h"
#include "catch.hpp"

#ifdef CUPDLP_GPU
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#endif

const bool dev_run = false;
const double kkt_tolerance = 1e-4;

TEST_CASE("hi-pdlp", "[pdlp]") {
  std::string model = "afiro";
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs h;
  REQUIRE(h.setOptionValue("output_flag", dev_run) == HighsStatus::kOk);
  REQUIRE(h.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(h.setOptionValue("solver", kHiPdlpString) == HighsStatus::kOk);
  REQUIRE(h.setOptionValue("kkt_tolerance", kkt_tolerance) ==
          HighsStatus::kOk);
  HighsStatus run_status = h.run();
  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const bool cupdlp_test = true;
  if (cupdlp_test) {
    REQUIRE(h.clearSolver() == HighsStatus::kOk);
    REQUIRE(h.setOptionValue("solver", kPdlpString) == HighsStatus::kOk);
    run_status = h.run();
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  }
  h.resetGlobalScheduler(true);
}

TEST_CASE("hi-pdlp-with-pslp-presolve", "[pdlp][pslp]") {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/afiro.mps";

  ExperimentalHiPdlpPslpRun none_result;
  ExperimentalHiPdlpPslpRun highs_result;
  ExperimentalHiPdlpPslpRun pslp_result;
  std::string error_message;

  REQUIRE(runExperimentalHiPdlpPslpBenchmark(model_file, "none", none_result,
                                             error_message) ==
          HighsStatus::kOk);
  REQUIRE(error_message.empty());

  error_message.clear();
  REQUIRE(runExperimentalHiPdlpPslpBenchmark(model_file, "highs", highs_result,
                                             error_message) ==
          HighsStatus::kOk);
  REQUIRE(error_message.empty());

  error_message.clear();
  REQUIRE(runExperimentalHiPdlpPslpBenchmark(model_file, "pslp", pslp_result,
                                             error_message) ==
          HighsStatus::kOk);
  REQUIRE(error_message.empty());

  std::cout << experimentalHiPdlpPslpCsvLine(none_result) << std::endl;
  std::cout << experimentalHiPdlpPslpCsvLine(highs_result) << std::endl;
  std::cout << experimentalHiPdlpPslpCsvLine(pslp_result) << std::endl;
  std::cout << "diagnostics,none_iterations=" << none_result.iterations
            << ",highs_iterations=" << highs_result.iterations
            << ",pslp_iterations=" << pslp_result.iterations
            << ",none_total_time=" << none_result.total_time
            << ",highs_total_time=" << highs_result.total_time
            << ",pslp_total_time=" << pslp_result.total_time
            << ",highs_presolve_time=" << highs_result.presolve_time
            << ",pslp_presolve_time=" << pslp_result.presolve_time
            << ",objective_delta_pslp_vs_none="
            << (pslp_result.objective - none_result.objective)
            << ",nnz_reduction="
            << (none_result.orig_nnz - pslp_result.reduced_nnz) << std::endl;

  REQUIRE(none_result.model_status == HighsModelStatus::kOptimal);
  REQUIRE(highs_result.model_status == HighsModelStatus::kOptimal);
  REQUIRE(pslp_result.model_status == HighsModelStatus::kOptimal);
  REQUIRE(std::isfinite(none_result.objective));
  REQUIRE(std::isfinite(highs_result.objective));
  REQUIRE(std::isfinite(pslp_result.objective));
  REQUIRE(none_result.iterations >= 0);
  REQUIRE(highs_result.iterations >= 0);
  REQUIRE(pslp_result.iterations >= 0);
  REQUIRE(highs_result.reduced_rows >= 0);
  REQUIRE(highs_result.reduced_cols >= 0);
  REQUIRE(highs_result.reduced_nnz >= 0);
  REQUIRE(pslp_result.reduced_rows >= 0);
  REQUIRE(pslp_result.reduced_cols >= 0);
  REQUIRE(pslp_result.reduced_nnz >= 0);
}
