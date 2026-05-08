#include <iostream>

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

  ExperimentalHiPdlpPslpRun plain_result;
  ExperimentalHiPdlpPslpRun pslp_result;
  std::string error_message;

  REQUIRE(runExperimentalHiPdlpPslpBenchmark(model_file, "plain", plain_result,
                                             error_message) ==
          HighsStatus::kOk);
  REQUIRE(error_message.empty());

  error_message.clear();
  REQUIRE(runExperimentalHiPdlpPslpBenchmark(model_file, "pslp", pslp_result,
                                             error_message) ==
          HighsStatus::kOk);
  REQUIRE(error_message.empty());

  std::cout << experimentalHiPdlpPslpCsvLine(plain_result) << std::endl;
  std::cout << experimentalHiPdlpPslpCsvLine(pslp_result) << std::endl;
  std::cout << "diagnostics,plain_iterations=" << plain_result.iterations
            << ",pslp_iterations=" << pslp_result.iterations
            << ",plain_total_time=" << plain_result.total_time
            << ",pslp_total_time=" << pslp_result.total_time
            << ",nnz_reduction="
            << (plain_result.orig_nnz - pslp_result.reduced_nnz) << std::endl;

  REQUIRE(plain_result.model_status == HighsModelStatus::kOptimal);
  REQUIRE(pslp_result.model_status == HighsModelStatus::kOptimal);
  REQUIRE(experimentalObjectivesClose(plain_result.objective,
                                      pslp_result.objective));
  REQUIRE(plain_result.iterations >= 0);
  REQUIRE(pslp_result.iterations >= 0);
  REQUIRE(pslp_result.reduced_rows >= 0);
  REQUIRE(pslp_result.reduced_cols >= 0);
  REQUIRE(pslp_result.reduced_nnz >= 0);
}
