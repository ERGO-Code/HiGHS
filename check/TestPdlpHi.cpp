#include <chrono>

#include "HCheckConfig.h"
#include "HConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
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
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(model_file) != HighsStatus::kError);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = h.run();
  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

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
