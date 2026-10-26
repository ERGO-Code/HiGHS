#include <chrono>
#include <cmath>

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
  h.resetThreadScheduler(true);
}

TEST_CASE("hi-pdlp-maximize", "[pdlp]") {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/3459.mps";

  // Reference: minimize with simplex
  Highs ref;
  ref.setOptionValue("output_flag", dev_run);
  REQUIRE(ref.readModel(model_file) != HighsStatus::kError);
  ref.setOptionValue("solver", kSimplexString);
  ref.setOptionValue("presolve", kHighsOffString);
  REQUIRE(ref.run() == HighsStatus::kOk);
  const double min_objective = ref.getInfo().objective_function_value;

  // Negate all costs and maximize: the optimal objective is -min_objective
  HighsLp lp = ref.getLp();
  for (double& cost : lp.col_cost_) cost = -cost;
  lp.sense_ = ObjSense::kMaximize;

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("presolve", kHighsOffString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  REQUIRE(h.run() == HighsStatus::kOk);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(std::fabs(h.getInfo().objective_function_value + min_objective) <
          1e-3 * (1.0 + std::fabs(min_objective)));

  // The objective offset must also be reported on the original (maximize)
  // scale
  const double offset = 5.0;
  lp.offset_ = offset;
  h.clearSolver();
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  REQUIRE(h.run() == HighsStatus::kOk);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(std::fabs(h.getInfo().objective_function_value -
                    (-min_objective + offset)) <
          1e-3 * (1.0 + std::fabs(min_objective)));

  h.resetGlobalScheduler(true);
}
