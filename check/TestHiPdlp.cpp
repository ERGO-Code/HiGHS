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

TEST_CASE("hi-pdlp", "[hipdlp]") {
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

TEST_CASE("hi-pdlp-maximize", "[hipdlp]") {
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

TEST_CASE("hi-pdlp-free-rows", "[hipdlp]") {
  // min x0 + 2 x1 + 3 x2
  // row 0:  x0 + x1       >= 1
  // row 1:  free (x0 - x2)
  // row 2:  x1 + x2       >= 2
  // row 3:  x0 + x1 + x2  =  4
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 4;
  lp.col_cost_ = {1.0, 2.0, 3.0};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {kHighsInf, kHighsInf, kHighsInf};
  lp.row_lower_ = {1.0, -kHighsInf, 2.0, 4.0};
  lp.row_upper_ = {kHighsInf, kHighsInf, kHighsInf, 4.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.start_ = {0, 3, 6, 9};
  lp.a_matrix_.index_ = {0, 1, 3, 0, 2, 3, 1, 2, 3};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 1, 1, -1, 1, 1};

  for (const ObjSense sense : {ObjSense::kMinimize, ObjSense::kMaximize}) {
    lp.sense_ = sense;

    Highs ref;
    ref.setOptionValue("output_flag", dev_run);
    REQUIRE(ref.passModel(lp) == HighsStatus::kOk);
    ref.setOptionValue("solver", kSimplexString);
    ref.setOptionValue("presolve", kHighsOffString);
    REQUIRE(ref.run() == HighsStatus::kOk);
    const double ref_objective = ref.getInfo().objective_function_value;

    Highs h;
    h.setOptionValue("output_flag", dev_run);
    REQUIRE(h.passModel(lp) == HighsStatus::kOk);
    h.setOptionValue("solver", kHiPdlpString);
    h.setOptionValue("presolve", kHighsOffString);
    h.setOptionValue("kkt_tolerance", kkt_tolerance);
    REQUIRE(h.run() == HighsStatus::kOk);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
    REQUIRE(std::fabs(h.getInfo().objective_function_value - ref_objective) <
            1e-3 * (1.0 + std::fabs(ref_objective)));

    const HighsSolution& solution = h.getSolution();
    REQUIRE(solution.row_dual.size() == 4);
    REQUIRE(solution.row_value.size() == 4);
    REQUIRE(solution.row_dual[1] == 0.0);  // free row
    h.resetGlobalScheduler(true);
  }
}

TEST_CASE("hi-pdlp-all-free-rows", "[hipdlp]") {
  // min x0 + 2 x1
  // 0 <= x0 <= 3, 1 <= x1 <= 2
  // two free rows, so there are no constraints left after they are removed:
  // the optimum is x = (0, 1), objective 2
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {1.0, 2.0};
  lp.col_lower_ = {0.0, 1.0};
  lp.col_upper_ = {3.0, 2.0};
  lp.row_lower_ = {-kHighsInf, -kHighsInf};
  lp.row_upper_ = {kHighsInf, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1.0, 1.0, 1.0, -1.0};

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("presolve", kHighsOffString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  REQUIRE(h.run() == HighsStatus::kOk);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(std::fabs(h.getInfo().objective_function_value - 2.0) < 1e-3);

  const HighsSolution& solution = h.getSolution();
  REQUIRE(solution.row_dual.size() == 2);
  REQUIRE(solution.row_dual[0] == 0.0);
  REQUIRE(solution.row_dual[1] == 0.0);

  h.resetGlobalScheduler(true);
}
