#include <cmath>
#include <iostream>
#include <vector>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/Filereader.h"
#include "ipm/hipo/ipm/Solver.h"
#include "ipm/hipo/ipm/Status.h"
#include "lp_data/HighsCallback.h"
#include "parallel/HighsParallel.h"

const bool dev_run = false;

void runHipoTest(Highs& highs, const std::string& model,
                 const double expected_obj) {
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("solver", kHipoString);
  highs.setOptionValue("timeless_log", kHighsOnString);

  std::string filename = std::string(HIGHS_DIR) + "/check/instances/" + model;
  highs.readModel(filename);

  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  const double actual_obj = highs.getObjectiveValue();
  REQUIRE(std::abs(actual_obj - expected_obj) / std::abs(expected_obj) < 1e-4);
}

TEST_CASE("test-hipo-afiro", "[highs_hipo]") {
  Highs highs;
  runHipoTest(highs, "afiro.mps", -464.753);
  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-deterministic", "[highs_hipo]") {
  // Test that hipo finds the exact same solution if run twice

  std::string model = "80bau3b.mps";
  const double expected_obj = 9.8722e5;

  HighsInt iter_1, iter_2;
  HighsSolution solution_1, solution_2;

  Highs highs;
  highs.setOptionValue(kRunCrossoverString, kHighsOffString);

  runHipoTest(highs, model, expected_obj);
  solution_1 = highs.getSolution();
  iter_1 = highs.getInfo().ipm_iteration_count;

  runHipoTest(highs, model, expected_obj);
  solution_2 = highs.getSolution();
  iter_2 = highs.getInfo().ipm_iteration_count;

  REQUIRE(iter_1 == iter_2);
  REQUIRE(solution_1.value_valid == solution_2.value_valid);
  REQUIRE(solution_1.dual_valid == solution_2.dual_valid);
  REQUIRE(solution_1.col_value == solution_2.col_value);
  REQUIRE(solution_1.row_value == solution_2.row_value);
  REQUIRE(solution_1.col_dual == solution_2.col_dual);
  REQUIRE(solution_1.row_dual == solution_2.row_dual);
}

TEST_CASE("test-hipo-orderings", "[highs_hipo]") {
  // Test that hipo orderings work correctly

  std::string model = "adlittle.mps";
  const double expected_obj = 2.2549e5;
  Highs highs;

  // metis
  highs.setOptionValue(kHipoOrderingString, kHipoMetisString);
  runHipoTest(highs, model, expected_obj);

  // amd
  highs.setOptionValue(kHipoOrderingString, kHipoAmdString);
  runHipoTest(highs, model, expected_obj);

  // rcm
  highs.setOptionValue(kHipoOrderingString, kHipoRcmString);
  runHipoTest(highs, model, expected_obj);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-qp", "[highs_hipo]") {
  Highs highs;
  runHipoTest(highs, "qptestnw.lp", -6.4500);
  runHipoTest(highs, "qjh.lp", -5.2500);
  runHipoTest(highs, "primal1.mps", -3.501296e-2);
  highs.resetGlobalScheduler(true);
}