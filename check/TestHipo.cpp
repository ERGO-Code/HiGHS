#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/Filereader.h"
#include "ipm/hipo/ipm/Solver.h"
#include "ipm/hipo/ipm/Status.h"
#include "lp_data/HighsCallback.h"
#include "parallel/HighsParallel.h"

// Example for using HiPO from its C++ interface. The program solves the Netlib
// problem afiro.

// #include <unistd.h>

#include <cmath>
#include <iostream>
#include <vector>

const bool dev_run = false;

TEST_CASE("test-hipo-afiro", "[highs_hipo]") {
  // Test that hipo runs and finds correct solution for afiro

  std::string model = "afiro.mps";
  const double expected_obj = -464.753;

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("solver", kHipoString);
  highs.setOptionValue("timeless_log", kHighsOnString);

  std::string filename = std::string(HIGHS_DIR) + "/check/instances/" + model;
  highs.readModel(filename);

  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  const double actual_obj = highs.getObjectiveValue();
  REQUIRE(std::abs(actual_obj - expected_obj) < 0.001);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-hipo-deterministic", "[highs_hipo]") {
  // Test that hipo finds the exact same solution if run twice

  std::string model = "80bau3b.mps";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue(kSolverString, kHipoString);
  highs.setOptionValue(kParallelString, kHighsOnString);
  highs.setOptionValue(kRunCrossoverString, kHighsOffString);

  std::string filename = std::string(HIGHS_DIR) + "/check/instances/" + model;
  highs.readModel(filename);

  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  const HighsSolution solution_1 = highs.getSolution();

  highs.run();
  const HighsSolution solution_2 = highs.getSolution();

  REQUIRE(solution_1.value_valid == solution_2.value_valid);
  REQUIRE(solution_1.dual_valid == solution_2.dual_valid);
  REQUIRE(solution_1.col_value == solution_2.col_value);
  REQUIRE(solution_1.row_value == solution_2.row_value);
  REQUIRE(solution_1.col_dual == solution_2.col_dual);
  REQUIRE(solution_1.row_dual == solution_2.row_dual);

  highs.resetGlobalScheduler(true);
}