#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;
const double inf = kHighsInf;

// No commas in test case name.
TEST_CASE("test-crossover", "[highs_crossover]") {
  HighsStatus return_status;
  std::string model;
  std::string model_file;
  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMaximize;
  lp.col_cost_ = {1, 2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {4};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1, 2};
  return_status = highs.passModel(lp);
  REQUIRE(return_status == HighsStatus::kOk);
  highs.run();
  const double require_optimal_objective = info.objective_function_value;
  highs.clearSolver();
  // Set solution to interior of optimal face
  HighsSolution solution;
  solution.col_value = {2, 1};
  return_status = highs.crossover(solution);
  REQUIRE(return_status == HighsStatus::kOk);
  highs.writeSolution("", 1);
  REQUIRE(require_optimal_objective == info.objective_function_value);
}
