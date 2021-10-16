#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const double inf = kHighsInf;
const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

HighsLp baseLp();

TEST_CASE("semi-variable-model", "[highs_test_semi_variables]") {
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  HighsStatus return_status;
  double optimal_objective_function_value;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp = baseLp();
  const HighsVarType continuous = HighsVarType::kContinuous;
  const HighsVarType semi_continuous = HighsVarType::kSemiContinuous;
  const HighsVarType semi_integer = HighsVarType::kSemiInteger;
  lp.integrality_ = {continuous, continuous, semi_continuous, continuous};
  const HighsInt semi_col = 2;
  const double save_semi_col_lower = lp.col_lower_[semi_col];
  const double save_semi_col_upper = lp.col_upper_[semi_col];
  lp.col_upper_[semi_col] = inf;
  return_status = highs.passModel(model);
  REQUIRE(return_status == HighsStatus::kError);
  lp.col_upper_[semi_col] = save_semi_col_upper;

  return_status = highs.passModel(model);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", true);
  optimal_objective_function_value = 6.83333;
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);
  // Remove the semi-condition and resolve
  highs.changeColIntegrality(semi_col, continuous);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", true);
  optimal_objective_function_value = 3.93333;
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);

  // Restore the semi-condition, change the cost and resolve
  highs.changeColIntegrality(semi_col, semi_continuous);
  highs.changeColCost(semi_col, -0.1);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", true);
  optimal_objective_function_value = 8.22333;
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);

  // Fix the variable at zero and resolve
  highs.changeColBounds(semi_col, 0, 0);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", true);
  optimal_objective_function_value = 6.83333;
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);

  // Change to sem-integer, restore the bounds and resolve
  highs.changeColIntegrality(semi_col, semi_integer);
  highs.changeColBounds(semi_col, save_semi_col_lower, save_semi_col_upper);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", true);
  optimal_objective_function_value = 8.13333;
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);
}

TEST_CASE("semi-variable-file", "[highs_test_semi_variables]") {
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  double optimal_objective_function_value;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  std::string model = "";
  std::string model_file;
  // Solve the same semi-continuous model from MPS and .lp files
  model = "semi-continuous";
  optimal_objective_function_value = 8.22333;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);
  // Solve the same semi-integer model from MPS and .lp files
  model = "semi-integer";
  optimal_objective_function_value = 8.13333;
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(fabs(info.objective_function_value -
               optimal_objective_function_value) < double_equal_tolerance);
}

HighsLp baseLp() {
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;
  lp.col_cost_ = {1, 2, -4, -3};
  lp.col_lower_ = {0, 0, 1.1, 0};
  lp.col_upper_ = {inf, inf, 10, inf};
  lp.row_lower_ = {-inf, 0, 0, 0.5};
  lp.row_upper_ = {5, inf, inf, inf};
  lp.a_start_ = {0, 3, 6, 7, 8};
  lp.a_index_ = {0, 1, 2, 0, 1, 2, 3, 3};
  lp.a_value_ = {1, 2, -1, 1, -1, 3, 1, 1};
  lp.format_ = MatrixFormat::kColwise;
  lp.sense_ = ObjSense::kMaximize;
  return lp;
}
