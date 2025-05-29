#include <cmath>

#include "Highs.h"
#include "catch.hpp"

#include "HCheckConfig.h"

const bool dev_run = false;
const double inf = kHighsInf;

void checkModelScaling(const HighsInt user_bound_scale,
                       const HighsInt user_cost_scale,
                       const HighsModel& unscaled_model,
                       const HighsModel& scaled_model);

void checkLpScaling(const HighsInt user_bound_scale,
                    const HighsInt user_cost_scale, const HighsLp& unscaled_lp,
                    const HighsLp& scaled_lp);

void checkSolutionScaling(const HighsInt user_bound_scale,
                          const HighsInt user_cost_scale,
                          const HighsSolution& unscaled_solution,
                          const HighsSolution& scaled_solution);

TEST_CASE("user-cost-scale-after-run", "[highs_user_scale]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);
  highs.run();
  HighsInfo unscaled_info = info;
  HighsSolution unscaled_solution = highs.getSolution();
  HighsLp unscaled_lp = highs.getLp();
  double max_primal_infeasibility = info.max_primal_infeasibility;
  double max_dual_infeasibility = info.max_dual_infeasibility;
  double sum_dual_infeasibilities = info.sum_dual_infeasibilities;
  double objective_function_value = info.objective_function_value;

  HighsInt user_bound_scale = 10;
  double user_bound_scale_value = std::pow(2, user_bound_scale);
  highs.setOptionValue("user_bound_scale", user_bound_scale);

  HighsInt user_cost_scale = 30;
  double user_cost_scale_value = std::pow(2, user_cost_scale);
  highs.setOptionValue("user_cost_scale", user_cost_scale);

  HighsLp scaled_lp = highs.getLp();
  HighsSolution scaled_solution = highs.getSolution();
  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);
  checkSolutionScaling(user_bound_scale, user_cost_scale, unscaled_solution,
                       scaled_solution);

  REQUIRE(highs.getModelStatus() == HighsModelStatus::kNotset);
  REQUIRE(info.dual_solution_status == kSolutionStatusInfeasible);
  REQUIRE(info.objective_function_value == user_cost_scale_value *
                                               user_bound_scale_value *
                                               objective_function_value);
  REQUIRE(info.num_dual_infeasibilities == kHighsIllegalInfeasibilityCount);
  REQUIRE(info.max_dual_infeasibility ==
          user_cost_scale_value * max_dual_infeasibility);
  REQUIRE(info.sum_dual_infeasibilities ==
          user_cost_scale_value * sum_dual_infeasibilities);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("user-cost-scale-after-load", "[highs_user_scale]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.readModel(filename);
  HighsLp unscaled_lp = highs.getLp();

  HighsInt user_bound_scale = 10;
  double user_bound_scale_value = std::pow(2, user_bound_scale);
  highs.setOptionValue("user_bound_scale", user_bound_scale);

  HighsInt user_cost_scale = 30;
  double user_cost_scale_value = std::pow(2, user_cost_scale);
  highs.setOptionValue("user_cost_scale", user_cost_scale);

  highs.readModel(filename);
  HighsLp scaled_lp = highs.getLp();

  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);
  //  checkSolutionScaling(user_bound_scale, user_cost_scale, unscaled_solution,
  //  scaled_solution);
  highs.run();

  highs.resetGlobalScheduler(true);
}

TEST_CASE("user-small-cost-scale", "[highs_user_scale]") {
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {10, 25};
  lp.sense_ = ObjSense::kMaximize;
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf, -inf};
  lp.row_upper_ = {80, 120};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 2, 4};
  highs.passModel(lp);
  highs.run();
  REQUIRE(solution.col_value[0] == 40);
  REQUIRE(solution.col_value[1] == 20);

  highs.setOptionValue("user_cost_scale", -30);
  highs.clearSolver();
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(solution.col_value[0] == 0);
  REQUIRE(solution.col_value[1] == 0);

  highs.setOptionValue("user_cost_scale", 0);

  highs.run();
  REQUIRE(solution.col_value[0] == 40);
  REQUIRE(solution.col_value[1] == 20);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("user-cost-scale-in-build", "[highs_user_scale]") {
  Highs unscaled_highs;
  Highs scaled_highs;
  unscaled_highs.setOptionValue("output_flag", dev_run);
  scaled_highs.setOptionValue("output_flag", dev_run);
  const HighsLp& unscaled_lp = unscaled_highs.getLp();
  const HighsLp& scaled_lp = scaled_highs.getLp();
  const HighsInfo& info = scaled_highs.getInfo();
  const HighsSolution& solution = scaled_highs.getSolution();
  const HighsInt user_cost_scale = -30;
  const HighsInt user_bound_scale = 10;
  const double unscaled_col0_cost = 1e14;
  unscaled_highs.addVar(0, inf);
  scaled_highs.addVar(0, inf);
  unscaled_highs.changeColCost(0, unscaled_col0_cost);
  scaled_highs.changeColCost(0, unscaled_col0_cost);

  scaled_highs.setOptionValue("user_cost_scale", user_cost_scale);
  scaled_highs.setOptionValue("user_bound_scale", user_bound_scale);
  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);

  const double unscaled_col1_cost = 1e12;
  unscaled_highs.addVar(1, inf);
  scaled_highs.addVar(1, inf);
  unscaled_highs.changeColCost(1, unscaled_col1_cost);
  scaled_highs.changeColCost(1, unscaled_col1_cost);
  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);

  std::vector<HighsInt> index = {0, 1};
  std::vector<double> value0 = {1, 2};
  std::vector<double> value1 = {1, 4};
  unscaled_highs.addRow(-inf, 120, 2, index.data(), value0.data());
  scaled_highs.addRow(-inf, 120, 2, index.data(), value0.data());
  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);

  unscaled_highs.addRow(-inf, 150, 2, index.data(), value1.data());
  scaled_highs.addRow(-inf, 150, 2, index.data(), value1.data());
  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);

  std::vector<double> cost = {0, 10};
  std::vector<double> lower = {2, 4};
  std::vector<double> upper = {inf, inf};
  std::vector<HighsInt> matrix_start = {0, 2};
  std::vector<HighsInt> matrix_index = {0, 1, 0, 1};
  std::vector<double> matrix_value = {1, 1, 2, 4};
  unscaled_highs.addCols(2, cost.data(), lower.data(), upper.data(), 4,
                         matrix_start.data(), matrix_index.data(),
                         matrix_value.data());
  scaled_highs.addCols(2, cost.data(), lower.data(), upper.data(), 4,
                       matrix_start.data(), matrix_index.data(),
                       matrix_value.data());
  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);

  lower = {-inf, 0};
  upper = {120, 150};
  matrix_start = {0, 2};
  matrix_index = {0, 2, 1, 3};
  matrix_value = {1, 1, 2, 4};
  unscaled_highs.addRows(2, lower.data(), upper.data(), 4, matrix_start.data(),
                         matrix_index.data(), matrix_value.data());
  scaled_highs.addRows(2, lower.data(), upper.data(), 4, matrix_start.data(),
                       matrix_index.data(), matrix_value.data());

  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);
}

void checkModelScaling(const HighsInt user_bound_scale,
                       const HighsInt user_cost_scale,
                       const HighsModel& unscaled_model,
                       const HighsModel& scaled_model) {
  checkLpScaling(user_bound_scale, user_cost_scale, unscaled_model.lp_,
                 scaled_model.lp_);
}

void checkLpScaling(const HighsInt user_bound_scale,
                    const HighsInt user_cost_scale, const HighsLp& unscaled_lp,
                    const HighsLp& scaled_lp) {
  const double user_bound_scale_value = std::pow(2, user_bound_scale);
  const double user_cost_scale_value = std::pow(2, user_cost_scale);
  REQUIRE(unscaled_lp.num_col_ == scaled_lp.num_col_);
  REQUIRE(unscaled_lp.num_row_ == scaled_lp.num_row_);
  for (HighsInt iCol = 0; iCol < unscaled_lp.num_col_; iCol++) {
    REQUIRE(scaled_lp.col_cost_[iCol] ==
            unscaled_lp.col_cost_[iCol] * user_cost_scale_value);
    if (unscaled_lp.col_lower_[iCol] > -inf)
      REQUIRE(scaled_lp.col_lower_[iCol] ==
              unscaled_lp.col_lower_[iCol] * user_bound_scale_value);
    if (unscaled_lp.col_upper_[iCol] < inf)
      REQUIRE(scaled_lp.col_upper_[iCol] ==
              unscaled_lp.col_upper_[iCol] * user_bound_scale_value);
  }
  for (HighsInt iRow = 0; iRow < unscaled_lp.num_row_; iRow++) {
    if (unscaled_lp.row_lower_[iRow] > -inf)
      REQUIRE(scaled_lp.row_lower_[iRow] ==
              unscaled_lp.row_lower_[iRow] * user_bound_scale_value);
    if (unscaled_lp.row_upper_[iRow] < inf)
      REQUIRE(scaled_lp.row_upper_[iRow] ==
              unscaled_lp.row_upper_[iRow] * user_bound_scale_value);
  }
}

void checkSolutionScaling(const HighsInt user_bound_scale,
                          const HighsInt user_cost_scale,
                          const HighsSolution& unscaled_solution,
                          const HighsSolution& scaled_solution) {
  const double user_bound_scale_value = std::pow(2, user_bound_scale);
  const double user_cost_scale_value = std::pow(2, user_cost_scale);
  for (HighsInt iCol = 0; iCol < HighsInt(unscaled_solution.col_value.size());
       iCol++) {
    REQUIRE(scaled_solution.col_value[iCol] ==
            unscaled_solution.col_value[iCol] * user_bound_scale_value);
    REQUIRE(scaled_solution.col_dual[iCol] ==
            unscaled_solution.col_dual[iCol] * user_cost_scale_value);
  }
  for (HighsInt iRow = 0; iRow < HighsInt(unscaled_solution.row_value.size());
       iRow++) {
    REQUIRE(scaled_solution.row_value[iRow] ==
            unscaled_solution.row_value[iRow] * user_bound_scale_value);
    REQUIRE(scaled_solution.row_dual[iRow] ==
            unscaled_solution.row_dual[iRow] * user_cost_scale_value);
  }
}
