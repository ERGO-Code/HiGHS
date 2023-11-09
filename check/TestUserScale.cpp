#include <cmath>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;
const double inf = kHighsInf;

TEST_CASE("user-cost-scale-after-run", "[highs_user_scale]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  const HighsLp& lp = highs.getLp();
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(filename);
  highs.run();
  HighsSolution unscaled_solution = solution;
  HighsInfo unscaled_info = info;
  HighsLp unscaled_lp = lp;
  double max_primal_infeasibility = info.max_primal_infeasibility;
  double max_dual_infeasibility = info.max_dual_infeasibility;
  double sum_dual_infeasibilities = info.sum_dual_infeasibilities;
  printf("Max primal infeasibility = %g\n", max_primal_infeasibility);
  printf("Max dual infeasibility = %g\n", max_dual_infeasibility);
  printf("Sum dual infeasibility = %g\n", sum_dual_infeasibilities);
  double objective_function_value = info.objective_function_value;

  HighsInt user_bound_scale = 10;
  double user_bound_scale_value = std::pow(2, user_bound_scale);
  highs.setOptionValue("user_bound_scale", user_bound_scale);

  HighsInt user_cost_scale = 30;
  double user_cost_scale_value = std::pow(2, user_cost_scale);
  highs.setOptionValue("user_cost_scale", user_cost_scale);

  REQUIRE(highs.getModelStatus() == HighsModelStatus::kNotset);
  REQUIRE(info.dual_solution_status == kSolutionStatusInfeasible);
  REQUIRE(info.objective_function_value ==
          user_cost_scale_value * user_bound_scale_value * objective_function_value);
  REQUIRE(info.num_dual_infeasibilities == kHighsIllegalInfeasibilityCount);
  REQUIRE(info.max_dual_infeasibility ==
          user_cost_scale_value * max_dual_infeasibility);
  REQUIRE(info.sum_dual_infeasibilities ==
          user_cost_scale_value * sum_dual_infeasibilities);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    REQUIRE(solution.col_value[iCol] == unscaled_solution.col_value[iCol] * user_bound_scale_value);
    REQUIRE(solution.col_dual[iCol] == unscaled_solution.col_dual[iCol] * user_cost_scale_value);
    REQUIRE(lp.col_cost_[iCol] == unscaled_lp.col_cost_[iCol] * user_cost_scale_value);
    if (lp.col_lower_[iCol] > -inf) 
      REQUIRE(lp.col_lower_[iCol] == unscaled_lp.col_lower_[iCol] * user_bound_scale_value);
    if (lp.col_upper_[iCol] < inf) 
      REQUIRE(lp.col_upper_[iCol] == unscaled_lp.col_upper_[iCol] * user_bound_scale_value);
  }
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    REQUIRE(solution.row_value[iRow] == unscaled_solution.row_value[iRow] * user_bound_scale_value);
    REQUIRE(solution.row_dual[iRow] == unscaled_solution.row_dual[iRow] * user_cost_scale_value);
    if (lp.row_lower_[iRow] > -inf) 
      REQUIRE(lp.row_lower_[iRow] == unscaled_lp.row_lower_[iRow] * user_bound_scale_value);
    if (lp.row_upper_[iRow] < inf) 
      REQUIRE(lp.row_upper_[iRow] == unscaled_lp.row_upper_[iRow] * user_bound_scale_value);
  }
}

TEST_CASE("user-cost-scale-before-run", "[highs_user_scale]") {
    Highs highs;
  const HighsLp& lp = highs.getLp();
  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();
  const HighsInt user_cost_scale = -30;
  const double user_cost_scale_value = std::pow(2, user_cost_scale);
  const double unscaled_col0_cost = 1e14;
  highs.addVar(0, 1);
  highs.changeColCost(0, unscaled_col0_cost);

  highs.setOptionValue("user_cost_scale", user_cost_scale);
  REQUIRE(lp.col_cost_[0] == unscaled_col0_cost *  user_cost_scale_value);

  const double unscaled_col1_cost = 1e12;
  highs.addVar(0, 1);
  highs.changeColCost(1, unscaled_col1_cost);
  REQUIRE(lp.col_cost_[1] == unscaled_col1_cost *  user_cost_scale_value);
}
