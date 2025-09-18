#include <cmath>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;  // false;
const double inf = kHighsInf;

void checkModelScaling(const HighsInt user_bound_scale,
                       const HighsInt user_cost_scale,
                       const HighsModel& unscaled_model,
                       const HighsModel& scaled_model);

void checkLpScaling(const HighsInt user_bound_scale,
                    const HighsInt user_cost_scale, const HighsLp& unscaled_lp,
                    const HighsLp& scaled_lp);

TEST_CASE("user-cost-scale-after-run", "[highs_user_scale]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);
  for (HighsInt k = 0; k < 2; k++) {
    highs.readModel(filename);
    highs.run();
    HighsInfo unscaled_info = info;
    HighsSolution unscaled_solution = highs.getSolution();
    HighsLp unscaled_lp = highs.getLp();

    HighsInt user_bound_scale = 10;
    REQUIRE(highs.setOptionValue("user_bound_scale", user_bound_scale) == HighsStatus::kOk);

    HighsInt user_cost_scale = 30;
    REQUIRE(highs.setOptionValue("user_cost_scale", user_cost_scale) == HighsStatus::kOk);

    HighsLp scaled_lp = highs.getLp();
    checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);

    REQUIRE(highs.getModelStatus() == HighsModelStatus::kNotset);
    REQUIRE(info.primal_solution_status == kSolutionStatusNone);
    REQUIRE(info.dual_solution_status == kSolutionStatusNone);
    REQUIRE(info.num_primal_infeasibilities == kHighsIllegalInfeasibilityCount);
    REQUIRE(info.max_primal_infeasibility == kHighsIllegalInfeasibilityMeasure);
    REQUIRE(info.sum_primal_infeasibilities == kHighsIllegalInfeasibilityMeasure);
    REQUIRE(info.num_dual_infeasibilities == kHighsIllegalInfeasibilityCount);
    REQUIRE(info.max_dual_infeasibility == kHighsIllegalInfeasibilityMeasure);
    REQUIRE(info.sum_dual_infeasibilities == kHighsIllegalInfeasibilityMeasure);

    filename =
      std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
    REQUIRE(highs.setOptionValue("user_bound_scale", 0) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("user_cost_scale", 0) == HighsStatus::kOk);
  }

  highs.resetGlobalScheduler(true);
}

TEST_CASE("user-cost-scale-after-load", "[highs_user_scale]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  for (HighsInt k = 0; k < 2; k++) {
    highs.readModel(filename);
    HighsLp unscaled_lp = highs.getLp();

    HighsInt user_bound_scale = 10;
    double user_bound_scale_value = std::pow(2, user_bound_scale);
    REQUIRE(highs.setOptionValue("user_bound_scale", user_bound_scale) == HighsStatus::kOk);

    HighsInt user_cost_scale = 30;
    double user_cost_scale_value = std::pow(2, user_cost_scale);
    REQUIRE(highs.setOptionValue("user_cost_scale", user_cost_scale) == HighsStatus::kOk);

    REQUIRE(highs.readModel(filename) == HighsStatus::kOk);
    HighsLp scaled_lp = highs.getLp();

    checkLpScaling(user_bound_scale, user_cost_scale, unscaled_lp, scaled_lp);

    filename =
      std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
    REQUIRE(highs.setOptionValue("user_bound_scale", 0) == HighsStatus::kOk);
    REQUIRE(highs.setOptionValue("user_cost_scale", 0) == HighsStatus::kOk);
  }

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

  REQUIRE(highs.setOptionValue("user_cost_scale", -30) == HighsStatus::kOk);
  highs.clearSolver();
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(solution.col_value[0] == 0);
  REQUIRE(solution.col_value[1] == 0);

  REQUIRE(highs.setOptionValue("user_cost_scale", 0) == HighsStatus::kOk);

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
  REQUIRE(scaled_lp.user_cost_scale_ == user_cost_scale);
  REQUIRE(scaled_lp.user_bound_scale_ == user_bound_scale);
  const bool has_integrality = scaled_lp.integrality_.size() > 0;
  const HighsSparseMatrix& unscaled_matrix = unscaled_lp.a_matrix_;
  REQUIRE(unscaled_matrix.isColwise());
  for (HighsInt iCol = 0; iCol < unscaled_lp.num_col_; iCol++) {
    bool continuous = !has_integrality || scaled_lp.integrality_[iCol] == HighsVarType::kContinuous;
    double value = unscaled_lp.col_cost_[iCol] * user_cost_scale_value;
    if (!continuous) value *= user_bound_scale_value;
    REQUIRE(scaled_lp.col_cost_[iCol] == value);
    if (unscaled_lp.col_lower_[iCol] > -inf) {
      value = unscaled_lp.col_lower_[iCol];
      if (continuous) value *= user_bound_scale_value;
      REQUIRE(scaled_lp.col_lower_[iCol] == value);
    }
    if (unscaled_lp.col_upper_[iCol] < inf) {
      value = unscaled_lp.col_upper_[iCol];
      if (continuous) value *= user_bound_scale_value;
      REQUIRE(scaled_lp.col_upper_[iCol] == value);
    }
    for (HighsInt iEl = unscaled_matrix.start_[iCol]; iEl < unscaled_matrix.start_[iCol+1]; iEl++) {
      value = unscaled_matrix.value_[iEl];
      if (!continuous) value *= user_bound_scale_value;
      REQUIRE(scaled_lp.a_matrix_.value_[iEl] == value);
    }
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

