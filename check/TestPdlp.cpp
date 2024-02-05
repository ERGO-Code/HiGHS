#include "HCheckConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = true;
const double double_equal_tolerance = 1e-3;

TEST_CASE("pdlp-distillation-lp", "[pdlp]") {
  SpecialLps special_lps;
  HighsLp lp;

  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(lp, require_model_status, optimal_objective);

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();
  const HighsOptions& options = highs.getOptions();
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("primal_feasibility_tolerance", 1e-4);
  highs.setOptionValue("dual_feasibility_tolerance", 1e-4);
  REQUIRE(highs.run() == HighsStatus::kOk);
  highs.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
}

TEST_CASE("pdlp-3d-lp", "[pdlp]") {
  SpecialLps special_lps;
  HighsLp lp;

  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.ThreeDLp(lp, require_model_status, optimal_objective);

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();
  const HighsOptions& options = highs.getOptions();
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("primal_feasibility_tolerance", 1e-4);
  highs.setOptionValue("dual_feasibility_tolerance", 1e-4);
  REQUIRE(highs.run() == HighsStatus::kOk);
  highs.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
}

TEST_CASE("pdlp-boxed-row-lp", "[pdlp]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {-1, -2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, 6};
  lp.row_lower_ = {3, -4};
  lp.row_upper_ = {10, 2};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 1, -1};
  double optimal_objective = -16;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  REQUIRE(highs.run() == HighsStatus::kOk);
  highs.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
}

TEST_CASE("pdlp-infeasible-lp", "[pdlp]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {-1, -2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {-1};
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  REQUIRE(highs.run() == HighsStatus::kOk);
  highs.writeSolution("", 1);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnboundedOrInfeasible);
}

TEST_CASE("pdlp-unbounded-lp", "[pdlp]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {-1, -2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {1};
  lp.row_upper_ = {inf};
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  REQUIRE(highs.run() == HighsStatus::kOk);
  highs.writeSolution("", 1);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnboundedOrInfeasible);
}
