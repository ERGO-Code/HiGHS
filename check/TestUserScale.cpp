#include <cmath>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;  // false;
const double inf = kHighsInf;

bool doubleEqual0(const double v0, const double v1) {
  double rel_difference = std::fabs(v0 - v1) / std::max(1.0, std::fabs(v0));
  bool ok = rel_difference < 1e-12;
  if (dev_run && !ok)
    printf("UserScaleDoubleEqual: %g and %g have relative difference = %g\n",
           v0, v1, rel_difference);
  return ok;
}

TEST_CASE("user-scale-after-run", "[highs_user_scale]") {
  const std::string mip_model = "flugpl";  //"rgn";//
  std::string model = "avgas";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);
  //    REQUIRE(highs.setOptionValue("presolve", kHighsOffString) ==
  //    HighsStatus::kOk);
  HighsInt num_k = 2;
  if (num_k == 1) model = mip_model;
  for (HighsInt k = 0; k < num_k; k++) {
    std::string filename =
        std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    highs.readModel(filename);
    HighsLp unscaled_lp = highs.getLp();
    const bool is_lp = !unscaled_lp.isMip();
    //    highs.writeModel("unscaled.mps");

    highs.run();
    double unscaled_objective = highs.getInfo().objective_function_value;

    HighsInt user_bound_scale = 1;
    double user_bound_scale_value = std::pow(2, user_bound_scale);
    REQUIRE(highs.setOptionValue("user_bound_scale", user_bound_scale) ==
            HighsStatus::kOk);

    HighsInt user_cost_scale = 4;
    double user_cost_scale_value = std::pow(2, user_cost_scale);
    REQUIRE(highs.setOptionValue("user_cost_scale", user_cost_scale) ==
            HighsStatus::kOk);

    highs.run();

    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
    REQUIRE(doubleEqual0(highs.getInfo().objective_function_value,
                         unscaled_objective));

    model = mip_model;
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
