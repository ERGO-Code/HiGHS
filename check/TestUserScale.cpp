#include <cmath>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const double inf = kHighsInf;

bool doubleEqual0(const double v0, const double v1) {
  double rel_difference = std::fabs(v0 - v1) / std::max(1.0, std::fabs(v0));
  bool ok = rel_difference < 1e-12;
  if (
      //dev_run &&
      !ok)
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

    HighsInt user_objective_scale = 4;
    double user_objective_scale_value = std::pow(2, user_objective_scale);
    REQUIRE(highs.setOptionValue("user_cost_scale", user_objective_scale) ==
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
  //  highs.setOptionValue("output_flag", dev_run);
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

  HighsInt suggested_objective_scale;
  HighsInt suggested_bound_scale;
  highs.getObjectiveBoundScaling(suggested_objective_scale, suggested_bound_scale);

  highs.passModel(lp);

  highs.getObjectiveBoundScaling(suggested_objective_scale, suggested_bound_scale);

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

  std::string model = "flugpl";
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  highs.readModel(filename);

  REQUIRE(highs.setOptionValue("user_cost_scale", -30) == HighsStatus::kOk);

  highs.run();

  highs.resetGlobalScheduler(true);
}

HighsLp lp0(const double cost, const double bound) {
  // This LP is unbounded and causes assert in presolve!
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {cost, -2*cost};
  lp.col_lower_ = {0, 1e-8};
  lp.col_upper_ = {bound, bound};
  lp.row_lower_ = {-kHighsInf, bound};
  lp.row_upper_ = {bound, kHighsInf};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 1, -1};
  return lp;
}

HighsLp mip0(const double cost, const double bound) {
  HighsLp lp = lp0(cost, bound);
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kContinuous};
  return lp;
}

HighsLp lp1(const double cost, const double col_lower, const double bound) {
  HighsLp lp;
  // Set up the LP
  //
  // min -4C x -7C y
  //
  // st x + y <= 6B; -2B <= x-y
  //
  // b <= [x, y] <= 10*B
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {-4*cost, -7*cost};
  lp.col_lower_ = {col_lower, col_lower};
  lp.col_upper_ = {10*bound, 10*bound};
  lp.row_lower_ = {-kHighsInf, -2*bound};
  lp.row_upper_ = {6*bound, kHighsInf};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 1, -1};
  return lp;
}

HighsLp mip1(const double cost, const double col_lower, const double bound) {
  HighsLp lp = lp1(cost, col_lower, bound);
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kContinuous};
  lp.col_lower_[0] = 0;
  return lp;
}

HighsHessian hessian(const double value) {
  HighsHessian hessian;
  hessian.dim_ = 2;
  hessian.start_ = {0, 1, 2};
  hessian.index_ = {0, 1};
  hessian.value_ = {value, 2*value};
  return hessian;
}

void testUserScale(Highs& h) {
  printf("\nWithout user scaling\n");
  h.writeModel("");
  h.run();
  h.writeSolution("", 1);
  double unscaled_objective_value = h.getInfo().objective_function_value;
  HighsInt suggested_objective_scale;
  HighsInt suggested_bound_scale;
  h.getObjectiveBoundScaling(suggested_objective_scale, suggested_bound_scale);
  const bool has_suggested_scaling = suggested_objective_scale || suggested_bound_scale;
  if (!has_suggested_scaling) {
    suggested_objective_scale = 2;
    suggested_bound_scale = 1;
  }

  h.setOptionValue("user_cost_scale", suggested_objective_scale);
  h.setOptionValue("user_bound_scale", suggested_bound_scale);
  printf("\nWith user scaling\n");
  h.run();
  h.writeSolution("", 1);
  REQUIRE(doubleEqual0(unscaled_objective_value, h.getInfo().objective_function_value));
}


TEST_CASE("ill-scaled-model", "[highs_user_scale]") {
  Highs h;
  const HighsInfo& info = h.getInfo();
  const HighsSolution& solution = h.getSolution();
  //  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("qp_regularization_value", 0);
  // Preolve on triggers assert
  const bool expose_presolve_bug = false;
  if (expose_presolve_bug) {
    h.setOptionValue("presolve", kHighsOffString);
    h.setOptionValue("presolve_reductions", 0);
    HighsLp lp = mip0(1.0, kHighsInf);
    h.passModel(lp);
    h.run();
  }

  const bool lp_test = true;
  if (lp_test) {
    HighsLp lp = lp1(1.0, 0.0, 1.0);
    h.passModel(lp);
    testUserScale(h);
  }

  const bool mip_test = false;
  if (mip_test) {
    HighsLp lp = mip1(1.0, 0.0, 1.0);
    h.passModel(lp);
    testUserScale(h);
  }

  const bool qp_test = true;
  if (qp_test) {
    HighsModel model;
    //    model.lp_ = lp1(1e-8, 1.0);
    //    model.hessian_ = hessian(1e-8);
    model.lp_ = lp1(1.0, 0.0, 1.0);
    model.hessian_ = hessian(1.0);
    h.passModel(model);
    testUserScale(h);
  }

  h.resetGlobalScheduler(true);

}

