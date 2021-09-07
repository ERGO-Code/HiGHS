#include "Highs.h"
#include "catch.hpp"

const double inf = kHighsInf;
const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

void detailedOutput(Highs& highs);
void dualiseTest(Highs& highs);
void simpleTest(Highs& highs);
void fixedColumnTest(Highs& highs);
void freeColumnTest(Highs& highs);
void distillationTest(Highs& highs);
void afiroTest(Highs& highs);

TEST_CASE("Dualise", "[highs_test_dualise]") {

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  //simpleTest(highs);
  //fixedColumnTest(highs);
  freeColumnTest(highs);
  //distillationTest(highs);
  //afiroTest(highs);
}

void dualiseTest(Highs& highs) {
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("simplex_dualise_strategy", kHighsOptionOff);
  highs.setBasis();
  highs.run();
  //  highs.writeSolution("", true);
  double primal_objective = info.objective_function_value;
  highs.setOptionValue("simplex_dualise_strategy", kHighsOptionOn);
  highs.setBasis();
  //  detailedOutput(highs);
  highs.run();
  //  highs.writeSolution("", true);
  double dual_objective = info.objective_function_value;
  double dl = fabs(primal_objective-dual_objective);
  REQUIRE(dl < double_equal_tolerance);
}

void detailedOutput(Highs& highs) {
  if (!dev_run) return;
  highs.setOptionValue("output_flag", true);
  highs.setOptionValue("log_dev_level", 1);
  highs.setOptionValue("highs_debug_level", 2);
}

void simpleTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {2, 1};
  lp.col_lower_ = {1, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {2};
  lp.row_upper_ = {inf};
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  highs.passModel(model);
  dualiseTest(highs);
  highs.clear();
}

void fixedColumnTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {8, 10};
  lp.col_lower_ = {1, 0};
  lp.col_upper_ = {1, inf};
  lp.row_lower_ = {7, 12, 6};
  lp.row_upper_ = {inf, inf, inf};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {2, 3, 2, 2, 4, 1};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  highs.passModel(model);
  dualiseTest(highs);
  highs.clear();
}

void freeColumnTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {8, 10};
  lp.col_lower_ = {0, -inf};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {7, 12, 6};
  lp.row_upper_ = {inf, inf, inf};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {2, 3, 2, 2, 4, 1};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  highs.passModel(model);
  dualiseTest(highs);
  highs.clear();
}

void distillationTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {8, 10};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {7, 12, 6};
  lp.row_upper_ = {inf, inf, inf};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {2, 3, 2, 2, 4, 1};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  highs.passModel(model);
  dualiseTest(highs);

  double x0_lower = 3;
  if (dev_run) printf("\nGive a lower bound on x0 of %g\n", x0_lower);
  highs.changeColBounds(0, x0_lower, inf);
  dualiseTest(highs);

  double x1_upper = 0.5;
  if (dev_run) printf("\nGive an upper bound on x1 of %g\n", x1_upper);
  //  optimal solution
  highs.changeColBounds(1, -inf, x1_upper);
  dualiseTest(highs);

  highs.clear();
}

void afiroTest(Highs& highs) {
  std::string model_file =
    std::string(HIGHS_DIR) + "/check/instances/afiro.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  dualiseTest(highs);
}
