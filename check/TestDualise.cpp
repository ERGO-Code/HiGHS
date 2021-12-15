#include "Highs.h"
#include "catch.hpp"

const double inf = kHighsInf;
const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

void detailedOutput(Highs& highs);
void dualiseTest(Highs& highs);
void simpleTest(Highs& highs);
void fixedColumnTest(Highs& highs);
void freeColumnTest(Highs& highs);
void colUpperBoundTest(Highs& highs);
void rowUpperBoundTest(Highs& highs);
void distillationTest(Highs& highs);
HighsLp distillationLp();
void instanceTest(Highs& highs, const std::string model_name);

TEST_CASE("Dualise", "[highs_test_dualise]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  // simpleTest(highs);
  // distillationTest(highs);
  // freeColumnTest(highs);
  // fixedColumnTest(highs);
  //  colUpperBoundTest(highs);
  //  rowUpperBoundTest(highs);
  //  instanceTest(highs, "avgas");
  //  instanceTest(highs, "afiro");
  //  instanceTest(highs, "adlittle");
  instanceTest(highs, "25fv47");
}

void dualiseTest(Highs& highs) {
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("presolve", "off");
  highs.setOptionValue("simplex_dualise_strategy", kHighsOptionOff);
  highs.setBasis();
  highs.run();
  //  if (dev_run) highs.writeSolution("", true);
  double primal_objective = info.objective_function_value;
  highs.setOptionValue("simplex_dualise_strategy", kHighsOptionOn);
  highs.setBasis();
  //  detailedOutput(highs);
  highs.run();
  // if (dev_run) highs.writeSolution("", true);
  double dual_objective = info.objective_function_value;
  double dl = fabs(primal_objective - dual_objective);
  REQUIRE(dl < double_equal_tolerance);
}

HighsLp distillationLp() {
  HighsLp lp;
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
  return lp;
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

void distillationTest(Highs& highs) {
  HighsModel model;
  model.lp_ = distillationLp();
  highs.passModel(model);
  dualiseTest(highs);

  double x0_lower = 3;
  if (dev_run) printf("\nGive a lower bound on x0 of %g\n", x0_lower);
  highs.changeColBounds(0, x0_lower, inf);
  dualiseTest(highs);

  double x1_upper = 0.5;
  if (dev_run) printf("\nGive an upper bound on x1 of %g\n", x1_upper);
  highs.changeColBounds(1, -inf, x1_upper);
  dualiseTest(highs);

  highs.clear();
}

void freeColumnTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp = distillationLp();
  if (dev_run) printf("\nFree column 1 of distillation\n");
  lp.col_lower_[1] = -inf;
  highs.passModel(model);
  dualiseTest(highs);
  highs.clear();
}

void fixedColumnTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp = distillationLp();
  double x0_fixed = 1;
  if (dev_run) printf("\nFix column 0 of distillation to be %g\n", x0_fixed);
  lp.col_lower_[0] = x0_fixed;
  lp.col_upper_[0] = x0_fixed;
  highs.passModel(model);
  dualiseTest(highs);
  highs.clear();
}

void colUpperBoundTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp = distillationLp();
  double col1_upper = 1;
  if (dev_run) printf("\nGive an upper bound on col 1 of %g\n", col1_upper);
  lp.col_upper_[1] = col1_upper;
  // Needs reduced lower bound for feasiblilty
  //  double col2_lower = 5.7; lp.col_lower_[2] = col2_lower;
  highs.passModel(model);
  dualiseTest(highs);
}

void rowUpperBoundTest(Highs& highs) {
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp = distillationLp();
  double row0_upper = 7.1;
  if (dev_run) printf("\nGive an upper bound on row 0 of %g\n", row0_upper);
  lp.row_upper_[0] = row0_upper;
  // Needs reduced lower bound for feasiblilty
  double row2_lower = 5.7;
  lp.row_lower_[2] = row2_lower;
  highs.passModel(model);
  dualiseTest(highs);
}

void instanceTest(Highs& highs, const std::string model_name) {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model_name + ".mps";
  if (dev_run) printf("\nSolving model %s\n", model_name.c_str());
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  dualiseTest(highs);
}
