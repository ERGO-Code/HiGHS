#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const double inf = kHighsInf;
const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

TEST_CASE("pwl-identity", "[highs_test_pwl]") {
  // f(x) = x with breakpoints at (0,0) and (10,10)
  // min y s.t. y = f(x), x in [0,10], y in [0,10]
  // Optimal: x=0, y=0
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(0.0, 10.0);  // y (col 1)

  // min y
  highs.changeColCost(1, 1.0);

  double xbp[] = {0.0, 10.0};
  double ybp[] = {0.0, 10.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 1, 2, xbp, ybp) ==
          HighsStatus::kOk);
  REQUIRE(highs.getNumPwlConstraints() == 1);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 0.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 0.0) < double_equal_tolerance);
  REQUIRE(fabs(info.objective_function_value - 0.0) < double_equal_tolerance);
}

TEST_CASE("pwl-single-segment", "[highs_test_pwl]") {
  // f(x) = 2x + 1 with breakpoints at (0,1) and (5,11)
  // min y s.t. y = f(x), x in [0,5], y free
  // Optimal: x=0, y=1
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 5.0);      // x (col 0)
  highs.addVar(-inf, inf);      // y (col 1)

  highs.changeColCost(1, 1.0);

  double xbp[] = {0.0, 5.0};
  double ybp[] = {1.0, 11.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 1, 2, xbp, ybp) ==
          HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 0.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 1.0) < double_equal_tolerance);
  REQUIRE(fabs(info.objective_function_value - 1.0) < double_equal_tolerance);
}

TEST_CASE("pwl-two-segment", "[highs_test_pwl]") {
  // f(x) with breakpoints (0,0), (2,4), (4,5)
  // Slopes: s0 = (4-0)/(2-0) = 2, s1 = (5-4)/(4-2) = 0.5
  // min y s.t. y = f(x), x = 3 (fixed), y free
  // At x=3: first segment is fully used (d0=2), second segment d1=1
  // y = 0 + 2*2 + 0.5*1 = 4.5
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(3.0, 3.0);      // x (col 0), fixed at 3
  highs.addVar(-inf, inf);      // y (col 1)

  highs.changeColCost(1, 1.0);

  double xbp[] = {0.0, 2.0, 4.0};
  double ybp[] = {0.0, 4.0, 5.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 1, 3, xbp, ybp) ==
          HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 3.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 4.5) < double_equal_tolerance);
  REQUIRE(fabs(info.objective_function_value - 4.5) < double_equal_tolerance);
}

TEST_CASE("pwl-inequality", "[highs_test_pwl]") {
  // Model y >= f(x) via auxiliary w = f(x) and constraint y >= w
  // f(x) with breakpoints (0,0), (5,10) => f(x)=2x
  // min y s.t. w = f(x), y >= w, x=3, y free, w free
  // => w = 6, y >= 6 => y = 6
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(3.0, 3.0);      // x (col 0)
  highs.addVar(-inf, inf);      // y (col 1)
  highs.addVar(-inf, inf);      // w (col 2, auxiliary)

  highs.changeColCost(1, 1.0);  // min y

  // PWL: w = f(x)
  double xbp[] = {0.0, 5.0};
  double ybp[] = {0.0, 10.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 2, 2, xbp, ybp) ==
          HighsStatus::kOk);

  // Constraint: y - w >= 0 => y >= w
  HighsInt indices[] = {1, 2};
  double values[] = {1.0, -1.0};
  highs.addRow(0.0, inf, 2, indices, values);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 3.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 6.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[2] - 6.0) < double_equal_tolerance);
}

TEST_CASE("pwl-multi", "[highs_test_pwl]") {
  // Two PWL constraints:
  // y1 = f1(x) with breakpoints (0,0), (10,10)  => f1(x) = x
  // y2 = f2(x) with breakpoints (0,10), (10,0)  => f2(x) = 10 - x
  // min y1 + y2, x in [0,10]
  // y1 + y2 = x + (10-x) = 10 for all x
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);     // x  (col 0)
  highs.addVar(-inf, inf);     // y1 (col 1)
  highs.addVar(-inf, inf);     // y2 (col 2)

  highs.changeColCost(1, 1.0);
  highs.changeColCost(2, 1.0);

  double xbp1[] = {0.0, 10.0};
  double ybp1[] = {0.0, 10.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 1, 2, xbp1, ybp1) ==
          HighsStatus::kOk);

  double xbp2[] = {0.0, 10.0};
  double ybp2[] = {10.0, 0.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 2, 2, xbp2, ybp2) ==
          HighsStatus::kOk);
  REQUIRE(highs.getNumPwlConstraints() == 2);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution& solution = highs.getSolution();
  // y1 + y2 should equal 10
  REQUIRE(fabs(solution.col_value[1] + solution.col_value[2] - 10.0) <
          double_equal_tolerance);
  REQUIRE(fabs(info.objective_function_value - 10.0) < double_equal_tolerance);
}

TEST_CASE("pwl-validation", "[highs_test_pwl]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // col 0
  highs.addVar(0.0, 10.0);  // col 1

  double xbp[] = {0.0, 5.0, 10.0};
  double ybp[] = {0.0, 5.0, 10.0};

  // input_col out of range
  REQUIRE(highs.addPiecewiseLinearConstraint(-1, 1, 3, xbp, ybp) ==
          HighsStatus::kError);
  REQUIRE(highs.addPiecewiseLinearConstraint(2, 1, 3, xbp, ybp) ==
          HighsStatus::kError);

  // output_col out of range
  REQUIRE(highs.addPiecewiseLinearConstraint(0, -1, 3, xbp, ybp) ==
          HighsStatus::kError);
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 2, 3, xbp, ybp) ==
          HighsStatus::kError);

  // input == output
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 0, 3, xbp, ybp) ==
          HighsStatus::kError);

  // Too few breakpoints
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 1, 1, xbp, ybp) ==
          HighsStatus::kError);

  // Non-increasing x breakpoints
  double bad_xbp[] = {0.0, 10.0, 5.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 1, 3, bad_xbp, ybp) ==
          HighsStatus::kError);

  // No constraints should have been added
  REQUIRE(highs.getNumPwlConstraints() == 0);
}

TEST_CASE("pwl-three-segment-optimize", "[highs_test_pwl]") {
  // f(x) with breakpoints (0,10), (3,1), (6,4), (10,20)
  // Slopes: s0 = (1-10)/3 = -3, s1 = (4-1)/3 = 1, s2 = (20-4)/4 = 4
  // min y s.t. y = f(x), x in [0,10], y free
  // Minimum of f is at x=3 where f(3)=1
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(-inf, inf);  // y (col 1)

  highs.changeColCost(1, 1.0);

  double xbp[] = {0.0, 3.0, 6.0, 10.0};
  double ybp[] = {10.0, 1.0, 4.0, 20.0};
  REQUIRE(highs.addPiecewiseLinearConstraint(0, 1, 4, xbp, ybp) ==
          HighsStatus::kOk);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const HighsSolution& solution = highs.getSolution();
  REQUIRE(fabs(solution.col_value[0] - 3.0) < double_equal_tolerance);
  REQUIRE(fabs(solution.col_value[1] - 1.0) < double_equal_tolerance);
  REQUIRE(fabs(info.objective_function_value - 1.0) < double_equal_tolerance);
}
