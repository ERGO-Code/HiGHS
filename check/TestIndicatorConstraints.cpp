#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const double inf = kHighsInf;
const bool dev_run = true;  // false;//
const double double_equal_tolerance = 1e-5;

void solveWriteReadSolve(Highs& highs, const double objective_value,
                         const std::vector<double> col_value,
                         const bool have_names = false);

void solveWriteReadSolve(Highs& highs, const double objective_value,
                         const std::vector<double> col_value,
                         const bool have_names) {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  std::string filename_mps = test_name + ".mps";
  std::string filename_lp = test_name + ".lp";

  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();

  const HighsInt to_k = 2;
  for (HighsInt k = 0; k < to_k; k++) {
    if (dev_run) printf("\nRun with k = %d\n==============\n\n", int(k));
    if (k == 1) {
      REQUIRE(highs.readModel(filename_mps) == HighsStatus::kOk);
    } else if (k == 2) {
      REQUIRE(highs.readModel(filename_lp) == HighsStatus::kOk);
    }
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

    if (dev_run)
      printf("Solution is (%g, %g) with objective %g\n", solution.col_value[0],
             solution.col_value[1], info.objective_function_value);
    if (dev_run)
      printf("Testing  is (%g, %g) with objective %g\n", col_value[0],
             col_value[1], objective_value);
    REQUIRE(fabs(solution.col_value[0] - col_value[0]) <
            double_equal_tolerance);
    REQUIRE(fabs(solution.col_value[1] - col_value[1]) <
            double_equal_tolerance);
    REQUIRE(fabs(info.objective_function_value - objective_value) <
            double_equal_tolerance);

    if (k == 0) {
      REQUIRE(highs.writeModel(filename_lp) ==
              (have_names ? HighsStatus::kOk : HighsStatus::kWarning));
      REQUIRE(highs.writeModel(filename_mps) == HighsStatus::kOk);
    }
  }
  //  std::remove(filename_lp.c_str());
  // std::remove(filename_mps.c_str());
}

TEST_CASE("indicator-simple-v1", "[highs_test_indicator]") {
  // min x
  // s.t. z=1 -> x >= 5
  //      z binary, x in [0, 10]
  // Optimal: z=0, x=0, obj=0
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(0.0, 1.0);   // z (col 1)

  // Costs: min x
  highs.changeColCost(0, 1.0);
  highs.changeColCost(1, 0.0);

  // z is integer (binary)
  highs.changeColIntegrality(1, HighsVarType::kInteger);

  // Indicator: z=1 -> x >= 5
  HighsInt indices[] = {0};
  double values[] = {1.0};
  REQUIRE(highs.addIndicatorConstraint(1, 1, 1, indices, values, 5.0, inf) ==
          HighsStatus::kOk);

  REQUIRE(highs.getNumIndicatorConstraints() == 1);

  const std::vector<double> col_value = {0, 0};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-simple-v0", "[highs_test_indicator]") {
  // min x
  // s.t. z=0 -> x >= 5
  //      z binary, x in [0, 10]
  // Optimal: z=1, x=0, obj=0
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(0.0, 1.0);   // z (col 1)

  highs.changeColCost(0, 1.0);
  highs.changeColCost(1, 0.0);
  highs.changeColIntegrality(1, HighsVarType::kInteger);

  // Indicator: z=0 -> x >= 5
  HighsInt indices[] = {0};
  double values[] = {1.0};
  REQUIRE(highs.addIndicatorConstraint(1, 0, 1, indices, values, 5.0, inf) ==
          HighsStatus::kOk);

  const std::vector<double> col_value = {0, 1};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-max-big-m", "[highs_test_indicator]") {}

TEST_CASE("indicator-range", "[highs_test_indicator]") {
  // min x
  // s.t. z=1 -> 3 <= x <= 7
  //      z binary, x in [0, 10]
  // Optimal: z=0, x=0, obj=0
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(0.0, 1.0);   // z (col 1)

  highs.changeColCost(0, 1.0);
  highs.changeColCost(1, 0.0);
  highs.changeColIntegrality(1, HighsVarType::kInteger);

  // Indicator: z=1 -> 3 <= x <= 7
  HighsInt indices[] = {0};
  double values[] = {1.0};
  REQUIRE(highs.addIndicatorConstraint(1, 1, 1, indices, values, 3.0, 7.0) ==
          HighsStatus::kOk);

  const std::vector<double> col_value = {0, 0};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-forced", "[highs_test_indicator]") {
  // min x
  // s.t. z == 1 (fixed)
  //      z=1 -> x >= 5
  //      x in [0, 10]
  // Optimal: x=5, obj=5
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(1.0, 1.0);   // z (col 1, fixed to 1)

  highs.changeColCost(0, 1.0);
  highs.changeColCost(1, 0.0);
  highs.changeColIntegrality(1, HighsVarType::kInteger);

  HighsInt indices[] = {0};
  double values[] = {1.0};
  REQUIRE(highs.addIndicatorConstraint(1, 1, 1, indices, values, 5.0, inf) ==
          HighsStatus::kOk);

  const std::vector<double> col_value = {5, 1};
  const double objective_value = 5;

  solveWriteReadSolve(highs, objective_value, col_value);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-multiple", "[highs_test_indicator]") {
  // min x + y
  // s.t. z=1 -> x >= 3
  //      z=1 -> y >= 4
  //      z binary, x in [0,10], y in [0,10]
  // Optimal: z=0, x=0, y=0, obj=0
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(0.0, 10.0);  // y (col 1)
  highs.addVar(0.0, 1.0);   // z (col 2)

  highs.changeColCost(0, 1.0);
  highs.changeColCost(1, 1.0);
  highs.changeColCost(2, 0.0);
  highs.changeColIntegrality(2, HighsVarType::kInteger);

  // z=1 -> x >= 3
  HighsInt idx0[] = {0};
  double val0[] = {1.0};
  REQUIRE(highs.addIndicatorConstraint(2, 1, 1, idx0, val0, 3.0, inf) ==
          HighsStatus::kOk);

  // z=1 -> y >= 4
  HighsInt idx1[] = {1};
  double val1[] = {1.0};
  REQUIRE(highs.addIndicatorConstraint(2, 1, 1, idx1, val1, 4.0, inf) ==
          HighsStatus::kOk);

  REQUIRE(highs.getNumIndicatorConstraints() == 2);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const std::vector<double> col_value = {0, 0, 0};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-validation", "[highs_test_indicator]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  highs.addVar(0.0, 10.0);  // x (col 0)
  highs.addVar(0.0, 1.0);   // z (col 1)
  highs.changeColIntegrality(1, HighsVarType::kInteger);

  HighsInt indices[] = {0};
  double values[] = {1.0};

  // Invalid binary_col
  REQUIRE(highs.addIndicatorConstraint(-1, 1, 1, indices, values, 0.0, inf) ==
          HighsStatus::kError);
  REQUIRE(highs.addIndicatorConstraint(5, 1, 1, indices, values, 0.0, inf) ==
          HighsStatus::kError);

  // Invalid binary_value
  REQUIRE(highs.addIndicatorConstraint(1, 2, 1, indices, values, 0.0, inf) ==
          HighsStatus::kError);

  // Non-integer variable
  REQUIRE(highs.addIndicatorConstraint(0, 1, 1, indices, values, 0.0, inf) ==
          HighsStatus::kError);

  // Valid call should succeed
  REQUIRE(highs.addIndicatorConstraint(1, 1, 1, indices, values, 5.0, inf) ==
          HighsStatus::kOk);
  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-mps", "[highs_test_indicator]") {
  // Test reading an MPS file with INDICATORS section
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/indicator1.mps";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);
  HighsStatus status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getNumIndicatorConstraints() == 1);

  const std::vector<double> col_value = {0, 0};
  const double objective_value = 0;

  const bool have_names = true;
  solveWriteReadSolve(highs, objective_value, col_value, have_names);

  highs.resetGlobalScheduler(true);
}
