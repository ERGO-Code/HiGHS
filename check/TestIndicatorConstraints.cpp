#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const double inf = kHighsInf;
const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

void solveWriteReadSolve(Highs& highs, const double objective_value,
                         const double* col_value,
                         const bool have_names = false);

void solveWriteReadSolve(Highs& highs, const double objective_value,
                         const double* col_value, const bool have_names) {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  std::string filename_mps = test_name + ".mps";
  std::string filename_lp = test_name + ".lp";

  const HighsInfo& info = highs.getInfo();
  const HighsSolution& solution = highs.getSolution();

  const HighsInt to_k = 3;
  for (HighsInt k = 0; k < to_k; k++) {
    if (dev_run) printf("\nRun with k = %d\n==============\n\n", int(k));
    if (k == 1) {
      REQUIRE(highs.readModel(filename_mps) == HighsStatus::kOk);
    } else if (k == 2) {
      REQUIRE(highs.readModel(filename_lp) == HighsStatus::kOk);
    }
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

    if (col_value) {
      if (dev_run)
        printf("Solution is (%g, %g) with objective %g\n",
               solution.col_value[0], solution.col_value[1],
               info.objective_function_value);
      if (dev_run)
        printf("Testing  is (%g, %g) with objective %g\n", col_value[0],
               col_value[1], objective_value);
      REQUIRE(fabs(solution.col_value[0] - col_value[0]) <
              double_equal_tolerance);
      REQUIRE(fabs(solution.col_value[1] - col_value[1]) <
              double_equal_tolerance);
    }
    REQUIRE(fabs(info.objective_function_value - objective_value) <
            double_equal_tolerance);

    if (k == 0) {
      REQUIRE(highs.writeModel("") ==
              (have_names ? HighsStatus::kOk : HighsStatus::kWarning));
      REQUIRE(highs.writeModel(filename_lp) == HighsStatus::kOk);
      REQUIRE(highs.writeModel(filename_mps) == HighsStatus::kOk);
    }
  }
  std::remove(filename_lp.c_str());
  std::remove(filename_mps.c_str());
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
  REQUIRE(highs.addIndicatorConstraint(1, 1, 5.0, inf, 1, indices, values) ==
          HighsStatus::kOk);

  REQUIRE(highs.getNumIndicatorConstraints() == 1);

  // min
  //  obj: +1 c0
  // st
  //  r0:c1 = 1 ->  +1 c0 >= +5
  // bounds
  //  c0 <= 10
  //  c1 <= 1
  // bin
  //  c1
  // end

  const std::vector<double> col_value = {0, 0};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value.data());

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
  REQUIRE(highs.addIndicatorConstraint(1, 0, 5.0, inf, 1, indices, values) ==
          HighsStatus::kOk);

  const std::vector<double> col_value = {0, 1};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value.data());

  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-max-big-m", "[highs_test_indicator]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  HighsInt x = 0;
  HighsInt y = 1;
  HighsInt z = 2;
  highs.addVar(0.0, 1.0);  // x (col 0)
  highs.addVar(0.0, 1.0);  // y (col 1)
  highs.addVar(0.0, 1.0);  // z (col 2)

  highs.changeColIntegrality(2, HighsVarType::kInteger);

  // Loop through the cases that yield just infinite upper big-M; just
  // infinite lower big-M and both
  //
  // Indicator: z=0 -> 3 <= -x + y <= 7
  HighsInt indices[] = {x, y};
  double values[] = {-1.0, 1.0};
  for (HighsInt k = 0; k < 3; k++) {
    if (k == 0) {
      // Infinite upper bound on x => infinite lower big-M
      highs.changeColBounds(x, 0.0, kHighsInf);
      highs.changeColBounds(y, 0.0, 1.0);
    } else if (k == 1) {
      // Infinite upper bound on y => infinite upper big-M
      highs.changeColBounds(x, 0.0, 1.0);
      highs.changeColBounds(y, 0.0, kHighsInf);
    } else {
      // Infinite upper bounds on x and y => infinite lower big-M and
      // upper big-M
      highs.changeColBounds(x, 0.0, kHighsInf);
      highs.changeColBounds(y, 0.0, kHighsInf);
    }
    REQUIRE(highs.addIndicatorConstraint(z, 0, 3.0, 7.0, 2, indices, values) ==
            HighsStatus::kOk);

    REQUIRE(highs.run() == HighsStatus::kError);
    highs.clearIndicatorConstraints();
  }

  highs.resetGlobalScheduler(true);
}

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
  REQUIRE(highs.addIndicatorConstraint(1, 1, 3.0, 7.0, 1, indices, values) ==
          HighsStatus::kOk);

  const std::vector<double> col_value = {0, 0};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value.data());

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
  REQUIRE(highs.addIndicatorConstraint(1, 1, 5.0, inf, 1, indices, values) ==
          HighsStatus::kOk);

  const std::vector<double> col_value = {5, 1};
  const double objective_value = 5;

  solveWriteReadSolve(highs, objective_value, col_value.data());

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
  REQUIRE(highs.addIndicatorConstraint(2, 1, 3.0, inf, 1, idx0, val0) ==
          HighsStatus::kOk);

  // z=1 -> y >= 4
  HighsInt idx1[] = {1};
  double val1[] = {1.0};
  REQUIRE(highs.addIndicatorConstraint(2, 1, 4.0, inf, 1, idx1, val1) ==
          HighsStatus::kOk);

  REQUIRE(highs.getNumIndicatorConstraints() == 2);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  const std::vector<double> col_value = {0, 0, 0};
  const double objective_value = 0;

  solveWriteReadSolve(highs, objective_value, col_value.data());

  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-validation", "[highs_test_indicator]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsIndicatorConstraints& ic = highs.getIndicatorConstraints();
  const std::vector<HighsInt>& ic_start = ic.matrix.start_;
  const std::vector<HighsInt>& ic_index = ic.matrix.index_;
  const std::vector<double>& ic_value = ic.matrix.value_;

  highs.addVar(0.0, 1.0);  // w (col 0)
  highs.addVar(0.0, 2.0);  // x (col 1)
  highs.addVar(0.0, 3.0);  // y (col 2)
  highs.addVar(0.0, 1.0);  // z (col 3)

  const HighsInt w = 0;
  const HighsInt x = 1;
  const HighsInt y = 2;
  const HighsInt z = 3;
  const HighsInt invalid_var0 = -1;
  const HighsInt invalid_var1 = 5;
  const HighsInt valid_binary_value = 1;
  const HighsInt invalid_binary_value = 2;

  highs.changeColIntegrality(z, HighsVarType::kInteger);

  HighsInt index[] = {w, x, y};
  double value[] = {1.0, 1.0, 1.0};
  HighsInt nnz = 3;

  HighsInt num_ic = 0;

  // Invalid binary_col
  REQUIRE(highs.addIndicatorConstraint(invalid_var0, valid_binary_value, 0.0, inf, nnz, index, value) ==
          HighsStatus::kError);
  REQUIRE(highs.getNumIndicatorConstraints() == num_ic);

  REQUIRE(highs.addIndicatorConstraint(invalid_var1, valid_binary_value, 0.0, inf, nnz, index, value) ==
          HighsStatus::kError);
  REQUIRE(ic.numIndicatorConstraints() == num_ic);

  // Invalid binary_value
  REQUIRE(highs.addIndicatorConstraint(z, invalid_binary_value, 0.0, inf, nnz, index, value) ==
          HighsStatus::kError);
  REQUIRE(ic.numIndicatorConstraints() == num_ic);

  // Non-integer variable
  REQUIRE(highs.addIndicatorConstraint(x, 1, 0.0, inf, nnz, index, value) ==
          HighsStatus::kError);
  REQUIRE(ic.numIndicatorConstraints() == num_ic);

  // Free indicator constraint
  REQUIRE(highs.addIndicatorConstraint(z, valid_binary_value, -inf, inf, nnz, index, value) ==
          HighsStatus::kError);
  REQUIRE(ic.numIndicatorConstraints() == num_ic);

  // Valid call should succeed
  REQUIRE(highs.addIndicatorConstraint(z, valid_binary_value, 5.0, inf, nnz, index, value) ==
          HighsStatus::kOk);
  num_ic++;
  REQUIRE(ic.numIndicatorConstraints() == num_ic);

  // Now pass a name
  std::string name = "IC0";

  // Small constraint coefficient should lead to warning
  value[1] = 1e-15;
  REQUIRE(highs.addIndicatorConstraint(z, valid_binary_value, 5.0, inf, nnz, index, value, name) ==
          HighsStatus::kWarning);
  num_ic++;
  REQUIRE(ic.numIndicatorConstraints() == num_ic);
  REQUIRE(ic_start[num_ic]-ic_start[num_ic-1] == 2);

  name = "IC1";
  // Large constraint coefficient should lead to error
  value[1] = 1e15;
  REQUIRE(highs.addIndicatorConstraint(z, valid_binary_value, 5.0, inf, nnz, index, value, name) ==
          HighsStatus::kError);
  REQUIRE(ic.numIndicatorConstraints() == num_ic);
  value[1] = 1;

  // Illegal index should lead to error
  index[1] = invalid_var0;
  REQUIRE(highs.addIndicatorConstraint(z, valid_binary_value, 5.0, inf, nnz, index, value, name) ==
          HighsStatus::kError);
  REQUIRE(ic.numIndicatorConstraints() == num_ic);

  // Duplicated index should lead to warning
  index[1] = w;
  REQUIRE(highs.addIndicatorConstraint(z, valid_binary_value, 5.0, inf, nnz, index, value, name) ==
          HighsStatus::kWarning);
  num_ic++;
  REQUIRE(ic.numIndicatorConstraints() == num_ic);
  REQUIRE(ic_start[num_ic]-ic_start[num_ic-1] == 2);
  HighsInt iEl = ic_start[num_ic];
  REQUIRE(ic_index[iEl] == w);
  REQUIRE(ic_value[iEl] == value[0] + value[1]);
  highs.resetGlobalScheduler(true);
}

TEST_CASE("indicator-mps", "[highs_test_indicator]") {
  // Test reading an MPS file with INDICATORS section
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/water.mps";
  Highs highs;
  const HighsInfo& info = highs.getInfo();
  highs.setOptionValue("output_flag", dev_run);
  HighsStatus status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);
  REQUIRE(highs.getNumIndicatorConstraints() == 280);

  const double objective_value = 20;

  const bool have_names = true;
  solveWriteReadSolve(highs, objective_value, nullptr, have_names);

  highs.resetGlobalScheduler(true);
}
