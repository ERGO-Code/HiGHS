#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;  // true;  //
const bool write_model = false;

const double inf = kHighsInf;

void testMps(std::string& model, const HighsInt iis_strategy,
             const HighsModelStatus require_model_status =
                 HighsModelStatus::kInfeasible);

void testFeasibilityRelaxation(
    std::string& model, const double lower_penalty, const double upper_penalty,
    const double rhs_penalty,
    const double require_feasibility_objective_function_value);

TEST_CASE("lp-incompatible-bounds", "[iis]") {
  // LP has row0 and col2 with inconsistent bounds.
  //
  // When prioritising rows, row0 and its constituent columns (1, 2) should be
  // found
  //
  // When prioritising columns, col2 and its constituent rows (0, 1) should be
  // found
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {0, 0, 0};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {1, 1, -1};
  lp.row_lower_ = {1, 0};
  lp.row_upper_ = {0, 1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {1, 2, 0, 2};
  lp.a_matrix_.value_ = {1, 1, 1, 1};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  // Perform the default (light) IIS check
  HighsIis iis;
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  // Perform full IIS
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  highs.setOptionValue("iis_strategy", kIisStrategyFromLpRowPriority);
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == 0);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusBoxed);

  highs.setOptionValue("iis_strategy", kIisStrategyFromLpColPriority);
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(iis.col_index_.size() == 1);
  REQUIRE(iis.row_index_.size() == 0);
  REQUIRE(iis.col_index_[0] == 2);
  REQUIRE(iis.col_bound_[0] == kIisBoundStatusBoxed);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("lp-empty-infeasible-row", "[iis]") {
  // Second row is empty, with bounds of [1, 2]
  const HighsInt empty_row = 1;
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf, 1, -inf};
  lp.row_upper_ = {8, 2, 9};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {2, 1, 1, 3};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);

  HighsIis iis;
  // Get IIS for empty row with positive lower bound
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  if (dev_run && write_model) {
    highs.writeModel("");
    highs.writeIisModel("");
  }
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == empty_row);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusLower);

  // Get IIS for empty row with negative upper bound
  double new_lower = -2;
  double new_upper = -1;
  assert(new_upper < 0);
  REQUIRE(highs.changeRowBounds(empty_row, new_lower, new_upper) ==
          HighsStatus::kOk);
  lp.row_lower_[empty_row] = new_lower;
  lp.row_upper_[empty_row] = new_upper;
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == empty_row);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusUpper);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("lp-get-iis-light", "[iis]") {
  HighsLp lp;
  lp.model_name_ = "lp-get-iis-light";
  lp.num_col_ = 4;
  lp.num_row_ = 3;
  lp.col_cost_ = {0, 0, 0, 0};
  lp.col_lower_ = {10, 10, 10, 10};
  lp.col_upper_ = {20, 20, 20, 20};
  lp.row_lower_ = {-inf, -10, -34};
  lp.row_upper_ = {30, 15, inf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 3, 7, 10};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2, 3, 1, 2, 3};
  lp.a_matrix_.value_ = {1.5, 2, 1, 4, -2, 1, 2, -2, -1.5, -1};
  //
  // Each of the three constraints constitutes an IIS. Find each by
  // freeing the upper bounds in turn.
  //
  // 1.5w + 2x + y <= 30
  //
  // -10 <= 4w -2x + y + 2z <= 15
  //
  // -2x -1.5y -z >= -34
  //
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  HighsIis iis;

  for (int l = 0; l < 3; l++) {
    for (int k = 0; k < 2; k++) {
      REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
      REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
      if (dev_run && write_model) {
        highs.writeModel("");
        highs.writeIisModel("");
        highs.writeIisModel(lp.model_name_ + ".lp");
      }
      if (k == 0) {
        // Now flip to column-wise for code coverage
        lp.a_matrix_.ensureColwise();
        highs.passModel(lp);
      }
    }
    // Prevent the first and then the second constraints from being an
    // IIS
    if (l == 0) {
      lp.row_upper_[0] = inf;
    } else if (l == 1) {
      lp.row_upper_[1] = inf;
    }
    highs.passModel(lp);
  }
  highs.resetGlobalScheduler(true);
}

TEST_CASE("lp-get-iis", "[iis]") {
  HighsLp lp;
  lp.model_name_ = "lp-get-iis";
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf, -inf, -inf};
  lp.row_upper_ = {8, 9, -2};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 2, 4, 6};
  lp.a_matrix_.index_ = {0, 1, 0, 1, 0, 1};
  lp.a_matrix_.value_ = {2, 1, 1, 3, 1, 1};
  //
  // 2x + y <= 8
  //
  // x + 3y <= 9
  //
  // x +  y <= -2
  //
  // x, y \in [0, inf)
  //
  //  lp.col_name_ = {"Col0", "Col1"};
  //  lp.row_name_ = {"Row0", "Row1", "Row2"};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  HighsIis iis;
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(iis.col_index_.size() == 2);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.col_index_[0] == 0);
  REQUIRE(iis.col_index_[1] == 1);
  REQUIRE(iis.row_index_[0] == 2);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("lp-get-iis-woodinfe", "[iis]") {
  std::string model = "woodinfe";
  testMps(model, kIisStrategyFromLpRowPriority);
  //  testMps(model, kIisStrategyFromRayRowPriority);
}

TEST_CASE("lp-get-iis-galenet", "[iis]") {
  // Dual ray corresponds to constraints
  //
  // r0: 0  <= c0 + c1      - c3 - c4 <=0
  //
  // r1: 20 <=           c2 + c3
  //
  // r2: 30 <=                     c4
  //
  // Where
  //
  // 0 <= c0 <= 10
  //
  // 0 <= c1 <= 10
  //
  // 0 <= c2 <=  2
  //
  // 0 <= c3 <= 20
  //
  // 0 <= c4 <= 30
  //
  // This is infeasible since c4 >= 30 and c4 <= 30 fices c4 = 30,
  // then c0 + c1 >= c3 + c4 >= 30 cannot be satisfied due to the
  // upper bounds of 10 on these variables
  //
  // r1 can be removed and infeasibility is retained, but not r0 or r2
  //
  // The upper bound on r0 can be removed
  //
  // The lower bounds on c0 and c1 can be removed, but not their upper
  // bounds
  //
  // c2 can be removed, as it is empty once r1 is removed
  //
  // c3 can be removed, as the value of c4 is sufficient to make r0
  // infeasible
  //
  // The bounds on c4 can be removed, since it's the lower bound from
  // r2 that makes r0 infeasible
  //
  // Hence only empty columns can be removed
  std::string model = "galenet";
  testMps(model, kIisStrategyFromLpRowPriority);
  //  testMps(model, kIisStrategyFromRayRowPriority);
}

TEST_CASE("lp-get-iis-avgas", "[iis]") {
  std::string model = "avgas";
  // For the whole LP calculation the elasticity filter only
  // identified feasibility, so the model status is not set
  testMps(model, kIisStrategyFromLpRowPriority, HighsModelStatus::kNotset);
  // For the ray calculation the model is solved, so its status is
  // known
  //  testMps(model, kIisStrategyFromRayRowPriority,
  //  HighsModelStatus::kOptimal);
}

TEST_CASE("lp-feasibility-relaxation", "[iis]") {
  // Using infeasible LP from AMPL documentation
  HighsLp lp;
  lp.model_name_ = "ampl_infeas";
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {1, -2};
  lp.col_lower_ = {5, -inf};
  lp.col_upper_ = {inf, inf};
  lp.col_names_ = {"X", "Y"};
  lp.row_lower_ = {2, -inf, -inf};
  lp.row_upper_ = {inf, 1, 20};
  lp.row_names_ = {"R0", "R1", "R2"};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {-1, -3, 20, 21, 2, 1};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  const HighsSolution& solution = h.getSolution();
  h.passModel(lp);
  //  h.run();

  const bool all_tests = false;
  const bool test0 = false || all_tests;
  const bool test1 = false || all_tests;
  const bool test2 = false || all_tests;
  const bool test3 = false || all_tests;
  if (test0) {
    // Vanilla feasibility relaxation
    if (dev_run)
      printf(
          "==============================\n"
          "Vanilla feasibility relaxation\n"
          "==============================\n");
    REQUIRE(h.feasibilityRelaxation(1, 1, 1) == HighsStatus::kOk);
    // Should get solution (1, 1)
    h.writeSolution("", 1);
    REQUIRE(solution.col_value[0] == 1);
    REQUIRE(solution.col_value[1] == 1);
  }

  if (test1) {
    // Now we want to force all variable lower bounds to be respected
    if (dev_run)
      printf(
          "========================\n"
          "Respect all lower bounds\n"
          "========================\n");
    h.feasibilityRelaxation(-1, 1, 1);
    // Should get solution (5, 1)
    h.writeSolution("", 1);
    REQUIRE(solution.col_value[0] == 5);
    REQUIRE(solution.col_value[1] == 1);
  }

  if (test2) {
    if (dev_run)
      printf(
          "===============================\n"
          "Local penalties RHS {1, -1, 10}\n"
          "===============================\n");
    // Now test local penalties
    //
    // constraint 0: normal weight
    //
    // constraint 1: cannot be violated
    //
    // constraint 2: rather not violate
    //
    std::vector<double> local_rhs_penalty = {1, -1, 10};
    h.feasibilityRelaxation(1, 1, 0, nullptr, nullptr,
                            local_rhs_penalty.data());
    // Should get slacks (-3, 4, 0)
    h.writeSolution("", 1);
    REQUIRE(solution.row_value[0] == lp.row_lower_[0] - 3);
    REQUIRE(solution.row_value[1] == lp.row_upper_[1] - 4);
    REQUIRE(solution.row_value[2] == lp.row_upper_[2]);
  }

  if (test3) {
    if (dev_run)
      printf(
          "==============================\n"
          "Local penalties RHS {10, 1, 1}\n"
          "==============================\n");
    //
    // constraint 0: rather not violate
    //
    // constraint 1: normal weight
    //
    // constraint 2: normal weight
    //
    std::vector<double> local_rhs_penalty = {10, 1, 1};
    h.feasibilityRelaxation(1, 1, 0, nullptr, nullptr,
                            local_rhs_penalty.data());
    // Should get slacks (18, 2, -1)
    REQUIRE(solution.row_value[0] == lp.row_lower_[0] + 18);
    REQUIRE(solution.row_value[1] == lp.row_upper_[1] - 2);
    REQUIRE(solution.row_value[2] == lp.row_upper_[2] + 1);
  }
  std::string model = "galenet";
  testFeasibilityRelaxation(model, 1, 1, 1, 28.0);
  model = "woodinfe";
  testFeasibilityRelaxation(model, 1, 1, 1, 15.0);
  model = "avgas";
  testFeasibilityRelaxation(model, 1, 1, 1, 0);

  h.resetGlobalScheduler(true);
}

void testMps(std::string& model, const HighsInt iis_strategy,
             const HighsModelStatus require_model_status) {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  //  if (iis_strategy == kIisStrategyFromRayRowPriority ||
  //      iis_strategy == kIisStrategyFromRayColPriority) {
  //    // For a ray strategy, solve the LP first
  //    REQUIRE(highs.run() == HighsStatus::kOk);
  //  }
  highs.setOptionValue("iis_strategy", iis_strategy);
  HighsIis iis;
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  if (dev_run && write_model) {
    highs.writeModel("");
    highs.writeIisModel("");
  }
  HighsInt num_iis_col = iis.col_index_.size();
  HighsInt num_iis_row = iis.row_index_.size();
  HighsModelStatus model_status = highs.getModelStatus();
  REQUIRE(model_status == require_model_status);
  if (model_status == HighsModelStatus::kInfeasible) {
    REQUIRE(num_iis_col > 0);
    REQUIRE(num_iis_row > 0);
    if (dev_run)
      printf("Model %s has IIS with %d columns and %d rows\n", model.c_str(),
             int(num_iis_col), int(num_iis_row));
  } else {
    REQUIRE(num_iis_col == 0);
    REQUIRE(num_iis_row == 0);
  }
}

void testFeasibilityRelaxation(
    std::string& model, const double lower_penalty, const double upper_penalty,
    const double rhs_penalty,
    const double require_feasibility_objective_function_value) {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.readModel(model_file);
  REQUIRE(h.feasibilityRelaxation(lower_penalty, upper_penalty, rhs_penalty) ==
          HighsStatus::kOk);
  REQUIRE(h.getInfo().objective_function_value ==
          require_feasibility_objective_function_value);
}
