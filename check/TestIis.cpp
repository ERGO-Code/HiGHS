#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const double inf = kHighsInf;

void testIis(const std::string& model, const HighsIis& iis);

void testMps(std::string& model, const HighsInt iis_strategy,
             const HighsModelStatus require_model_status =
                 HighsModelStatus::kInfeasible);

void testFeasibilityRelaxation(
    std::string& model, const double lower_penalty, const double upper_penalty,
    const double rhs_penalty,
    const double require_feasibility_objective_function_value);

void checkIisLp(HighsLp& lp, const HighsIis& iis, const HighsLp& iis_lp);
		
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
  //  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  // Perform the light IIS check
  highs.setOptionValue("iis_strategy", kIisStrategyLight);
  HighsLp iis_lp;
  REQUIRE(highs.getIisLp(iis_lp) == HighsStatus::kOk);
  HighsIis iis;
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  checkIisLp(lp, iis, iis_lp);  
  /*
  // Perform full IIS
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  highs.setOptionValue("iis_strategy", kIisStrategyFromLpRowPriority);
  HighsIis iis;
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
  */
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
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  HighsIis iis;
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == empty_row);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusLower);
  REQUIRE(highs.changeRowBounds(empty_row, -2, -1) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == empty_row);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusUpper);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("lp-get-iis", "[iis]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf, -inf, -inf};
  lp.row_upper_ = {8, 9, -2};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 4, 6};
  lp.a_matrix_.index_ = {0, 1, 0, 1, 0, 1};
  lp.a_matrix_.value_ = {2, 1, 1, 3, 1, 1};
  //  lp.col_name_ = {"Col0", "Col1"};
  //  lp.row_name_ = {"Row0", "Row1", "Row2"};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  highs.setOptionValue("iis_strategy", kIisStrategyFromLpRowPriority);
  HighsIis iis;
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
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

void testIis(const std::string& model, const HighsIis& iis) {
  HighsModelStatus model_status;
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  HighsInt num_iis_col = iis.col_index_.size();
  HighsInt num_iis_row = iis.row_index_.size();

  Highs highs;
  highs.setOptionValue("output_flag", false);
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  const HighsLp& incumbent_lp = highs.getLp();
  HighsLp lp = highs.getLp();
  // Zero the objective
  lp.col_cost_.assign(lp.num_col_, 0);
  REQUIRE(highs.changeColsCost(0, lp.num_col_ - 1, lp.col_cost_.data()) ==
          HighsStatus::kOk);

  // Save the bounds
  std::vector<double> iis_col_lower;
  std::vector<double> iis_col_upper;
  std::vector<double> iis_row_lower;
  std::vector<double> iis_row_upper;
  for (HighsInt iisCol = 0; iisCol < num_iis_col; iisCol++) {
    HighsInt iCol = iis.col_index_[iisCol];
    iis_col_lower.push_back(lp.col_lower_[iCol]);
    iis_col_upper.push_back(lp.col_upper_[iCol]);
  }
  for (HighsInt iisRow = 0; iisRow < num_iis_row; iisRow++) {
    HighsInt iRow = iis.row_index_[iisRow];
    iis_row_lower.push_back(lp.row_lower_[iRow]);
    iis_row_upper.push_back(lp.row_upper_[iRow]);
  }

  // Free all the columns and rows
  lp.col_lower_.assign(lp.num_col_, -inf);
  lp.col_upper_.assign(lp.num_col_, inf);
  lp.row_lower_.assign(lp.num_row_, -inf);
  lp.row_upper_.assign(lp.num_row_, inf);
  // Restore the bounds for the IIS columns and rows
  for (HighsInt iisCol = 0; iisCol < num_iis_col; iisCol++) {
    HighsInt iCol = iis.col_index_[iisCol];
    lp.col_lower_[iCol] = iis_col_lower[iisCol];
    lp.col_upper_[iCol] = iis_col_upper[iisCol];
  }
  for (HighsInt iisRow = 0; iisRow < num_iis_row; iisRow++) {
    HighsInt iRow = iis.row_index_[iisRow];
    lp.row_lower_[iRow] = iis_row_lower[iisRow];
    lp.row_upper_[iRow] = iis_row_upper[iisRow];
  }

  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  for (HighsInt iisCol = 0; iisCol < num_iis_col; iisCol++) {
    HighsInt iCol = iis.col_index_[iisCol];
    HighsInt iis_bound = iis.col_bound_[iisCol];
    const double lower = lp.col_lower_[iCol];
    const double upper = lp.col_upper_[iCol];
    double to_lower = lower;
    double to_upper = upper;
    REQUIRE(iis_bound != kIisBoundStatusDropped);
    REQUIRE(iis_bound != kIisBoundStatusNull);
    REQUIRE(iis_bound != kIisBoundStatusBoxed);
    if (iis_bound == kIisBoundStatusLower) {
      to_lower = -inf;
    } else if (iis_bound == kIisBoundStatusUpper) {
      to_upper = inf;
    } else if (iis_bound == kIisBoundStatusFree) {
      if (dev_run)
        printf("IIS Col %2d (LP col %6d) status %s\n", int(iisCol), int(iCol),
               iis.iisBoundStatusToString(iis_bound).c_str());
      continue;
    }
    REQUIRE(highs.changeColBounds(iCol, to_lower, to_upper) ==
            HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    model_status = highs.getModelStatus();
    if (dev_run)
      printf(
          "IIS Col %2d (LP col %6d) status %s: removal yields model status "
          "%s\n",
          int(iisCol), int(iCol), iis.iisBoundStatusToString(iis_bound).c_str(),
          highs.modelStatusToString(model_status).c_str());
    REQUIRE(model_status == HighsModelStatus::kOptimal);
    REQUIRE(highs.changeColBounds(iCol, lower, upper) == HighsStatus::kOk);
  }
  for (HighsInt iisRow = 0; iisRow < num_iis_row; iisRow++) {
    HighsInt iRow = iis.row_index_[iisRow];
    HighsInt iis_bound = iis.row_bound_[iisRow];
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
    double to_lower = lower;
    double to_upper = upper;
    REQUIRE(iis_bound != kIisBoundStatusDropped);
    REQUIRE(iis_bound != kIisBoundStatusNull);
    REQUIRE(iis_bound != kIisBoundStatusFree);
    REQUIRE(iis_bound != kIisBoundStatusBoxed);
    if (iis_bound == kIisBoundStatusLower) {
      to_lower = -inf;
    } else if (iis_bound == kIisBoundStatusUpper) {
      to_upper = inf;
    }
    REQUIRE(highs.changeRowBounds(iRow, to_lower, to_upper) ==
            HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    model_status = highs.getModelStatus();
    if (dev_run)
      printf(
          "IIS Row %2d (LP row %6d) status %s: removal yields model status "
          "%s\n",
          int(iisRow), int(iRow), iis.iisBoundStatusToString(iis_bound).c_str(),
          highs.modelStatusToString(model_status).c_str());
    //    if (model_status != HighsModelStatus::kOptimal)
    //    highs.writeSolution("", kSolutionStylePretty);
    REQUIRE(model_status == HighsModelStatus::kOptimal);
    REQUIRE(highs.changeRowBounds(iRow, lower, upper) == HighsStatus::kOk);
  }

  highs.resetGlobalScheduler(true);
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
    testIis(model, iis);
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

void checkIisLp(HighsLp& lp, const HighsIis& iis, const HighsLp& iis_lp) {
  HighsInt iis_num_col = iis.col_index_.size();
  HighsInt iis_num_row = iis.row_index_.size();
  REQUIRE(iis_lp.num_col_ == iis_num_col);
  REQUIRE(iis_lp.num_row_ == iis_num_row);

  lp.a_matrix_.ensureColwise();
  std::vector<HighsInt> iis_row;
  iis_row.assign(lp.num_row_, -1);
  for (HighsInt iisRow = 0; iisRow < iis_num_row; iisRow++) {
    HighsInt iRow = iis.row_index_[iisRow];
    iis_row[iRow] = iisRow;
    REQUIRE(iis_lp.row_lower_[iisRow] == lp.row_lower_[iRow]);
    REQUIRE(iis_lp.row_upper_[iisRow] == lp.row_upper_[iRow]);
  }
    
  for (HighsInt iisCol = 0; iisCol < iis_num_col; iisCol++) {
    HighsInt iCol = iis.col_index_[iisCol];
    REQUIRE(iis_lp.col_cost_[iisCol] == lp.col_cost_[iCol]);
    REQUIRE(iis_lp.col_lower_[iisCol] == lp.col_lower_[iCol]);
    REQUIRE(iis_lp.col_upper_[iisCol] == lp.col_upper_[iCol]);
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
	 iEl < lp.a_matrix_.start_[iCol+1]; iisCol++) {
      HighsInt iRow = lp.a_matrix_.index_[iEl];
      HighsInt iisRow = iis_row[iRow];
      if (iisRow >= 0) {
	REQUIRE(iis_lp.a_matrix_.index_[iisCol] == iisRow);
	REQUIRE(iis_lp.a_matrix_.value_[iisCol] == lp.a_matrix_.value_[iEl]);
      }
    }
  }

}

