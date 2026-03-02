#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;  // true;//
const bool write_model = false;

const double inf = kHighsInf;

const HighsInt kIisStrategyFromRayColPriority =
    kIisStrategyFromRay + kIisStrategyColPriority;
const HighsInt kIisStrategyFromLpColPriority =
    kIisStrategyFromLp + kIisStrategyColPriority;

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
  REQUIRE(iis.valid_ == true);
  REQUIRE(iis.status_ == kIisModelStatusIrreducible);

  // Perform full IIS
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);

  highs.setOptionValue("iis_strategy", kIisStrategyFromLp);
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(iis.valid_ == true);
  REQUIRE(iis.status_ == kIisModelStatusIrreducible);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == 0);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusBoxed);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    REQUIRE(iis.col_status_[iCol] == kIisStatusNotInConflict);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (iRow == iis.row_index_[0]) {
      REQUIRE(iis.row_status_[iRow] == kIisStatusInConflict);
    } else {
      REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);
    }
  }

  highs.setOptionValue("iis_strategy", kIisStrategyFromLpColPriority);
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(iis.valid_ == true);
  REQUIRE(iis.status_ == kIisModelStatusIrreducible);
  REQUIRE(iis.col_index_.size() == 1);
  REQUIRE(iis.row_index_.size() == 0);
  REQUIRE(iis.col_index_[0] == 2);
  REQUIRE(iis.col_bound_[0] == kIisBoundStatusBoxed);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (iCol == iis.col_index_[0]) {
      REQUIRE(iis.col_status_[iCol] == kIisStatusInConflict);
    } else {
      REQUIRE(iis.col_status_[iCol] == kIisStatusNotInConflict);
    }
  }
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
    REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);

  // Now give two columns and two rows incompatible bounds, and ensure
  // that just one of each (the first encountered) is found
  lp.col_upper_[0] = -1;
  lp.row_upper_[1] = -1;
  const bool two_inconsistent_rows = lp.row_upper_[0] < lp.row_lower_[0] &&
                                     lp.row_upper_[1] < lp.row_lower_[1];
  const bool two_inconsistent_cols = lp.col_upper_[0] < lp.col_lower_[0] &&
                                     lp.col_upper_[2] < lp.col_lower_[2];
  REQUIRE(two_inconsistent_cols);
  REQUIRE(two_inconsistent_rows);

  highs.passModel(lp);
  highs.setOptionValue("iis_strategy", kIisStrategyFromLp);
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(iis.valid_ == true);
  REQUIRE(iis.status_ == kIisModelStatusIrreducible);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == 0);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusBoxed);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    REQUIRE(iis.col_status_[iCol] == kIisStatusNotInConflict);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (iRow == iis.row_index_[0]) {
      REQUIRE(iis.row_status_[iRow] == kIisStatusInConflict);
    } else {
      REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);
    }
  }

  highs.setOptionValue("iis_strategy", kIisStrategyFromLpColPriority);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
  REQUIRE(iis.valid_ == true);
  REQUIRE(iis.status_ == kIisModelStatusIrreducible);
  REQUIRE(iis.col_index_.size() == 1);
  REQUIRE(iis.row_index_.size() == 0);
  REQUIRE(iis.col_index_[0] == 0);
  REQUIRE(iis.col_bound_[0] == kIisBoundStatusBoxed);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (iCol == iis.col_index_[0]) {
      REQUIRE(iis.col_status_[iCol] == kIisStatusInConflict);
    } else {
      REQUIRE(iis.col_status_[iCol] == kIisStatusNotInConflict);
    }
  }
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
    REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);

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
  REQUIRE(iis.valid_ == true);
  REQUIRE(iis.status_ == kIisModelStatusIrreducible);
  if (dev_run && write_model) {
    highs.writeModel("");
    highs.writeIisModel("");
  }
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == empty_row);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusLower);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    REQUIRE(iis.col_status_[iCol] == kIisStatusNotInConflict);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (iRow == iis.row_index_[0]) {
      REQUIRE(iis.row_status_[iRow] == kIisStatusInConflict);
    } else {
      REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);
    }
  }

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
  REQUIRE(iis.valid_ == true);
  REQUIRE(iis.status_ == kIisModelStatusIrreducible);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 1);
  REQUIRE(iis.row_index_[0] == empty_row);
  REQUIRE(iis.row_bound_[0] == kIisBoundStatusUpper);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    REQUIRE(iis.col_status_[iCol] == kIisStatusNotInConflict);
  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
    if (iRow == iis.row_index_[0]) {
      REQUIRE(iis.row_status_[iRow] == kIisStatusInConflict);
    } else {
      REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);
    }
  }

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
  //      1.5w + 2x +   y      <= 30
  //
  // -10 <= 4w - 2x +   y + 2z <= 15
  //
  // -34 <=    - 2x -1.5y -  z
  //
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  HighsIis iis;

  for (int l = 0; l < 3; l++) {
    for (int k = 0; k < 2; k++) {
      REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
      REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
      REQUIRE(iis.valid_ == true);
      REQUIRE(iis.status_ == kIisModelStatusIrreducible);
      REQUIRE(iis.row_index_.size() == 1);
      HighsInt iis_row = iis.row_index_[0];
      if (lp.a_matrix_.isColwise()) {
        for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
          for (HighsInt iEl = lp.a_matrix_.start_[iCol];
               iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
            if (lp.a_matrix_.index_[iEl] == iis_row) {
              REQUIRE(iis.col_status_[iCol] == kIisStatusInConflict);
              break;
            }
          }
        }
      } else {
        for (HighsInt iEl = lp.a_matrix_.start_[iis_row];
             iEl < lp.a_matrix_.start_[iis_row + 1]; iEl++) {
          HighsInt iCol = lp.a_matrix_.index_[iEl];
          REQUIRE(iis.col_status_[iCol] == kIisStatusInConflict);
        }
      }
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
        if (iRow == iis.row_index_[0]) {
          REQUIRE(iis.row_status_[iRow] == kIisStatusInConflict);
        } else {
          REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);
        }
      }
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
  const HighsLp& highs_lp = highs.getLp();
  // First pass with incumbent matrix colwise; second with it
  // rowwise
  highs.ensureColwise();
  REQUIRE(highs_lp.a_matrix_.isColwise());
  for (HighsInt k = 0; k < 2; k++) {
    REQUIRE(highs.getIis(iis) == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
    REQUIRE(iis.valid_ == true);
    REQUIRE(iis.status_ == kIisModelStatusIrreducible);
    REQUIRE(iis.col_index_.size() == 2);
    REQUIRE(iis.row_index_.size() == 1);
    REQUIRE(iis.col_index_[0] == 0);
    REQUIRE(iis.col_index_[1] == 1);
    REQUIRE(iis.row_index_[0] == 2);

    HighsInt iis_row = iis.row_index_[0];
    if (lp.a_matrix_.isColwise()) {
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
        for (HighsInt iEl = lp.a_matrix_.start_[iCol];
             iEl < lp.a_matrix_.start_[iCol + 1]; iEl++) {
          if (lp.a_matrix_.index_[iEl] == iis_row) {
            REQUIRE(iis.col_status_[iCol] == kIisStatusInConflict);
            break;
          }
        }
      }
    } else {
      for (HighsInt iEl = lp.a_matrix_.start_[iis_row];
           iEl < lp.a_matrix_.start_[iis_row + 1]; iEl++) {
        HighsInt iCol = lp.a_matrix_.index_[iEl];
        REQUIRE(iis.col_status_[iCol] == kIisStatusInConflict);
      }
    }
    for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++) {
      if (iRow == iis.row_index_[0]) {
        REQUIRE(iis.row_status_[iRow] == kIisStatusInConflict);
      } else {
        REQUIRE(iis.row_status_[iRow] == kIisStatusNotInConflict);
      }
    }
    highs.clearSolver();
    highs.ensureRowwise();
    REQUIRE(highs_lp.a_matrix_.isRowwise());
  }
  highs.resetGlobalScheduler(true);
}

TEST_CASE("lp-get-iis-woodinfe", "[iis]") {
  std::string model = "woodinfe";
  testMps(model, kIisStrategyLight);
  testMps(model, kIisStrategyFromLp);
  //  testMps(model, kIisStrategyFromRay);
  //
  // No need for a +kIisStrategyIrreducible test, since kIisStrategyFromLp
  // yields IIS
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
  // This is infeasible since c4 >= 30 and c4 <= 30 fixes c4 = 30,
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
  //
  // If the elasticity filter is used, then it identifies the
  // following infeasibility system
  //
  // r4:  0 <= c2 + c3      - c6 - c7
  //
  // r6: 20 <=           c5 + c6
  //
  // r7: 30 <=                     c7
  //
  // This is infeasible since c7 >= 30 gives 30 <= c2 + c3, but c2 and
  // c3 have upper bounds of 10
  //
  // Hence the IIS does not require r6 or c5, and consists of r4 and r7
  // (>=0) with c2 <= 10; c3 <= 10; c6 free; c7 free

  std::string model = "galenet";
  testMps(model, kIisStrategyLight, HighsModelStatus::kNotset);
  testMps(model, kIisStrategyFromLp);
  testMps(model, kIisStrategyFromLp + kIisStrategyIrreducible);
}

TEST_CASE("lp-get-iis-avgas", "[iis]") {
  std::string model = "avgas";
  // For the whole LP calculation the elasticity filter only
  // identified feasibility, so the model status is not set
  testMps(model, kIisStrategyFromLp, HighsModelStatus::kNotset);
  // For the ray calculation the model is solved, so its status is
  // known
  //  testMps(model, kIisStrategyFromRay,
  //  HighsModelStatus::kOptimal);
}

TEST_CASE("lp-feasibility-relaxation", "[iis]") {
  // Using infeasible MIP from AMPL documentation
  //
  // https://mp.ampl.com/features-guide.html#feasibility-relaxation
  //
  // min x-2y
  //
  // 2 <= -x + 21y
  //
  //     -3x +  2y <= 1
  //
  //     20x +   y <= 20
  //
  //  5 <= x;    y free
  //
  //  Vanilla feasibility relaxation
  //
  // min e0 + e1 + e2 + e3
  //
  // 2 <= e1 - x + 21y
  //
  //     -3x +  2y - e2 <= 1
  //
  //     20x +   y - e3 <= 20
  //
  //  5 <= e0 + x
  //
  //  x, y free
  //
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

  const bool all_tests = true;
  const bool test0 = true || all_tests;
  const bool test1 = true || all_tests;
  const bool test2 = true || all_tests;
  const bool test3 = true || all_tests;
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
    // Now test local penalties, allowing lower bounds to be violated,
    // but not upper bounds
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
    // AMPL says: Should get slacks (-3, 4, 0) corresponding to (x, y)
    // = (1, 0), giving objective = 4 + 3 = 7
    //
    // However, (x, y) = (0, 0) gives slacks (-2, 1, 20) and also
    // objective = 5 + 2 = 7
    //
    // The 64-bit integer build (11/11/25) gives (x, y) = (1, 0), and
    // the 32-bit integer build gives (x, y) = (0, 0). Since this may
    // vary randomly in future, the test below is dependent on whether
    // x = 0 or x = 1
    const bool solution0 = solution.col_value[0] == 0;
    double r0_slack = solution0 ? -2 : -3;
    double r1_slack = solution0 ? 1 : 4;
    double r2_slack = solution0 ? 20 : 0;
    h.writeSolution("", 1);
    REQUIRE(solution.row_value[0] == lp.row_lower_[0] + r0_slack);
    REQUIRE(solution.row_value[1] == lp.row_upper_[1] - r1_slack);
    REQUIRE(solution.row_value[2] == lp.row_upper_[2] - r2_slack);
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
    h.writeSolution("", 1);
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
  //  if (iis_strategy == kIisStrategyFromRay ||
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
    REQUIRE(iis.valid_ == true);
    const bool find_irreducible = kIisStrategyIrreducible & iis_strategy;
    if (find_irreducible) REQUIRE(iis.status_ == kIisModelStatusIrreducible);
    const HighsInt in_iis_status = iis.status_ == kIisModelStatusIrreducible
                                       ? kIisStatusInConflict
                                       : kIisStatusMaybeInConflict;
    for (HighsInt iX = 0; iX < num_iis_col; iX++)
      REQUIRE(iis.col_status_[iis.col_index_[iX]] == in_iis_status);
    for (HighsInt iX = 0; iX < num_iis_row; iX++)
      REQUIRE(iis.row_status_[iis.row_index_[iX]] == in_iis_status);
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

TEST_CASE("feasible-lp-iis", "[iis]") {
  HighsLp lp;
  lp.model_name_ = "chip";
  lp.sense_ = ObjSense::kMaximize;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {10, 25};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.col_names_ = {"Tables", "SetsOfChairs"};
  lp.row_lower_ = {-inf, -inf};
  lp.row_upper_ = {80, 120};
  lp.row_names_ = {"Assembly", "Finishng"};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 2, 4};
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.passModel(lp);
  HighsIis iis;
  // With kIisStrategyLight, feasibility of the LP is not determined
  // until it's been solved
  h.getIis(iis);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 0);
  REQUIRE(iis.col_status_[0] == kIisStatusMaybeInConflict);
  REQUIRE(iis.col_status_[1] == kIisStatusMaybeInConflict);
  REQUIRE(iis.row_status_[0] == kIisStatusMaybeInConflict);
  REQUIRE(iis.row_status_[1] == kIisStatusMaybeInConflict);

  h.run();

  h.getIis(iis);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 0);
  REQUIRE(iis.col_status_[0] == kIisStatusNotInConflict);
  REQUIRE(iis.col_status_[1] == kIisStatusNotInConflict);
  REQUIRE(iis.row_status_[0] == kIisStatusNotInConflict);
  REQUIRE(iis.row_status_[1] == kIisStatusNotInConflict);

  h.passModel(lp);
  // With kIisStrategyFromLp, feasibility of the LP is determined
  h.setOptionValue("iis_strategy", kIisStrategyFromLp);

  h.getIis(iis);
  REQUIRE(iis.col_index_.size() == 0);
  REQUIRE(iis.row_index_.size() == 0);
  REQUIRE(iis.col_status_[0] == kIisStatusNotInConflict);
  REQUIRE(iis.col_status_[1] == kIisStatusNotInConflict);
  REQUIRE(iis.row_status_[0] == kIisStatusNotInConflict);
  REQUIRE(iis.row_status_[1] == kIisStatusNotInConflict);

  h.resetGlobalScheduler(true);
}

TEST_CASE("write-iis_model-file", "[iis]") {
  // Reproduces #2635, and adds code coverage for writing IIS model
  // and solution at runtime
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string test_mps = test_name + ".mps";
  const std::string test_sol = test_name + ".sol";
  HighsLp lp;
  lp.model_name_ = "2635";
  lp.num_col_ = 5;
  lp.num_row_ = 4;
  lp.col_cost_ = {-1, -1, 0, 0, 0};
  lp.col_lower_ = {0, 0, 10, 0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf, 10, kHighsInf, kHighsInf};
  lp.row_lower_ = {0, 0, 0, 20};
  lp.row_upper_ = {0, 0, 0, kHighsInf};
  lp.a_matrix_.start_ = {0, 2, 4, 5, 7, 10};
  lp.a_matrix_.index_ = {1, 3, 2, 3, 1, 0, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {-1, 1, -1, 1, 1, 1, 1, -0.5, -1, 1};
  lp.col_names_ = {"x0", "x1", "x2", "x3", "x4"};
  lp.row_names_ = {"x3-.5x4=0", "x2-x4-x0=0", "x3+x4-x1=0", "x0+x1>=20"};
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("log_file", "highs.log");
  h.setOptionValue("write_iis_model_file", test_mps);
  h.setOptionValue("solution_file", test_sol);
  const HighsInt to_k = 2;
  for (HighsInt k = 0; k <= to_k; k++) {
    if (k == 0) {
      h.setOptionValue("iis_strategy", kIisStrategyLight);
    } else if (k == 1) {
      h.setOptionValue("iis_strategy", kIisStrategyFromLp);
    } else {
      assert(k == 2);
      h.setOptionValue("iis_strategy",
                       kIisStrategyFromLp + kIisStrategyIrreducible);
    }
    if (dev_run)
      printf(
          "\nPass %d with iis_strategy = %d\n============================\n\n",
          int(k), int(h.getOptions().iis_strategy));
    h.passModel(lp);
    h.run();
    h.readModel(test_mps);
    if (k == 0) {
      REQUIRE(h.getLp().num_col_ == 0);
      REQUIRE(h.getLp().num_row_ == 0);
    } else {
      // Use h.optimizeModel(); to avoid unnecessary MPS and solution
      // file
      h.optimizeModel();
      REQUIRE(h.getModelStatus() == HighsModelStatus::kInfeasible);
    }
  }
  std::remove(test_mps.c_str());
  std::remove(test_sol.c_str());
  h.resetGlobalScheduler(true);
}
