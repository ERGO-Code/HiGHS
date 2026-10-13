#include <cmath>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "util/HighsSolverSelect.h"

namespace {
const bool dev_run = false;

// minimize    f  =  x0 +  x1
// subject to              x1 <= 7
//             5 <=  x0 + 2x1 <= 15
//             6 <= 3x0 + 2x1
//             0 <= x0 <= 4; 1 <= x1
//
// Optimum: x0 = 0.5, x1 = 2.25, f = 2.75 (vertex where the second and
// third constraints are both tight).
HighsLp exampleLp() {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {1, 1};
  lp.col_lower_ = {0, 1};
  lp.col_upper_ = {4, kHighsInf};
  lp.row_lower_ = {-kHighsInf, 5, 6};
  lp.row_upper_ = {7, 15, kHighsInf};
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.start_ = {0, 2, 5};
  lp.a_matrix_.index_ = {1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {1, 3, 1, 2, 2};
  return lp;
}

}  // namespace

// Every value below was worked out by hand from exampleLp() (see the
// comments on each row/column), independently of computeLpFeatures's own
// logic, so this is a real check of the arithmetic - not just a
// change-detector pinned to whatever the code currently returns.
TEST_CASE("solver-select-features", "[solver_select]") {
  const HighsLp lp = exampleLp();
  const HighsLpFeatures f = computeLpFeatures(lp);

  REQUIRE(f.num_row == 3);
  REQUIRE(f.num_col == 2);
  REQUIRE(f.num_nz == 5);
  REQUIRE(f.num_integer_col == 0);

  // Rows: x1<=7 (inequality), 5<=x0+2x1<=15 (ranged), 3x0+2x1>=6 (inequality)
  REQUIRE(f.relative_num_equalities == 0.0);
  REQUIRE(std::fabs(f.relative_num_inequalities - 2.0 / 3.0) < 1e-12);
  REQUIRE(std::fabs(f.relative_num_ranged_rows - 1.0 / 3.0) < 1e-12);
  REQUIRE(f.relative_num_free_rows == 0.0);

  // Columns: x0 in [0,4] (boxed), x1 in [1,inf) (singly bounded)
  REQUIRE(f.relative_num_cols_without_upper == 0.5);
  REQUIRE(f.relative_num_cols_without_lower == 0.0);
  REQUIRE(f.relative_num_free_cols == 0.0);
  REQUIRE(f.relative_num_boxed_cols == 0.5);
  REQUIRE(f.relative_num_singly_bounded_cols == 0.5);
  REQUIRE(f.relative_num_fixed_cols == 0.0);

  // |a_ij| = {1,1,2,2,3}: the two 1s and two 2s cluster, the 3 doesn't
  REQUIRE(f.relative_num_almost_identical_nonzeros == 0.8);
  // c = {1,1}: identical
  REQUIRE(f.relative_num_identical_obj_values == 1.0);
  // representative rhs per row = {7,15,6}: all distinct
  REQUIRE(f.relative_num_identical_rhs_values == 0.0);

  REQUIRE(f.max_abs_matrix_coefficient == 3.0);
  REQUIRE(f.min_abs_matrix_coefficient == 1.0);
  REQUIRE(f.relative_max_matrix_coefficient == 3.0);
  REQUIRE(f.max_abs_obj_coefficient == 1.0);
  REQUIRE(f.min_abs_obj_coefficient == 1.0);
  REQUIRE(f.relative_max_obj_coefficient == 1.0);
  REQUIRE(f.max_abs_rhs_value == 15.0);
  REQUIRE(f.min_abs_rhs_value == 5.0);
  REQUIRE(f.relative_max_rhs_value == 3.0);

  REQUIRE(std::fabs(f.matrix_nonzero_density - 5.0 / 6.0) < 1e-12);
  REQUIRE(f.avg_nonzeros_per_col == 2.5);
  REQUIRE(std::fabs(f.avg_nonzeros_per_row - 5.0 / 3.0) < 1e-12);
  REQUIRE(f.relative_num_dense_rows == 0.0);
  REQUIRE(f.max_nonzeros_in_row == 2.0);
  REQUIRE(f.max_nonzeros_in_col == 3.0);
}

TEST_CASE("solver-select-strategy-option-bounds", "[solver_select]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  // Valid range is [0, 3]: kSimplexOnly, kNColsRowsNnzIntegrality, kFull,
  // kHeuristic
  for (const HighsInt value : {0, 1, 2, 3})
    REQUIRE(highs.setOptionValue("solver_select_strategy", value) ==
            HighsStatus::kOk);
  for (const HighsInt value : {-1, 4})
    REQUIRE(highs.setOptionValue("solver_select_strategy", value) ==
            HighsStatus::kError);
}

TEST_CASE("solver-select-require-basis-option", "[solver_select]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.setOptionValue("solver_select_require_basis", true) ==
          HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("solver_select_require_basis", false) ==
          HighsStatus::kOk);
}

// Every strategy x require_basis combination, via the real end-to-end
// solver="choose" dispatch in HighsSolve.cpp, solves exampleLp() optimally
TEST_CASE("solver-select-each-strategy-solves-optimally", "[solver_select]") {
  for (const HighsInt strategy : {0, 1, 2, 3}) {
    for (const bool require_basis : {false, true}) {
      INFO("strategy " << strategy << ", require_basis " << require_basis);
      Highs highs;
      highs.setOptionValue("output_flag", dev_run);
      highs.passModel(exampleLp());
      highs.setOptionValue("solver", "choose");
      highs.setOptionValue("solver_select_strategy", strategy);
      highs.setOptionValue("solver_select_require_basis", require_basis);

      REQUIRE(highs.run() == HighsStatus::kOk);
      REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
      REQUIRE(std::fabs(highs.getInfo().objective_function_value - 2.75) <
              1e-6);
    }
  }
}

TEST_CASE("solver-select-require-basis-yields-a-basis", "[solver_select]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(exampleLp());
  highs.setOptionValue("solver", "choose");
  highs.setOptionValue("solver_select_strategy", 3);
  highs.setOptionValue("solver_select_require_basis", true);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(highs.getBasis().valid);
}
