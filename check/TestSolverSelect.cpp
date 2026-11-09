#include <cmath>
#include <string>
#include <utility>
#include <vector>

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

// Every value below was worked out by hand from exampleLp()
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
  // Valid range is [0, 6]: kSimplexOnly, kColsRowsNnz,
  // kRankedFeatures, kHeuristic, kDualSimplex, kIpx, kHipo
  for (const HighsInt value : {0, 1, 2, 3, 4, 5, 6})
    REQUIRE(highs.setOptionValue("solver_select_strategy", value) ==
            HighsStatus::kOk);
  for (const HighsInt value : {-1, 7})
    REQUIRE(highs.setOptionValue("solver_select_strategy", value) ==
            HighsStatus::kError);
}

// Every strategy x require_basis combination, via the real end-to-end
// solver="choose" dispatch in HighsSolve.cpp, solves exampleLp() optimally
TEST_CASE("solver-select-each-strategy-solves-optimally", "[solver_select]") {
  for (const HighsInt strategy : {0, 1, 2, 3, 4, 5, 6}) {
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

// Strategies 4, 5 and 6 choose their solver whatever the LP, even one
// smaller than strategies 1 and 2 send to dual simplex
TEST_CASE("solver-select-forced-strategies", "[solver_select]") {
  const std::vector<std::pair<HighsInt, HighsSolverSelect>> forced = {
      {4, HighsSolverSelect::kDualSimplex},
      {5, HighsSolverSelect::kIpx},
      {6, HighsSolverSelect::kHipo}};
  for (const auto& strategy_solver : forced) {
    const HighsInt strategy = strategy_solver.first;
    const HighsSolverSelect solver = strategy_solver.second;
    INFO("strategy " << strategy);
    REQUIRE(selectSolver(exampleLp(), static_cast<SolverSelectStrategy>(
                                          strategy)) == solver);
    Highs highs;
    highs.setOptionValue("output_flag", dev_run);
    highs.passModel(exampleLp());
    highs.setOptionValue("solver", "choose");
    highs.setOptionValue("solver_select_strategy", strategy);
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
    const bool ipm_ran = highs.getInfo().ipm_iteration_count > 0;
    REQUIRE(ipm_ran == (solver != HighsSolverSelect::kDualSimplex));
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

// Highs::run skips presolve when there is a basis because it assumes
// that solver="choose" will use it, so selection must not discard it
TEST_CASE("solver-select-keeps-warm-start", "[solver_select]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(exampleLp());
  highs.setOptionValue("solver", "choose");
  highs.setOptionValue("presolve", kHighsOffString);
  // exampleLp() is dense, so strategy 3 chooses IPM
  highs.setOptionValue("solver_select_strategy", 3);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getInfo().ipm_iteration_count > 0);
  REQUIRE(highs.getBasis().valid);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(highs.getInfo().ipm_iteration_count == 0);
  REQUIRE(highs.getInfo().simplex_iteration_count == 0);
}

TEST_CASE("solver-select-basis-needs-crossover", "[solver_select]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(exampleLp());
  highs.setOptionValue("solver", "choose");
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("solver_select_strategy", 3);
  highs.setOptionValue("run_crossover", kHighsOffString);

  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(!highs.getBasis().valid);

  highs.clearSolver();
  highs.setOptionValue("solver_select_require_basis", true);
  REQUIRE(highs.run() == HighsStatus::kError);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kSolveError);
}

TEST_CASE("solver-select-count-integer-cols", "[solver_select]") {
  const std::vector<HighsVarType> integrality = {
      HighsVarType::kContinuous, HighsVarType::kInteger,
      HighsVarType::kSemiContinuous, HighsVarType::kSemiInteger,
      HighsVarType::kImplicitInteger};
  // Implicit integers count as continuous
  REQUIRE(countIntegerCols(integrality) == 3);
  REQUIRE(countIntegerCols({}) == 0);
}

TEST_CASE("solver-select-small-model-uses-dual-simplex", "[solver_select]") {
  HighsLpFeatures f;
  f.num_row = 10;
  f.num_col = 10;
  f.num_nz = 20;
  for (const SolverSelectStrategy strategy :
       {SolverSelectStrategy::kColsRowsNnz,
        SolverSelectStrategy::kRankedFeatures})
    REQUIRE(selectSolverByModel(f, strategy) ==
            HighsSolverSelect::kDualSimplex);
}

TEST_CASE("solver-select-heuristic", "[solver_select]") {
  HighsLpFeatures f;
  // Empty
  REQUIRE(selectSolverByHeuristic(f) == HighsSolverSelect::kDualSimplex);
  // Small and sparse
  f.num_row = 1000;
  f.num_col = 1000;
  f.num_nz = 5000;
  f.avg_nonzeros_per_row = 5;
  REQUIRE(selectSolverByHeuristic(f) == HighsSolverSelect::kDualSimplex);
  // Small and dense
  f.avg_nonzeros_per_row = 100;
  REQUIRE(selectSolverByHeuristic(f) == HighsSolverSelect::kIpx);
  // Very large and well conditioned
  f.num_nz = 6000000;
  REQUIRE(selectSolverByHeuristic(f) == HighsSolverSelect::kHipo);
  // Very large and ill conditioned
  f.relative_max_matrix_coefficient = 1e9;
  REQUIRE(selectSolverByHeuristic(f) == HighsSolverSelect::kIpx);
}

// With strategy 3 the root LP of these MIPs is dense enough to choose
// IPM, so this runs IPM + crossover and then the whole tree from that
// basis, including the sub-MIPs, which make their own selections
TEST_CASE("solver-select-mip-tree", "[solver_select]") {
  const std::vector<std::pair<std::string, double>> mips = {
      {"flugpl", 1201500}, {"egout", 568.1007}, {"bell5", 8966406.49}};
  for (const auto& mip : mips) {
    for (const HighsInt strategy : {0, 1, 2, 3, 4, 5, 6}) {
      INFO(mip.first << ", strategy " << strategy);
      const std::string model_file =
          std::string(HIGHS_DIR) + "/check/instances/" + mip.first + ".mps";
      Highs highs;
      highs.setOptionValue("output_flag", dev_run);
      REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
      highs.setOptionValue("solver_select_strategy", strategy);
      REQUIRE(highs.run() == HighsStatus::kOk);
      REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
      const HighsInfo& info = highs.getInfo();
      REQUIRE(std::fabs(info.objective_function_value - mip.second) <
              1e-4 * mip.second);
      REQUIRE(info.max_integrality_violation <= 1e-6);
      REQUIRE(info.max_primal_infeasibility <= 1e-6);
    }
  }
}

// The fitted models index the feature vector, so a change of order
// needs them to be refitted
TEST_CASE("solver-select-feature-order", "[solver_select]") {
  const std::vector<std::string> expected = {
      "num_row",
      "num_col",
      "num_nz",
      "num_integer_col",
      "relative_num_equalities",
      "relative_num_inequalities",
      "relative_num_ranged_rows",
      "relative_num_free_rows",
      "relative_num_cols_without_upper",
      "relative_num_cols_without_lower",
      "relative_num_free_cols",
      "relative_num_boxed_cols",
      "relative_num_singly_bounded_cols",
      "relative_num_fixed_cols",
      "relative_num_almost_identical_nonzeros",
      "relative_num_identical_obj_values",
      "relative_num_identical_rhs_values",
      "relative_max_matrix_coefficient",
      "relative_max_obj_coefficient",
      "relative_max_rhs_value",
      "max_abs_matrix_coefficient",
      "min_abs_matrix_coefficient",
      "max_abs_obj_coefficient",
      "min_abs_obj_coefficient",
      "max_abs_rhs_value",
      "min_abs_rhs_value",
      "matrix_nonzero_density",
      "avg_nonzeros_per_col",
      "avg_nonzeros_per_row",
      "relative_num_dense_rows",
      "max_nonzeros_in_row",
      "max_nonzeros_in_col"};
  REQUIRE(highsLpFeatureNames() == expected);
}

// Without presolve, greenbea is above the size threshold, so this walks
// each of the fitted trees
TEST_CASE("solver-select-fitted-models-solve", "[solver_select]") {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/greenbea.mps";
  for (const HighsInt strategy : {1, 2}) {
    for (const bool require_basis : {false, true}) {
      INFO("strategy " << strategy << ", require_basis " << require_basis);
      Highs highs;
      highs.setOptionValue("output_flag", dev_run);
      REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
      highs.setOptionValue("presolve", kHighsOffString);
      highs.setOptionValue("solver_select_strategy", strategy);
      highs.setOptionValue("solver_select_require_basis", require_basis);
      REQUIRE(highs.run() == HighsStatus::kOk);
      REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
      REQUIRE(std::fabs(highs.getInfo().objective_function_value +
                        7.2555248130e7) < 1);
    }
  }
}

// Presolve runs before selection, so must allow for IPM without
// crossover needing dual postsolve
TEST_CASE("solver-select-ipm-without-crossover-presolve", "[solver_select]") {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/etamacro.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  highs.setOptionValue("solver_select_strategy", 3);
  highs.setOptionValue("run_crossover", kHighsOffString);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
}

// IPM can end optimal without crossover having run, which leaves no
// basis. If one is required, simplex must then find it, otherwise the
// MIP solver would select IPM again for each LP
TEST_CASE("solver-select-basis-without-crossover", "[solver_select]") {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/etamacro.mps";
  for (const bool require_basis : {false, true}) {
    INFO("require_basis " << require_basis);
    Highs highs;
    highs.setOptionValue("output_flag", dev_run);
    REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
    highs.setOptionValue("solver_select_strategy", 3);
    // etamacro is dense enough for strategy 3 to choose IPM
    highs.setOptionValue("run_crossover", kHighsChooseString);
    highs.setOptionValue("solver_select_require_basis", require_basis);
    REQUIRE(highs.run() != HighsStatus::kError);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
    REQUIRE(highs.getInfo().ipm_iteration_count > 0);
    if (require_basis) REQUIRE(highs.getBasis().valid);
  }
}
