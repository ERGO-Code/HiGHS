#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsModelUtils.h"

const bool dev_run = false;
const double inf = kHighsInf;

void test_crossover(Highs& highs, HighsLp& lp);
void report(const std::string message, const HighsLp& lp,
            const HighsSolution& solution, const HighsBasis& basis);

// No commas in test case name.
TEST_CASE("test-crossover", "[highs_crossover]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMaximize;
  lp.col_cost_ = {1, 2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {4};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1, 2};
  HighsLp df_lp = lp;

  Highs highs;
  // From feasible non-vertex for feasible LP
  if (dev_run) printf("\nCrossover from feasible non-vertex for feasible LP\n");
  test_crossover(highs, lp);

  if (dev_run)
    printf("\nCrossover from infeasible non-vertex for feasible LP\n");
  lp = df_lp;
  lp.col_lower_[0] = 4;
  test_crossover(highs, lp);

  if (dev_run)
    printf("\nCrossover from infeasible non-vertex for infeasible LP\n");
  lp = df_lp;
  lp.col_lower_ = {4, 2};
  test_crossover(highs, lp);

  if (dev_run)
    printf("\nCrossover from feasible non-vertex for unbounded LP\n");
  lp = df_lp;
  lp.row_lower_ = {4};
  lp.row_upper_ = {inf};
  test_crossover(highs, lp);
}

void test_crossover(Highs& highs, HighsLp& lp) {
  HighsStatus return_status;
  highs.clear();
  highs.setOptionValue("output_flag", dev_run);
  //  highs.setOptionValue("output_flag", false);
  return_status = highs.passModel(lp);
  REQUIRE(return_status == HighsStatus::kOk);
  highs.run();
  const std::vector<double> from_primal_value = {2, 1};
  HighsModelStatus model_status = highs.getModelStatus();
  const HighsInfo& info = highs.getInfo();
  const double require_optimal_objective = info.objective_function_value;

  highs.clearSolver();
  // Set solution to interior of optimal face
  HighsSolution solution;
  solution.col_value = from_primal_value;
  HighsBasis basis;
  report("Before crossover", lp, solution, basis);
  return_status = highs.crossover(solution);
  if (model_status == HighsModelStatus::kOptimal) {
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == model_status);
    REQUIRE(require_optimal_objective == info.objective_function_value);
  } else if (model_status == HighsModelStatus::kInfeasible) {
    // Crossover returns "imprecise" => HighsStatus::kWarning and
    // HighsModelStatus::kUnknown
    REQUIRE(return_status == HighsStatus::kWarning);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);
  } else if (model_status == HighsModelStatus::kUnbounded) {
    // Crossover returns "imprecise" => HighsStatus::kWarning and
    // HighsModelStatus::kUnknown
    REQUIRE(return_status == HighsStatus::kWarning);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);
  } else {
    assert(1 == 0);
  }
  report("After crossover", lp, solution, basis);
  report("From Highs", lp, highs.getSolution(), highs.getBasis());
}
void report(const std::string message, const HighsLp& lp,
            const HighsSolution& solution, const HighsBasis& basis) {
  if (!dev_run) return;
  double objective = lp.offset_;
  printf("Solution: %s\n  Ix       Value   Status\n", message.c_str());
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const double value = solution.col_value[iCol];
    const std::string status =
        basis.valid ? utilBasisStatusToString(basis.col_status[iCol]) : "";
    objective += value * lp.col_cost_[iCol];
    printf("%4d %11.4g   %s\n", (int)iCol, value, status.c_str());
  }
  printf(" Obj %11.4g\n", objective);
}
