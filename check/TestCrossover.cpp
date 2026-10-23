#include "HCheckConfig.h"
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

  highs.resetGlobalScheduler(true);
}

void report(const std::string message, const HighsLp& lp,
            const HighsSolution& solution, const HighsBasis& basis) {
  if (!dev_run) return;
  double objective = lp.offset_;
  printf("Solution: %s\n  Ix       Value   Status\n", message.c_str());
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    const double value = solution.col_value[iCol];
    const std::string status =
        basis.valid ? utilBasisStatusToString(basis.col_status[iCol]).full_
                    : "";
    objective += value * lp.col_cost_[iCol];
    printf("%4d %11.4g   %s\n", (int)iCol, value, status.c_str());
  }
  printf(" Obj %11.4g\n", objective);
}

const std::vector<std::string> solvers{kPdlpString, kHiPdlpString, kIpxString
#ifdef HIPO
                                       ,
                                       kHipoString
#endif
};

TEST_CASE("issue-3459", "[highs_lp_solver]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  const bool use_3459_model = true;
  if (use_3459_model) {
    std::string model_file =
        std::string(HIGHS_DIR) + "/check/instances/3459.mps";
    h.readModel(model_file);
  } else {
    HighsLp lp;
    lp.num_col_ = 2;
    lp.num_row_ = 3;
    lp.col_cost_ = {1, -1};
    lp.col_lower_ = {0, 0};
    lp.col_upper_ = {kHighsInf, kHighsInf};
    lp.row_lower_ = {1, 2, 3};  // Infeasible or unbounded
    lp.row_upper_ = {2, 3, 4};  // Infeasible or unbounded
    lp.row_lower_ = {2, 0, 0};
    lp.row_upper_ = {3, 1, 1};
    lp.a_matrix_.start_ = {0, 3, 6};
    lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
    lp.a_matrix_.value_ = {1, -1, 2, 2, 1, -1};
    h.passModel(lp);
  }

  HighsStatus return_status;
  HighsSolution solution;
  h.setOptionValue(kPresolveString, kHighsOffString);
  h.setOptionValue("run_crossover", kHighsOffString);
  const double objective_value = 17.555148585;

  auto testSolvers = [&](const char* type) {
    for (auto& solver : solvers) {
      if (dev_run)
        printf("====================\n%s (%s)\n====================\n",
               solver.c_str(), type);
      h.setOptionValue(kSolverString, solver);
      return_status = h.run();
      REQUIRE(return_status == HighsStatus::kOk);
      REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
      REQUIRE(std::fabs(h.getObjectiveValue() - objective_value) < 1e-6);

      solution = h.getSolution();

      return_status = h.crossover(solution);
      REQUIRE(return_status == HighsStatus::kOk);
      REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
      REQUIRE(std::fabs(h.getObjectiveValue() - objective_value) < 1e-6);
    }
  };

  testSolvers("original");

  HighsInt free_row0 = 10;
  HighsInt free_row1 = 14;
  REQUIRE(std::fabs(h.getSolution().row_dual[free_row0]) <
          h.getOptions().dual_feasibility_tolerance);
  REQUIRE(std::fabs(h.getSolution().row_dual[free_row1]) <
          h.getOptions().dual_feasibility_tolerance);

  h.changeRowBounds(free_row0, -kHighsInf, kHighsInf);
  h.changeRowBounds(free_row1, -kHighsInf, kHighsInf);

  testSolvers("with two free rows");

  h.resetGlobalScheduler(true);
}
