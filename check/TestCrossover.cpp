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

  highs.resetThreadScheduler(true);
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

const std::vector<std::string> solvers{kPdlpString,
                                       // kHiPdlpString,
                                       kIpxString
#ifdef HIPO
                                       ,
                                       kHipoString
#endif
};

TEST_CASE("issue-3459", "[highs_crossover]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  const bool use_3459_model = false;
  double objective_value;
  HighsInt free_row0;
  HighsInt free_row1;
  if (use_3459_model) {
    std::string model_file =
        std::string(HIGHS_DIR) + "/check/instances/3459.mps";
    h.readModel(model_file);
    objective_value = 17.555148585;
    free_row0 = 10;
    free_row1 = 14;
  } else {
    // This LP has optimal solution
    //
    // Columns
    // Ix St  Lo  Up    Pr  Du  Name
    //  0 UB  20  30    30  -4  Tables
    //  1 LB   5   7     5   2  Seats
    //  2 BS  20  40    25   0  Sets_of_chairs
    //
    // Rows
    // Ix St  Lo  Up     Pr Du  Name
    //  0 BS  70  90  86.25  0  Assembly
    //  1 BS 100 120 116.25  0  Redundant
    //  2 UB 125 145 145    -6  Finishing
    //
    // So, in setting up for crossover there are boxed variables at
    // both lower and upper bounds, and boxed constraints testing the
    // addition of explicit slacks. Also, row 0 and row 1 can be made
    // free, so testing the handling of free rows - which are omitted
    // from IPX
    HighsLp lp;
    lp.sense_ = ObjSense::kMaximize;
    lp.num_col_ = 3;
    lp.num_row_ = 3;
    lp.col_cost_ = {-10, -16, -24};
    lp.col_lower_ = {20, 5, 20};
    lp.col_upper_ = {30, 7, 40};
    lp.row_lower_ = {70, 100, 125};
    lp.row_upper_ = {90, 120, 145};
    lp.a_matrix_.start_ = {0, 3, 6, 9};
    lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2, 0, 1, 2};
    lp.a_matrix_.value_ = {1, 1, 1, 1.25, 2.25, 3, 2, 3, 4};
    lp.col_names_ = {"Tables", "Seats", "Sets_of_chairs"};
    lp.row_names_ = {"Assembly", "Redundant", "Finishing"};
    objective_value = -980;
    if (lp.sense_ == ObjSense::kMaximize) {
      for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
        lp.col_cost_[iCol] *= -1;
      objective_value *= -1;
      ;
    }
    h.passModel(lp);
    free_row0 = 0;
    free_row1 = 1;
  }

  HighsStatus return_status;
  HighsSolution solution;
  h.setOptionValue(kPresolveString, kHighsOffString);
  h.setOptionValue("run_crossover", kHighsOffString);

  auto objectiveError = [&]() {
    return std::fabs(h.getObjectiveValue() - objective_value) /
           (1.0 + std::fabs(objective_value));
  };

  auto testSolvers = [&](const char* type) {
    for (auto& solver : solvers) {
      if (dev_run)
        printf("====================\n%s (%s)\n====================\n",
               solver.c_str(), type);
      h.setOptionValue(kSolverString, solver);
      return_status = h.run();
      REQUIRE(return_status == HighsStatus::kOk);
      REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
      REQUIRE(objectiveError() < 1e-6);
      //      h.writeSolution("", 1);
      solution = h.getSolution();
      return_status = h.crossover(solution);
      REQUIRE(return_status == HighsStatus::kOk);
      REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
      REQUIRE(objectiveError() < 1e-6);
      //      h.writeSolution("", 1);
    }
  };

  testSolvers("original");

  REQUIRE(std::fabs(h.getSolution().row_dual[free_row0]) <
          h.getOptions().dual_feasibility_tolerance);
  REQUIRE(std::fabs(h.getSolution().row_dual[free_row1]) <
          h.getOptions().dual_feasibility_tolerance);

  REQUIRE(h.changeRowBounds(free_row0, -kHighsInf, kHighsInf) ==
          HighsStatus::kOk);
  REQUIRE(h.changeRowBounds(free_row1, -kHighsInf, kHighsInf) ==
          HighsStatus::kOk);

  testSolvers("with free row");

  h.resetGlobalScheduler(true);
}
