#include "HCheckConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;

void presolveSolvePostsolve(const std::string& model_file,
                            const bool solve_relaxation = false);

TEST_CASE("presolve-solve-postsolve-lp", "[highs_test_presolve]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";
  presolveSolvePostsolve(model_file);
}

TEST_CASE("presolve-solve-postsolve-mip", "[highs_test_presolve]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
  presolveSolvePostsolve(model_file);
}

TEST_CASE("presolve-solve-postsolve-relaxation", "[highs_test_presolve]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
  presolveSolvePostsolve(model_file, true);
}

TEST_CASE("presolve", "[highs_test_presolve]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  // Make sure that an empty LP returns kNotReduced
  const HighsModel& presolved_model = highs.getPresolvedModel();
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kNotReduced);

  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  highs.readModel(model_file);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kReduced);
  REQUIRE(!presolved_model.isEmpty());

  model_file = std::string(HIGHS_DIR) + "/check/instances/gas11.mps";
  highs.readModel(model_file);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() ==
          HighsPresolveStatus::kUnboundedOrInfeasible);
  REQUIRE(presolved_model.isEmpty());

  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  SpecialLps special_lps;

  special_lps.scipLpi3Lp(lp, require_model_status);
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kInfeasible);
  REQUIRE(highs.getModelStatus() == require_model_status);
  REQUIRE(presolved_model.isEmpty());

  special_lps.distillationLp(lp, require_model_status, optimal_objective);
  // Have to set matrix dimensions to match presolved_model.lp_
  lp.setMatrixDimensions();
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(lp.equalButForNames(presolved_model.lp_));
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kNotReduced);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kNotset);
  REQUIRE(!presolved_model.isEmpty());

  special_lps.primalDualInfeasible1Lp(lp, require_model_status);
  highs.passModel(lp);
  REQUIRE(highs.presolve() == HighsStatus::kOk);
  REQUIRE(highs.getModelPresolveStatus() == HighsPresolveStatus::kInfeasible);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(presolved_model.isEmpty());
}

TEST_CASE("empty-row", "[highs_test_presolve]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.col_cost_ = {-7.0, -6.0, -5.0};
  lp.col_lower_ = {-73.0, -83.0, -94.0};
  lp.col_upper_ = {62.0, 96.0, 62.0};
  lp.row_lower_ = {-19.0};
  lp.row_upper_ = {11.0};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 0};
  // LP has empty constraint matrix so doesn't need to be presolved,
  // and shouldn't be since this would cause vacuous null pointer
  // operation in util/HighsMatrixSlice.h (see #1531)
  highs.passModel(lp);
  highs.run();
  const HighsSolution& solution = highs.getSolution();
  const HighsBasis& basis = highs.getBasis();
  REQUIRE(HighsInt(solution.row_value.size()) == lp.num_row_);
  REQUIRE(HighsInt(basis.row_status.size()) == lp.num_row_);
}

void presolveSolvePostsolve(const std::string& model_file,
                            const bool solve_relaxation) {
  Highs highs0;
  Highs highs1;
  highs0.setOptionValue("output_flag", dev_run);
  highs1.setOptionValue("output_flag", dev_run);
  HighsStatus return_status;
  highs0.readModel(model_file);
  if (solve_relaxation) highs0.setOptionValue("solver", kSimplexString);
  return_status = highs0.presolve();
  REQUIRE(return_status == HighsStatus::kOk);
  HighsPresolveStatus model_presolve_status = highs0.getModelPresolveStatus();
  if (model_presolve_status == HighsPresolveStatus::kTimeout) {
    if (dev_run)
      printf("Presolve timeout: return status = %d\n", (int)return_status);
  }
  HighsLp lp = highs0.getPresolvedLp();
  highs1.passModel(lp);
  if (solve_relaxation) highs1.setOptionValue("solver", kSimplexString);
  highs1.setOptionValue("presolve", kHighsOffString);
  highs1.run();
  HighsSolution solution = highs1.getSolution();
  const double objective_value = highs1.getInfo().objective_function_value;
  if (lp.isMip() && !solve_relaxation) {
    return_status = highs0.postsolve(solution);
    REQUIRE(return_status == HighsStatus::kWarning);
    HighsModelStatus model_status = highs0.getModelStatus();
    REQUIRE(model_status == HighsModelStatus::kUnknown);
    const double dl_objective_value =
        std::fabs(highs0.getInfo().objective_function_value - objective_value);
    REQUIRE(dl_objective_value < 1e-12);
    REQUIRE(highs0.getInfo().primal_solution_status == kSolutionStatusFeasible);
    double mip_feasibility_tolerance;
    highs0.getOptionValue("mip_feasibility_tolerance",
                          mip_feasibility_tolerance);
    REQUIRE(highs0.getInfo().max_integrality_violation <=
            mip_feasibility_tolerance);
  } else {
    HighsBasis basis = highs1.getBasis();
    return_status = highs0.postsolve(solution, basis);
    REQUIRE(return_status == HighsStatus::kOk);
    HighsModelStatus model_status = highs0.getModelStatus();
    REQUIRE(model_status == HighsModelStatus::kOptimal);
    REQUIRE(highs0.getInfo().simplex_iteration_count <= 0);
  }
}

HighsStatus zeroCostColSing() {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(1);
  lp.a_matrix_.start_.push_back(2);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.value_.push_back(0.5);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.value_.push_back(0.5);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(0.1);
  lp.row_upper_.push_back(0.9);

  lp.col_cost_.push_back(0);
  lp.col_cost_.push_back(1);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

// handled by doubleton equality
HighsStatus colSingDoubletonEquality() {
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 2;

  lp.a_matrix_.format_ = MatrixFormat::kColwise;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(2);
  lp.a_matrix_.start_.push_back(3);
  lp.a_matrix_.start_.push_back(4);
  lp.a_matrix_.start_.push_back(5);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(1);

  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(1);
  lp.row_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);
  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(1);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

HighsStatus colSingDoubletonInequality() {
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 2;

  lp.a_matrix_.format_ = MatrixFormat::kColwise;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(2);
  lp.a_matrix_.start_.push_back(3);
  lp.a_matrix_.start_.push_back(4);
  lp.a_matrix_.start_.push_back(5);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(1);

  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);
  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(1);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

// handled by doubleton equality
HighsStatus twoColSingDoubletonEquality() {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(1);
  lp.a_matrix_.start_.push_back(2);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(0);

  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(1);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

// handled by special case.
HighsStatus twoColSingDoubletonInequality() {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(1);
  lp.a_matrix_.start_.push_back(2);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(0);

  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  highs.run();
  status = highs.run();
  return status;
}

// No commas in test case name.
TEST_CASE("zero-cost [presolve-col-sing]") {
  std::cout << "Presolve 1." << std::endl;
  HighsStatus status = zeroCostColSing();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("col-sing-doubleton-eq [presolve-col-sing]") {
  std::cout << "Presolve 2." << std::endl;
  HighsStatus status = colSingDoubletonEquality();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("col-sing-doubleton-ineq [presolve-col-sing]") {
  std::cout << "Presolve 3." << std::endl;
  HighsStatus status = colSingDoubletonInequality();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("two-col-sing-doubleton-eq [presolve-col-sing]") {
  std::cout << "Presolve 4." << std::endl;
  HighsStatus status = twoColSingDoubletonEquality();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("two-col-sing-doubleton-ineq [presolve-col-sing]") {
  std::cout << "Presolve 5." << std::endl;
  HighsStatus status = twoColSingDoubletonInequality();
  std::string str = highsStatusToString(status);
  REQUIRE(str == "OK");
}

// test case failing
HighsStatus issue425() {
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(3);
  lp.a_matrix_.start_.push_back(5);
  lp.a_matrix_.start_.push_back(6);
  lp.a_matrix_.start_.push_back(7);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.index_.push_back(2);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.index_.push_back(3);
  lp.a_matrix_.value_.push_back(1);

  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.value_.push_back(2);
  lp.a_matrix_.index_.push_back(3);
  lp.a_matrix_.value_.push_back(1);

  lp.a_matrix_.index_.push_back(3);
  lp.a_matrix_.value_.push_back(1);

  lp.a_matrix_.index_.push_back(3);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, kHighsInf);

  std::vector<double> b{1, 2, 2, 4};
  lp.row_lower_ = b;
  lp.row_upper_ = b;

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

TEST_CASE("presolve-issue-425") {
  std::cout << std::endl;
  std::cout << "Presolve issue 425." << std::endl;
  HighsStatus status = issue425();
  REQUIRE(status == HighsStatus::kOk);
}
