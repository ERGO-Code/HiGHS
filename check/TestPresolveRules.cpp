#include <numeric>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "presolve/HPresolve.h"
#include "presolve/HighsPostsolveStack.h"

const bool dev_run = false;

void solveAndCheck(const std::string& message, const HighsLp& lp, Highs& h,
                   const std::string& solver, bool use_presolve,
                   const HighsInt require_presolved_model_num_col = -1,
                   const HighsInt require_presolved_model_num_row = -1,
                   const HighsInt require_presolved_model_num_nz = -1);

void presolveOffOn(const std::string& message, const HighsLp& lp, Highs& h,
                   const std::vector<std::string>& solvers,
                   const HighsInt require_presolved_model_num_col = -1,
                   const HighsInt require_presolved_model_num_row = -1,
                   const HighsInt require_presolved_model_num_nz = -1);

TEST_CASE("test-col-stuffing", "[highs_test_presolve_rules]") {
  HighsLp lp;

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("presolve_rule_test", kPresolveRuleColStuffing);
  REQUIRE(h.setOptionValue("presolve_rule_logging", true) == HighsStatus::kOk);
  // Initial sweep doesn't yield reductions, but switch it off for clarity
  REQUIRE(h.setOptionValue("presolve_rule_off",
                           1 << kPresolveRuleInitialSweep) == HighsStatus::kOk);
  const bool lp0 = true;
  const bool lp1 = true;
  const bool lp1a = true;
  const bool lp1b = true;

  if (lp0) {
    lp.num_col_ = 3;
    lp.num_row_ = 1;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_ = {1.8, 0.9, 1};
    lp.col_lower_.assign(lp.num_col_, 0);
    lp.col_upper_.assign(lp.num_col_, 1);
    lp.row_lower_ = {-kHighsInf};
    lp.row_upper_ = {4};
    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp.a_matrix_.start_ = {0, lp.num_col_};
    lp.a_matrix_.index_.resize(lp.num_col_);
    std::iota(lp.a_matrix_.index_.begin(), lp.a_matrix_.index_.end(), 0);
    lp.a_matrix_.value_ = {3, 2, 2};

    for (int k = 0; k < 2; k++) {
      if (dev_run) printf("\n3-variable knapsack: %s\n", k == 0 ? "LP" : "IP");
      REQUIRE(h.passModel(lp) == HighsStatus::kOk);
      h.setOptionValue("presolve_rule_test", kPresolveRuleColStuffing);
      h.run();
      if (dev_run) h.writeSolution("", 1);
      lp.integrality_.assign(lp.num_col_, HighsVarType::kInteger);
    }
    lp.clear();
  }

  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, 1);
  lp.row_lower_ = {2.0};
  lp.row_upper_ = {kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, lp.num_col_};
  lp.a_matrix_.index_.resize(lp.num_col_);
  std::iota(lp.a_matrix_.index_.begin(), lp.a_matrix_.index_.end(), 0);
  const std::vector<std::string> solvers = {kSimplexString, kIpmString,
                                            kHiPdlpString};
  if (lp1) {
    lp.col_cost_.assign(lp.num_col_, 1);
    lp.a_matrix_.value_.assign(lp.num_col_, 1);
    presolveOffOn("Capturing neos-787933 issue", lp, h, solvers);
  }
  if (lp1a) {
    lp.col_cost_ = {2, 1};
    lp.a_matrix_.value_.assign(lp.num_col_, 1);
    presolveOffOn("Variant A neos-787933 issue", lp, h, solvers);
  }
  if (lp1b) {
    lp.col_cost_ = {-2, -1};
    lp.a_matrix_.value_.assign(lp.num_col_, 1);
    presolveOffOn("Variant B neos-787933 issue", lp, h, solvers);
  }
  lp.clear();

  h.resetGlobalScheduler(true);
}

/*
TEST_CASE("test-weakly-dominated-col-upper", "[highs_test_presolve_rules]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.setOptionValue("presolve_rule_logging", true) == HighsStatus::kOk);
  // LP is
  //
  // min -y, subject to x+y <= 0, x >= 0; 0 <= x <= 1, y free
  //
  // Optimal solution is x = 1; y = -1, with x nonbasic with dual -1, and
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_lower_ = {-kHighsInf, -kHighsInf};
  lp.col_upper_ = {1,  kHighsInf};
  lp.row_lower_ = {-kHighsInf,         1};
  lp.row_upper_ = {         0, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 3};
  lp.a_matrix_.index_ = {0, 1, 0};
  lp.a_matrix_.value_ = {1, 1, 1};

  bool maximize_first = true;
  std::string sense_string = "";
  std::string test_string = "";

  for (HighsInt k = 0; k < 2; k++) {
    // Passes are minimize c^Tx and maximize -c^Tx according to
    // maximize_first
    if (maximize_first) {
      lp.sense_ = ObjSense::kMaximize;
      sense_string = "maximize";
      lp.col_cost_ = {0, 1};
    } else {
      lp.sense_ = ObjSense::kMinimize;
      sense_string = "minimize";
      lp.col_cost_ = {0, -1};
    }
    //  REQUIRE(h.setOptionValue("presolve_rule_test", 0) == HighsStatus::kOk);
    //  test_string = "vanilla-presolve-" + sense_string;
    //  presolveOffOn(test_string, lp, h);

    REQUIRE(h.setOptionValue("presolve_rule_test",
kPresolveRuleWeaklyDominatedColUpper) == HighsStatus::kOk);

    // test_string = "initial-sweep+test-weakly-dominated-col-upper-" +
sense_string;
    // presolveOffOn(test_string, lp, h, 1, 1, 1);

    REQUIRE(h.setOptionValue("presolve_rule_off", 1 <<
kPresolveRuleInitialSweep) == HighsStatus::kOk);

    test_string = "test-weakly-dominated-col-upper-" + sense_string;
    presolveOffOn(test_string, lp, h, 1, 2, 1);

    maximize_first = !maximize_first;
  }
  h.resetGlobalScheduler(true);
}
*/

TEST_CASE("test-parallel-rows-cut-ordering", "[highs_test_presolve_rules]") {
  // Rows 0 and 1 are parallel (both [1, 1]). Row 0 is marked as a
  // cut. detectParallelRowsAndCols must remove the cut row (0) and
  // keep the non-cut row (1), not the other way around.
  // Row 2 involves only col 0, breaking column parallelism.
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1, 2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {10, 10};
  lp.row_lower_ = {-kHighsInf, -kHighsInf, -kHighsInf};
  lp.row_upper_ = {5, 5, 3};
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.start_ = {0, 3, 5};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 1};

  HighsOptions options;
  options.presolve_rule_test = kPresolveRuleParallelRowsAndCols;
  options.presolve_rule_off = 1 << kPresolveRuleInitialSweep;
  options.output_flag = dev_run;

  HighsTimer timer;
  timer.start();

  presolve::HighsPostsolveStack postsolve_stack;
  postsolve_stack.initializeIndexMaps(lp.num_row_, lp.num_col_);
  // Mark parallel row 0 as a cut
  postsolve_stack.setRowType(0,
                             presolve::HighsPostsolveStack::OrigRowType::kCut);

  presolve::HPresolve presolve;
  presolve.setInput(lp, options, -1, &timer);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  timer.stop();
  REQUIRE(status == HighsModelStatus::kNotset);
  // One row must have been removed
  REQUIRE(lp.num_row_ == 1);
  // The surviving row must be original row 1 (non-cut), not row 0 (cut)
  REQUIRE(postsolve_stack.getOrigRowIndex(0) == 1);
  REQUIRE(!postsolve_stack.isCutRow(0));
}

TEST_CASE("test-effective-costs", "[highs_test_presolve]") {
  // Debugging ZeroObjSingletonContinuousCol for germanrr highlighted
  // the deficiency in computing the active_cost_norm when the
  // objective is f = z, with z = c^Tx and z free. In
  // HighsSolution.cpp is the method getEffectiveCosts that
  // substitutes all free column singletons into the objective to get
  // the "effective costs".
  Highs h;
  //  h.setOptionValue("output_flag", dev_run);
  bool test_all = true;
  bool test_lp0 = test_all;
  bool test_lp1 = test_all;
  bool test_lp2 = test_all;

  if (test_lp0) {
    HighsLp lp;
  // First LP is
  //
  // min 4z
  //
  // -1 <=    x + y - 2z <= 1
  //
  // -1 <= 201x + y      <= 1
  //
  // 0 <= x <= 1, y, z free
  //
  // where the bounds on the two constraints and non-unit coefficients
  // of z in the objective and first contraint give code coverage
  //
  // Aiming to minimize 4z, and bound is given by 2z >= x + y - 1, so
  // substitute z = (x+y-1)/2 into the objective to give
  //
  // min 2x + 2y - 2
  //
  // y is then minimized with bound is given by y >= -201x - 1, so
  // substitute y = -201x - 1 into the objective to give
  //
  // min 2x +(-402x-2) - 2 = -400x - 4
  //
  // This function is minimized when x = 1 to give y = -202 and z =
  // -101 with objective -404
  //
  // The optimal dual values are -400 for x, -2 for row 0 and 2 for
  // row 1. However, although this example tests code coverage on
  // identifying free column singletons and a double free column
  // singleton identified in getEffectiveCosts, the dual of -400 for
  // the only nonbasic column means that there are no active costs, so
  // active_cost_norm is zero (hence absolute and relative dual
  // infeasibility measures are identical).
  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {0, 0, 4};
  lp.col_lower_ = {0, -kHighsInf, -kHighsInf};
  lp.col_upper_ = {1, kHighsInf, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 3, 5};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1};
  lp.a_matrix_.value_ = {1, 1, -2, 201, 1};
  lp.row_lower_ = {-1, -1};
  lp.row_upper_ = {1, 1};
  h.passModel(lp);
  h.setOptionValue("log_dev_level", 1);
  h.setOptionValue("presolve_rule_logging", kHighsOnString);
  h.run();
  REQUIRE(h.getInfo().active_cost_norm == 0);
  printf("LP0: active_cost_norm = %g\n", h.getInfo().active_cost_norm);
  }
  
  if (test_lp1) {
  HighsLp lp;
  // Here's a simpler example that reflects the behaviour observed
  // with germanrr, where the cost row of the matrix introduced many
  // large costs. Hence the presolved model had a large value for
  // active_cost_norm but, after postsolve, the model had
  // active_cost_norm = 1.

  double cost = 1e5;
  double eps = 1e-4;
  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {0, 0, 1};
  lp.col_lower_ = {0, 0, -kHighsInf};
  lp.col_upper_ = {1, 1, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 3, 5};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1};
  lp.a_matrix_.value_ = {cost, cost - eps, 1, 1, 1, 1};
  lp.row_lower_ = {0, 1};
  lp.row_upper_ = {0, 1};
  h.passModel(lp);

  h.run();
  REQUIRE(h.getInfo().active_cost_norm == cost);
  printf("LP1: active_cost_norm = %g\n", h.getInfo().active_cost_norm);
  }
  test_lp2 = true;
  if (test_lp2) {
  // Finally gas11 has 61 free column singletons
  h.setOptionValue("output_flag", true);
  const std::string model = "gas11";
   std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  h.readModel(model_file);
  REQUIRE(h.setOptionValue(kPresolveString, kHighsOffString) == HighsStatus::kOk);

  HighsStatus return_status = h.run();
  REQUIRE(return_status == HighsStatus::kOk);

  printf("LP2: active_cost_norm = %g\n", h.getInfo().active_cost_norm);
  
  }
  h.resetGlobalScheduler(true);
}

void solveAndCheck(const std::string& message, const HighsLp& lp, Highs& h,
                   const std::string& solver, bool use_presolve,
                   const HighsInt require_presolved_model_num_col,
                   const HighsInt require_presolved_model_num_row,
                   const HighsInt require_presolved_model_num_nz) {
  const HighsRunData& run_data = h.getRunData();
  std::string run_crossover = kHighsOnString;
  bool basis_postsolve = true;
  if (solver == kIpmString) {
    run_crossover = kHighsOffString;
    basis_postsolve = false;
  } else if (solver == kHiPdlpString) {
    basis_postsolve = false;
  }
  std::string presolve = use_presolve ? kHighsOnString : kHighsOffString;
  h.setOptionValue(kPresolveString, presolve);
  h.setOptionValue(kRunCrossoverString, run_crossover);
  h.setOptionValue(kSolverString, solver);
  if (dev_run)
    printf("\n============\n%s: presolve = %s; solver = %s%s\n============\n\n",
           message.c_str(), presolve.c_str(), solver.c_str(),
           solver == kIpmString ? ("; run_crossover = " + run_crossover).c_str()
                                : "");
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.run();
  if (dev_run) h.writeSolution("", 1);
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(h.getInfo().num_primal_infeasibilities == 0);
  REQUIRE(h.getInfo().num_dual_infeasibilities == 0);
  if (use_presolve) {
    // Ensure that the model is reduced as expected
    if (require_presolved_model_num_col >= 0)
      REQUIRE(run_data.presolved_model_num_col ==
              require_presolved_model_num_col);
    if (require_presolved_model_num_row >= 0)
      REQUIRE(run_data.presolved_model_num_row ==
              require_presolved_model_num_row);
    if (require_presolved_model_num_nz >= 0)
      REQUIRE(run_data.presolved_model_num_nz ==
              require_presolved_model_num_nz);
    if (require_presolved_model_num_col == 0 &&
        require_presolved_model_num_row == 0)
      REQUIRE(h.getInfo().simplex_iteration_count == 0);
    // Ensure that any basis postsolve is correct
    if (basis_postsolve)
      REQUIRE(run_data.num_simplex_iterations_after_postsolve == 0);
  }
}

void presolveOffOn(const std::string& message, const HighsLp& lp, Highs& h,
                   const std::vector<std::string>& solvers,
                   const HighsInt require_presolved_model_num_col,
                   const HighsInt require_presolved_model_num_row,
                   const HighsInt require_presolved_model_num_nz) {
  // Presolve off - to get the optimal solution to debug presolve
  solveAndCheck(message, lp, h, kSimplexString, false);
  // Presolve on with each solver
  for (const std::string& solver : solvers) {
    solveAndCheck(message, lp, h, solver, true, require_presolved_model_num_col,
                  require_presolved_model_num_row,
                  require_presolved_model_num_nz);
  }
}
