#include <numeric>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "mip/HighsCliqueTable.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "parallel/HighsParallel.h"
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

TEST_CASE("test-implied-bound-aggregation", "[highs_test_presolve_rules]") {
  // Three VLB rows on a continuous variable z with binary variables
  // x0, x1, x2 in a clique (from assignment equation). Probing should
  // discover the clique and aggregation should merge the three VLBs
  // into one stronger row.
  //
  // min z
  // x0 + x1 + x2 = 1
  // z >= 10*x0  (stored as -z + 10*x0 <= 0)
  // z >= 20*x1  (stored as -z + 20*x1 <= 0)
  // z >= 30*x2  (stored as -z + 30*x2 <= 0)
  // z >= 0, x0,x1,x2 binary
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {kHighsInf, 1, 1, 1};
  lp.integrality_ = {HighsVarType::kContinuous, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {1, -kHighsInf, -kHighsInf, -kHighsInf};
  lp.row_upper_ = {1, 0, 0, 0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 4;
  lp.a_matrix_.start_ = {0, 3, 5, 7, 9};
  lp.a_matrix_.index_ = {1, 2, 3, 0, 1, 0, 2, 0, 3};
  lp.a_matrix_.value_ = {-1, -1, -1, 1, 10, 1, 20, 1, 30};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test", kPresolveRuleProbing);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  // Populate VLBs: z >= 10*x0 (row 1), z >= 20*x1 (row 2), z >= 30*x2 (row 3)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addVLB(0, 1, 10.0, 0.0, 1);
  implications.addVLB(0, 2, 20.0, 0.0, 2);
  implications.addVLB(0, 3, 30.0, 0.0, 3);

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // The three VLBs form a clique, so aggregation merges them into one
  // row: z - 10*x0 - 20*x1 - 30*x2 >= 0
  const HighsLp& presolved = *mipsolver.model_;
  HighsInt agg_row = -1;
  for (HighsInt i = 0; i < mipsolver.numRow(); i++) {
    if (postsolve_stack.getOrigRowIndex(i) == 0) continue;
    REQUIRE(agg_row == -1);
    agg_row = i;
  }
  REQUIRE(agg_row >= 0);
  REQUIRE(presolved.row_lower_[agg_row] == 0.0);
  REQUIRE(presolved.row_upper_[agg_row] == kHighsInf);

  std::vector<double> coeffs(4, 0.0);
  for (HighsInt j = 0; j < presolved.num_col_; j++) {
    for (HighsInt p = presolved.a_matrix_.start_[j];
         p < presolved.a_matrix_.start_[j + 1]; p++) {
      if (presolved.a_matrix_.index_[p] != agg_row) continue;
      coeffs[postsolve_stack.getOrigColIndex(j)] =
          presolved.a_matrix_.value_[p];
    }
  }
  REQUIRE(coeffs[0] == 1.0);
  REQUIRE(coeffs[1] == -10.0);
  REQUIRE(coeffs[2] == -20.0);
  REQUIRE(coeffs[3] == -30.0);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-implied-bound-aggregation-vub", "[highs_test_presolve_rules]") {
  // Three VUB rows on a continuous variable z with binary variables
  // x0, x1, x2 where x0+x1+x2 = 2 (so at most one xi can be 0).
  // The complement literals (1-x0), (1-x1), (1-x2) form a clique.
  // VUBs have positive coefficients so the activating literal is x=0
  // (val=0), exercising the complement offset path in mergeCliques.
  //
  // max z
  // x0 + x1 + x2 = 2
  // z <= 10*x0 + 20  (stored as z - 10*x0 <= 20)
  // z <= 20*x1 + 10  (stored as z - 20*x1 <= 10)
  // z <= 30*x2        (stored as z - 30*x2 <= 0)
  // 0 <= z <= 30, x0,x1,x2 binary
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;
  lp.sense_ = ObjSense::kMaximize;
  lp.col_cost_ = {1, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {30, 1, 1, 1};
  lp.integrality_ = {HighsVarType::kContinuous, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {2, -kHighsInf, -kHighsInf, -kHighsInf};
  lp.row_upper_ = {2, 20, 10, 0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 4;
  lp.a_matrix_.start_ = {0, 3, 5, 7, 9};
  lp.a_matrix_.index_ = {1, 2, 3, 0, 1, 0, 2, 0, 3};
  lp.a_matrix_.value_ = {1, 1, 1, 1, -10, 1, -20, 1, -30};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test", kPresolveRuleProbing);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  // Populate VUBs with positive coef (val=0 complement path):
  // z <= 10*x0 + 20 (row 1), z <= 20*x1 + 10 (row 2), z <= 30*x2 (row 3)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addVUB(0, 1, 10.0, 20.0, 30.0, false, 1);
  implications.addVUB(0, 2, 20.0, 10.0, 30.0, false, 2);
  implications.addVUB(0, 3, 30.0, 0.0, 30.0, false, 3);

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // The three VUBs form a clique, so aggregation merges them into one
  // row: z - 10*x0 - 20*x1 - 30*x2 <= -30
  const HighsLp& presolved = *mipsolver.model_;
  HighsInt agg_row = -1;
  for (HighsInt i = 0; i < mipsolver.numRow(); i++) {
    if (postsolve_stack.getOrigRowIndex(i) == 0) continue;
    REQUIRE(agg_row == -1);
    agg_row = i;
  }
  REQUIRE(agg_row >= 0);
  REQUIRE(presolved.row_lower_[agg_row] == -kHighsInf);
  REQUIRE(presolved.row_upper_[agg_row] == -30.0);

  std::vector<double> coeffs(4, 0.0);
  for (HighsInt j = 0; j < presolved.num_col_; j++) {
    for (HighsInt p = presolved.a_matrix_.start_[j];
         p < presolved.a_matrix_.start_[j + 1]; p++) {
      if (presolved.a_matrix_.index_[p] != agg_row) continue;
      coeffs[postsolve_stack.getOrigColIndex(j)] =
          presolved.a_matrix_.value_[p];
    }
  }
  REQUIRE(coeffs[0] == 1.0);
  REQUIRE(coeffs[1] == -10.0);
  REQUIRE(coeffs[2] == -20.0);
  REQUIRE(coeffs[3] == -30.0);

  HighsTaskExecutor::shutdown(true);
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
  h.setOptionValue("output_flag", dev_run);
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
    // where the bounds on the two constraints and non-unit
    // coefficients of z in the objective and first contraint give
    // code coverage
    //
    // Aiming to minimize 4z, and bound is given by 2z >= x + y - 1,
    // so substitute z = (x+y-1)/2 into the objective to give
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
    // the only nonbasic column means that there are no active costs,
    // so active_cost_norm is zero (hence absolute and relative dual
    // infeasibility measures are identical).
    lp.model_name_ = "LP0";
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
    lp.model_name_ = "LP1";
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
  }
  if (test_lp2) {
    // Finally gas11 has 61 free column singletons: 55 in the first
    // pass, and 6 in the second.
    const std::string model = "gas11";
    std::string model_file =
        std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    REQUIRE(h.readModel(model_file) == HighsStatus::kWarning);
    REQUIRE(h.setOptionValue(kPresolveString, kHighsOffString) ==
            HighsStatus::kOk);
    HighsStatus return_status = h.run();
    REQUIRE(return_status == HighsStatus::kOk);
    double active_cost_norm = 2.000000001e+7;
    REQUIRE(std::fabs(h.getInfo().active_cost_norm - active_cost_norm) <= 1e-8);
  }
}

TEST_CASE("test-fourier-motzkin", "[highs_test_presolve_rules]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("presolve_rule_test", kPresolveRuleFourierMotzkin);
  h.setOptionValue("presolve_rule_logging", true);
  h.setOptionValue("log_dev_level", 1);

  const bool lp0 = true;
  const bool lp1 = true;  // Makes eliminations marginal, and leaves x2=0
  const bool lp2 = true;

  // No PDLP due to numerical issues with FM postsolve
  const std::vector<std::string> solvers = {kSimplexString, kIpmString};

  // From "A novel linear optimization presolve technique based on
  // Fourier-Motzkin elimination", Zhang, Ploskas and Sahinidis,
  // Mathematical Programming Computation (2026) 18:345-378
  HighsLp lp;

  lp.num_col_ = 4;
  lp.num_row_ = 3;

  lp.col_cost_.assign(lp.num_col_, 0);
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, kHighsInf);
  lp.col_upper_[0] = 40.0;

  lp.row_lower_.assign(lp.num_row_, -kHighsInf);
  lp.row_upper_ = {-30, 50, 40};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 3, 6, 9};
  lp.a_matrix_.index_ = {0, 1, 3, 1, 2, 3, 1, 2, 3};
  lp.a_matrix_.value_ = {-1, 1, -1, 2, 1, 2, 3, -1, 3};

  if (lp0) {
    REQUIRE(h.passModel(lp) == HighsStatus::kOk);
    presolveOffOn("FM example from paper", lp, h, solvers);
  }

  lp.col_upper_[0] = 5.0;
  lp.row_upper_ = {-30, 75, 50};

  if (lp1) {
    REQUIRE(h.passModel(lp) == HighsStatus::kOk);
    presolveOffOn("FM example from paper - tightened", lp, h, solvers);
  }

  lp.col_cost_ = {1, 2, 3, 4};

  REQUIRE(h.passModel(lp) == HighsStatus::kOk);

  if (lp2) {
    // Objective reformulation is needed since all costs are nonzero
    h.setOptionValue("presolve_fm_obj_reformulation", true);
    presolveOffOn("FM example from paper - tightened and with costs", lp, h,
                  solvers, 1, 6, 6);
  }

  h.resetGlobalScheduler(true);
}

TEST_CASE("test-parallel-cols-merge-lp", "[highs_test_presolve_rules]") {
  // Example 8 (LP) from Gamrath et al. 2015: parallel column merge.
  //
  //   min  2x1 + 4x2 + x3
  //   s.t. -x1 - 2x2 - x3 <= -10
  //        0 <= x1 <= 3, 0 <= x2 <= 4, 0 <= x3 <= 5
  //
  // Columns 1 and 2 are parallel with lambda = 2, c2 = lambda*c1.
  // Merge y := x1 + 2x2 in [0, 11], cost 2y.
  // Presolved: min 2y + x3, -y - x3 <= -10, y in [0,11], x3 in [0,5].
  // Optimal x* = (0, 2.5, 5), obj = 15.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {2, 4, 1};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {3, 4, 5};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {-10};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 3};
  lp.a_matrix_.index_ = {0, 1, 2};
  lp.a_matrix_.value_ = {-1, -2, -1};

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("presolve_rule_test", kPresolveRuleParallelRowsAndCols);
  h.presolve();
  REQUIRE(h.getPresolvedLp().num_col_ == 2);
  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(h.getInfo().num_primal_infeasibilities == 0);
  REQUIRE(std::abs(h.getObjectiveValue() - 15) < 1e-8);

  h.resetGlobalScheduler(true);
}

TEST_CASE("test-parallel-cols-merge-ip", "[highs_test_presolve_rules]") {
  // Example 8 (IP) from Gamrath et al. 2015: parallel column merge.
  //
  //   min  2x1 + 4x2 + x3
  //   s.t. -x1 - 2x2 - x3 <= -10
  //        0 <= x1 <= 3, 0 <= x2 <= 4, 0 <= x3 <= 5
  //        x1, x2, x3 integer
  //
  // Columns 1 and 2 are parallel with lambda = 2, c2 = lambda*c1.
  // Merge y := x1 + 2x2 in {0, ..., 11}, cost 2y.
  // Presolved: min 2y + x3, -y - x3 <= -10, y in [0,11], x3 in [0,5].
  // Optimal x* = (1, 2, 5), obj = 15.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {2, 4, 1};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {3, 4, 5};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {-10};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 3};
  lp.a_matrix_.index_ = {0, 1, 2};
  lp.a_matrix_.value_ = {-1, -2, -1};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("presolve_rule_test", kPresolveRuleParallelRowsAndCols);
  h.presolve();
  REQUIRE(h.getPresolvedLp().num_col_ == 2);
  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(h.getInfo().num_primal_infeasibilities == 0);
  REQUIRE(std::abs(h.getObjectiveValue() - 15) < 1e-8);

  h.resetGlobalScheduler(true);
}

TEST_CASE("test-parallel-cols-merge-floor-rounding",
          "[highs_test_presolve_rules]") {
  // Exercises the floor branch in DuplicateColumn::undo postsolve rounding.
  // Mixed integer/continuous parallel columns with colLower = -inf so that
  // the initial postsolve decomposition pushes duplicateCol below its lower
  // bound. After clipping duplicateCol to its lower bound, col is recomputed
  // and floor-rounded.
  //
  //   min  x1 + 2x2 +  x3 + x4
  //   s.t. x1 + 2x2 + 3x3 + x4 >= 10
  //                    x3 + x4 <= 6
  //        x1 integer in (-inf, 5], x2 continuous in [3, 4],
  //        x3 integer in [0   , 5], x4 continuous in [0, 5]
  //
  // col = x1 (integer), duplicateCol = x2 (continuous), colScale = 2.
  // Merge y := x1 + 2x2 in (-inf, 13], cost y.
  // After merge: min y + x3 + x4, y + 3x3 + x4 >= 10, x3 + x4 <= 6.
  // Optimal: x3 = 5, x4 = 0, y = -5, obj = 0.
  // Postsolve: col = min(0, 5) = 0, duplicateCol = (-5-0)/2 = -2.5 < 3.
  // Clip duplicateCol to 3, recompute col = -5 - 2*3 = -11, floor(-11) = -11.
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 2;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1, 2, 1, 1};
  lp.col_lower_ = {-kHighsInf, 3, 0, 0};
  lp.col_upper_ = {5, 4, 5, 5};
  lp.row_lower_ = {10, -kHighsInf};
  lp.row_upper_ = {kHighsInf, 6};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 4, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 2, 3};
  lp.a_matrix_.value_ = {1, 2, 3, 1, 1, 1};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kContinuous,
                     HighsVarType::kInteger, HighsVarType::kContinuous};

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("presolve_rule_test", kPresolveRuleParallelRowsAndCols);
  h.presolve();
  REQUIRE(h.getPresolvedLp().num_col_ == 3);
  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(h.getInfo().num_primal_infeasibilities == 0);
  REQUIRE(std::abs(h.getObjectiveValue()) < 1e-8);

  h.resetGlobalScheduler(true);
}

TEST_CASE("test-parallel-cols-merge-ceil-rounding",
          "[highs_test_presolve_rules]") {
  // Exercises the ceil branch in DuplicateColumn::undo postsolve rounding.
  // Mixed integer/continuous parallel columns with finite colLower so that
  // the initial postsolve decomposition pushes duplicateCol above its upper
  // bound. After clipping duplicateCol to its upper bound, col is recomputed
  // and ceil-rounded.
  //
  //   min  -x1 - 2x2 +  x3 + x4
  //   s.t.  x1 + 2x2 + 3x3 + x4 >= 10
  //        2x1 + 4x2 +  x3 + x4 <= 21
  //                     x3 + x4 <= 6
  //        x1 integer in [0, 10], x2 continuous in [0, 2],
  //        x3 integer in [0,  5], x4 continuous in [0, 5]
  //
  // col = x1 (integer), duplicateCol = x2 (continuous), colScale = 2.
  // Merge y := x1 + 2x2 in [0, 14], cost -y.
  // After merge: min -y + x3 + x4, y + 3x3 + x4 >= 10, 2y + x3 + x4 <= 21.
  // Optimal: x3 = 0, x4 = 0, y = 10.5, obj = -10.5.
  // Postsolve: col = colLower = 0, duplicateCol = (10.5-0)/2 = 5.25 > 2.
  // Clip duplicateCol to 2, recompute col = 10.5 - 2*2 = 6.5, ceil(6.5) = 7.
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 3;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {-1, -2, 1, 1};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {10, 2, 5, 5};
  lp.row_lower_ = {10, -kHighsInf, -kHighsInf};
  lp.row_upper_ = {kHighsInf, 21, 6};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 4, 8, 10};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 0, 1, 2, 3, 2, 3};
  lp.a_matrix_.value_ = {1, 2, 3, 1, 2, 4, 1, 1, 1, 1};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kContinuous,
                     HighsVarType::kInteger, HighsVarType::kContinuous};

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("presolve_rule_test", kPresolveRuleParallelRowsAndCols);
  h.presolve();
  REQUIRE(h.getPresolvedLp().num_col_ == 3);
  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  REQUIRE(h.getInfo().num_primal_infeasibilities == 0);
  REQUIRE(std::abs(h.getObjectiveValue() + 10.5) < 1e-8);

  h.resetGlobalScheduler(true);
}

TEST_CASE("test-clique-extract-origin", "[highs_test_presolve_rules]") {
  // normaliseCliqueRows flips >= rows to <= form so that extractCliques
  // recognises them as set packing constraints with a row origin.
  // Clique merging then extends the 3-clique from row 0 with (x3,0)
  // and subsumes the size-2 cliques, deleting their origin rows.
  //   row 0: x0 + x1 + x2 <= 1       (set packing, 3-clique)
  //   row 1: -x0 + x3 >= 0           (x3 >= x0, implication)
  //   row 2: -x1 + x3 >= 0           (x3 >= x1, implication)
  //   row 3: -x2 + x3 >= 0           (x3 >= x2, implication)
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0, 1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0, 0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0, 1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf, 0.0, 0.0, 0.0};
  lp.row_upper_ = {1.0, kHighsInf, kHighsInf, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 4;
  lp.a_matrix_.start_ = {0, 2, 4, 6, 9};
  lp.a_matrix_.index_ = {0, 1, 0, 2, 0, 3, 1, 2, 3};
  lp.a_matrix_.value_ = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test", kPresolveRuleProbing);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  REQUIRE(mipsolver.numRow() <= 1);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-normalise-unequal-coeff", "[highs_test_presolve_rules]") {
  // Same as test-clique-extract-origin but row 0 has unequal coefficients,
  // exercising the normalisation path in normaliseCliqueRows.
  //   row 0: x0 + 2*x1 + 2*x2 <= 2  (unequal coeffs, clique covers all 3)
  //   row 1: -x0 + x3 >= 0
  //   row 2: -x1 + x3 >= 0
  //   row 3: -x2 + x3 >= 0
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0, 1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0, 0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0, 1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf, 0.0, 0.0, 0.0};
  lp.row_upper_ = {2.0, kHighsInf, kHighsInf, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 4;
  lp.a_matrix_.start_ = {0, 2, 4, 6, 9};
  lp.a_matrix_.index_ = {0, 1, 0, 2, 0, 3, 1, 2, 3};
  lp.a_matrix_.value_ = {1.0, -1.0, 2.0, -1.0, 2.0, -1.0, 1.0, 1.0, 1.0};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test", kPresolveRuleProbing);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  REQUIRE(mipsolver.numRow() <= 1);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-normalise-equation", "[highs_test_presolve_rules]") {
  // 2x0 + 2x1 + 2x2 = 2 is a valid set partitioning constraint
  // (coefficients == rhs). normaliseCliqueRows should normalise it
  // and extractCliques should find a 3-clique with all pairs.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {2.0};
  lp.row_upper_ = {2.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {2.0, 2.0, 2.0};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test", kPresolveRuleProbing);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  REQUIRE(cliquetable.haveCommonClique({0, 1}, {1, 1}));
  REQUIRE(cliquetable.haveCommonClique({0, 1}, {2, 1}));
  REQUIRE(cliquetable.haveCommonClique({1, 1}, {2, 1}));

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-normalise-trivial-fixing", "[highs_test_presolve_rules]") {
  // 3*x0 + x1 + x2 <= 2: coefficient of x0 exceeds rhs,
  // so normaliseCliqueRows fixes x0 to lower bound.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {2.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {3.0, 1.0, 1.0};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  auto result = presolve.normaliseCliqueRows(postsolve_stack);
  mipsolver.timer_.stop();
  REQUIRE(static_cast<int>(result) == 0);
  // x0 must have been fixed (recorded on the postsolve stack)
  REQUIRE(postsolve_stack.numReductions() == 1);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-normalise-complemented", "[highs_test_presolve_rules]") {
  // 3*x0 - x1 - x2 <= -1: after complementing x1, x2 the transformed
  // row is 3*x0 + (1-x1) + (1-x2) <= rhs=1. x0 has coefficient 3 > 1,
  // so it is trivially fixed to 0. remaining complemented variables are
  // normalised to -1 coefficients with row_upper = 1 - 2 = -1.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {-1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {3.0, -1.0, -1.0};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  auto result = presolve.normaliseCliqueRows(postsolve_stack);
  mipsolver.timer_.stop();
  REQUIRE(static_cast<int>(result) == 0);
  // x0 must have been fixed (recorded on the postsolve stack)
  REQUIRE(postsolve_stack.numReductions() == 1);
  // row should be normalised to -x1 - x2 <= -1
  REQUIRE(mipsolver.model_->row_upper_[0] == -1.0);
  REQUIRE(mipsolver.model_->row_lower_[0] == -kHighsInf);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-normalise-fix-to-upper", "[highs_test_presolve_rules]") {
  // -3*x0 - x1 - x2 <= -4: after complementing all variables the transformed
  // row is 3*(1-x0) + (1-x1) + (1-x2) <= rhs=1. x0 has complemented
  // coefficient 3 > 1 and complementation -1, so it is fixed to upper bound.
  // Remaining variables are normalised to -1 coefficients with
  // row_upper = 1 - 2 = -1.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {-4.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {-3.0, -1.0, -1.0};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  auto result = presolve.normaliseCliqueRows(postsolve_stack);
  mipsolver.timer_.stop();
  REQUIRE(static_cast<int>(result) == 0);
  // x0 must have been fixed to upper bound
  REQUIRE(postsolve_stack.numReductions() == 1);
  // row should be normalised to -x1 - x2 <= -1
  REQUIRE(mipsolver.model_->row_upper_[0] == -1.0);
  REQUIRE(mipsolver.model_->row_lower_[0] == -kHighsInf);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-normalise-non-integral-rhs", "[highs_test_presolve_rules]") {
  // 2*x0 + 2*x1 <= 3.7: integralScale = 0.5 gives x0 + x1 <= 1.85.
  // floor(1.85) = 1, so the normalised row is x0 + x1 <= 1 (a clique).
  // With std::round the rhs would become 2, losing the clique.
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {3.7};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 2;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {2.0, 2.0};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsProfiling profiling;

  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.timer_.start();
  profiling.initialize(mipsolver.timer_, true, true);
  mipsolver.setProfiling(&profiling);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->init();
  mipsolver.mipdata_->setupDomainPropagation();

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  auto result = presolve.normaliseCliqueRows(postsolve_stack);
  mipsolver.timer_.stop();
  REQUIRE(static_cast<int>(result) == 0);
  // row should be normalised to x0 + x1 <= 1
  REQUIRE(mipsolver.model_->row_upper_[0] == 1.0);
  REQUIRE(mipsolver.model_->row_lower_[0] == -kHighsInf);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-clique-no-delete-ranged-row", "[highs_test_presolve_rules]") {
  // Negative test: ranged rows must not be deleted by clique merging because
  // the extracted clique is a relaxation (only captures one side).
  //   row 0: x0 + x1 + x2 <= 1       (set packing, 3-clique)
  //   row 1: 0 <= -x0 + x3 <= 1      (ranged)
  //   row 2: 0 <= -x1 + x3 <= 1      (ranged)
  //   row 3: 0 <= -x2 + x3 <= 1      (ranged)
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0, 1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0, 0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0, 1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf, 0.0, 0.0, 0.0};
  lp.row_upper_ = {1.0, 1.0, 1.0, 1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 4;
  lp.a_matrix_.start_ = {0, 2, 4, 6, 9};
  lp.a_matrix_.index_ = {0, 1, 0, 2, 0, 3, 1, 2, 3};
  lp.a_matrix_.value_ = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0};

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->feastol = 1e-6;
  mipsolver.mipdata_->postSolveStack.initializeIndexMaps(4, 4);
  mipsolver.mipdata_->setupDomainPropagation();

  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  HighsDomain& domain = mipsolver.mipdata_->getDomain();

  cliquetable.extractCliques(mipsolver);
  cliquetable.runCliqueMerging(domain);

  const std::vector<HighsInt>& deleted = cliquetable.getDeletedRows();
  REQUIRE(deleted.empty());

  highs.resetGlobalScheduler(true);
}

TEST_CASE("test-clique-implied-equality", "[highs_test_presolve_rules]") {
  // 0.5 <= x0 + x1 <= 1: after rounding, lhs = ceil(0.5) = 1 and
  // rhs = floor(1) = 1, so this is an equality clique (set partitioning).
  // Without the fix, lhs was not computed and the clique was stored as
  // non-equality (set packing), losing the lower bound if the row is deleted.
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {1.0, 1.0};
  lp.col_lower_ = {0.0, 0.0};
  lp.col_upper_ = {1.0, 1.0};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {0.5};
  lp.row_upper_ = {1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 2;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1.0, 1.0};

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);

  HighsCallback callback(&highs);
  const HighsOptions& options = highs.getOptions();
  HighsSolution solution;
  HighsMipSolver mipsolver(callback, options, lp, solution);
  mipsolver.mipdata_ =
      std::unique_ptr<HighsMipSolverData>(new HighsMipSolverData(mipsolver));
  mipsolver.mipdata_->feastol = 1e-6;
  mipsolver.mipdata_->postSolveStack.initializeIndexMaps(1, 2);
  mipsolver.mipdata_->setupDomainPropagation();

  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  HighsDomain& domain = mipsolver.mipdata_->getDomain();

  cliquetable.extractCliques(mipsolver);

  // fix x0 = 0 and propagate; for an equality clique x0 + x1 = 1,
  // this forces x1 = 1. For a non-equality clique x0 + x1 <= 1,
  // fixing x0 = 0 does not constrain x1.
  domain.fixCol(0, 0.0);
  domain.propagate();
  REQUIRE(!domain.infeasible());
  REQUIRE(domain.isFixed(1));
  REQUIRE(domain.col_lower_[1] == 1.0);

  highs.resetGlobalScheduler(true);
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
