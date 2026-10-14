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
  bool test_lp1 = test_all;
  bool test_lp2 = test_all;

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

TEST_CASE("test-impl-aware-clique-extraction", "[highs_test_presolve_rules]") {
  // two-column clique discovered via implied bounds on non-binaries.
  //
  // row 0: 3*x0 + 3*x1 + y1 + y2 <= 12
  // row 1: y1 >= 5*x0    (stored as -5*x0 + y1 >= 0)
  // row 2: y2 >= 5*x1    (stored as -5*x1 + y2 >= 0)
  // x0, x1 binary, y1, y2 continuous in [0, 10]
  //
  // standard clique extraction on row 0 fails: binary coefficients 3+3=6 < 12.
  // on row 0 alone, x0=1 gives y1 >= 5, but 3 + 3*x1 + 5 + y2 <= 12 does not
  // fix x1 (probing would, via y2 <= 4 and row 2). implAwareConstrPropagation
  // uses both implications on row 0 simultaneously:
  // combined weight = 3 + 3 + max(5,0) + max(0,5) = 16 > 12.
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 3;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 10, 10};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kContinuous, HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf, 0, 0};
  lp.row_upper_ = {12, kHighsInf, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 3;
  lp.a_matrix_.start_ = {0, 2, 4, 6, 8};
  lp.a_matrix_.index_ = {0, 1, 0, 2, 0, 1, 0, 2};
  lp.a_matrix_.value_ = {3, -5, 3, -5, 1, 1, 1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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
  // cliques are only added once cliques have been extracted (by probing)
  mipsolver.mipdata_->cliquesExtracted = true;

  // manually populate VLBs and implications (normally done by probing)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  // x0=1 implies y1 >= 5 (VLB: y1 >= 5*x0)
  implications.addVLB(2, 0, 5.0, 0.0, 1);
  implications.addImplication(0, 1, 2,
                              HighsImplications::Implication{5.0, kHighsInf});
  // x1=1 implies y2 >= 5 (VLB: y2 >= 5*x1)
  implications.addVLB(3, 1, 5.0, 0.0, 2);
  implications.addImplication(1, 1, 3,
                              HighsImplications::Implication{5.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // implAwareConstrPropagation should have found that x0=1 and x1=1
  // can't coexist (combined implied activity exceeds row 0's upper bound)
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  REQUIRE(cliquetable.haveCommonClique({0, 1}, {1, 1}));

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-conflict-clique", "[highs_test_presolve_rules]") {
  // two-column clique from contradicting implied bounds on a non-binary.
  //
  // row 0: 4*x0 + 4*x1 + y <= 12
  // x0, x1 binary; y, z continuous in [0, 10]; z is not in any row
  // VLBs: y >= 3*x0, y >= 3*x1
  // implications: x0=1 => z >= 6, x1=1 => z <= 4
  //
  // threshold = 12 and the implied activity of x0=1 and x1=1 is only
  // 4 + 4 + max(3, 3) = 11 <= 12, but they imply contradicting bounds on z
  // (z >= 6 and z <= 4), so they form a clique independently of the row.
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 10, 10};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kContinuous, HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {12};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {4, 4, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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
  // cliques are only added once cliques have been extracted (by probing)
  mipsolver.mipdata_->cliquesExtracted = true;

  // VLBs: y >= 3*x0, y >= 3*x1
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addVLB(2, 0, 3.0, 0.0);
  implications.addVLB(2, 1, 3.0, 0.0);
  // x0=1 => z >= 6, x1=1 => z <= 4
  implications.addImplication(0, 1, 3,
                              HighsImplications::Implication{6.0, kHighsInf});
  implications.addImplication(1, 1, 3,
                              HighsImplications::Implication{-kHighsInf, 4.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0=1 and x1=1 can't coexist
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  REQUIRE(cliquetable.haveCommonClique({0, 1}, {1, 1}));

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-nonbinary-tightening",
          "[highs_test_presolve_rules]") {
  // non-binary integer bound tightening via piecewise linear walk.
  //
  // row 0: 3*x0 + y <= 8
  // x0 binary, y integer in [0, 10]
  // implication: x0=0 => y >= 9
  //
  // standard propagation: min_activity=0, threshold=8, so y <= 8.
  // VI-aware: for 5 < y < 9, the implication makes x0=0 infeasible, and
  // x0=1 violates the row (3 + y > 8). for y >= 9 the row is violated as
  // well. hence y <= 5, which the piecewise walk from the upper endpoint
  // finds. the implication also lifts x0.weightLower to 9 > 8, so x0 is fixed
  // to 1 and removed.
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {1, 10};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {8};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 2;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {3, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // manually add implication: x0=0 => y >= 9 (no VLB needed)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 1,
                              HighsImplications::Implication{9.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0 is fixed and removed, so y is at presolved index 0. the piecewise
  // walk should tighten y's upper bound from 10 to 5
  REQUIRE(mipsolver.model_->num_col_ == 1);
  REQUIRE(mipsolver.model_->col_upper_[0] <= 5.0 + 1e-6);
  REQUIRE(mipsolver.model_->col_upper_[0] >= 5.0 - 1e-6);

  // x0 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {5.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-paper-example-3-5", "[highs_test_presolve_rules]") {
  // Chen et al. 2026, Example 3: binary fixing via an implied bound and
  // cliques on a binary with zero coefficient.
  //
  // row 0: x2 + 0.9*x3 + 0.5*x5 <= 2
  // x1(col0), x2(col1), x3(col2), x4(col3) binary; x5(col4) continuous [0, 3]
  // x1 and x4 have zero coefficient in the row.
  //
  // VIs (24a)-(24h); binary-binary VIs are cliques:
  //   (24a) x1=0 => x2>=1        clique {x1=0, x2=0}
  //   (24b) x1=0 => x3>=1        clique {x1=0, x3=0}
  //   (24c)-(24d) x1=0 => 0.4<=x5<=0.5
  //   (24e) x1=1 => x4>=1        clique {x1=1, x4=0}
  //   (24f) x5>0.5 => x2=0       x2=1 => x5<=0.5
  //   (24g) x5<1 => x2=1         x2=0 => x5>=1
  //   (24h) x5>2 => x3=1         x3=0 => x5<=2
  //
  // the implied bound x5>=0.4 lifts x1.weightLower by 0.5*0.4 = 0.2, clique
  // propagation adds |1| + |0.9| = 1.9, so x1.weightLower = 2.1 > 2 and x1 is
  // fixed to 1. the paper's bound x5 <= 2.2 (Examples 4-5) is not derived:
  // bounds of continuous columns are not tightened, only fixed. the piecewise
  // walk is tested in test-impl-aware-piecewise-walk. fixing x1 = 1 also
  // fixes x4 = 1 via the clique from (24e), as in the paper. the clique table
  // only changes x4's bounds; the fixed column would be removed by later
  // presolve, which the rule test does not run.
  HighsLp lp;
  lp.num_col_ = 5;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 3};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {2};
  // x1, x4 not in row; x2(col1)=1, x3(col2)=0.9, x5(col4)=0.5
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 5;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 0, 1, 2, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {1, 0.9, 0.5};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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
  // fixings are only propagated through the clique table once cliques have
  // been extracted (by probing)
  mipsolver.mipdata_->cliquesExtracted = true;

  // cliques from (24a), (24b), (24e)
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  HighsCliqueTable::CliqueVar clq1[] = {{0, 0}, {1, 0}};
  cliquetable.doAddClique(clq1, 2);
  HighsCliqueTable::CliqueVar clq2[] = {{0, 0}, {2, 0}};
  cliquetable.doAddClique(clq2, 2);
  HighsCliqueTable::CliqueVar clq3[] = {{0, 1}, {3, 0}};
  cliquetable.doAddClique(clq3, 2);

  // implications on x5 from (24c)-(24d), (24f)-(24h)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 4,
                              HighsImplications::Implication{0.4, 0.5});
  implications.addImplication(1, 1, 4,
                              HighsImplications::Implication{-kHighsInf, 0.5});
  implications.addImplication(1, 0, 4,
                              HighsImplications::Implication{1.0, kHighsInf});
  implications.addImplication(2, 0, 4,
                              HighsImplications::Implication{-kHighsInf, 2.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x1 is fixed and removed; presolved cols are x2(0), x3(1), x4(2), x5(3)
  REQUIRE(mipsolver.model_->num_col_ == 4);

  // x4 is fixed to 1
  REQUIRE(mipsolver.model_->col_lower_[2] == 1.0);
  REQUIRE(mipsolver.model_->col_upper_[2] == 1.0);

  // x1 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {0.0, 0.0, 1.0, 1.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-paper-example-3-5-geq",
          "[highs_test_presolve_rules]") {
  // Chen et al. 2026, Example 3 with the row negated into a >= row, so that
  // the row is propagated in the negated direction.
  //
  // row 0: -x2 - 0.9*x3 - 0.5*x5 >= -2
  // x1(col0), x2(col1), x3(col2), x4(col3) binary; x5(col4) continuous [0, 3]
  // x1 and x4 have zero coefficient in the row.
  //
  // cliques and implications as in test-impl-aware-paper-example-3-5.
  // expected result is the same: x1 = 1 and x4 = 1.
  HighsLp lp;
  lp.num_col_ = 5;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 3};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kContinuous};
  lp.row_lower_ = {-2};
  lp.row_upper_ = {kHighsInf};
  // x1, x4 not in row; x2(col1)=-1, x3(col2)=-0.9, x5(col4)=-0.5
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 5;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 0, 1, 2, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {-1, -0.9, -0.5};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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
  // fixings are only propagated through the clique table once cliques have
  // been extracted (by probing)
  mipsolver.mipdata_->cliquesExtracted = true;

  // cliques from (24a), (24b), (24e)
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  HighsCliqueTable::CliqueVar clq1[] = {{0, 0}, {1, 0}};
  cliquetable.doAddClique(clq1, 2);
  HighsCliqueTable::CliqueVar clq2[] = {{0, 0}, {2, 0}};
  cliquetable.doAddClique(clq2, 2);
  HighsCliqueTable::CliqueVar clq3[] = {{0, 1}, {3, 0}};
  cliquetable.doAddClique(clq3, 2);

  // implications on x5 from (24c)-(24d), (24f)-(24h)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 4,
                              HighsImplications::Implication{0.4, 0.5});
  implications.addImplication(1, 1, 4,
                              HighsImplications::Implication{-kHighsInf, 0.5});
  implications.addImplication(1, 0, 4,
                              HighsImplications::Implication{1.0, kHighsInf});
  implications.addImplication(2, 0, 4,
                              HighsImplications::Implication{-kHighsInf, 2.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x1 is fixed and removed; presolved cols are x2(0), x3(1), x4(2), x5(3)
  REQUIRE(mipsolver.model_->num_col_ == 4);

  // x4 is fixed to 1
  REQUIRE(mipsolver.model_->col_lower_[2] == 1.0);
  REQUIRE(mipsolver.model_->col_upper_[2] == 1.0);

  // x1 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {0.0, 0.0, 1.0, 1.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-fixed-col", "[highs_test_presolve_rules]") {
  // fixed column that has not been removed before the impl-aware propagation
  // runs (rule test mode skips the initial row and column presolve).
  //
  // row 0: 3*x0 + y + 2*z <= 10
  // x0 binary, y integer in [0, 10], z integer fixed at 1
  // implication: x0=0 => y >= 9
  //
  // the fixed column contributes 2 to the minimum activity, so the threshold
  // is 8 and the result equals test-impl-aware-nonbinary-tightening: x0 is
  // fixed to 1 (x0.weightLower = 9 > 8) and y <= 5. ignoring z's
  // contribution would give threshold 10: x0 would not be fixed and y would
  // not be tightened (y = 10 with x0 = 0 has weight 10 <= 10). z itself must
  // stay fixed at 1.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0};
  lp.col_lower_ = {0, 0, 1};
  lp.col_upper_ = {1, 10, 1};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {10};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {3, 1, 2};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // implication: x0=0 => y >= 9
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 1,
                              HighsImplications::Implication{9.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0 is fixed and removed, so y and z are at presolved indices 0 and 1
  REQUIRE(mipsolver.model_->num_col_ == 2);
  REQUIRE(mipsolver.model_->col_upper_[0] <= 5.0 + 1e-6);
  REQUIRE(mipsolver.model_->col_upper_[0] >= 5.0 - 1e-6);
  REQUIRE(mipsolver.model_->col_lower_[1] == 1.0);
  REQUIRE(mipsolver.model_->col_upper_[1] == 1.0);

  // x0 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {5.0, 1.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-objective-cutoff", "[highs_test_presolve_rules]") {
  // impl-aware propagation of the objective cutoff, including the objective
  // offset and a fixed column that has not been removed yet.
  //
  // min 3*x0 + y + 2*z + 4
  // row 0: x0 + y <= 100 (redundant)
  // x0 binary, y integer in [0, 10], z integer fixed at 1
  // implication: x0=0 => y >= 9
  //
  // the objective bound 14 gives upper_limit = 14 in the original frame, so
  // 3*x0 + y + 2*z <= 14 - 4 = 10. z contributes 2 to the minimum, so the
  // threshold is 8 and, as in test-impl-aware-nonbinary-tightening, x0 is
  // fixed to 1 and y <= 5. ignoring the offset (threshold 12) or z
  // (threshold 10) leaves x0 unfixed and y not tightened.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 4;
  lp.col_cost_ = {3, 1, 2};
  lp.col_lower_ = {0, 0, 1};
  lp.col_upper_ = {1, 10, 1};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {100};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
  highs.setOptionValue("objective_bound", 14.0);
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
  REQUIRE(mipsolver.mipdata_->upper_limit == 14.0);

  // implication: x0=0 => y >= 9
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 1,
                              HighsImplications::Implication{9.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0 is fixed and removed, so y and z are at presolved indices 0 and 1
  REQUIRE(mipsolver.model_->num_col_ == 2);
  REQUIRE(mipsolver.model_->col_upper_[0] <= 5.0 + 1e-6);
  REQUIRE(mipsolver.model_->col_upper_[0] >= 5.0 - 1e-6);
  REQUIRE(mipsolver.model_->col_lower_[1] == 1.0);
  REQUIRE(mipsolver.model_->col_upper_[1] == 1.0);

  // x0 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {5.0, 1.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-binary-infeasible", "[highs_test_presolve_rules]") {
  // infeasibility detection when neither value of a binary satisfies the row.
  //
  // row 0: 3*x0 + y + z <= 8
  // x0 binary, y and z integer in [0, 10]
  // VLBs: y >= -9*x0 + 9 (x0=0 => y >= 9), z >= 6*x0 (x0=1 => z >= 6)
  //
  // x0.weightLower = 9 > 8 and x0.weightUpper = 3 + 6 = 9 > 8, so the
  // problem is infeasible.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {1, 10, 10};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {8};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {3, 1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // VLBs: y >= -9*x0 + 9, z >= 6*x0
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addVLB(1, 0, -9.0, 9.0);
  implications.addVLB(2, 0, 6.0, 0.0);

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kInfeasible);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-infeasible-literal", "[highs_test_presolve_rules]") {
  // binary fixing from a literal whose implied bounds contradict the bounds of
  // a non-binary, independently of the row's threshold.
  //
  // row 0: x0 + y <= 10
  // x0 binary; y continuous in [0, 10]; z continuous in [0, 4], not in any row
  // implication: x0=0 => z >= 6
  //
  // the implication does not lift x0 since z has a zero coefficient, so
  // x0.weightLower = 0 <= 10. but z >= 6 contradicts z <= 4, so x0=0 is
  // infeasible and x0 is fixed to 1.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {1, 10, 4};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kContinuous,
                     HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {10};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // x0=0 => z >= 6
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 2,
                              HighsImplications::Implication{6.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0 is fixed and removed
  REQUIRE(mipsolver.model_->num_col_ == 2);

  // x0 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value.assign(mipsolver.model_->num_col_, 0.0);
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-paper-example-6", "[highs_test_presolve_rules]") {
  // Chen et al. 2026, Example 6: non-binary variable not in the row has its
  // lower bound tightened via VIs to binaries in the row.
  //
  // row 0: x1 + x2 + x3 + x4 + 0.1*x5 <= 2   (54a)
  // x1-x5(col0-4) binary; x6(col5) integer [0, 4], not in row
  // x6 is continuous in the paper, but continuous bounds are only fixed, not
  // tightened. the bound lies on a breakpoint, so the result is the same.
  // VIs (55a)-(55e) as implications: x1=0 => x6>=3, x2=0 => x6>=3,
  // x3=0 => x6>=2, x4=0 => x6>=2, x5=0 => x6<=1
  //
  // weight of the row as a function of x6 = d: 4 on [0, 1], 4.1 on (1, 2),
  // 2.1 on [2, 3), 0.1 at 3. the smallest d with weight <= 2 is 3, so
  // x6 >= 3 (57). without (55e) the weight at d = 2 is 2, giving x6 >= 2 as
  // in Achterberg et al. (2013) (56).
  HighsLp lp;
  lp.num_col_ = 6;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 1, 4};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {2};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 6;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 4, 5, 5};
  lp.a_matrix_.index_ = {0, 0, 0, 0, 0};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 0.1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // x1=0 => x6>=3, x2=0 => x6>=3, x3=0 => x6>=2, x4=0 => x6>=2,
  // x5=0 => x6<=1
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 5,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(1, 0, 5,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(2, 0, 5,
                              HighsImplications::Implication{2.0, kHighsInf});
  implications.addImplication(3, 0, 5,
                              HighsImplications::Implication{2.0, kHighsInf});
  implications.addImplication(4, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 1.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x6 lower bound should be tightened from 0 to 3
  REQUIRE(mipsolver.model_->col_lower_[5] >= 3.0 - 1e-6);
  REQUIRE(mipsolver.model_->col_lower_[5] <= 3.0 + 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-paper-example-7", "[highs_test_presolve_rules]") {
  // Chen et al. 2026, Example 7: non-binary lower bound tightening
  // with multiple breakpoints at the same value.
  //
  // row 0: x1 + x2 + x3 + x4 + 0.1*x5 + 0.2*x6 <= 2.2   (58)
  // x1-x5(col0-4) binary, x6(col5) integer [0, 4]
  // x6 is continuous in the paper, but continuous bounds are only fixed, not
  // tightened. the bound lies on a breakpoint, so the result is the same.
  // VIs (55a)-(55d) as implications: x1=0 => x6>=3, x2=0 => x6>=3,
  //                                  x3=0 => x6>=2, x4=0 => x6>=2
  //
  // x6.weightLower = 4 (sum of excesses from 4 lower-type breakpoints)
  // threshold = 2.2, weightLower > threshold -> tighten lower bound.
  // walk from lb=0: breakpoints at {2,2,3,3} deactivate weight, w(d) =
  // 2 + 0.2*d > 2.2 on [2, 3) and w(3) = 0.6. paper result (59): x6 >= 3.
  // relaxing x6 from the row as in Achterberg et al. (2013) gives x6 >= 2.
  HighsLp lp;
  lp.num_col_ = 6;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 1, 4};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {2.2};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 6;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 4, 5, 6};
  lp.a_matrix_.index_ = {0, 0, 0, 0, 0, 0};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 0.1, 0.2};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // x1=0 => x6>=3, x2=0 => x6>=3, x3=0 => x6>=2, x4=0 => x6>=2
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 5,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(1, 0, 5,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(2, 0, 5,
                              HighsImplications::Implication{2.0, kHighsInf});
  implications.addImplication(3, 0, 5,
                              HighsImplications::Implication{2.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x6 lower bound should be tightened from 0 to 3
  REQUIRE(mipsolver.model_->col_lower_[5] >= 3.0 - 1e-6);
  REQUIRE(mipsolver.model_->col_lower_[5] <= 3.0 + 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-clique-extraction-complemented",
          "[highs_test_presolve_rules]") {
  // test-impl-aware-clique-extraction with y1 complemented: y1' = 10 - y1.
  // exercises a negative non-binary coefficient and a VUB (instead of a VLB).
  //
  // row 0: 3*x0 + 3*x1 - y1' + y2 <= 2
  // row 1: -5*x0 - y1' >= -10   (from y1 >= 5*x0, substituting y1 = 10 - y1')
  // row 2: -5*x1 + y2 >= 0      (unchanged)
  // x0, x1 binary, y1', y2 continuous in [0, 10]
  //
  // VUB on y1': y1' <= -5*x0 + 10, so x0=1 => y1' <= 5 (i.e. y1 >= 5)
  // VLB on y2: y2 >= 5*x1 (unchanged)
  // threshold = 2 - (-10) = 12 and the combined implied weight is 16 as in
  // the uncomplemented test, so {x0=1, x1=1} is a clique.
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 3;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 10, 10};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kContinuous, HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf, -10, 0};
  lp.row_upper_ = {2, kHighsInf, kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 3;
  lp.a_matrix_.start_ = {0, 2, 4, 6, 8};
  lp.a_matrix_.index_ = {0, 1, 0, 2, 0, 1, 0, 2};
  lp.a_matrix_.value_ = {3, -5, 3, -5, -1, -1, 1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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
  // cliques are only added once cliques have been extracted (by probing)
  mipsolver.mipdata_->cliquesExtracted = true;

  // manually populate VUB, VLB and implications (normally done by probing)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  // x0=1 implies y1' <= 5 (VUB: y1' <= -5*x0 + 10, from y1 >= 5*x0)
  implications.addVUB(2, 0, -5.0, 10.0, 1);
  implications.addImplication(0, 1, 2,
                              HighsImplications::Implication{-kHighsInf, 5.0});
  // x1=1 implies y2 >= 5 (VLB: y2 >= 5*x1, unchanged)
  implications.addVLB(3, 1, 5.0, 0.0, 2);
  implications.addImplication(1, 1, 3,
                              HighsImplications::Implication{5.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0=1 and x1=1 can't coexist (same clique as original, via VUB path)
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  REQUIRE(cliquetable.haveCommonClique({0, 1}, {1, 1}));

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-conflict-clique-complemented",
          "[highs_test_presolve_rules]") {
  // test-impl-aware-conflict-clique with x1 complemented: x1 here is 1 - x1
  // of the original test.
  //
  // row 0: 4*x0 - 4*x1 + y <= 8
  // x0, x1 binary; y, z continuous in [0, 10]; z is not in any row
  // VLBs: y >= 3*x0, y >= -3*x1 + 3
  // implications: x0=1 => z >= 6, x1=0 => z <= 4
  //
  // threshold = 8 - (-4) = 12, x0.weightUpper = x1.weightLower = 7.
  // implied activity of the pair (x0=1, x1=0) is 4 + 4 + 3 = 11 <= 12, but
  // the implied bounds on z contradict, so (x0=1, x1=0) form a clique.
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 10, 10};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kContinuous, HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {8};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {4, -4, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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
  // cliques are only added once cliques have been extracted (by probing)
  mipsolver.mipdata_->cliquesExtracted = true;

  // VLBs: y >= 3*x0, y >= -3*x1 + 3
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addVLB(2, 0, 3.0, 0.0);
  implications.addVLB(2, 1, -3.0, 3.0);
  // x0=1 => z >= 6, x1=0 => z <= 4
  implications.addImplication(0, 1, 3,
                              HighsImplications::Implication{6.0, kHighsInf});
  implications.addImplication(1, 0, 3,
                              HighsImplications::Implication{-kHighsInf, 4.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0=1 and x1=0 can't coexist
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  REQUIRE(cliquetable.haveCommonClique({0, 1}, {1, 0}));

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-nonbinary-tightening-complemented",
          "[highs_test_presolve_rules]") {
  // test-impl-aware-nonbinary-tightening with y complemented: y' = 10 - y.
  // exercises a negative non-binary coefficient and lower bound tightening.
  //
  // row 0: 3*x0 - y' <= -2
  // x0 binary, y' integer in [0, 10]
  // implication: x0=0 => y' <= 1  (from x0=0 => y >= 9)
  //
  // standard propagation: min_activity=-10, threshold=8, so y' >= 2.
  // VI-aware: for 1 < y' < 5, the implication makes x0=0 infeasible, and
  // x0=1 violates the row (3 - y' > -2). for y' <= 1 the row is violated as
  // well. hence y' >= 5 (original: y <= 5), which the piecewise walk from the
  // lower endpoint finds. the implication also lifts x0.weightLower to
  // 9 > 8, so x0 is fixed to 1 and removed.
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {1, 10};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {-2};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 2;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {3, -1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // x0=0 => y' <= 1 (from x0=0 => y >= 9, complemented)
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 1,
                              HighsImplications::Implication{-kHighsInf, 1.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x0 is fixed and removed, so y' is at presolved index 0. y' lower bound
  // tightened from 0 to 5 (original: y upper bound 10 -> 5)
  REQUIRE(mipsolver.model_->num_col_ == 1);
  REQUIRE(mipsolver.model_->col_lower_[0] >= 5.0 - 1e-6);
  REQUIRE(mipsolver.model_->col_lower_[0] <= 5.0 + 1e-6);

  // x0 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {5.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-paper-example-3-5-complemented",
          "[highs_test_presolve_rules]") {
  // Chen et al. 2026, Example 3 with x5 complemented: x5' = 3 - x5. exercises
  // a negative non-binary coefficient in the lift.
  //
  // row 0: x2 + 0.9*x3 - 0.5*x5' <= 0.5
  // x1(col0), x2(col1), x3(col2), x4(col3) binary; x5'(col4) continuous [0, 3]
  // x1 and x4 have zero coefficient in the row.
  //
  // cliques as in test-impl-aware-paper-example-3-5, implications on x5'
  // complemented: x1=0 => 2.5<=x5'<=2.6, x2=1 => x5'>=2.5, x2=0 => x5'<=2,
  // x3=0 => x5'>=1. x5'<=2.6 lifts x1.weightLower by -0.5*(2.6-3) = 0.2, as
  // in the uncomplemented example. expected result is the same: x1 = 1 and x4
  // = 1.
  HighsLp lp;
  lp.num_col_ = 5;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 3};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {0.5};
  // x1, x4 not in row; x2(col1)=1, x3(col2)=0.9, x5'(col4)=-0.5
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 5;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 0, 1, 2, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {1, 0.9, -0.5};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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
  // fixings are only propagated through the clique table once cliques have
  // been extracted (by probing)
  mipsolver.mipdata_->cliquesExtracted = true;

  // cliques from (24a), (24b), (24e)
  HighsCliqueTable& cliquetable = mipsolver.mipdata_->cliquetable;
  HighsCliqueTable::CliqueVar clq1[] = {{0, 0}, {1, 0}};
  cliquetable.doAddClique(clq1, 2);
  HighsCliqueTable::CliqueVar clq2[] = {{0, 0}, {2, 0}};
  cliquetable.doAddClique(clq2, 2);
  HighsCliqueTable::CliqueVar clq3[] = {{0, 1}, {3, 0}};
  cliquetable.doAddClique(clq3, 2);

  // implications on x5' from (24c)-(24d), (24f)-(24h), complemented
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 4,
                              HighsImplications::Implication{2.5, 2.6});
  implications.addImplication(1, 1, 4,
                              HighsImplications::Implication{2.5, kHighsInf});
  implications.addImplication(1, 0, 4,
                              HighsImplications::Implication{-kHighsInf, 2.0});
  implications.addImplication(2, 0, 4,
                              HighsImplications::Implication{1.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x1 is fixed and removed; presolved cols are x2(0), x3(1), x4(2), x5'(3)
  REQUIRE(mipsolver.model_->num_col_ == 4);

  // x4 is fixed to 1
  REQUIRE(mipsolver.model_->col_lower_[2] == 1.0);
  REQUIRE(mipsolver.model_->col_upper_[2] == 1.0);

  // x1 is restored to 1 by postsolve
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {0.0, 0.0, 1.0, 2.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[0] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-piecewise-walk", "[highs_test_presolve_rules]") {
  // piecewise walk on an integer non-binary with a lower-type and an
  // upper-type breakpoint, built from VIs (24g) and (24h) of Chen et al. 2026.
  //
  // row 0: x2 + 0.9*x3 + 0.5*x5 <= 2
  // x2(col0), x3(col1) binary; x5(col2) integer [0, 3]
  // implications: x2=0 => x5>=1  (lower-type breakpoint at 1, excess 1)
  //               x3=0 => x5<=2  (upper-type breakpoint at 2, excess 0.9)
  //
  // implied activity: w(x5) = 0.5*x5 + [x5 < 1]*1 + [x5 > 2]*0.9, threshold 2.
  // w(3) = 2.4 > 2, so the upper walk starts at 3 and crosses the threshold
  // mid-segment: 0.5*x5 + 0.9 <= 2 gives x5 <= 2.2, floored to x5 <= 2.
  // w(0) = 1 <= 2, so the lower bound is not tightened. standard propagation
  // only gives x5 <= 4, and no binary is fixed.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {1, 1, 3};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {2};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {1, 0.9, 0.5};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 2,
                              HighsImplications::Implication{1.0, kHighsInf});
  implications.addImplication(1, 0, 2,
                              HighsImplications::Implication{-kHighsInf, 2.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  REQUIRE(mipsolver.model_->num_col_ == 3);
  REQUIRE(mipsolver.model_->col_lower_[2] == 0.0);
  REQUIRE(mipsolver.model_->col_upper_[2] == 2.0);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-piecewise-walk-complemented",
          "[highs_test_presolve_rules]") {
  // test-impl-aware-piecewise-walk with x5 complemented: x5' = 3 - x5.
  // exercises the lower walk with a negative coefficient.
  //
  // row 0: x2 + 0.9*x3 - 0.5*x5' <= 0.5
  // x2(col0), x3(col1) binary; x5'(col2) integer [0, 3]
  // implications: x2=0 => x5'<=2  (upper-type breakpoint at 2, excess 1)
  //               x3=0 => x5'>=1  (lower-type breakpoint at 1, excess 0.9)
  //
  // implied activity above the minimum: w(x5') = 0.5*(3-x5') + [x5' < 1]*0.9 +
  // [x5' > 2]*1, threshold 2. w(0) = 2.4 > 2, so the lower walk starts at 0
  // and crosses the threshold mid-segment: 2.4 - 0.5*x5' <= 2 gives
  // x5' >= 0.8, rounded up to x5' >= 1 (original: x5 <= 2). w(3) = 1 <= 2, so
  // the upper bound is not tightened.
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {1, 1, 3};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {0.5};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 3;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {1, 0.9, -0.5};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 2,
                              HighsImplications::Implication{-kHighsInf, 2.0});
  implications.addImplication(1, 0, 2,
                              HighsImplications::Implication{1.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  REQUIRE(mipsolver.model_->num_col_ == 3);
  REQUIRE(mipsolver.model_->col_lower_[2] == 1.0);
  REQUIRE(mipsolver.model_->col_upper_[2] == 3.0);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-piecewise-walk-tie", "[highs_test_presolve_rules]") {
  // lower-type and upper-type breakpoints at the same value. neither type is
  // active at the breakpoint itself, so it is an isolated feasible point that
  // both walks only find if the breakpoint that deactivates is processed
  // before the one that activates (sort tie-break: lower-type first).
  //
  // row 0: x0 + x1 + x2 + x3 <= 1.5
  // x0..x3(col0..col3) binary; y(col4) integer [0, 6], not in row
  // implications: x0=0 => y>=3, x1=0 => y>=3  (lower-type breakpoints at 3)
  //               x2=0 => y<=3, x3=0 => y<=3  (upper-type breakpoints at 3)
  //
  // implied activity: w(y) = [y < 3]*2 + [y > 3]*2, threshold 1.5.
  // w(0) = 2 > 1.5: walking right, both lower-types deactivate at 3 and
  // w(3) = 0, so y >= 3. w(6) = 2 > 1.5: walking left, both upper-types
  // deactivate at 3, so y <= 3. y is fixed to 3. with upper-types first,
  // the walk right would see 2 -> 3 -> 4 -> 3 -> 2 at y = 3 and miss the
  // feasible point.
  HighsLp lp;
  lp.num_col_ = 5;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 6};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {1.5};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 5;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 4, 4};
  lp.a_matrix_.index_ = {0, 0, 0, 0};
  lp.a_matrix_.value_ = {1, 1, 1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 4,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(1, 0, 4,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(2, 0, 4,
                              HighsImplications::Implication{-kHighsInf, 3.0});
  implications.addImplication(3, 0, 4,
                              HighsImplications::Implication{-kHighsInf, 3.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // y is fixed and removed; the binaries remain
  REQUIRE(mipsolver.model_->num_col_ == 4);

  // postsolve restores y = 3
  HighsSolution sol;
  sol.value_valid = true;
  sol.col_value = {0.0, 0.0, 0.0, 0.0};
  postsolve_stack.undoPrimal(options, sol);
  REQUIRE(sol.col_value[4] == 3.0);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-piecewise-walk-infeasible",
          "[highs_test_presolve_rules]") {
  // infeasibility detection when the implied activity of a non-binary exceeds
  // the threshold on its whole domain.
  //
  // row 0: x0 + x1 + x2 + x3 <= 1.5
  // x0..x3(col0..col3) binary; y(col4) integer [0, 6], not in row
  // implications: x0=0 => y>=3, x1=0 => y>=3  (lower-type breakpoints at 3)
  //               x2=0 => y<=2, x3=0 => y<=2  (upper-type breakpoints at 2)
  //
  // implied activity: w(y) = [y < 3]*2 + [y > 2]*2, threshold 1.5.
  // w(y) >= 2 > 1.5 for all y in [0, 6], so neither walk crosses the
  // threshold and the problem is infeasible. no binary is fixed, since y has a
  // zero coefficient and does not lift any binary.
  HighsLp lp;
  lp.num_col_ = 5;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 6};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {1.5};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 5;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 4, 4};
  lp.a_matrix_.index_ = {0, 0, 0, 0};
  lp.a_matrix_.value_ = {1, 1, 1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 4,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(1, 0, 4,
                              HighsImplications::Implication{3.0, kHighsInf});
  implications.addImplication(2, 0, 4,
                              HighsImplications::Implication{-kHighsInf, 2.0});
  implications.addImplication(3, 0, 4,
                              HighsImplications::Implication{-kHighsInf, 2.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kInfeasible);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-piecewise-walk-infeasible-upper",
          "[highs_test_presolve_rules]") {
  // infeasibility detection in the upper walk. the lower walk does not start,
  // since an upper-type breakpoint below lb is neither counted in w(lb) nor
  // activated by the lower walk. the upper walk keeps it active.
  //
  // row 0: x0 + x1 + x2 <= 1.5
  // x0..x2(col0..col2) binary; y(col3) integer [1, 3], not in row
  // implications: x0=0 => y<=0  (infeasible literal, upper-type at 0)
  //               x1=0 => y>=2  (lower-type breakpoint at 2)
  //               x2=0 => y<=1  (upper-type breakpoint at 1)
  //
  // x0=0 contradicts y >= 1, so x0 is fixed to 1. implied activity:
  // w(y) = 1 + [y < 2]*1 + [y > 1]*1 = 2 > 1.5 for all y in [1, 3], so the
  // problem is infeasible. w(lb) without x0 is 1 <= 1.5, so the lower walk is
  // skipped. w(ub) = 2 > 1.5, and the upper walk does not cross the threshold.
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 1};
  lp.col_upper_ = {1, 1, 1, 3};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {1.5};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 4;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 3};
  lp.a_matrix_.index_ = {0, 0, 0};
  lp.a_matrix_.value_ = {1, 1, 1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 3,
                              HighsImplications::Implication{-kHighsInf, 0.0});
  implications.addImplication(1, 0, 3,
                              HighsImplications::Implication{2.0, kHighsInf});
  implications.addImplication(2, 0, 3,
                              HighsImplications::Implication{-kHighsInf, 1.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kInfeasible);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-paper-example-6-complemented",
          "[highs_test_presolve_rules]") {
  // Chen et al. 2026, Example 6 with x6 complemented: x6' = 4 - x6. exercises
  // upper bound tightening of a non-binary variable not in the row.
  //
  // row 0: x1 + x2 + x3 + x4 + 0.1*x5 <= 2   (unchanged, x6' not in row)
  // x1-x5(col0-4) binary; x6'(col5) integer [0, 4], not in row
  // x6 is continuous in the paper, but continuous bounds are only fixed, not
  // tightened. the bound lies on a breakpoint, so the result is the same.
  // VIs complemented: x1=0 => x6'<=1, x2=0 => x6'<=1, x3=0 => x6'<=2,
  // x4=0 => x6'<=2, x5=0 => x6'>=3
  //
  // walk from ub=4 mirrors test-impl-aware-paper-example-6.
  // result: x6' <= 1 (original: x6 >= 3).
  HighsLp lp;
  lp.num_col_ = 6;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 1, 4};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {2};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 6;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 4, 5, 5};
  lp.a_matrix_.index_ = {0, 0, 0, 0, 0};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 0.1};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // x1=0 => x6'<=1, x2=0 => x6'<=1, x3=0 => x6'<=2, x4=0 => x6'<=2,
  // x5=0 => x6'>=3
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 1.0});
  implications.addImplication(1, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 1.0});
  implications.addImplication(2, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 2.0});
  implications.addImplication(3, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 2.0});
  implications.addImplication(4, 0, 5,
                              HighsImplications::Implication{3.0, kHighsInf});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x6' upper bound tightened from 4 to 1 (original: x6 lower bound 0 -> 3)
  REQUIRE(mipsolver.model_->col_upper_[5] <= 1.0 + 1e-6);
  REQUIRE(mipsolver.model_->col_upper_[5] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("test-impl-aware-paper-example-7-complemented",
          "[highs_test_presolve_rules]") {
  // Chen et al. 2026, Example 7 with x6 complemented: x6' = 4 - x6. exercises
  // a negative non-binary coefficient with multiple upper-type breakpoints at
  // the same value.
  //
  // row 0: x1 + x2 + x3 + x4 + 0.1*x5 - 0.2*x6' <= 1.4
  // x1-x5(col0-4) binary, x6'(col5) integer [0, 4]
  // x6 is continuous in the paper, but continuous bounds are only fixed, not
  // tightened. the bound lies on a breakpoint, so the result is the same.
  // VIs complemented: x1=0 => x6'<=1, x2=0 => x6'<=1,
  //                   x3=0 => x6'<=2, x4=0 => x6'<=2
  //
  // x6'.weightUpper = 4 (sum of excesses from 4 upper-type breakpoints)
  // threshold = 1.4 - (-0.8) = 2.2, weightUpper > threshold -> tighten upper
  // bound. walk from ub=4 mirrors test-impl-aware-paper-example-7.
  // result: x6' <= 1 (original: x6 >= 3).
  HighsLp lp;
  lp.num_col_ = 6;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {0, 0, 0, 0, 0, 0};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 1, 4};
  lp.integrality_ = {HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger,
                     HighsVarType::kInteger, HighsVarType::kInteger};
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {1.4};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  lp.a_matrix_.num_col_ = 6;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.start_ = {0, 1, 2, 3, 4, 5, 6};
  lp.a_matrix_.index_ = {0, 0, 0, 0, 0, 0};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 0.1, -0.2};

  highs::parallel::initialize_scheduler(1);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve_rule_test",
                       kPresolveRuleImplAwareConstrPropagation);
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

  // x1=0 => x6'<=1, x2=0 => x6'<=1, x3=0 => x6'<=2, x4=0 => x6'<=2
  HighsImplications& implications = mipsolver.mipdata_->implications;
  implications.addImplication(0, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 1.0});
  implications.addImplication(1, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 1.0});
  implications.addImplication(2, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 2.0});
  implications.addImplication(3, 0, 5,
                              HighsImplications::Implication{-kHighsInf, 2.0});

  presolve::HighsPostsolveStack& postsolve_stack =
      mipsolver.mipdata_->postSolveStack;

  presolve::HPresolve presolve;
  presolve.setInput(mipsolver, -1);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  REQUIRE(status == HighsModelStatus::kNotset);

  // x6' upper bound tightened from 4 to 1 (original: x6 lower bound 0 -> 3)
  REQUIRE(mipsolver.model_->col_upper_[5] <= 1.0 + 1e-6);
  REQUIRE(mipsolver.model_->col_upper_[5] >= 1.0 - 1e-6);

  HighsTaskExecutor::shutdown(true);
}

TEST_CASE("issue-3364", "[highs_test_presolve_rules]") {
  // Row i (2) is an equation with two parallel candidates in its bucket:
  // A (0), which carries two integer singletons, and B (1), which is exactly
  // parallel. A and B are never compared with each other: neither is an
  // equation, so the singleton rules skip the pair. Visiting A adds row i into
  // A, which leaves the integer inequality 2 s1 + 3 s2 <= 7, and rowPresolve
  // on it calls storeRow(A). If the scan then went on to B, row i would be
  // deleted with A's positions, leaving i's own nonzeros linked to a deleted
  // row, which shrinkProblem would renumber to row -1.
  // Rows D (3) and E (4) keep x0, x1 and x3 from being parallel columns; s1
  // and s2 are integers whose coefficients have a fractional ratio, so they
  // are not merged as parallel columns either.
  //
  //   cols: x0 x1 x3 s1 s2 x4 x5
  //   A:  x0 + 2 x1 + 3 x3 + 2 s1 + 3 s2  <= 13
  //   B:  x0 + 2 x1 + 3 x3                 <= 8
  //   i:  x0 + 2 x1 + 3 x3                  = 6
  //   D:  x0                     + x4      <= 7
  //   E:       x1                     + x5 <= 7
  HighsLp lp;
  lp.num_col_ = 7;
  lp.num_row_ = 5;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_ = {-1, -1, -1, 1, 1, -1, -1};
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, 10);
  lp.integrality_ = {HighsVarType::kContinuous, HighsVarType::kContinuous,
                     HighsVarType::kContinuous, HighsVarType::kInteger,
                     HighsVarType::kInteger,    HighsVarType::kContinuous,
                     HighsVarType::kContinuous};
  lp.row_lower_ = {-kHighsInf, -kHighsInf, 6, -kHighsInf, -kHighsInf};
  lp.row_upper_ = {13, 8, 6, 7, 7};
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 4, 8, 11, 12, 13, 14, 15};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 0, 1, 2, 4, 0, 1, 2, 0, 0, 3, 4};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 2, 2, 2, 1, 3, 3, 3, 2, 3, 1, 1};
  lp.col_names_ = {"x0", "x1", "x3", "s1", "s2", "x4", "x5"};
  lp.row_names_ = {"A", "B", "i", "D", "E"};

  HighsOptions options;
  options.presolve_rule_test = kPresolveRuleParallelRowsAndCols;
  // Rows with singletons are only considered when LP presolve need not
  // maintain a basis for postsolve
  options.solver = kIpmString;
  options.run_crossover = kHighsOffString;
  options.output_flag = dev_run;

  HighsTimer timer;
  timer.start();

  presolve::HighsPostsolveStack postsolve_stack;
  postsolve_stack.initializeIndexMaps(lp.num_row_, lp.num_col_);

  presolve::HPresolve presolve;
  presolve.setInput(lp, options, -1, &timer);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  //  timer.stop();
  REQUIRE(status == HighsModelStatus::kNotset);
  // Every nonzero of the presolved model must belong to a row of it
  REQUIRE(lp.a_matrix_.isColwise());
  for (const HighsInt row : lp.a_matrix_.index_) {
    REQUIRE(row >= 0);
    REQUIRE(row < lp.num_row_);
  }
}

TEST_CASE("test-parallel-rows-singleton-count", "[highs_test_presolve_rules]") {
  // rows C (0), D (1) and i (2) are parallel except for one singleton each,
  // and row i is an equation. adding row i into its first candidate (C or D)
  // leaves u - s <= -10 (or t - s <= -10), which is a forcing row, so s is
  // fixed at its upper bound and removed from row i. row i is then the
  // equation x + 2 y = -6 without singletons, and adding it into the other
  // candidate leaves a singleton row. this second reduction is only found if
  // the scan continues and the number of singletons in row i is up to date.
  // rows F (3) and G (4) keep x and y from being parallel columns.
  //
  //   cols: x y s u t w1 w2
  //   C:  x + 2 y     + u             <= -6
  //   D:  x + 2 y         + t         <= -6
  //   i:  x + 2 y + s                  = 4
  //   F:  x                   + w1    <= 7
  //   G:      y                   + w2 <= 7
  HighsLp lp;
  lp.num_col_ = 7;
  lp.num_row_ = 5;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_cost_.assign(lp.num_col_, 1);
  lp.col_lower_ = {-10, -10, 0, 0, 0, 0, 0};
  lp.col_upper_.assign(lp.num_col_, 10);
  lp.row_lower_ = {-kHighsInf, -kHighsInf, 4, -kHighsInf, -kHighsInf};
  lp.row_upper_ = {-6, -6, 4, 7, 7};
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 4, 8, 9, 10, 11, 12, 13};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 0, 1, 2, 4, 2, 0, 1, 3, 4};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1};

  HighsOptions options;
  options.presolve_rule_test = kPresolveRuleParallelRowsAndCols;
  // rows with singletons are only considered when LP presolve need not
  // maintain a basis for postsolve
  options.solver = kIpmString;
  options.run_crossover = kHighsOffString;
  options.output_flag = dev_run;

  HighsTimer timer;
  timer.start();

  presolve::HighsPostsolveStack postsolve_stack;
  postsolve_stack.initializeIndexMaps(lp.num_row_, lp.num_col_);

  presolve::HPresolve presolve;
  presolve.setInput(lp, options, -1, &timer);
  REQUIRE(presolve.okSetupPresolveDataStructures());
  HighsModelStatus status = presolve.run(postsolve_stack);
  timer.stop();
  REQUIRE(status == HighsModelStatus::kNotset);
  // rows C and D must have been removed, leaving rows i, F and G
  REQUIRE(lp.num_row_ == 3);
  REQUIRE(postsolve_stack.getOrigRowIndex(0) == 2);
  REQUIRE(postsolve_stack.getOrigRowIndex(1) == 3);
  REQUIRE(postsolve_stack.getOrigRowIndex(2) == 4);
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
