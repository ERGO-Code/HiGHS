#include <numeric>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
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
