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
  REQUIRE(lp.row_upper_[0] == -1.0);
  REQUIRE(lp.row_lower_[0] == -kHighsInf);

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
