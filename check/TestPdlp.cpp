#include "HCheckConfig.h"
#include "HConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;
const double double_equal_tolerance = 1e-3;
const double kkt_tolerance = 1e-4;
#ifdef CUPDLP_CPU
TEST_CASE("pdlp-distillation-lp", "[pdlp]") {
  SpecialLps special_lps;
  HighsLp lp;

  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(lp, require_model_status, optimal_objective);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsInfo& info = highs.getInfo();
  const HighsOptions& options = highs.getOptions();
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = HighsStatus::kOk;
  bool optimal = true;

  const HighsInt pdlp_iteration_count_optimal = 160;
  run_status = highs.run();
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
  if (optimal) {
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  } else {
    REQUIRE(run_status == HighsStatus::kWarning);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);
  }

  HighsInt pdlp_iteration_count = highs.getInfo().pdlp_iteration_count;
  REQUIRE(pdlp_iteration_count > 0);
  REQUIRE(pdlp_iteration_count == pdlp_iteration_count_optimal);

  // Now run with half the iteration count as the limit to test
  // iteration limit termination
  //
  // Now that PDLP hot starts, have to clear the solution
  //
  highs.clearSolver();

  highs.setOptionValue("pdlp_iteration_limit",
                       pdlp_iteration_count_optimal / 2);
  run_status = highs.run();

  REQUIRE(run_status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kIterationLimit);
  pdlp_iteration_count = highs.getInfo().pdlp_iteration_count;
  REQUIRE(pdlp_iteration_count > 0);
  REQUIRE(pdlp_iteration_count == (pdlp_iteration_count_optimal / 2) - 1);

  highs.resetGlobalScheduler(true);
}
#else
// CUPDLP_GPU
TEST_CASE("pdlp-distillation-lp", "[pdlp]") {
  SpecialLps special_lps;
  HighsLp lp;

  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.distillationLp(lp, require_model_status, optimal_objective);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsInfo& info = highs.getInfo();
  const HighsOptions& options = highs.getOptions();
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = HighsStatus::kOk;
  bool optimal = true;

  run_status = highs.run();
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
  if (optimal) {
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  } else {
    REQUIRE(run_status == HighsStatus::kWarning);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);
  }

  // IG todo add iteration count
  // HighsInt pdlp_iteration_count = highs.getInfo().pdlp_iteration_count;
  // REQUIRE(pdlp_iteration_count > 0);
  // REQUIRE(pdlp_iteration_count == 160);

  // Now run with half the iteration count as the limit to test
  // iteration limit termination

  // highs.setOptionValue("pdlp_iteration_limit", pdlp_iteration_count / 2);
  // run_status = highs.run();

  // REQUIRE(run_status == HighsStatus::kWarning);
  // REQUIRE(highs.getModelStatus() == HighsModelStatus::kIterationLimit);
  // pdlp_iteration_count = highs.getInfo().pdlp_iteration_count;
  // REQUIRE(pdlp_iteration_count > 0);
  // REQUIRE(pdlp_iteration_count == 79);

  highs.resetGlobalScheduler(true);
}
#endif

TEST_CASE("pdlp-3d-lp", "[pdlp]") {
  SpecialLps special_lps;
  HighsLp lp;

  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.ThreeDLp(lp, require_model_status, optimal_objective);

  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsInfo& info = highs.getInfo();
  const HighsOptions& options = highs.getOptions();
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = highs.run();
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
  const bool optimal = true;
  if (optimal) {
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  } else {
    REQUIRE(run_status == HighsStatus::kWarning);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);
  }

  highs.resetGlobalScheduler(true);
}

TEST_CASE("pdlp-boxed-row-lp", "[pdlp]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 2;
  lp.col_cost_ = {-1, -2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, 6};
  lp.row_lower_ = {3, -4};
  lp.row_upper_ = {10, 2};
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {0, 1, 0, 1};
  lp.a_matrix_.value_ = {1, 1, 1, -1};
  double optimal_objective = -16;
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsInfo& info = highs.getInfo();
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = highs.run();
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
  const bool optimal = true;
  if (optimal) {
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
  } else {
    REQUIRE(run_status == HighsStatus::kWarning);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnknown);
  }

  highs.resetGlobalScheduler(true);
}

TEST_CASE("pdlp-infeasible-lp", "[pdlp]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {-1, -2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {-1};
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("kkt_tolerance", kkt_tolerance);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnboundedOrInfeasible);

  highs.resetGlobalScheduler(true);
}

TEST_CASE("pdlp-unbounded-lp", "[pdlp]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {-1, -2};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {1};
  lp.row_upper_ = {inf};
  lp.a_matrix_.start_ = {0, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1, 1};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
  highs.setOptionValue("solver", kPdlpString);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("kkt_tolerance", kkt_tolerance);
  REQUIRE(highs.run() == HighsStatus::kOk);
  if (dev_run) highs.writeSolution("", 1);
  const bool not_unbounded = false;
  if (not_unbounded) {
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnboundedOrInfeasible);
  } else {
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
  }

  highs.resetGlobalScheduler(true);
}

void pdlpRestart(const std::string& model) {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  REQUIRE(h.readModel(model_file) == HighsStatus::kOk);
  h.setOptionValue("solver", kPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = h.run();
  const bool was_optimal = h.getModelStatus() == HighsModelStatus::kOptimal;
  h.setOptionValue("presolve", kHighsOffString);
  run_status = h.run();
}

TEST_CASE("pdlp-restart", "[pdlp]") {
  pdlpRestart("adlittle");
  //  pdlpRestart("shell");
  //  pdlpRestart("25fv47");
}

TEST_CASE("pdlp-restart-lp", "[pdlp]") {
  // LP with examples of, respectively, LB, EQ, BX and UB constraints
  // to test conversion of the incumbent solution into the cuPDLP-C
  // problem space
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 4;
  lp.col_cost_ = {1, 3, 5};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {kHighsInf, kHighsInf, kHighsInf};
  lp.row_lower_ = {1, 3, 2, -kHighsInf};
  lp.row_upper_ = {kHighsInf, 3, 10, 5};
  lp.sense_ = ObjSense::kMaximize;
  lp.a_matrix_.start_ = {0, 4, 8, 12};
  lp.a_matrix_.index_ = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  lp.a_matrix_.value_ = {1, 1, 1, 1, 2, 1, 2, 2, 4, 3, 2, 3};
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("solver", kPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = h.run();

  h.setOptionValue("presolve", kHighsOffString);
  run_status = h.run();
}

TEST_CASE("pdlp-restart-add-row", "[pdlp]") {
  SpecialLps special_lps;
  HighsLp lp;

  HighsModelStatus require_model_status;
  double optimal_objective;
  special_lps.ThreeDLp(lp, require_model_status, optimal_objective);

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  const HighsInfo& info = h.getInfo();
  const HighsOptions& options = h.getOptions();
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("solver", kPdlpString);
  h.setOptionValue("presolve", kHighsOffString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  HighsStatus run_status = h.run();
  if (dev_run) h.writeSolution("", 1);
  REQUIRE(std::abs(info.objective_function_value - optimal_objective) <
          double_equal_tolerance);
  const bool optimal = true;
  if (optimal) {
    REQUIRE(run_status == HighsStatus::kOk);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
  } else {
    REQUIRE(run_status == HighsStatus::kWarning);
    REQUIRE(h.getModelStatus() == HighsModelStatus::kUnknown);
  }

  if (dev_run) h.writeSolution("", 1);
  HighsSolution solution = h.getSolution();

  std::vector<HighsInt> index = {0, 1, 2};
  std::vector<double> value = {1, 1, 1};
  h.addRow(-kHighsInf, 2, 3, index.data(), value.data());

  solution.row_dual.push_back(0);

  h.setSolution(solution);
  run_status = h.run();
  if (dev_run) h.writeSolution("", 1);
}

TEST_CASE("hi-pdlp", "[pdlp]") {
  std::string model = "adlittle";//"avgas";
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs h;
  //  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(model_file) == HighsStatus::kOk);
  h.setOptionValue("solver", kHiPdlpString);
  h.setOptionValue("kkt_tolerance", kkt_tolerance);
  h.setOptionValue("pdlp_iteration_limit", 10000);
  h.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);
  HighsStatus run_status = h.run();
  //  REQUIRE(run_status == HighsStatus::kOk);
  //  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);
}
