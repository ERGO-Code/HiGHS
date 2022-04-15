#include <cstdio>

#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = false;

void runWriteReadCheckSolution(Highs& highs, const std::string model,
                               const HighsModelStatus require_model_status);

void runSetLpSolution(const std::string model);

TEST_CASE("check-solution", "[highs_check_solution]") {
  std::string model = "";
  std::string model_file;
  HighsStatus read_status;
  HighsStatus require_read_status;
  HighsModelStatus require_model_status;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  //  const HighsInfo& info = highs.getInfo();

  const bool test_st_test23 = false;
  if (test_st_test23) {
    model = "st-test23";
    model_file = "st-test23.lp";
    require_read_status = HighsStatus::kWarning;
  } else {
    model = "avgas";  // 25fv47";
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    require_read_status = HighsStatus::kOk;
  }

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == require_read_status);

  require_model_status = HighsModelStatus::kOptimal;
  runWriteReadCheckSolution(highs, model, require_model_status);
  SpecialLps special_lps;
  HighsLp lp;
  double optimal_objective;

  model = "distillation";
  special_lps.distillationMip(lp, require_model_status, optimal_objective);
  highs.passModel(lp);
  runWriteReadCheckSolution(highs, model, require_model_status);

  lp.clear();
  model = "primalDualInfeasible1Lp";
  special_lps.primalDualInfeasible1Lp(lp, require_model_status);
  highs.passModel(lp);
  runWriteReadCheckSolution(highs, model, require_model_status);
}

TEST_CASE("check-set-mip-solution", "[highs_check_solution]") {
  HighsStatus return_status;
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  //  const HighsInfo& info = highs.getInfo();
  if (dev_run) printf("\nSolving from scratch\n");
  highs.setOptionValue("output_flag", dev_run);

  highs.readModel(model_file);
  highs.run();
  HighsSolution solution = highs.getSolution();
  highs.clear();
  if (dev_run) printf("\nSolving from saved solution\n");
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);

  return_status = highs.setSolution(solution);
  REQUIRE(return_status == HighsStatus::kOk);
  highs.run();
}

TEST_CASE("check-set-lp-solution", "[highs_check_solution]") {
  //  runSetLpSolution("avgas");
  runSetLpSolution("adlittle");
  runSetLpSolution("shell");
  runSetLpSolution("stair");
}

void runWriteReadCheckSolution(Highs& highs, const std::string model,
                               const HighsModelStatus require_model_status) {
  HighsStatus run_status;
  HighsStatus return_status;
  std::string solution_file;
  HighsModelStatus status = HighsModelStatus::kNotset;
  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  status = highs.getModelStatus();
  REQUIRE(status == require_model_status);

  solution_file = model + ".sol";
  if (dev_run) return_status = highs.writeSolution("", kSolutionStyleRaw);
  return_status = highs.writeSolution(solution_file, kSolutionStyleRaw);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.readSolution(solution_file, kSolutionStyleRaw);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.checkSolutionFeasibility();
  REQUIRE(return_status == HighsStatus::kOk);
  std::remove(solution_file.c_str());
}

void runSetLpSolution(const std::string model) {
  HighsStatus return_status;
  Highs highs;
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  const HighsInfo& info = highs.getInfo();
  if (dev_run) printf("\nSolving from scratch\n");
  highs.setOptionValue("output_flag", dev_run);
  if (dev_run) highs.setOptionValue("log_dev_level", 1);

  highs.readModel(model_file);
  highs.run();
  HighsInt simplex_iteration_count0 = info.simplex_iteration_count;
  HighsSolution solution = highs.getSolution();
  highs.clear();
  if (dev_run) printf("\nSolving from saved solution\n");
  highs.setOptionValue("output_flag", dev_run);
  if (dev_run) highs.setOptionValue("log_dev_level", 1);
  highs.readModel(model_file);

  // solution.col_value.assign(highs.getNumCol(), 0);

  return_status = highs.setSolution(solution);
  REQUIRE(return_status == HighsStatus::kOk);

  highs.run();
  // Use a reduction in iteration count as a anity check that starting
  // from the optimal solution has worked
  HighsInt simplex_iteration_count1 = info.simplex_iteration_count;
  if (dev_run)
    printf(
        "For model %s: iteration counts are %d from logical basis and %d from "
        "optimal solution\n",
        model.c_str(), (int)simplex_iteration_count0,
        (int)simplex_iteration_count1);
  REQUIRE(simplex_iteration_count1 < simplex_iteration_count0);
}
