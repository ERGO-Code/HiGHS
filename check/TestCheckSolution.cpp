#include <cstdio>

#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"

const bool dev_run = true;

void runWriteReadCheckSolution(Highs& highs, const std::string model,
                               const HighsModelStatus require_model_status);

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

TEST_CASE("check-set-solution", "[highs_check_solution]") {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  if (dev_run) printf("\nSolving from scratch\n");
  highs.setOptionValue("output_flag", dev_run);

  highs.readModel(model_file);
  highs.run();
  HighsSolution solution = highs.getSolution();
  highs.clear();
  if (dev_run) printf("\nSolving from saved solution\n");
  highs.setOptionValue("output_flag", dev_run);
  highs.readModel(model_file);

  highs.setSolution(solution);
  highs.run();
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
