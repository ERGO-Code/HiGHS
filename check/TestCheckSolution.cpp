#include <cstdio>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("check-solution", "[highs_check_solution]") {
  std::string model = "";
  std::string model_file;
  std::string solution_file;
  HighsStatus run_status;
  HighsStatus read_status;
  HighsStatus require_read_status;
  HighsStatus return_status;
  Highs highs;
  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
  }
  //  const HighsInfo& info = highs.getInfo();

  const bool test_st_test23 = false;
  if (test_st_test23) {
    model = "st-test23";
    model_file = "st-test23.lp";
    require_read_status = HighsStatus::kWarning;
  } else {
    model = "25fv47";
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    require_read_status = HighsStatus::kOk;
  }

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == require_read_status);

  const bool solve_ml = true;
  HighsModelStatus status = HighsModelStatus::kNotset;
  if (solve_ml) {
    run_status = highs.run();
    REQUIRE(run_status == HighsStatus::kOk);

    status = highs.getModelStatus();
    REQUIRE(status == HighsModelStatus::kOptimal);

    solution_file = model + ".sol";
    return_status = highs.writeSolution(solution_file, kWriteSolutionStyleRaw);
    REQUIRE(return_status == HighsStatus::kOk);
  }

  return_status = highs.readSolution(solution_file, kWriteSolutionStyleRaw);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.checkSolutionFeasibility();
  REQUIRE(return_status == HighsStatus::kOk);
  std::remove(solution_file.c_str());
}
