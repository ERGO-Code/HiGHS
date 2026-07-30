#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/HMPSIO.h"

const bool dev_run = false;

TEST_CASE("run-data-md", "[highs_run_data]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  // Use this name so that it can be copied to docs and provides code
  // coverage
  const std::string run_data_file = "HighsRunData.md";
  REQUIRE(h.writeRunData(run_data_file) == HighsStatus::kOk);
  if (!dev_run) std::remove(run_data_file.c_str());
}

TEST_CASE("highs-run-data", "[highs_run_data]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string highs_run_data_file = test_name + ".run_data";

  Highs h;
  if (!dev_run) h.setOptionValue("output_flag", false);
  const HighsRunData& highs_run_data = h.getRunData();

  auto testRunData = [&](const std::string& filename) {
    HighsStatus return_status = h.readModel(filename);
    REQUIRE(return_status == HighsStatus::kOk);

    // Cannot write run_data since not valid before run()
    return_status = h.writeRunData("");
    REQUIRE(return_status == HighsStatus::kWarning);

    HighsRunDataType highs_run_data_type;
    return_status = h.getRunDataType("presolved_num_col", highs_run_data_type);
    REQUIRE(return_status == HighsStatus::kError);
    return_status =
        h.getRunDataType("presolved_model_num_col", highs_run_data_type);
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(highs_run_data_type == HighsRunDataType::kInt);

    return_status = h.getRunDataType("presolving_time", highs_run_data_type);
    REQUIRE(return_status == HighsStatus::kError);
    return_status = h.getRunDataType("presolve_time", highs_run_data_type);
    REQUIRE(return_status == HighsStatus::kOk);
    REQUIRE(highs_run_data_type == HighsRunDataType::kDouble);

    // Run data not valid before run()
    HighsInt presolved_model_num_col;
    return_status =
        h.getRunDataValue("presolved_model_num_col", presolved_model_num_col);
    REQUIRE(return_status == HighsStatus::kWarning);

    return_status = h.run();
    REQUIRE(return_status == HighsStatus::kOk);

    if (dev_run) {
      return_status = h.writeRunData("");
      REQUIRE(return_status == HighsStatus::kOk);
    }

    return_status = h.writeRunData(highs_run_data_file);
    REQUIRE(return_status == HighsStatus::kOk);

    // Wrong name for objective
    return_status =
        h.getRunDataValue("presolved_num_col", presolved_model_num_col);
    REQUIRE(return_status == HighsStatus::kError);

    // Right name for objective
    return_status =
        h.getRunDataValue("presolved_model_num_col", presolved_model_num_col);
    REQUIRE(return_status == HighsStatus::kOk);

    if (dev_run)
      printf("From getRunDataValue: presolved_model_num_col = %d\n",
             int(presolved_model_num_col));

    double presolve_time;
    // Wrong name for simplex iteration count
    return_status = h.getRunDataValue("presolving_time", presolve_time);
    REQUIRE(return_status == HighsStatus::kError);

    // Right name for presolve time
    return_status = h.getRunDataValue("presolve_time", presolve_time);
    REQUIRE(return_status == HighsStatus::kOk);

    const HighsModelStatus model_status = h.getModelStatus();
    if (dev_run) {
      printf("From getModelStatus: model_status = %s\n",
             h.modelStatusToString(model_status).c_str());
      printf("From getRunData: presolved_model_num_col = %d\n",
             int(highs_run_data.presolved_model_num_col));
      printf("From getRunData: presolved_model_num_row = %d\n",
             int(highs_run_data.presolved_model_num_row));
      printf("From getRunData: presolved_model_num_nz  = %d\n",
             int(highs_run_data.presolved_model_num_nz));
      if (!h.getLp().isMip())
        printf(
            "From getRunData: num_simplex_iterations_after_postsolve  = %d\n",
            int(highs_run_data.num_simplex_iterations_after_postsolve));
      printf("From getRunData:  presolve_time = %g\n",
             highs_run_data.presolve_time);
      printf("From getRunData:     solve_time = %g\n",
             highs_run_data.solve_time);
      printf("From getRunData: postsolve_time = %g\n",
             highs_run_data.postsolve_time);
    }
    REQUIRE(highs_run_data.presolved_model_num_col >= 0);
    REQUIRE(highs_run_data.presolved_model_num_row >= 0);
    REQUIRE(highs_run_data.presolved_model_num_nz >= 0);
    if (!h.getLp().isMip())
      REQUIRE(highs_run_data.num_simplex_iterations_after_postsolve == 0);
    REQUIRE(highs_run_data.presolve_time >= 0);
    REQUIRE(highs_run_data.solve_time >= 0);
    REQUIRE(highs_run_data.postsolve_time >= 0);
  };

  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  testRunData(filename);

  filename = std::string(HIGHS_DIR) + "/check/instances/egout-ac.mps";
  testRunData(filename);

  // Doesn't work for MIPs yet, but wait until profiling is merged in
  // to avoid conflicts
  //
  //  filename = std::string(HIGHS_DIR) + "/check/instances/flugpl.mps";
  //  testRunData(filename);

  if (!dev_run) std::remove(highs_run_data_file.c_str());

  h.resetGlobalScheduler(true);
}
