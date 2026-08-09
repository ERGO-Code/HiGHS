#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/HMPSIO.h"

const bool dev_run = false;  // true;//

const std::vector<std::string> solvers{
    //    kHighsChooseString
    kSimplexString,
    //  kIpxString,
    kHipoString
    //    kQpAsmString
    //    kHiPdlpString
};

void testRunData(Highs& h, const bool irreducible, const bool reduces_to_empty,
                 const std::string& run_data_file);

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
  // Doesn't work for MIPs yet, but wait until profiling is merged in
  // to avoid conflicts
  const std::vector<std::string> models{
      "adlittle", "egout-ac"
      //    "flugpl"
  };
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string run_data_file = test_name + ".run_data";

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  const bool irreducible = false;
  for (auto& model : models) {
    std::string filename =
        std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    REQUIRE(h.readModel(filename) == HighsStatus::kOk);
    HighsStatus return_status = h.readModel(filename);
    REQUIRE(return_status == HighsStatus::kOk);
    const bool reduces_to_empty = model == "egout-ac" ? true : false;

    for (auto& solver : solvers)
      testRunData(h, irreducible, reduces_to_empty, run_data_file);
  }
  if (!dev_run) std::remove(run_data_file.c_str());

  h.resetGlobalScheduler(true);
}

TEST_CASE("highs-run-data-presolve", "[highs_run_data]") {
  const std::vector<std::string> models{"adlittle", "flugpl"};
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string run_data_file = test_name + ".run_data";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  const HighsRunData& run_data = h.getRunData();
  const HighsLp& lp = h.getLp();
  for (auto& model : models) {
    std::string filename =
        std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    REQUIRE(h.readModel(filename) == HighsStatus::kOk);
    const bool irreducible = true;
    const bool reduces_to_empty = false;
    for (auto& solver : solvers) {
      h.setOptionValue("solver", solver);
      //      if (dev_run)
      printf("\n!>>>>%s-%s<<<<\n", model.c_str(), solver.c_str());

      REQUIRE(h.presolve() == HighsStatus::kOk);
      HighsLp presolved_lp = h.getPresolvedLp();

      h.passModel(presolved_lp);
      h.setOptionValue("solve_relaxation", true);
      h.setOptionValue(kPresolveString, kHighsOffString);
      testRunData(h, irreducible, reduces_to_empty, run_data_file);
    }
  }

  h.resetGlobalScheduler(true);
}

void testRunData(Highs& h, const bool irreducible, const bool reduces_to_empty,
                 const std::string& run_data_file) {
  assert(irreducible != reduces_to_empty);
  const HighsRunData& run_data = h.getRunData();
  const HighsLp& lp = h.getLp();

  std::string presolve;
  h.getOptionValue(kPresolveString, presolve);
  const bool run_presolve = presolve != kHighsOffString;

  // Cannot write run_data since not valid before run()
  HighsStatus return_status = h.writeRunData("");
  REQUIRE(return_status == HighsStatus::kWarning);

  HighsRunDataType run_data_type;
  return_status = h.getRunDataType("presolved_num_col", run_data_type);
  REQUIRE(return_status == HighsStatus::kError);
  return_status = h.getRunDataType("presolved_model_num_col", run_data_type);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(run_data_type == HighsRunDataType::kInt);

  return_status = h.getRunDataType("presolving_time", run_data_type);
  REQUIRE(return_status == HighsStatus::kError);
  return_status = h.getRunDataType("presolve_time", run_data_type);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(run_data_type == HighsRunDataType::kDouble);

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

  return_status = h.writeRunData(run_data_file);
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
  // Wrong name for presolve_time
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
           int(run_data.presolved_model_num_col));
    printf("From getRunData: presolved_model_num_row = %d\n",
           int(run_data.presolved_model_num_row));
    printf("From getRunData: presolved_model_num_nz  = %d\n",
           int(run_data.presolved_model_num_nz));
    if (!h.getLp().isMip())
      printf("From getRunData: num_simplex_iterations_after_postsolve  = %d\n",
             int(run_data.num_simplex_iterations_after_postsolve));
    printf("From getRunData:  presolve_time = %g\n", run_data.presolve_time);
    printf("From getRunData:     solve_time = %g\n", run_data.solve_time);
    printf("From getRunData: postsolve_time = %g\n", run_data.postsolve_time);
  }
  if (run_presolve) {
    REQUIRE(run_data.presolve_time >= 0);
    REQUIRE(run_data.presolve_time < kHighsInf);
    REQUIRE(run_data.presolved_model_num_col >= 0);
    REQUIRE(run_data.presolved_model_num_row >= 0);
    REQUIRE(run_data.presolved_model_num_nz >= 0);
    if (!irreducible) {
      REQUIRE(run_data.presolved_model_num_col < lp.num_col_);
      REQUIRE(run_data.presolved_model_num_row < lp.num_row_);
      REQUIRE(run_data.presolved_model_num_nz < lp.a_matrix_.numNz());
    }
    if (reduces_to_empty) {
      REQUIRE(run_data.presolved_model_num_col == 0);
      REQUIRE(run_data.presolved_model_num_row == 0);
      REQUIRE(run_data.presolved_model_num_nz == 0);
    }
    REQUIRE(run_data.postsolve_time >= 0);
    REQUIRE(run_data.postsolve_time < kHighsInf);
    if (!h.getLp().isMip()) {
      REQUIRE(run_data.num_simplex_iterations_after_postsolve == 0);
    }
  } else {
    REQUIRE(run_data.presolve_time == kHighsIllegalDoubleMeasure);
    REQUIRE(run_data.presolved_model_num_col == kHighsIllegalIntMeasure);
    REQUIRE(run_data.presolved_model_num_row == kHighsIllegalIntMeasure);
    REQUIRE(run_data.presolved_model_num_nz == kHighsIllegalIntMeasure);
    REQUIRE(run_data.postsolve_time == kHighsIllegalDoubleMeasure);
    REQUIRE(run_data.num_simplex_iterations_after_postsolve ==
            kHighsIllegalIntMeasure);
  }
  REQUIRE(run_data.solve_time >= 0);
  REQUIRE(run_data.solve_time < kHighsInf);
  h.clearSolver();
}
