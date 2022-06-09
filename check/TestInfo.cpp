#include <cstdio>

#include "FilereaderEms.h"
#include "HMPSIO.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("highs-info", "[highs_info]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  //  filename = std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& highs_info = highs.getInfo();

  HighsStatus return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  // Cannot write info since not valid before run()
  return_status = highs.writeInfo("");
  REQUIRE(return_status == HighsStatus::kWarning);

  // Only use IPX if and using 32-bit arithmetic
  bool use_ipx = false;
#ifndef HIGHSINT64
  use_ipx = true;
#endif
  if (use_ipx) {
    return_status = highs.setOptionValue("solver", "ipm");
    REQUIRE(return_status == HighsStatus::kOk);
  }

  // Info not valid before run()
  double objective_function_value;
  return_status =
      highs.getInfoValue("objective_function_value", objective_function_value);
  REQUIRE(return_status == HighsStatus::kWarning);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    return_status = highs.writeInfo("");
    REQUIRE(return_status == HighsStatus::kOk);
  }

  std::string highs_info_file = "Highs.info";
  return_status = highs.writeInfo(highs_info_file);
  REQUIRE(return_status == HighsStatus::kOk);

  // Wrong name for objective
  return_status =
      highs.getInfoValue("objective_value", objective_function_value);
  REQUIRE(return_status == HighsStatus::kError);

  // Right name for objective
  return_status =
      highs.getInfoValue("objective_function_value", objective_function_value);
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run)
    printf("From getInfoValue: objective_function_value = %g\n",
           objective_function_value);

  HighsInt simplex_iteration_count;
  // Wrong name for simplex iteration count
  return_status =
      highs.getInfoValue("iteration_count", simplex_iteration_count);
  REQUIRE(return_status == HighsStatus::kError);

  // Right name for simplex iteration count
  return_status =
      highs.getInfoValue("simplex_iteration_count", simplex_iteration_count);
  REQUIRE(return_status == HighsStatus::kOk);

  const HighsModelStatus model_status = highs.getModelStatus();
  if (dev_run) {
    printf("From getModelStatus: model_status = %s\n",
           highs.modelStatusToString(model_status).c_str());
    printf("From getInfo: objective_function_value = %g\n",
           highs_info.objective_function_value);
    if (use_ipx) {
      printf("From getInfo: ipm_iteration_count = %" HIGHSINT_FORMAT "\n",
             highs_info.ipm_iteration_count);
    } else {
      printf("From getInfo: simplex_iteration_count = %" HIGHSINT_FORMAT "\n",
             highs_info.simplex_iteration_count);
    }
  }
  std::remove(highs_info_file.c_str());
}
