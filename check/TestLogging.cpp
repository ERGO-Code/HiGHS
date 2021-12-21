#include <cstdio>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

TEST_CASE("logging", "[highs_logging]") {
  std::string model;
  std::string model_file;
  std::string log_file;
  HighsStatus return_status;

  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  log_file = "temp.log";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  // By default, initial output is just to to console
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kOk);

  // Setting log_file to a non-empty string opens the file
  highs.setOptionValue("log_file", log_file);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (!dev_run) std::remove(log_file.c_str());
}
