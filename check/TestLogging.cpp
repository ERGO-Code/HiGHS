#include <cstdio>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("logging", "[highs_logging]") {
  std::string model;
  std::string model_file;
  std::string log_file;
  HighsStatus return_status;

  model = "avgas";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  log_file = "temp.log";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  // By default, initial output is just to to console
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kOk);

  // Setting log_file to a non-empty string opens the file
  highs.setOptionValue(kLogFileString, log_file);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  // Setting log_file to an empty string closes the file and prevents
  // further logging to file
  highs.setOptionValue(kLogFileString, "");

  // Setting log_to_console false suppresses logging to console
  highs.setOptionValue("log_to_console", false);
  // Writing out the the info should puroduce no output to console or file
  highs.run();

  // Setting log_to_console true restores logging to console
  highs.setOptionValue("log_to_console", true);
  if (dev_run) printf("After setting log_to_console = true\n");
  // Writing out the the info should puroduce output to console
  highs.run();

  if (!dev_run) std::remove(log_file.c_str());
}

TEST_CASE("no-logging", "[highs_logging]") {
  std::string model;
  std::string model_file;
  std::string log_file;
  HighsStatus return_status;

  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
}
