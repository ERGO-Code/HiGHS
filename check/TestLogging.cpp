#include <cstdio>
#include <fstream>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("logging", "[highs_logging]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string log_file = test_name + ".log";
  std::string model;
  std::string model_file;
  HighsStatus return_status;

  model = "avgas";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

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

  highs.resetThreadScheduler(true);
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

  highs.resetThreadScheduler(true);
}

// Number of lines of a file that contain the log header
static HighsInt numLogHeaders(const std::string& file_name) {
  std::ifstream file(file_name);
  HighsInt num_log_headers = 0;
  for (std::string line; std::getline(file, line);)
    if (line.find("Running HiGHS") != std::string::npos) num_log_headers++;
  return num_log_headers;
}

TEST_CASE("log-header", "[highs_logging]") {
  const std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  HighsCallbackFunctionType callback = [](int, const std::string&,
                                          const HighsCallbackOutput*,
                                          HighsCallbackInput*, void*) {};
  // The C API sets callbacks through this overload of setCallback
  HighsCCallbackType c_callback = [](int, const char*,
                                     const HighsCallbackDataOut*,
                                     HighsCallbackDataIn*, void*) {};
  for (HighsInt instance = 0; instance < 2; instance++) {
    // Each instance writes the log header to its own log file
    const std::string log_file =
        "log-header" + std::to_string(instance) + ".log";
    Highs highs;
    highs.setOptionValue("log_to_console", dev_run);
    highs.setOptionValue(kLogFileString, log_file);
    // Setting, starting and stopping callbacks does not write the
    // header, so that a log callback can receive it
    REQUIRE(highs.setCallback(c_callback, nullptr) == HighsStatus::kOk);
    REQUIRE(highs.setCallback(callback) == HighsStatus::kOk);
    const int logging = kCallbackLogging;
    REQUIRE(highs.startCallback(logging) == HighsStatus::kOk);
    REQUIRE(highs.stopCallback(logging) == HighsStatus::kOk);
    REQUIRE(highs.startCallback(kCallbackLogging) == HighsStatus::kOk);
    REQUIRE(highs.stopCallback(kCallbackLogging) == HighsStatus::kOk);
    REQUIRE(numLogHeaders(log_file) == 0);
    // The header is written once, when the instance is first used
    REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    // Close the log file
    highs.setOptionValue(kLogFileString, "");
    REQUIRE(numLogHeaders(log_file) == 1);
    std::remove(log_file.c_str());
  }

  Highs::resetThreadScheduler(true);
}
