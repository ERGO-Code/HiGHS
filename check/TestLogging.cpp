#include <cstdio>

//#include "HMPSIO.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

TEST_CASE("logging", "[highs_logging]") {
  std::string model;
  std::string model_file;
  HighsStatus return_status;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  
  // Read mps
  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kOk);

  highs.setOptionValue("log_file", "temp.log");

}
