#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

// No commas in test case name.
TEST_CASE("HighsModel", "[highs_model]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  HighsStatus status;
  Highs highs;
  highs.readModel("filename");
  const HighsLp& lp = highs.getLp();
  HighsModel model;
  model.lp_ = lp;
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
}
