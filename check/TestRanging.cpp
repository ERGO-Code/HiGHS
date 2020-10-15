#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

TEST_CASE("Ranging", "[highs_test_ranging]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }

  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;

  model_file = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  require_model_status = HighsModelStatus::OPTIMAL;
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);
  lp = highs.getLp();

  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);
  //  double optimal_objective = highs.getObjectiveValue();

  HighsRanging ranging;
  REQUIRE(highs.getRanging(ranging) == HighsStatus::OK);
  
}
