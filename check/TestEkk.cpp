#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

TEST_CASE("Ekk", "[highs_test_ekk]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;

  const bool from_file = true;
  if (from_file) {
    std::string model_file = std::string(HIGHS_DIR) +
                             // "/check/instances/stair.mps";
                             // "/check/instances/adlittle.mps";
                             "/check/instances/cycle.mps";
    REQUIRE(highs.readModel(model_file) == HighsStatus::OK);
  } else {
    SpecialLps special_lps;
    //    special_lps.blendingLp(lp, require_model_status, optimal_objective);
    special_lps.distillationLp(lp, require_model_status, optimal_objective);
    highs.passModel(lp);
  }
  REQUIRE(highs.setHighsOptionValue("simplex_strategy", SIMPLEX_STRATEGY_EKK) ==
          HighsStatus::OK);
  highs.setHighsOptionValue("message_level", 6);
  REQUIRE(highs.run() == HighsStatus::OK);
}
