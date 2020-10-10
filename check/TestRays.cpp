#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

TEST_CASE("Dual-ray", "[highs_test_rays]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  HighsLp lp;
  HighsModelStatus require_model_status;
  SpecialLps special_lps;
  special_lps.issue285Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
}
