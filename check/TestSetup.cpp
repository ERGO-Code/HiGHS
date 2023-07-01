#include "HCheckConfig.h"
#include "catch.hpp"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsStatus.h"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("correct-print-status", "[highs_data]") {
  std::string str = highsStatusToString(HighsStatus::kOk);
  REQUIRE(str == "OK");
}
