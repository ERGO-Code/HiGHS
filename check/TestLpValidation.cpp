#include "catch.hpp"
#include "HighsLp.h"
#include "HighsOptions.h"
#include "HighsStatus.h"
#include "HighsLpUtils.h"

// No commas in test case name.
TEST_CASE("LP-validation", "[highs_data]") {
  HighsLp lp;
  HighsOptions options;
  HighsStatus return_status;
  return_status = assessLp(lp, options);
  REQUIRE(return_status == HighsStatus::OK);
}
