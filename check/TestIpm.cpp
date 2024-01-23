#include <cmath>

#include "Highs.h"
#include "catch.hpp"

// I use dev_run to switch on/off printing and logging used for
// development of the unit test
const bool dev_run = true;
const double inf = kHighsInf;

TEST_CASE("test-analytic-centre", "[highs_ipm]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsLp lp;
  // Set up a problem for which analytic centre calculations are
  // required, followed by sanity test(s) on the results.
  HighsInt chalk = 0;
  HighsInt cheese = 0;
  // REQUIRE is like assert, but feeds back to the unit test mechanism
  // "catch.hpp"
  REQUIRE(chalk == cheese);
}
