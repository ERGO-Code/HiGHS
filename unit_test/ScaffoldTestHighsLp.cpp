// HiGHS Scaffold example of a unit test file.

// Catch2 has been somewhat updated so this uses the new version 
// todo: when moving the pre-release check/ unit tests make sure to copy external/catch/ subdir

#include "Highs.h"
#include "catch.hpp"

TEST_CASE("HighsLp", "init_highslp_scaffold") {
  // This test belongs to the IsPrimeTest test case.
  HighsLp lp;
  REQUIRE(lp.numCol_ == 0);
}
