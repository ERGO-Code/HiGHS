// HiGHS Scaffold example of a unit test file.

// Catch2 has been somewhat updated since we added it to HiGHS so when moving
// the pre-release check/ unit tests make sure to copy external/catch/ subdir
// too and use that version of catch.

// This tells Catch to provide a main() -
// only do this in one cpp file per directory of unit tests defining its own
// executable target in the CMakeLists.txt file.
#define CATCH_CONFIG_MAIN

#include "Highs.h"
#include "catch.hpp"

TEST_CASE("HighsLp", "init_highslp_scaffold") {
  // This test belongs to the IsPrimeTest test case.
  HighsLp lp;
  REQUIRE(lp.numCol_ == 0);
}
