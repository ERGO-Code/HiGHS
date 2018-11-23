#include "catch.hpp"
#include "HighsLp.h"

// No commas in test case name.
TEST_CASE("correct-print", "[highs_data]") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str == "OK.");
}
