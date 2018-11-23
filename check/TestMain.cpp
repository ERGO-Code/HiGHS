#define CATCH_CONFIG_MAIN
#include "HighsLp.h"
#include "catch.hpp"

// No commas in test case name.
TEST_CASE("correct-print", "[highs_data]") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str == "OK");
}