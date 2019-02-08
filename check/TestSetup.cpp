#include "catch.hpp"
#include "HighsLp.h"
#include "Filereader.h"

// No commas in test case name.
TEST_CASE("correct-print-input-status", "[highs_data]") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str == "OK");
}
