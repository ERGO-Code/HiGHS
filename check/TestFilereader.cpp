#include "catch.hpp"
#include "HighsLp.h"
#include "LoadProblem.h"
#include "HighsSetup.h"

// No commas in test case name.
TEST_CASE("read-mps-ems", "[highs_filereader]") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str == "OK");
}
