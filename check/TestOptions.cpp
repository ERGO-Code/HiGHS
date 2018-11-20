#include "catch.hpp"
#include "HighsLp.h"

TEST_CASE("Correct status print", "[highs_data]") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str == "OK");
}