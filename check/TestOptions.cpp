#include "catch.hpp"
#include "HighsLp.h"

TEST_CASE("wrong-print", "[highs_data]") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str != "OK");
}