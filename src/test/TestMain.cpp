#define CATCH_CONFIG_MAIN
#include "HighsLp.h"
#include "catch.hpp"

TEST_CASE("status") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str == "OK");
}