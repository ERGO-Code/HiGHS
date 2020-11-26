#include "HighsLp.h"
#include "HighsStatus.h"
#include "catch.hpp"

const bool dev_run = false;
const bool use_ekk = true;

// No commas in test case name.
TEST_CASE("correct-print-status", "[highs_data]") {
  std::string str = HighsStatusToString(HighsStatus::OK);
  REQUIRE(str == "OK");
}
