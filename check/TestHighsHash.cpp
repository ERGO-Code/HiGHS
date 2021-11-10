
#include "catch.hpp"
#include "util/HighsHash.h"

TEST_CASE("Highs_log2i", "[util]") {
  // test 32 bit and 64 bit values whoes log2 value should be floored
  uint32_t x = 12345;
  uint64_t y = 1234567891011;
  REQUIRE(HighsHashHelpers::log2i(x) == 13);
  REQUIRE(HighsHashHelpers::log2i(y) == 40);

  // test 2^13 and 2^40 which should yield the same result
  x = 8192;
  y = 1099511627776;
  REQUIRE(HighsHashHelpers::log2i(x) == 13);
  REQUIRE(HighsHashHelpers::log2i(y) == 40);
}
