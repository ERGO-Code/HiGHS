
#include "../extern/filereaderlp/reader.hpp"
#include "Highs.h"
#include "catch.hpp"

TEST_CASE(
    "highs can load .lp files that contain quadratic terms without a space") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/quadratic.lp";
  Highs highs;
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);

  // todo: read the coefficients/variables off the highs model
}
