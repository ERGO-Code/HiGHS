#include "Highs.h"
#include "HCheckConfig.h"
#include "catch.hpp"
#include "filereaderlp/reader.hpp"

const bool dev_run = false;

TEST_CASE("lp-file-format-quad-no-space", "[LpFileFormat]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/quadratic.lp";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel(filename) == HighsStatus::kOk);

  // todo: read the coefficients/variables off the highs model
}
