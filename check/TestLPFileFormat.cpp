#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "filereaderlp/reader.hpp"

const bool dev_run = false;

TEST_CASE("lp-file-format-quad-no-space", "[LpFileFormat]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/qcqp.lp";

  // HiGHS cannot handle quadratic constraints as there are in qcqp.lp

  // Test that filereaderlp does not throw an exception
  REQUIRE_NOTHROW(readinstance(filename));

  // Test that HiGHS returns an error when passes a QCQP problem
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.readModel(filename) == HighsStatus::kError);

  // todo: read the coefficients/variables off the highs model
}
