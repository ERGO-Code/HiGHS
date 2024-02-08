#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

TEST_CASE("tsp-p01", "[highs_test_tsp_solver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/p01.mps";
  const double optimal_obective_value = 263;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  highs.run();
  REQUIRE(highs.getObjectiveValue() == optimal_obective_value);
}
