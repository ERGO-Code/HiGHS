#include <cstdio>

#include "FilereaderLp.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("qpsolver", "[qpsolver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/qptestnw.lp";

  Highs highs;

  HighsStatus return_status = highs.readModel(filename);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  double objval = highs.getObjectiveValue();
  REQUIRE(fabs(objval + 6.45) < 10E-5);
  REQUIRE(return_status == HighsStatus::kOk);

  HighsSolution sol = highs.getSolution();
  REQUIRE(fabs(sol.col_value[0] - 1.4) < 10E-5);
  REQUIRE(fabs(sol.col_value[1] - 1.7) < 10E-5);
}
