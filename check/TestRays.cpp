#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

TEST_CASE("Dual-ray", "[highs_test_rays]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  SpecialLps special_lps;
  special_lps.issue285Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  bool has_dual_ray;
  vector<double> dual_ray_values;
  dual_ray_values.resize(lp.numCol_+lp.numRow_);
  REQUIRE(highs.getDualRay(has_dual_ray, &dual_ray_values[0]) ==
          HighsStatus::OK);
  //  REQUIRE(has_dual_ray==true);

  // Check that there is no ray for this LP
  special_lps.issue272Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::OK);
  REQUIRE(has_dual_ray == false);
}
