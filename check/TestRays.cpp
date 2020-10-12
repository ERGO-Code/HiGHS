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
  //  special_lps.issue285Lp(lp, require_model_status);
  REQUIRE(highs.setHighsOptionValue("presolve", "off") == HighsStatus::OK);

  special_lps.scipLpi3Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  bool has_dual_ray;
  vector<double> dual_ray_values;
  dual_ray_values.resize(lp.numRow_);
  REQUIRE(highs.getDualRay(has_dual_ray, &dual_ray_values[0]) ==
          HighsStatus::OK);
  REQUIRE(has_dual_ray == true);
  vector<double> exp_dualray = {0.5, -1};
  if (dev_run) {
    printf("Dual ray:\nRow    computed    expected\n");
    for (int iRow = 0; iRow < lp.numRow_; iRow++)
      printf("%3d %11.4g %11.4g\n", iRow, dual_ray_values[iRow],
             exp_dualray[iRow]);
  }

  // Check that there is no ray for this LP
  special_lps.issue272Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::OK);
  REQUIRE(has_dual_ray == false);

  special_lps.scipLpi2Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);
  highs.writeSolution("", true);

  bool has_primal_ray;
  vector<double> primal_ray_values;
  primal_ray_values.resize(lp.numRow_);
  REQUIRE(highs.getPrimalRay(has_primal_ray, &primal_ray_values[0]) ==
          HighsStatus::OK);
  REQUIRE(has_primal_ray == true);
  vector<double> exp_primalray = {0.5, -1};
  if (dev_run) {
    printf("Primal ray:\nRow    computed    expected\n");
    for (int iRow = 0; iRow < lp.numRow_; iRow++)
      printf("%3d %11.4g %11.4g\n", iRow, primal_ray_values[iRow],
             exp_primalray[iRow]);
  }
}
