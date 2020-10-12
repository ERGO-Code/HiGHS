#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

void checkPrimalRay(Highs& highs,
		    const vector<double>& primal_ray_value) {
  int numCol = highs.getLp().numCol_;
  bool ray_error = false;
  double primal_ray_slope = 0;
  const vector<double>& colCost = highs.getLp().colCost_;
  const vector<double>& colLower = highs.getLp().colLower_;
  const vector<double>& colUpper = highs.getLp().colUpper_;
  const ObjSense sense = highs.getLp().sense_;
  double dual_feasibility_tolerance;
  highs.getHighsOptionValue("dual_feasibility_tolerance", dual_feasibility_tolerance);
  for (int iCol=0; iCol<numCol; iCol++) {
    primal_ray_slope += primal_ray_value[iCol]*colCost[iCol];
    if (primal_ray_value[iCol] > 0) {
      // Upper bound must be infinite
      if (colUpper[iCol] < HIGHS_CONST_INF) {
	ray_error = true;
	printf("Column %d has primal ray value %g and finite upper bound of %g\n",
	       iCol, primal_ray_value[iCol], colUpper[iCol]);
      }
    }
    else if (primal_ray_value[iCol] < 0) {
      // Lower bound must be infinite
      if (colLower[iCol] > -HIGHS_CONST_INF) {
	ray_error = true;
	printf("Column %d has primal ray value %g and finite lower bound of %g\n",
	       iCol, primal_ray_value[iCol], colLower[iCol]);
      }
    }
  }
  primal_ray_slope *= (int)sense;
  REQUIRE(!ray_error);
  REQUIRE(primal_ray_slope > dual_feasibility_tolerance);
}

TEST_CASE("Rays", "[highs_test_rays]") {
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

  // Test dual ray for infeasible LP
  special_lps.scipLpi3Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  bool has_dual_ray;
  vector<double> dual_ray_value;
  dual_ray_value.resize(lp.numRow_);
  REQUIRE(highs.getDualRay(has_dual_ray, &dual_ray_value[0]) ==
          HighsStatus::OK);
  REQUIRE(has_dual_ray == true);
  vector<double> exp_dualray = {0.5, -1};
  if (dev_run) {
    printf("Dual ray:\nRow    computed    expected\n");
    for (int iRow = 0; iRow < lp.numRow_; iRow++)
      printf("%3d %11.4g %11.4g\n", iRow, dual_ray_value[iRow],
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

  // Test primal ray for unbounded LP
  special_lps.scipLpi2Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);
  highs.writeSolution("", true);

  bool has_primal_ray;
  vector<double> primal_ray_value;
  primal_ray_value.resize(lp.numRow_);
  REQUIRE(highs.getPrimalRay(has_primal_ray, &primal_ray_value[0]) ==
          HighsStatus::OK);
  REQUIRE(has_primal_ray == true);
  vector<double> exp_primalray = {0.5, -1};
  if (dev_run) {
    printf("Primal ray:\nRow    computed    expected\n");
    for (int iRow = 0; iRow < lp.numRow_; iRow++)
      printf("%3d %11.4g %11.4g\n", iRow, primal_ray_value[iRow],
             exp_primalray[iRow]);
  }

  checkPrimalRay(highs, primal_ray_value);
}
