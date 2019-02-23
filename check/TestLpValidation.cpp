#include "catch.hpp"
#include "HighsLp.h"
#include "HighsOptions.h"
#include "HighsStatus.h"
#include "HighsTimer.h"
#include "HighsLpUtils.h"
#include "HighsSimplexInterface.h"

// No commas in test case name.
TEST_CASE("LP-validation", "[highs_data]") {
  HighsLp lp;
  HighsOptions options;
  HighsTimer timer;
  HighsStatus return_status;
  return_status = checkLp(lp);
  REQUIRE(return_status == HighsStatus::OK);
  bool normalise = true;
  return_status = assessLp(lp, options, normalise);
  REQUIRE(return_status == HighsStatus::OK);
  reportLp(lp);
  int XnumNewCol = 1;
  int XnumNewNZ = 0;
  vector<double> XcolCost; XcolCost.resize(XnumNewCol);
  vector<double> XcolLower; XcolLower.resize(XnumNewCol);
  vector<double> XcolUpper; XcolUpper.resize(XnumNewCol);
  vector<int> XAstart; XAstart.resize(XnumNewCol);
  vector<int> XAindex; 
  vector<double> XAvalue;
  bool force = false;
  HighsModelObject hmo(lp, options, timer);
  HighsSimplexInterface hsi(hmo);
  return_status = hsi.util_add_cols(XnumNewCol, (double*)&XcolCost[0], (double*)&XcolLower[0], (double*)&XcolUpper[0],
				    XnumNewNZ, (int*)&XAstart[0], (int*)&XAindex[0], (double*)&XAindex[0],
				    force);
  REQUIRE(return_status == HighsStatus::OK);
  reportLp(lp);
}
