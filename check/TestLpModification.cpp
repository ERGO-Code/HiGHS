#include "catch.hpp"
#include "Highs.h"
#include "Avgas.h"

// No commas in test case name.
TEST_CASE("LP-modification", "[highs_data]") {
  // Create an empty LP
  HighsLp lp;
  HighsOptions options;
  HighsSetMessagelevel(ML_ALWAYS);

  Avgas avgas;
  int num_row;
  vector<double> rowLower;
  vector<double> rowUpper;
  avgas.rows(num_row, rowLower, rowUpper);

  int num_col = 0;
  int num_nz = 0;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  for (int col = 0; col < 8; col++) {
    avgas.col(col, num_col, num_nz, colCost, colLower, colUpper, Astart, Aindex, Avalue);
  }
  Highs highs(options);
  HighsStatus return_status;
  return_status = highs.initializeLp(lp);
  //  printf("initializeLp: return_status = %s\n", HighsStatusToString(return_status).c_str());
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  HighsStatusReport("highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::LpEmpty);

  //  highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], &Astart[0], num_nz, &Aindex[0], &Avalue[0], false);
}

