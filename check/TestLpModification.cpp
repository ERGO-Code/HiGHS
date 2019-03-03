#include "catch.hpp"
#include "Highs.h"
#include "lp_data/HighsLpUtils.h"
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
  const HighsLp &reference_lp = highs.getLp();
  HighsSetMessagelevel(ML_ALWAYS);
  reportLp(reference_lp, 2);
  HighsSetMessagelevel(ML_NONE);

  // Test the extraction of intervals from interval
  bool interval = false;
  int from_col = 0;
  int to_col = 0;
  bool set = false;
  int num_set_entries = 4;
  int col_set[] = {1, 4, 5, 9};
  int current_set_entry;
  bool mask = false;
  int col_mask[] = {0,1,0,0,1,1,0,0,1,0};

  int col_dim = num_row;
  int delete_from_col;
  int delete_to_col;
  int keep_from_col;
  int keep_to_col;

  interval = true;
  from_col = 3;
  to_col = 7;

  keep_to_col = 0;
  printf("With index interval [%d, %d) in [%d, %d)\n", from_col, to_col, 0, col_dim);
  update_delete_keep_ix(col_dim,
			interval, from_col, to_col,
			set, 0, NULL,
			mask, NULL,
			delete_from_col, delete_to_col,
			keep_from_col, keep_to_col,
			current_set_entry);
  printf("Keep [%d, %d); delete [%d, %d); keep [%d, %d)\n", 0, delete_from_col, delete_from_col, delete_to_col, keep_from_col, keep_to_col);
  interval = false;

  set = true;
  current_set_entry = 0;
  keep_to_col = 0;
  printf("With index set\n");
  for (int set = 0; set < num_set_entries; set++) {
    printf("%2d: %d\n", set, col_set[set]);
  }
  update_delete_keep_ix(col_dim,
			interval, 0, 0,
			set, num_set_entries, col_set,
			mask, NULL,
			delete_from_col, delete_to_col,
			keep_from_col, keep_to_col,
			current_set_entry);
  printf("Keep   [%2d, %2d)\n", 0, col_set[0]);
  for (int set = 0; set < num_set_entries; set++) {
    update_delete_keep_ix(col_dim,
			  interval, 0, 0,
			  set, num_set_entries, col_set,
			  mask, NULL,
			  delete_from_col, delete_to_col,
			  keep_from_col, keep_to_col,
			  current_set_entry);
    printf("Delete [%2d, %2d); keep [%2d, %2d)\n", delete_from_col, delete_to_col, keep_from_col, keep_to_col);
    if (keep_to_col == col_dim) break;
  }

  //  highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], &Astart[0], num_nz, &Aindex[0], &Avalue[0]);
}

