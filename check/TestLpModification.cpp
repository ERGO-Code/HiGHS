#include "catch.hpp"
#include "Highs.h"
#include "lp_data/HighsLpUtils.h"
#include "Avgas.h"

void test_delete_keep(const int row_dim,
		      const bool interval, const int from_row, const int to_row,
		      const bool set, const int num_set_entries, const int* row_set,
		      const bool mask, const int* row_mask) {
  int delete_from_row;
  int delete_to_row;
  int keep_from_row;
  int keep_to_row;
  int current_set_entry;
  if (interval) {
    printf("With index interval [%d, %d) in [%d, %d)\n", from_row, to_row, 0, row_dim);
  } else if (set) {
    printf("With index set\n");
    for (int set = 0; set < num_set_entries; set++) printf(" %2d", set); printf("\n");
    for (int set = 0; set < num_set_entries; set++) printf(" %2d", row_set[set]); printf("\n");
  } else {
    printf("With index mask\n");
    for (int row = 0; row < row_dim; row++) printf(" %2d", row); printf("\n");
    for (int row = 0; row < row_dim; row++) printf(" %2d", row_mask[row]); printf("\n");
  }
  
  keep_from_row = 0;
  if (interval) {
    keep_to_row = from_row;
  } else if (set) {
    current_set_entry = 0;
    keep_to_row = row_set[0];
  } else {
    keep_to_row = row_dim;
    for (int row = 0; row < row_dim; row++) {
      if (row_mask[row]) {
	keep_to_row = row;
	break;	
      }
    }
  }
  printf("Keep   [%2d, %2d)\n", 0, keep_to_row);
  if (keep_to_row >= row_dim) return;
  for (int k = 0; k < row_dim; k++) {
    update_delete_keep_ix(row_dim,
			  interval, from_row, to_row,
			  set, num_set_entries, row_set,
			  mask, row_mask,
			  delete_from_row, delete_to_row,
			  keep_from_row, keep_to_row,
			  current_set_entry);
    printf("Delete [%2d, %2d); keep [%2d, %2d)\n", delete_from_row, delete_to_row, keep_from_row, keep_to_row);
    if (delete_to_row == row_dim || keep_to_row == row_dim) break;
  }
}

bool test_all_delete_keep(int num_row) {
  // Test the extraction of intervals from interval, set and mask
  bool interval = false;
  bool set = false;
  bool mask = false;
  int row_dim = num_row;
  
  int from_row = 3;
  int to_row = 7;
  int num_set_entries = 4;
  int row_set[] = {1, 4, 5, 8};
  int current_set_entry;
  int row_mask[] = {0,1,0,0,1,1,0,0,1,0};
  for (int pass = 0; pass < 3; pass++) {
    printf("\nTesting delete-keep: pass %d\n", pass);
    if (pass == 1) {
      from_row = 0;
      row_set[0] = 0;
      row_mask[0] = 1;
    } else if (pass == 2) {
      from_row = 3;
      to_row = 10;
      row_set[0] = 1;
      row_set[3] = 9;
      row_mask[0] = 0;
      row_mask[9] = 1;
    }
    
    interval = true;
    test_delete_keep(row_dim,
		     interval, from_row, to_row,
		     set, num_set_entries, row_set,
		     mask, row_mask);
    interval = false;
    set = true;
    test_delete_keep(row_dim,
		     interval, from_row, to_row,
		     set, num_set_entries, row_set,
		     mask, row_mask);
    set = false;
    mask = true;
    test_delete_keep(row_dim,
		     interval, from_row, to_row,
		     set, num_set_entries, row_set,
		     mask, row_mask);
  }
  return true;
}


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

  //  test_all_delete_keep(num_row);

  // Add column vectors and matrix to model with no rows returns an error
  bool return_bool;
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], &Astart[0], num_nz, &Aindex[0], &Avalue[0]);
  REQUIRE(!return_bool);

  // Add column vectors to model with no rows returns OK
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], NULL, 0, NULL, NULL);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportLp(reference_lp, 2);
  HighsSetMessagelevel(ML_NONE);
}

