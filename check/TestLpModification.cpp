#include "catch.hpp"
#include "Highs.h"
#include "lp_data/HighsLpUtils.h"
#include "Avgas.h"

bool areLpColEqual(
		   const int num_col0, const double* colCost0, const double* colLower0, const double* colUpper0,
		   const int num_nz0, const int* Astart0, const int* Aindex0, const double* Avalue0,
		   const int num_col1, const double* colCost1, const double* colLower1, const double* colUpper1,
		   const int num_nz1, const int* Astart1, const int* Aindex1, const double* Avalue1,
		   const double infinite_bound) {
  if (num_col0 != num_col1) {
    printf("areLpColEqual: %d = num_col0 != num_col1 = %d\n", num_col0, num_col1);
    return false;
  }
  if (!num_col0) return true;
  int num_col = num_col0;
  for (int col = 0; col < num_col; col++) {
    if (colCost0[col] != colCost1[col]) {
      printf("areLpColEqual: %g = colCost0[%d] != colCost1[%d] = %g\n", colCost0[col], col, col, colCost1[col]);
      return false;
    }
  }
  for (int col = 0; col < num_col; col++) {
    if (colLower0[col] <= -infinite_bound && colLower1[col] <= -infinite_bound) continue;
    if (colLower0[col] != colLower1[col]) {
      printf("areLpColEqual: %g = colLower0[%d] != colLower1[%d] = %g\n", colLower0[col], col, col, colLower1[col]);
      return false;
    }
    if (colUpper0[col] >= infinite_bound && colUpper1[col] >= infinite_bound) continue;
    if (colUpper0[col] != colUpper1[col]) {
      printf("areLpColEqual: %g = colUpper0[%d] != colUpper1[%d] = %g\n", colUpper0[col], col, col, colUpper1[col]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    printf("areLpColEqual: %d = num_nz0 != num_nz1 = %d\n", num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (int col = 0; col < num_col; col++) {
    if (Astart0[col] != Astart1[col]) {
      printf("areLpColEqual: %d = Astart0[%d] != Astart1[%d] = %d\n", Astart0[col], col, col, Astart1[col]);
      return false;
    }
  }
  int num_nz = num_nz0;
  for (int nz = 0; nz < num_nz; nz++) {
    if (Aindex0[nz] != Aindex1[nz]) {
      printf("areLpColEqual: %d = Aindex0[%d] != Aindex1[%d] = %d\n", Aindex0[nz], nz, nz, Aindex1[nz]);
      return false;
    }
    if (Avalue0[nz] != Avalue1[nz]) {
      printf("areLpColEqual: %g = Avalue0[%d] != Avalue1[%d] = %g\n", Avalue0[nz], nz, nz, Avalue1[nz]);
      return false;
    }
  }
  return true;		      
}

bool areLpRowEqual(
		   const int num_row0,const double* rowLower0, const double* rowUpper0,
		   const int num_nz0, const int* ARstart0, const int* ARindex0, const double* ARvalue0,
		   const int num_row1, const double* rowLower1, const double* rowUpper1,
		   const int num_nz1, const int* ARstart1, const int* ARindex1, const double* ARvalue1,
		   const double infinite_bound) {
  if (num_row0 != num_row1) {
    printf("areLpRowEqual: %d = num_row0 != num_row1 = %d\n", num_row0, num_row1);
    return false;
  }
  if (!num_row0) return true;
  int num_row = num_row0;
  for (int row = 0; row < num_row; row++) {
    if (rowLower0[row] <= -infinite_bound && rowLower1[row] <= -infinite_bound) continue;
    if (rowLower0[row] != rowLower1[row]) {
      printf("areLpRowEqual: %g = rowLower0[%d] != rowLower1[%d] = %g\n", rowLower0[row], row, row, rowLower1[row]);
      return false;
    }
    if (rowUpper0[row] >= infinite_bound && rowUpper1[row] >= infinite_bound) continue;
    if (rowUpper0[row] != rowUpper1[row]) {
      printf("areLpRowEqual: %g = rowUpper0[%d] != rowUpper1[%d] = %g\n", rowUpper0[row], row, row, rowUpper1[row]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    printf("areLpRowEqual: %d = num_nz0 != num_nz1 = %d\n", num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (int row = 0; row < num_row; row++) {
    if (ARstart0[row] != ARstart1[row]) {
      printf("areLpRowEqual: %d = ARstart0[%d] != ARstart1[%d] = %d\n", ARstart0[row], row, row, ARstart1[row]);
      return false;
    }
  }
  int num_nz = num_nz0;
  for (int nz = 0; nz < num_nz; nz++) {
    if (ARindex0[nz] != ARindex1[nz]) {
      printf("areLpRowEqual: %d = ARindex0[%d] != ARindex1[%d] = %d\n", ARindex0[nz], nz, nz, ARindex1[nz]);
      return false;
    }
    if (ARvalue0[nz] != ARvalue1[nz]) {
      printf("areLpRowEqual: %g = ARvalue0[%d] != ARvalue1[%d] = %g\n", ARvalue0[nz], nz, nz, ARvalue1[nz]);
      return false;
    }
  }
  return true;		      
}


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
    update_out_in_ix(row_dim,
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
  const int avgas_num_col = 8;
  const int avgas_num_row = 10;
  int num_row = 0;
  int num_row_nz = 0;
  vector<double> rowLower;
  vector<double> rowUpper;
  vector<int> ARstart;
  vector<int> ARindex;
  vector<double> ARvalue;
  
  for (int row = 0; row < avgas_num_row; row++) {
    avgas.row(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex, ARvalue);
  }

  int num_col = 0;
  int num_col_nz = 0;
  vector<double> colCost;
  vector<double> colLower;
  vector<double> colUpper;
  vector<int> Astart;
  vector<int> Aindex;
  vector<double> Avalue;
  for (int col = 0; col < avgas_num_col; col++) {
    avgas.col(col, num_col, num_col_nz, colCost, colLower, colUpper, Astart, Aindex, Avalue);
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

  // Adding column vectors and matrix to model with no rows returns an error
  bool return_bool;
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], &Astart[0], num_col_nz, &Aindex[0], &Avalue[0]);
  REQUIRE(!return_bool);

  // Adding column vectors to model with no rows returns OK
  return_bool = highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], NULL, 0, NULL, NULL);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportLp(reference_lp, 2);
  HighsSetMessagelevel(ML_NONE);

  // Adding row vectors and matrix to model with columns returns OK
  return_bool = highs.addRows(num_row, &rowLower[0], &rowUpper[0], &ARstart[0], num_row_nz, &ARindex[0], &ARvalue[0]);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportLp(reference_lp, 2);
  HighsSetMessagelevel(ML_NONE);

  return_status = highs.run();
  HighsStatusReport("highs.run()", return_status);
  REQUIRE(return_status == HighsStatus::Optimal);
  
  HighsSetMessagelevel(ML_ALWAYS);
  reportLp(reference_lp, 2);
  HighsSetMessagelevel(ML_NONE);

  // Getting columns from the LP is OK
  int col1357_col_mask[] = {0, 1, 0, 1, 0, 1, 0, 1};
  int col1357_col_set[] = {1, 3, 5, 7};
  int col1357_num_ix = 4;
  int col1357_num_col;
  int col1357_num_nz;
  double *col1357_cost = (double *)malloc(sizeof(double) * col1357_num_ix);
  double *col1357_lower = (double *)malloc(sizeof(double) * col1357_num_ix);
  double *col1357_upper = (double *)malloc(sizeof(double) * col1357_num_ix);
  int *col1357_start = (int *)malloc(sizeof(int) * col1357_num_ix);
  int *col1357_index = (int *)malloc(sizeof(int) * num_col_nz);
  double *col1357_value = (double *)malloc(sizeof(double) * num_col_nz);

    
  return_bool = highs.getCols(col1357_num_ix, col1357_col_set, col1357_num_col, col1357_cost, col1357_lower, col1357_upper,
			      col1357_num_nz, col1357_start, col1357_index, col1357_value);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportMtx("Get by set\nColumn", col1357_num_col, col1357_num_nz, col1357_start, col1357_index, col1357_value);
  HighsSetMessagelevel(ML_NONE);
  
  return_bool = highs.getCols(3, 7, col1357_num_col, col1357_cost, col1357_lower, col1357_upper,
			      col1357_num_nz, col1357_start, col1357_index, col1357_value);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportMtx("Get by interval\nColumn", col1357_num_col, col1357_num_nz, col1357_start, col1357_index, col1357_value);
  HighsSetMessagelevel(ML_NONE);
  
  return_bool = highs.getCols(col1357_col_mask, col1357_num_col, col1357_cost, col1357_lower, col1357_upper,
			      col1357_num_nz, col1357_start, col1357_index, col1357_value);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportMtx("Get by mask\nColumn", col1357_num_col, col1357_num_nz, col1357_start, col1357_index, col1357_value);
  HighsSetMessagelevel(ML_NONE);
  
  return_bool = highs.deleteCols(col1357_num_ix, col1357_col_set);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportLp(reference_lp, 2);
  HighsSetMessagelevel(ML_NONE);

  // Getting columns from the LP is OK
  int row1357_row_set[] = {1, 3, 5, 7};
  int row1357_num_ix = 4;
  int row1357_num_row;
  int row1357_num_nz;
  double *row1357_lower = (double *)malloc(sizeof(double) * row1357_num_ix);
  double *row1357_upper = (double *)malloc(sizeof(double) * row1357_num_ix);
  int *row1357_start = (int *)malloc(sizeof(int) * row1357_num_ix);
  int *row1357_index = (int *)malloc(sizeof(int) * num_row_nz);
  double *row1357_value = (double *)malloc(sizeof(double) * num_row_nz);

    
  return_bool = highs.getRows(row1357_num_ix, row1357_row_set, row1357_num_row, row1357_lower, row1357_upper,
			      row1357_num_nz, row1357_start, row1357_index, row1357_value);
  REQUIRE(return_bool);

  HighsSetMessagelevel(ML_ALWAYS);
  reportMtx("Row   ", row1357_num_row, row1357_num_nz, row1357_start, row1357_index, row1357_value);
  HighsSetMessagelevel(ML_NONE);
  
}

