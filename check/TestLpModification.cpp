#include "Avgas.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"

const bool dev_run = true;

void HighsStatusReport(const HighsLogOptions& log_options, std::string message,
                       HighsStatus status) {
  highsLogUser(log_options, HighsLogType::INFO, "%s: HighsStatus = %d - %s\n",
               message.c_str(), (int)status,
               HighsStatusToString(status).c_str());
}

void callRun(Highs& highs, const HighsLogOptions& log_options,
             std::string message, const HighsStatus require_return_status) {
  HighsStatus return_status = highs.run();
  HighsStatusReport(log_options, message, return_status);
  REQUIRE(return_status == require_return_status);
#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis(message);
#endif
}

bool areLpColEqual(const int num_col0, const double* colCost0,
                   const double* colLower0, const double* colUpper0,
                   const int num_nz0, const int* Astart0, const int* Aindex0,
                   const double* Avalue0, const int num_col1,
                   const double* colCost1, const double* colLower1,
                   const double* colUpper1, const int num_nz1,
                   const int* Astart1, const int* Aindex1,
                   const double* Avalue1, const double infinite_bound) {
  if (num_col0 != num_col1) {
    if (dev_run)
      printf("areLpColEqual: %d = num_col0 != num_col1 = %d\n", num_col0,
             num_col1);
    return false;
  }
  if (!num_col0) return true;
  int num_col = num_col0;
  for (int col = 0; col < num_col; col++) {
    if (colCost0[col] != colCost1[col]) {
      if (dev_run)
        printf("areLpColEqual: %g = colCost0[%d] != colCost1[%d] = %g\n",
               colCost0[col], col, col, colCost1[col]);
      return false;
    }
  }
  for (int col = 0; col < num_col; col++) {
    if (colLower0[col] <= -infinite_bound && colLower1[col] <= -infinite_bound)
      continue;
    if (colLower0[col] != colLower1[col]) {
      if (dev_run)
        printf("areLpColEqual: %g = colLower0[%d] != colLower1[%d] = %g\n",
               colLower0[col], col, col, colLower1[col]);
      return false;
    }
    if (colUpper0[col] >= infinite_bound && colUpper1[col] >= infinite_bound)
      continue;
    if (colUpper0[col] != colUpper1[col]) {
      if (dev_run)
        printf("areLpColEqual: %g = colUpper0[%d] != colUpper1[%d] = %g\n",
               colUpper0[col], col, col, colUpper1[col]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    if (dev_run)
      printf("areLpColEqual: %d = num_nz0 != num_nz1 = %d\n", num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (int col = 0; col < num_col; col++) {
    if (Astart0[col] != Astart1[col]) {
      if (dev_run)
        printf("areLpColEqual: %d = Astart0[%d] != Astart1[%d] = %d\n",
               Astart0[col], col, col, Astart1[col]);
      return false;
    }
  }
  int num_nz = num_nz0;
  for (int nz = 0; nz < num_nz; nz++) {
    if (Aindex0[nz] != Aindex1[nz]) {
      if (dev_run)
        printf("areLpColEqual: %d = Aindex0[%d] != Aindex1[%d] = %d\n",
               Aindex0[nz], nz, nz, Aindex1[nz]);
      return false;
    }
    if (Avalue0[nz] != Avalue1[nz]) {
      if (dev_run)
        printf("areLpColEqual: %g = Avalue0[%d] != Avalue1[%d] = %g\n",
               Avalue0[nz], nz, nz, Avalue1[nz]);
      return false;
    }
  }
  return true;
}

bool areLpRowEqual(const int num_row0, const double* rowLower0,
                   const double* rowUpper0, const int num_nz0,
                   const int* ARstart0, const int* ARindex0,
                   const double* ARvalue0, const int num_row1,
                   const double* rowLower1, const double* rowUpper1,
                   const int num_nz1, const int* ARstart1, const int* ARindex1,
                   const double* ARvalue1, const double infinite_bound) {
  if (num_row0 != num_row1) {
    if (dev_run)
      printf("areLpRowEqual: %d = num_row0 != num_row1 = %d\n", num_row0,
             num_row1);
    return false;
  }
  if (!num_row0) return true;
  int num_row = num_row0;
  for (int row = 0; row < num_row; row++) {
    if (rowLower0[row] <= -infinite_bound && rowLower1[row] <= -infinite_bound)
      continue;
    if (rowLower0[row] != rowLower1[row]) {
      if (dev_run)
        printf("areLpRowEqual: %g = rowLower0[%d] != rowLower1[%d] = %g\n",
               rowLower0[row], row, row, rowLower1[row]);
      return false;
    }
    if (rowUpper0[row] >= infinite_bound && rowUpper1[row] >= infinite_bound)
      continue;
    if (rowUpper0[row] != rowUpper1[row]) {
      if (dev_run)
        printf("areLpRowEqual: %g = rowUpper0[%d] != rowUpper1[%d] = %g\n",
               rowUpper0[row], row, row, rowUpper1[row]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    if (dev_run)
      printf("areLpRowEqual: %d = num_nz0 != num_nz1 = %d\n", num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (int row = 0; row < num_row; row++) {
    if (ARstart0[row] != ARstart1[row]) {
      if (dev_run)
        printf("areLpRowEqual: %d = ARstart0[%d] != ARstart1[%d] = %d\n",
               ARstart0[row], row, row, ARstart1[row]);
      return false;
    }
  }
  int num_nz = num_nz0;
  for (int nz = 0; nz < num_nz; nz++) {
    if (ARindex0[nz] != ARindex1[nz]) {
      if (dev_run)
        printf("areLpRowEqual: %d = ARindex0[%d] != ARindex1[%d] = %d\n",
               ARindex0[nz], nz, nz, ARindex1[nz]);
      return false;
    }
    if (ARvalue0[nz] != ARvalue1[nz]) {
      if (dev_run)
        printf("areLpRowEqual: %g = ARvalue0[%d] != ARvalue1[%d] = %g\n",
               ARvalue0[nz], nz, nz, ARvalue1[nz]);
      return false;
    }
  }
  return true;
}

bool areLpEqual(const HighsLp lp0, const HighsLp lp1,
                const double infinite_bound) {
  bool return_bool;
  if (lp0.numCol_ > 0 && lp1.numCol_ > 0) {
    int lp0_num_nz = lp0.Astart_[lp0.numCol_];
    int lp1_num_nz = lp1.Astart_[lp1.numCol_];
    return_bool = areLpColEqual(
        lp0.numCol_, &lp0.colCost_[0], &lp0.colLower_[0], &lp0.colUpper_[0],
        lp0_num_nz, &lp0.Astart_[0], &lp0.Aindex_[0], &lp0.Avalue_[0],
        lp1.numCol_, &lp1.colCost_[0], &lp1.colLower_[0], &lp1.colUpper_[0],
        lp1_num_nz, &lp1.Astart_[0], &lp1.Aindex_[0], &lp1.Avalue_[0],
        infinite_bound);
    if (!return_bool) return return_bool;
  }
  if (lp0.numRow_ > 0 && lp1.numRow_ > 0) {
    int lp0_num_nz = 0;
    int lp1_num_nz = 0;
    return_bool = areLpRowEqual(
        lp0.numRow_, &lp0.rowLower_[0], &lp0.rowUpper_[0], lp0_num_nz, NULL,
        NULL, NULL, lp1.numRow_, &lp1.rowLower_[0], &lp1.rowUpper_[0],
        lp1_num_nz, NULL, NULL, NULL, infinite_bound);
  }
  return return_bool;
}

void testDeleteKeep(const HighsIndexCollection& index_collection) {
  int delete_from_index;
  int delete_to_index;
  int keep_from_index;
  int keep_to_index;
  int current_set_entry;
  const int* set = index_collection.set_;
  const int* mask = index_collection.mask_;
  const int dimension = index_collection.dimension_;
  if (dev_run) {
    if (index_collection.is_interval_) {
      printf("With index interval [%d, %d] in [%d, %d]\n",
             index_collection.from_, index_collection.to_, 0, dimension - 1);
    } else if (index_collection.is_set_) {
      printf("With index set\n");
      for (int entry = 0; entry < index_collection.set_num_entries_; entry++)
        printf(" %2d", entry);
      printf("\n");
      for (int entry = 0; entry < index_collection.set_num_entries_; entry++)
        printf(" %2d", set[entry]);
      printf("\n");
    } else {
      printf("With index mask\n");
      for (int index = 0; index < dimension; index++) printf(" %2d", index);
      printf("\n");
      for (int index = 0; index < dimension; index++)
        printf(" %2d", mask[index]);
      printf("\n");
    }
  }

  keep_from_index = 0;
  if (index_collection.is_interval_) {
    keep_to_index = index_collection.from_ - 1;
  } else if (index_collection.is_set_) {
    current_set_entry = 0;
    keep_to_index = set[0] - 1;
  } else {
    keep_to_index = dimension;
    for (int index = 0; index < dimension; index++) {
      if (mask[index]) {
        keep_to_index = index - 1;
        break;
      }
    }
  }
  if (dev_run) printf("Keep   [%2d, %2d]\n", 0, keep_to_index);
  if (keep_to_index >= dimension - 1) return;
  for (int k = 0; k < dimension; k++) {
    updateIndexCollectionOutInIndex(index_collection, delete_from_index,
                                    delete_to_index, keep_from_index,
                                    keep_to_index, current_set_entry);
    if (dev_run)
      printf("Delete [%2d, %2d]; keep [%2d, %2d]\n", delete_from_index,
             delete_to_index, keep_from_index, keep_to_index);
    if (delete_to_index >= dimension - 1 || keep_to_index >= dimension - 1)
      break;
  }
}

bool testAllDeleteKeep(int num_row) {
  // Test the extraction of intervals from index collections
  int set[] = {1, 4, 5, 8};
  int mask[] = {0, 1, 0, 0, 1, 1, 0, 0, 1, 0};

  HighsIndexCollection index_collection;
  index_collection.dimension_ = num_row;
  index_collection.is_interval_ = false;
  index_collection.from_ = 3;
  index_collection.to_ = 6;
  index_collection.is_set_ = false;
  index_collection.set_num_entries_ = 4;
  index_collection.set_ = &set[0];
  index_collection.is_mask_ = false;
  index_collection.mask_ = &mask[0];

  int save_from = index_collection.from_;
  int save_set_0 = set[0];
  int save_mask_0 = mask[0];

  int to_pass = 2;  // 2
  for (int pass = 0; pass <= to_pass; pass++) {
    if (dev_run) printf("\nTesting delete-keep: pass %d\n", pass);
    if (pass == 1) {
      // Mods to test LH limit behaviour
      index_collection.from_ = 0;
      set[0] = 0;
      mask[0] = 1;
    } else if (pass == 2) {
      // Mods to test RH limit behaviour
      index_collection.from_ = save_from;
      index_collection.to_ = 9;
      set[0] = save_set_0;
      set[3] = 9;
      mask[0] = save_mask_0;
      mask[9] = 1;
    }

    index_collection.is_interval_ = true;
    testDeleteKeep(index_collection);
    index_collection.is_interval_ = false;

    index_collection.is_set_ = true;
    testDeleteKeep(index_collection);
    index_collection.is_set_ = false;

    index_collection.is_mask_ = true;
    testDeleteKeep(index_collection);
  }
  return true;
}

void messageReportLp(const char* message, const HighsLp& lp) {
  HighsLogOptions log_options;
  bool output_flag;
  bool log_to_console;
  int log_dev_level;
  output_flag = dev_run;
  log_to_console = true;
  log_dev_level = LOG_DEV_LEVEL_VERBOSE;
  log_options.output_flag = &output_flag;
  log_options.log_file_stream = NULL;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  highsLogDev(log_options, HighsLogType::VERBOSE, "\nReporting LP: %s\n",
              message);
  reportLp(log_options, lp, HighsLogType::VERBOSE);
}

void messageReportMatrix(const char* message, const int num_col,
                         const int num_nz, const int* start, const int* index,
                         const double* value) {
  HighsLogOptions log_options;
  bool output_flag = true;
  bool log_to_console = false;
  int log_dev_level = LOG_DEV_LEVEL_INFO;
  log_options.log_file_stream = stdout;
  log_options.output_flag = &output_flag;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  highsLogDev(log_options, HighsLogType::VERBOSE, "\nReporting Matrix: %s\n",
              message);
  reportMatrix(log_options, message, num_col, num_nz, start, index, value);
}

// No commas in test case name.
TEST_CASE("LP-modification", "[highs_data]") {
  if (dev_run) printf("testAllDeleteKeep\n");
  testAllDeleteKeep(10);

  HighsOptions options;
  options.log_dev_level = LOG_DEV_LEVEL_VERBOSE;

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
    avgas.row(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex,
              ARvalue);
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
    avgas.col(col, num_col, num_col_nz, colCost, colLower, colUpper, Astart,
              Aindex, Avalue);
  }

  HighsStatus return_status;
  HighsModelStatus model_status;

  // Create two empty LPs: one to be initialised as AVGAS by adding
  // all the columns and rows separately, the other to be built by
  // adding piecemeal.
  HighsLp avgas_lp;
  HighsLp lp;

  Highs avgas_highs(options);
  if (!dev_run) {
    avgas_highs.setHighsOptionValue("output_flag", false);
  }
  return_status = avgas_highs.passModel(avgas_lp);
  HighsStatusReport(options.log_options, "avgas_highs.passModel(avgas_lp)",
                    return_status);
  REQUIRE(return_status == HighsStatus::OK);

  REQUIRE(avgas_highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              0, NULL, NULL, NULL));
  REQUIRE(avgas_highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                              &ARstart[0], &ARindex[0], &ARvalue[0]));

  return_status = avgas_highs.writeModel("");
  HighsStatusReport(options.log_options, "avgas_highs.writeModel(\"\")",
                    return_status);
  REQUIRE(return_status == HighsStatus::OK);

  Highs highs(options);
  if (!dev_run) {
    highs.setHighsOptionValue("output_flag", false);
  }
  return_status = highs.setHighsOptionValue("highs_debug_level", 2);
  HighsStatusReport(options.log_options, "\"highs_debug_level\", 2",
                    return_status);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.passModel(lp);
  HighsStatusReport(options.log_options, "highs.passModel(lp)", return_status);
  REQUIRE(return_status == HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::NOTSET);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::MODEL_EMPTY);

  // Adding column vectors and matrix to model with no rows returns an error
  REQUIRE(!highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                         num_col_nz, &Astart[0], &Aindex[0], &Avalue[0]));

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  return_status = highs.writeModel("");
  HighsStatusReport(options.log_options, "highs.writeModel(\"\")",
                    return_status);
  REQUIRE(return_status == HighsStatus::OK);

  // Adding row vectors and matrix to model with columns returns OK
  REQUIRE(highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                        &ARstart[0], &ARindex[0], &ARvalue[0]));

  return_status = highs.writeModel("");
  HighsStatusReport(options.log_options, "highs.writeModel(\"\")",
                    return_status);
  REQUIRE(return_status == HighsStatus::OK);

  REQUIRE(
      areLpEqual(highs.getLp(), avgas_highs.getLp(), options.infinite_bound));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  double avgas_optimal_objective_value;
  highs.getHighsInfoValue("objective_function_value",
                          avgas_optimal_objective_value);
  double optimal_objective_value;

  // Getting columns from the LP is OK
  int col1357_col_mask[] = {0, 1, 0, 1, 0, 1, 0, 1};
  int col1357_col_set[] = {1, 3, 5, 7};
  int col1357_illegal_col_set[] = {3, 7, 1, 5};
  int col1357_num_ix = 4;
  int col1357_num_col;
  int col1357_num_nz;
  double* col1357_cost = (double*)malloc(sizeof(double) * col1357_num_ix);
  double* col1357_lower = (double*)malloc(sizeof(double) * col1357_num_ix);
  double* col1357_upper = (double*)malloc(sizeof(double) * col1357_num_ix);
  int* col1357_start = (int*)malloc(sizeof(int) * col1357_num_ix);
  int* col1357_index = (int*)malloc(sizeof(int) * num_col_nz);
  double* col1357_value = (double*)malloc(sizeof(double) * num_col_nz);

  REQUIRE(highs.getCols(3, 6, col1357_num_col, col1357_cost, col1357_lower,
                        col1357_upper, col1357_num_nz, col1357_start,
                        col1357_index, col1357_value));

  REQUIRE(!highs.getCols(col1357_num_ix, col1357_illegal_col_set,
                         col1357_num_col, col1357_cost, col1357_lower,
                         col1357_upper, col1357_num_nz, col1357_start,
                         col1357_index, col1357_value));

  REQUIRE(highs.getCols(col1357_num_ix, col1357_col_set, col1357_num_col,
                        col1357_cost, col1357_lower, col1357_upper,
                        col1357_num_nz, col1357_start, col1357_index,
                        col1357_value));

  REQUIRE(highs.getCols(col1357_col_mask, col1357_num_col, col1357_cost,
                        col1357_lower, col1357_upper, col1357_num_nz,
                        col1357_start, col1357_index, col1357_value));

  // Try to delete an empty range of cols: OK
  REQUIRE(highs.deleteCols(0, -1));

  // Try to delete more cols than there are: ERROR
  REQUIRE(!highs.deleteCols(0, num_col + 1));

  REQUIRE(highs.deleteCols(col1357_num_ix, col1357_col_set));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After deleting columns 1, 3, 5, 7");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.addCols(col1357_num_col, col1357_cost, col1357_lower,
                        col1357_upper, col1357_num_nz, col1357_start,
                        col1357_index, col1357_value));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After restoring columns 1, 3, 5, 7\n");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  highs.getHighsInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After re-solving");
#endif

  // Delete all the columns: OK
  REQUIRE(highs.deleteCols(0, num_col - 1));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After deleting all columns");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Delete all the rows: OK
  REQUIRE(highs.deleteRows(0, num_row - 1));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After deleteRows(0, num_row - 1)");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Adding row vectors and matrix to model with columns returns OK
  REQUIRE(highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                        &ARstart[0], &ARindex[0], &ARvalue[0]));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("With columns but and rows");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Getting rows from the LP is OK
  int from_row_ix = 0;
  int to_row_ix = 3;
  int row0135789_row_set[] = {0, 1, 3, 5, 7, 8, 9};
  int row0135789_row_mask[] = {1, 1, 0, 1, 0, 1, 0, 1, 1, 1};
  int row0135789_num_ix = 7;
  int row0135789_num_row;
  int row0135789_num_nz;
  double* row0135789_lower =
      (double*)malloc(sizeof(double) * row0135789_num_ix);
  double* row0135789_upper =
      (double*)malloc(sizeof(double) * row0135789_num_ix);
  int* row0135789_start = (int*)malloc(sizeof(int) * row0135789_num_ix);
  int* row0135789_index = (int*)malloc(sizeof(int) * num_row_nz);
  double* row0135789_value = (double*)malloc(sizeof(double) * num_row_nz);

  REQUIRE(highs.getRows(from_row_ix, to_row_ix, row0135789_num_row,
                        row0135789_lower, row0135789_upper, row0135789_num_nz,
                        row0135789_start, row0135789_index, row0135789_value));

  REQUIRE(highs.getRows(row0135789_num_ix, row0135789_row_set,
                        row0135789_num_row, row0135789_lower, row0135789_upper,
                        row0135789_num_nz, row0135789_start, row0135789_index,
                        row0135789_value));

  REQUIRE(highs.getRows(row0135789_row_mask, row0135789_num_row,
                        row0135789_lower, row0135789_upper, row0135789_num_nz,
                        row0135789_start, row0135789_index, row0135789_value));

  REQUIRE(highs.getRows(row0135789_num_ix, row0135789_row_set,
                        row0135789_num_row, row0135789_lower, row0135789_upper,
                        row0135789_num_nz, row0135789_start, row0135789_index,
                        row0135789_value));

  REQUIRE(highs.deleteRows(row0135789_num_ix, row0135789_row_set));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After deleting rows 0-1, 3, 5, 7-9");
#endif

  int row012_row_set[] = {0, 1, 2};
  int row012_row_mask[] = {1, 1, 1};
  int row012_num_ix = 3;
  int row012_num_row;
  int row012_num_nz;
  double* row012_lower = (double*)malloc(sizeof(double) * row012_num_ix);
  double* row012_upper = (double*)malloc(sizeof(double) * row012_num_ix);
  int* row012_start = (int*)malloc(sizeof(int) * row012_num_ix);
  int* row012_index = (int*)malloc(sizeof(int) * num_row_nz);
  double* row012_value = (double*)malloc(sizeof(double) * num_row_nz);

  REQUIRE(highs.getRows(row012_num_ix, row012_row_set, row012_num_row,
                        row012_lower, row012_upper, row012_num_nz, row012_start,
                        row012_index, row012_value));

  REQUIRE(highs.deleteRows(row012_row_mask));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After deleting rows 0-2");
#endif

  // Delete all the columns: OK
  REQUIRE(highs.deleteCols(0, num_col - 1));

#ifdef HiGHSDEV
  messageReportLp("After deleting all columns", highs.getLp());
  highs.reportModelStatusSolutionBasis("After deleting all columns");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Can't add rows with no columns
  REQUIRE(!highs.addRows(row0135789_num_row, row0135789_lower, row0135789_upper,
                         row0135789_num_nz, row0135789_start, row0135789_index,
                         row0135789_value));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.addRows(row0135789_num_row, row0135789_lower, row0135789_upper,
                        row0135789_num_nz, row0135789_start, row0135789_index,
                        row0135789_value));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.addRows(row012_num_row, row012_lower, row012_upper,
                        row012_num_nz, row012_start, row012_index,
                        row012_value));

#ifdef HiGHSDEV
  messageReportLp("After restoring all rows", highs.getLp());
  highs.reportModelStatusSolutionBasis("After restoring all rows");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  highs.getHighsInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After resolve");
#endif

  // Try to delete an empty range of rows: OK
  REQUIRE(highs.deleteRows(0, -1));

  // Try to delete more rows than there are: ERROR
  REQUIRE(!highs.deleteRows(0, num_row));

  REQUIRE(highs.getCols(col1357_col_mask, col1357_num_col, col1357_cost,
                        col1357_lower, col1357_upper, col1357_num_nz,
                        col1357_start, col1357_index, col1357_value));

  REQUIRE(highs.deleteCols(col1357_num_ix, col1357_col_set));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  int col0123_col_mask[] = {1, 1, 1, 1};
  //  int col0123_col_set[] = {0, 1, 2, 3};
  int col0123_num_ix = 4;
  int col0123_num_col;
  int col0123_num_nz;
  double* col0123_cost = (double*)malloc(sizeof(double) * col0123_num_ix);
  double* col0123_lower = (double*)malloc(sizeof(double) * col0123_num_ix);
  double* col0123_upper = (double*)malloc(sizeof(double) * col0123_num_ix);
  int* col0123_start = (int*)malloc(sizeof(int) * col0123_num_ix);
  int* col0123_index = (int*)malloc(sizeof(int) * num_col_nz);
  double* col0123_value = (double*)malloc(sizeof(double) * num_col_nz);

  REQUIRE(highs.getCols(col0123_col_mask, col0123_num_col, col0123_cost,
                        col0123_lower, col0123_upper, col0123_num_nz,
                        col0123_start, col0123_index, col0123_value));
  //  messageReportMatrix("Get col1357 by mask\nRow   ", col1357_num_col,
  //  col1357_num_nz, col1357_start, col1357_index, col1357_value);
  //  messageReportMatrix("Get col0123 by mask\nRow   ", col0123_num_col,
  //  col0123_num_nz, col0123_start, col0123_index, col0123_value);

  REQUIRE(highs.deleteRows(0, num_row - 1));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.deleteCols(col0123_col_mask));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After deleting all rows and columns");
#endif

  // Adding row vectors to model with no columns returns OK
  REQUIRE(highs.addRows(row0135789_num_row, row0135789_lower, row0135789_upper,
                        0, NULL, NULL, NULL));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After restoring 7 rows");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.addRows(row012_num_row, row012_lower, row012_upper, 0,
                        row012_start, row012_index, row012_value));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After restoring all rows");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.addCols(col1357_num_col, col1357_cost, col1357_lower,
                        col1357_upper, col1357_num_nz, col1357_start,
                        col1357_index, col1357_value));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After restoring columns 1, 3, 5, 7");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis(
      "After solving after restoring all rows and columns 1, 3, 5, 7");
#endif

  REQUIRE(highs.addCols(col0123_num_col, col0123_cost, col0123_lower,
                        col0123_upper, col0123_num_nz, col0123_start,
                        col0123_index, col0123_value));

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After restoring columns 0-3");
#endif

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  highs.getHighsInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

  // Fix columns 1, 3, 5, 7 to check resetting of their nonbasic status
  col1357_lower[0] = 0;
  col1357_lower[1] = 0;
  col1357_lower[2] = 0;
  col1357_lower[3] = 0;
  col1357_upper[0] = 0;
  col1357_upper[1] = 0;
  col1357_upper[2] = 0;
  col1357_upper[3] = 0;

  REQUIRE(highs.changeColsBounds(col1357_num_ix, col1357_col_set, col1357_lower,
                                 col1357_upper));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Now restore the upper bounds to check resetting of their nonbasic status
  col1357_upper[0] = 1;
  col1357_upper[1] = 1;
  col1357_upper[2] = 1;
  col1357_upper[3] = 1;

  REQUIRE(highs.changeColsBounds(col1357_num_ix, col1357_col_set, col1357_lower,
                                 col1357_upper));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  highs.getHighsInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

  const HighsLp& local_lp = highs.getLp();
  row0135789_lower[0] = local_lp.rowLower_[0];
  row0135789_lower[1] = local_lp.rowLower_[1];
  row0135789_lower[2] = local_lp.rowLower_[3];
  row0135789_lower[3] = local_lp.rowLower_[5];
  row0135789_lower[4] = local_lp.rowLower_[7];
  row0135789_lower[5] = local_lp.rowLower_[8];
  row0135789_lower[6] = local_lp.rowLower_[9];
  row0135789_upper[0] = local_lp.rowLower_[0];
  row0135789_upper[1] = local_lp.rowLower_[1];
  row0135789_upper[2] = local_lp.rowLower_[3];
  row0135789_upper[3] = local_lp.rowLower_[5];
  row0135789_upper[4] = local_lp.rowLower_[7];
  row0135789_upper[5] = local_lp.rowLower_[8];
  row0135789_upper[6] = local_lp.rowLower_[9];

  REQUIRE(highs.changeRowsBounds(row0135789_num_ix, row0135789_row_set,
                                 row0135789_lower, row0135789_upper));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  row0135789_upper[0] = local_lp.rowUpper_[0];
  row0135789_upper[1] = local_lp.rowUpper_[1];
  row0135789_upper[2] = local_lp.rowUpper_[3];
  row0135789_upper[3] = local_lp.rowUpper_[5];
  row0135789_upper[4] = local_lp.rowUpper_[7];
  row0135789_upper[5] = local_lp.rowUpper_[8];
  row0135789_upper[6] = local_lp.rowUpper_[9];

  REQUIRE(highs.changeRowsBounds(row0135789_num_ix, row0135789_row_set,
                                 row0135789_lower, row0135789_upper));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.deleteRows(0, num_row - 1));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.deleteCols(0, num_col - 1));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("After deleteing all rows and columns");
#endif

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

#ifdef HiGHSDEV
  highs.reportModelStatusSolutionBasis("With columns but no rows");
#endif

  // Adding row vectors and matrix to model with columns returns OK
  REQUIRE(highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                        &ARstart[0], &ARindex[0], &ARvalue[0]));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  col1357_cost[0] = 2.01;
  col1357_cost[1] = 2.31;
  col1357_cost[2] = 2.51;
  col1357_cost[3] = 2.71;
  col1357_lower[0] = 0.01;
  col1357_lower[1] = 0.31;
  col1357_lower[2] = 0.51;
  col1357_lower[3] = 0.71;
  col1357_upper[0] = 1.01;
  col1357_upper[1] = 1.31;
  col1357_upper[2] = 1.51;
  col1357_upper[3] = 1.71;

  row0135789_lower[0] = -0.01;
  row0135789_lower[1] = -0.11;
  row0135789_lower[2] = -0.31;
  row0135789_lower[3] = -0.51;
  row0135789_lower[4] = -0.71;
  row0135789_lower[5] = -0.81;
  row0135789_lower[6] = -0.91;
  row0135789_upper[0] = 3.01;
  row0135789_upper[1] = 3.11;
  row0135789_upper[2] = 3.31;
  row0135789_upper[3] = 3.51;
  row0135789_upper[4] = 3.71;
  row0135789_upper[5] = 3.81;
  row0135789_upper[6] = 3.91;

  // Attempting to set a cost to infinity may return error
  REQUIRE(highs.changeColCost(7, HIGHS_CONST_INF) == allow_infinite_costs);

  // Attempting to set a cost to a finite value returns OK
  REQUIRE(highs.changeColCost(7, 77));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.changeColsCost(col1357_num_ix, col1357_col_set, col1357_cost));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Attempting to set row bounds with infinite lower bound returns error
  REQUIRE(!highs.changeRowBounds(2, HIGHS_CONST_INF, 3.21));

  REQUIRE(highs.changeRowBounds(2, -HIGHS_CONST_INF, 3.21));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  // Attempting to set col bounds with -infinite upper bound returns error
  REQUIRE(!highs.changeColBounds(2, 0.21, -HIGHS_CONST_INF));

  REQUIRE(highs.changeColBounds(2, 0.21, HIGHS_CONST_INF));

  REQUIRE(highs.changeRowsBounds(row0135789_num_ix, row0135789_row_set,
                                 row0135789_lower, row0135789_upper));

  REQUIRE(highs.changeColsBounds(col1357_num_ix, col1357_col_set, col1357_lower,
                                 col1357_upper));

  // Return the LP to its original state with a mask
  REQUIRE(highs.changeColsCost(col1357_col_mask, &colCost[0]));

  REQUIRE(highs.changeColBounds(2, colLower[2], colUpper[2]));

  REQUIRE(highs.changeColsBounds(col1357_col_mask, &colLower[0], &colUpper[0]));

  REQUIRE(
      highs.changeRowsBounds(row0135789_row_mask, &rowLower[0], &rowUpper[0]));

  REQUIRE(highs.changeRowBounds(2, rowLower[2], rowUpper[2]));

  REQUIRE(
      areLpEqual(avgas_highs.getLp(), highs.getLp(), options.infinite_bound));

  int before_num_col;
  int after_num_col;
  int rm_col;
  int before_num_row;
  int after_num_row;
  int rm_row;

  before_num_col = highs.getNumCols();
  rm_col = 0;
  REQUIRE(highs.deleteCols(rm_col, rm_col));
  after_num_col = highs.getNumCols();
  if (dev_run)
    printf("After removing col %d / %d have %d cols\n", rm_col, before_num_col,
           after_num_col);
  REQUIRE(after_num_col == before_num_col - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  before_num_row = highs.getNumRows();
  rm_row = 0;
  REQUIRE(highs.deleteRows(rm_row, rm_row));
  after_num_row = highs.getNumRows();
  if (dev_run)
    printf("After removing row %d / %d have %d rows\n", rm_row, before_num_row,
           after_num_row);
  REQUIRE(after_num_row == before_num_row - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  before_num_col = highs.getNumCols();
  rm_col = before_num_col - 1;
  REQUIRE(highs.deleteCols(rm_col, rm_col));
  after_num_col = highs.getNumCols();
  if (dev_run)
    printf("After removing col %d / %d have %d cols\n", rm_col, before_num_col,
           after_num_col);
  REQUIRE(after_num_col == before_num_col - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  before_num_row = highs.getNumRows();
  rm_row = before_num_row - 1;
  REQUIRE(highs.deleteRows(rm_row, rm_row));
  after_num_row = highs.getNumRows();
  if (dev_run)
    printf("After removing row %d / %d have %d rows\n", rm_row, before_num_row,
           after_num_row);
  REQUIRE(after_num_row == before_num_row - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(!highs.scaleCol(-1, 2.0));

  REQUIRE(!highs.scaleCol(highs.getNumCols(), 2.0));

  REQUIRE(!highs.scaleCol(0, 0));

  REQUIRE(highs.scaleCol(highs.getNumCols() - 1, 2.0));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.scaleCol(0, -2.0));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(!highs.scaleRow(-1, 2.0));

  REQUIRE(!highs.scaleRow(highs.getNumRows(), 2.0));

  REQUIRE(!highs.scaleRow(0, 0));

  REQUIRE(highs.scaleRow(0, 2.0));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  REQUIRE(highs.scaleRow(highs.getNumRows() - 1, -2.0));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);
}

TEST_CASE("LP-getcols", "[highs_data]") {
  Highs highs;
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  int aindex[2] = {0, 1};
  double avalue[2] = {1.0, -1.0};
  highs.addRow(0.0, 0.0, 2, aindex, avalue);
  int num_cols;
  int num_nz;
  int matrix_start[2] = {-1, -1};
  highs.getCols(0, 1, num_cols, NULL, NULL, NULL, num_nz, matrix_start, NULL,
                NULL);
  REQUIRE(num_cols == 2);
  REQUIRE(num_nz == 2);
  REQUIRE(matrix_start[0] == 0);
  REQUIRE(matrix_start[1] == 1);
  int matrix_indices[2] = {-1, -1};
  double matrix_values[2] = {0.0, 0.0};
  highs.getCols(0, 1, num_cols, NULL, NULL, NULL, num_nz, matrix_start,
                matrix_indices, matrix_values);
  REQUIRE(matrix_indices[0] == 0);
  REQUIRE(matrix_indices[1] == 0);
  REQUIRE(matrix_values[0] == 1.0);
  REQUIRE(matrix_values[1] == -1.0);
}

TEST_CASE("LP-getrows", "[highs_data]") {
  Highs highs;
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  int aindex = 0;
  double avalue = 1.0;
  highs.addRow(0.0, 0.0, 1, &aindex, &avalue);
  aindex = 1;
  avalue = -2.0;
  highs.addRow(0.0, 0.0, 1, &aindex, &avalue);
  int num_rows;
  int num_nz;
  int matrix_start[2] = {-1, -1};
  highs.getRows(0, 1, num_rows, NULL, NULL, num_nz, matrix_start, NULL, NULL);
  REQUIRE(num_rows == 2);
  REQUIRE(num_nz == 2);
  REQUIRE(matrix_start[0] == 0);
  REQUIRE(matrix_start[1] == 1);
  int matrix_indices[2] = {-1, -1};
  double matrix_values[2] = {0.0, 0.0};
  highs.getRows(0, 1, num_rows, NULL, NULL, num_nz, matrix_start,
                matrix_indices, matrix_values);
  REQUIRE(matrix_indices[0] == 0);
  REQUIRE(matrix_indices[1] == 1);
  REQUIRE(matrix_values[0] == 1.0);
  REQUIRE(matrix_values[1] == -2.0);
}

TEST_CASE("LP-interval-changes", "[highs_data]") {

  HighsStatus run_status;
  HighsStatus return_status;

  Highs highs;
  const HighsOptions& options = highs.getHighsOptions();
  const HighsInfo& info = highs.getHighsInfo();

  highs.setHighsOptionValue("output_flag", dev_run);
  highs.setHighsOptionValue("log_to_console", true);
  highs.setHighsOptionValue("log_dev_level", LOG_DEV_LEVEL_VERBOSE);

  std::string model_file = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);

  const HighsLp& lp = highs.getLp();

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  double avgas_optimal_objective_function_value = info.objective_function_value;

  REQUIRE(info.objective_function_value == avgas_optimal_objective_function_value);
  //  messageReportLp("LP-interval-changes", lp);

  // Change an interval of column costs
  const int from_col = 2;
  const int to_col = 5;
  int get_num_col;
  int get_num_nz;
  vector<double> og_col2345_cost;
  vector<double> set_col2345_cost;
  vector<double> get_col2345_cost;
  og_col2345_cost.resize(lp.numCol_);
  set_col2345_cost.resize(lp.numCol_);
  get_col2345_cost.resize(lp.numCol_);
  set_col2345_cost[0] = 2.0;
  set_col2345_cost[1] = 3.0;
  set_col2345_cost[2] = 4.0;
  set_col2345_cost[3] = 5.0;
  REQUIRE(highs.getCols(from_col, to_col, get_num_col, &og_col2345_cost[0], NULL, NULL, get_num_nz, NULL, NULL, NULL));
  REQUIRE(highs.changeColsCost(from_col, to_col, &set_col2345_cost[0]));
  REQUIRE(highs.getCols(from_col, to_col, get_num_col, &get_col2345_cost[0], NULL, NULL, get_num_nz, NULL, NULL, NULL));
  for (int iCol = from_col; iCol < to_col + 1; iCol++)
    REQUIRE(get_col2345_cost[iCol] == set_col2345_cost[iCol]);
  REQUIRE(highs.changeColsCost(from_col, to_col, &og_col2345_cost[0]));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::OK);

  double optimal_objective_function_value;
  highs.getHighsInfoValue("objective_function_value", optimal_objective_function_value);
  REQUIRE(optimal_objective_function_value == avgas_optimal_objective_function_value);
}
