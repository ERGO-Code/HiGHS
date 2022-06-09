#include "Avgas.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsRandom.h"
#include "util/HighsUtils.h"

const bool dev_run = false;
const double inf = kHighsInf;
const double double_equal_tolerance = 1e-5;
void HighsStatusReport(const HighsLogOptions& log_options, std::string message,
                       HighsStatus status);

void callRun(Highs& highs, const HighsLogOptions& log_options,
             std::string message, const HighsStatus require_return_status);

bool areLpColEqual(const HighsInt num_col0, const double* colCost0,
                   const double* colLower0, const double* colUpper0,
                   const HighsInt num_nz0, const HighsInt* Astart0,
                   const HighsInt* Aindex0, const double* Avalue0,
                   const HighsInt num_col1, const double* colCost1,
                   const double* colLower1, const double* colUpper1,
                   const HighsInt num_nz1, const HighsInt* Astart1,
                   const HighsInt* Aindex1, const double* Avalue1,
                   const double infinite_bound);

bool areLpRowEqual(const HighsInt num_row0, const double* rowLower0,
                   const double* rowUpper0, const HighsInt num_nz0,
                   const HighsInt* ARstart0, const HighsInt* ARindex0,
                   const double* ARvalue0, const HighsInt num_row1,
                   const double* rowLower1, const double* rowUpper1,
                   const HighsInt num_nz1, const HighsInt* ARstart1,
                   const HighsInt* ARindex1, const double* ARvalue1,
                   const double infinite_bound);

bool areLpEqual(const HighsLp lp0, const HighsLp lp1,
                const double infinite_bound);

void testDeleteKeep(const HighsIndexCollection& index_collection);

bool testAllDeleteKeep(HighsInt num_row);

void messageReportLp(const char* message, const HighsLp& lp);

void messageReportMatrix(const char* message, const HighsInt num_col,
                         const HighsInt num_nz, const HighsInt* start,
                         const HighsInt* index, const double* value);

// No commas in test case name.
TEST_CASE("LP-717-od", "[highs_data]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  REQUIRE(highs.addCol(0.0, -inf, inf, 0, nullptr, nullptr) ==
          HighsStatus::kOk);
  std::vector<HighsInt> index = {0};
  std::vector<double> value = {1.0};
  REQUIRE(highs.addRow(2.0, inf, 1, &index[0], &value[0]) == HighsStatus::kOk);
  REQUIRE(highs.addCol(0.0, -inf, inf, 0, nullptr, nullptr) ==
          HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
}

TEST_CASE("LP-717-full0", "[highs_data]") {
  // Add columns to an LP with a small matrix held row-wise so that
  // the orientation is flipped
  HighsInt row_block_num_col = 2;
  HighsInt row_block_num_row = 3;
  HighsInt col_block_num_col = 3;
  HighsInt col_block_num_row = 3;

  HighsLp lp;
  lp.num_col_ = row_block_num_col + col_block_num_col;
  lp.num_row_ = row_block_num_row;
  lp.col_cost_ = {-2, -1, -2, -3, -3};
  lp.col_lower_ = {0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 1};
  lp.row_lower_ = {-inf, -inf, -inf};
  lp.row_upper_ = {-2, 2, 1};
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 3, 8, 10};
  lp.a_matrix_.index_ = {0, 2, 4, 0, 1, 2, 3, 4, 1, 3};
  lp.a_matrix_.value_ = {1, -1, -3, 1, 1, 1, -2, 3, 1, 2};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsLp& highs_lp = highs.getLp();
  highs.passModel(lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  HighsInfo info0 = highs.getInfo();
  HighsSolution solution0 = highs.getSolution();
  highs.clear();
  if (!dev_run) highs.setOptionValue("output_flag", false);
  std::vector<double> row_block_col_cost;
  std::vector<double> row_block_col_lower;
  std::vector<double> row_block_col_upper;
  std::vector<double> row_block_row_lower;
  std::vector<double> row_block_row_upper;

  std::vector<double> col_block_col_cost;
  std::vector<double> col_block_col_lower;
  std::vector<double> col_block_col_upper;
  std::vector<double> col_block_row_lower;
  std::vector<double> col_block_row_upper;

  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (iCol < row_block_num_col) {
      row_block_col_cost.push_back(lp.col_cost_[iCol]);
      row_block_col_lower.push_back(lp.col_lower_[iCol]);
      row_block_col_upper.push_back(lp.col_upper_[iCol]);
    } else {
      col_block_col_cost.push_back(lp.col_cost_[iCol]);
      col_block_col_lower.push_back(lp.col_lower_[iCol]);
      col_block_col_upper.push_back(lp.col_upper_[iCol]);
    }
  }
  row_block_row_lower = lp.row_lower_;
  row_block_row_upper = lp.row_upper_;

  HighsInt row_block_format = (HighsInt)MatrixFormat::kRowwise;
  HighsInt row_block_num_nz;
  std::vector<HighsInt> row_block_start;
  std::vector<HighsInt> row_block_index;
  std::vector<double> row_block_value;

  for (HighsInt iRow = 0; iRow < row_block_num_row; iRow++) {
    row_block_start.push_back(row_block_index.size());
    for (HighsInt iEl = lp.a_matrix_.start_[iRow];
         iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
      HighsInt iCol = lp.a_matrix_.index_[iEl];
      if (iCol < row_block_num_col) {
        row_block_index.push_back(iCol);
        row_block_value.push_back(lp.a_matrix_.value_[iEl]);
      }
    }
  }
  row_block_num_nz = row_block_index.size();

  REQUIRE(highs.addCols(row_block_num_col, &row_block_col_cost[0],
                        &row_block_col_lower[0], &row_block_col_upper[0], 0,
                        nullptr, nullptr, nullptr) == HighsStatus::kOk);

  REQUIRE(highs.addRows(row_block_num_row, &row_block_row_lower[0],
                        &row_block_row_upper[0], row_block_num_nz,
                        &row_block_start[0], &row_block_index[0],
                        &row_block_value[0]) == HighsStatus::kOk);

  if (dev_run)
    printf("After adding a row-wise matrix, LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);

  // Column block
  HighsInt col_block_format = (HighsInt)MatrixFormat::kColwise;
  HighsInt col_block_num_nz = 6;
  std::vector<HighsInt> col_block_start = {0, 2, 4};
  std::vector<HighsInt> col_block_index = {0, 1, 1, 2, 0, 1};
  std::vector<double> col_block_value = {-1, 1, -2, 2, -3, 3};
  REQUIRE(highs.addCols(col_block_num_col, &col_block_col_cost[0],
                        &col_block_col_lower[0], &col_block_col_upper[0],
                        col_block_num_nz, &col_block_start[0],
                        &col_block_index[0],
                        &col_block_value[0]) == HighsStatus::kOk);
  if (dev_run)
    printf("After adding a column-wise matrix, LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);

  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  if (dev_run)
    printf("After run() LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);
}

TEST_CASE("LP-717-full1", "[highs_data]") {
  // Add columns to an LP with a larger matrix held row-wise so that
  // the orientation is not flipped
  HighsInt row_block_num_col = 5;
  HighsInt row_block_num_row = 3;
  HighsInt col_block_num_col = 3;
  HighsInt col_block_num_row = 3;

  HighsLp lp;
  lp.num_col_ = row_block_num_col + col_block_num_col;
  lp.num_row_ = row_block_num_row;
  lp.col_cost_ = {-1, -1, -1, -1, -2, -2, -3, -3};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 1, 1, 1, 1};
  lp.row_lower_ = {-inf, -inf, -inf};
  lp.row_upper_ = {-5, 5, 1};
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 5, 13, 16};
  lp.a_matrix_.index_ = {0, 2, 4, 5, 7, 0, 1, 2, 3, 4, 5, 6, 7, 1, 3, 6};
  lp.a_matrix_.value_ = {1, -1, 1, -1, -3, 1, 1, 1, 1, 1, 1, -2, 3, 1, -1, 2};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsLp& highs_lp = highs.getLp();
  highs.passModel(lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  HighsInfo info0 = highs.getInfo();
  HighsSolution solution0 = highs.getSolution();
  highs.clear();
  if (!dev_run) highs.setOptionValue("output_flag", false);
  std::vector<double> row_block_col_cost;
  std::vector<double> row_block_col_lower;
  std::vector<double> row_block_col_upper;
  std::vector<double> row_block_row_lower;
  std::vector<double> row_block_row_upper;

  std::vector<double> col_block_col_cost;
  std::vector<double> col_block_col_lower;
  std::vector<double> col_block_col_upper;
  std::vector<double> col_block_row_lower;
  std::vector<double> col_block_row_upper;

  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (iCol < row_block_num_col) {
      row_block_col_cost.push_back(lp.col_cost_[iCol]);
      row_block_col_lower.push_back(lp.col_lower_[iCol]);
      row_block_col_upper.push_back(lp.col_upper_[iCol]);
    } else {
      col_block_col_cost.push_back(lp.col_cost_[iCol]);
      col_block_col_lower.push_back(lp.col_lower_[iCol]);
      col_block_col_upper.push_back(lp.col_upper_[iCol]);
    }
  }
  row_block_row_lower = lp.row_lower_;
  row_block_row_upper = lp.row_upper_;

  HighsInt row_block_format = (HighsInt)MatrixFormat::kRowwise;
  HighsInt row_block_num_nz;
  std::vector<HighsInt> row_block_start;
  std::vector<HighsInt> row_block_index;
  std::vector<double> row_block_value;

  for (HighsInt iRow = 0; iRow < row_block_num_row; iRow++) {
    row_block_start.push_back(row_block_index.size());
    for (HighsInt iEl = lp.a_matrix_.start_[iRow];
         iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
      HighsInt iCol = lp.a_matrix_.index_[iEl];
      if (iCol < row_block_num_col) {
        row_block_index.push_back(iCol);
        row_block_value.push_back(lp.a_matrix_.value_[iEl]);
      }
    }
  }
  row_block_num_nz = row_block_index.size();

  REQUIRE(highs.addCols(row_block_num_col, &row_block_col_cost[0],
                        &row_block_col_lower[0], &row_block_col_upper[0], 0,
                        nullptr, nullptr, nullptr) == HighsStatus::kOk);

  REQUIRE(highs.addRows(row_block_num_row, &row_block_row_lower[0],
                        &row_block_row_upper[0], row_block_num_nz,
                        &row_block_start[0], &row_block_index[0],
                        &row_block_value[0]) == HighsStatus::kOk);

  if (dev_run)
    printf("After adding a row-wise matrix, LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);

  // Column block
  HighsInt col_block_format = (HighsInt)MatrixFormat::kColwise;
  HighsInt col_block_num_nz = 6;
  std::vector<HighsInt> col_block_start = {0, 2, 4};
  std::vector<HighsInt> col_block_index = {0, 1, 1, 2, 0, 1};
  std::vector<double> col_block_value = {-1, 1, -2, 2, -3, 3};
  REQUIRE(highs.addCols(col_block_num_col, &col_block_col_cost[0],
                        &col_block_col_lower[0], &col_block_col_upper[0],
                        col_block_num_nz, &col_block_start[0],
                        &col_block_index[0],
                        &col_block_value[0]) == HighsStatus::kOk);
  if (dev_run)
    printf("After adding a column-wise matrix, LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);

  const bool equal_lp = lp == highs_lp;
  REQUIRE(equal_lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  if (dev_run)
    printf("After run() LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);
}

TEST_CASE("LP-717-full2", "[highs_data]") {
  // Add columns and then rows to an LP with a larger matrix held
  // row-wise so that the orientation is not flipped
  HighsInt row_block_num_col = 5;
  HighsInt row_block_num_row = 3;
  HighsInt col_block_num_col = 3;
  HighsInt col_block_num_row = 3;

  HighsLp lp;
  lp.num_col_ = row_block_num_col + col_block_num_col;
  lp.num_row_ = 2 * row_block_num_row;
  lp.col_cost_ = {-1, -1, -1, -1, -2, -2, -3, -3};
  lp.col_lower_ = {0, 0, 0, 0, 0, 0, 0, 0};
  lp.col_upper_ = {1, 1, 1, 1, 1, 1, 1, 1};
  lp.row_lower_ = {-inf, -inf, -inf};
  lp.row_upper_ = {-1, 6, 2};
  for (HighsInt iRow = 0; iRow < row_block_num_row; iRow++) {
    lp.row_lower_.push_back(lp.row_lower_[iRow]);
    lp.row_upper_.push_back(lp.row_upper_[iRow]);
  }
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = lp.num_row_;
  lp.a_matrix_.start_ = {0, 5, 13, 16, 19, 24, 26};
  lp.a_matrix_.index_ = {0, 2, 4, 5, 7, 0, 1, 2, 3, 4, 5, 6, 7,
                         1, 3, 6, 0, 2, 4, 0, 1, 2, 3, 4, 1, 3};
  lp.a_matrix_.value_ = {1, -1, 1, -1, -3, 1, 1, 1, 1, 1, 1, -2, 3,
                         1, -1, 2, 1,  -1, 1, 1, 1, 1, 1, 1, 1,  -1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsLp& highs_lp = highs.getLp();
  highs.passModel(lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  HighsInfo info0 = highs.getInfo();
  HighsSolution solution0 = highs.getSolution();
  highs.clear();
  if (!dev_run) highs.setOptionValue("output_flag", false);
  std::vector<double> row_block_col_cost;
  std::vector<double> row_block_col_lower;
  std::vector<double> row_block_col_upper;
  std::vector<double> row_block_row_lower;
  std::vector<double> row_block_row_upper;

  std::vector<double> col_block_col_cost;
  std::vector<double> col_block_col_lower;
  std::vector<double> col_block_col_upper;
  std::vector<double> col_block_row_lower;
  std::vector<double> col_block_row_upper;

  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    if (iCol < row_block_num_col) {
      row_block_col_cost.push_back(lp.col_cost_[iCol]);
      row_block_col_lower.push_back(lp.col_lower_[iCol]);
      row_block_col_upper.push_back(lp.col_upper_[iCol]);
    } else {
      col_block_col_cost.push_back(lp.col_cost_[iCol]);
      col_block_col_lower.push_back(lp.col_lower_[iCol]);
      col_block_col_upper.push_back(lp.col_upper_[iCol]);
    }
  }
  row_block_row_lower = lp.row_lower_;
  row_block_row_upper = lp.row_upper_;

  HighsInt row_block_format = (HighsInt)MatrixFormat::kRowwise;
  HighsInt row_block_num_nz;
  std::vector<HighsInt> row_block_start;
  std::vector<HighsInt> row_block_index;
  std::vector<double> row_block_value;

  for (HighsInt iRow = 0; iRow < row_block_num_row; iRow++) {
    row_block_start.push_back(row_block_index.size());
    for (HighsInt iEl = lp.a_matrix_.start_[iRow];
         iEl < lp.a_matrix_.start_[iRow + 1]; iEl++) {
      HighsInt iCol = lp.a_matrix_.index_[iEl];
      if (iCol < row_block_num_col) {
        row_block_index.push_back(iCol);
        row_block_value.push_back(lp.a_matrix_.value_[iEl]);
      }
    }
  }
  row_block_num_nz = row_block_index.size();

  REQUIRE(highs.addCols(row_block_num_col, &row_block_col_cost[0],
                        &row_block_col_lower[0], &row_block_col_upper[0], 0,
                        nullptr, nullptr, nullptr) == HighsStatus::kOk);

  REQUIRE(highs.addRows(row_block_num_row, &row_block_row_lower[0],
                        &row_block_row_upper[0], row_block_num_nz,
                        &row_block_start[0], &row_block_index[0],
                        &row_block_value[0]) == HighsStatus::kOk);

  if (dev_run)
    printf("After adding a row-wise matrix, LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);

  // Column block
  HighsInt col_block_format = (HighsInt)MatrixFormat::kColwise;
  HighsInt col_block_num_nz = 6;
  std::vector<HighsInt> col_block_start = {0, 2, 4};
  std::vector<HighsInt> col_block_index = {0, 1, 1, 2, 0, 1};
  std::vector<double> col_block_value = {-1, 1, -2, 2, -3, 3};
  REQUIRE(highs.addCols(col_block_num_col, &col_block_col_cost[0],
                        &col_block_col_lower[0], &col_block_col_upper[0],
                        col_block_num_nz, &col_block_start[0],
                        &col_block_index[0],
                        &col_block_value[0]) == HighsStatus::kOk);
  if (dev_run)
    printf("After adding a column-wise matrix, LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);

  REQUIRE(highs.addRows(row_block_num_row, &row_block_row_lower[0],
                        &row_block_row_upper[0], row_block_num_nz,
                        &row_block_start[0], &row_block_index[0],
                        &row_block_value[0]) == HighsStatus::kOk);

  if (dev_run)
    printf("After adding a row-wise matrix, LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);

  const bool equal_lp = lp == highs_lp;
  REQUIRE(equal_lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  if (dev_run)
    printf("After run() LP matrix has format %d\n",
           (int)highs_lp.a_matrix_.format_);
}

TEST_CASE("LP-modification", "[highs_data]") {
  if (dev_run) printf("testAllDeleteKeep\n");
  testAllDeleteKeep(10);

  HighsOptions options;
  //  options.log_dev_level = kHighsLogDevLevelVerbose;

  Avgas avgas;
  const HighsInt avgas_num_col = 8;
  const HighsInt avgas_num_row = 10;
  HighsInt num_row = 0;
  HighsInt num_row_nz = 0;
  std::vector<double> rowLower;
  std::vector<double> rowUpper;
  std::vector<HighsInt> ARstart;
  std::vector<HighsInt> ARindex;
  std::vector<double> ARvalue;

  for (HighsInt row = 0; row < avgas_num_row; row++) {
    avgas.row(row, num_row, num_row_nz, rowLower, rowUpper, ARstart, ARindex,
              ARvalue);
  }

  HighsInt num_col = 0;
  HighsInt num_col_nz = 0;
  std::vector<double> colCost;
  std::vector<double> colLower;
  std::vector<double> colUpper;
  std::vector<HighsInt> Astart;
  std::vector<HighsInt> Aindex;
  std::vector<double> Avalue;
  for (HighsInt col = 0; col < avgas_num_col; col++) {
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

  Highs avgas_highs;
  avgas_highs.passOptions(options);
  if (!dev_run) avgas_highs.setOptionValue("output_flag", false);
  return_status = avgas_highs.passModel(avgas_lp);
  HighsStatusReport(options.log_options, "avgas_highs.passModel(avgas_lp)",
                    return_status);
  REQUIRE(return_status == HighsStatus::kOk);

  REQUIRE(avgas_highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                              0, NULL, NULL, NULL) == HighsStatus::kOk);
  REQUIRE(avgas_highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                              &ARstart[0], &ARindex[0],
                              &ARvalue[0]) == HighsStatus::kOk);

  //  return_status = avgas_highs.writeModel("");

  Highs highs;
  highs.passOptions(options);
  if (!dev_run) highs.setOptionValue("output_flag", false);
  return_status = highs.setOptionValue("highs_debug_level", 3);
  REQUIRE(return_status == HighsStatus::kOk);

  lp.model_name_ = "Building avgas";
  return_status = highs.passModel(lp);
  REQUIRE(return_status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kNotset);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kModelEmpty);

  // Adding column vectors and matrix to model with no rows returns an error
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0],
                        num_col_nz, &Astart[0], &Aindex[0],
                        &Avalue[0]) == HighsStatus::kError);

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  //  return_status = highs.writeModel("");

  // Adding row vectors and matrix to model with columns returns OK
  REQUIRE(highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                        &ARstart[0], &ARindex[0],
                        &ARvalue[0]) == HighsStatus::kOk);

  //  return_status = highs.writeModel("");

  REQUIRE(
      areLpEqual(highs.getLp(), avgas_highs.getLp(), options.infinite_bound));

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  double avgas_optimal_objective_value;
  highs.getInfoValue("objective_function_value", avgas_optimal_objective_value);
  double optimal_objective_value;

  // Getting columns from the LP is OK
  HighsInt col1357_col_mask[] = {0, 1, 0, 1, 0, 1, 0, 1};
  HighsInt col1357_col_set[] = {1, 3, 5, 7};
  HighsInt col1357_illegal_col_set[] = {3, 7, 1, 5};
  HighsInt col1357_num_ix = 4;
  HighsInt col1357_num_col;
  HighsInt col1357_num_nz;
  double* col1357_cost = (double*)malloc(sizeof(double) * col1357_num_ix);
  double* col1357_lower = (double*)malloc(sizeof(double) * col1357_num_ix);
  double* col1357_upper = (double*)malloc(sizeof(double) * col1357_num_ix);
  HighsInt* col1357_start =
      (HighsInt*)malloc(sizeof(HighsInt) * col1357_num_ix);
  HighsInt* col1357_index = (HighsInt*)malloc(sizeof(HighsInt) * num_col_nz);
  double* col1357_value = (double*)malloc(sizeof(double) * num_col_nz);

  REQUIRE(highs.getCols(3, 6, col1357_num_col, col1357_cost, col1357_lower,
                        col1357_upper, col1357_num_nz, col1357_start,
                        col1357_index, col1357_value) == HighsStatus::kOk);

  REQUIRE(highs.getCols(col1357_num_ix, col1357_illegal_col_set,
                        col1357_num_col, col1357_cost, col1357_lower,
                        col1357_upper, col1357_num_nz, col1357_start,
                        col1357_index, col1357_value) == HighsStatus::kError);

  REQUIRE(highs.getCols(col1357_num_ix, col1357_col_set, col1357_num_col,
                        col1357_cost, col1357_lower, col1357_upper,
                        col1357_num_nz, col1357_start, col1357_index,
                        col1357_value) == HighsStatus::kOk);

  REQUIRE(highs.getCols(col1357_col_mask, col1357_num_col, col1357_cost,
                        col1357_lower, col1357_upper, col1357_num_nz,
                        col1357_start, col1357_index,
                        col1357_value) == HighsStatus::kOk);

  // Try to delete an empty range of cols: OK
  REQUIRE(highs.deleteCols(0, -1) == HighsStatus::kOk);

  // Try to delete more cols than there are: ERROR
  REQUIRE(highs.deleteCols(0, num_col + 1) == HighsStatus::kError);

  REQUIRE(highs.deleteCols(col1357_num_ix, col1357_col_set) ==
          HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.addCols(col1357_num_col, col1357_cost, col1357_lower,
                        col1357_upper, col1357_num_nz, col1357_start,
                        col1357_index, col1357_value) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  highs.getInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

  // Delete all the columns: OK
  REQUIRE(highs.deleteCols(0, num_col - 1) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Delete all the rows: OK
  REQUIRE(highs.deleteRows(0, num_row - 1) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Adding row vectors and matrix to model with columns returns OK
  REQUIRE(highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                        &ARstart[0], &ARindex[0],
                        &ARvalue[0]) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Getting rows from the LP is OK
  HighsInt from_row_ix = 0;
  HighsInt to_row_ix = 3;
  HighsInt row0135789_row_set[] = {0, 1, 3, 5, 7, 8, 9};
  HighsInt row0135789_row_mask[] = {1, 1, 0, 1, 0, 1, 0, 1, 1, 1};
  HighsInt row0135789_num_ix = 7;
  HighsInt row0135789_num_row;
  HighsInt row0135789_num_nz;
  double* row0135789_lower =
      (double*)malloc(sizeof(double) * row0135789_num_ix);
  double* row0135789_upper =
      (double*)malloc(sizeof(double) * row0135789_num_ix);
  HighsInt* row0135789_start =
      (HighsInt*)malloc(sizeof(HighsInt) * row0135789_num_ix);
  HighsInt* row0135789_index = (HighsInt*)malloc(sizeof(HighsInt) * num_row_nz);
  double* row0135789_value = (double*)malloc(sizeof(double) * num_row_nz);

  REQUIRE(highs.getRows(from_row_ix, to_row_ix, row0135789_num_row,
                        row0135789_lower, row0135789_upper, row0135789_num_nz,
                        row0135789_start, row0135789_index,
                        row0135789_value) == HighsStatus::kOk);

  REQUIRE(highs.getRows(row0135789_num_ix, row0135789_row_set,
                        row0135789_num_row, row0135789_lower, row0135789_upper,
                        row0135789_num_nz, row0135789_start, row0135789_index,
                        row0135789_value) == HighsStatus::kOk);

  REQUIRE(highs.getRows(row0135789_row_mask, row0135789_num_row,
                        row0135789_lower, row0135789_upper, row0135789_num_nz,
                        row0135789_start, row0135789_index,
                        row0135789_value) == HighsStatus::kOk);

  REQUIRE(highs.getRows(row0135789_num_ix, row0135789_row_set,
                        row0135789_num_row, row0135789_lower, row0135789_upper,
                        row0135789_num_nz, row0135789_start, row0135789_index,
                        row0135789_value) == HighsStatus::kOk);

  REQUIRE(highs.deleteRows(row0135789_num_ix, row0135789_row_set) ==
          HighsStatus::kOk);

  HighsInt row012_row_set[] = {0, 1, 2};
  HighsInt row012_row_mask[] = {1, 1, 1};
  HighsInt row012_num_ix = 3;
  HighsInt row012_num_row;
  HighsInt row012_num_nz;
  double* row012_lower = (double*)malloc(sizeof(double) * row012_num_ix);
  double* row012_upper = (double*)malloc(sizeof(double) * row012_num_ix);
  HighsInt* row012_start = (HighsInt*)malloc(sizeof(HighsInt) * row012_num_ix);
  HighsInt* row012_index = (HighsInt*)malloc(sizeof(HighsInt) * num_row_nz);
  double* row012_value = (double*)malloc(sizeof(double) * num_row_nz);

  REQUIRE(highs.getRows(row012_num_ix, row012_row_set, row012_num_row,
                        row012_lower, row012_upper, row012_num_nz, row012_start,
                        row012_index, row012_value) == HighsStatus::kOk);

  REQUIRE(highs.deleteRows(row012_row_mask) == HighsStatus::kOk);

  // Delete all the columns: OK
  REQUIRE(highs.deleteCols(0, num_col - 1) == HighsStatus::kOk);

  messageReportLp("After deleting all columns", highs.getLp());

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Can't add rows with no columns
  REQUIRE(highs.addRows(row0135789_num_row, row0135789_lower, row0135789_upper,
                        row0135789_num_nz, row0135789_start, row0135789_index,
                        row0135789_value) == HighsStatus::kError);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.addRows(row0135789_num_row, row0135789_lower, row0135789_upper,
                        row0135789_num_nz, row0135789_start, row0135789_index,
                        row0135789_value) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  //  highs.setOptionValue("log_dev_level", 2);
  //  highs.setOptionValue("highs_debug_level", 3);
  //  highs.setOptionValue("output_flag", true);
  REQUIRE(highs.addRows(row012_num_row, row012_lower, row012_upper,
                        row012_num_nz, row012_start, row012_index,
                        row012_value) == HighsStatus::kOk);

  //  messageReportLp("After restoring all rows", highs.getLp());

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  //  highs.setOptionValue("log_dev_level", 0);
  //  highs.setOptionValue("highs_debug_level", 0);
  //  highs.setOptionValue("output_flag", false);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  highs.getInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(std::fabs(optimal_objective_value - avgas_optimal_objective_value) <
          double_equal_tolerance);

  // Try to delete an empty range of rows: OK
  REQUIRE(highs.deleteRows(0, -1) == HighsStatus::kOk);

  // Try to delete more rows than there are: ERROR
  REQUIRE(highs.deleteRows(0, num_row) == HighsStatus::kError);

  REQUIRE(highs.getCols(col1357_col_mask, col1357_num_col, col1357_cost,
                        col1357_lower, col1357_upper, col1357_num_nz,
                        col1357_start, col1357_index,
                        col1357_value) == HighsStatus::kOk);

  REQUIRE(highs.deleteCols(col1357_num_ix, col1357_col_set) ==
          HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  HighsInt col0123_col_mask[] = {1, 1, 1, 1};
  //  HighsInt col0123_col_set[] = {0, 1, 2, 3};
  HighsInt col0123_num_ix = 4;
  HighsInt col0123_num_col;
  HighsInt col0123_num_nz;
  double* col0123_cost = (double*)malloc(sizeof(double) * col0123_num_ix);
  double* col0123_lower = (double*)malloc(sizeof(double) * col0123_num_ix);
  double* col0123_upper = (double*)malloc(sizeof(double) * col0123_num_ix);
  HighsInt* col0123_start =
      (HighsInt*)malloc(sizeof(HighsInt) * col0123_num_ix);
  HighsInt* col0123_index = (HighsInt*)malloc(sizeof(HighsInt) * num_col_nz);
  double* col0123_value = (double*)malloc(sizeof(double) * num_col_nz);

  REQUIRE(highs.getCols(col0123_col_mask, col0123_num_col, col0123_cost,
                        col0123_lower, col0123_upper, col0123_num_nz,
                        col0123_start, col0123_index,
                        col0123_value) == HighsStatus::kOk);
  //  messageReportMatrix("Get col1357 by mask\nRow   ", col1357_num_col,
  //  col1357_num_nz, col1357_start, col1357_index, col1357_value);
  //  messageReportMatrix("Get col0123 by mask\nRow   ", col0123_num_col,
  //  col0123_num_nz, col0123_start, col0123_index, col0123_value);

  REQUIRE(highs.deleteRows(0, num_row - 1) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.deleteCols(col0123_col_mask) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Adding row vectors to model with no columns returns OK
  REQUIRE(highs.addRows(row0135789_num_row, row0135789_lower, row0135789_upper,
                        0, NULL, NULL, NULL) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.addRows(row012_num_row, row012_lower, row012_upper, 0,
                        row012_start, row012_index,
                        row012_value) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.addCols(col1357_num_col, col1357_cost, col1357_lower,
                        col1357_upper, col1357_num_nz, col1357_start,
                        col1357_index, col1357_value) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  REQUIRE(highs.addCols(col0123_num_col, col0123_cost, col0123_lower,
                        col0123_upper, col0123_num_nz, col0123_start,
                        col0123_index, col0123_value) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  highs.getInfoValue("objective_function_value", optimal_objective_value);
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
                                 col1357_upper) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Now restore the upper bounds to check resetting of their nonbasic status
  col1357_upper[0] = 1;
  col1357_upper[1] = 1;
  col1357_upper[2] = 1;
  col1357_upper[3] = 1;

  REQUIRE(highs.changeColsBounds(col1357_num_ix, col1357_col_set, col1357_lower,
                                 col1357_upper) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  highs.getInfoValue("objective_function_value", optimal_objective_value);
  REQUIRE(optimal_objective_value == avgas_optimal_objective_value);

  const HighsLp& local_lp = highs.getLp();
  row0135789_lower[0] = local_lp.row_lower_[0];
  row0135789_lower[1] = local_lp.row_lower_[1];
  row0135789_lower[2] = local_lp.row_lower_[3];
  row0135789_lower[3] = local_lp.row_lower_[5];
  row0135789_lower[4] = local_lp.row_lower_[7];
  row0135789_lower[5] = local_lp.row_lower_[8];
  row0135789_lower[6] = local_lp.row_lower_[9];
  row0135789_upper[0] = local_lp.row_lower_[0];
  row0135789_upper[1] = local_lp.row_lower_[1];
  row0135789_upper[2] = local_lp.row_lower_[3];
  row0135789_upper[3] = local_lp.row_lower_[5];
  row0135789_upper[4] = local_lp.row_lower_[7];
  row0135789_upper[5] = local_lp.row_lower_[8];
  row0135789_upper[6] = local_lp.row_lower_[9];

  REQUIRE(highs.changeRowsBounds(row0135789_num_ix, row0135789_row_set,
                                 row0135789_lower,
                                 row0135789_upper) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  row0135789_upper[0] = local_lp.row_upper_[0];
  row0135789_upper[1] = local_lp.row_upper_[1];
  row0135789_upper[2] = local_lp.row_upper_[3];
  row0135789_upper[3] = local_lp.row_upper_[5];
  row0135789_upper[4] = local_lp.row_upper_[7];
  row0135789_upper[5] = local_lp.row_upper_[8];
  row0135789_upper[6] = local_lp.row_upper_[9];

  REQUIRE(highs.changeRowsBounds(row0135789_num_ix, row0135789_row_set,
                                 row0135789_lower,
                                 row0135789_upper) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.deleteRows(0, num_row - 1) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.deleteCols(0, num_col - 1) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Adding column vectors to model with no rows returns OK
  REQUIRE(highs.addCols(num_col, &colCost[0], &colLower[0], &colUpper[0], 0,
                        NULL, NULL, NULL) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Adding row vectors and matrix to model with columns returns OK
  REQUIRE(highs.addRows(num_row, &rowLower[0], &rowUpper[0], num_row_nz,
                        &ARstart[0], &ARindex[0],
                        &ARvalue[0]) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

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
  return_status = highs.changeColCost(7, kHighsInf);
  if (kHighsAllowInfiniteCosts) {
    REQUIRE(return_status == HighsStatus::kOk);
  } else {
    REQUIRE(return_status == HighsStatus::kError);
  }

  // Attempting to set a cost to a finite value returns OK
  REQUIRE(highs.changeColCost(7, 77) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.changeColsCost(col1357_num_ix, col1357_col_set, col1357_cost) ==
          HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Attempting to set row bounds with infinite lower bound returns error
  REQUIRE(highs.changeRowBounds(2, kHighsInf, 3.21) == HighsStatus::kError);

  REQUIRE(highs.changeRowBounds(2, -kHighsInf, 3.21) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  // Attempting to set col bounds with -infinite upper bound returns error
  REQUIRE(highs.changeColBounds(2, 0.21, -kHighsInf) == HighsStatus::kError);

  REQUIRE(highs.changeColBounds(2, 0.21, kHighsInf) == HighsStatus::kOk);

  REQUIRE(highs.changeRowsBounds(row0135789_num_ix, row0135789_row_set,
                                 row0135789_lower,
                                 row0135789_upper) == HighsStatus::kOk);

  REQUIRE(highs.changeColsBounds(col1357_num_ix, col1357_col_set, col1357_lower,
                                 col1357_upper) == HighsStatus::kOk);

  // Return the LP to its original state with a mask
  REQUIRE(highs.changeColsCost(col1357_col_mask, &colCost[0]) ==
          HighsStatus::kOk);

  REQUIRE(highs.changeColBounds(2, colLower[2], colUpper[2]) ==
          HighsStatus::kOk);

  REQUIRE(highs.changeColsBounds(col1357_col_mask, &colLower[0],
                                 &colUpper[0]) == HighsStatus::kOk);

  REQUIRE(highs.changeRowsBounds(row0135789_row_mask, &rowLower[0],
                                 &rowUpper[0]) == HighsStatus::kOk);

  REQUIRE(highs.changeRowBounds(2, rowLower[2], rowUpper[2]) ==
          HighsStatus::kOk);

  avgas_highs.setMatrixFormat(MatrixFormat::kColwise);
  REQUIRE(
      areLpEqual(avgas_highs.getLp(), highs.getLp(), options.infinite_bound));

  HighsInt before_num_col;
  HighsInt after_num_col;
  HighsInt rm_col;
  HighsInt before_num_row;
  HighsInt after_num_row;
  HighsInt rm_row;

  before_num_col = highs.getNumCol();
  rm_col = 0;
  REQUIRE(highs.deleteCols(rm_col, rm_col) == HighsStatus::kOk);
  after_num_col = highs.getNumCol();
  if (dev_run)
    printf("After removing col %" HIGHSINT_FORMAT " / %" HIGHSINT_FORMAT
           " have %" HIGHSINT_FORMAT " cols\n",
           rm_col, before_num_col, after_num_col);
  REQUIRE(after_num_col == before_num_col - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  before_num_row = highs.getNumRow();
  rm_row = 0;
  REQUIRE(highs.deleteRows(rm_row, rm_row) == HighsStatus::kOk);
  after_num_row = highs.getNumRow();
  if (dev_run)
    printf("After removing row %" HIGHSINT_FORMAT " / %" HIGHSINT_FORMAT
           " have %" HIGHSINT_FORMAT " rows\n",
           rm_row, before_num_row, after_num_row);
  REQUIRE(after_num_row == before_num_row - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  before_num_col = highs.getNumCol();
  rm_col = before_num_col - 1;
  REQUIRE(highs.deleteCols(rm_col, rm_col) == HighsStatus::kOk);
  after_num_col = highs.getNumCol();
  if (dev_run)
    printf("After removing col %" HIGHSINT_FORMAT " / %" HIGHSINT_FORMAT
           " have %" HIGHSINT_FORMAT " cols\n",
           rm_col, before_num_col, after_num_col);
  REQUIRE(after_num_col == before_num_col - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  before_num_row = highs.getNumRow();
  rm_row = before_num_row - 1;
  REQUIRE(highs.deleteRows(rm_row, rm_row) == HighsStatus::kOk);
  after_num_row = highs.getNumRow();
  if (dev_run)
    printf("After removing row %" HIGHSINT_FORMAT " / %" HIGHSINT_FORMAT
           " have %" HIGHSINT_FORMAT " rows\n",
           rm_row, before_num_row, after_num_row);
  REQUIRE(after_num_row == before_num_row - 1);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.scaleCol(-1, 2.0) == HighsStatus::kError);

  REQUIRE(highs.scaleCol(highs.getNumCol(), 2.0) == HighsStatus::kError);

  REQUIRE(highs.scaleCol(0, 0) == HighsStatus::kError);

  REQUIRE(highs.scaleCol(highs.getNumCol() - 1, 2.0) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.scaleCol(0, -2.0) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.scaleRow(-1, 2.0) == HighsStatus::kError);

  REQUIRE(highs.scaleRow(highs.getNumRow(), 2.0) == HighsStatus::kError);

  REQUIRE(highs.scaleRow(0, 0) == HighsStatus::kError);

  REQUIRE(highs.scaleRow(0, 2.0) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  REQUIRE(highs.scaleRow(highs.getNumRow() - 1, -2.0) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);
}

TEST_CASE("LP-getcols", "[highs_data]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  HighsInt aindex[2] = {0, 1};
  double avalue[2] = {1.0, -1.0};
  highs.addRow(0.0, 0.0, 2, aindex, avalue);
  HighsInt num_cols;
  HighsInt num_nz;
  HighsInt matrix_start[2] = {-1, -1};
  highs.getCols(0, 1, num_cols, NULL, NULL, NULL, num_nz, matrix_start, NULL,
                NULL);
  REQUIRE(num_cols == 2);
  REQUIRE(num_nz == 2);
  REQUIRE(matrix_start[0] == 0);
  REQUIRE(matrix_start[1] == 1);
  HighsInt matrix_indices[2] = {-1, -1};
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
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  highs.addCol(-1.0, 0.0, 1.0, 0, NULL, NULL);
  HighsInt aindex = 0;
  double avalue = 1.0;
  highs.addRow(0.0, 0.0, 1, &aindex, &avalue);
  aindex = 1;
  avalue = -2.0;
  highs.addRow(0.0, 0.0, 1, &aindex, &avalue);
  HighsInt num_rows;
  HighsInt num_nz;
  HighsInt matrix_start[2] = {-1, -1};
  highs.getRows(0, 1, num_rows, NULL, NULL, num_nz, matrix_start, NULL, NULL);
  REQUIRE(num_rows == 2);
  REQUIRE(num_nz == 2);
  REQUIRE(matrix_start[0] == 0);
  REQUIRE(matrix_start[1] == 1);
  HighsInt matrix_indices[2] = {-1, -1};
  double matrix_values[2] = {0.0, 0.0};
  highs.getRows(0, 1, num_rows, NULL, NULL, num_nz, matrix_start,
                matrix_indices, matrix_values);
  REQUIRE(matrix_indices[0] == 0);
  REQUIRE(matrix_indices[1] == 1);
  REQUIRE(matrix_values[0] == 1.0);
  REQUIRE(matrix_values[1] == -2.0);
}

TEST_CASE("LP-interval-changes", "[highs_data]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  const HighsOptions& options = highs.getOptions();
  const HighsInfo& info = highs.getInfo();

  highs.setOptionValue("log_to_console", true);
  highs.setOptionValue("log_dev_level", kHighsLogDevLevelVerbose);

  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);

  const HighsLp& lp = highs.getLp();

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  double avgas_optimal_objective_function_value = info.objective_function_value;

  REQUIRE(info.objective_function_value ==
          avgas_optimal_objective_function_value);
  //  messageReportLp("LP-interval-changes", lp);

  // Change an interval of column costs
  HighsInt from_col = 2;
  HighsInt to_col = 5;
  HighsInt set_num_col = to_col - from_col + 1;
  HighsInt get_num_col;
  HighsInt get_num_nz;
  std::vector<double> og_col2345_cost;
  std::vector<double> set_col2345_cost;
  std::vector<double> get_col2345_cost;
  og_col2345_cost.resize(lp.num_col_);
  set_col2345_cost.resize(set_num_col);
  get_col2345_cost.resize(lp.num_col_);
  set_col2345_cost[0] = 2.0;
  set_col2345_cost[1] = 3.0;
  set_col2345_cost[2] = 4.0;
  set_col2345_cost[3] = 5.0;
  REQUIRE(highs.getCols(from_col, to_col, get_num_col, &og_col2345_cost[0],
                        NULL, NULL, get_num_nz, NULL, NULL,
                        NULL) == HighsStatus::kOk);
  REQUIRE(highs.changeColsCost(from_col, to_col, &set_col2345_cost[0]) ==
          HighsStatus::kOk);
  REQUIRE(highs.getCols(from_col, to_col, get_num_col, &get_col2345_cost[0],
                        NULL, NULL, get_num_nz, NULL, NULL,
                        NULL) == HighsStatus::kOk);
  REQUIRE(get_num_col == set_num_col);
  for (HighsInt usr_col = 0; usr_col < get_num_col; usr_col++)
    REQUIRE(get_col2345_cost[usr_col] == set_col2345_cost[usr_col]);
  REQUIRE(highs.changeColsCost(from_col, to_col, &og_col2345_cost[0]) ==
          HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  double optimal_objective_function_value;
  highs.getInfoValue("objective_function_value",
                     optimal_objective_function_value);
  REQUIRE(optimal_objective_function_value ==
          avgas_optimal_objective_function_value);

  from_col = 0;
  to_col = 4;
  set_num_col = to_col - from_col + 1;
  std::vector<double> og_col01234_lower;
  std::vector<double> og_col01234_upper;
  std::vector<double> set_col01234_lower;
  std::vector<double> get_col01234_lower;
  og_col01234_lower.resize(lp.num_col_);
  og_col01234_upper.resize(lp.num_col_);
  set_col01234_lower.resize(set_num_col);
  get_col01234_lower.resize(lp.num_col_);
  set_col01234_lower[0] = 0.0;
  set_col01234_lower[1] = 1.0;
  set_col01234_lower[2] = 2.0;
  set_col01234_lower[3] = 3.0;
  set_col01234_lower[4] = 4.0;
  REQUIRE(highs.getCols(from_col, to_col, get_num_col, NULL,
                        &og_col01234_lower[0], &og_col01234_upper[0],
                        get_num_nz, NULL, NULL, NULL) == HighsStatus::kOk);
  REQUIRE(highs.changeColsBounds(from_col, to_col, &set_col01234_lower[0],
                                 &og_col01234_upper[0]) == HighsStatus::kOk);
  REQUIRE(highs.getCols(from_col, to_col, get_num_col, NULL,
                        &get_col01234_lower[0], &og_col01234_upper[0],
                        get_num_nz, NULL, NULL, NULL) == HighsStatus::kOk);
  REQUIRE(get_num_col == set_num_col);
  for (HighsInt usr_col = 0; usr_col < get_num_col; usr_col++)
    REQUIRE(get_col01234_lower[usr_col] == set_col01234_lower[usr_col]);
  REQUIRE(highs.changeColsBounds(from_col, to_col, &og_col01234_lower[0],
                                 &og_col01234_upper[0]) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  highs.getInfoValue("objective_function_value",
                     optimal_objective_function_value);
  REQUIRE(optimal_objective_function_value ==
          avgas_optimal_objective_function_value);

  HighsInt from_row = 5;
  HighsInt to_row = 9;
  HighsInt set_num_row = to_row - from_row + 1;
  HighsInt get_num_row;
  std::vector<double> og_row56789_lower;
  std::vector<double> og_row56789_upper;
  std::vector<double> set_row56789_lower;
  std::vector<double> get_row56789_lower;
  og_row56789_lower.resize(lp.num_row_);
  og_row56789_upper.resize(lp.num_row_);
  set_row56789_lower.resize(set_num_row);
  get_row56789_lower.resize(lp.num_row_);
  set_row56789_lower[0] = 5.0;
  set_row56789_lower[1] = 6.0;
  set_row56789_lower[2] = 7.0;
  set_row56789_lower[3] = 8.0;
  set_row56789_lower[4] = 9.0;
  REQUIRE(highs.getRows(from_row, to_row, get_num_row, &og_row56789_lower[0],
                        &og_row56789_upper[0], get_num_nz, NULL, NULL,
                        NULL) == HighsStatus::kOk);
  REQUIRE(highs.changeRowsBounds(from_row, to_row, &set_row56789_lower[0],
                                 &og_row56789_upper[0]) == HighsStatus::kOk);
  REQUIRE(highs.getRows(from_row, to_row, get_num_row, &get_row56789_lower[0],
                        &og_row56789_upper[0], get_num_nz, NULL, NULL,
                        NULL) == HighsStatus::kOk);
  REQUIRE(get_num_row == set_num_row);
  for (HighsInt usr_row = 0; usr_row < get_num_row; usr_row++)
    REQUIRE(get_row56789_lower[usr_row] == set_row56789_lower[usr_row]);
  REQUIRE(highs.changeRowsBounds(from_row, to_row, &og_row56789_lower[0],
                                 &og_row56789_upper[0]) == HighsStatus::kOk);

  callRun(highs, options.log_options, "highs.run()", HighsStatus::kOk);

  highs.getInfoValue("objective_function_value",
                     optimal_objective_function_value);
  REQUIRE(optimal_objective_function_value ==
          avgas_optimal_objective_function_value);
}
TEST_CASE("LP-delete", "[highs_data]") {
  // Rather better testing of deleteCols() and deleteRows()
  Highs highs;
  HighsOptions options;
  HighsLogOptions& log_options = options.log_options;

  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
    options.output_flag = false;
  }

  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);

  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);

  const HighsLp& lp = highs.getLp();

  callRun(highs, log_options, "highs.run()", HighsStatus::kOk);

  double adlittle_objective_function_value;
  highs.getInfoValue("objective_function_value",
                     adlittle_objective_function_value);

  HighsRandom random(0);
  double objective_function_value;
  HighsInt num_nz = lp.a_matrix_.numNz();
  std::vector<HighsInt> mask;
  std::vector<HighsInt> mask_check;
  HighsInt get_num_nz;
  std::vector<HighsInt> get_start;
  std::vector<HighsInt> get_index;
  std::vector<double> get_cost;
  std::vector<double> get_lower;
  std::vector<double> get_upper;
  std::vector<double> get_value;

  // Test deleteCols
  HighsInt num_col = lp.num_col_;
  HighsInt rm_num_col = num_col / 5;
  assert(rm_num_col >= 10);
  mask.assign(num_col, 0);
  mask_check.assign(num_col, 0);
  HighsInt num_col_k = 0;
  for (;;) {
    HighsInt iCol = random.integer(num_col);
    if (mask[iCol]) continue;
    mask[iCol] = 1;
    num_col_k++;
    if (num_col_k >= rm_num_col) break;
  }
  HighsInt new_col_index = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    if (!mask[iCol]) {
      mask_check[iCol] = new_col_index;
      new_col_index++;
    } else {
      mask_check[iCol] = -1;
    }
  }
  HighsInt get_num_col;
  get_cost.resize(rm_num_col);
  get_lower.resize(rm_num_col);
  get_upper.resize(rm_num_col);
  get_start.resize(rm_num_col);
  get_index.resize(num_nz);
  get_value.resize(num_nz);

  // Get the set of cols to be removed - so that they can be reintroduced
  REQUIRE(highs.getCols(&mask[0], get_num_col, &get_cost[0], &get_lower[0],
                        &get_upper[0], get_num_nz, &get_start[0], &get_index[0],
                        &get_value[0]) == HighsStatus::kOk);
  REQUIRE(get_num_col == rm_num_col);
  get_index.resize(get_num_nz);
  get_value.resize(get_num_nz);

  // Remove the set of cols
  REQUIRE(highs.deleteCols(&mask[0]) == HighsStatus::kOk);
  REQUIRE(mask == mask_check);
  REQUIRE(lp.num_col_ == num_col - rm_num_col);

  // Replace the set of cols
  REQUIRE(highs.addCols(get_num_col, &get_cost[0], &get_lower[0], &get_upper[0],
                        get_num_nz, &get_start[0], &get_index[0],
                        &get_value[0]) == HighsStatus::kOk);
  REQUIRE(lp.num_col_ == num_col);

  callRun(highs, log_options, "highs.run()", HighsStatus::kOk);

  highs.getInfoValue("objective_function_value", objective_function_value);

  REQUIRE(
      std::fabs(objective_function_value - adlittle_objective_function_value) <
      double_equal_tolerance);

  // Test deleteRows
  HighsInt num_row = lp.num_row_;
  HighsInt rm_num_row = num_row / 5;
  assert(rm_num_row >= 10);
  mask.assign(num_row, 0);
  mask_check.assign(num_row, 0);
  HighsInt num_row_k = 0;
  for (;;) {
    HighsInt iRow = random.integer(num_row);
    if (mask[iRow]) continue;
    mask[iRow] = 1;
    num_row_k++;
    if (num_row_k >= rm_num_row) break;
  }
  HighsInt new_row_index = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    if (!mask[iRow]) {
      mask_check[iRow] = new_row_index;
      new_row_index++;
    } else {
      mask_check[iRow] = -1;
    }
  }
  HighsInt get_num_row;
  get_lower.resize(rm_num_row);
  get_upper.resize(rm_num_row);
  get_start.resize(rm_num_row);
  get_index.resize(num_nz);
  get_value.resize(num_nz);

  // Get the set of rows to be removed - so that they can be reintroduced
  REQUIRE(highs.getRows(&mask[0], get_num_row, &get_lower[0], &get_upper[0],
                        get_num_nz, &get_start[0], &get_index[0],
                        &get_value[0]) == HighsStatus::kOk);
  REQUIRE(get_num_row == rm_num_row);
  get_index.resize(get_num_nz);
  get_value.resize(get_num_nz);

  // Remove the set of rows
  REQUIRE(highs.deleteRows(&mask[0]) == HighsStatus::kOk);
  REQUIRE(mask == mask_check);
  REQUIRE(lp.num_row_ == num_row - rm_num_row);

  // Replace the set of rows
  REQUIRE(highs.addRows(get_num_row, &get_lower[0], &get_upper[0], get_num_nz,
                        &get_start[0], &get_index[0],
                        &get_value[0]) == HighsStatus::kOk);
  REQUIRE(lp.num_row_ == num_row);

  callRun(highs, log_options, "highs.run()", HighsStatus::kOk);

  highs.getInfoValue("objective_function_value", objective_function_value);

  REQUIRE(
      std::fabs(objective_function_value - adlittle_objective_function_value) <
      double_equal_tolerance);
}

void HighsStatusReport(const HighsLogOptions& log_options, std::string message,
                       HighsStatus status) {
  if (!dev_run) return;
  highsLogUser(log_options, HighsLogType::kInfo,
               "%s: HighsStatus = %" HIGHSINT_FORMAT " - %s\n", message.c_str(),
               (int)status, highsStatusToString(status).c_str());
}

void callRun(Highs& highs, const HighsLogOptions& log_options,
             std::string message, const HighsStatus require_return_status) {
  HighsStatus return_status = highs.run();
  HighsStatusReport(log_options, message, return_status);
  REQUIRE(return_status == require_return_status);
}

bool areLpColEqual(const HighsInt num_col0, const double* colCost0,
                   const double* colLower0, const double* colUpper0,
                   const HighsInt num_nz0, const HighsInt* Astart0,
                   const HighsInt* Aindex0, const double* Avalue0,
                   const HighsInt num_col1, const double* colCost1,
                   const double* colLower1, const double* colUpper1,
                   const HighsInt num_nz1, const HighsInt* Astart1,
                   const HighsInt* Aindex1, const double* Avalue1,
                   const double infinite_bound) {
  if (num_col0 != num_col1) {
    if (dev_run)
      printf("areLpColEqual: %" HIGHSINT_FORMAT
             " = num_col0 != num_col1 = %" HIGHSINT_FORMAT "\n",
             num_col0, num_col1);
    return false;
  }
  if (!num_col0) return true;
  HighsInt num_col = num_col0;
  for (HighsInt col = 0; col < num_col; col++) {
    if (colCost0[col] != colCost1[col]) {
      if (dev_run)
        printf("areLpColEqual: %g = colCost0[%" HIGHSINT_FORMAT
               "] != colCost1[%" HIGHSINT_FORMAT "] = %g\n",
               colCost0[col], col, col, colCost1[col]);
      return false;
    }
  }
  for (HighsInt col = 0; col < num_col; col++) {
    if (colLower0[col] <= -infinite_bound && colLower1[col] <= -infinite_bound)
      continue;
    if (colLower0[col] != colLower1[col]) {
      if (dev_run)
        printf("areLpColEqual: %g = colLower0[%" HIGHSINT_FORMAT
               "] != colLower1[%" HIGHSINT_FORMAT "] = %g\n",
               colLower0[col], col, col, colLower1[col]);
      return false;
    }
    if (colUpper0[col] >= infinite_bound && colUpper1[col] >= infinite_bound)
      continue;
    if (colUpper0[col] != colUpper1[col]) {
      if (dev_run)
        printf("areLpColEqual: %g = colUpper0[%" HIGHSINT_FORMAT
               "] != colUpper1[%" HIGHSINT_FORMAT "] = %g\n",
               colUpper0[col], col, col, colUpper1[col]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    if (dev_run)
      printf("areLpColEqual: %" HIGHSINT_FORMAT
             " = num_nz0 != num_nz1 = %" HIGHSINT_FORMAT "\n",
             num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (HighsInt col = 0; col < num_col; col++) {
    if (Astart0[col] != Astart1[col]) {
      if (dev_run)
        printf("areLpColEqual: %" HIGHSINT_FORMAT " = Astart0[%" HIGHSINT_FORMAT
               "] != Astart1[%" HIGHSINT_FORMAT "] = %" HIGHSINT_FORMAT "\n",
               Astart0[col], col, col, Astart1[col]);
      return false;
    }
  }
  HighsInt num_nz = num_nz0;
  for (HighsInt nz = 0; nz < num_nz; nz++) {
    if (Aindex0[nz] != Aindex1[nz]) {
      if (dev_run)
        printf("areLpColEqual: %" HIGHSINT_FORMAT " = Aindex0[%" HIGHSINT_FORMAT
               "] != Aindex1[%" HIGHSINT_FORMAT "] = %" HIGHSINT_FORMAT "\n",
               Aindex0[nz], nz, nz, Aindex1[nz]);
      return false;
    }
    if (Avalue0[nz] != Avalue1[nz]) {
      if (dev_run)
        printf("areLpColEqual: %g = Avalue0[%" HIGHSINT_FORMAT
               "] != Avalue1[%" HIGHSINT_FORMAT "] = %g\n",
               Avalue0[nz], nz, nz, Avalue1[nz]);
      return false;
    }
  }
  return true;
}

bool areLpRowEqual(const HighsInt num_row0, const double* rowLower0,
                   const double* rowUpper0, const HighsInt num_nz0,
                   const HighsInt* ARstart0, const HighsInt* ARindex0,
                   const double* ARvalue0, const HighsInt num_row1,
                   const double* rowLower1, const double* rowUpper1,
                   const HighsInt num_nz1, const HighsInt* ARstart1,
                   const HighsInt* ARindex1, const double* ARvalue1,
                   const double infinite_bound) {
  if (num_row0 != num_row1) {
    if (dev_run)
      printf("areLpRowEqual: %" HIGHSINT_FORMAT
             " = num_row0 != num_row1 = %" HIGHSINT_FORMAT "\n",
             num_row0, num_row1);
    return false;
  }
  if (!num_row0) return true;
  HighsInt num_row = num_row0;
  for (HighsInt row = 0; row < num_row; row++) {
    if (rowLower0[row] <= -infinite_bound && rowLower1[row] <= -infinite_bound)
      continue;
    if (rowLower0[row] != rowLower1[row]) {
      if (dev_run)
        printf("areLpRowEqual: %g = rowLower0[%" HIGHSINT_FORMAT
               "] != rowLower1[%" HIGHSINT_FORMAT "] = %g\n",
               rowLower0[row], row, row, rowLower1[row]);
      return false;
    }
    if (rowUpper0[row] >= infinite_bound && rowUpper1[row] >= infinite_bound)
      continue;
    if (rowUpper0[row] != rowUpper1[row]) {
      if (dev_run)
        printf("areLpRowEqual: %g = rowUpper0[%" HIGHSINT_FORMAT
               "] != rowUpper1[%" HIGHSINT_FORMAT "] = %g\n",
               rowUpper0[row], row, row, rowUpper1[row]);
      return false;
    }
  }
  if (num_nz0 != num_nz1) {
    if (dev_run)
      printf("areLpRowEqual: %" HIGHSINT_FORMAT
             " = num_nz0 != num_nz1 = %" HIGHSINT_FORMAT "\n",
             num_nz0, num_nz1);
    return false;
  }
  if (!num_nz0) return true;
  for (HighsInt row = 0; row < num_row; row++) {
    if (ARstart0[row] != ARstart1[row]) {
      if (dev_run)
        printf("areLpRowEqual: %" HIGHSINT_FORMAT
               " = ARstart0[%" HIGHSINT_FORMAT "] != ARstart1[%" HIGHSINT_FORMAT
               "] = %" HIGHSINT_FORMAT "\n",
               ARstart0[row], row, row, ARstart1[row]);
      return false;
    }
  }
  HighsInt num_nz = num_nz0;
  for (HighsInt nz = 0; nz < num_nz; nz++) {
    if (ARindex0[nz] != ARindex1[nz]) {
      if (dev_run)
        printf("areLpRowEqual: %" HIGHSINT_FORMAT
               " = ARindex0[%" HIGHSINT_FORMAT "] != ARindex1[%" HIGHSINT_FORMAT
               "] = %" HIGHSINT_FORMAT "\n",
               ARindex0[nz], nz, nz, ARindex1[nz]);
      return false;
    }
    if (ARvalue0[nz] != ARvalue1[nz]) {
      if (dev_run)
        printf("areLpRowEqual: %g = ARvalue0[%" HIGHSINT_FORMAT
               "] != ARvalue1[%" HIGHSINT_FORMAT "] = %g\n",
               ARvalue0[nz], nz, nz, ARvalue1[nz]);
      return false;
    }
  }
  return true;
}

bool areLpEqual(const HighsLp lp0, const HighsLp lp1,
                const double infinite_bound) {
  bool return_bool;
  if (lp0.a_matrix_.format_ != lp1.a_matrix_.format_) return false;
  if (lp0.num_col_ > 0 && lp1.num_col_ > 0) {
    HighsInt lp0_num_nz = lp0.a_matrix_.start_[lp0.num_col_];
    HighsInt lp1_num_nz = lp1.a_matrix_.start_[lp1.num_col_];
    return_bool = areLpColEqual(
        lp0.num_col_, &lp0.col_cost_[0], &lp0.col_lower_[0], &lp0.col_upper_[0],
        lp0_num_nz, &lp0.a_matrix_.start_[0], &lp0.a_matrix_.index_[0],
        &lp0.a_matrix_.value_[0], lp1.num_col_, &lp1.col_cost_[0],
        &lp1.col_lower_[0], &lp1.col_upper_[0], lp1_num_nz,
        &lp1.a_matrix_.start_[0], &lp1.a_matrix_.index_[0],
        &lp1.a_matrix_.value_[0], infinite_bound);
    if (!return_bool) return return_bool;
  }
  if (lp0.num_row_ > 0 && lp1.num_row_ > 0) {
    HighsInt lp0_num_nz = 0;
    HighsInt lp1_num_nz = 0;
    return_bool = areLpRowEqual(
        lp0.num_row_, &lp0.row_lower_[0], &lp0.row_upper_[0], lp0_num_nz, NULL,
        NULL, NULL, lp1.num_row_, &lp1.row_lower_[0], &lp1.row_upper_[0],
        lp1_num_nz, NULL, NULL, NULL, infinite_bound);
  }
  return return_bool;
}

void testDeleteKeep(const HighsIndexCollection& index_collection) {
  HighsInt delete_from_index;
  HighsInt delete_to_index;
  HighsInt keep_from_index;
  HighsInt keep_to_index;
  HighsInt current_set_entry;
  const std::vector<HighsInt>& set = index_collection.set_;
  const std::vector<HighsInt>& mask = index_collection.mask_;
  const HighsInt dimension = index_collection.dimension_;
  if (dev_run) {
    if (index_collection.is_interval_) {
      printf("With index interval [%" HIGHSINT_FORMAT ", %" HIGHSINT_FORMAT
             "] in [%d, %" HIGHSINT_FORMAT "]\n",
             index_collection.from_, index_collection.to_, 0, dimension - 1);
    } else if (index_collection.is_set_) {
      printf("With index set\n");
      for (HighsInt entry = 0; entry < index_collection.set_num_entries_;
           entry++)
        printf(" %2" HIGHSINT_FORMAT "", entry);
      printf("\n");
      for (HighsInt entry = 0; entry < index_collection.set_num_entries_;
           entry++)
        printf(" %2" HIGHSINT_FORMAT "", set[entry]);
      printf("\n");
    } else {
      printf("With index mask\n");
      for (HighsInt index = 0; index < dimension; index++)
        printf(" %2" HIGHSINT_FORMAT "", index);
      printf("\n");
      for (HighsInt index = 0; index < dimension; index++)
        printf(" %2" HIGHSINT_FORMAT "", mask[index]);
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
    for (HighsInt index = 0; index < dimension; index++) {
      if (mask[index]) {
        keep_to_index = index - 1;
        break;
      }
    }
  }
  if (dev_run)
    printf("Keep   [%2d, %2" HIGHSINT_FORMAT "]\n", 0, keep_to_index);
  if (keep_to_index >= dimension - 1) return;
  for (HighsInt k = 0; k < dimension; k++) {
    updateOutInIndex(index_collection, delete_from_index, delete_to_index,
                     keep_from_index, keep_to_index, current_set_entry);
    if (dev_run)
      printf("Delete [%2" HIGHSINT_FORMAT ", %2" HIGHSINT_FORMAT
             "]; keep [%2" HIGHSINT_FORMAT ", %2" HIGHSINT_FORMAT "]\n",
             delete_from_index, delete_to_index, keep_from_index,
             keep_to_index);
    if (delete_to_index >= dimension - 1 || keep_to_index >= dimension - 1)
      break;
  }
}

bool testAllDeleteKeep(HighsInt num_row) {
  // Test the extraction of intervals from index collections
  std::vector<HighsInt> set = {1, 4, 5, 8};
  std::vector<HighsInt> mask = {0, 1, 0, 0, 1, 1, 0, 0, 1, 0};

  HighsIndexCollection index_collection;
  index_collection.dimension_ = num_row;
  index_collection.is_interval_ = false;
  index_collection.from_ = 3;
  index_collection.to_ = 6;
  index_collection.is_set_ = false;
  index_collection.set_num_entries_ = 4;
  index_collection.set_ = set;
  index_collection.is_mask_ = false;
  index_collection.mask_ = mask;

  HighsInt save_from = index_collection.from_;
  HighsInt save_set_0 = set[0];
  HighsInt save_mask_0 = mask[0];

  HighsInt to_pass = 2;  // 2
  for (HighsInt pass = 0; pass <= to_pass; pass++) {
    if (dev_run)
      printf("\nTesting delete-keep: pass %" HIGHSINT_FORMAT "\n", pass);
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
  HighsInt log_dev_level;
  output_flag = dev_run;
  log_to_console = true;
  log_dev_level = kHighsLogDevLevelVerbose;
  log_options.output_flag = &output_flag;
  log_options.log_file_stream = NULL;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  highsLogDev(log_options, HighsLogType::kVerbose, "\nReporting LP: %s\n",
              message);
  reportLp(log_options, lp, HighsLogType::kVerbose);
}

void messageReportMatrix(const char* message, const HighsInt num_col,
                         const HighsInt num_nz, const HighsInt* start,
                         const HighsInt* index, const double* value) {
  HighsLogOptions log_options;
  bool output_flag = true;
  bool log_to_console = false;
  HighsInt log_dev_level = kHighsLogDevLevelInfo;
  log_options.log_file_stream = stdout;
  log_options.output_flag = &output_flag;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  highsLogDev(log_options, HighsLogType::kVerbose, "\nReporting Matrix: %s\n",
              message);
  reportMatrix(log_options, message, num_col, num_nz, start, index, value);
}
