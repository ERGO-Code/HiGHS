#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
//#include "io/FilereaderLp.h"

const bool dev_run = true;
const double inf = kHighsInf;

void testIis(std::string& model,
	     HighsInt& num_iis_col, HighsInt& num_iis_row,
	     HighsInt* iis_col_index,
	     HighsInt* iis_row_index) {
  std::string model_file =
    std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs highs;
  highs.setOptionValue("output_flag", false);
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  HighsLp lp = highs.getLp();
  // Zero the objective
  lp.col_cost_.assign(lp.num_col_, 0);
  REQUIRE(highs.changeColsCost(0, lp.num_col_-1, lp.col_cost_.data()) == HighsStatus::kOk);

  for (HighsInt iisRow = 0; iisRow < num_iis_row; iisRow++) {
    HighsInt iRow = iis_row_index[iisRow];
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
    REQUIRE(highs.changeRowBounds(iRow, -kHighsInf, kHighsInf) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
    REQUIRE(highs.changeRowBounds(iRow, lower, upper) == HighsStatus::kOk);
    if (dev_run)
      printf("Removing IIS Row %d (LP row %d) yields optimality\n", int(iisRow), int(iRow));
  }
  for (HighsInt iisCol = 0; iisCol < num_iis_col; iisCol++) {
    HighsInt iCol = iis_col_index[iisCol];
    const double lower = lp.col_lower_[iCol];
    const double upper = lp.col_upper_[iCol];
    REQUIRE(highs.changeColBounds(iCol, -kHighsInf, kHighsInf) == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
    REQUIRE(highs.changeColBounds(iCol, lower, upper) == HighsStatus::kOk);
    if (dev_run)
      printf("Removing IIS Col %d (LP row %d) yields optimality\n", int(iisCol), int(iCol));
  }
}

void testMps(std::string& model) {
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);

  HighsInt num_iis_col;
  HighsInt num_iis_row;

  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row) == HighsStatus::kOk);
  REQUIRE(num_iis_col > 0);
  REQUIRE(num_iis_row > 0);
  std::vector<HighsInt> iis_col_index(num_iis_col);
  std::vector<HighsInt> iis_row_index(num_iis_row);
  std::vector<HighsInt> iis_col_bound(num_iis_col);
  std::vector<HighsInt> iis_row_bound(num_iis_row);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row, iis_col_index.data(),
                       iis_row_index.data(),
                       // iis_col_bound.data(), iis_row_bound.data()
                       nullptr, nullptr) == HighsStatus::kOk);
  if (dev_run)
    printf("Model %s has IIS with %d columns and %d rows\n", model.c_str(), int(num_iis_col), int(num_iis_row));
  testIis(model, num_iis_col, num_iis_row, iis_col_index.data(), iis_row_index.data());
}

TEST_CASE("lp-incompatible-bounds", "[iis]") {
  // LP has row0 and col2 with inconsistent bounds.
  //
  // When prioritising rows, row0 and its constituent columns (1, 2) should be found
  //
  // When prioritising columns, col2 and its constituent rows (0, 1) should be found
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 2;
  lp.col_cost_ = {0, 0, 0};
  lp.col_lower_ = {0, 0, 0};
  lp.col_upper_ = {1, 1, -1};
  lp.row_lower_ = {1, 0};
  lp.row_upper_ = {0, 1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 4};
  lp.a_matrix_.index_ = {1, 2, 0, 2};
  lp.a_matrix_.value_ = {1, 1, 1, 1};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  HighsInt num_iis_col;
  HighsInt num_iis_row;
  std::vector<HighsInt> iis_col_index;
  std::vector<HighsInt> iis_row_index;
  std::vector<HighsInt> iis_col_bound;
  std::vector<HighsInt> iis_row_bound;
  highs.setOptionValue("iis_strategy", kIisStrategyFromLpRowPriority);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row) == HighsStatus::kOk);
  REQUIRE(num_iis_col == 2);
  REQUIRE(num_iis_row == 1);
  iis_col_index.resize(num_iis_col);
  iis_row_index.resize(num_iis_row);
  iis_col_bound.resize(num_iis_col);
  iis_row_bound.resize(num_iis_row);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row, iis_col_index.data(),
                       iis_row_index.data(),
                       // iis_col_bound.data(), iis_row_bound.data()
                       nullptr, nullptr) == HighsStatus::kOk);
  REQUIRE(iis_col_index[0] == 1);
  REQUIRE(iis_col_index[1] == 2);
  REQUIRE(iis_row_index[0] == 0);
  highs.setOptionValue("iis_strategy", kIisStrategyFromLpColPriority);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row) == HighsStatus::kOk);
  REQUIRE(num_iis_col == 1);
  REQUIRE(num_iis_row == 2);
  iis_col_index.resize(num_iis_col);
  iis_row_index.resize(num_iis_row);
  iis_col_bound.resize(num_iis_col);
  iis_row_bound.resize(num_iis_row);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row, iis_col_index.data(),
                       iis_row_index.data(),
                       // iis_col_bound.data(), iis_row_bound.data()
                       nullptr, nullptr) == HighsStatus::kOk);
  REQUIRE(iis_col_index[0] == 2);
  REQUIRE(iis_row_index[0] == 0);
  REQUIRE(iis_row_index[1] == 1);
}

TEST_CASE("lp-get-iis", "[iis]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf, -inf, -inf};
  lp.row_upper_ = {8, 9, -2};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2, 4, 6};
  lp.a_matrix_.index_ = {0, 1, 0, 1, 0, 1};
  lp.a_matrix_.value_ = {2, 1, 1, 3, 1, 1};
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.passModel(lp);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  HighsInt num_iis_col;
  HighsInt num_iis_row;
  std::vector<HighsInt> iis_col_index;
  std::vector<HighsInt> iis_row_index;
  std::vector<HighsInt> iis_col_bound;
  std::vector<HighsInt> iis_row_bound;
  highs.setOptionValue("iis_strategy", kIisStrategyFromLpRowPriority);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row) == HighsStatus::kOk);
  REQUIRE(num_iis_col == 2);
  REQUIRE(num_iis_row == 1);
  iis_col_index.resize(num_iis_col);
  iis_row_index.resize(num_iis_row);
  iis_col_bound.resize(num_iis_col);
  iis_row_bound.resize(num_iis_row);
  REQUIRE(highs.getIis(num_iis_col, num_iis_row, iis_col_index.data(),
                       iis_row_index.data(),
                       // iis_col_bound.data(), iis_row_bound.data()
                       nullptr, nullptr) == HighsStatus::kOk);
  REQUIRE(iis_col_index[0] == 0);
  REQUIRE(iis_col_index[1] == 1);
  REQUIRE(iis_row_index[0] == 2);

}

TEST_CASE("lp-get-iis-galenet", "[iis]") {
  std::string model = "galenet";
  testMps(model);
}
