#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
//#include "io/FilereaderLp.h"

const bool dev_run = false;
const double inf = kHighsInf;

TEST_CASE("lp-incompatible-bounds", "[iis]") {
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
  //  highs.setOptionValue("output_flag", dev_run);
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
  std::string model = "galenet";
  std::string model_file =
      std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs highs;
  //  highs.setOptionValue("output_flag", dev_run);

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
}
