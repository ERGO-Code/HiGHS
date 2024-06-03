#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
//#include "io/FilereaderLp.h"

const bool dev_run = false;
const double inf = kHighsInf;

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
