#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const double inf = kHighsInf;
const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

HighsLp baseLp();

TEST_CASE("semi-continuous", "[highs_test_semi_variables]") {
  Highs highs;
  HighsStatus return_status;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  HighsModel model;
  HighsLp& lp = model.lp_;
  lp = baseLp();
  const HighsVarType continuous = HighsVarType::kContinuous;
  const HighsVarType semi_continuous = HighsVarType::kSemiContinuous;
  lp.integrality_ = {continuous, continuous, semi_continuous, continuous};
  const double save_upper2 = lp.col_upper_[2];
  lp.col_upper_[2] = inf;
  return_status = highs.passModel(model);
  REQUIRE(return_status == HighsStatus::kError);
  lp.col_upper_[2] = save_upper2;

  return_status = highs.passModel(model);
  REQUIRE(return_status == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
}

HighsLp baseLp() {
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;
  lp.col_cost_ = {1, 2, -4, -3};
  lp.col_lower_ = {0, 0, 1.1, 0};
  lp.col_upper_ = {inf, inf, 10, inf};
  lp.row_lower_ = {-inf, 0, 0, 0.5};
  lp.row_upper_ = {5, inf, inf, inf};
  lp.a_start_ = {0, 3, 6, 7, 8};
  lp.a_index_ = {0, 1,  2, 0,  1, 2, 3, 3};
  lp.a_value_ = {1, 2, -1, 1, -1, 3, 1, 1};
  lp.format_ = MatrixFormat::kColwise;
  lp.sense_ = ObjSense::kMaximize;
  return lp;
}
