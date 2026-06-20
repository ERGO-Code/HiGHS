#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

#include <numeric>

const bool dev_run = false;

TEST_CASE("test-col-stuffing", "[highs_test_presolve_rules]") {
  HighsLp lp;
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMaximize;
  lp.col_cost_ = {1.8, 0.9, 1};
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, 1);
  //  lp.integrality_.assign(lp.num_col_, HighsVarType::kInteger);
  lp.row_lower_ = {-kHighsInf};
  lp.row_upper_ = {4};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, lp.num_col_};
  lp.a_matrix_.index_.resize( lp.num_col_);
  std::iota(lp.a_matrix_.index_.begin(), lp.a_matrix_.index_.end(), 0);
  lp.a_matrix_.value_ = {3, 2, 2};

  Highs h;
  //  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  h.setOptionValue("presolve_rule_test", kPresolveRuleColStuffing);
  h.run();
  h.writeSolution("", 1);

  h.resetGlobalScheduler(true);
  
}
