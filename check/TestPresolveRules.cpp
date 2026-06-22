#include <numeric>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("test-col-stuffing", "[highs_test_presolve_rules]") {
  HighsLp lp;

  Highs h;
  //  h.setOptionValue("output_flag", dev_run);
  const HighsLp& presolved_lp = h.getPresolvedLp();
  h.setOptionValue("presolve_rule_test", kPresolveRuleColStuffing);
  const bool lp0 = false;
  const bool lp1 = false;
  const bool lp1a = true;
  const bool lp1b = false;

  if (lp0) {
    lp.num_col_ = 3;
    lp.num_row_ = 1;
    lp.sense_ = ObjSense::kMaximize;
    lp.col_cost_ = {1.8, 0.9, 1};
    lp.col_lower_.assign(lp.num_col_, 0);
    lp.col_upper_.assign(lp.num_col_, 1);
    lp.row_lower_ = {-kHighsInf};
    lp.row_upper_ = {4};
    lp.a_matrix_.format_ = MatrixFormat::kRowwise;
    lp.a_matrix_.start_ = {0, lp.num_col_};
    lp.a_matrix_.index_.resize(lp.num_col_);
    std::iota(lp.a_matrix_.index_.begin(), lp.a_matrix_.index_.end(), 0);
    lp.a_matrix_.value_ = {3, 2, 2};

    for (int k = 0; k < 2; k++) {
      printf("\n3-variable knapsack: %s\n", k == 0 ? "LP" : "IP");
      REQUIRE(h.passModel(lp) == HighsStatus::kOk);
      h.setOptionValue("presolve_rule_test", kPresolveRuleColStuffing);
      h.run();
      h.writeSolution("", 1);
      lp.integrality_.assign(lp.num_col_, HighsVarType::kInteger);
    }
    lp.clear();
  }

  auto presolveOffOn = [&](const std::string& message) {
    h.setOptionValue(kPresolveString, kHighsOffString);
    for (int k = 0; k < 2; k++) {
      printf("\n============\n%s: presolve = %s\n============\n\n",
             message.c_str(), k == 0 ? "off" : "on");
      REQUIRE(h.passModel(lp) == HighsStatus::kOk);
      h.run();
      h.writeSolution("", 1);
      if (k == 1) REQUIRE(h.getInfo().simplex_iteration_count == 0);
      h.setOptionValue(kPresolveString, kHighsOnString);
    }
  };

  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.sense_ = ObjSense::kMinimize;
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, 1);
  lp.row_lower_ = {2.0};
  lp.row_upper_ = {kHighsInf};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, lp.num_col_};
  lp.a_matrix_.index_.resize(lp.num_col_);
  std::iota(lp.a_matrix_.index_.begin(), lp.a_matrix_.index_.end(), 0);
  if (lp1) {
    lp.col_cost_.assign(lp.num_col_, 1);
    lp.a_matrix_.value_.assign(lp.num_col_, 1);
    presolveOffOn("Capturing neos-787933 issue");
  }
  if (lp1a) {
    lp.col_cost_ = {2, 1};
    lp.a_matrix_.value_.assign(lp.num_col_, 1);
    presolveOffOn("Variant A neos-787933 issue");
  }
  if (lp1b) {
    lp.col_cost_ = {-2, -1};
    lp.a_matrix_.value_.assign(lp.num_col_, 1);
    presolveOffOn("Variant B neos-787933 issue");
  }
  lp.clear();

  h.resetGlobalScheduler(true);
}
