#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include "Highs.h"

// test case failing
HighsStatus issue425() {
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 4;

  lp.a_start_.push_back(0);
  lp.a_start_.push_back(3);
  lp.a_start_.push_back(5);
  lp.a_start_.push_back(6);
  lp.a_start_.push_back(7);

  lp.a_index_.push_back(0);
  lp.a_value_.push_back(1);
  lp.a_index_.push_back(2);
  lp.a_value_.push_back(1);
  lp.a_index_.push_back(3);
  lp.a_value_.push_back(1);

  lp.a_index_.push_back(1);
  lp.a_value_.push_back(2);
  lp.a_index_.push_back(3);
  lp.a_value_.push_back(1);

  lp.a_index_.push_back(3);
  lp.a_value_.push_back(1);

  lp.a_index_.push_back(3);
  lp.a_value_.push_back(1);

  lp.colLower_.assign(lp.num_col_, 0);
  lp.colUpper_.assign(lp.num_col_, kHighsInf);
  
  std::vector<double> b{1, 2, 2, 4};
  lp.rowLower_ = b;
  lp.rowUpper_ = b;

  lp.colCost_.push_back(1);
  lp.colCost_.push_back(1);
  lp.colCost_.push_back(1);
  lp.colCost_.push_back(2);

  lp.format_ = MatrixFormat::kColwise;
  
  Highs highs;
  HighsStatus status = highs.passModel(lp);
   assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}


TEST_CASE("presolve-issue-425") {
  std::cout << std::endl;
  std::cout << "Presolve issue 425." << std::endl;
  HighsStatus status = issue425();
  std::string str = HighsStatusToString(status);
  REQUIRE(str == "OK");
}
