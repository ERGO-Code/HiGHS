#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include "Highs.h"

// test case failing
HighsStatus issue425() {
  HighsLp lp;
  lp.numCol_ = 4;
  lp.numRow_ = 4;

  lp.Astart_.push_back(0);
  lp.Astart_.push_back(3);
  lp.Astart_.push_back(5);
  lp.Astart_.push_back(6);
  lp.Astart_.push_back(7);

  lp.Aindex_.push_back(0);
  lp.Avalue_.push_back(1);
  lp.Aindex_.push_back(2);
  lp.Avalue_.push_back(1);
  lp.Aindex_.push_back(3);
  lp.Avalue_.push_back(1);

  lp.Aindex_.push_back(1);
  lp.Avalue_.push_back(2);
  lp.Aindex_.push_back(3);
  lp.Avalue_.push_back(1);

  lp.Aindex_.push_back(3);
  lp.Avalue_.push_back(1);

  lp.Aindex_.push_back(3);
  lp.Avalue_.push_back(1);

  lp.colLower_.assign(lp.numCol_, 0);
  lp.colUpper_.assign(lp.numCol_, HIGHS_CONST_INF);
  
  std::vector<double> b{1, 2, 2, 4};
  lp.rowLower_ = b;
  lp.rowUpper_ = b;

  lp.colCost_.push_back(1);
  lp.colCost_.push_back(1);
  lp.colCost_.push_back(1);
  lp.colCost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
   assert(status == HighsStatus::OK);

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