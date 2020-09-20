#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest.h>
#include "Highs.h"

int factorial(int number) { return number <= 1 ? number : factorial(number - 1) * number; }

TEST_CASE("testing the factorial function") {
    CHECK(factorial(1) == 1);
    CHECK(factorial(2) == 2);
    CHECK(factorial(3) == 6);
    CHECK(factorial(10) == 3628800);
}

// New tests, for each of the special cases.
// - col sing doubleton equality 
// - col sing doubleton inequality 
// - test zero cost col sing

// test zero cost col sing
HighsStatus zeroCostColSing() {
  HighsLp lp;
  lp.numCol_ = 2;
  lp.numRow_ = 1;

  lp.Astart_.push_back(0);
  lp.Astart_.push_back(1);
  lp.Astart_.push_back(2);

  lp.Aindex_.push_back(0);
  lp.Avalue_.push_back(0.5);

  lp.Aindex_.push_back(0);
  lp.Avalue_.push_back(0.5);

  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);
  
  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);

  lp.rowLower_.push_back(0.1);
  lp.rowUpper_.push_back(0.9);

  lp.colCost_.push_back(0);
  lp.colCost_.push_back(1);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::OK);

  status = highs.run();
  return status;
}

// handled by doubleton equality
HighsStatus colSingDoubletonEquality() 
{
  HighsLp lp;
  lp.numCol_ = 4;
  lp.numRow_ = 2;

  lp.Astart_.push_back(0);
  lp.Astart_.push_back(2);
  lp.Astart_.push_back(3);
  lp.Astart_.push_back(4);
  lp.Astart_.push_back(5);

  lp.Aindex_.push_back(0);
  lp.Aindex_.push_back(1);
  lp.Aindex_.push_back(0);
  lp.Aindex_.push_back(1);
  lp.Aindex_.push_back(1);

  lp.Avalue_.push_back(0.5);
  lp.Avalue_.push_back(0.5);
  lp.Avalue_.push_back(1);
  lp.Avalue_.push_back(1);
  lp.Avalue_.push_back(1);

  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);
  
  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);

  lp.colLower_.push_back(0);
  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);
  lp.colUpper_.push_back(1);

  lp.rowLower_.push_back(1);
  lp.rowUpper_.push_back(1);

  lp.rowLower_.push_back(0);
  lp.rowUpper_.push_back(1);

  lp.colCost_.push_back(1);
  lp.colCost_.push_back(2);
  lp.colCost_.push_back(1);
  lp.colCost_.push_back(1);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::OK);

  status = highs.run();
  return status;
}

HighsStatus colSingDoubletonInequality() 
{
  HighsLp lp;
  lp.numCol_ = 4;
  lp.numRow_ = 2;

  lp.Astart_.push_back(0);
  lp.Astart_.push_back(2);
  lp.Astart_.push_back(3);
  lp.Astart_.push_back(4);
  lp.Astart_.push_back(5);

  lp.Aindex_.push_back(0);
  lp.Aindex_.push_back(1);
  lp.Aindex_.push_back(0);
  lp.Aindex_.push_back(1);
  lp.Aindex_.push_back(1);

  lp.Avalue_.push_back(0.5);
  lp.Avalue_.push_back(0.5);
  lp.Avalue_.push_back(1);
  lp.Avalue_.push_back(1);
  lp.Avalue_.push_back(1);

  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);
  
  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);

  lp.colLower_.push_back(0);
  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);
  lp.colUpper_.push_back(1);

  lp.rowLower_.push_back(0);
  lp.rowUpper_.push_back(1);

  lp.rowLower_.push_back(0);
  lp.rowUpper_.push_back(1);

  lp.colCost_.push_back(1);
  lp.colCost_.push_back(2);
  lp.colCost_.push_back(1);
  lp.colCost_.push_back(1);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::OK);

  status = highs.run();
  return status;
}

// handled by doubleton equality
HighsStatus twoColSingDoubletonEquality() 
{
  HighsLp lp;
  lp.numCol_ = 2;
  lp.numRow_ = 1;

  lp.Astart_.push_back(0);
  lp.Astart_.push_back(1);
  lp.Astart_.push_back(2);

  lp.Aindex_.push_back(0);
  lp.Aindex_.push_back(0);

  lp.Avalue_.push_back(1);
  lp.Avalue_.push_back(1);

  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);
  
  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);

  lp.rowLower_.push_back(1);
  lp.rowUpper_.push_back(1);

  lp.colCost_.push_back(1);
  lp.colCost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::OK);

  status = highs.run();
  return status;
}

// handled by special case.
HighsStatus twoColSingDoubletonInequality() 
{
  HighsLp lp;
  lp.numCol_ = 2;
  lp.numRow_ = 1;

  lp.Astart_.push_back(0);
  lp.Astart_.push_back(1);
  lp.Astart_.push_back(2);

  lp.Aindex_.push_back(0);
  lp.Aindex_.push_back(0);

  lp.Avalue_.push_back(1);
  lp.Avalue_.push_back(1);

  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);
  
  lp.colLower_.push_back(0);
  lp.colUpper_.push_back(1);

  lp.rowLower_.push_back(0);
  lp.rowUpper_.push_back(1);

  lp.colCost_.push_back(1);
  lp.colCost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::OK);

  highs.run();
  status = highs.run();
  return status;
}

// No commas in test case name.
TEST_CASE("zero-cost [presolve-col-sing]") {
  std::cout << "Presolve 1." << std::endl;
  HighsStatus status = zeroCostColSing();
  std::string str = HighsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("col-sing-doubleton-eq [presolve-col-sing]") {
  std::cout << "Presolve 2." << std::endl;
  HighsStatus status =  colSingDoubletonEquality();
  std::string str = HighsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("col-sing-doubleton-ineq [presolve-col-sing]") {
  std::cout << "Presolve 3." << std::endl;
  HighsStatus status =  colSingDoubletonInequality();
  std::string str = HighsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("two-col-sing-doubleton-eq [presolve-col-sing]") {
  std::cout << "Presolve 4." << std::endl;
  HighsStatus status =  twoColSingDoubletonEquality();
  std::string str = HighsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("two-col-sing-doubleton-ineq [presolve-col-sing]") {
  std::cout << "Presolve 5." << std::endl;
  HighsStatus status =  twoColSingDoubletonInequality();
  std::string str = HighsStatusToString(status);
  REQUIRE(str == "OK");
}