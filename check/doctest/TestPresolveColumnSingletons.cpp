#include <doctest.h>
#include "Highs.h"

HighsInt factorial(HighsInt number) { return number <= 1 ? number : factorial(number - 1) * number; }

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
  lp.num_col_ = 2;
  lp.num_row_ = 1;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(1);
  lp.a_matrix_.start_.push_back(2);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.value_.push_back(0.5);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.value_.push_back(0.5);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(0.1);
  lp.row_upper_.push_back(0.9);

  lp.col_cost_.push_back(0);
  lp.col_cost_.push_back(1);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

// handled by doubleton equality
HighsStatus colSingDoubletonEquality() 
{
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 2;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(2);
  lp.a_matrix_.start_.push_back(3);
  lp.a_matrix_.start_.push_back(4);
  lp.a_matrix_.start_.push_back(5);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(1);

  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(1);
  lp.row_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);
  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(1);

  lp.format_ = MatrixFormat::kColwise;
  
  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

HighsStatus colSingDoubletonInequality() 
{
  HighsLp lp;
  lp.num_col_ = 4;
  lp.num_row_ = 2;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(2);
  lp.a_matrix_.start_.push_back(3);
  lp.a_matrix_.start_.push_back(4);
  lp.a_matrix_.start_.push_back(5);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(1);
  lp.a_matrix_.index_.push_back(1);

  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(0.5);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);
  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(1);

  lp.format_ = MatrixFormat::kColwise;
  
  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

// handled by doubleton equality
HighsStatus twoColSingDoubletonEquality() 
{
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(1);
  lp.a_matrix_.start_.push_back(2);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(0);

  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(1);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  status = highs.run();
  return status;
}

// handled by special case.
HighsStatus twoColSingDoubletonInequality() 
{
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 1;

  lp.a_matrix_.start_.push_back(0);
  lp.a_matrix_.start_.push_back(1);
  lp.a_matrix_.start_.push_back(2);

  lp.a_matrix_.index_.push_back(0);
  lp.a_matrix_.index_.push_back(0);

  lp.a_matrix_.value_.push_back(1);
  lp.a_matrix_.value_.push_back(1);

  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);
  
  lp.col_lower_.push_back(0);
  lp.col_upper_.push_back(1);

  lp.row_lower_.push_back(0);
  lp.row_upper_.push_back(1);

  lp.col_cost_.push_back(1);
  lp.col_cost_.push_back(2);

  Highs highs;
  HighsStatus status = highs.passModel(lp);
  assert(status == HighsStatus::kOk);

  highs.run();
  status = highs.run();
  return status;
}

// No commas in test case name.
TEST_CASE("zero-cost [presolve-col-sing]") {
  std::cout << "Presolve 1." << std::endl;
  HighsStatus status = zeroCostColSing();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("col-sing-doubleton-eq [presolve-col-sing]") {
  std::cout << "Presolve 2." << std::endl;
  HighsStatus status =  colSingDoubletonEquality();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("col-sing-doubleton-ineq [presolve-col-sing]") {
  std::cout << "Presolve 3." << std::endl;
  HighsStatus status =  colSingDoubletonInequality();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("two-col-sing-doubleton-eq [presolve-col-sing]") {
  std::cout << "Presolve 4." << std::endl;
  HighsStatus status =  twoColSingDoubletonEquality();
  std::string str = highsStatusToString(status);
  CHECK(str == "OK");
}

TEST_CASE("two-col-sing-doubleton-ineq [presolve-col-sing]") {
  std::cout << "Presolve 5." << std::endl;
  HighsStatus status =  twoColSingDoubletonInequality();
  std::string str = highsStatusToString(status);
  REQUIRE(str == "OK");
}

