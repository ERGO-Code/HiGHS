#define CATCH_CONFIG_MAIN

#include "catch.hpp"

// Compile with
// g++ -std=c++11 -Wall -I$.../HiGHS/src/external/catch.hpp -c 000-CatchMain.cpp

bool isPrime(int n) {
 return (n==2 || n==3);
}

TEST_CASE("Correct status print" , "[highs_data]") {
  REQUIRE(isPrime(2) == true);
  REQUIRE(isPrime(2) == false);
}

TEST_CASE("Correct status print", "[highs_data]") {
  std::string str = HighsInputStatusToString(HighsInputStatus::OK);
  REQUIRE(str == "OK");
}