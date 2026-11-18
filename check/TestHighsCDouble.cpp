#include "HCheckConfig.h"
#include "catch.hpp"
#include "util/HighsCDouble.h"
#include "util/HighsRandom.h"

void testCeil(HighsCDouble x) {
  double ceil_x;
  double double_x;
  ceil_x = double(ceil(x));
  double_x = double(x);
  REQUIRE(ceil_x >= double_x);
  REQUIRE(ceil(x) >= x);
}

void testFloor(HighsCDouble x) {
  double floor_x;
  double double_x;
  floor_x = double(floor(x));
  double_x = double(x);
  REQUIRE(floor_x <= double_x);
  REQUIRE(floor(x) <= x);
}

TEST_CASE("HighsCDouble-ceil", "[util]") {
  HighsCDouble x;
  x = -1e-34;
  testCeil(x);
  x = -1e-32;
  testCeil(x);
  x = -1e-30;
  testCeil(x);
  x = -1e-23;
  testCeil(x);
  x = -1e-12;
  testCeil(x);
  x = -1e-1;
  testCeil(x);
  x = -0.99;
  testCeil(x);

  x = 0.99;
  testCeil(x);
  x = 1e-1;
  testCeil(x);
  x = 1e-12;
  testCeil(x);
  // This and rest failed in #2041
  x = 1e-23;
  testCeil(x);
  x = 1e-30;
  testCeil(x);
  x = 1e-32;
  testCeil(x);
  x = 1e-34;
  testCeil(x);

  HighsRandom rand;
  for (HighsInt k = 0; k < 1000; k++) {
    double man = rand.fraction();
    HighsInt power = 2 - rand.integer(5);
    double exp = std::pow(10, power);
    x = man * exp;
    testCeil(x);
  }
}

TEST_CASE("HighsCDouble-division", "[util]") {
  // see issue #3604: the division did not normalize its operands, so that the
  // gap between the high and the low part grew with each division until all
  // precision was lost

  // returns a HighsCDouble with the given value (|value| < 1) whose high part
  // is zero, i.e. that is not normalized. this happens with the lazily
  // normalized addition and subtraction.
  auto unnormalized = [](double value) {
    HighsCDouble x = 1e16;
    // value is smaller than half an ulp of 1e16 and ends up in the low part
    x += value;
    x -= 1e16;
    return x;
  };

  // repeated division of an unnormalized number by a double
  const double value = 0.9167;
  HighsCDouble x = unnormalized(value);
  REQUIRE(double(x) == value);
  for (HighsInt k = 0; k < 60; k++) {
    x /= 1.0;
    REQUIRE(std::abs(double(x) - value) <= 1e-15);
  }

  // division of an unnormalized number by an unnormalized number
  const double dividend = 0.375;
  HighsCDouble y = unnormalized(dividend) / unnormalized(value);
  REQUIRE(std::abs(double(y) - dividend / value) <= 1e-15);
}

TEST_CASE("HighsCDouble-floor", "[util]") {
  HighsCDouble x;

  x = 1e-34;
  testFloor(x);
  x = 1e-32;
  testFloor(x);
  x = 1e-30;
  testFloor(x);
  x = 1e-23;
  testFloor(x);
  x = 1e-12;
  testFloor(x);
  x = 1e-1;
  testFloor(x);
  x = 0.99;
  testFloor(x);

  x = -0.99;
  testFloor(x);
  x = -1e-1;
  testFloor(x);
  x = -1e-12;
  testFloor(x);
  // This and rest failed in #2041
  x = -1e-23;
  testFloor(x);
  x = -1e-30;
  testFloor(x);
  x = -1e-32;
  testFloor(x);
  x = -1e-34;
  testFloor(x);

  HighsRandom rand;
  for (HighsInt k = 0; k < 1000; k++) {
    double man = rand.fraction();
    HighsInt power = 2 - rand.integer(5);
    double exp = std::pow(10, power);
    x = man * exp;
    testFloor(x);
  }
}
