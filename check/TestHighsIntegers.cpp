#include "HCheckConfig.h"
#include "catch.hpp"
#include "util/HighsCDouble.h"
#include "util/HighsIntegers.h"
#include "util/HighsRandom.h"

const bool dev_run = false;

TEST_CASE("HighsIntegers", "[util]") {
  double x1 = 6.4700675;
  double x2 = 0.27425;
  double x3 = 5.68625;

  // simple test case, that requires denominators above 1000, but not
  // well behaved ones as they contain multiple powers of two
  // this should get found by multiplying them by 600 before running
  // the continued fraction algorithm
  std::vector<double> tmp{x1, x2, x3};
  double integralscalar = HighsIntegers::integralScale(tmp, 1e-6, 1e-9);
  REQUIRE(integralscalar == 400000);

  if (dev_run) printf("integral scalar is %g\n", integralscalar);

  // stress test the algorithm with a constructed case:
  // add 6 fractions with prime number denominators just below 1000.
  // This will blow up the common denominator to a value around 9e17
  // with this magnitude the results are still representable
  // in an 64bit integers, but not in 53bits double precision.
  // The double precision error is already far above 1.0
  // so for computing the correct fraction it is necessary to use HighsCDouble
  // which represents a roughly quad precision number as the unevaluated sum of
  // two double precision numbers, otherwise the algorithm will fail.
  int64_t primes[] = {967, 971, 977, 983};
  std::vector<double> values{1.0 / primes[0], 2.0 / primes[1], 3.0 / primes[2],
                             4.0 / primes[3]};

  integralscalar = HighsIntegers::integralScale(values, 1e-6, 1e-6);

  REQUIRE(integralscalar == primes[0] * primes[1] * primes[2] * primes[3]);

  if (dev_run) printf("integral scalar is %g\n", integralscalar);
}

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
