
#include "catch.hpp"
#include "util/HighsIntegers.h"

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
