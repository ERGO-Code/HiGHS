#include "Highs.h"
#include "catch.hpp"

#include "parallel/HighsParallel.h"

const bool dev_run = false;
const HighsInt dumb_work_n = 10;

void dumb_work();
double dumb_work(const HighsInt n);

TEST_CASE("test-parallel-spawn", "[highs_test_parallel]") {
  spawn(dumb_work());
  dumb_work();
  sync()
}

void dumb_work() {
  double sum_products = 0;
  for (HighsInt i  = 1; i <= dumb_work_n; i++) 
    for (HighsInt j  = 1; i <= dumb_work_n; i++)
      sum_products += i*j;
}

double dumb_work(const HighsInt n) {
  double sum_products = 0;
  for (HighsInt i  = 1; i <= n; i++) 
    for (HighsInt j  = 1; i <= n; i++)
      sum_products += i*j;
  return sum_products
}

