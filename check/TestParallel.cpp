#include "catch.hpp"
#include "parallel/HighsParallel.h"

using namespace highs;

const HighsInt numThreads = (std::thread::hardware_concurrency() + 1) / 2;
const bool dev_run = true;
const HighsInt dumb_work_n = 10;

void dumb_work();
double dumb_work(const HighsInt n);

TEST_CASE("test-parallel-spawn", "[highs_test_parallel]") {
  if (dev_run)
    printf("std::thread::hardware_concurrency() returns %d\n",
           std::thread::hardware_concurrency());
  if (dev_run) printf("Using %d threads\n", int(numThreads));
  HighsTaskExecutor::shutdown();
  parallel::initialize_scheduler(numThreads);
  // Why is the following "invalid use of void expression"
  //  parallel::spawn(dumb_work());
  //  dumb_work();
  //  parallel::sync();
  double v10;
  parallel::spawn([&]() { v10 = dumb_work(10); });
  const double v11 = dumb_work(11);
  printf("Before sync(): v10 = %g; v11 = %g\n", v10, v11);
  parallel::sync();
  printf("After  sync(): v10 = %g; v11 = %g\n", v10, v11);
}

void dumb_work() {
  double sum_products = 0;
  for (HighsInt i = 1; i <= dumb_work_n; i++)
    for (HighsInt j = 1; i <= dumb_work_n; i++) sum_products += i * j;
}

double dumb_work(const HighsInt n) {
  double sum_products = 0;
  for (HighsInt i = 1; i <= n; i++)
    for (HighsInt j = 1; i <= n; i++) sum_products += i * j;
  return sum_products;
}
