#include "catch.hpp"
#include "parallel/HighsParallel.h"

using namespace highs;

const HighsInt numThreads = (std::thread::hardware_concurrency() + 1) / 2;
const bool dev_run = true;

void dumbWork();
double dumbWork(const HighsInt n);

void noReturnSpawn();
void returnSpawn();

TEST_CASE("test-parallel", "[highs_test_parallel]") {
  if (dev_run) {
    printf("std::thread::hardware_concurrency() returns %d\n",
           std::thread::hardware_concurrency());
    printf("Using %d threads\n", int(numThreads));
  }
  HighsTaskExecutor::shutdown();
  parallel::initialize_scheduler(numThreads);
  noReturnSpawn();
  returnSpawn();
}

void noReturnSpawn() {
  // Why is the following "invalid use of void expression"
  //  parallel::spawn(dumbWork());
  //  dumbWork();
  //  parallel::sync();
}

void returnSpawn() {
  double v10;
  parallel::spawn([&]() { v10 = dumbWork(10); });
  const double v11 = dumbWork(11);
  if (dev_run) printf("Before sync(): v10 = %g; v11 = %g\n", v10, v11);
  parallel::sync();
  if (dev_run) printf("After  sync(): v10 = %g; v11 = %g\n", v10, v11);
}

void dumbWork() {
  const HighsInt dumb_work_n = 10;
  double sum_products = 0;
  for (HighsInt i = 1; i <= dumb_work_n; i++)
    for (HighsInt j = 1; i <= dumb_work_n; i++) sum_products += i * j;
}

double dumbWork(const HighsInt n) {
  double sum_products = 0;
  for (HighsInt i = 1; i <= n; i++)
    for (HighsInt j = 1; i <= n; i++) sum_products += i * j;
  return sum_products;
}
