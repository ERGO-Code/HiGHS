#include <cstdio>

#include "Highs.h"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"
#include "io/LoadProblem.h"
#include "TestInstanceLoad.h"
#include "catch.hpp"

#ifdef __linux__
#include <unistd.h>
#elif _WIN32
#define NOGDI
#include <windows.h>
#else

#endif

const double kOptimalQap04 = 32;

// No commas in test case name.
TEST_CASE("ff-qap04", "[highs_presolve]") {
  HighsOptions options;
  std::string dir = GetCurrentWorkingDir();

  std::cout << dir << std::endl;

  // For debugging use the latter.
  options.model_file = dir + "/../../check/instances/qap04.mps";
  //options.model_file = dir + "/check/instances/qap04.mps";

  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  REQUIRE(read_status == HighsStatus::OK);

  Highs highs;
  options.icrash = true;
  options.icrash_starting_weight = 10;

  highs.options_ = options;
  HighsStatus init_status = highs.initializeLp(lp);
  REQUIRE(init_status == HighsStatus::OK);

  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::OK);

  const HighsSolution& x = highs.getSolution();

  // assert objective is close enough to the optimal one
  // c'x
  double ctx = highs.getObjectiveValue();
  double difference = std::fabs(kOptimalQap04 - ctx);
  REQUIRE(difference < 0.18);

  // assert residual is below threshold
  // r
  std::vector<double> residual;
  HighsSolution solution = highs.getSolution();
  HighsStatus result = calculateResidual(highs.getLp(), solution, residual);
  REQUIRE(result == HighsStatus::OK);

  double r = getNorm2(residual);
  REQUIRE(r < 1e-08);
}