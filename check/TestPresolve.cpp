#include <cstdio>

#include "Highs.h"
#include "TestInstanceLoad.h"
#include "catch.hpp"
#include "io/LoadProblem.h"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"

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
  // options.model_file = dir + "/check/instances/qap04.mps";

  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  REQUIRE(read_status == HighsStatus::OK);

  Highs highs;
  options.icrash = true;
  options.icrash_starting_weight = 10;
  options.icrash_approximate_minimization_iterations = 100;

  highs.options_ = options;
  HighsStatus init_status = highs.initializeLp(lp);
  REQUIRE(init_status == HighsStatus::OK);

  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::OK);

  ICrashInfo info = highs.getICrashInfo();

  // assert objective is close enough to the optimal one
  // c'x
  REQUIRE(std::fabs(info.final_lp_objective - kOptimalQap04) < 0.18);

  // assert residual is below threshold
  REQUIRE(info.final_residual_norm_2 >= -1e08);
  REQUIRE(info.final_residual_norm_2 < 1e08);
}