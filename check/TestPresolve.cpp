#include <cstdio>

#include "Highs.h"
#include "catch.hpp"
#include "io/LoadProblem.h"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"

const double kOptimalQap04 = 32;

// No commas in test case name.
TEST_CASE("irash-qap04", "[highs_presolve]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/qap04.mps";

  Highs highs;
  HighsStatus highs_status = highs.readModel(filename);
  REQUIRE(highs_status==HighsStatus::OK);

  HighsOptions options;
  options.icrash = true;
  options.icrash_starting_weight = 10;
  options.icrash_approximate_minimization_iterations = 100;

  // highs.options_ is now private!
  // highs.options_ = options;
  highs_status = highs.passHighsOptions(options);
  REQUIRE(highs_status==HighsStatus::OK);

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
