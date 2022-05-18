#include <cstdio>

#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsLpUtils.h"
#include "util/HighsUtils.h"

const double kOptimalQap04 = 32;

// No commas in test case name.
TEST_CASE("icrash-qap04", "[highs_presolve]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/qap04.mps";

  Highs highs;
  HighsStatus highs_status = highs.readModel(filename);
  REQUIRE(highs_status == HighsStatus::kOk);

  HighsOptions options;

  options.icrash = true;
  options.icrash_starting_weight = 10;
  options.icrash_approx_iter = 100;
  options.simplex_strategy = kSimplexStrategyPrimal;

  highs_status = highs.passHighsOptions(options);
  REQUIRE(highs_status == HighsStatus::kOk);

  HighsStatus run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  ICrashInfo info = highs.getICrashInfo();

  // assert objective is close enough to the optimal one
  // c'x
  REQUIRE(std::fabs(info.final_lp_objective - kOptimalQap04) < 0.18);

  // assert residual is below threshold
  REQUIRE(info.final_residual_norm_2 >= -1e08);
  REQUIRE(info.final_residual_norm_2 < 1e08);
}
