#include "Highs.h"
#include "catch.hpp"
//#include "parallel/HighsParallel.h"
#include "parallel/HighsRaceTimer.h"

const bool dev_run = true;

void testSolver(Highs& highs, const std::string solver,
                const HighsInt simplex_strategy, const double time_limit) {
  highs.setBasis();
  highs.setSolution();
  highs.setOptionValue("solver", solver);
  if (simplex_strategy == kSimplexStrategyDualPlain)
    highs.setOptionValue("parallel", kHighsOffString);
  if (simplex_strategy >= 0)
    highs.setOptionValue("simplex_strategy", simplex_strategy);
  double run_time = -highs.getRunTime();
  HighsStatus return_status = highs.run();
  run_time += highs.getRunTime();
  HighsModelStatus model_status = highs.getModelStatus();
  printf(
      "Running for %11.4gs (limit %11.4g) with \"%s\" (strategy %d) returns "
      "\"%s\" and model status \"%s\"\n",
      run_time, time_limit, solver.c_str(), int(simplex_strategy),
      highsStatusToString(return_status).c_str(),
      highs.modelStatusToString(model_status).c_str());
  REQUIRE(model_status == HighsModelStatus::kRaceTimerStop);
}

TEST_CASE("test-race-timer", "[highs_test_race_timer]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";

  highs.setOptionValue("presolve", kHighsOffString);
  highs.setOptionValue("solver", kIpmString);
  highs.readModel(model_file);
  double time_limit = -highs.getRunTime();
  highs.run();
  time_limit += highs.getRunTime();
  time_limit /= 10;
  HighsLp lp = highs.getLp();

  HighsRaceTimer<double> race_timer;
  highs.passRaceTimer(&race_timer);

  race_timer.decreaseLimit(time_limit);

  // IPM
  testSolver(highs, kIpmString, -1, time_limit);
  // Serial dual simplex
  testSolver(highs, kSimplexString, 1, time_limit);
  // Parallel dual simplex (SIP)
  testSolver(highs, kSimplexString, 2, time_limit);
  // Parallel dual simplex (PAMI)
  testSolver(highs, kSimplexString, 3, time_limit);
  // Serial primal simplex
  testSolver(highs, kSimplexString, 4, time_limit);
}
