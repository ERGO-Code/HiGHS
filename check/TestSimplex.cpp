#include "Avgas.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/LoadProblem.h"
//#include "lp_data/HighsLpUtils.h"

void testSolver(Highs& highs, const std::string solver,
                const int simplex_strategy, int& simplex_iteration_count) {
  const int adlittle_default_simplex_iteration_count = 86;
  double default_time_limit;
  int default_simplex_iteration_limit;
  HighsModelStatus model_status;
  HighsStatus return_status;

  const HighsInfo& info = highs.getHighsInfo();

  return_status = highs.setHighsOptionValue("solver", solver);
  REQUIRE(return_status == HighsStatus::OK);

  return_status =
      highs.setHighsOptionValue("simplex_strategy", simplex_strategy);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.getHighsOptionValue("time_limit", default_time_limit);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.getHighsOptionValue("simplex_iteration_limit",
                                            default_simplex_iteration_limit);
  REQUIRE(return_status == HighsStatus::OK);

  // Vanilla solve: get solution time to calibrate time limit test
  double run_time = highs.getHighsRunTime();
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);
  const double single_solve_run_time = highs.getHighsRunTime() - run_time;

  simplex_iteration_count += adlittle_default_simplex_iteration_count;
  REQUIRE(info.simplex_iteration_count == simplex_iteration_count);

  // Only perform the time limit test if the solve time is large enough
  const double min_run_time_for_test = 0.001;
  if (single_solve_run_time > min_run_time_for_test) {
    const int ideal_num_solve = 10;
    const double local_time_limit = ideal_num_solve * single_solve_run_time;

    // Solve with time limit
    run_time = highs.getHighsRunTime();
    printf("Current run time is %g\n", run_time);

    double use_time_limit = run_time + local_time_limit;
    return_status = highs.setHighsOptionValue("time_limit", use_time_limit);
    REQUIRE(return_status == HighsStatus::OK);

    const int max_num_solve = 10 * ideal_num_solve;
    int num_solve;
    for (num_solve = 0; num_solve < max_num_solve; num_solve++) {
      return_status = highs.setBasis();
      return_status = highs.run();
      if (highs.getModelStatus() == HighsModelStatus::REACHED_TIME_LIMIT) break;
    }
    REQUIRE(num_solve < max_num_solve);
    run_time = highs.getHighsRunTime();
    printf("Current run time is %g: time limit is %g (difference = %g)\n",
           run_time, use_time_limit, run_time - use_time_limit);

    printf("Required %d solves (ideally %d - max %d)\n", num_solve,
           ideal_num_solve, max_num_solve);
  } else {
    printf(
        "Not performed the time limit test since solve time is %g <= %g = "
        "min_run_time_for_test\n",
        single_solve_run_time, min_run_time_for_test);
  }
  return_status = highs.setHighsOptionValue("time_limit", default_time_limit);
  REQUIRE(return_status == HighsStatus::OK);

  // Solve with iteration limit
  /*
  simplex_iteration_count = info.simplex_iteration_count;
  const int further_simplex_iterations = 10;
  int use_simplex_iteration_limit;
  // First of all check that no iterations are performed if the
  // iteration limit is the current number of iterations

  use_simplex_iteration_limit = simplex_iteration_count;

  return_status = highs.setHighsOptionValue("simplex_iteration_limit",
                                            use_simplex_iteration_limit);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.setBasis();
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  model_status = highs.getModelStatus();
  printf("Returns status = %d; model status = %s\n", (int)return_status,
         highs.highsModelStatusToString(model_status).c_str());
  REQUIRE(return_status == HighsStatus::Warning);
  REQUIRE(model_status == HighsModelStatus::REACHED_ITERATION_LIMIT);
  REQUIRE(simplex_iteration_count == info.simplex_iteration_count);
  simplex_iteration_count = info.simplex_iteration_count;

  // Now check that it stops after 10 iterations
  use_simplex_iteration_limit =
      simplex_iteration_count + further_simplex_iterations;
  return_status = highs.setHighsOptionValue("simplex_iteration_limit",
                                            use_simplex_iteration_limit);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.setBasis();
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::Warning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::REACHED_ITERATION_LIMIT);

  simplex_iteration_count = use_simplex_iteration_limit;
  REQUIRE(simplex_iteration_count == info.simplex_iteration_count);
  simplex_iteration_count = info.simplex_iteration_count;

  return_status = highs.setHighsOptionValue("simplex_iteration_limit",
                                            default_simplex_iteration_limit);
  REQUIRE(return_status == HighsStatus::OK);
  */
}

// No commas in test case name.
TEST_CASE("LP-simplex", "[highs_simplex]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string model = "adlittle";
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  int simplex_iteration_count = 0;

  HighsOptions options;
  options.model_file = filename;

  // Read mps.
  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  REQUIRE(read_status == HighsStatus::OK);

  Highs highs(options);
  HighsStatus return_status = highs.passModel(lp);
  REQUIRE(return_status == HighsStatus::OK);

  std::string solver = "simplex";
  int from_i = (int)SimplexStrategy::SIMPLEX_STRATEGY_MIN;
  int to_i = 1 + from_i;  // (int)SimplexStrategy::SIMPLEX_STRATEGY_MAX;
  for (int i = from_i; i < to_i; i++) {
    SimplexStrategy simplex_strategy = static_cast<SimplexStrategy>(i);
    printf("Simplex strategy %d\n", (int)simplex_strategy);
    testSolver(highs, solver, simplex_strategy, simplex_iteration_count);
  }
  /*
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  options.model_file = filename;
  read_status = loadLpFromFile(options, lp);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.passModel(lp);
  REQUIRE(return_status == HighsStatus::OK);

  */
  /*
  return_status = highs.getHighsInfoValue("simplex_iteration_count",
  simplex_iteration_count); REQUIRE(return_status == HighsStatus::OK);
  REQUIRE(simplex_iteration_count == 86);
  */
}
