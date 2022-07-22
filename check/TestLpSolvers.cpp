#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

struct IterationCount {
  HighsInt simplex;
  HighsInt ipm;
  HighsInt crossover;
};

void testSolver(Highs& highs, const std::string solver,
                IterationCount& default_iteration_count,
                const HighsInt int_simplex_strategy = 0) {
  double default_time_limit;
  HighsInt default_simplex_iteration_limit;
  HighsInt default_ipm_iteration_limit;
  HighsModelStatus model_status;
  HighsStatus return_status;
  const bool perform_timeout_test = false;  // true;  //
  bool use_simplex = solver == "simplex";
  const HighsInfo& info = highs.getInfo();

  if (!dev_run) highs.setOptionValue("output_flag", false);
  return_status = highs.setOptionValue("solver", solver);
  REQUIRE(return_status == HighsStatus::kOk);

  if (use_simplex) {
    SimplexStrategy simplex_strategy =
        static_cast<SimplexStrategy>(int_simplex_strategy);
    if (simplex_strategy == SimplexStrategy::kSimplexStrategyDualTasks) return;
    if (dev_run)
      printf("Simplex strategy %" HIGHSINT_FORMAT "\n", int_simplex_strategy);
    return_status = highs.setOptionValue("simplex_strategy", simplex_strategy);
    REQUIRE(return_status == HighsStatus::kOk);
  }

  return_status = highs.getOptionValue("time_limit", default_time_limit);
  REQUIRE(return_status == HighsStatus::kOk);

  if (use_simplex) {
    return_status = highs.getOptionValue("simplex_iteration_limit",
                                         default_simplex_iteration_limit);
    REQUIRE(return_status == HighsStatus::kOk);
    // Clear the solver information - necessary if this is the second
    // or subsequent call to testSolver
    return_status = highs.clearSolver();
    REQUIRE(return_status == HighsStatus::kOk);
  } else {
    return_status = highs.getOptionValue("ipm_iteration_limit",
                                         default_ipm_iteration_limit);
    REQUIRE(return_status == HighsStatus::kOk);
  }

  // Vanilla solve: get solution time to calibrate time limit test
  double run_time = highs.getRunTime();
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);
  const double single_solve_run_time = highs.getRunTime() - run_time;

  if (use_simplex) {
    REQUIRE(info.simplex_iteration_count == default_iteration_count.simplex);
  } else {
    if (dev_run)
      printf("IPM: %" HIGHSINT_FORMAT "; Crossover: %" HIGHSINT_FORMAT "\n",
             info.ipm_iteration_count, info.crossover_iteration_count);
    REQUIRE(info.ipm_iteration_count == default_iteration_count.ipm);
    REQUIRE(info.crossover_iteration_count ==
            default_iteration_count.crossover);
  }

  // Only perform the time limit test if the solve time is large enough
  const double min_run_time_for_test = 0.001;
  if (perform_timeout_test && single_solve_run_time > min_run_time_for_test) {
    const HighsInt ideal_num_solve = 10;
    const double local_time_limit = ideal_num_solve * single_solve_run_time;

    // Solve with time limit
    run_time = highs.getRunTime();
    if (dev_run) printf("Current run time is %g\n", run_time);

    double use_time_limit = run_time + local_time_limit;
    return_status = highs.setOptionValue("time_limit", use_time_limit);
    REQUIRE(return_status == HighsStatus::kOk);

    const HighsInt max_num_solve = 10 * ideal_num_solve;
    HighsInt num_solve;
    for (num_solve = 0; num_solve < max_num_solve; num_solve++) {
      if (use_simplex) return_status = highs.setBasis();
      return_status = highs.run();
      if (highs.getModelStatus() == HighsModelStatus::kTimeLimit) break;
    }
    REQUIRE(num_solve < max_num_solve);
    run_time = highs.getRunTime();
    if (dev_run)
      printf("Current run time is %g: time limit is %g (difference = %g)\n",
             run_time, use_time_limit, run_time - use_time_limit);

    if (dev_run)
      printf("Required %" HIGHSINT_FORMAT " solves (ideally %" HIGHSINT_FORMAT
             " - max %" HIGHSINT_FORMAT ")\n",
             num_solve, ideal_num_solve, max_num_solve);
  } else {
    if (dev_run)
      printf(
          "Not performed the time limit test since solve time is %g <= %g = "
          "min_run_time_for_test\n",
          single_solve_run_time, min_run_time_for_test);
  }
  return_status = highs.setOptionValue("time_limit", default_time_limit);
  REQUIRE(return_status == HighsStatus::kOk);
  if (!use_simplex) {
    if (dev_run)
      printf("IPM: %" HIGHSINT_FORMAT "; Crossover: %" HIGHSINT_FORMAT "\n",
             info.ipm_iteration_count, info.crossover_iteration_count);
  }
  // Solve with iteration limit
  // First of all check that no iterations are performed if the
  // iteration limit is zero
  if (use_simplex) {
    return_status = highs.setOptionValue("simplex_iteration_limit", 0);
    REQUIRE(return_status == HighsStatus::kOk);

    return_status = highs.setBasis();
    REQUIRE(return_status == HighsStatus::kOk);
  } else {
    return_status = highs.setOptionValue("ipm_iteration_limit", 0);
    REQUIRE(return_status == HighsStatus::kOk);
  }

  return_status = highs.run();
  model_status = highs.getModelStatus();
  if (dev_run)
    printf("Returns status = %" HIGHSINT_FORMAT "; model status = %s\n",
           (HighsInt)return_status,
           highs.modelStatusToString(model_status).c_str());
  REQUIRE(return_status == HighsStatus::kWarning);
  REQUIRE(model_status == HighsModelStatus::kIterationLimit);

  if (use_simplex) {
    REQUIRE(info.simplex_iteration_count == 0);
  } else {
    REQUIRE(info.ipm_iteration_count == 0);
  }

  // Now check that simplex/IPM stops after 10/5 iterations
  const HighsInt further_simplex_iterations = 10;
  const HighsInt further_ipm_iterations = 5;
  if (use_simplex) {
    if (dev_run)
      printf("Setting simplex_iteration_limit = %" HIGHSINT_FORMAT "\n",
             further_simplex_iterations);
    return_status = highs.setOptionValue("simplex_iteration_limit",
                                         further_simplex_iterations);
    REQUIRE(return_status == HighsStatus::kOk);
    return_status = highs.clearSolver();
    REQUIRE(return_status == HighsStatus::kOk);
  } else {
    if (dev_run)
      printf("Setting ipm_iteration_limit = %" HIGHSINT_FORMAT "\n",
             further_ipm_iterations);
    return_status =
        highs.setOptionValue("ipm_iteration_limit", further_ipm_iterations);
    REQUIRE(return_status == HighsStatus::kOk);
  }

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kIterationLimit);

  if (use_simplex) {
    REQUIRE(info.simplex_iteration_count == further_simplex_iterations);
    return_status = highs.setOptionValue("simplex_iteration_limit",
                                         default_simplex_iteration_limit);
    REQUIRE(return_status == HighsStatus::kOk);
  } else {
    REQUIRE(info.ipm_iteration_count == further_ipm_iterations);
    return_status = highs.setOptionValue("ipm_iteration_limit",
                                         default_ipm_iteration_limit);
    REQUIRE(return_status == HighsStatus::kOk);
  }
}

void testSolversSetup(const std::string model,
                      IterationCount& model_iteration_count,
                      vector<HighsInt>& simplex_strategy_iteration_count) {
  if (model.compare("adlittle") == 0) {
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::kSimplexStrategyChoose] = 87;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::kSimplexStrategyDualPlain] = 87;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::kSimplexStrategyDualTasks] = 72;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::kSimplexStrategyDualMulti] = 73;
    simplex_strategy_iteration_count[(
        int)SimplexStrategy::kSimplexStrategyPrimal] = 94;
    model_iteration_count.ipm = 13;
    model_iteration_count.crossover = 2;
  }
}

void testSolvers(Highs& highs, IterationCount& model_iteration_count,
                 const vector<HighsInt>& simplex_strategy_iteration_count) {
  bool have_omp = true;

  /*
  HighsInt i = (HighsInt)SimplexStrategy::kSimplexStrategyPrimal;
  model_iteration_count.simplex = simplex_strategy_iteration_count[i];
  testSolver(highs, "simplex", model_iteration_count, i);
  */

  HighsInt from_i = (HighsInt)SimplexStrategy::kSimplexStrategyMin;
  HighsInt to_i =
      (HighsInt)SimplexStrategy::kSimplexStrategyDualMulti;  // PRIMAL;  // NUM;
  for (HighsInt i = from_i; i < to_i; i++) {
    if (!have_omp) {
      if (i == (HighsInt)SimplexStrategy::kSimplexStrategyDualTasks) continue;
      if (i == (HighsInt)SimplexStrategy::kSimplexStrategyDualMulti) continue;
    }
    model_iteration_count.simplex = simplex_strategy_iteration_count[i];
    testSolver(highs, "simplex", model_iteration_count, i);
  }
  // Only use IPX with 32-bit arithmetic
  // ToDo This is no longer true
#ifndef HIGHSINT64
  testSolver(highs, "ipm", model_iteration_count);
#endif
}

// No commas in test case name.
TEST_CASE("LP-solver", "[highs_lp_solver]") {
  std::string model;
  std::string model_file;
  IterationCount model_iteration_count;
  vector<HighsInt> simplex_strategy_iteration_count;
  simplex_strategy_iteration_count.resize(
      (HighsInt)SimplexStrategy::kSimplexStrategyNum);

  HighsLp lp;
  //  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  // Read mps
  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  testSolversSetup(model, model_iteration_count,
                   simplex_strategy_iteration_count);

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::kOk);

  return_status = highs.setBasis();
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  testSolvers(highs, model_iteration_count, simplex_strategy_iteration_count);

  // Now check that we can change model within the same Highs instance
  // First reset all the options to their default values
  return_status = highs.resetOptions();
  REQUIRE(return_status == HighsStatus::kOk);

  if (!dev_run) highs.setOptionValue("output_flag", false);

  model_file = std::string(HIGHS_DIR) + "/check/instances/etamacro.mps";
  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::kOk);

  return_status = highs.setBasis();
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  const HighsInfo& info = highs.getInfo();
  REQUIRE(info.num_dual_infeasibilities == 0);

  REQUIRE(info.simplex_iteration_count == 476);  // 444);

  HighsModelStatus model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  // Test the solver without scaling
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("simplex_scale_strategy", 0) ==
          HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  REQUIRE(info.simplex_iteration_count == 621);  // 584);  //
}

TEST_CASE("dual-objective-upper-bound", "[highs_lp_solver]") {
  std::string filename;
  HighsStatus status;
  HighsModelStatus model_status;
  bool bool_status;
  const double min_objective_function_value = -11.6389290663705;
  const double max_objective_function_value = 111.650960689315;
  const double smaller_min_objective_bound = -110.0;
  const double larger_min_objective_bound = -45.876;
  const double use_max_objective_bound = 150.0;
  double save_objective_bound;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();

  //  status = highs.setOptionValue("log_dev_level",
  //  kHighsLogDevLevelVerbose);

  double error;
  filename = std::string(HIGHS_DIR) + "/check/instances/e226.mps";
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  // Solve vanilla
  if (dev_run) printf("\nSolving vanilla LP\n");
  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  error = fabs((info.objective_function_value - min_objective_function_value) /
               min_objective_function_value);
  if (dev_run) printf("\nOptimal objective value error = %g\n", error);
  REQUIRE(error < 1e-14);

  // Set dual objective value upper bound after saving the default value
  status = highs.getOptionValue("objective_bound", save_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.setOptionValue("objective_bound", larger_min_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  // Solve again
  if (dev_run)
    printf(
        "\nSolving LP with presolve and dual objective value upper bound of "
        "%g\n",
        larger_min_objective_bound);
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  // Switch off presolve
  status = highs.setOptionValue("presolve", "off");
  REQUIRE(status == HighsStatus::kOk);

  // Solve again
  // This larger dual objective value upper bound is satisfied during phase 2
  if (dev_run)
    printf(
        "\nSolving LP without presolve and larger dual objective value upper "
        "bound of %g\n",
        larger_min_objective_bound);
  status = highs.clearSolver();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kObjectiveBound);

  // Solve again
  // This smaller dual objective value upper bound is satisfied at the start of
  // phase 2
  if (dev_run)
    printf(
        "\nSolving LP without presolve and smaller dual objective value upper "
        "bound of %g\n",
        smaller_min_objective_bound);
  status = highs.setOptionValue("objective_bound", smaller_min_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.setBasis();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kObjectiveBound);

  // Solve as maximization and ensure that the dual objective value upper bound
  // isn't used
  bool_status =
      highs.changeObjectiveSense(ObjSense::kMaximize) == HighsStatus::kOk;
  REQUIRE(bool_status);

  status = highs.setOptionValue("objective_bound", use_max_objective_bound);
  REQUIRE(status == HighsStatus::kOk);

  // Solve again
  if (dev_run)
    printf(
        "\nSolving LP as maximization without presolve and dual objective "
        "value "
        "upper bound of %g\n",
        use_max_objective_bound);
  status = highs.setBasis();
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  model_status = highs.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::kOptimal);

  error = fabs((info.objective_function_value - max_objective_function_value) /
               max_objective_function_value);
  if (dev_run) printf("\nOptimal objective value error = %g\n", error);
  REQUIRE(error < 1e-10);
}
