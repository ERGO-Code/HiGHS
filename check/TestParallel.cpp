#include "Highs.h"
#include "catch.hpp"
#include "parallel/HighsParallel.h"
#include "parallel/HighsRaceTimer.h"

using namespace highs;

const HighsInt numThreads = (std::thread::hardware_concurrency() + 1) / 2;
const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

bool equalObjective(const double obj0, const double obj1);

void dumbWork();
double dumbWork(const HighsInt n);

void noReturnSpawn();
void returnSpawn();
double concurrentLpSolve(const bool use_race_timer = true);

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
  double total_time = 0;
  total_time += concurrentLpSolve(false);
  printf("Total time is %g\n", total_time);
}

void noReturnSpawn() {
  // @Leona Why is the following "invalid use of void expression"
  //
  // this works, as it creates a lambda that calls the overload of dumbWork()
  // without an argument and overload resolution works normally
  parallel::spawn([]() { dumbWork(); });

  // this would work if there was no overload "double dumbWork(const HighsInt
  // n);"
  // parallel::spawn(dumbWork);

  // Another version that works is the following.
  // This resolves the overload resolution to store the function pointer for the
  // void dumbWork() overload void (*f)() = dumbWork; with the function pointer
  // you can use it similarly to the version "parallel::spawn(dumbWork)"
  // parallel::spawn(f);

  dumbWork();
  parallel::sync();
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

double concurrentLpSolve(const bool use_race_timer) {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  std::string model_file;
  model_file = std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";

  highs.readModel(model_file);
  HighsLp lp = highs.getLp();
  // Run parallel primal simplex, dual simplex and IPM
  vector<HighsInt> lp_solvers;
  lp_solvers.push_back(kLpSolverPrimalSimplex);
  lp_solvers.push_back(kLpSolverDualSimplex);
  lp_solvers.push_back(kLpSolverInteriorPoint);
  HighsInt num_lp_solvers = lp_solvers.size();
  // @Leona How can I form a vector of Highs instances?
  //
  // this should work:
  vector<Highs> parallel_highs;
  parallel_highs.resize(num_lp_solvers);
  for (HighsInt lp_solver = 0; lp_solver < num_lp_solvers; lp_solver++) {
    parallel_highs[lp_solver].passModel(lp);
  }

  // Alternatively this should work to:
  // vector<std::unique_ptr<Highs>> parallel_highs;
  // for (HighsInt lp_solver = 0; lp_solver< num_lp_solvers; lp_solver++) {
  //  parallel_highs.emplace_back(new Highs());
  //  parallel_highs[lp_solver]->passModel(lp);
  //}

  // when you have this as a vector you can also use the parallel::for_each call
  // parallel::for_each(0, num_lp_solvers, [&](HighsInt start, HighsInt end){
  //  for(HighsInt lp_solver = start; lp_solver < end; ++lp_solver)
  //  {
  //    parallel_highs[lp_solver].setOptionValue("output_flag", false);
  //    parallel_highs[lp_solver].passModel(lp);
  //    // todo, set pointer to race timer, store return status
  //    parallel_highs[lp_solver].run();
  //
  //    // todo call race_timer decreaseLimit() function if solved successfully
  //  }
  //});

  // the additional for-loop inside the for_each call would not really be
  // necessary with grainSize = 1, which is the default value of the fourth
  // argument to for_each in that case end-start should always be equal to 1 and
  // you can remove the for-loop and change "[&](HighsInt start, HighsInt end){"
  // to "[&](HighsInt lp_solver, HighsInt){"

  //
  // vector<Highs*> parallel_highs;
  // for (HighsInt lp_solver = 0; lp_solver< num_lp_solvers; lp_solver++) {
  // parallel_highs.push_back(new Highs());
  // parallel_highs[lp_solver]->passModel(lp);
  // }
  // parallel_highs[kLpSolverPrimalSimplex]->setOptionValue("simplex_strategy",
  // kSimplexStrategyPrimal);

  Highs highs_primal_simplex;
  Highs highs_dual_simplex;
  Highs highs_ipm;
  highs_primal_simplex.setOptionValue("output_flag", false);
  highs_dual_simplex.setOptionValue("output_flag", false);
  highs_ipm.setOptionValue("output_flag", false);

  HighsStatus primal_simplex_status;
  HighsStatus dual_simplex_status;
  HighsStatus ipm_status;
  double primal_simplex_time;
  double dual_simplex_time;
  double ipm_time;
  double primal_simplex_objective;
  double dual_simplex_objective;
  double ipm_objective;
  //  parallel::spawn([&]() { primal_simplex_status =
  //  highs_primal_simplex.passModel(lp); }); parallel::sync();

  primal_simplex_status = highs_primal_simplex.passModel(lp);
  REQUIRE(primal_simplex_status == HighsStatus::kOk);

  dual_simplex_status = highs_dual_simplex.passModel(lp);
  REQUIRE(dual_simplex_status == HighsStatus::kOk);

  ipm_status = highs_ipm.passModel(lp);
  REQUIRE(ipm_status == HighsStatus::kOk);

  highs_primal_simplex.setOptionValue("simplex_strategy",
                                      kSimplexStrategyPrimal);
  highs_dual_simplex.setOptionValue("simplex_strategy", kSimplexStrategyDual);
  highs_ipm.setOptionValue("solver", kIpmString);

  primal_simplex_status = HighsStatus::kError;
  dual_simplex_status = HighsStatus::kError;
  ipm_status = HighsStatus::kError;

  if (use_race_timer) {
    // @Leona I realise that I need to specify a type - double in this
    // case - but I don't know how, and the following gives "error:
    // missing template arguments before ‘race_timer’"
    HighsRaceTimer<double> race_timer;

    parallel::TaskGroup tg;
    tg.spawn([&]() {
      // todo: somehow pass a pointer to the race_timer into the solver. The
      // solver should check if such a pointer was passed and call the
      // limitReached(currentTime) function to determine whether a limit was
      // reached. If limitReached returns true then the solver should stop.
      // e.g. something like:
      // if (race_timer_ptr && race_timer_ptr->limitReached(currentTime))
      //  modelstatus = HighsModelStatus::kIterationLimit;

      primal_simplex_status = highs_primal_simplex.run();
      // this should check for the status returned when the race timer limit was
      // reached and only call decrease limit if it was not reached
      if (highs_dual_simplex.getModelStatus() !=
          HighsModelStatus::kIterationLimit)
        race_timer.decreaseLimit(highs_primal_simplex.getHighsRunTime());
    });
    tg.spawn([&]() {
      dual_simplex_status = highs_dual_simplex.run();
      if (highs_dual_simplex.getModelStatus() !=
          HighsModelStatus::kIterationLimit)
        race_timer.decreaseLimit(highs_primal_simplex.getHighsRunTime());
    });
    tg.spawn([&]() {
      ipm_status = highs_ipm.run();
      // this should check for the status returned when the race timer limit was
      // reached and only call decrease limit if it was not reached
      if (highs_dual_simplex.getModelStatus() !=
          HighsModelStatus::kIterationLimit)
        race_timer.decreaseLimit(highs_primal_simplex.getHighsRunTime());
    });

    tg.taskWait();
  } else {
    parallel::TaskGroup tg;
    tg.spawn([&]() { primal_simplex_status = highs_primal_simplex.run(); });
    tg.spawn([&]() { dual_simplex_status = highs_dual_simplex.run(); });
    tg.spawn([&]() { ipm_status = highs_ipm.run(); });

    tg.sync();
    tg.sync();
    tg.sync();

    // alternatively  instead of the 3 sync() calls you can call: tg.taskWait();
  }

  primal_simplex_time = highs_primal_simplex.getRunTime();
  dual_simplex_time = highs_dual_simplex.getRunTime();
  ipm_time = highs_ipm.getRunTime();

  primal_simplex_objective =
      highs_primal_simplex.getInfo().objective_function_value;
  dual_simplex_objective =
      highs_dual_simplex.getInfo().objective_function_value;
  ipm_objective = highs_ipm.getInfo().objective_function_value;

  printf("Status values: primal = %s; dual = %s; ipm = %s\n",
         highsStatusToString(primal_simplex_status).c_str(),
         highsStatusToString(dual_simplex_status).c_str(),
         highsStatusToString(ipm_status).c_str());

  printf("Solve times:  primal = %f6.4; dual = %f6.4; ipm = %f6.4\n",
         primal_simplex_time, dual_simplex_time, ipm_time);

  REQUIRE(equalObjective(dual_simplex_objective, primal_simplex_objective));
  REQUIRE(equalObjective(dual_simplex_objective, ipm_objective));

  return primal_simplex_time + dual_simplex_time + ipm_time;
}

bool equalObjective(const double obj0, const double obj1) {
  return std::fabs((obj0 - obj1) / std::max(1.0, fabs(obj0))) <
         double_equal_tolerance;
}
