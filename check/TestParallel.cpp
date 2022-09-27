#include "Highs.h"
#include "catch.hpp"
#include "parallel/HighsParallel.h"
#include "parallel/HighsRaceTimer.h"

using namespace highs;

const HighsInt numThreads = (std::thread::hardware_concurrency() + 1) / 2;
const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

bool equalObjective(const double obj0, const double obj1);

void noReturnDumbWork();
double dumbWork(const HighsInt n);

void noReturnSpawn();
void returnSpawn();
double concurrentLpSolve(const bool use_highs_run = true);

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
  total_time += concurrentLpSolve(true);
  //  total_time += concurrentLpSolve(false);
  printf("Total time is %g\n", total_time);
}

void noReturnSpawn() {
  // this worked when there was an overload of dumbWork(), as it
  // creates a lambda that calls the overload of dumbWork() without an
  // argument and overload resolution works normally
  // parallel::spawn([]() { dumbWork(); });

  // this works now that there is no overload "double dumbWork(const
  // HighsInt n);"
  //
  parallel::spawn(&noReturnDumbWork);
  // @Leona. It compiles (and runs OK), but with the warning, ah okay. Nevermind
  // the warning.
  // @Julian The warning goes away if you explicitly pass it as a function
  // pointer. Never tested it with function pointers since it does not make
  // sense to use it with just a function pointer.
  //         Having only a function pointer means there is no return and no
  //         arguments, so by definition nothing can happen. If you use it with
  //         a lambda expression that captures the arguments to a function and
  //         calls it, then this will not occur.
  //
  //   In file included from
  //   /home/jajhall/HiGHS/src/parallel/HighsSplitDeque.h:28,
  //                  from
  //                  /home/jajhall/HiGHS/src/parallel/HighsTaskExecutor.h:25,
  //                  from /home/jajhall/HiGHS/src/parallel/HighsMutex.h:18,
  //                  from /home/jajhall/HiGHS/src/parallel/HighsParallel.h:16,
  //                  from /home/jajhall/HiGHS/check/TestParallel.cpp:3:
  // /home/jajhall/HiGHS/src/parallel/HighsTask.h: In instantiation of ‘void
  // HighsTask::setTaskData(F&&) [with F = void (&)()]’:
  // /home/jajhall/HiGHS/src/parallel/HighsSplitDeque.h:294:5:   required from
  // ‘void HighsSplitDeque::push(F&&) [with F = void (&)()]’
  // /home/jajhall/HiGHS/src/parallel/HighsParallel.h:41:3:   required from
  // ‘void highs::parallel::spawn(HighsSplitDeque*, F&&) [with F = void (&)()]’
  // /home/jajhall/HiGHS/src/parallel/HighsParallel.h:46:8:   required from
  // ‘void highs::parallel::spawn(F&&) [with F = void (&)()]’
  // /home/jajhall/HiGHS/check/TestParallel.cpp:44:35:   required from here
  // /home/jajhall/HiGHS/src/parallel/HighsTask.h:99:19: warning: invalid
  // application of ‘sizeof’ to a function type [-Wpointer-arith]
  //    99 |     static_assert(sizeof(F) <= sizeof(taskData),
  //       |                   ^~~~~~~~~

  // Another version that works is the following.
  // This resolves the overload resolution to store the function pointer for the
  // void dumbWork() overload void (*f)() = dumbWork; with the function pointer
  // you can use it similarly to the version "parallel::spawn(dumbWork)"
  // parallel::spawn(f);

  noReturnDumbWork();
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

void noReturnDumbWork() {
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

double concurrentLpSolve(const bool use_highs_run) {
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

  const HighsInt primal_simplex_index = 0;
  const HighsInt dual_simplex_index = 1;
  const HighsInt ipm_index = 2;

  assert(lp_solvers[primal_simplex_index] == kLpSolverPrimalSimplex);
  assert(lp_solvers[dual_simplex_index] == kLpSolverDualSimplex);
  assert(lp_solvers[ipm_index] == kLpSolverInteriorPoint);

  // Set up a vector of pointers to Highs instances
  vector<std::unique_ptr<Highs>> parallel_highs;
  for (HighsInt lp_solver = 0; lp_solver < num_lp_solvers; lp_solver++)
    parallel_highs.emplace_back(new Highs());

  std::vector<double> run_time(num_lp_solvers);
  std::vector<double> objective_function_value(num_lp_solvers);
  std::vector<HighsStatus> run_status(num_lp_solvers);
  std::vector<HighsModelStatus> model_status(num_lp_solvers);

  // Force the primal simplex instance to use primal simplex
  parallel_highs[primal_simplex_index]->setOptionValue("simplex_strategy",
                                                       kSimplexStrategyDual);
  // Force the dual simplex instance to use dual simplex
  parallel_highs[dual_simplex_index]->setOptionValue("simplex_strategy",
                                                     kSimplexStrategyPrimal);
  // Force the IPM instance to use IPM
  parallel_highs[ipm_index]->setOptionValue("solver", kIpmString);

  run_status.assign(num_lp_solvers, HighsStatus::kError);
  // when you have this as a vector you can also use the parallel::for_each call
  parallel::for_each(0, num_lp_solvers, [&](HighsInt start, HighsInt end) {
    for (HighsInt lp_solver = start; lp_solver < end; ++lp_solver) {
      parallel_highs[lp_solver]->setOptionValue("output_flag", false);
      parallel_highs[lp_solver]->passModel(lp);
      // todo, set pointer to race timer
      run_status[lp_solver] = parallel_highs[lp_solver]->run();
      run_time[lp_solver] = parallel_highs[lp_solver]->getRunTime();
      objective_function_value[lp_solver] =
          parallel_highs[lp_solver]->getInfo().objective_function_value;
      parallel_highs[lp_solver]->setBasis();
      parallel_highs[lp_solver]->setSolution();

      // todo call race_timer decreaseLimit() function if solved successfully
    }
  });
  printf("Status values: primal = %s; dual = %s; ipm = %s\n",
         highsStatusToString(run_status[0]).c_str(),
         highsStatusToString(run_status[1]).c_str(),
         highsStatusToString(run_status[2]).c_str());

  printf("Solve times:  primal = %f6.4; dual = %f6.4; ipm = %f6.4\n",
         run_time[0], run_time[1], run_time[2]);

  // the additional for-loop inside the for_each call would not really be
  // necessary with grainSize = 1, which is the default value of the fourth
  // argument to for_each in that case end-start should always be equal to 1 and
  // you can remove the for-loop and change "[&](HighsInt start, HighsInt end){"
  // to "[&](HighsInt lp_solver, HighsInt){"

  HighsRaceTimer<double> race_timer;

  if (use_highs_run) {
    parallel::TaskGroup tg;
    for (HighsInt lp_solver = 0; lp_solver < num_lp_solvers; lp_solver++) {
      tg.spawn([&parallel_highs, &run_status, &model_status, &race_timer, lp_solver]() {
        // todo: somehow pass a pointer to the race_timer into the solver. The
        // solver should check if such a pointer was passed and call the
        // limitReached(currentTime) function to determine whether a limit was
        // reached. If limitReached returns true then the solver should stop.
        // e.g. something like:
        // if (race_timer_ptr && race_timer_ptr->limitReached(currentTime))
        //  modelstatus = HighsModelStatus::kIterationLimit;
        parallel_highs[lp_solver]->passRaceTimer(&race_timer);
        run_status[lp_solver] = parallel_highs[lp_solver]->run();
	model_status[lp_solver] = parallel_highs[lp_solver]->getModelStatus();
        // this should check for the status returned when the race timer limit
        // was reached and only call decrease limit if it was not reached
	printf("LP solver %d returns status %s and model status %s\n",
	       int(lp_solver),
	       highsStatusToString(run_status[lp_solver]).c_str(),
	       parallel_highs[lp_solver]->modelStatusToString(model_status[lp_solver]).c_str()); fflush(stdout);
        if (model_status[lp_solver] != HighsModelStatus::kRaceTimerStop)
          race_timer.decreaseLimit(parallel_highs[lp_solver]->getRunTime());
      });
    }
    tg.taskWait();
  } else {
    printf("Setting up parallel HighsLpSolverObject\n"); fflush(stdout);
    vector<std::unique_ptr<HighsLpSolverObject>> parallel_lp_solver_object;
    vector<std::unique_ptr<HighsBasis>> parallel_basis;
    vector<std::unique_ptr<HighsSolution>> parallel_solution;
    vector<std::unique_ptr<HighsInfo>> parallel_info;
    vector<std::unique_ptr<HEkk>> parallel_ekk_instance;
    vector<std::unique_ptr<HighsOptions>> parallel_options;
    vector<std::unique_ptr<HighsTimer>> parallel_timer;
    for (HighsInt lp_solver = 0; lp_solver < num_lp_solvers; lp_solver++) {
      parallel_basis.emplace_back(new HighsBasis());
      parallel_solution.emplace_back(new HighsSolution());
      parallel_info.emplace_back(new HighsInfo());
      parallel_ekk_instance.emplace_back(new HEkk());
      parallel_options.emplace_back(new HighsOptions());
      parallel_timer.emplace_back(new HighsTimer());
      parallel_lp_solver_object.emplace_back(new HighsLpSolverObject(lp,
								     *parallel_basis[lp_solver],
								     *parallel_solution[lp_solver],
								     *parallel_info[lp_solver],
								     *parallel_ekk_instance[lp_solver],
								     *parallel_options[lp_solver],
								     *parallel_timer[lp_solver]));
      parallel_ekk_instance[lp_solver]->race_timer_ = &race_timer;
    }
    printf("Set up parallel HighsLpSolverObject\n"); 
  }
  double total_time = 0;
  for (HighsInt lp_solver = 0; lp_solver < num_lp_solvers; lp_solver++) {
    objective_function_value[lp_solver] =
        parallel_highs[lp_solver]->getInfo().objective_function_value;
    run_time[lp_solver] = parallel_highs[lp_solver]->getRunTime();
    total_time += run_time[lp_solver];
  }

  printf("Status values: primal = %s; dual = %s; ipm = %s\n",
         highsStatusToString(run_status[primal_simplex_index]).c_str(),
         highsStatusToString(run_status[dual_simplex_index]).c_str(),
         highsStatusToString(run_status[ipm_index]).c_str());

  printf("Solve times:  primal = %f6.4; dual = %f6.4; ipm = %f6.4\n",
         run_time[0], run_time[1], run_time[2]);

  REQUIRE(equalObjective(objective_function_value[dual_simplex_index],
                         objective_function_value[primal_simplex_index]));
  REQUIRE(equalObjective(objective_function_value[dual_simplex_index],
                         objective_function_value[ipm_index]));

  return total_time;
}

bool equalObjective(const double obj0, const double obj1) {
  return std::fabs((obj0 - obj1) / std::max(1.0, fabs(obj0))) <
         double_equal_tolerance;
}
