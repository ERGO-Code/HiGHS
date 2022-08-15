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
  // this worked when there was an overload of dumbWork(), as it
  // creates a lambda that calls the overload of dumbWork() without an
  // argument and overload resolution works normally
  // parallel::spawn([]() { dumbWork(); });

  // this works now that there is no overload "double dumbWork(const
  // HighsInt n);"
  //
  parallel::spawn(&noReturnDumbWork);
  // @Leona. It compiles (and runs OK), but with the warning, ah okay. Nevermind the warning.
  // @Julian The warning goes away if you explicitly pass it as a function pointer. Never tested it with function pointers since it does not make sense to use it with just a function pointer.
  //         Having only a function pointer means there is no return and no arguments, so by definition nothing can happen. If you use it with a lambda expression that captures the arguments to a function and calls it, then this will not occur.
  //
  //   In file included from /home/jajhall/HiGHS/src/parallel/HighsSplitDeque.h:28,
  //                  from /home/jajhall/HiGHS/src/parallel/HighsTaskExecutor.h:25,
  //                  from /home/jajhall/HiGHS/src/parallel/HighsMutex.h:18,
  //                  from /home/jajhall/HiGHS/src/parallel/HighsParallel.h:16,
  //                  from /home/jajhall/HiGHS/check/TestParallel.cpp:3:
// /home/jajhall/HiGHS/src/parallel/HighsTask.h: In instantiation of ‘void HighsTask::setTaskData(F&&) [with F = void (&)()]’:
// /home/jajhall/HiGHS/src/parallel/HighsSplitDeque.h:294:5:   required from ‘void HighsSplitDeque::push(F&&) [with F = void (&)()]’
// /home/jajhall/HiGHS/src/parallel/HighsParallel.h:41:3:   required from ‘void highs::parallel::spawn(HighsSplitDeque*, F&&) [with F = void (&)()]’
// /home/jajhall/HiGHS/src/parallel/HighsParallel.h:46:8:   required from ‘void highs::parallel::spawn(F&&) [with F = void (&)()]’
// /home/jajhall/HiGHS/check/TestParallel.cpp:44:35:   required from here
// /home/jajhall/HiGHS/src/parallel/HighsTask.h:99:19: warning: invalid application of ‘sizeof’ to a function type [-Wpointer-arith]
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
  // @Leona this (in /*...*/) should work:
  //
  // @Leona this doesn't compile: it gives
  // @Julian, ah I see, then use the unique_ptr version. This requires the HiGHS class to get a move constructor and/or a copy constructor. Wrapping it in a unique_ptr works around that.
  // adding a move/copy constructor should not be too hard, but there are some members, e.g. the FILE* for input/output handling that might need special handling, so the unique_ptr version is simpler than adding a move contructor.
  // In file included from /usr/include/c++/9/vector:66,
  //                  from /home/jajhall/HiGHS/src/lp_data/HighsLpUtils.h:19,
  //                  from /home/jajhall/HiGHS/src/Highs.h:21,
  //                  from /home/jajhall/HiGHS/check/TestParallel.cpp:1:
// /usr/include/c++/9/bits/stl_uninitialized.h: In instantiation of ‘_ForwardIterator std::uninitialized_copy(_InputIterator, _In// putIterator, _ForwardIterator) [with _InputIterator = std::move_iterator<Highs*>; _ForwardIterator = Highs*]’:
// /usr/include/c++/9/bits/stl_uninitialized.h:307:37:   required from ‘_ForwardIterator std::__uninitialized_copy_a(_InputIterator, _InputIterator, _ForwardIterator, std::allocator<_Tp>&) [with _InputIterator = std::move_iterator<Highs*>; _ForwardIterator = Highs*; _Tp = Highs]’
// /usr/include/c++/9/bits/stl_uninitialized.h:329:2:   required from ‘_ForwardIterator std::__uninitialized_move_if_noexcept_a(_InputIterator, _InputIterator, _ForwardIterator, _Allocator&) [with _InputIterator = Highs*; _ForwardIterator = Highs*; _Allocator = std::allocator<Highs>]’
// /usr/include/c++/9/bits/vector.tcc:659:48:   required from ‘void std::vector<_Tp, _Alloc>::_M_default_append(std::vector<_Tp, _Alloc>::size_type) [with _Tp = Highs; _Alloc = std::allocator<Highs>; std::vector<_Tp, _Alloc>::size_type = long unsigned int]’
// /usr/include/c++/9/bits/stl_vector.h:937:4:   required from ‘void std::vector<_Tp, _Alloc>::resize(std::vector<_Tp, _Alloc>::size_type) [with _Tp = Highs; _Alloc = std::allocator<Highs>; std::vector<_Tp, _Alloc>::size_type = long unsigned int]’
// /home/jajhall/HiGHS/check/TestParallel.cpp:98:39:   required from here
// /usr/include/c++/9/bits/stl_uninitialized.h:127:72: error: static assertion failed: result type must be constructible from value type of input range
//   127 |       static_assert(is_constructible<_ValueType2, decltype(*__first)>::value,
//       |
  
  /*
  std::vector<Highs> parallel_highs;
  parallel_highs.resize(num_lp_solvers);
  for (HighsInt lp_solver = 0; lp_solver < num_lp_solvers; lp_solver++) {
    parallel_highs[lp_solver].passModel(lp);
  }
  */
  
  // Alternatively this should work too:
  //
  // @Leona: Yes, it does
  vector<std::unique_ptr<Highs>> parallel_highs;
  for (HighsInt lp_solver = 0; lp_solver< num_lp_solvers; lp_solver++) {
    parallel_highs.emplace_back(new Highs());
    parallel_highs[lp_solver]->passModel(lp);
  }

   std::vector<double> run_time(num_lp_solvers);
   std::vector<double> objective_function_value(num_lp_solvers);
   std::vector<HighsStatus> run_status(num_lp_solvers);
   run_status.assign(num_lp_solvers, HighsStatus::kError);
  // when you have this as a vector you can also use the parallel::for_each call
   parallel::for_each(0, num_lp_solvers, [&](HighsInt start, HighsInt end){
    for (HighsInt lp_solver = start; lp_solver < end; ++lp_solver)
    {
      parallel_highs[lp_solver]->setOptionValue("output_flag", false);
      parallel_highs[lp_solver]->passModel(lp);
      // todo, set pointer to race timer
      run_status[lp_solver] = parallel_highs[lp_solver]->run();
      run_time[lp_solver] = parallel_highs[lp_solver]->getRunTime();
      objective_function_value[lp_solver] = parallel_highs[lp_solver]->getInfo().objective_function_value;
      
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
