/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2022 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Leona Gottwald and Michael    */
/*    Feldmeier                                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsSolve.cpp
 * @brief Class-independent utilities for HiGHS
 */

#include "Highs.h"
#include "ipm/IpxWrapper.h"
#include "lp_data/HighsSolutionDebug.h"
#include "parallel/HighsParallel.h"
#include "parallel/HighsRaceTimer.h"
#include "simplex/HApp.h"

using namespace highs;

bool boundTermination(const HighsModelStatus& model_status) {
   switch (model_status) {
   case HighsModelStatus::kObjectiveBound:
   case HighsModelStatus::kObjectiveTarget:
     return true;
     break;
   case HighsModelStatus::kLoadError:
   case HighsModelStatus::kModelError:
   case HighsModelStatus::kModelEmpty:
   case HighsModelStatus::kOptimal:
   case HighsModelStatus::kInfeasible:
   case HighsModelStatus::kUnboundedOrInfeasible:
   case HighsModelStatus::kUnbounded:
   case HighsModelStatus::kNotset:
   case HighsModelStatus::kPresolveError:
   case HighsModelStatus::kSolveError:
   case HighsModelStatus::kPostsolveError:
   case HighsModelStatus::kTimeLimit:
   case HighsModelStatus::kIterationLimit:
   case HighsModelStatus::kUnknown:
   case HighsModelStatus::kInterrupted:
   case HighsModelStatus::kRaceTimerStop:
     return false;
     break;
   default:
     // All cases should have been considered so assert on reaching here
     assert(1 == 0);
   }
}

bool positiveModelStatus(const HighsModelStatus& model_status) {
   switch (model_status) {
   case HighsModelStatus::kLoadError:
   case HighsModelStatus::kModelError:
   case HighsModelStatus::kModelEmpty:
   case HighsModelStatus::kOptimal:
   case HighsModelStatus::kInfeasible:
   case HighsModelStatus::kUnboundedOrInfeasible:
   case HighsModelStatus::kUnbounded:
   case HighsModelStatus::kObjectiveBound:
   case HighsModelStatus::kObjectiveTarget:
     return true;
     break;
   case HighsModelStatus::kNotset:
   case HighsModelStatus::kPresolveError:
   case HighsModelStatus::kSolveError:
   case HighsModelStatus::kPostsolveError:
   case HighsModelStatus::kTimeLimit:
   case HighsModelStatus::kIterationLimit:
   case HighsModelStatus::kUnknown:
   case HighsModelStatus::kInterrupted:
   case HighsModelStatus::kRaceTimerStop:
     return false;
     break;
   default:
     // All cases should have been considered so assert on reaching here
     assert(1 == 0);
   }
}

void logLocalSolverOutcome(const HighsStatus return_status,
			     const HighsLpSolverObject& solver_object) {
  highsLogUser(solver_object.options_.log_options, HighsLogType::kInfo,
               "Local   solver %2d (time = %11.4g) returns status %s and model status %s\n",
	 int(solver_object.spawn_id_),
	 solver_object.run_time_,
	 highsStatusToString(return_status).c_str(),
	 utilModelStatusToString(solver_object.model_status_).c_str());
}
HighsStatus solveLpReturn(const HighsStatus return_status,
                          HighsLpSolverObject& solver_object,
                          const string message) {
  if (return_status != HighsStatus::kError) {
    // Analyse the HiGHS (basic) solution
    if (debugHighsLpSolution(message, solver_object) ==
        HighsDebugStatus::kLogicalError)
      return HighsStatus::kError;
  }
  solver_object.run_time_ = solver_object.timer_.readRunHighsClock() - solver_object.run_time_;
  logLocalSolverOutcome(return_status, solver_object);
  return return_status;
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus solveLp(HighsLpSolverObject& solver_object, const string message) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  HighsOptions& options = solver_object.options_;
  // Reset unscaled model status and solution params - except for
  // iteration counts
  resetModelStatusAndHighsInfo(solver_object);
  highsLogUser(options.log_options, HighsLogType::kInfo,
               (message + "\n").c_str());
  if (options.highs_debug_level > kHighsDebugLevelMin) {
    // Shouldn't have to check validity of the LP since this is done when it is
    // loaded or modified
    call_status = assessLp(solver_object.lp_, options);
    // If any errors have been found or normalisation carried out,
    // call_status will be ERROR or WARNING, so only valid return is OK.
    assert(call_status == HighsStatus::kOk);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "assessLp");
    if (return_status == HighsStatus::kError) return return_status;
  }
  // Make sure that the solver option is OK
  assert(solverOptionOk(options.solver));
  if (!solver_object.lp_.num_row_) {
    // Unconstrained LP so solve directly
    call_status = solveUnconstrainedLp(solver_object);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "solveUnconstrainedLp");
    return solveLpReturn(return_status, solver_object, message);
  }
  // LP has positive number of constraints
  //
  // Determine whether to use concurrent LP solvers
  //
  HighsInt max_concurrent_solvers = 1;
  const HighsInt max_threads = highs::parallel::num_threads();
  if (!solver_object.basis_.valid && !solver_object.solution_.value_valid) {
    // No hot start, so possibly use concurrent LP solvers
    max_concurrent_solvers = max_threads;
    // Prevent concurrent LP solvers unless parallel option is on
    if (options.parallel != kHighsOnString) max_concurrent_solvers = 1;
    // Prevent concurrent LP solvers if IPM is to be used
    if (options.solver == kIpmString) max_concurrent_solvers = 1;
    // Prevent concurrent LP solvers if the simplex strategy is
    // anything other than choose
    if (options.simplex_strategy != kSimplexStrategyChoose)
      max_concurrent_solvers = 1;
  }
  // Start with (at most) serial dual, serial primal and parallel dual
  // running concurrently
  const HighsInt max_spawnable_solvers = 3;
  const HighsInt spawnable_threads = max_threads - 1;
  std::string use_solver = options.solver;
  HighsInt spawnable_solvers = max_concurrent_solvers == 1 ? 0 : max_spawnable_solvers;
  // Save copies of simplex_strategy and simplex_max_concurrency
  const HighsInt simplex_strategy = options.simplex_strategy;
  const HighsInt simplex_max_concurrency = options.simplex_max_concurrency;

  assert(max_concurrent_solvers >= 1);
  // Ensure that concurrent LP solvers aren't used if IPM is forced
  assert(max_concurrent_solvers == 1 || use_solver != kIpmString);

  if (max_concurrent_solvers > 1) {
    // If using concurrent LP solvers, determine which solver is to be
    // used by this Highs instance, and how many to be used
    if (use_solver == kHighsChooseString) {
      // Can use interior point for LP, so this Highs instance uses
      // IPM
      use_solver = kIpmString;
    } else {
      assert(use_solver == kSimplexString);
      // Can only use simplex solvers, so this Highs instance uses
      // serial dual simplex
      options.simplex_strategy = kSimplexStrategyDualPlain;
      options.simplex_max_concurrency = 1;
      // One of the concurrent simplex solvers is running on this
      // Highs instance, so reduce the number of solvers that can be
      // spawned
      spawnable_solvers--;
    }
    // If there is only one thread avaiable for parallel dual simplex,
    // then reduce the number of simplex that can be spawned -
    // preventing it from being used
    if (spawnable_solvers == spawnable_threads)
      spawnable_solvers--;
    printf(
        "Spawning %d LP solvers using (upto) %d threads: \n",
        int(spawnable_solvers), int(spawnable_threads));
    assert(spawnable_solvers > 0);
  }
  // Determine which solvers are to be spawned
  vector<HighsInt> spawned_solver;
  if (spawnable_solvers) {
    HighsInt available_solvers = spawnable_solvers;
    HighsInt available_threads = spawnable_threads;
    if (use_solver == kIpmString) {
      // This thread is running IPM, so need a thread running serial
      // dual simplex
      assert(available_solvers);
      assert(available_threads);
      spawned_solver.push_back(kSimplexStrategyDualPlain);
      available_solvers--;
      available_threads--;
    } else {
      // This thread must be using serial dual simplex
      assert(options.simplex_strategy == kSimplexStrategyDualPlain);
      assert(options.simplex_max_concurrency == 1);
    }
    if (available_solvers) {
      // Have a thread running primal simplex
      assert(available_solvers);
      assert(available_threads);
      spawned_solver.push_back(kSimplexStrategyPrimal);
      available_solvers--;
      available_threads--;
    }
  }
  HighsInt num_spawned_solver = spawned_solver.size();
  // Set up a vector of pointers to Highs instances
  vector<std::unique_ptr<Highs>> parallel_highs;
  for (HighsInt ix = 0; ix < num_spawned_solver; ix++) {
    parallel_highs.emplace_back(new Highs());
    parallel_highs[ix]->passSpawnId(ix);
    // Concurrent instances run simplex silently, without (further) scaling
    parallel_highs[ix]->passOptions(options);
    const bool debug_log = true;
    if (debug_log) {
      std::string log_file = "HiGHS" + std::to_string(ix) + ".log";
      parallel_highs[ix]->setOptionValue("log_to_console", false);
      parallel_highs[ix]->setOptionValue("log_file", log_file);
      parallel_highs[ix]->setOptionValue("log_dev_level", 1);
    } else {
      parallel_highs[ix]->setOptionValue("output_flag", false);
    }
    parallel_highs[ix]->setOptionValue("solver", kSimplexString);
    parallel_highs[ix]->setOptionValue("simplex_strategy", spawned_solver[ix]);
    //    parallel_highs[ix]->setOptionValue("simplex_scale_strategy", kSimplexScaleStrategyOff);
      parallel_highs[ix]->setOptionValue("parallel", kHighsOffString);
    if (spawned_solver[ix] == kSimplexStrategyDualTasks) {
      // Using PAMI, so force parallel to be on
      assert(111==222);
    }
    parallel_highs[ix]->passModel(solver_object.lp_);
  }

  std::vector<double> run_time;
  std::vector<double> objective_function_value;
  std::vector<HighsStatus> run_status;
  std::vector<HighsModelStatus> model_status;
  if (num_spawned_solver) {
    run_time.resize(num_spawned_solver);
    objective_function_value.resize(num_spawned_solver);
    model_status.resize(num_spawned_solver);
    run_status.assign(num_spawned_solver, HighsStatus::kError);
  }
  printf("Spawn Id = %2d: num_spawned_solver = %1d; use_solver = %7s; strategy = %d; parallel = %3s; max concurrency = %d\n",
	 int(solver_object.spawn_id_), int(num_spawned_solver),
	 use_solver.c_str(), int(options.simplex_strategy),
	 options.parallel.c_str(), int(options.simplex_max_concurrency));
  HighsRaceTimer<double> race_timer;
  parallel::TaskGroup tg;
  for (HighsInt solver = 0; solver < num_spawned_solver; solver++) {
    tg.spawn([&parallel_highs, &run_status, &model_status, &race_timer,
	      &run_time, solver]() {
	       run_time[solver] = parallel_highs[solver]->getRunTime();
	       parallel_highs[solver]->passRaceTimer(&race_timer);
	       run_status[solver] = parallel_highs[solver]->run();
	       model_status[solver] = parallel_highs[solver]->getModelStatus();
	       // this should check for the status returned when the race timer limit
	       // was reached and only call decrease limit if it was not reached
	       printf("Spawned solver %2d (time = %11.4g) returns status %s and model status %s\n",
		      int(solver),
		      parallel_highs[solver]->getRunTime() - run_time[solver],
		      highsStatusToString(run_status[solver]).c_str(),
		      parallel_highs[solver]
		      ->modelStatusToString(model_status[solver])
		      .c_str());
	       fflush(stdout);
	       if (positiveModelStatus(model_status[solver])) {
		 run_time[solver] =
		   parallel_highs[solver]->getRunTime() - run_time[solver];
		 race_timer.decreaseLimit(run_time[solver]);
	       }
	     });
  }
  solver_object.run_time_ = solver_object.timer_.readRunHighsClock();
  if (use_solver == kIpmString) {
    // Use IPM for this Highs instance
    try {
      call_status = solveLpIpx(solver_object);
    } catch (const std::exception& exception) {
      highsLogDev(options.log_options, HighsLogType::kError,
                  "Exception %s in solveLpIpx\n", exception.what());
      call_status = HighsStatus::kError;
    }
    // Non-error return requires a primal solution
    assert(solver_object.solution_.value_valid);
    // Set the return_status, model status and, for completeness, scaled
    // model status
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "solveLpIpx");
    if (return_status == HighsStatus::kError) {
      // Can only allow this return if this is a spawned solver
      if (solver_object.spawn_id_ >=0) return solveLpReturn(return_status, solver_object, message);
      logLocalSolverOutcome(return_status, solver_object);
      assert(111==333);
    }
    // Get the objective and any KKT failures
    solver_object.highs_info_.objective_function_value =
        solver_object.lp_.objectiveValue(solver_object.solution_.col_value);
    getLpKktFailures(options, solver_object.lp_, solver_object.solution_,
                     solver_object.basis_, solver_object.highs_info_);
    // Seting the IPM-specific values of (highs_)info_ has been done in
    // solveLpIpx
    if ((solver_object.model_status_ == HighsModelStatus::kUnknown ||
         (solver_object.model_status_ ==
              HighsModelStatus::kUnboundedOrInfeasible &&
          !options.allow_unbounded_or_infeasible)) &&
        options.run_crossover) {
      // IPX has returned a model status that HiGHS would rather
      // avoid, so perform simplex clean-up if crossover was allowed.
      //
      // This is an unusual situation, and the cost will usually be
      // acceptable. Worst case is if crossover wasn't run, in which
      // case there's no basis to start simplex
      //
      // ToDo: Check whether simplex can exploit the primal solution returned by
      // IPX
      highsLogUser(options.log_options, HighsLogType::kWarning,
                   "Imprecise solution returned from IPX, so use simplex to "
                   "clean up\n");
      // Reset the return status since it will now be determined by
      // the outcome of the simplex solve
      return_status = HighsStatus::kOk;
      // ToDo: Need to ensure that the clean-up uses serial simplex
      options.simplex_max_concurrency = 1;
      call_status = solveLpSimplex(solver_object);
      // Restore simplex_max_concurrency
      options.simplex_max_concurrency = simplex_max_concurrency;
      return_status = interpretCallStatus(options.log_options, call_status,
                                          return_status, "solveLpSimplex");
      if (return_status == HighsStatus::kError) {
	// Can only allow this return if this is a spawned solver
	if (solver_object.spawn_id_ >= 0) return solveLpReturn(return_status, solver_object, message);
	logLocalSolverOutcome(return_status, solver_object);
	assert(111==444);
      }
    }
  } else {
    assert(use_solver == kSimplexString || use_solver == kHighsChooseString);
    // If there are concurrent simplex solvers, this thread must be
    // serial dual simplex
    assert(spawnable_solvers == 0 ||
	   (options.simplex_max_concurrency == 1 &&
	    options.simplex_strategy == kSimplexStrategyDualPlain));
    // Use Simplex
    call_status = solveLpSimplex(solver_object);
    return_status = interpretCallStatus(options.log_options, call_status,
                                        return_status, "solveLpSimplex");
    // Restore simplex_max_concurrency and simplex_strategy
    options.simplex_max_concurrency = simplex_max_concurrency;
    options.simplex_strategy = simplex_strategy;
    if (return_status == HighsStatus::kError) {
      // Can only allow this return if this is a spawned solver
      if (solver_object.spawn_id_ >=0) return solveLpReturn(return_status, solver_object, message);
      logLocalSolverOutcome(return_status, solver_object);
      assert(111==555);
    }
  }
  // Check for solution consistency
  if (!isSolutionRightSize(solver_object.lp_, solver_object.solution_)) {
    return_status = HighsStatus::kError;
    highsLogUser(options.log_options, HighsLogType::kError,
		 "Inconsistent solution returned from solver\n");
    // Can only allow this return if this is a spawned solver
    if (solver_object.spawn_id_ >= 0) solveLpReturn(return_status, solver_object, message); 
    logLocalSolverOutcome(return_status, solver_object);
    assert(111==666);
  }
  double solver_run_time = solver_object.timer_.readRunHighsClock() - solver_object.run_time_;
  if (positiveModelStatus(solver_object.model_status_))
    race_timer.decreaseLimit(solver_run_time);
  if (solver_object.spawn_id_ < 0) {
    double save_run_time = solver_object.run_time_;
    solver_object.run_time_ = solver_run_time;
    if (num_spawned_solver) logLocalSolverOutcome(return_status, solver_object);
    solver_object.run_time_ = save_run_time;
  }
  tg.taskWait();
  if (num_spawned_solver) {
    HighsModelStatus have_model_status = solver_object.model_status_;
    bool have_positive_model_status = positiveModelStatus(have_model_status);
    bool have_bound_termination = boundTermination(have_model_status);
    // Look for positive model status elsewhere and either take it or
    // check that it's the same as what's known
    for (HighsInt solver = 0; solver < num_spawned_solver; solver++) {
      const HighsModelStatus this_model_status = model_status[solver];
      const bool this_bound_termination = boundTermination(this_model_status);
      if (have_positive_model_status) {
	if (positiveModelStatus(this_model_status)) {
	  // Here's another positive model status: should be the
	  // same, unless one is a bound termination
	  if (!this_bound_termination || !have_bound_termination) {
	    if (this_model_status != have_model_status) {
	      highsLogUser(options.log_options, HighsLogType::kError,
			   "Inconsistent positive model status from concurrent solvers: (%s vs %s)\n",
			   utilModelStatusToString(this_model_status),
			   utilModelStatusToString(have_model_status));
	      return solveLpReturn(HighsStatus::kError, solver_object, message);
	    }	    
	  }
	}
      } else {
	// Here's the first positive model status: extract solution
	// and basis if available
	solver_object.basis_ = parallel_highs[solver]->getBasis();
	solver_object.solution_ = parallel_highs[solver]->getSolution();
	solver_object.highs_info_ = parallel_highs[solver]->getInfo();
	//	solver_object.ekk_instance_ = parallel_highs[solver]->getEkk();
	have_model_status = this_model_status;
	have_bound_termination = this_bound_termination;
	highsLogUser(options.log_options, HighsLogType::kInfo,
			   "Solution obtained from spawned %s solver\n",
		     simplexStrategyToString(spawned_solver[solver]).c_str());
	assert(111==000);
      }
    }
  }
  return solveLpReturn(return_status, solver_object, message);
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsInfo
HighsStatus solveUnconstrainedLp(HighsLpSolverObject& solver_object) {
  return (solveUnconstrainedLp(solver_object.options_, solver_object.lp_,
                               solver_object.model_status_,
                               solver_object.highs_info_,
                               solver_object.solution_, solver_object.basis_));
}

// Solves an unconstrained LP without scaling, setting HighsBasis, HighsSolution
// and HighsInfo
HighsStatus solveUnconstrainedLp(const HighsOptions& options, const HighsLp& lp,
                                 HighsModelStatus& model_status,
                                 HighsInfo& highs_info, HighsSolution& solution,
                                 HighsBasis& basis) {
  // Aliase to model status and solution parameters
  resetModelStatusAndHighsInfo(model_status, highs_info);

  // Check that the LP really is unconstrained!
  assert(lp.num_row_ == 0);
  if (lp.num_row_ != 0) return HighsStatus::kError;

  highsLogUser(options.log_options, HighsLogType::kInfo,
               "Solving an unconstrained LP with %" HIGHSINT_FORMAT
               " columns\n",
               lp.num_col_);

  solution.col_value.assign(lp.num_col_, 0);
  solution.col_dual.assign(lp.num_col_, 0);
  basis.col_status.assign(lp.num_col_, HighsBasisStatus::kNonbasic);
  // No rows for primal solution, dual solution or basis
  solution.row_value.clear();
  solution.row_dual.clear();
  basis.row_status.clear();

  double primal_feasibility_tolerance = options.primal_feasibility_tolerance;
  double dual_feasibility_tolerance = options.dual_feasibility_tolerance;

  // Initialise the objective value calculation. Done using
  // HighsSolution so offset is vanilla
  double objective = lp.offset_;
  bool infeasible = false;
  bool unbounded = false;

  highs_info.num_primal_infeasibilities = 0;
  highs_info.max_primal_infeasibility = 0;
  highs_info.sum_primal_infeasibilities = 0;
  highs_info.num_dual_infeasibilities = 0;
  highs_info.max_dual_infeasibility = 0;
  highs_info.sum_dual_infeasibilities = 0;

  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    double cost = lp.col_cost_[iCol];
    double dual = (HighsInt)lp.sense_ * cost;
    double lower = lp.col_lower_[iCol];
    double upper = lp.col_upper_[iCol];
    double value;
    double primal_infeasibility = 0;
    double dual_infeasibility = -1;
    HighsBasisStatus status = HighsBasisStatus::kNonbasic;
    if (lower > upper) {
      // Inconsistent bounds, so set the variable to lower bound,
      // unless it's infinite. Otherwise set the variable to upper
      // bound, unless it's infinite. Otherwise set the variable to
      // zero.
      if (highs_isInfinity(lower)) {
        // Lower bound of +inf
        if (highs_isInfinity(-upper)) {
          // Upper bound of -inf
          value = 0;
          status = HighsBasisStatus::kZero;
          primal_infeasibility = kHighsInf;
          dual_infeasibility = std::fabs(dual);
        } else {
          // Finite upper bound - since lower exceeds it
          value = upper;
          status = HighsBasisStatus::kUpper;
          primal_infeasibility = lower - value;
          dual_infeasibility = std::max(dual, 0.);
        }
      } else {
        // Finite lower bound
        value = lower;
        status = HighsBasisStatus::kLower;
        primal_infeasibility = value - upper;
        dual_infeasibility = std::max(-dual, 0.);
      }
    } else if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free column: set to zero and record dual infeasiblility
      value = 0;
      status = HighsBasisStatus::kZero;
      dual_infeasibility = std::fabs(dual);
    } else if (dual >= dual_feasibility_tolerance) {
      // Column with sufficiently positive dual
      if (!highs_isInfinity(-lower)) {
        // Set to this finite lower bound
        value = lower;
        status = HighsBasisStatus::kLower;
        dual_infeasibility = 0;
      } else {
        // Infinite lower bound so set to upper bound and record dual
        // infeasiblility
        value = upper;
        status = HighsBasisStatus::kUpper;
        dual_infeasibility = dual;
      }
    } else if (dual <= -dual_feasibility_tolerance) {
      // Column with sufficiently negative dual
      if (!highs_isInfinity(upper)) {
        // Set to this finite upper bound
        value = upper;
        status = HighsBasisStatus::kUpper;
        dual_infeasibility = 0;
      } else {
        // Infinite upper bound so set to lower bound and record dual
        // infeasiblility
        value = lower;
        status = HighsBasisStatus::kLower;
        dual_infeasibility = -dual;
      }
    } else {
      // Column with sufficiently small dual: set to lower bound (if
      // finite) otherwise upper bound
      if (highs_isInfinity(-lower)) {
        value = upper;
        status = HighsBasisStatus::kUpper;
      } else {
        value = lower;
        status = HighsBasisStatus::kLower;
      }
      dual_infeasibility = std::fabs(dual);
    }
    assert(status != HighsBasisStatus::kNonbasic);
    assert(dual_infeasibility >= 0);
    solution.col_value[iCol] = value;
    solution.col_dual[iCol] = (HighsInt)lp.sense_ * dual;
    basis.col_status[iCol] = status;
    objective += value * cost;
    if (primal_infeasibility > primal_feasibility_tolerance)
      highs_info.num_primal_infeasibilities++;
    highs_info.sum_primal_infeasibilities += primal_infeasibility;
    highs_info.max_primal_infeasibility =
        std::max(primal_infeasibility, highs_info.max_primal_infeasibility);
    if (dual_infeasibility > dual_feasibility_tolerance)
      highs_info.num_dual_infeasibilities++;
    highs_info.sum_dual_infeasibilities += dual_infeasibility;
    highs_info.max_dual_infeasibility =
        std::max(dual_infeasibility, highs_info.max_dual_infeasibility);
  }
  highs_info.objective_function_value = objective;
  solution.value_valid = true;
  solution.dual_valid = true;
  basis.valid = true;
  highs_info.basis_validity = kBasisValidityValid;
  setSolutionStatus(highs_info);
  if (highs_info.num_primal_infeasibilities) {
    // Primal infeasible
    model_status = HighsModelStatus::kInfeasible;
  } else if (highs_info.num_dual_infeasibilities) {
    // Dual infeasible => primal unbounded for unconstrained LP
    model_status = HighsModelStatus::kUnbounded;
  } else {
    model_status = HighsModelStatus::kOptimal;
  }

  return HighsStatus::kOk;
}
