/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Highs.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Highs.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>

#include "HConfig.h"
#include "io/Filereader.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsStatus.h"
#include "presolve/FindFeasibility.h"
#include "presolve/Presolve.h"
#include "simplex/HApp.h"
#include "simplex/HighsSimplexInterface.h"

// until add_row_.. functions are moved to HighsLpUtils.h
#include "simplex/HSimplex.h"

Highs::Highs() {
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  simplex_has_run_ = false;
}

  HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                  const std::string& value) {
    OptionStatus status;
    HighsOptionType type;
    status = getOptionType(option, type);
    switch( type )
      {
      case HighsOptionType::BOOL:
	printf("ERROR: No method to set option %s of type bool\n", option.c_str());
	return HighsStatus::Error;
      case HighsOptionType::INT:
	status = setOptionValue(options_, option, atoi(value.c_str()));
	break;
      case HighsOptionType::DOUBLE:
	status = setOptionValue(options_, option, atof(value.c_str()));
	break;
      case HighsOptionType::STRING:
	status = setOptionValue(options_, option, value);
	break;
      default:
	printf("ERROR: No method to set option %s of unknown type %d\n", option.c_str(), (int)type);
	return HighsStatus::Error;
      }
    if (status == OptionStatus::OK) return HighsStatus::OK;
    return HighsStatus::Error;
  }

  HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                  const double& value) {
    OptionStatus status = setOptionValue(options_, option, value);
    if (status == OptionStatus::OK) return HighsStatus::OK;
    return HighsStatus::Error;
  }

  HighsStatus Highs::setHighsOptionValue(const std::string& option,
                                  const int& value) {
    OptionStatus status = setOptionValue(options_, option, value);
    if (status == OptionStatus::OK) return HighsStatus::OK;
    return HighsStatus::Error;
  }

  HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                  std::string& value) {
    OptionStatus status;
    HighsOptionType type;
    std::stringstream stringstream_value;
    status = getOptionType(option, type);
    switch( type )
      {
      case HighsOptionType::BOOL:
	printf("ERROR: No method to get option %s of type bool\n", option.c_str());
	return HighsStatus::Error;
      case HighsOptionType::INT:
	int int_value;
	status = getOptionValue(options_, option, int_value);
	if (status == OptionStatus::OK) {
	  stringstream_value << int_value;
	  value = stringstream_value.str();
	}
	break;
      case HighsOptionType::DOUBLE:
	double double_value;
	status = getOptionValue(options_, option, double_value);
	if (status == OptionStatus::OK) {
	  stringstream_value << double_value;
	  value = stringstream_value.str();
	}
	break;
      case HighsOptionType::STRING:
	status = getOptionValue(options_, option, value);
	break;
      default:
	printf("ERROR: No method to get option %s of unknown type %d\n", option.c_str(), (int)type);
	return HighsStatus::Error;
      }
    if (status == OptionStatus::OK) return HighsStatus::OK;
    return HighsStatus::Error;
  }
  HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                  double& value) {
    OptionStatus status = getOptionValue(options_, option, value);
    printf("getHighsOptionValue(%s, %g)\n", option.c_str(), value);
    if (status == OptionStatus::OK) return HighsStatus::OK;
    return HighsStatus::Error;
  }

  HighsStatus Highs::getHighsOptionValue(const std::string& option,
                                  int& value) {
    OptionStatus status = getOptionValue(options_, option, value);
    if (status == OptionStatus::OK) return HighsStatus::OK;
    return HighsStatus::Error;
  }


HighsStatus Highs::initializeLp(const HighsLp& lp) {
  // todo:(julian) add code to check that LP is valid.
  lp_ = lp;

  // hmos_[0] is the HighsModelObject corresponding to the original LP
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  simplex_has_run_ = false;
  return HighsStatus::OK;
}

HighsStatus Highs::initializeFromFile(const std::string filename) {
  Filereader* reader = Filereader::getFilereader(filename.c_str());
  HighsLp model;
  this->options_.filename = filename;

  FilereaderRetcode retcode = reader->readModelFromFile(this->options_, model);
  if (retcode != FilereaderRetcode::OK) {
    return HighsStatus::Error;
  }

  return this->initializeLp(model);
}

HighsStatus Highs::writeToFile(const std::string filename) {
  HighsLp model = this->lp_;

  Filereader* writer = Filereader::getFilereader(filename.c_str());
  FilewriterRetcode retcode = writer->writeModelToFile(filename.c_str(), model);
  if (retcode != FilewriterRetcode::OK) {
    return HighsStatus::Error;
  }

  return HighsStatus::OK;
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run() {
  // If running as hsol, reset any changed options
  if (options_.run_as_hsol) setHsolOptions(options_);
  // Assess the LP, normalising its values
  bool normalise = true;
  HighsStatus return_status = assessLp(lp_, options_, normalise);
  if (return_status == HighsStatus::Error) return return_status;

  // For the moment runFeasibility as standalone.
  if (options_.find_feasibility) {
    // use when you do something with solution depending on whether we have
    // dualized or not.
    // HighsSolution& solution = solution_;

    // options_.messageLevel = HighsPrintMessageLevel::ML_DETAILED;
    // HighsSetIO(options_);

    // Add slacks and make sure a minimization problem is passed to
    // runFeasibility.
    HighsLp primal = transformIntoEqualityProblem(lp_);
    if (options_.feasibility_strategy_dualize) {
      // Add slacks & dualize.
      HighsLp dual = dualizeEqualityProblem(primal);
      // dualizeEqualityProblem returns a minimization problem.
      initializeLp(dual);
    } else {
      // If maximization, minimize before calling runFeasibility.
      if (primal.sense_ != OBJSENSE_MINIMIZE) {
        for (int col = 0; col < primal.numCol_; col++)
          primal.colCost_[col] = -primal.colCost_[col];
      }
      initializeLp(primal);
    }

    if (options_.feasibility_strategy ==
        FeasibilityStrategy::kApproxComponentWise)
      return runFeasibility(lp_, solution_, MinimizationType::kComponentWise);
    else if (options_.feasibility_strategy == FeasibilityStrategy::kApproxExact)
      return runFeasibility(lp_, solution_, MinimizationType::kExact);
    else if (options_.feasibility_strategy ==
             FeasibilityStrategy::kDirectSolve) {
      // Proceed to normal exection of run().
      // If dualize has been called replace LP is replaced with dual in code
      // above.
    }
  }

  // Return immediately if the LP has no columns
  if (!lp_.numCol_) return HighsStatus::LpEmpty;

  // todo: check options.
  HighsSetIO(options_);

  reportOptionsValue(options_, 0);
  HighsPrintMessage(ML_VERBOSE, "Solving %s", lp_.model_name_.c_str());
  if (options_.mip) return runBnb();

  // Running as LP solver: start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  // Record the initial time and zero the overall iteration count
  double initial_time = timer_.readRunHighsClock();
  int solve_iteration_count = 0;
  int postsolve_iteration_count = 0;
  // todo: make sure it should remain Init for calls of run() after
  // simplex_has_run_ is valid.
  HighsStatus solve_status = HighsStatus::Init;
  // Define identifiers to refer to the HMO of the original LP
  // (0) and the HMO created when using presolve (1)
  const int original_hmo = 0;
  const int presolve_hmo = 1;
  // Keep track of the hmo that is the most recently solved. By default it's the
  // original LP
  int solved_hmo = original_hmo;

  // Initial solve. Presolve, choose solver (simplex, ipx), postsolve.
  int iteration_count;
  if (!simplex_has_run_) {
    // Presolve. runPresolve handles the level of presolving (0 = don't
    // presolve).
    timer_.start(timer_.presolve_clock);
    PresolveInfo presolve_info(options_.presolve_option, lp_);
    HighsPresolveStatus presolve_status = runPresolve(presolve_info);
    timer_.stop(timer_.presolve_clock);

    // Run solver.
    switch (presolve_status) {
      case HighsPresolveStatus::NotPresolved: {
	hmos_[solved_hmo].lp_.lp_name_ = "Original LP";
        solve_status = callRunSolver(hmos_[solved_hmo], iteration_count,
                                     "Not presolved: solving the LP");
        solve_iteration_count += iteration_count;
        break;
      }
      case HighsPresolveStatus::NotReduced: {
	hmos_[solved_hmo].lp_.lp_name_ = "Unreduced LP";
        solve_status = callRunSolver(hmos_[solved_hmo], iteration_count,
				     "Problem not reduced by presolve: solving the LP");
        solve_iteration_count += iteration_count;
        break;
      }
      case HighsPresolveStatus::Reduced: {
        HighsLp& reduced_lp = presolve_info.getReducedProblem();
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.
        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        // Log the presolve reductions
        logPresolveReductions(hmos_[original_hmo].lp_, hmos_[presolve_hmo].lp_);
        // Record the HMO to be solved
        solved_hmo = presolve_hmo;
	hmos_[solved_hmo].lp_.lp_name_ = "Presolved LP";
        solve_status = callRunSolver(hmos_[solved_hmo], iteration_count,
                                     "Solving the presolved LP");
        solve_iteration_count += iteration_count;
        break;
      }
      case HighsPresolveStatus::ReducedToEmpty: {
        // Proceed to postsolve.
        break;
      }
      case HighsPresolveStatus::Infeasible:
      case HighsPresolveStatus::Unbounded: {
        HighsStatus result =
            (presolve_status == HighsPresolveStatus::Infeasible)
                ? HighsStatus::Infeasible
                : HighsStatus::Unbounded;
        HighsPrintMessage(ML_ALWAYS,
                          "Problem status detected on presolve: %s\n",
                          HighsStatusToString(result).c_str());

        // Report this way for the moment. May modify after merge with
        // OSIinterface branch which has new way of setting up a
        // HighsModelObject and can support multiple calls to run(). Stop and
        // read the HiGHS clock, then work out time for this call
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
        double lp_solve_final_time = timer_.readRunHighsClock();

        std::stringstream message_not_opt;
        message_not_opt << std::endl;
        message_not_opt << "Run status : " << HighsStatusToString(result)
                        << std::endl;
        message_not_opt << "Time       : " << std::fixed << std::setprecision(3)
                        << lp_solve_final_time - initial_time << std::endl;

        message_not_opt << std::endl;

        HighsPrintMessage(ML_MINIMAL, message_not_opt.str().c_str());
        return result;
      }
      default: {
        // case HighsPresolveStatus::Error
        HighsPrintMessage(ML_ALWAYS, "Presolve failed.");
        if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
        return HighsStatus::PresolveError;
      }
    }
    // Postsolve. Does nothing if there were no reductions during presolve.
    if (solve_status == HighsStatus::Optimal) {
      if (presolve_status == HighsPresolveStatus::Reduced) {
        // If presolve is nontrivial, extract the optimal solution
        // and basis for the presolved problem in order to generate
        // the solution and basis for postsolve to use to generate a
        // solution(?) and basis that is, hopefully, optimal. This is
        // confirmed or corrected by hot-starting the simplex solver
        presolve_info.reduced_solution_ = hmos_[solved_hmo].solution_;
        presolve_info.presolve_[0].setBasisInfo(
            hmos_[solved_hmo].basis_.col_status,
            hmos_[solved_hmo].basis_.row_status);
        // Run postsolve
        timer_.start(timer_.postsolve_clock);
        HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
        timer_.stop(timer_.postsolve_clock);
        if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
          HighsPrintMessage(ML_VERBOSE, "Postsolve finished.");
          // Set solution(?) and basis to hot-start the simplex solver
          // for the original_hmo
          hmos_[original_hmo].solution_ = presolve_info.recovered_solution_;

          hmos_[original_hmo].basis_.col_status =
              presolve_info.presolve_[0].getColStatus();
          hmos_[original_hmo].basis_.row_status =
              presolve_info.presolve_[0].getRowStatus();
          hmos_[original_hmo].basis_.valid_ = true;
          // Analyse the Highs basic solution returned from postsolve
          HighsSimplexInterface simplex_interface(hmos_[original_hmo]);
          int report_level = -1;
#ifdef HiGHSDEV
          report_level = 1;
#endif
          simplex_interface.analyseHighsSolutionAndBasis(
              report_level, "after returning from postsolve");
          // Now hot-start the simplex solver for the original_hmo
          solved_hmo = original_hmo;
          // Save the options to allow the best simplex strategy to
          // be used
          HighsOptions save_options = hmos_[solved_hmo].options_;
          HighsOptions& options = hmos_[solved_hmo].options_;
          options.simplex_strategy = SimplexStrategy::CHOOSE;
          // Set the message level to ML_ALWAYS so that data for
          // individual iterations are reported
          bool full_iteration_logging = false;
          if (full_iteration_logging) HighsSetMessagelevel(ML_ALWAYS);
	  hmos_[solved_hmo].lp_.lp_name_ = "Postsolve LP";
          solve_status = callRunSolver(hmos_[solved_hmo], iteration_count,
				       "Solving the original LP from the solution after postsolve");
          postsolve_iteration_count = iteration_count;
          solve_iteration_count += iteration_count;
          // Recover the options
          options = save_options;
          // Reset the message level
          if (full_iteration_logging)
            HighsSetMessagelevel(options_.messageLevel);
        }
      }
    }
  } else {
    // The problem has been solved before so we ignore presolve/postsolve/ipx.
    solved_hmo = original_hmo;
    hmos_[solved_hmo].lp_.lp_name_ = "Re-solved LP";
    solve_status = callRunSolver(hmos_[solved_hmo], iteration_count,
				 "Re-solving the LP");
    solve_iteration_count += iteration_count;
  }
  // else if (reduced problem failed to solve) {
  //   todo: handle case when presolved problem failed to solve. Try to solve
  //   again with no presolve.
  // }

  assert(hmos_.size() > 0);
  solution_ = hmos_[original_hmo].solution_;
  basis_ = hmos_[original_hmo].basis_;

  // Report times
  if (hmos_[original_hmo].report_model_operations_clock) {
    std::vector<int> clockList{
        timer_.presolve_clock,
	timer_.solve_clock,
	timer_.postsolve_clock};
    timer_.report("ModelOperations", clockList);
  }
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
  double lp_solve_final_time = timer_.readRunHighsClock();

  std::stringstream message;
  message << std::endl;
  message << "Run status : " << HighsStatusToString(solve_status) << std::endl;
  message
      << "Iterations : "
      << solve_iteration_count  // hmos_[solved_hmo].simplex_info_.iteration_count
      << std::endl;

  if (solve_status == HighsStatus::Optimal)
    message << "Objective  : " << std::scientific
            << hmos_[original_hmo].simplex_info_.dual_objective_value
            << std::endl;

  message << "Time       : " << std::fixed << std::setprecision(3)
          << lp_solve_final_time - initial_time << std::endl;

  message << "Postsolve  : " << postsolve_iteration_count << std::endl;

  message << std::endl;

  HighsPrintMessage(ML_MINIMAL, message.str().c_str());

  return solve_status;
}

const HighsLp& Highs::getLp() const { return lp_; }

const HighsSolution& Highs::getSolution() const { return solution_; }

const HighsBasis& Highs::getBasis() const { return basis_; }

double Highs::getObjectiveValue() const {
  if (hmos_.size() > 0) {
    return hmos_[0].simplex_info_.dual_objective_value;
  } else {
    // todo: ipx case
    // todo: error/warning message
  }
  return 0;
}

int Highs::getIterationCount() const {
  if (hmos_.size() == 0) return 0;
  return hmos_[0].simplex_info_.iteration_count;
}

HighsStatus Highs::setSolution(const HighsSolution& solution) {
  underDevelopmentLogMessage("setSolution");
  // Check if solution is valid.
  assert((int)solution_.col_value.size() != 0 ||
         (int)solution_.col_value.size() != lp_.numCol_);
  assert((int)solution.col_dual.size() == 0 ||
         (int)solution.col_dual.size() == lp_.numCol_);
  assert((int)solution.row_dual.size() == 0 ||
         (int)solution.row_dual.size() == lp_.numRow_);

  if (solution.col_value.size()) solution_.col_value = solution.col_value;
  if (solution.col_dual.size()) solution_.col_dual = solution.col_dual;
  if (solution.row_dual.size()) solution_.row_dual = solution.row_dual;

  HighsStatus result_values = HighsStatus::NotSet;
  HighsStatus result_duals = HighsStatus::NotSet;

  if (solution.col_value.size() > 0)
    result_values = calculateRowValues(lp_, solution_);
  if (solution.row_dual.size() > 0)
    result_duals = calculateColDuals(lp_, solution_);

  if (result_values == HighsStatus::Error || result_duals == HighsStatus::Error)
    return HighsStatus::Error;

  return HighsStatus::OK;
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  underDevelopmentLogMessage("setBasis");
  basis_ = basis;
  return HighsStatus::OK;
}

void Highs::reportSolution() {
  reportModelBoundSol(true, lp_.numCol_, lp_.colLower_, lp_.colUpper_,
                      lp_.col_names_, solution_.col_value, solution_.col_dual,
                      basis_.col_status);
  reportModelBoundSol(false, lp_.numRow_, lp_.rowLower_, lp_.rowUpper_,
                      lp_.row_names_, solution_.row_value, solution_.row_dual,
                      basis_.row_status);
}

bool Highs::addRow(const double lower_bound, const double upper_bound,
                   const int num_new_nz, const int* indices,
                   const double* values) {
  int starts = 0;
  return addRows(1, &lower_bound, &upper_bound, num_new_nz, &starts, indices,
                 values);
}

bool Highs::addRows(const int num_new_row, const double* lower_bounds,
                    const double* upper_bounds, const int num_new_nz,
                    const int* starts, const int* indices,
                    const double* values) {
  underDevelopmentLogMessage("addRows");
  HighsStatus return_status = HighsStatus::NotSet;
  // if simplex has not solved already
  if (!simplex_has_run_) {
    return_status = addLpRows(lp_, num_new_row, lower_bounds, upper_bounds,
                              num_new_nz, starts, indices, values, options_);
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    // todo: change to take int return value
    return_status = interface.addRows(num_new_row, lower_bounds, upper_bounds,
                                      num_new_nz, starts, indices, values);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::addCol(const double cost, const double lower_bound,
                   const double upper_bound, const int num_new_nz,
                   const int* indices, const double* values) {
  int starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, num_new_nz, &starts,
                 indices, values);
}

bool Highs::addCols(const int num_new_col, const double* costs,
                    const double* lower_bounds, const double* upper_bounds,
                    const int num_new_nz, const int* starts, const int* indices,
                    const double* values) {
  underDevelopmentLogMessage("addCols");
  HighsStatus return_status = HighsStatus::NotSet;
  // if simplex has not solved already
  if (!simplex_has_run_) {
    return_status =
        addLpCols(lp_, num_new_col, costs, lower_bounds, upper_bounds,
                  num_new_nz, starts, indices, values, options_);
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    // todo: change to take int return value
    return_status =
        interface.addCols(num_new_col, costs, lower_bounds, upper_bounds,
                          num_new_nz, starts, indices, values);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeObjectiveSense(const int sense) {
  underDevelopmentLogMessage("changeObjectiveSense");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    this->lp_.sense_ = sense;
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.changeObjectiveSense(sense);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeColCost(const int col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const int num_set_entries, const int* set,
                           const double* cost) {
  underDevelopmentLogMessage("changeColsCost");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i = 0; i < num_set_entries; i++) {
      this->lp_.colCost_[set[i]] = cost[i];
    }
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeCosts(num_set_entries, set, cost);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeColsCost(const int* mask, const double* cost) {
  underDevelopmentLogMessage("changeColsCost");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i = 0; i < this->lp_.numCol_; i++) {
      if (mask[i]) this->lp_.colCost_[i] = cost[i];
    }
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeCosts(mask, cost);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeColBounds(const int col, const double lower,
                            const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  underDevelopmentLogMessage("changeColsBounds");
  HighsStatus return_status = HighsStatus::NotSet;
  /*
  if (!simplex_has_run_) {
    printf("changeColsBounds: Simplex has not run\n");
    for (int i = 0; i < num_set_entries; i++) {
      this->lp_.colLower_[set[i]] = lower[i];
      this->lp_.colUpper_[set[i]] = upper[i];
    }

  } else {
  */
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status =
        interface.changeColBounds(num_set_entries, set, lower, upper);
    //  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeColsBounds(const int from_col, const int to_col,
                             const double* lower, const double* upper) {
  underDevelopmentLogMessage("changeColsBounds");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    int num_changed_bounds = to_col - from_col;
    for (int i = 0; i < num_changed_bounds; i++) {
      this->lp_.colLower_[from_col + i] = lower[i];
      this->lp_.colUpper_[from_col + i] = upper[i];
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeColBounds(from_col, to_col, lower, upper);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeColsBounds(const int* mask, const double* lower,
                             const double* upper) {
  underDevelopmentLogMessage("changeColsBounds");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i = 0; i < this->lp_.numCol_; i++) {
      if (mask[i]) {
        this->lp_.colLower_[i] = lower[i];
        this->lp_.colUpper_[i] = upper[i];
      }
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeColBounds(mask, lower, upper);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeRowBounds(const int row, const double lower,
                            const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const int num_set_entries, const int* set,
                             const double* lower, const double* upper) {
  underDevelopmentLogMessage("changeRowsBounds");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i = 0; i < num_set_entries; i++) {
      this->lp_.rowLower_[set[i]] = lower[i];
      this->lp_.rowUpper_[set[i]] = upper[i];
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status =
        interface.changeRowBounds(num_set_entries, set, lower, upper);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeRowsBounds(const int* mask, const double* lower,
                             const double* upper) {
  underDevelopmentLogMessage("changeRowsBounds");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i = 0; i < this->lp_.numRow_; i++) {
      if (mask[i]) {
        this->lp_.rowLower_[i] = lower[i];
        this->lp_.rowUpper_[i] = upper[i];
      }
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeRowBounds(mask, lower, upper);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::changeCoeff(const int row, const int col, const double value) {
  underDevelopmentLogMessage("changeCoeff");
  HighsStatus return_status = HighsStatus::NotSet;
  assert(hmos_.size() > 0);
  HighsSimplexInterface interface(hmos_[0]);

  return_status = interface.changeCoefficient(row, col, value);
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
  
}

bool Highs::getCols(const int from_col, const int to_col, int& num_col,
                    double* costs, double* lower, double* upper, int& num_nz,
                    int* start, int* index, double* value) {
  underDevelopmentLogMessage("getCols");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getCols(from_col, to_col, num_col, costs, lower,
                                      upper, num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::getCols(const int n, const int* set, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getCols");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getCols(n, set, num_col, costs, lower, upper,
                                      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::getCols(const int* col_mask, int& num_col, double* costs,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getCols");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getCols(col_mask, num_col, costs, lower, upper,
                                      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::getRows(const int from_row, const int to_row, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getRows");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getRows(from_row, to_row, num_row, lower, upper,
                                      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::getRows(const int num_set_entries, const int* set, int& num_row,
                    double* lower, double* upper, int& num_nz, int* start,
                    int* index, double* value) {
  underDevelopmentLogMessage("getRows");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getRows(num_set_entries, set, num_row, lower,
                                      upper, num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::getRows(const int* mask, int& num_row, double* lower, double* upper,
                    int& num_nz, int* start, int* index, double* value) {
  underDevelopmentLogMessage("getRows");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getRows(mask, num_row, lower, upper, num_nz,
                                      start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::deleteCols(const int from_col, const int to_col) {
  underDevelopmentLogMessage("deleteCols");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.deleteCols(from_col, to_col);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::deleteCols(const int num_set_entries, const int* set) {
  underDevelopmentLogMessage("deleteCols");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.deleteCols(num_set_entries, set);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::deleteCols(int* mask) {
  underDevelopmentLogMessage("deleteCols");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.deleteCols(mask);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::deleteRows(const int from_row, const int to_row) {
  underDevelopmentLogMessage("deleteRows");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.deleteRows(from_row, to_row);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::deleteRows(const int num_set_entries, const int* set) {
  underDevelopmentLogMessage("deleteRows");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.deleteRows(num_set_entries, set);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

bool Highs::deleteRows(int* mask) {
  underDevelopmentLogMessage("deleteRows");
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.deleteRows(mask);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet)
    return false;
  return true;
}

// Private methods
HighsPresolveStatus Highs::runPresolve(PresolveInfo& info) {
  if (options_.presolve_option != PresolveOption::ON)
    return HighsPresolveStatus::NotPresolved;

  if (info.lp_ == nullptr) return HighsPresolveStatus::NullError;

  if (info.presolve_.size() == 0) return HighsPresolveStatus::NotReduced;

  info.presolve_[0].load(*(info.lp_));

  // Initialize a new presolve class instance for the LP given in presolve info
  return info.presolve_[0].presolve();
}

HighsPostsolveStatus Highs::runPostsolve(PresolveInfo& info) {
  if (info.presolve_.size() != 0) {
    bool solution_ok =
        isSolutionConsistent(info.getReducedProblem(), info.reduced_solution_);
    if (!solution_ok)
      return HighsPostsolveStatus::ReducedSolutionDimenionsError;

    // todo: error handling + see todo in run()
    info.presolve_[0].postsolve(info.reduced_solution_,
                                info.recovered_solution_);

    return HighsPostsolveStatus::SolutionRecovered;
  } else {
    return HighsPostsolveStatus::NoPostsolve;
  }
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus Highs::callRunSolver(HighsModelObject& model, int& iteration_count,
                                 const string message) {
  HighsLogMessage(HighsMessageType::INFO, message.c_str());
  // Handle the case of unconstrained LPs here
  if (!model.lp_.numRow_) {
    HighsSimplexInterface simplex_interface(model);
    iteration_count = 0;
    return simplex_interface.lpStatusToHighsStatus(solveUnconstrainedLp(model));
  }
  int initial_iteration_count = model.simplex_info_.iteration_count;
  HighsStatus solve_status = runSolver(model);
  int final_iteration_count = model.simplex_info_.iteration_count;
  iteration_count = final_iteration_count - initial_iteration_count;
  return solve_status;
}

// The method below runs simplex or ipx solver on the lp.
HighsStatus Highs::runSolver(HighsModelObject& model) {
  bool normalise = true;
  HighsStatus return_status = assessLp(model.lp_, model.options_, normalise);
  if (return_status == HighsStatus::Error) return return_status;

  HighsStatus status = HighsStatus::Init;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = solveModelSimplex(model);
  simplex_has_run_ = true;
#else
  // IPX
  // todo:Check options for simplex-specific options
  // use model.lp_, model.solution_
  // status = runIpxSolver(options_, lp_, solution_);
  // If ipx crossover did not find optimality set up simplex.

#endif

  if (status != HighsStatus::Optimal) return status;

  // Check.
  if (!isSolutionConsistent(model.lp_, model.solution_)) {
    std::cout << "Error: Inconsistent solution returned from solver.\n";
  }

  // todo:
  // assert(KktSatisfied(lp, solution));

  return status;
}

// Solve a mixed integer problem using branch and bound.
HighsStatus Highs::runBnb() {
  HighsPrintMessage(ML_ALWAYS, "Using branch and bound solver\n");

  // Need to start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  double mip_solve_initial_time = timer_.readRunHighsClock();

  // Start tree by making root node.
  std::unique_ptr<Node> root = std::unique_ptr<Node>(new Node(-1, 0, 0));

  root->integer_variables = lp_.integrality_;
  root->col_lower_bound = lp_.colLower_;
  root->col_upper_bound = lp_.colUpper_;

  HighsStatus status = solveRootNode(*(root.get()));
  if (status != HighsStatus::Optimal) {
    HighsPrintMessage(ML_ALWAYS,
                      "Root note not solved to optimality. Status: %s\n",
                      HighsStatusToString(status).c_str());
    return status;
  }

  // The method branch(...) below calls chooseBranchingVariable(..) which
  // currently returns the first violated one. If a branching variable is found
  // children are added to the stack. If there are no more violated integrality
  // constraints we have a feasible solution, if it is best than current best,
  // the current best is updated.
  int message_level = ML_DETAILED | ML_VERBOSE;
  options_.messageLevel = message_level;
  Tree tree(*(root.get()));
  tree.branch(*(root.get()));

  // While stack not empty.
  //   Solve node.
  //   Branch.
  while (!tree.empty()) {
    Node& node = tree.next();
    HighsStatus status = solveNode(node);
    tree.pop();

    if (status == HighsStatus::Infeasible) continue;

    options_.messageLevel = message_level;
    tree.branch(node);
  }

  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
  double mip_solve_final_time = timer_.readRunHighsClock();

  if (tree.getBestSolution().size() > 0) {
    std::stringstream message;
    message << std::endl;
    message << "Optimal solution found.";
    message << std::endl;
    message << "Run status : " << HighsStatusToString(HighsStatus::Optimal)
            << std::endl;
    message << "Objective  : " << std::scientific << tree.getBestObjective()
            << std::endl;
    message << "Time       : " << std::fixed << std::setprecision(3)
            << mip_solve_final_time - mip_solve_initial_time << std::endl;
    message << std::endl;

    HighsPrintMessage(ML_MINIMAL, message.str().c_str());
  } else {
    HighsPrintMessage(ML_ALWAYS, "No feasible solution found.\n");
  }

  return HighsStatus::OK;
}

HighsStatus Highs::solveNode(Node& node) {
  // Apply column bounds from node to LP.
  const bool check_call = false;
  const bool call_changeColsBounds = true;
  if (call_changeColsBounds) {
    changeColsBounds(0, lp_.numCol_, &node.col_lower_bound[0],
                     &node.col_upper_bound[0]);
  } else {
    // Change the LP directly and invalidate the simplex information
    lp_.colLower_ = node.col_lower_bound;
    lp_.colUpper_ = node.col_upper_bound;
    hmos_[0].simplex_lp_status_.valid = false;
  }

  // Call warm start.
  //  HighsStatus status = run();
  // call works but simply calling run() should be enough and will call hot
  // start in the same way as a user would call it from the outside

  int iteration_count0;
  int iteration_count1;
  int solve0_iteration_count;
  int solve1_iteration_count;
  double solve0_objective_value;
  double solve1_objective_value;
  int solve0_status;
  int solve1_status;

  iteration_count0 = hmos_[0].simplex_info_.iteration_count;

  HighsStatus status = solveModelSimplex(hmos_[0]);
  simplex_has_run_ = true;

  iteration_count1 = hmos_[0].simplex_info_.iteration_count;
  solve0_iteration_count = iteration_count1 - iteration_count0;
  solve0_objective_value = hmos_[0].simplex_info_.dual_objective_value;
  solve0_status = (int)status;
  printf("Solve0: Obj = %12g; Iter =%6d; Status =%2d\n", solve0_objective_value,
         solve0_iteration_count, solve0_status);

  if (check_call) {
    // Generate a fresh model object for the LP at this node
    hmos_[0].simplex_lp_status_.has_basis = false;
    hmos_[0].basis_.valid_ = false;
    iteration_count0 = hmos_[0].simplex_info_.iteration_count;
    HighsStatus status = solveModelSimplex(hmos_[0]);
    iteration_count1 = hmos_[0].simplex_info_.iteration_count;
    solve1_iteration_count = iteration_count1 - iteration_count0;
    solve1_objective_value = hmos_[0].simplex_info_.dual_objective_value;
    solve1_status = (int)status;
    printf("Solve1: Obj = %12g; Iter =%6d; Status =%2d\n",
           solve1_objective_value, solve1_iteration_count, solve1_status);
    double rlv_objective_value_difference =
        fabs(solve1_objective_value - solve0_objective_value) /
        max(1.0, fabs(solve1_objective_value));
    if (solve0_status != solve1_status) {
      // Look for unequal status
      printf(
          "!! NodeSolveInequality: Status difference: Status0=%2d; Status1=%2d "
          "!!\n",
          solve0_status, solve1_status);
    } else if (solve0_status != (int)HighsStatus::Infeasible) {
      // Unless infeasible, look for unequal objective
      if (rlv_objective_value_difference > 1e-12)
        printf(
            "!! NodeSolveInequality: Relative objective difference = %12g !!\n",
            rlv_objective_value_difference);
    }
  }

  // Set solution.
  if (status == HighsStatus::Optimal) {
    node.primal_solution = hmos_[0].solution_.col_value;
    node.objective_value = hmos_[0].simplex_info_.dual_objective_value;
  }

  // JAJH(8519) Need to understand why simplex_has_run_ is false for lp_
  // changeColsBounds(0, lp_.numCol_, &node.col_lower_bound[0],
  // &node.col_upper_bound[0]);

  // Solve with a new hmo (replace with code above)
  // initializeLp(lp_);
  // lp_.colLower_ = node.col_lower_bound;
  // lp_.colUpper_ = node.col_upper_bound;

  // HighsStatus status = solveModelSimplex(hmos_[0]);

  // // Set solution.
  // if (status == HighsStatus::Optimal) {
  //   node.primal_solution = hmos_[0].solution_.col_value;
  //   node.objective_value = hmos_[0].simplex_info_.dual_objective_value;
  // }

  return status;
}

HighsStatus Highs::solveRootNode(Node& root) {
  // No presolve for the moment.
  options_.messageLevel = ML_NONE;
  // HighsStatus status = run();
  // call works but simply calling run() should be enough.
  HighsStatus status = solveModelSimplex(hmos_[0]);
  simplex_has_run_ = true;

  if (status == HighsStatus::Optimal) {
    root.primal_solution = hmos_[0].solution_.col_value;
    root.objective_value = hmos_[0].simplex_info_.dual_objective_value;
  }

  return status;
}

void Highs::underDevelopmentLogMessage(const string method_name) {
  HighsLogMessage(HighsMessageType::WARNING, "Method %s is still under development and behaviour may be unpredictable", method_name.c_str());
}

