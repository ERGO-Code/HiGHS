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
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsStatus.h"
#include "presolve/Presolve.h"
#include "presolve/FindFeasibility.h"
#include "simplex/HApp.h"
#include "simplex/HighsSimplexInterface.h"

// until add_row_.. functions are moved to HighsLpUtils.h
#include "simplex/HSimplex.h"
HighsStatus Highs::initializeLp(const HighsLp &lp) {
  // todo:(julian) add code to check that LP is valid.
  lp_ = lp;

  // hmos_[0] is the HighsModelObject corresponding to the original LP
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  simplex_has_run_ = false;
  return HighsStatus::OK;
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run() {

  //  HighsSetMessagelevel(HighsPrintMessageLevel::ML_ALWAYS); reportLp(lp_, 1);
  bool normalise = true;
  HighsStatus return_status = assessLp(lp_, options_, normalise);
  if (return_status == HighsStatus::Error) return return_status;

  // For the moment runFeasibility as standalone.
  if (options_.find_feasibility) {
    // use when you do something with solution depending on whether we have dualized or not.
    HighsSolution& solution = solution_;

    //options_.messageLevel = HighsPrintMessageLevel::ML_DETAILED;
    //HighsSetIO(options_);

    // Add slacks and make sure a minimization problem is passed to runFeasibility.
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


    if (options_.feasibility_strategy == FeasibilityStrategy::kApproxComponentWise)
      return runFeasibility(lp_, solution_, MinimizationType::kComponentWise);
    else if (options_.feasibility_strategy == FeasibilityStrategy::kApproxExact)
      return runFeasibility(lp_, solution_, MinimizationType::kExact);
    else if (options_.feasibility_strategy == FeasibilityStrategy::kDirectSolve) {
      // Proceed to normal exection of run().
      // If dualize has been called replace LP is replaced with dual in code above.
    }
  }


  // Return immediately if the LP has no columns
  if (!lp_.numCol_) return HighsStatus::LpEmpty;

  // todo: check options.
  HighsSetIO(options_);

  HighsPrintMessage(ML_VERBOSE, "Solving %s", lp_.model_name_.c_str());
  if (options_.mip)
    return runBnb();

  // Running as LP solver: start the HiGHS clock unless it's already running
  bool run_highs_clock_already_running = timer_.runningRunHighsClock();
  if (!run_highs_clock_already_running) timer_.startRunHighsClock();
  double tt0 = timer_.readRunHighsClock();
  // todo: make sure it should remain Init for calls of run() after
  // simplex_has_run_ is valid.
  HighsStatus solve_status = HighsStatus::Init;
  // Define identifiers to refer to the HMO of the original LP
  // (0) and the HMO created when using presolve (1)
  const int original_hmo = 0;
  const int presolve_hmo = 1;
  // Keep track of the hmo that is the most recently solved. By default it's the original LP
  int solved_hmo = original_hmo;

  // Initial solve. Presolve, choose solver (simplex, ipx), postsolve.
  if (!simplex_has_run_) {
    // Presolve. runPresolve handles the level of presolving (0 = don't presolve).
    timer_.start(timer_.presolve_clock);
    PresolveInfo presolve_info(options_.presolve_option, lp_);
    HighsPresolveStatus presolve_status = runPresolve(presolve_info);
    timer_.stop(timer_.presolve_clock);
    
    // Run solver.
    switch (presolve_status) {
    case HighsPresolveStatus::NotPresolved: {
      solve_status = runSolver(hmos_[solved_hmo]);
      break;
    }
    case HighsPresolveStatus::NotReduced: {
      printf("Problem not reduced\n");
      solve_status = runSolver(hmos_[solved_hmo]);
      break;
    }
    case HighsPresolveStatus::Reduced: {
      HighsLp& reduced_lp = presolve_info.getReducedProblem();
      // Add reduced lp object to vector of HighsModelObject,
      // so the last one in lp_ is the presolved one.
      hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
      solved_hmo = presolve_hmo;
      solve_status = runSolver(hmos_[solved_hmo]);
      break;
    }
    case HighsPresolveStatus::ReducedToEmpty: {
      // Proceed to postsolve.
      break;
    }
    case HighsPresolveStatus::Infeasible:
    case HighsPresolveStatus::Unbounded: {
      HighsStatus result = (presolve_status == HighsPresolveStatus::Infeasible) ?
	HighsStatus::Infeasible : HighsStatus::Unbounded;
      HighsPrintMessage(ML_ALWAYS, "Problem status detected on presolve: %s\n",
			HighsStatusToString(result).c_str());
      
      // Report this way for the moment. May modify after merge with OSIinterface
      // branch which has new way of setting up a HighsModelObject and can support
      // multiple calls to run().
      // Stop and read the HiGHS clock, then work out time for this call
      if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
      double tt1 = timer_.readRunHighsClock();
      
      std::stringstream message_not_opt;
      message_not_opt << std::endl;
      message_not_opt << "Run status : " << HighsStatusToString(result)
		      << std::endl;
      message_not_opt << "Time       : " << std::fixed << std::setprecision(3)
		      << tt1-tt0 << std::endl;
      
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
    bool run_postsolve = false;//true;
    if (run_postsolve) {
      // Postsolve. Does nothing if there were no reductions during presolve.
      if (solve_status == HighsStatus::Optimal) {
	if (presolve_status == HighsPresolveStatus::Reduced) {
	  // If presolve is nontrivial, extract the optimal solution
	  // and basis for the presolved problem in order to generate
	  // the solution and basis for postsolve to use to generate a
	  // solution(?) and basis that is, hopefully, optimal. This is
	  // confirmed or corrected by hot-starting the simplex solver
	  presolve_info.reduced_solution_ = hmos_[solved_hmo].solution_;
	  presolve_info.presolve_[0].setBasisInfo(hmos_[solved_hmo].simplex_basis_.basicIndex_,
						  hmos_[solved_hmo].simplex_basis_.nonbasicFlag_,
						  hmos_[solved_hmo].simplex_basis_.nonbasicMove_);
	  // Run postsolve
	  timer_.start(timer_.postsolve_clock);
	  HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
	  timer_.stop(timer_.postsolve_clock);
	  if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
	    HighsPrintMessage(ML_VERBOSE, "Postsolve finished.");
	    // Set solution(?) and basis to hot-start the simplex solver
	    // for the original_hmo
	    hmos_[original_hmo].simplex_basis_.valid_ = true;
	    hmos_[original_hmo].simplex_basis_.basicIndex_ = presolve_info.presolve_[0].getBasisIndex();
	    hmos_[original_hmo].simplex_basis_.nonbasicFlag_ = presolve_info.presolve_[0].getNonbasicFlag();
	    hmos_[original_hmo].simplex_basis_.nonbasicMove_ = presolve_info.presolve_[0].getNonbasicMove();
	    options_.clean_up = true;
	    // Now hot-start the simplex solver for the original_hmo
	    solved_hmo = original_hmo;
	    solve_status = runSolver(hmos_[solved_hmo]);
	  }
	}
      }
    } else {
      // Hack to get data for reporting when bypassing postsolve
      if (solved_hmo == presolve_hmo) {
	hmos_[original_hmo].simplex_info_.iteration_count =
	  hmos_[solved_hmo].simplex_info_.iteration_count;
	hmos_[original_hmo].simplex_info_.dualObjectiveValue =
	  hmos_[solved_hmo].simplex_info_.dualObjectiveValue;
      }
    }
  } else {
    // The problem has been solved before so we ignore presolve/postsolve/ipx.
    solved_hmo = original_hmo;
    solve_status = runSolver(hmos_[solved_hmo]);
  }
  // else if (reduced problem failed to solve) {
  //   todo: handle case when presolved problem failed to solve. Try to solve again
  //   with no presolve.
  // }
  
  assert(hmos_.size() > 0);
  solution_ = hmos_[original_hmo].solution_;
  basis_ = hmos_[original_hmo].basis_;
  
  // Report times
  if (hmos_[original_hmo].reportModelOperationsClock) {
    std::vector<int> clockList{timer_.presolve_clock, timer_.scale_clock,
                               timer_.crash_clock, timer_.solve_clock,
                               timer_.postsolve_clock};
    timer_.report("ModelOperations", clockList);
  }
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
  double tt1 = timer_.readRunHighsClock();

  std::stringstream message;
  message << std::endl;
  message << "Run status : " << HighsStatusToString(solve_status)
          << std::endl;
  message << "Iterations : " << hmos_[original_hmo].simplex_info_.iteration_count
          << std::endl;

  if (solve_status == HighsStatus::Optimal)
    message << "Objective  : " << std::scientific
            << hmos_[original_hmo].simplex_info_.dualObjectiveValue << std::endl;

  message << "Time       : " << std::fixed << std::setprecision(3)
          << tt1-tt0 << std::endl;

  message << std::endl;

  HighsPrintMessage(ML_MINIMAL, message.str().c_str());

  return solve_status;
}

const HighsLp &Highs::getLp() const {
  return lp_;
}

const HighsSolution &Highs::getSolution() const {
  return solution_;
}

const HighsBasis &Highs::getBasis() const {
  return basis_;
}

double Highs::getObjectiveValue() const {
  if (hmos_.size() > 0) {
    return hmos_[0].simplex_info_.dualObjectiveValue;
  } else {
    // todo: ipx case
    // todo: error/warning message
  }
  return 0;
}

const int Highs::getIterationCount() const {
  if (hmos_.size() == 0) return 0;
  return hmos_[0].simplex_info_.iteration_count;
}

HighsStatus Highs::setSolution(const HighsSolution &solution) {
  // Check if solution is valid.
  assert(solution_.col_value.size() != 0 ||
         solution_.col_value.size() != lp_.numCol_);
  assert(solution.col_dual.size() == 0 ||
         solution.col_dual.size() == lp_.numCol_);
  assert(solution.row_dual.size() == 0 ||
         solution.row_dual.size() == lp_.numRow_);

  if (solution.col_value.size()) solution_.col_value = solution.col_value;
  if (solution.col_dual.size()) solution_.col_dual = solution.col_dual;
  if (solution.row_dual.size()) solution_.row_dual = solution.row_dual;

  HighsStatus result_values = HighsStatus::NotSet;
  HighsStatus result_duals = HighsStatus::NotSet;

  if (solution.col_value.size() > 0)
    result_values = calculateRowValues(lp_, solution_);
  if (solution.row_dual.size() > 0)
    result_duals = calculateColDuals(lp_, solution_);

  if (result_values == HighsStatus::Error ||
      result_duals == HighsStatus::Error);
    return HighsStatus::Error;

  return HighsStatus::OK;
}

HighsStatus Highs::setBasis(const HighsBasis& basis) {
  basis_ = basis;
  return HighsStatus::OK;
}

void Highs::reportSolution() {
  reportModelBoundSol(true, lp_.numCol_,
		      lp_.colLower_, lp_.colUpper_,
		      lp_.col_names_,
		      solution_.col_value, solution_.col_dual,
		      basis_.col_status);
  reportModelBoundSol(false, lp_.numRow_,
		      lp_.rowLower_, lp_.rowUpper_,
		      lp_.row_names_,
		      solution_.row_value, solution_.row_dual,
		      basis_.row_status);
}

bool Highs::addRow(const double lower_bound, const double upper_bound,
                   const int num_new_nz, const int *indices, const double *values) {
  int starts = 0;
  return addRows(1, &lower_bound, &upper_bound, 
                 num_new_nz, &starts, indices, values);
}

bool Highs::addRows(const int num_new_row,
		    const double *lower_bounds, const double *upper_bounds, 
                    const int num_new_nz,
		    const int *starts, const int *indices, const double *values) {
  HighsStatus return_status = HighsStatus::NotSet;
  // if simplex has not solved already
  if (!simplex_has_run_) {
    return_status = add_lp_rows(lp_, num_new_row, lower_bounds, upper_bounds,
				num_new_nz, starts, indices, values, options_);
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    // todo: change to take int return value
    return_status = interface.util_add_rows(num_new_row, lower_bounds, upper_bounds,
					    num_new_nz, starts, indices, values);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::addCol(const double cost, const double lower_bound, const double upper_bound,
		   const int num_new_nz, const int *indices, const double *values) {
  int starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, 
                 num_new_nz, &starts, indices, values);
}


bool Highs::addCols(const int num_new_col, 
                    const double *costs, const double *lower_bounds, const double *upper_bounds,
                    const int num_new_nz,
		    const int *starts, const int *indices, const double *values) {
  HighsStatus return_status = HighsStatus::NotSet;
  // if simplex has not solved already
  if (!simplex_has_run_) {
    return_status = add_lp_cols(lp_, num_new_col,
				costs, lower_bounds, upper_bounds,
				num_new_nz, starts, indices, values, options_);
  } else
  {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    // todo: change to take int return value
    return_status = interface.util_add_cols(num_new_col,
					    costs, lower_bounds, upper_bounds,
					    num_new_nz, starts, indices, values);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::changeObjectiveSense(const int sense) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    this->lp_.sense_ = sense;
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.change_ObjSense(sense);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::changeColCost(const int col, const double cost) {
  return changeColsCost(1, &col, &cost);
}

bool Highs::changeColsCost(const int num_set_entries, const int* set, const double* cost) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i=0; i<num_set_entries; i++) {
      this->lp_.colCost_[set[i]] = cost[i];
    }
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.change_costs(num_set_entries, set, cost);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;  
}

bool Highs::changeColsCost(const int* mask, const double* cost) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i=0; i<this->lp_.numCol_; i++) {
      if (mask[i]) this->lp_.colCost_[i] = cost[i];
    }
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.change_costs(mask, cost);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;  
}

bool Highs::changeColBounds(const int col, const double lower, const double upper) {
  return changeColsBounds(1, &col, &lower, &upper);
}

bool Highs::changeColsBounds(const int num_set_entries, const int *set, const double *lower, const double *upper) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i = 0; i < num_set_entries; i++) {
      this->lp_.colLower_[set[i]] = lower[i];
      this->lp_.colUpper_[set[i]] = upper[i];
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeColBounds(num_set_entries, set, lower, upper);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::changeColsBounds(const int from_col, const int to_col, const double *lower, const double *upper) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    int num_changed_bounds = to_col-from_col;
    for (int i = 0; i < num_changed_bounds; i++) {
      this->lp_.colLower_[from_col+i] = lower[i];
      this->lp_.colUpper_[from_col+i] = upper[i];
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeColBounds(from_col, to_col, lower, upper);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::changeColsBounds(const int *mask, const double *lower, const double *upper) {
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
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::changeRowBounds(const int row,
			    const double lower, const double upper) {
  return changeRowsBounds(1, &row, &lower, &upper);
}

bool Highs::changeRowsBounds(const int num_set_entries,
			     const int *set, const double *lower, const double *upper) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    for (int i = 0; i < num_set_entries; i++) {
      this->lp_.rowLower_[set[i]] = lower[i];
      this->lp_.rowUpper_[set[i]] = upper[i];
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.changeRowBounds(num_set_entries, set, lower, upper);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::changeRowsBounds(const int *mask, const double *lower, const double *upper) {
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
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::getCols(const int from_col, const int to_col,
		    int &num_col, double *costs, double *lower, double *upper,
		    int &num_nz, int *start, int *index, double *value) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getCols(from_col, to_col,
				      num_col, costs, lower, upper,
				      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::getCols(const int n, const int *set,
		    int &num_col, double *costs, double *lower, double *upper,
		    int &num_nz, int *start, int *index, double *value) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getCols(n, set,
				      num_col, costs, lower, upper,
				      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::getCols(const int *col_mask,
		    int &num_col, double *costs, double *lower, double *upper,
		    int &num_nz, int *start, int *index, double *value) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getCols(col_mask,
		      num_col, costs, lower, upper,
		      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::getRows(const int from_row, const int to_row,
		    int &num_row, double *lower, double *upper,
		    int &num_nz, int *start, int *index, double *value) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getRows(from_row, to_row,
		      num_row, lower, upper,
		      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::getRows(const int num_set_entries, const int *set,
		    int &num_row, double *lower, double *upper,
		    int &num_nz, int *start, int *index, double *value) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getRows(num_set_entries, set,
		      num_row, lower, upper,
		      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::getRows(const int *mask,
		    int &num_row, double *lower, double *upper,
		    int &num_nz, int *start, int *index, double *value) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    return_status = interface.getRows(mask,
		      num_row, lower, upper,
		      num_nz, start, index, value);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::deleteCols(const int from_col, const int to_col) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.delete_cols(from_col, to_col);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::deleteCols(const int num_set_entries, const int *set) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.delete_cols(num_set_entries, set);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::deleteCols(int *mask) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.delete_cols(mask);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::deleteRows(const int from_row, const int to_row) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.delete_rows(from_row, to_row);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::deleteRows(const int num_set_entries, const int *set) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.delete_rows(num_set_entries, set);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::deleteRows(int *mask) {
  HighsStatus return_status = HighsStatus::NotSet;
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    return_status = interface.delete_rows(mask);
  }
  if (return_status == HighsStatus::Error ||
      return_status == HighsStatus::NotSet) return false;
  return true;
}

bool Highs::writeMPS(const char* filename) {
  return writeLpAsMPS(filename, lp_);
}

// Private methods
HighsPresolveStatus Highs::runPresolve(PresolveInfo &info) {
  if (options_.presolve_option != PresolveOption::ON)
    return HighsPresolveStatus::NotPresolved;

  if (info.lp_ == nullptr) return HighsPresolveStatus::NullError;

  if (info.presolve_.size() == 0) return HighsPresolveStatus::NotReduced;

  info.presolve_[0].load(*(info.lp_));

  // Initialize a new presolve class instance for the LP given in presolve info
  return info.presolve_[0].presolve();
}

HighsPostsolveStatus Highs::runPostsolve(PresolveInfo &info) {
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
HighsStatus Highs::runSolver(HighsModelObject &model) {
  // trim bounds to get OSI unit test to pass: delete when checkLp
  // disappears.
  HighsLp& lp = model.lp_;
  for (int i = 0; i < lp.numRow_; i++) {
      if (lp.rowLower_[i] < -HIGHS_CONST_INF)
        lp.rowLower_[i] = -HIGHS_CONST_INF;
      if (lp.rowUpper_[i] > HIGHS_CONST_INF)
        lp.rowUpper_[i] = HIGHS_CONST_INF;
    }

  for (int j = 0; j < lp.numCol_; j++) {
      if (lp.colLower_[j] < -HIGHS_CONST_INF)
        lp.colLower_[j] = -HIGHS_CONST_INF;
      if (lp.colUpper_[j] > HIGHS_CONST_INF)
        lp.colUpper_[j] = HIGHS_CONST_INF;
    }

  assert(checkLp(model.lp_) == HighsStatus::OK);
  //  assert(assessLp(lp, options) == HighsStatus::OK);
  
  HighsStatus status = HighsStatus::Init;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = runSimplexSolver(options_, model);
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
  double tt0 = timer_.readRunHighsClock();

  // Start tree by making root node.
   std::unique_ptr<Node> root =std::unique_ptr<Node>(new Node(-1, 0, 0));

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

    if (status == HighsStatus::Infeasible)
      continue;

    options_.messageLevel = message_level;
    tree.branch(node);
  }
  
  // Stop and read the HiGHS clock, then work out time for this call
  if (!run_highs_clock_already_running) timer_.stopRunHighsClock();
  double tt1 = timer_.readRunHighsClock();

  if (tree.getBestSolution().size() > 0) {
    std::stringstream message;
    message << std::endl;
    message << "Optimal solution found.";
    message << std::endl;
    message << "Run status : " << HighsStatusToString(HighsStatus::Optimal)
            << std::endl;
    message << "Objective  : " << std::scientific
            << tree.getBestObjective() << std::endl;
    message << "Time       : " << std::fixed << std::setprecision(3)
            << tt1-tt0 << std::endl;
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
    changeColsBounds(0, lp_.numCol_, &node.col_lower_bound[0], &node.col_upper_bound[0]);
  } else {
    // Change the LP directly and ivalidate the simplex information
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

  HighsStatus status = runSimplexSolver(options_, hmos_[0]);
  simplex_has_run_ = true;

  iteration_count1 = hmos_[0].simplex_info_.iteration_count;
  solve0_iteration_count = iteration_count1 - iteration_count0;
  solve0_objective_value = hmos_[0].simplex_info_.dualObjectiveValue;
  solve0_status = (int)status;
  printf("Solve0: Obj = %12g; Iter =%6d; Status =%2d\n", solve0_objective_value, solve0_iteration_count, solve0_status);

  if (check_call) {
    // Generate a fresh model object for the LP at this node
    hmos_[0].simplex_lp_status_.valid = false;
    hmos_[0].basis_.valid_ = false;
    hmos_[0].simplex_basis_.valid_ = false;
    iteration_count1 = hmos_[0].simplex_info_.iteration_count;
    HighsStatus status = runSimplexSolver(options_, hmos_[0]);
    iteration_count1 = hmos_[0].simplex_info_.iteration_count;
    solve1_iteration_count = iteration_count1 - iteration_count0;
    solve1_iteration_count = iteration_count1 - iteration_count0;
    solve1_objective_value = hmos_[0].simplex_info_.dualObjectiveValue;
    solve1_status = (int)status;
    printf("Solve1: Obj = %12g; Iter =%6d; Status =%2d\n", solve1_objective_value, solve1_iteration_count, solve1_status);
    double rlv_objective_value_difference = fabs(solve1_objective_value - solve0_objective_value)/max(1.0, fabs(solve1_objective_value));
    if (solve0_status != solve1_status) {
      // Look for unequal status
      printf("!! NodeSolveInequality: Status difference: Status0=%2d; Status1=%2d !!\n", solve0_status, solve1_status);
    } else if (solve0_status != (int)HighsStatus::Infeasible) {
      // Unless infeasible, look for unequal objective
      if (rlv_objective_value_difference > 1e-12)
	printf("!! NodeSolveInequality: Relative objective difference = %12g !!\n", rlv_objective_value_difference);
    }
  }

  // Set solution.
  if (status == HighsStatus::Optimal) {
    node.primal_solution = hmos_[0].solution_.col_value;
    node.objective_value = hmos_[0].simplex_info_.dualObjectiveValue;
  }

  // JAJH(8519) Need to understand why simplex_has_run_ is false for lp_

  // Solve with a new hmo (replace with code above)
  // initializeLp(lp_);
  // lp_.colLower_ = node.col_lower_bound;
  // lp_.colUpper_ = node.col_upper_bound;

  // HighsStatus status = runSimplexSolver(options_, hmos_[0]);

  // // Set solution.
  // if (status == HighsStatus::Optimal) {
  //   node.primal_solution = hmos_[0].solution_.col_value;
  //   node.objective_value = hmos_[0].simplex_info_.dualObjectiveValue;
  // }

  return status;
}

HighsStatus Highs::solveRootNode(Node& root) {
  // No presolve for the moment.
  options_.messageLevel = ML_NONE;
  //HighsStatus status = run();
  // call works but simply calling run() should be enough.
  HighsStatus status = runSimplexSolver(options_, hmos_[0]);
  simplex_has_run_ = true;
  
  if (status == HighsStatus::Optimal) {
    root.primal_solution = hmos_[0].solution_.col_value;
    root.objective_value = hmos_[0].simplex_info_.dualObjectiveValue;
  }

  return status;
}
