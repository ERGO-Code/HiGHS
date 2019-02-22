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

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsStatus.h"
#include "presolve/Presolve.h"
#include "simplex/HApp.h"
#include "simplex/HighsSimplexInterface.h"

// until add_row_.. functions are moved to HighsLpUtils.h
#include "simplex/HSimplex.h"
HighsStatus Highs::initializeLp(const HighsLp &lp) {
  // todo:(julian) add code to check that LP is valid.
  lp_ = lp;

  // For the moment hmos_[0] is the HighsModelObject corresponding to the
  // original lp.
  hmos_.clear();
  hmos_.push_back(HighsModelObject(lp_, options_, timer_));
  simplex_has_run_ = false;
  return HighsStatus::OK;
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run() {
  // todo: check if lp is empty.

  // todo: check options.
  HighsSetIO(options_);

  HighsPrintMessage(HighsMessageType::INFO, "Solving %s",
                    lp_.model_name_.c_str());

  timer_.startRunHighsClock();
  // todo: make sure it should remain Init for calls of run() after
  // simplex_has_run_ is valid.
  HighsStatus solve_status = HighsStatus::Init;

  // Initial solve. Presolve, choose solver (simplex, ipx), postsolve.
  if (!simplex_has_run_) {
    // Presolve. runPresolve handles the level of presolving (0 = don't presolve).
    timer_.start(timer_.presolve_clock);
    PresolveInfo presolve_info(options_.presolve_option, lp_);
    HighsPresolveStatus presolve_status = runPresolve(presolve_info);
    timer_.stop(timer_.presolve_clock);

    // Run solver.
    switch (presolve_status) {
      case HighsPresolveStatus::NotReduced: {
        solve_status = runSolver(hmos_[0]);
        break;
      }
      case HighsPresolveStatus::Reduced: {
        HighsLp &reduced_lp = presolve_info.getReducedProblem();
        // Add reduced lp object to vector of HighsModelObject,
        // so the last one in lp_ is the presolved one.
        hmos_.push_back(HighsModelObject(reduced_lp, options_, timer_));
        solve_status = runSolver(hmos_[1]);
        break;
      }
      case HighsPresolveStatus::ReducedToEmpty: {
        // Proceed to postsolve.
        break;
      }
      case HighsPresolveStatus::Infeasible:
      case HighsPresolveStatus::Unbounded: {
        HighsStatus result = (presolve_status == HighsPresolveStatus::Infeasible)
                                ? HighsStatus::Infeasible
                                : HighsStatus::Unbounded;
        HighsPrintMessage(ML_ALWAYS, "Problem status detected on presolve: %s\n",
                          HighsStatusToString(result).c_str());
        // for tests
        HighsPrintMessage(ML_ALWAYS, "Run: NOT-OPT\n");
        return result;
      }
      default: {
        // case HighsPresolveStatus::Error
        HighsPrintMessage(HighsMessageType::ERROR, "Presolve failed.");
        return HighsStatus::PresolveError;
      }
    }

    // Postsolve. Does nothing if there were no reductions during presolve.
    if (solve_status == HighsStatus::Optimal) {
      if (presolve_status == HighsPresolveStatus::Reduced) {
        presolve_info.reduced_solution_ = hmos_[1].solution_;
        presolve_info.presolve_[0].setBasisInfo(hmos_[1].basis_.basicIndex_,
                                                hmos_[1].basis_.nonbasicFlag_,
                                                hmos_[1].basis_.nonbasicMove_);
      }

      timer_.start(timer_.postsolve_clock);
      HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
      timer_.stop(timer_.postsolve_clock);
      if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
        HighsPrintMessage(HighsMessageType::INFO, "Postsolve finished.");

        // Set solution and basis info for simplex clean up.
        // Original LP is in lp_[0] so we set the basis information there.
        hmos_[0].basis_.basicIndex_ = presolve_info.presolve_[0].getBasisIndex();
        hmos_[0].basis_.nonbasicFlag_ =
            presolve_info.presolve_[0].getNonbasicFlag();
        hmos_[0].basis_.nonbasicMove_ =
            presolve_info.presolve_[0].getNonbasicMove();

        options_.clean_up = true;

        solve_status = runSolver(hmos_[0]);
      }
    }
  } else {
    // The problem has been solved before so we ignore presolve/postsolve/ipx.
    solve_status = runSolver(hmos_[0]);
  }

  // Report status.
  HighsSimplexInterface simplex_interface(hmos_[0]);
  switch (solve_status) {
    case HighsStatus::Optimal:
    case HighsStatus::ReachedDualObjectiveUpperBound: {
      solution_ = hmos_[0].solution_;
      basis_ = getHighsBasis(hmos_[0].basis_);
      break;
    }
    case HighsStatus::Infeasible:
    case HighsStatus::Unbounded: {
      if (options_.presolve_option == PresolveOption::ON) {
        // todo: handle case. Try to solve again with no presolve.
        HighsPrintMessage(ML_ALWAYS, "Reduced problem status not optiaml: %s\n",
                          HighsStatusToString(solve_status));
        break;
      }
    }
  }

  // Report in old way so tests pass. Change to better output and modify tests.
  simplex_interface.report_simplex_outcome("Run");

  // Report times
  if (hmos_[0].reportModelOperationsClock) {
    std::vector<int> clockList{timer_.presolve_clock, timer_.scale_clock,
                               timer_.crash_clock, timer_.solve_clock,
                               timer_.postsolve_clock};
    timer_.report("ModelOperations", clockList);
  }
#ifdef HiGHSDEV
/* todo: do elsewhere once timing is added.
    bool rpBnchmk = false;
    if (rpBnchmk) {
      int numCol = highs_model.lp_.numCol_;
      int numRow = highs_model.lp_.numRow_;
      printf(
          "\nBnchmkHsol99,hsol,%3d,%16s,Presolve %s,"
          "Crash %s,EdWt %s,Price %s,%d,%d,%10.3f,%10.3f,"
          "%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,"
          "%20.10e,%10d,%10.3f,"
          "%d\n",
          model.getPrStatus(), highs_model.lp_->model_name_.c_str(),
   Presolve_ArgV, Crash_ArgV, EdWt_ArgV, Price_ArgV, numRow, numCol, setupTime,
          presolve1Time, crashTime, crossoverTime, presolve2Time, solveTime,
          postsolveTime, simplex_info_.dualObjectiveValue,
   simplex_info_.iteration_count, model.totalTime, solver.n_wg_DSE_wt); cout <<
   flush;
    }
*/
#endif

  timer_.stopRunHighsClock();

  return HighsStatus::OK;
}

HighsPresolveStatus Highs::runPresolve(PresolveInfo &info) {
  if (options_.presolve_option != PresolveOption::ON)
    return HighsPresolveStatus::NotReduced;

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
  assert(checkLp(model.lp_) == HighsStatus::OK);

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

bool Highs::addRow(const double lower_bound, const double upper_bound,
                   const int num_new_nz, const int *columns,
                   const double *values, const bool force) {
  int row_starts = 0;
  return addRows(1, &lower_bound, &upper_bound, &row_starts, num_new_nz,
                 columns, values);
}

bool Highs::addRows(const int num_new_rows, const double *lower_bounds,
                    const double *upper_bounds, const int *row_starts,
                    const int num_new_nz, const int *columns,
                    const double *values, const bool force) {
  // if simplex has not solved already
  if (!simplex_has_run_) {
    add_lp_rows(lp_, num_new_rows, lower_bounds, upper_bounds, num_new_nz,
                row_starts, columns, values, options_);
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    // todo: change to take int return value
    interface.util_add_rows(num_new_rows, lower_bounds, upper_bounds,
                            num_new_nz, row_starts, columns, values);
  }

  return 0;
}

bool Highs::addCol(const double cost, const double lower_bound,
                   const double upper_bound, const int num_new_nz,
                   const int *rows, const double *values, const bool force) {
  int col_starts = 0;
  return addCols(1, &cost, &lower_bound, &upper_bound, &col_starts, num_new_nz,
                 rows, values);
}

bool Highs::addCols(const int num_new_cols, const double *column_costs,
                    const double *lower_bounds, const double *upper_bounds,
                    const int *col_starts, const int num_new_nz,
                    const int *rows, const double *values, const bool force) {
  // if simplex has not solved already
  if (!simplex_has_run_) {
    add_lp_cols(lp_, num_new_cols, column_costs, lower_bounds, upper_bounds,
                num_new_nz, col_starts, rows, values, options_);
  }
  // else (if simplex has solved already)
  {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    // todo: change to take int return value
    interface.util_add_rows(num_new_cols, lower_bounds, upper_bounds,
                            num_new_nz, col_starts, rows, values);
  }

  return 0;
}

double Highs::getObjectiveValue() const {
  if (hmos_.size() > 0) {
    int last = hmos_.size() - 1;
    return hmos_[last].simplex_info_.dualObjectiveValue;
  } else {
    // todo: ipx case
    // todo: error/warning message
  }
  return 0;
}

const HighsLp &Highs::getLp() const { return lp_; }

const HighsSolution &Highs::getSolution() const { return solution_; }

const HighsBasis_new &Highs::getBasis() const { return basis_; }

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

bool Highs::changeObjectiveSense(int sense) {
  if (!simplex_has_run_) {
    this->lp_.sense_ = sense;
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);
    interface.change_ObjSense(sense);
  }
  return true;
}

bool Highs::changeRowBounds(int index, double lower, double higher) {
  if (!simplex_has_run_) {
    this->lp_.rowLower_[index] = lower;
    this->lp_.rowUpper_[index] = higher;
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.change_row_bounds_set(1, &index, &lower, &higher, true);
  }
  return true;
}

bool Highs::changeColBounds(int index, double lower, double higher) {
  if (!simplex_has_run_) {
    this->lp_.colLower_[index] = lower;
    this->lp_.colUpper_[index] = higher;
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.change_col_bounds_set(1, &index, &lower, &higher, true);
  }
  return true;
}

bool Highs::changeRowsBounds(int n, int *index, double *lower, double *higher) {
  if (!simplex_has_run_) {
    for (int i = 0; i < n; i++) {
      this->lp_.rowLower_[index[i]] = lower[i];
      this->lp_.rowUpper_[index[i]] = higher[i];
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.change_row_bounds_set(n, index, lower, higher, true);
  }
  return true;
}

bool Highs::changeColsBounds(int n, int *index, double *lower, double *higher) {
  if (!simplex_has_run_) {
    for (int i = 0; i < n; i++) {
      this->lp_.colLower_[index[i]] = lower[i];
      this->lp_.colUpper_[index[i]] = higher[i];
    }

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.change_col_bounds_set(n, index, lower, higher, true);
  }
  return true;
}

bool Highs::changeObjCoef(int index, double coef) {
  if (!simplex_has_run_) {
    this->lp_.colCost_[index] = coef;
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.change_costs_set(1, &index, &coef);
  }
  return true;
}

bool Highs::deleteRows(const int n, const int *indices) {
  std::vector<int> rows;
  rows.assign(indices, indices + n);

  if (!simplex_has_run_) {
    // TODO: modify local lp

  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.util_delete_row_set(rows);
  }
  return true;
}

bool Highs::deleteCols(const int n, const int *indices) {
  std::vector<int> cols;
  cols.assign(indices, indices + n);
  if (!simplex_has_run_) {
    // TODO: modify local lp
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.util_delete_col_set(cols);
  }
  return true;
}

bool Highs::changeObjCoefs(int n, int* index, double* coef) {
  if (!simplex_has_run_) {
    for (int i=0; i<n; i++) {
      this->lp_.colCost_[index[i]] = coef[i];
    }
  } else {
    assert(hmos_.size() > 0);
    HighsSimplexInterface interface(hmos_[0]);

    interface.change_costs_set(n, index, coef);
  }
  return true;  
}