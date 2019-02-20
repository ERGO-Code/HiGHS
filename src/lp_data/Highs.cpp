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
#include "simplex/HApp.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelObject.h"
#include "simplex/HighsSimplexInterface.h"
#include "lp_data/HighsStatus.h"
#include "presolve/Presolve.h"

// until add_row_.. functions are moved to HighsLpUtils.h
#include "simplex/HSimplex.h"
HighsStatus Highs::initializeLp(HighsLp &lp) {
  // todo:(julian) add code to check that LP is valid.
  lp_ = lp;
  return HighsStatus::OK;
}

int Highs::HighsAddVariable(double obj, double lo, double hi) {
  if (this->runSuccessful) {
    HighsSimplexInterface simplex_interface(this->hmos_[0]);
    simplex_interface.util_add_cols(1, &obj, &lo, &hi, 0, NULL, NULL, NULL);
    return 0; //TODO
    
  } else {
    // build initial model using HighsModelBuilder
    HighsVar* newVariable;
    this->builder.HighsCreateVar(NULL, lo, hi, obj, HighsVarType::CONT, &newVariable);
    return this->builder.getNumberOfVariables();
  }
}

// Checks the options calls presolve and postsolve if needed. Solvers are called
// with runSolver(..)
HighsStatus Highs::run() {
  HighsPrintMessage(HighsMessageType::INFO, "Solving %s", lp_.model_name_.c_str());
  // Not solved before, so create an instance of HighsModelObject.
  hmos_.push_back(HighsModelObject(lp_, timer));

  // Options for HighsPrintMessage and HighsPrintMessage
  options_.logfile = stdout;//fopen("HiGHS.log", "w");
  options_.output = stdout;
  options_.messageLevel = ML_DEFAULT;
  HighsSetIO(options_);

  //Define clocks
  HighsTimer &timer = hmos_[0].timer_;
  timer.startRunHighsClock();

  // Presolve. runPresolve handles the level of presolving (0 = don't presolve).
  timer.start(timer.presolve_clock);
  PresolveInfo presolve_info(options_.presolve_option, lp_);
  HighsPresolveStatus presolve_status = runPresolve(presolve_info);
  timer.stop(timer.presolve_clock);
 
  // Run solver.
  HighsStatus solve_status = HighsStatus::Init;
  switch (presolve_status) {
    case HighsPresolveStatus::NotReduced: {
      solve_status = runSolver(hmos_[0]);
      break;
    }
    case HighsPresolveStatus::Reduced: {
      HighsLp& reduced_lp = presolve_info.getReducedProblem();
      // Add reduced lp object to vector of HighsModelObject,
      // so the last one in lp_ is the presolved one.
      hmos_.push_back(HighsModelObject(reduced_lp, timer));
      solve_status = runSolver(hmos_[1]);
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
      presolve_info.presolve_[0].setBasisInfo(
          hmos_[1].basis_.basicIndex_, hmos_[1].basis_.nonbasicFlag_,
          hmos_[1].basis_.nonbasicMove_);
    }

    timer.start(timer.postsolve_clock);
    HighsPostsolveStatus postsolve_status = runPostsolve(presolve_info);
    timer.stop(timer.postsolve_clock);
    if (postsolve_status == HighsPostsolveStatus::SolutionRecovered) {
       HighsPrintMessage(HighsMessageType::INFO, "Postsolve finished.");

      // Set solution and basis info for simplex clean up.
      // Original LP is in lp_[0] so we set the basis information there.
      hmos_[0].basis_.basicIndex_ =
          presolve_info.presolve_[0].getBasisIndex();
      hmos_[0].basis_.nonbasicFlag_ =
          presolve_info.presolve_[0].getNonbasicFlag();
      hmos_[0].basis_.nonbasicMove_ =
          presolve_info.presolve_[0].getNonbasicMove();

      options_.clean_up = true;

      solve_status = runSolver(hmos_[0]);
    }
  }
  
  assert(hmos_.size() > 0);
  int last = hmos_.size() - 1;
  solution_ = hmos_[last].solution_;

  HighsSimplexInterface simplex_interface(hmos_[0]);
  if (solve_status != HighsStatus::Optimal) {
    if (solve_status == HighsStatus::Infeasible ||
        solve_status == HighsStatus::Unbounded) {
      if (options_.presolve_option == PresolveOption::ON) {
        std::stringstream ss;
        ss << "Reduced problem status: " << HighsStatusToString(solve_status) << ".";
        HighsPrintMessage(HighsMessageType::ERROR,ss.str().c_str());
        // todo: handle case. Try to solve again with no presolve?
        return HighsStatus::NotImplemented;
      } else {
        std::cout << "Solver terminated with a non-optimal status: "
                  << HighsStatusToString(solve_status) << std::endl;
        simplex_interface.report_simplex_outcome("Run");
      }
    }
  } else {
    // Report in old way so tests pass.
    simplex_interface.report_simplex_outcome("Run");
  }

  if (hmos_[0].reportModelOperationsClock) {
    // Report times
    std::vector<int> clockList{timer.presolve_clock, timer.scale_clock, timer.crash_clock, timer.solve_clock, timer.postsolve_clock};
    timer.report("ModelOperations", clockList);
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
          model.getPrStatus(), highs_model.lp_->model_name_.c_str(), Presolve_ArgV,
          Crash_ArgV, EdWt_ArgV, Price_ArgV, numRow, numCol, setupTime,
          presolve1Time, crashTime, crossoverTime, presolve2Time, solveTime,
          postsolveTime, simplex_info_.dualObjectiveValue, simplex_info_.iteration_count,
          model.totalTime, solver.n_wg_DSE_wt);
      cout << flush;
    }
*/
#endif

  timer.stopRunHighsClock();

  return HighsStatus::OK;
}

HighsPresolveStatus Highs::runPresolve(PresolveInfo& info) {
  if (options_.presolve_option != PresolveOption::ON)
    return HighsPresolveStatus::NotReduced;

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
HighsStatus Highs::runSolver(HighsModelObject& model) {
  assert(checkLp(model.lp_) == HighsStatus::OK);

  HighsStatus status = HighsStatus::Init;
#ifndef IPX
  // HiGHS
  // todo: Without the presolve part, so will be
  //     = solve_simplex(options, reduced_lp, reduced_solution)
  status = runSimplexSolver(options_, model);
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
            const int num_new_nz,
            const int *columns, const double *values,
            const bool force) {
  int row_starts = 0;
  return addRows(1, &lower_bound, &upper_bound, &row_starts,
                 num_new_nz, columns, values);
}

bool Highs::addRows(const int num_new_rows,
             const double *lower_bounds, const double *upper_bounds,
             const int *row_starts,
             const int num_new_nz,
             const int *columns, const double *values,
             const bool force) {
  // if simplex has not solved already
  if (hmos_.size() == 0) {
    // add_lp_rows(lp, ..., options);

    // the above currently being
    // add_rows_to_lp_vectors
    // add_rows_to_lp_matrix
  }
  // else (if simplex has solved already)
  {
    assert(hmos_.size() > 0);
    int last = hmos_.size() - 1;
    HighsSimplexInterface interface(hmos_[last]);

    // todo: change to take int return value
    interface.util_add_rows(num_new_rows, lower_bounds, upper_bounds,
                            num_new_nz, row_starts, columns, values);
  }

  return 0;
}

bool Highs::addCol(const double lower_bound, const double upper_bound,
            const int num_new_nz,
            const int *rows, const double *values,
            const bool force) {
  int col_starts = 0;
  return addCols(1, &lower_bound, &upper_bound, &col_starts,
                 num_new_nz, rows, values);
}

bool Highs::addCols(const int num_new_cols,
             const double *lower_bounds, const double *upper_bounds,
             const int *col_starts,
             const int num_new_nz,
             const int *rows, const double *values,
             const bool force) {
  // if simplex has not solved already
  if (hmos_.size() == 0) {
    // add_lp_cols(lp, ..., options);  
    
    // the above currently being
    // add_cols_to_lp_vectors
    // add_cols_to_lp_matrix
  }
  // else (if simplex has solved already)
  {
    assert(hmos_.size() > 0);
    int last = hmos_.size() - 1;
    HighsSimplexInterface interface(hmos_[last]);

    // todo: change to take int return value
    interface.util_add_rows(num_new_cols, lower_bounds, upper_bounds,
                            num_new_nz, col_starts, rows, values);
  }
  
  return 0;
}

  double Highs::getObjectiveValue() {
    if (hmos_.size() > 0) {
      int last = hmos_.size() - 1;
      return hmos_[last].simplex_info_.dualObjectiveValue;
    } else {
      // todo: ipx case
      // todo: error/warning message
    }
    return 0;
  }