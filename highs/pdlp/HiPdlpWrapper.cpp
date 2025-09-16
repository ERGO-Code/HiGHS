/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/HiPdlpWrapper.cpp
 * @brief
 * @author Julian Hall
 */
#include "pdlp/HiPdlpWrapper.h"

#include "pdlp/hipdlp/logger.hpp"
#include "pdlp/hipdlp/pdhg.hpp"
#include "pdlp/hipdlp/restart.hpp"

HighsStatus solveLpHiPdlp(HighsLpSolverObject& solver_object) {
  return solveLpHiPdlp(solver_object.options_, solver_object.timer_,
                       solver_object.lp_, solver_object.basis_,
                       solver_object.solution_, solver_object.model_status_,
                       solver_object.highs_info_, solver_object.callback_);
}

HighsStatus solveLpHiPdlp(const HighsOptions& options, HighsTimer& timer,
                          const HighsLp& lp, HighsBasis& highs_basis,
                          HighsSolution& highs_solution,
                          HighsModelStatus& model_status, HighsInfo& highs_info,
                          HighsCallback& callback) {
  // Indicate that no imprecise solution has (yet) been found
  resetModelStatusAndHighsInfo(model_status, highs_info);

  /*** Order of operations
   * Preprocess with HiPdlp
   * Scale with HiPdlp
   * Solve with HiPdlp
   * Unscale with HiPdlp
   * Postprocess with HiPDLP
   * ***/
  PDLPSolver pdlp;

  // 0. Set up logger and params
  pdlp.setup(options, timer);

  // 1. Pass the LP to be solved
  pdlp.passLp(&lp);

  // 2. Preprocess with HiPdlp
  pdlp.preprocessLp();

  // 3. Scale with HiPdlp
  pdlp.scaleProblem();

  // 4. Solve with HiPdlp
  std::vector<double> x, y;
  pdlp.solve(x, y);

  // 5. Unscale with HiPdlp
  pdlp.unscaleSolution(x, y);

  // 6. Postprocess with HiPDLP
  HighsSolution pdlp_solution;
  PostSolveRetcode retcode = pdlp.postprocess(pdlp_solution);
  if (retcode != PostSolveRetcode::OK) {
    assert(111 == 444);
    return HighsStatus::kError;
  }

  // --- Print Summary ---
  pdlp.logSummary();

  highs_info.pdlp_iteration_count = pdlp.getIterationCount();

  model_status = HighsModelStatus::kUnknown;
  highs_solution.clear();
  highs_basis.valid = false;

  const TerminationStatus termination_status = pdlp.getTerminationCode();
  switch (termination_status) {
    case TerminationStatus::OPTIMAL: {
      model_status = HighsModelStatus::kOptimal;
      break;
    }
    case TerminationStatus::INFEASIBLE: {
      assert(111 == 222);
      model_status = HighsModelStatus::kInfeasible;
      return HighsStatus::kOk;
      break;
    }
    case TerminationStatus::UNBOUNDED: {
      assert(111 == 333);
      model_status = HighsModelStatus::kUnbounded;
      return HighsStatus::kOk;
      break;
    }
    case TerminationStatus::TIMEOUT: {
      // ToDo IterationLimit termination needs to be handled separately
      model_status =
          highs_info.pdlp_iteration_count >= options.pdlp_iteration_limit
              ? HighsModelStatus::kIterationLimit
              : HighsModelStatus::kTimeLimit;
      break;
    }
    case TerminationStatus::WARNING:
    case TerminationStatus::FEASIBLE: {
      assert(111 == 555);
      model_status = HighsModelStatus::kUnknown;
      return HighsStatus::kError;
    }
    default:
      assert(termination_status == TerminationStatus::ERROR);
      return HighsStatus::kError;
  }
  assert(termination_status == TerminationStatus::OPTIMAL ||
         termination_status == TerminationStatus::TIMEOUT);
  highs_solution.col_value = x;
  highs_solution.col_value.resize(lp.num_col_);
  highs_solution.row_dual = y;
  lp.a_matrix_.product(highs_solution.row_value, highs_solution.col_value);
  lp.a_matrix_.productTranspose(highs_solution.col_dual,
                                highs_solution.row_dual);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    highs_solution.col_dual[iCol] =
        lp.col_cost_[iCol] - highs_solution.col_dual[iCol];
  highs_solution.value_valid = true;
  highs_solution.dual_valid = true;
  return HighsStatus::kOk;
}
