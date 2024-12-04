/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <climits>

#include "lp_data/HighsModelUtils.h"  // for typeToString Remove 1423
#include "mip/HighsMipSolverData.h"
#include "mip/feasibilityjump.hh"

void HighsMipSolverData::feasibilityJump() {
  // This is the (presolved) model being solved
  const HighsLp* model = this->mipsolver.model_;
  const HighsLogOptions& log_options = mipsolver.options_mip_->log_options;
  printf("HighsMipSolverData::feasibilityJump called with primal bound of %g\n",
         lower_bound);

  const HighsInt MAX_TOTAL_EFFORT = 1e6;
  const HighsInt MAX_EFFORT_SINCE_LAST_IMPROVEMENT = 1e3;

  bool found_integer_feasible_solution = false;
  std::vector<double> col_value(model->num_col_, 0.0);
  double objective_function_value;

  printf("DEBUG: configuring feasibility jump\n");

  // Configure Feasibility Jump and pass it the problem
  external_feasibilityjump::equalityTolerance = epsilon;
  external_feasibilityjump::violationTolerance = feastol;

  auto solver = external_feasibilityjump::FeasibilityJumpSolver(0, 1);

  printf("DEBUG: adding cols (vars)\n");

  for (int i = 0; i < model->num_col_; ++i) {
    external_feasibilityjump::VarType fjVarType;
    if (model->integrality_[i] == HighsVarType::kContinuous) {
      fjVarType = external_feasibilityjump::VarType::Continuous;
    } else if (model->integrality_[i] == HighsVarType::kInteger) {
      fjVarType = external_feasibilityjump::VarType::Integer;
    } else {
      printf(
          "Feasibility Jump only supports continuous and integer variables, "
          "but integrality_[%d] is %s. "
          "Skipping Feasibility Jump...\n",
          int(i), typeToString(model->integrality_[i]).c_str());
      return;
    }
    solver.addVar(fjVarType, model->col_lower_[i], model->col_upper_[i],
                  model->col_cost_[i]);
  }

  printf("DEBUG: adding rows (constraints)\n");

  HighsInt numNzThisRow;
  HighsInt* rowIndexBuffer = new HighsInt[model->num_col_];
  double* rowValueBuffer = new double[model->num_col_];

  for (int i = 0; i < model->num_row_; ++i) {
    bool hasFiniteLower = std::isfinite(model->row_lower_[i]);
    bool hasFiniteUpper = std::isfinite(model->row_upper_[i]);
    if (hasFiniteLower || hasFiniteUpper) {
      model->a_matrix_.getRow(i, numNzThisRow, rowIndexBuffer, rowValueBuffer);
      // TODO: check if relax_continuous == 1 really is the right thing to pass
      if (hasFiniteLower) {
        solver.addConstraint(external_feasibilityjump::RowType::Gte,
                             model->row_lower_[i], numNzThisRow, rowIndexBuffer,
                             rowValueBuffer, 1);
      }
      if (hasFiniteUpper) {
        solver.addConstraint(external_feasibilityjump::RowType::Lte,
                             model->row_upper_[i], numNzThisRow, rowIndexBuffer,
                             rowValueBuffer, 1);
      }
    }
  }

  delete[] rowValueBuffer;
  delete[] rowIndexBuffer;

  printf("DEBUG: calling Feasibility Jump\n");

  auto fjControlCallback =
      [=, &col_value, &found_integer_feasible_solution,
       &objective_function_value](external_feasibilityjump::FJStatus status)
      -> external_feasibilityjump::CallbackControlFlow {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "From Feasibility Jump callback\n");
    if (status.solution != nullptr) {
      found_integer_feasible_solution = true;
      col_value = std::vector<double>(status.solution,
                                      status.solution + status.numVars);
      objective_function_value = status.solutionObjectiveValue;
    }
    if (status.effortSinceLastImprovement > MAX_EFFORT_SINCE_LAST_IMPROVEMENT ||
        status.totalEffort > MAX_TOTAL_EFFORT) {
      return external_feasibilityjump::CallbackControlFlow::Terminate;
    } else {
      return external_feasibilityjump::CallbackControlFlow::Continue;
    }
  };

  solver.solve(nullptr, fjControlCallback);

  if (found_integer_feasible_solution) {
    // Feasibility jump has found a solution, so call addIncumbent to
    // (possibly) update the incumbent
    //    highsLogUser(log_options, HighsLogType::kInfo,
    printf(
        "DEBUG: Feasibility Jump has found an integer feasible solution with "
        "objective value %g\n",
        objective_function_value);
    printf("DEBUG: Solution: [");
    for(HighsInt iCol = 0; iCol < std::min(10, int(col_value.size())); iCol++)
      printf(" %g", col_value[iCol]);
    printf("]\n");
    if (!trySolution(col_value, kSolutionSourceFeasibilityJump))
      printf("DEBUG: Feasibility Jump solution was not integer feasible\n");
  } else {
    //    highsLogUser(log_options, HighsLogType::kInfo,
    printf(
        "DEBUG: Feasibility Jump has not found an integer feasible solution\n");
  }
}
