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

#include "mip/HighsMipSolverData.h"
#include "mip/feasibilityjump.hh"

void HighsMipSolverData::feasibilityJump() {
  // This is the (presolved) model being solved
  const HighsLp* model = this->mipsolver.model_;
  const HighsLogOptions& log_options = mipsolver.options_mip_->log_options;
  highsLogUser(
      log_options, HighsLogType::kInfo,
      "HighsMipSolverData::feasibilityJump called with primal bound of %g\n",
      lower_bound);
#ifdef HIGHSINT64
  highsLogUser(log_options, HighsLogType::kInfo,
               "Feasibility Jump code uses 'int' so isn't currently compatible "
               "with a 64-bit HighsInt. Skipping Feasibility Jump\n",
               lower_bound);
  return;
#else
  const HighsInt kMaxTotalEffort = 1e6;
  const HighsInt kMaxEffortSinceLastImprovement = 1e3;

  bool found_integer_feasible_solution = false;
  std::vector<double> col_value(model->num_col_, 0.0);
  double objective_function_value;

  // Configure Feasibility Jump and pass it the problem
  external_feasibilityjump::equalityTolerance = epsilon;
  external_feasibilityjump::violationTolerance = feastol;

  auto solver = external_feasibilityjump::FeasibilityJumpSolver(0, 1);

  for (int i = 0; i < model->num_col_; ++i) {
    assert(model->integrality_[i] == HighsVarType::kContinuous ||
           model->integrality_[i] == HighsVarType::kInteger ||
           model->integrality_[i] == HighsVarType::kImplicitInteger);
    external_feasibilityjump::VarType fjVarType;
    if (model->integrality_[i] == HighsVarType::kContinuous) {
      fjVarType = external_feasibilityjump::VarType::Continuous;
    } else {
      fjVarType = external_feasibilityjump::VarType::Integer;
    }
    solver.addVar(fjVarType, model->col_lower_[i], model->col_upper_[i],
                  model->col_cost_[i]);
  }

  HighsInt row_num_nz;
  HighsInt* row_index_buffer = new HighsInt[model->num_col_];
  double* row_value_buffer = new double[model->num_col_];

  for (int i = 0; i < model->num_row_; ++i) {
    bool hasFiniteLower = std::isfinite(model->row_lower_[i]);
    bool hasFiniteUpper = std::isfinite(model->row_upper_[i]);
    if (hasFiniteLower || hasFiniteUpper) {
      model->a_matrix_.getRow(i, row_num_nz, row_index_buffer,
                              row_value_buffer);
      // TODO: check if relax_continuous == 1 really is the right thing to
      // pass
      if (hasFiniteLower) {
        solver.addConstraint(external_feasibilityjump::RowType::Gte,
                             model->row_lower_[i], row_num_nz, row_index_buffer,
                             row_value_buffer, 1);
      }
      if (hasFiniteUpper) {
        solver.addConstraint(external_feasibilityjump::RowType::Lte,
                             model->row_upper_[i], row_num_nz, row_index_buffer,
                             row_value_buffer, 1);
      }
    }
  }

  delete[] row_value_buffer;
  delete[] row_index_buffer;

  auto fjControlCallback =
      [=, &col_value, &found_integer_feasible_solution,
       &objective_function_value](external_feasibilityjump::FJStatus status)
      -> external_feasibilityjump::CallbackControlFlow {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "From Feasibility Jump callback\n");
    highsLogUser(log_options, HighsLogType::kInfo, "Total effort: %d\n",
                 status.totalEffort);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Effort since last improvement: %d\n",
                 status.effortSinceLastImprovement);
    highsLogUser(log_options, HighsLogType::kInfo, "Number of variables: %d\n",
                 status.numVars);
    highsLogUser(log_options, HighsLogType::kInfo,
                 "Solution objective value: %0.2f\n",
                 status.solutionObjectiveValue);
    highsLogUser(log_options, HighsLogType::kInfo, "Solution found: %d\n",
                 status.solution != nullptr);
    if (status.solution != nullptr) {
      found_integer_feasible_solution = true;
      col_value = std::vector<double>(status.solution,
                                      status.solution + status.numVars);
      objective_function_value = status.solutionObjectiveValue;
    }
    if (status.effortSinceLastImprovement > kMaxEffortSinceLastImprovement ||
        status.totalEffort > kMaxTotalEffort) {
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
    for (HighsInt iCol = 0; iCol < std::min(10, int(col_value.size())); iCol++)
      printf(" %g", col_value[iCol]);
    printf("]\n");
    if (!trySolution(col_value, kSolutionSourceFeasibilityJump))
      printf("DEBUG: Feasibility Jump solution was not integer feasible\n");
  } else {
    //    highsLogUser(log_options, HighsLogType::kInfo,
    printf(
        "DEBUG: Feasibility Jump has not found an integer feasible "
        "solution\n");
  }
#endif
}
