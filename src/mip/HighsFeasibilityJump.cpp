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
#include "util/HighsSparseMatrix.h"

void HighsMipSolverData::feasibilityJump() {
  // This is the (presolved) model being solved
  const HighsLp* model = this->mipsolver.model_;
  const HighsLogOptions& log_options = mipsolver.options_mip_->log_options;
  highsLogUser(
      log_options, HighsLogType::kInfo,
      "HighsMipSolverData::feasibilityJump called with primal bound of %g\n",
      lower_bound);
#ifdef HIGHSINT64
  // TODO(BenChampion): make FJ work with 64-bit HighsInt
  highsLogUser(log_options, HighsLogType::kInfo,
               "Feasibility Jump code uses 'int' so isn't currently compatible "
               "with a 64-bit HighsInt. Skipping Feasibility Jump.\n");
  return;
#else
  // TODO(BenChampion): pick more thoughtful values for these
  const HighsInt kMaxTotalEffort = 1e6;
  const HighsInt kMaxEffortSinceLastImprovement = 1e3;

  bool found_integer_feasible_solution = false;
  std::vector<double> col_value(model->num_col_, 0.0);
  double objective_function_value;

  // Configure Feasibility Jump and pass it the problem
  external_feasibilityjump::equalityTolerance = epsilon;
  external_feasibilityjump::violationTolerance = feastol;

  auto solver = external_feasibilityjump::FeasibilityJumpSolver(0, 1);

  for (HighsInt col = 0; col < model->num_col_; ++col) {
    assert(model->integrality_[col] == HighsVarType::kContinuous ||
           model->integrality_[col] == HighsVarType::kInteger ||
           model->integrality_[col] == HighsVarType::kImplicitInteger);
    external_feasibilityjump::VarType fjVarType;
    if (model->integrality_[col] == HighsVarType::kContinuous) {
      fjVarType = external_feasibilityjump::VarType::Continuous;
    } else {
      fjVarType = external_feasibilityjump::VarType::Integer;
    }
    // TODO(BenChampion): do we handle sense of objective correctly?
    solver.addVar(fjVarType, model->col_lower_[col], model->col_upper_[col],
                  model->col_cost_[col]);
    // TODO(BenChampion): what about other infeasibilities/unboundedness?
    if (model->col_lower_[col] > model->col_upper_[col]) {
      highsLogUser(
          log_options, HighsLogType::kInfo,
          "Detected infeasible column bounds. Skipping Feasibility Jump");
      return;
    }
    // TODO(BenChampion): any other cases where infinite bounds are problematic?
    double initial_assignment = 0;
    if (std::isfinite(model->col_lower_[col])) {
      initial_assignment = model->col_lower_[col];
    } else if (std::isfinite(model->col_upper_[col])) {
      initial_assignment = model->col_upper_[col];
    }
    col_value[col] = initial_assignment;
  }

  HighsSparseMatrix a_matrix;
  a_matrix.createRowwise(model->a_matrix_);

  for (HighsInt row = 0; row < model->num_row_; ++row) {
    bool hasFiniteLower = std::isfinite(model->row_lower_[row]);
    bool hasFiniteUpper = std::isfinite(model->row_upper_[row]);
    if (hasFiniteLower || hasFiniteUpper) {
      HighsInt row_num_nz = a_matrix.start_[row + 1] - a_matrix.start_[row];
      auto row_index = a_matrix.index_.data() + a_matrix.start_[row];
      auto row_value = a_matrix.value_.data() + a_matrix.start_[row];
      if (hasFiniteLower) {
        solver.addConstraint(external_feasibilityjump::RowType::Gte,
                             model->row_lower_[row], row_num_nz, row_index,
                             row_value, /* relax_continuous = */ 0);
      }
      if (hasFiniteUpper) {
        solver.addConstraint(external_feasibilityjump::RowType::Lte,
                             model->row_upper_[row], row_num_nz, row_index,
                             row_value, /* relax_continuous = */ 0);
      }
    }
  }

  auto fjControlCallback =
      [=, &col_value, &found_integer_feasible_solution,
       &objective_function_value](external_feasibilityjump::FJStatus status)
      -> external_feasibilityjump::CallbackControlFlow {
    // TODO(BenChampion): these are really debug, remove/reduce
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
    // TODO(BenChampion): which solution to pick? First? Last? Best?
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

  solver.solve(col_value.data(), fjControlCallback);
  // TODO(BenChampion): remove all DEBUG
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
    if (!trySolution(col_value, kSolutionSourceFeasibilityJump)) {
      printf("DEBUG: Feasibility Jump solution was not integer feasible\n");
    } else {
      printf("DEBUG: Feasibility Jump solution is integer feasible.\n");
    }

  } else {
    //    highsLogUser(log_options, HighsLogType::kInfo,
    printf(
        "DEBUG: Feasibility Jump has not found an integer feasible "
        "solution\n");
  }
#endif
}
