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

HighsModelStatus HighsMipSolverData::feasibilityJump() {
  // This is the (presolved) model being solved
  const HighsLp* model = this->mipsolver.model_;
  const HighsLogOptions& log_options = mipsolver.options_mip_->log_options;
  double sense_multiplier = static_cast<double>(model->sense_);

#ifdef HIGHSINT64
  // TODO(BenChampion,9999-12-31): make FJ work with 64-bit HighsInt
  highsLogUser(log_options, HighsLogType::kInfo,
               "Feasibility Jump code isn't currently compatible "
               "with a 64-bit HighsInt: skipping Feasibility Jump\n");
  return HighsModelStatus::kNotset;
#else
  const size_t kMaxTotalEffort = std::pow(2, 30);  // Originally 1e6
  const size_t kMaxEffortSinceLastImprovement =
      std::pow(2, 25);  // Originally 1e3

  bool found_integer_feasible_solution = false;
  std::vector<double> col_value(model->num_col_, 0.0);
  double objective_function_value;

  // Configure Feasibility Jump and pass it the problem
  external_feasibilityjump::equalityTolerance = epsilon;
  external_feasibilityjump::violationTolerance = feastol;

  int verbosity = mipsolver.submip ? 0 : mipsolver.options_mip_->log_dev_level;
  auto solver = external_feasibilityjump::FeasibilityJumpSolver(
      /* seed = */ 0, /* verbosity = */ verbosity);

  for (HighsInt col = 0; col < model->num_col_; ++col) {
    double lower = model->col_lower_[col];
    double upper = model->col_upper_[col];

    assert(model->integrality_[col] == HighsVarType::kContinuous ||
           model->integrality_[col] == HighsVarType::kInteger ||
           model->integrality_[col] == HighsVarType::kImplicitInteger);
    external_feasibilityjump::VarType fjVarType;
    if (model->integrality_[col] == HighsVarType::kContinuous) {
      fjVarType = external_feasibilityjump::VarType::Continuous;
    } else {
      fjVarType = external_feasibilityjump::VarType::Integer;
      lower = std::ceil(lower);
      upper = std::floor(upper);
    }

    const bool legal_bounds = lower <= upper && lower < kHighsInf &&
                              upper > -kHighsInf && !std::isnan(lower) &&
                              !std::isnan(upper);
    if (!legal_bounds) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "HighsMipSolverData::feasibilityJump() has detected "
                   "infeasible/illegal bounds [%g, %g] for "
                   "column %d: MIP is infeasible\n",
                   lower, upper, int(col));
      assert(legal_bounds);
      return HighsModelStatus::kInfeasible;
    }
    solver.addVar(fjVarType, lower, upper,
                  sense_multiplier * model->col_cost_[col]);

    double initial_assignment = 0.0;
    if (std::isfinite(lower)) {
      initial_assignment = lower;
    } else if (std::isfinite(upper)) {
      initial_assignment = upper;
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
    if (status.solution != nullptr) {
      found_integer_feasible_solution = true;
      col_value = std::vector<double>(status.solution,
                                      status.solution + status.numVars);
      objective_function_value =
          model->offset_ + sense_multiplier * status.solutionObjectiveValue;
      if (verbosity > 0) {
        printf("Feasibility Jump has found a solution");
        if (model->num_col_ < 10) {
          printf(" [");
          for (HighsInt col = 0; col < std::min(10, model->num_col_); ++col)
            printf(" %g", col_value[col]);
          printf("]");
        }
        printf(" with objective %g\n", objective_function_value);
      }
    }
    if (status.effortSinceLastImprovement > kMaxEffortSinceLastImprovement ||
        status.totalEffort > kMaxTotalEffort) {
      return external_feasibilityjump::CallbackControlFlow::Terminate;
    } else {
      return external_feasibilityjump::CallbackControlFlow::Continue;
    }
  };

  if (verbosity > 0)
    printf(
        "Feasibility Jump: kMaxTotalEffort = %zd; "
        "kMaxEffortSinceLastImprovement = %zd\n",
        kMaxTotalEffort, kMaxEffortSinceLastImprovement);
  solver.solve(col_value.data(), fjControlCallback);

  if (found_integer_feasible_solution) {
    // Initial assignments that violate integrality or column bounds can lead to
    // infeasible results. Even if those initial assignments should not occur,
    // use trySolution rather than addIncumbent for an explicit check.
    bool is_really_feasible =
        trySolution(col_value, kSolutionSourceFeasibilityJump);
    if (!is_really_feasible) {
      highsLogUser(log_options, HighsLogType::kInfo,
                   "Discarding infeasible result from Feasibility Jump\n");
    }
    assert(is_really_feasible);
  }
  return HighsModelStatus::kNotset;
#endif
}
