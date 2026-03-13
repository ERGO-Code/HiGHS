/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipWorker.h"

#include "mip/HighsMipSolverData.h"
#include "mip/MipTimer.h"

HighsMipWorker::HighsMipWorker(const HighsMipSolver& mipsolver,
                               HighsLpRelaxation* lp, HighsDomain* domain,
                               HighsCutPool* cutpool,
                               HighsConflictPool* conflictpool,
                               HighsPseudocost* pseudocost)
    : mipsolver_(mipsolver),
      mipdata_(*mipsolver_.mipdata_),
      lp_(lp),
      globaldom_(domain),
      cutpool_(cutpool),
      conflictpool_(conflictpool),
      pseudocost_(pseudocost),
      randgen(mipsolver.options_mip_->random_seed) {
  upper_bound = mipdata_.upper_bound;
  upper_limit = mipdata_.upper_limit;
  optimality_limit = mipdata_.optimality_limit;
  search_ptr_ =
      std::unique_ptr<HighsSearch>(new HighsSearch(*this, getPseudocost()));
  sepa_ptr_ = std::unique_ptr<HighsSeparation>(new HighsSeparation(*this));
  search_ptr_->setLpRelaxation(lp_);
  sepa_ptr_->setLpRelaxation(lp_);
}

const HighsMipSolver& HighsMipWorker::getMipSolver() const {
  return mipsolver_;
}

void HighsMipWorker::resetSearch() {
  search_ptr_.reset();
  search_ptr_ =
      std::unique_ptr<HighsSearch>(new HighsSearch(*this, getPseudocost()));
  search_ptr_->setLpRelaxation(lp_);
}

void HighsMipWorker::resetSepa() {
  sepa_ptr_.reset();
  sepa_ptr_ = std::unique_ptr<HighsSeparation>(new HighsSeparation(*this));
  sepa_ptr_->setLpRelaxation(lp_);
}

bool HighsMipWorker::addIncumbent(const std::vector<double>& sol, double solobj,
                                  int solution_source) {
  if (solobj < upper_bound) {
    // Get the transformed objective and solution if required
    const std::pair<bool, double> transformed_solobj =
        transformNewIntegerFeasibleSolution(sol);
    if (transformed_solobj.first && transformed_solobj.second < upper_bound) {
      upper_bound = transformed_solobj.second;
      double new_upper_limit =
          mipdata_.computeNewUpperLimit(upper_bound, 0.0, 0.0);
      if (new_upper_limit < upper_limit) {
        upper_limit = new_upper_limit;
        optimality_limit = mipdata_.computeNewUpperLimit(
            upper_bound, mipsolver_.options_mip_->mip_abs_gap,
            mipsolver_.options_mip_->mip_rel_gap);
      }
    }
    // Can't repair solutions locally, so also buffer infeasible ones
    solutions_.emplace_back(sol, solobj, solution_source);
  }
  return true;
}

std::pair<bool, double> HighsMipWorker::transformNewIntegerFeasibleSolution(
    const std::vector<double>& sol) const {
  HighsSolution solution;
  solution.col_value = sol;
  solution.value_valid = true;

  // Perform primal postsolve to get the original column values
  mipsolver_.mipdata_->postSolveStack.undoPrimal(*mipsolver_.options_mip_,
                                                 solution, -1, true);

  // Determine the row values, as they aren't computed in primal
  // postsolve
  HighsStatus return_status =
      calculateRowValuesQuad(*mipsolver_.orig_model_, solution);
  if (kAllowDeveloperAssert) assert(return_status == HighsStatus::kOk);

  // compute the objective value in the original space
  double bound_violation_ = 0;
  double row_violation_ = 0;
  double integrality_violation_ = 0;

  HighsCDouble mipsolver_quad_objective_value = 0;

  bool feasible = mipsolver_.solutionFeasible(
      mipsolver_.orig_model_, solution.col_value, &solution.row_value,
      bound_violation_, row_violation_, integrality_violation_,
      mipsolver_quad_objective_value);

  const double transformed_solobj = static_cast<double>(
      static_cast<HighsInt>(mipsolver_.orig_model_->sense_) *
          mipsolver_quad_objective_value -
      mipsolver_.model_->offset_);

  return std::make_pair(feasible, transformed_solobj);
}

bool HighsMipWorker::trySolution(const std::vector<double>& solution,
                                 const int solution_source) {
  if (static_cast<int>(solution.size()) != mipsolver_.model_->num_col_)
    return false;

  HighsCDouble obj = 0;

  for (HighsInt i = 0; i != mipsolver_.model_->num_col_; ++i) {
    if (solution[i] < mipsolver_.model_->col_lower_[i] - mipdata_.feastol)
      return false;
    if (solution[i] > mipsolver_.model_->col_upper_[i] + mipdata_.feastol)
      return false;
    if (mipsolver_.variableType(i) == HighsVarType::kInteger &&
        fractionality(solution[i]) > mipdata_.feastol)
      return false;

    obj += mipsolver_.colCost(i) * solution[i];
  }

  for (HighsInt i = 0; i != mipsolver_.model_->num_row_; ++i) {
    double rowactivity = 0.0;

    HighsInt start = mipdata_.ARstart_[i];
    HighsInt end = mipdata_.ARstart_[i + 1];

    for (HighsInt j = start; j != end; ++j)
      rowactivity += solution[mipdata_.ARindex_[j]] * mipdata_.ARvalue_[j];

    if (rowactivity > mipsolver_.rowUpper(i) + mipdata_.feastol) return false;
    if (rowactivity < mipsolver_.rowLower(i) - mipdata_.feastol) return false;
  }

  return addIncumbent(solution, static_cast<double>(obj), solution_source);
}