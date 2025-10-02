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

HighsMipWorker::HighsMipWorker(const HighsMipSolver& mipsolver__,
                               HighsLpRelaxation* lprelax_, HighsDomain* domain,
                               HighsCutPool* cutpool,
                               HighsConflictPool* conflictpool)
    : mipsolver_(mipsolver__),
      mipdata_(*mipsolver_.mipdata_.get()),
      lprelaxation_(lprelax_),
      pseudocost_(mipsolver__),
      globaldom_(domain),
      cutpool_(cutpool),
      conflictpool_(conflictpool) {
  upper_bound = mipdata_.upper_bound;
  search_ptr_ =
      std::unique_ptr<HighsSearch>(new HighsSearch(*this, pseudocost_));
  sepa_ptr_ = std::unique_ptr<HighsSeparation>(new HighsSeparation(*this));

  // add local cutpool
  search_ptr_->getLocalDomain().addCutpool(*cutpool_);
  search_ptr_->getLocalDomain().addConflictPool(*conflictpool_);

  search_ptr_->setLpRelaxation(lprelaxation_);
  sepa_ptr_->setLpRelaxation(lprelaxation_);
}

const HighsMipSolver& HighsMipWorker::getMipSolver() { return mipsolver_; }

void HighsMipWorker::resetSearch() {
  search_ptr_.reset();
  search_ptr_ =
      std::unique_ptr<HighsSearch>(new HighsSearch(*this, pseudocost_));
  search_ptr_->setLpRelaxation(lprelaxation_);
}

void HighsMipWorker::resetSepa() {
  sepa_ptr_.reset();
  sepa_ptr_ = std::unique_ptr<HighsSeparation>(new HighsSeparation(*this));
  sepa_ptr_->setLpRelaxation(lprelaxation_);
}

bool HighsMipWorker::addIncumbent(const std::vector<double>& sol, double solobj,
                                  int solution_source) {
  if (solobj < upper_bound) {
    // Get the transformed objective and solution if required
    const std::pair<bool, double> transformed_solobj =
        transformNewIntegerFeasibleSolution(sol);
    if (transformed_solobj.first && transformed_solobj.second < upper_bound) {
      upper_bound = transformed_solobj.second;
    }
    // Can't repair solutions locally, so also buffer infeasible ones
    solutions_.emplace_back(sol, solobj, solution_source);
  }
  return true;
}

std::pair<bool, double> HighsMipWorker::transformNewIntegerFeasibleSolution(
    const std::vector<double>& sol) {
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