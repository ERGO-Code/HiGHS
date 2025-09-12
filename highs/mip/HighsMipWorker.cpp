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
                               HighsLpRelaxation& lprelax_)
    : mipsolver_(mipsolver__),
      mipdata_(*mipsolver_.mipdata_.get()),
      lprelaxation_(lprelax_),
      pseudocost_(mipsolver__),
      cutpool_(mipsolver_.numCol(), mipsolver_.options_mip_->mip_pool_age_limit,
               mipsolver_.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver_.options_mip_->mip_pool_age_limit,
                    mipsolver_.options_mip_->mip_pool_soft_limit),
      upper_bound(kHighsInf) {
  // std::cout << mipdata_.domain.changedcolsflags_.size() << std::endl;
  search_ptr_ =
      std::unique_ptr<HighsSearch>(new HighsSearch(*this, pseudocost_));

  // Register cutpool and conflict pool in local search domain.
  // Add global cutpool.
  search_ptr_->getLocalDomain().addCutpool(mipsolver_.mipdata_->cutpool);
  search_ptr_->getLocalDomain().addConflictPool(
      mipsolver_.mipdata_->conflictPool);

  // cutpool_.matrix_.AheadNeg_.assign(mipsolver__.numCol(), -1);
  // cutpool_.matrix_.AheadPos_.assign(mipsolver__.numCol(), -1);

  // std::vector<HighsInt> AheadPos_;
  // std::vector<HighsInt> AheadNeg_;

  // add local cutpool
  search_ptr_->getLocalDomain().addCutpool(cutpool_);
  search_ptr_->getLocalDomain().addConflictPool(conflictpool_);

  // printf(
  //     "lprelax_ parameter address in constructor of mipworker %p, %d columns,
  //     " "and "
  //     "%d rows\n",
  //     (void*)&lprelax_, int(lprelax_.getLpSolver().getNumCol()),
  //     int(lprelax_.getLpSolver().getNumRow()));

  // printf(
  //     "lprelaxation_ address in constructor of mipworker %p, %d columns, and
  //     "
  //     "%d rows\n",
  //     (void*)&lprelaxation_, int(lprelaxation_.getLpSolver().getNumCol()),
  //     int(lprelaxation_.getLpSolver().getNumRow()));

  // HighsSearch has its own relaxation initialized no nullptr.

  search_ptr_->setLpRelaxation(&lprelaxation_);

  // printf(
  //     "Search has lp member in constructor of mipworker with address %p, %d "
  //     "columns, and %d rows\n",
  //     (void*)&search_ptr_->lp,
  //     int(search_ptr_->lp->getLpSolver().getNumCol()),
  //     int(search_ptr_->lp->getLpSolver().getNumRow()));
}

const HighsMipSolver& HighsMipWorker::getMipSolver() { return mipsolver_; }

void HighsMipWorker::resetSearchDomain() {
  search_ptr_.reset();
  search_ptr_ =
      std::unique_ptr<HighsSearch>(new HighsSearch(*this, pseudocost_));
  // search_ptr_->getLocalDomain().addCutpool(mipsolver_.mipdata_->cutpool);
  // search_ptr_->getLocalDomain().addConflictPool(
  //     mipsolver_.mipdata_->conflictPool);
  cutpool_ = HighsCutPool(mipsolver_.numCol(),
                          mipsolver_.options_mip_->mip_pool_age_limit,
                          mipsolver_.options_mip_->mip_pool_soft_limit);
  conflictpool_ =
      HighsConflictPool(5 * mipsolver_.options_mip_->mip_pool_age_limit,
                        mipsolver_.options_mip_->mip_pool_soft_limit);
  search_ptr_->getLocalDomain().addCutpool(cutpool_);
  search_ptr_->getLocalDomain().addConflictPool(conflictpool_);
  search_ptr_->setLpRelaxation(&lprelaxation_);
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