/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipWorker.h"

#include "mip/HighsMipSolverData.h"

HighsMipWorker::HighsMipWorker(const HighsMipSolver& mipsolver__,
                               const HighsLpRelaxation& lprelax_)
    : mipsolver_(mipsolver__),
      mipdata_(*mipsolver_.mipdata_.get()),
      lprelaxation_(lprelax_),
      cutpool_(mipsolver_.numCol(),
               mipsolver_.options_mip_->mip_pool_age_limit,
               mipsolver_.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver_.options_mip_->mip_pool_age_limit,
                    mipsolver_.options_mip_->mip_pool_soft_limit),
      pseudocost_(mipsolver__)
      {

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
  search_ptr_->setLpRelaxation(&lprelaxation_);

  printf(
      "lprelax_ parameter address in constructor of mipworker %p, %d columns, and "
      "%d rows\n",
      (void*)&lprelax_, int(lprelax_.getLpSolver().getNumCol()),
      int(lprelax_.getLpSolver().getNumRow()));

 printf(
      "lprelaxation_ address in constructor of mipworker %p, %d columns, and "
      "%d rows\n",
      (void*)&lprelaxation_, int(lprelaxation_.getLpSolver().getNumCol()),
      int(lprelaxation_.getLpSolver().getNumRow()));

  // HighsSearch has its own relaxation initialized no nullptr.
  search_ptr_->setLpRelaxation(&lprelaxation_);

  printf(
      "Search has lp member in constructor of mipworker with address %p, %d "
      "columns, and %d rows\n",
      (void*)&search_ptr_->lp, int(search_ptr_->lp->getLpSolver().getNumCol()),
      int(search_ptr_->lp->getLpSolver().getNumRow()));
}

const HighsMipSolver& HighsMipWorker::getMipSolver() { return mipsolver_; }