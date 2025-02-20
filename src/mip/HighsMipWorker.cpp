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
      cutpool_(mipsolver__.numCol(),
               mipsolver__.options_mip_->mip_pool_age_limit,
               mipsolver__.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver__.options_mip_->mip_pool_age_limit,
                    mipsolver__.options_mip_->mip_pool_soft_limit),
    //   cliquetable_(mipsolver__.numCol()),
      pseudocost_(mipsolver__),
      pscostinit_(pseudocost_, 1),
    //   clqtableinit_(mipsolver_.numCol()),
      implicinit_(mipsolver_),
      pscostinit(pscostinit_),
      implicinit(implicinit_)
    //   clqtableinit(clqtableinit_) 
      {

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
      "lprelaxation_ address in constructor of mipworker %p, %d columns, and "
      "%d rows\n",
      (void*)&lprelaxation_, int(lprelaxation_.getLpSolver().getNumCol()),
      int(lprelaxation_.getLpSolver().getNumRow()));

  printf(
      "Search has lp member in constructor of mipworker with address %p, %d "
      "columns, and %d rows\n",
      (void*)&search_ptr_->lp, int(search_ptr_->lp->getLpSolver().getNumCol()),
      int(search_ptr_->lp->getLpSolver().getNumRow()));
}

//  HighsMipWorker::~HighsMipWorker() {
//     delete search_ptr;
//   };

const HighsMipSolver& HighsMipWorker::getMipSolver() { return mipsolver_; }

// HighsSearch& HighsMipWorker::getSearch() { return (*search_ptr); }