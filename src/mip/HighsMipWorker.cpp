/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipWorker.h"

HighsMipWorker::HighsMipWorker(const HighsMipSolver& mipsolver__)
    : mipsolver_(mipsolver__),
      lprelaxation_(mipsolver__),
      cutpool_(mipsolver__.numCol(),
               mipsolver__.options_mip_->mip_pool_age_limit,
               mipsolver__.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver__.options_mip_->mip_pool_age_limit,
                    mipsolver__.options_mip_->mip_pool_soft_limit),
      cliquetable_(mipsolver__.numCol()),
      // mipsolver(mipsolver__),
      pseudocost_(mipsolver__),
      // search_(*this, pseudocost_),
      search_(mipsolver_, pseudocost_),

      pscostinit_(pseudocost_, 1),
      clqtableinit_(mipsolver_.numCol()),
      implicinit_(mipsolver_),

      pscostinit(pscostinit_),
      implicinit(implicinit_),
      clqtableinit(clqtableinit_) {
  // Register cutpool and conflict pool in local search domain.
  search_.getLocalDomain().addCutpool(cutpool_);
  search_.getLocalDomain().addConflictPool(conflictpool_);

  // Initialize mipdata_.
  // mipdata_ = decltype(mipdata_)(new HighsMipSolverData(mipsolver__));
  // mipdata_->init();
}

const HighsMipSolver& HighsMipWorker::getMipSolver() { return mipsolver_; }

HighsSearch& HighsMipWorker::getSearch() { return search_; }