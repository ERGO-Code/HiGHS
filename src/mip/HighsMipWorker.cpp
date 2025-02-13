/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipWorker.h"

#include "mip/HighsSearch.h"
#include "mip/HighsMipSolverData.h"

HighsMipWorker::HighsMipWorker(const HighsMipSolver& mipsolver__, const HighsLpRelaxation& lprelax_)
    : mipsolver_(mipsolver__),
      // mipsolver_worker_(mipsolver__),
      // lprelaxation_(mipsolver__),
      // lprelaxation_(mipsolver__) required setLpRelaxation to be called after, but here
      // we use the local relaxation so we can initialize it in the constructor
      lprelaxation_(lprelax_),
      cutpool_(mipsolver__.numCol(),
               mipsolver__.options_mip_->mip_pool_age_limit,
               mipsolver__.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver__.options_mip_->mip_pool_age_limit,
                    mipsolver__.options_mip_->mip_pool_soft_limit),
      cliquetable_(mipsolver__.numCol()),
      // mipsolver(mipsolver__),
      pseudocost_(mipsolver__),
      // search_(mipsolver_, pseudocost_),

      pscostinit_(pseudocost_, 1),
      clqtableinit_(mipsolver_.numCol()),
      implicinit_(mipsolver_),

      pscostinit(pscostinit_),
      implicinit(implicinit_),
      clqtableinit(clqtableinit_) {
  // Register cutpool and conflict pool in local search domain.

  // search_ptr_= std::unique_ptr<HighsSearch>(new HighsSearch(mipsolver_, pseudocost_));

  search_ptr_= std::unique_ptr<HighsSearch>(new HighsSearch(*this, pseudocost_));
  // search_ptr_shared_ = std::shared_ptr<HighsSearch>(new HighsSearch(*this, pseudocost_));
  // search_ptr = new HighsSearch(*this, pseudocost_);

  // add global cutpool 
  // search_ptr_->getLocalDomain().addCutpool(mipsolver__.mipdata_->cutpool);
  // search_ptr_->getLocalDomain().addConflictPool(mipsolver_.mipdata_->conflictPool);

  // cutpool_.matrix_.AheadNeg_.assign(mipsolver__.numCol(), -1);
  // cutpool_.matrix_.AheadPos_.assign(mipsolver__.numCol(), -1);

  // std::vector<HighsInt> AheadPos_;
  // std::vector<HighsInt> AheadNeg_;

  // add local cutpool 
  search_ptr_->getLocalDomain().addCutpool(cutpool_);
  search_ptr_->getLocalDomain().addConflictPool(conflictpool_);
  search_ptr_->setLpRelaxation(&lprelaxation_);

  // search_ptr_shared_->getLocalDomain().addCutpool(cutpool_);
  // search_ptr_shared_->getLocalDomain().addConflictPool(conflictpool_);
  // search_ptr_shared_->setLpRelaxation(&lprelaxation_);

  // search_ptr->getLocalDomain().addCutpool(cutpool_);
  // search_ptr->getLocalDomain().addConflictPool(conflictpool_);
  // search_ptr->setLpRelaxation(&lprelaxation_);


      printf("lprelaxation_ address in constructor of mipworker %p, %d columns, and %d rows\n",
	     (void*)&lprelaxation_,
	     int(lprelaxation_.getLpSolver().getNumCol()),
	     int(lprelaxation_.getLpSolver().getNumRow()));

      printf("Search has lp member in constructor of mipworker with address %p, %d columns, and %d rows\n",
	     (void*)&search_ptr_->lp,
	     int(search_ptr_->lp->getLpSolver().getNumCol()),
	     int(search_ptr_->lp->getLpSolver().getNumRow()));

      // printf("Search has lp member in constructor of mipworker with address %p, %d columns, and %d rows\n",
	    //  (void*)&search_ptr_shared_->lp,
	    //  int(search_ptr_shared_->lp->getLpSolver().getNumCol()),
	    //  int(search_ptr_shared_->lp->getLpSolver().getNumRow()));

      // printf("Search has lp member in constructor of mipworker with address %p, %d columns, and %d rows\n",
	    //  (void*)search_ptr->lp,
	    //  int(search_ptr->lp->getLpSolver().getNumCol()),
	    //  int(search_ptr->lp->getLpSolver().getNumRow()));

  // Initialize mipdata_.
  // mipdata_ = decltype(mipdata_)(new HighsMipSolverData(mipsolver__));
  // mipdata_->init();
}

//  HighsMipWorker::~HighsMipWorker() {
//     delete search_ptr;
//   };

const HighsMipSolver& HighsMipWorker::getMipSolver() { return mipsolver_; }

HighsSearch& HighsMipWorker::getSearch() { return *search_ptr_; }
// HighsSearch& HighsMipWorker::getSearch() { return (*search_ptr); }
// HighsSearch& HighsMipWorker::getSearch() { return (*search_ptr_shared_); }