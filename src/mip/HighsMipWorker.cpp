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
      cutpool_(mipsolver__.numCol(), mipsolver__.options_mip_->mip_pool_age_limit,
               mipsolver__.options_mip_->mip_pool_soft_limit),
      conflictpool_(5 * mipsolver__.options_mip_->mip_pool_age_limit,
                    mipsolver__.options_mip_->mip_pool_soft_limit),
      cliquetable_(mipsolver__.numCol()),
      mipsolver(mipsolver__),
      pseudocost(mipsolver),
      search_(mipsolver, pseudocost) {}


HighsMipSolver& HighsMipWorker::getMipSolver() {
  return mipsolver;
}