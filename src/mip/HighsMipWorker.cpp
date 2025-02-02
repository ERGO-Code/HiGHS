/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsMipWorker.h"

HighsMipWorker::HighsMipWorker(HighsMipSolver& mipsolver)
  : mipsolver_(mipsolver)//,
    //    cutpool_(mipsolver.numCol(), mipsolver.options_mip_->mip_pool_age_limit,
    //                mipsolver.options_mip_->mip_pool_soft_limit),
    //    conflictPool_(5 * mipsolver.options_mip_->mip_pool_age_limit,
    //		 mipsolver.options_mip_->mip_pool_soft_limit),
    //    cliquetable_(mipsolver.numCol()),
    //    lprelaxation_(mipsolver)
{}

