/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available Hias open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_MIP_WORKER_H_
#define HIGHS_MIP_WORKER_H_

#include "mip/HighsCliqueTable.h"
#include "mip/HighsConflictPool.h"
#include "mip/HighsCutPool.h"
//#include "mip/HighsDomain.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsSearch.h"
//#include "mip/HighsNodeQueue.h"
//#include "mip/HighsPseudocost.h"
//#include "mip/HighsSeparation.h"
//#include "presolve/HighsSymmetry.h"
//#include "util/HighsHash.h"

class HighsMipWorker {
  const HighsMipSolver& mipsolver_;
  HighsLpRelaxation lprelaxation_;
  HighsCutPool cutpool_;
  HighsConflictPool conflictpool_;
  HighsCliqueTable cliquetable_;

  // Not sure if this should be here or elsewhere.
  HighsMipSolver mipsolver;

  // Not sure if this should be here or in HighsSearch.
  HighsPseudocost pseudocost;
  
  HighsSearch search_;

public:
  HighsMipWorker(const HighsMipSolver& mipsolver);

  HighsMipSolver& getMipSolver();
};

#endif
