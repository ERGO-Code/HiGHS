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

// #include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsSearch.h"

// #include "mip/HighsNodeQueue.h"
// #include "mip/HighsPseudocost.h"
// #include "mip/HighsSeparation.h"
// #include "presolve/HighsSymmetry.h"
// #include "util/HighsHash.h"

class HighsMipWorker {
  const HighsMipSolver& mipsolver_;
  HighsLpRelaxation lprelaxation_;
  HighsCutPool cutpool_;
  HighsConflictPool conflictpool_;
  HighsCliqueTable cliquetable_;

  // Not sure if this should be here or elsewhere.
  // HighsMipSolver mipsolver;

  // Not sure if this should be here or in HighsSearch.
  HighsPseudocost pseudocost_;

  HighsSearch search_;

 public:
  HighsMipWorker(const HighsMipSolver& mipsolver__);

  const HighsMipSolver& getMipSolver();

  HighsSearch& getSearch();

  // members for worker threads.
  HighsPseudocostInitialization pscostinit_;
  HighsCliqueTable clqtableinit_;
  HighsImplications implicinit_;

  // References to members, initialized to local objects for worker threads,
  // modify to mip solver for main worker.
  HighsPseudocostInitialization& pscostinit;
  HighsCliqueTable& clqtableinit;
  HighsImplications& implicinit;

  // std::unique_ptr<HighsMipSolverData> mipdata_;
};

#endif
