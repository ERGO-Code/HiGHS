/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available Hias open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_MIP_WORKER_H_
#define HIGHS_MIP_WORKER_H_

// #include "mip/HighsCliqueTable.h"
#include "mip/HighsConflictPool.h"
#include "mip/HighsCutPool.h"

// #include "mip/HighsDomain.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"

// #include "mip/HighsNodeQueue.h"
#include "mip/HighsPseudocost.h"
// #include "mip/HighsSeparation.h"
// #include "presolve/HighsSymmetry.h"
// #include "util/HighsHash.h"

class HighsSearch;

class HighsMipWorker {
 public: 

  const HighsMipSolver& mipsolver_;
  const HighsMipSolverData& mipdata_;

  // HighsCliqueTable cliquetable_;

  HighsPseudocost pseudocost_;

  std::unique_ptr<HighsSearch> search_ptr_;

  // HighsMipWorker(const HighsMipSolver& mipsolver__);
  HighsMipWorker(const HighsMipSolver& mipsolver__, const HighsLpRelaxation& lprelax_);
 
  // HighsMipWorker(const HighsMipWorker& mipworker);

  const HighsMipSolver& getMipSolver();

  HighsLpRelaxation lprelaxation_;

  HighsCutPool cutpool_;
  HighsConflictPool conflictpool_;

  // members for worker threads.
  HighsPseudocost pscost_;
  // HighsCliqueTable clqtableinit_;
  /// HighsImplications implicinit_;

  struct Solution {
    double row_violation_;
    double bound_violation_;
    double integrality_violation_;
    std::vector<double> solution_;
    double solution_objective_;
  };

  const bool checkLimits(int64_t nodeOffset = 0) const;

  // ... implement necessary methods for HighsSearch
 
  ~HighsMipWorker() {
      // search_ptr_.release();
      search_ptr_.reset();
  }

  HighsPrimalHeuristics::Statistics heur_stats;

  // todo: 
  // timer_
  // sync too
  // or name times differently for workers in the same timer instance in mipsolver.

};

#endif
