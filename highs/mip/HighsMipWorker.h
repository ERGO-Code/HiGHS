/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available Hias open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_MIP_WORKER_H_
#define HIGHS_MIP_WORKER_H_

#include "mip/HighsConflictPool.h"
#include "mip/HighsCutPool.h"
#include "mip/HighsImplications.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsPrimalHeuristics.h"
#include "mip/HighsPseudocost.h"
// #include "mip/HighsSeparation.h"

class HighsSearch;

class HighsMipWorker {
 public:
  const HighsMipSolver& mipsolver_;
  const HighsMipSolverData& mipdata_;

  HighsPseudocost pseudocost_;
  HighsLpRelaxation& lprelaxation_;
  HighsDomain& globaldom_;
  HighsCutPool* cutpool_;
  HighsConflictPool& conflictpool_;

  std::unique_ptr<HighsSearch> search_ptr_;

  const HighsMipSolver& getMipSolver();

  double upper_bound;

  std::vector<std::tuple<std::vector<double>, double, int>> solutions_;

  HighsPrimalHeuristics::Statistics heur_stats;

  HighsRandom randgen;

  // HighsMipWorker(const HighsMipSolver& mipsolver__);
  HighsMipWorker(const HighsMipSolver& mipsolver__,
                 HighsLpRelaxation& lprelax_,
                 HighsDomain& domain,
                 HighsCutPool* cutpool,
                 HighsConflictPool& conflictpool);

  ~HighsMipWorker() {
    // search_ptr_.release();
    search_ptr_.reset();
  }

  const bool checkLimits(int64_t nodeOffset = 0) const;

  void resetSearch();

  // bool addIncumbent(const std::vector<double>& sol, double solobj,
  //                   const int solution_source,
  //                   const bool print_display_line = true);

  bool addIncumbent(const std::vector<double>& sol, double solobj, int solution_source);

  std::pair<bool, double> transformNewIntegerFeasibleSolution(
      const std::vector<double>& sol);

  // todo:
  // timer_
  // sync too
  // or name times differently for workers in the same timer instance in
  // mipsolver.
};

#endif
