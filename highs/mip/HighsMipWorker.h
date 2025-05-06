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
#include "mip/HighsNodeQueue.h"
#include "mip/HighsPrimalHeuristics.h"
#include "mip/HighsPseudocost.h"
// #include "mip/HighsSeparation.h"

class HighsSearch;

class HighsMipWorker {
 public:
  const HighsMipSolver& mipsolver_;
  const HighsMipSolverData& mipdata_;

  HighsPseudocost pseudocost_;

  std::unique_ptr<HighsSearch> search_ptr_;

  const HighsMipSolver& getMipSolver();

  HighsLpRelaxation& lprelaxation_;

  HighsCutPool cutpool_;
  HighsConflictPool conflictpool_;

  HighsNodeQueue nodequeue_;

  struct Solution {
    double row_violation_;
    double bound_violation_;
    double integrality_violation_;
    std::vector<double> solution_;
    double solution_objective_;
  };

  double upper_bound;

  Solution solution_;

  HighsPrimalHeuristics::Statistics heur_stats;

  HighsRandom randgen;

  // HighsMipWorker(const HighsMipSolver& mipsolver__);
  HighsMipWorker(const HighsMipSolver& mipsolver__,
                 HighsLpRelaxation& lprelax_);

  ~HighsMipWorker() {
    // search_ptr_.release();
    search_ptr_.reset();
  }

  const bool checkLimits(int64_t nodeOffset = 0) const;

  bool addIncumbent(const std::vector<double>& sol, double solobj,
                    const int solution_source,
                    const bool print_display_line = true);

  double transformNewIntegerFeasibleSolution(
      const std::vector<double>& sol,
      const bool possibly_store_as_new_incumbent = true);

  bool syncSolution(double& row_violation,
    double& bound_violation,
    double& integrality_violation,
    std::vector<double>& solution,
    double& solution_objective) {
      row_violation = solution_.row_violation_;
      bound_violation = solution_.bound_violation_;
      integrality_violation = solution_.integrality_violation_;
      solution = solution_.solution_;
      solution_objective = solution_.solution_objective_;

      return true;
    }

  bool syncNodeQueue(HighsNodeQueue& nodequeue);

  // todo:
  // timer_
  // sync too
  // or name times differently for workers in the same timer instance in
  // mipsolver.
};

#endif
