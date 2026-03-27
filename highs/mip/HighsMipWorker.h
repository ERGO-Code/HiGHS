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
#include "mip/HighsNodeQueue.h"
#include "mip/HighsPrimalHeuristics.h"
#include "mip/HighsPseudocost.h"
#include "mip/HighsSeparation.h"

class HighsSearch;

class HighsMipWorker {
 public:
  struct SepaStatistics {
    SepaStatistics() : numNeighbourhoodQueries(0), sepa_lp_iterations(0) {}

    int64_t numNeighbourhoodQueries;
    int64_t sepa_lp_iterations;
  };
  struct HeurStatistics {
    HeurStatistics()
        : total_repair_lp(0),
          total_repair_lp_feasible(0),
          total_repair_lp_iterations(0),
          lp_iterations(0) {
      successObservations = 0;
      numSuccessObservations = 0;
      infeasObservations = 0;
      numInfeasObservations = 0;
      max_submip_level = 0;
      termination_status_ = HighsModelStatus::kNotset;
    }

    int64_t total_repair_lp;
    int64_t total_repair_lp_feasible;
    int64_t total_repair_lp_iterations;
    int64_t lp_iterations;

    double successObservations;
    HighsInt numSuccessObservations;
    double infeasObservations;
    HighsInt numInfeasObservations;
    HighsInt max_submip_level;
    HighsModelStatus termination_status_;
  };
  const HighsMipSolver& mipsolver_;
  const HighsMipSolverData& mipdata_;

  HighsLpRelaxation* lp_;
  HighsDomain* globaldom_;
  HighsCutPool* cutpool_;
  HighsConflictPool* conflictpool_;
  HighsPseudocost* pseudocost_;

  std::unique_ptr<HighsSearch> search_ptr_;
  std::unique_ptr<HighsSeparation> sepa_ptr_;

  HighsNodeQueue nodequeue;

  const HighsMipSolver& getMipSolver() const;

  double upper_bound;
  double upper_limit;
  double optimality_limit;

  std::vector<std::tuple<std::vector<double>, double, int>> solutions_;

  bool heuristics_allowed;

  HeurStatistics heur_stats;
  SepaStatistics sepa_stats;

  HighsRandom randgen;

  HighsMipWorker(const HighsMipSolver& mipsolver, HighsLpRelaxation* lp,
                 HighsDomain* domain, HighsCutPool* cutpool,
                 HighsConflictPool* conflictpool, HighsPseudocost* pseudocost);

  ~HighsMipWorker() {
    search_ptr_.reset();
    sepa_ptr_.reset();
  }

  const bool checkLimits(int64_t nodeOffset = 0) const;

  void resetSearch();

  void resetSepa();

  HighsDomain& getGlobalDomain() const { return *globaldom_; };

  HighsPseudocost& getPseudocost() const { return *pseudocost_; };

  HighsConflictPool& getConflictPool() const { return *conflictpool_; };

  HighsCutPool& getCutPool() const { return *cutpool_; };

  HighsLpRelaxation& getLpRelaxation() const { return *lp_; };

  bool addIncumbent(const std::vector<double>& sol, double solobj,
                    int solution_source);

  std::pair<bool, double> transformNewIntegerFeasibleSolution(
      const std::vector<double>& sol) const;

  bool trySolution(const std::vector<double>& solution,
                   const int solution_source);

  void setAllowHeuristics(const bool allowed) { heuristics_allowed = allowed; }

  bool getAllowHeuristics() const { return heuristics_allowed; }
};

#endif
