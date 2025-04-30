/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_PRIMAL_HEURISTICS_H_
#define HIGHS_PRIMAL_HEURISTICS_H_

#include <vector>

#include "lp_data/HStruct.h"
#include "lp_data/HighsLp.h"
#include "util/HighsRandom.h"

class HighsMipSolver;
class HighsMipWorker;

class HighsPrimalHeuristics {
 private:
  const HighsMipSolver& mipsolver;

  // HighsMipWorker& mipworker;
  // const HighsMipSolver& mipsolver;

  std::vector<HighsInt> intcols;

 public:
  struct Statistics {
    Statistics()
        : total_repair_lp(0),
          total_repair_lp_feasible(0),
          total_repair_lp_iterations(0),
          lp_iterations(0) {
      successObservations = 0;
      numSuccessObservations = 0;
      infeasObservations = 0;
      numInfeasObservations = 0;
    }

    size_t total_repair_lp;
    size_t total_repair_lp_feasible;
    size_t total_repair_lp_iterations;
    size_t lp_iterations;

    double successObservations;
    HighsInt numSuccessObservations;
    double infeasObservations;
    HighsInt numInfeasObservations;

    // still need to create in the mipworker
    // probably keep them separate

    // HighsRandom randgen;
  };

  HighsPrimalHeuristics(HighsMipSolver& mipsolver);
  // HighsPrimalHeuristics(HighsMipSolver& mipsolver);

  void setupIntCols();

  bool solveSubMip(HighsMipWorker& worker, const HighsLp& lp,
                   const HighsBasis& basis, double fixingRate,
                   std::vector<double> colLower, std::vector<double> colUpper,
                   HighsInt maxleaves, HighsInt maxnodes, HighsInt stallnodes);

  double determineTargetFixingRate(HighsMipWorker& worker);

  void rootReducedCost(HighsMipWorker& worker);

  // void RENS(const std::vector<double>& relaxationsol);
  void RENS(HighsMipWorker& worker, const std::vector<double>& relaxationsol);

  // void RINS(const std::vector<double>& relaxationsol);
  void RINS(HighsMipWorker& worker, const std::vector<double>& relaxationsol);

  void feasibilityPump(HighsMipWorker& worker);

  void centralRounding(HighsMipWorker& worker);

  void flushStatistics(HighsMipWorker& worker);

  bool tryRoundedPoint(HighsMipWorker& worker, const std::vector<double>& point,
                       const int solution_source);

  bool linesearchRounding(HighsMipWorker& worker,
                          const std::vector<double>& point1,
                          const std::vector<double>& point2,
                          const int solution_source);

  void randomizedRounding(HighsMipWorker& worker,
                          const std::vector<double>& relaxationsol);

  void shifting(HighsMipWorker& worker,
                const std::vector<double>& relaxationsol);

  void ziRound(HighsMipWorker& worker,
               const std::vector<double>& relaxationsol);
};

#endif
