/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_PRIMAL_HEURISTICS_H_
#define HIGHS_PRIMAL_HEURISTICS_H_

#include <vector>

#include "Highs.h"
//#include "lp_data/HStruct.h"
//#include "lp_data/HighsLp.h"
//#include "util/HighsRandom.h"

struct TrivialHeuristicRecord {
  HighsInt not_run = 0;
  HighsInt cannot_run = 0;
  HighsInt fail = 0;
  HighsInt feasible = 0;
  HighsInt improvement = 0;
};

enum TrivialHeuristic {
  kTrivialHeuristicZero = 0,
  kTrivialHeuristicUpper,
  kTrivialHeuristicLower,
  kTrivialHeuristicLock,
  kTrivialHeuristicCount
};

struct TrivialHeuristicData {
  std::vector<TrivialHeuristicRecord> record;
};

class HighsMipSolver;

class HighsPrimalHeuristics {
 private:
  HighsMipSolver& mipsolver;
  size_t lp_iterations;

  double successObservations;
  HighsInt numSuccessObservations;
  double infeasObservations;
  HighsInt numInfeasObservations;

  HighsRandom randgen;

  std::vector<HighsInt> intcols;
  TrivialHeuristicData trivial_heuristics_statistics_;

 public:
  HighsPrimalHeuristics(HighsMipSolver& mipsolver);

  void setupIntCols();

  bool solveSubMip(const HighsLp& lp, const HighsBasis& basis,
                   double fixingRate, std::vector<double> colLower,
                   std::vector<double> colUpper, HighsInt maxleaves,
                   HighsInt maxnodes, HighsInt stallnodes);

  double determineTargetFixingRate();

  void rootReducedCost();

  void RENS(const std::vector<double>& relaxationsol);

  void RINS(const std::vector<double>& relaxationsol);

  void feasibilityPump();

  void centralRounding();

  void flushStatistics();

  bool tryRoundedPoint(const std::vector<double>& point, char source);

  bool linesearchRounding(const std::vector<double>& point1,
                          const std::vector<double>& point2, char source);

  void randomizedRounding(const std::vector<double>& relaxationsol);

  void trivial();
  void runTrivial(Highs* highs = nullptr);
  void initialiseLocalTrivialHeuristicsStatistics();
  void downCopyLocalTrivialHeuristicsStatistics(
      const TrivialHeuristicData& from_statistics);
  void upCopyLocalTrivialHeuristicsStatistics(
      TrivialHeuristicData& to_statistics);
};

void initialiseTrivialHeuristicsStatistics(TrivialHeuristicData& statistics);
void copyTrivialHeuristicsStatistics(
    const TrivialHeuristicData& from_statistics,
    TrivialHeuristicData& to_statistics);
void reportTrivialHeuristicsStatistics(
    const HighsLogOptions& log_options,
    const TrivialHeuristicData& to_statistics);

#endif
