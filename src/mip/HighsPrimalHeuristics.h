/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_PRIMAL_HEURISTICS_H_
#define HIGHS_PRIMAL_HEURISTICS_H_

#include <random>
#include <vector>

class HighsMipSolver;

class HighsPrimalHeuristics {
 private:
  HighsMipSolver& mipsolver;
  size_t lp_iterations;
  std::mt19937 randgen;

 public:
  HighsPrimalHeuristics(HighsMipSolver& mipsolver);

  void solveSubMip(std::vector<double> colLower, std::vector<double> colUpper,
                   int maxleaves, int maxnodes);

  void RENS(const std::vector<double>& relaxationsol);

  void RINS(const std::vector<double>& relaxationsol);

  void feasibilityPump();

  void centralRounding();

  void flushStatistics();

  void clique();

  bool tryRoundedPoint(const std::vector<double>& point, char source);

  bool linesearchRounding(const std::vector<double>& point1,
                          const std::vector<double>& point2, char source);

  void randomizedRounding(const std::vector<double>& relaxationsol);
};

#endif
