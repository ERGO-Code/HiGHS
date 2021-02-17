/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsSeparator.h
 * @brief Base class for separators
 *
 * @author Leona Gottwald
 */

#ifndef MIP_HIGHS_SEPARATOR_H_
#define MIP_HIGHS_SEPARATOR_H_

#include <cstddef>
#include <cstdint>

class HighsLpRelaxation;
class HighsTransformedLp;
class HighsCutPool;
class HighsLpAggregator;
class HighsMipSolver;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsSeparator {
 private:
  std::size_t numCutsFound;
  std::size_t numCalls;
  int clockIndex;

 public:
  HighsSeparator(const HighsMipSolver& mipsolver, const char* name,
                 const char* ch3_name);

  virtual void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                                  HighsLpAggregator& lpAggregator,
                                  HighsTransformedLp& transLp,
                                  HighsCutPool& cutpool) = 0;

  void run(HighsLpRelaxation& lpRelaxation, HighsLpAggregator& lpAggregator,
           HighsTransformedLp& transLp, HighsCutPool& cutpool);

  std::size_t getNumCutsFound() const { return numCutsFound; }

  std::size_t getNumCalls() const { return numCalls; }

  int getClockIndex() const { return clockIndex; }

  virtual ~HighsSeparator() {}
};

#endif