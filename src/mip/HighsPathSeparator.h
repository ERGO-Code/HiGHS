/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsPathSeparator.h
 * @brief Class for separating cuts from heuristically aggregating rows from the
 * LP to indetify path's in a network
 *
 * @author Leona Gottwald
 */

#ifndef MIP_HIGHS_PATH_SEPARATOR_H_
#define MIP_HIGHS_PATH_SEPARATOR_H_

#include <random>

#include "mip/HighsSeparator.h"

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsPathSeparator : public HighsSeparator {
 private:
  std::mt19937 randgen;

 public:
  void separateLpSolution(HighsLpRelaxation& lpRelaxation,
                          HighsLpAggregator& lpAggregator,
                          HighsTransformedLp& transLp,
                          HighsCutPool& cutpool) override;

  HighsPathSeparator(const HighsMipSolver& mipsolver)
      : HighsSeparator(mipsolver, "PathAggr sepa", "Agg") {}
};

#endif