/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsSeparator.h"

#include <string>

#include "mip/HighsCutPool.h"
#include "mip/HighsLpRelaxation.h"
#include "mip/HighsMipSolver.h"

HighsSeparator::HighsSeparator(const HighsMipSolver& mipsolver,
                               const char* name)
    : numCutsFound(0), numCalls(0) {
  clockIndex = mipsolver.timer_.clock_def(name);
}

void HighsSeparator::run(HighsLpRelaxation& lpRelaxation,
                         HighsLpAggregator& lpAggregator,
                         HighsTransformedLp& transLp, HighsCutPool& cutpool) {
  ++numCalls;
  HighsInt currNumCuts = cutpool.getNumCuts();

  lpRelaxation.getMipSolver().timer_.start(clockIndex);
  separateLpSolution(lpRelaxation, lpAggregator, transLp, cutpool);
  lpRelaxation.getMipSolver().timer_.stop(clockIndex);

  numCutsFound += cutpool.getNumCuts() - currNumCuts;
}
