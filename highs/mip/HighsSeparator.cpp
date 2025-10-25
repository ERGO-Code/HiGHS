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
                               const std::string& name)
    : numCutsFound(0), numCalls(0) {
  this->clockIndex = -1;
  // Don't get the clock index when analyse_mip_time is false - as
  // will generally be the case, and always so for sub-MIPs
  if (mipsolver.analysis_.analyse_mip_time) {
    this->clockIndex = mipsolver.analysis_.getSepaClockIndex(name);
    assert(this->clockIndex > 0);
  }
}

void HighsSeparator::run(HighsLpRelaxation& lpRelaxation,
                         HighsLpAggregator& lpAggregator,
                         HighsTransformedLp& transLp, HighsCutPool& cutpool) {
  ++numCalls;
  HighsInt currNumCuts = cutpool.getNumCuts();

  lpRelaxation.getMipSolver().analysis_.mipTimerStart(clockIndex);
  separateLpSolution(lpRelaxation, lpAggregator, transLp, cutpool);
  lpRelaxation.getMipSolver().analysis_.mipTimerStop(clockIndex);

  numCutsFound += cutpool.getNumCuts() - currNumCuts;
}
