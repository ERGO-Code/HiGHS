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
  this->analyse_mip_time = mipsolver.analysis_.analyse_mip_time;
  this->clockIndex = -1;
  // Don't define clocks when analyse_mip_time is false - as will
  // generally be the case, and always so for sub-MIPs
  if (this->analyse_mip_time) 
    this->clockIndex = mipsolver.timer_.clock_def(name);
  if (this->analyse_mip_time) {
    printf("Defined clock %2d for %s\n", int(this->clockIndex), name);
  }
}

void HighsSeparator::run(HighsLpRelaxation& lpRelaxation,
                         HighsLpAggregator& lpAggregator,
                         HighsTransformedLp& transLp, HighsCutPool& cutpool) {
  ++numCalls;
  HighsInt currNumCuts = cutpool.getNumCuts();

  // Don't start/stop clocks when analyse_mip_time is false
  if (this->analyse_mip_time) 
    lpRelaxation.getMipSolver().timer_.start(clockIndex);
  separateLpSolution(lpRelaxation, lpAggregator, transLp, cutpool);
  if (this->analyse_mip_time) 
    lpRelaxation.getMipSolver().timer_.stop(clockIndex);

  numCutsFound += cutpool.getNumCuts() - currNumCuts;
}
