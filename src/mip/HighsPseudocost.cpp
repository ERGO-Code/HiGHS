#include "mip/HighsPseudocost.h"

#include "mip/HighsMipSolverData.h"

HighsPseudocost::HighsPseudocost(const HighsMipSolver& mipsolver)
    : pseudocostup(mipsolver.numCol()),
      pseudocostdown(mipsolver.numCol()),
      nsamplesup(mipsolver.numCol()),
      nsamplesdown(mipsolver.numCol()),
      inferencesup(mipsolver.numCol()),
      inferencesdown(mipsolver.numCol()),
      ninferencesup(mipsolver.numCol()),
      ninferencesdown(mipsolver.numCol()),
      ncutoffsup(mipsolver.numCol()),
      ncutoffsdown(mipsolver.numCol()),
      cost_total(0),
      inferences_total(0),
      nsamplestotal(0),
      ninferencestotal(0),
      ncutoffstotal(0),
      minreliable(mipsolver.options_mip_->mip_pscost_minreliable) {
  if (mipsolver.pscostinit != nullptr) {
    if (mipsolver.pscostinit->nsamplestotal != 0) {
      cost_total = mipsolver.pscostinit->cost_total;
      nsamplestotal = 1;
    }
    if (mipsolver.pscostinit->ninferencestotal != 0) {
      inferences_total = mipsolver.pscostinit->ninferencestotal;
      ninferencestotal = 1;
    }
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      int origCol = mipsolver.mipdata_->postSolveStack.getOrigColIndex(i);

      if (mipsolver.pscostinit->nsamplesup[origCol] != 0) {
        pseudocostup[i] = mipsolver.pscostinit->pseudocostup[origCol];
        nsamplesup[i] = 1;
      }

      if (mipsolver.pscostinit->nsamplesdown[origCol] != 0) {
        pseudocostdown[i] = mipsolver.pscostinit->pseudocostdown[origCol];
        nsamplesdown[i] = 1;
      }

      if (mipsolver.pscostinit->ninferencesup[origCol] != 0) {
        inferencesup[i] = mipsolver.pscostinit->inferencesup[origCol];
        ninferencesup[i] = 1;
      }

      if (mipsolver.pscostinit->ninferencesdown[origCol] != 0) {
        inferencesdown[i] = mipsolver.pscostinit->inferencesdown[origCol];
        ninferencesdown[i] = 1;
      }
    }
  }
}
