#include "mip/HighsPseudocost.h"

#include "mip/HighsMipSolverData.h"

HighsPseudocost::HighsPseudocost(const HighsMipSolver& mipsolver)
    : pseudocostup(mipsolver.numCol()),
      pseudocostdown(mipsolver.numCol()),
      nsamplesup(mipsolver.numCol()),
      nsamplesdown(mipsolver.numCol()),
      ncutoffsup(mipsolver.numCol()),
      ncutoffsdown(mipsolver.numCol()),
      cost_total(0),
      nsamplestotal(0),
      ncutoffstotal(0),
      minreliable(mipsolver.options_mip_->mip_pscost_minreliable) {
  if (mipsolver.pscostinit != nullptr) {
    if (mipsolver.pscostinit->nsamplestotal != 0) {
      cost_total = mipsolver.pscostinit->cost_total /
                   mipsolver.pscostinit->nsamplestotal;
      nsamplestotal = 1;
    }
    ncutoffstotal = mipsolver.pscostinit->ncutoffstotal;
    for (int i = 0; i != mipsolver.numCol(); ++i) {
      int origCol = mipsolver.mipdata_->postSolveStack.getOrigColIndex(i);
      if (mipsolver.pscostinit->nsamplesup[origCol] != 0) {
        pseudocostup[i] = mipsolver.pscostinit->pseudocostup[origCol] /
                          mipsolver.pscostinit->nsamplesup[origCol];
        nsamplesup[i] = 1;
      }
      if (mipsolver.pscostinit->nsamplesdown[origCol] != 0) {
        pseudocostdown[i] = mipsolver.pscostinit->pseudocostdown[origCol] /
                            mipsolver.pscostinit->nsamplesdown[origCol];
        nsamplesdown[i] = 1;
      }

      ncutoffsup[i] = mipsolver.pscostinit->ncutoffsup[origCol];
      ncutoffsdown[i] = mipsolver.pscostinit->ncutoffsdown[origCol];
    }
  }
}