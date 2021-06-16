/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
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
    cost_total = mipsolver.pscostinit->cost_total;
    nsamplestotal = mipsolver.pscostinit->nsamplestotal;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      HighsInt origCol = mipsolver.mipdata_->postSolveStack.getOrigColIndex(i);

      pseudocostup[i] = mipsolver.pscostinit->pseudocostup[origCol];
      nsamplesup[i] = mipsolver.pscostinit->nsamplesup[origCol];
      pseudocostdown[i] = mipsolver.pscostinit->pseudocostdown[origCol];
      nsamplesdown[i] = mipsolver.pscostinit->nsamplesdown[origCol];
    }
  }
}

HighsPseudocostInitialization::HighsPseudocostInitialization(
    const HighsPseudocost& pscost, HighsInt maxCount)
    : pseudocostup(pscost.pseudocostup),
      pseudocostdown(pscost.pseudocostdown),
      nsamplesup(pscost.nsamplesup),
      nsamplesdown(pscost.nsamplesdown),
      cost_total(pscost.cost_total),
      nsamplestotal(std::min(int64_t{1}, pscost.nsamplestotal)) {
  HighsInt ncol = pseudocostup.size();
  for (HighsInt i = 0; i != ncol; ++i) {
    nsamplesup[i] = std::min(nsamplesup[i], maxCount);
    nsamplesdown[i] = std::min(nsamplesdown[i], maxCount);
  }
}

HighsPseudocostInitialization::HighsPseudocostInitialization(
    const HighsPseudocost& pscost, HighsInt maxCount,
    const presolve::HighsPostsolveStack& postsolveStack)
    : cost_total(pscost.cost_total),
      nsamplestotal(std::min(int64_t{1}, pscost.nsamplestotal)) {
  pseudocostup.resize(postsolveStack.getOrigNumCol());
  pseudocostdown.resize(postsolveStack.getOrigNumCol());
  nsamplesup.resize(postsolveStack.getOrigNumCol());
  nsamplesdown.resize(postsolveStack.getOrigNumCol());

  HighsInt ncols = pscost.pseudocostup.size();

  for (HighsInt i = 0; i != ncols; ++i) {
    pseudocostup[postsolveStack.getOrigColIndex(i)] = pscost.pseudocostup[i];
    pseudocostdown[postsolveStack.getOrigColIndex(i)] =
        pscost.pseudocostdown[i];
    nsamplesup[postsolveStack.getOrigColIndex(i)] =
        std::min(maxCount, pscost.nsamplesup[i]);
    nsamplesdown[postsolveStack.getOrigColIndex(i)] =
        std::min(maxCount, pscost.nsamplesdown[i]);
  }
}
