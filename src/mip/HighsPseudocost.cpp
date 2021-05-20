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
    inferences_total = mipsolver.pscostinit->inferences_total;
    nsamplestotal = mipsolver.pscostinit->nsamplestotal;
    ninferencestotal = mipsolver.pscostinit->ninferencestotal;
    for (HighsInt i = 0; i != mipsolver.numCol(); ++i) {
      HighsInt origCol = mipsolver.mipdata_->postSolveStack.getOrigColIndex(i);

      pseudocostup[i] = mipsolver.pscostinit->pseudocostup[origCol];
      nsamplesup[i] = mipsolver.pscostinit->nsamplesup[origCol];
      pseudocostdown[i] = mipsolver.pscostinit->pseudocostdown[origCol];
      nsamplesdown[i] = mipsolver.pscostinit->nsamplesdown[origCol];
      inferencesup[i] = mipsolver.pscostinit->inferencesup[origCol];
      ninferencesup[i] = mipsolver.pscostinit->ninferencesup[origCol];
      inferencesdown[i] = mipsolver.pscostinit->inferencesdown[origCol];
      ninferencesdown[i] = mipsolver.pscostinit->ninferencesdown[origCol];
    }
  }
}

HighsPseudocostInitialization::HighsPseudocostInitialization(
    const HighsPseudocost& pscost, HighsInt maxCount)
    : pseudocostup(pscost.pseudocostup),
      pseudocostdown(pscost.pseudocostdown),
      nsamplesup(pscost.nsamplesup),
      nsamplesdown(pscost.nsamplesdown),
      inferencesup(pscost.inferencesup),
      inferencesdown(pscost.inferencesdown),
      ninferencesup(pscost.ninferencesup),
      ninferencesdown(pscost.ninferencesdown),
      cost_total(pscost.cost_total),
      inferences_total(pscost.inferences_total),
      nsamplestotal(std::min(int64_t{1}, pscost.nsamplestotal)),
      ninferencestotal(std::min(int64_t{1}, pscost.ninferencestotal)) {
  HighsInt ncol = pseudocostup.size();
  for (HighsInt i = 0; i != ncol; ++i) {
    nsamplesup[i] = std::min(nsamplesup[i], maxCount);
    nsamplesdown[i] = std::min(nsamplesdown[i], maxCount);
    ninferencesup[i] = std::min(ninferencesup[i], HighsInt{1});
    ninferencesdown[i] = std::min(ninferencesdown[i], HighsInt{1});
  }
}

HighsPseudocostInitialization::HighsPseudocostInitialization(
    const HighsPseudocost& pscost, HighsInt maxCount,
    const presolve::HighsPostsolveStack& postsolveStack)
    : cost_total(pscost.cost_total),
      inferences_total(pscost.inferences_total),
      nsamplestotal(std::min(int64_t{1}, pscost.nsamplestotal)),
      ninferencestotal(std::min(int64_t{1}, pscost.ninferencestotal)) {
  pseudocostup.resize(postsolveStack.getOrigNumCol());
  pseudocostdown.resize(postsolveStack.getOrigNumCol());
  nsamplesup.resize(postsolveStack.getOrigNumCol());
  nsamplesdown.resize(postsolveStack.getOrigNumCol());
  inferencesup.resize(postsolveStack.getOrigNumCol());
  inferencesdown.resize(postsolveStack.getOrigNumCol());
  ninferencesup.resize(postsolveStack.getOrigNumCol());
  ninferencesdown.resize(postsolveStack.getOrigNumCol());

  HighsInt ncols = pscost.pseudocostup.size();

  for (HighsInt i = 0; i != ncols; ++i) {
    pseudocostup[postsolveStack.getOrigColIndex(i)] = pscost.pseudocostup[i];
    pseudocostdown[postsolveStack.getOrigColIndex(i)] =
        pscost.pseudocostdown[i];
    nsamplesup[postsolveStack.getOrigColIndex(i)] =
        std::min(maxCount, pscost.nsamplesup[i]);
    nsamplesdown[postsolveStack.getOrigColIndex(i)] =
        std::min(maxCount, pscost.nsamplesdown[i]);
    inferencesup[postsolveStack.getOrigColIndex(i)] = pscost.inferencesup[i];
    inferencesdown[postsolveStack.getOrigColIndex(i)] =
        pscost.inferencesdown[i];
    ninferencesup[postsolveStack.getOrigColIndex(i)] = 1;
    ninferencesdown[postsolveStack.getOrigColIndex(i)] = 1;
  }
}
