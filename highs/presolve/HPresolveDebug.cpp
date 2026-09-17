/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HPresolve.h"

// #include <algorithm>
// #include <atomic>
// #include <cmath>
// #include <limits>
//
// #include "../extern/pdqsort/pdqsort.h"
#include "Highs.h"
// #include "io/HighsIO.h"
// #include "lp_data/HConst.h"
// #include "lp_data/HStruct.h"
// #include "lp_data/HighsLpUtils.h"
// #include "lp_data/HighsModelUtils.h"
// #include "lp_data/HighsSolution.h"
// #include "mip/HighsCliqueTable.h"
// #include "mip/HighsImplications.h"
// #include "mip/HighsMipSolverData.h"
// #include "mip/HighsObjectiveFunction.h"
// #include "mip/MipTimer.h"
#include "presolve/HighsPostsolveStack.h"
// #include "presolve/PresolveTimer.h"
#include "test_kkt/DevKkt.h"
// #include "util/HFactor.h"
// #include "util/HighsCDouble.h"
// #include "util/HighsIntegers.h"
// #include "util/HighsLinearSumBounds.h"
// #include "util/HighsMemoryAllocation.h"
// #include "util/HighsSplay.h"
// #include "util/HighsUtils.h"

namespace presolve {
/*
void HPresolve::debug(const HighsLp& lp, const HighsOptions& options) {
HighsSolution reducedsol;
HighsBasis reducedbasis;

HighsSolution sol;
HighsBasis basis;

HighsLp model = lp;
model.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);

HighsPostsolveStack postsolve_stack;
postsolve_stack.initializeIndexMaps(lp.num_row_, lp.num_col_);
{
  HPresolve presolve;
  presolve.setInput(model, options, options.presolve_reduction_limit);
  if (presolve.run(postsolve_stack) != HighsModelStatus::kNotset) return;
  Highs highs;
  highs.passModel(model);
  highs.passOptions(options);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.run();
  if (highs.getModelStatus() != HighsModelStatus::kOptimal) return;
  reducedsol = highs.getSolution();
  reducedbasis = highs.getBasis();
}
model = lp;
sol = reducedsol;
basis = reducedbasis;
postsolve_stack.undo(options, sol, basis);
refineBasis(lp, sol, basis);
calculateRowValuesQuad(model, sol);
#if 0
Highs highs;
highs.passModel(model);
highs.passOptions(options);
highs.setSolution(sol);
basis.debug_origin_name = "HPresolve::debug";
highs.setBasis(basis);
highs.run();
return;
#endif
std::vector<HighsInt> flagCol(lp.num_col_, 1);
std::vector<HighsInt> flagRow(lp.num_row_, 1);
std::vector<HighsInt> Aend;
std::vector<HighsInt> ARstart;
std::vector<HighsInt> ARindex;
std::vector<double> ARvalue;
dev_kkt_check::KktInfo kktinfo = dev_kkt_check::initInfo();
Aend.assign(model.a_matrix_.start_.begin() + 1, model.a_matrix_.start_.end());
highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                     model.a_matrix_.index_, model.a_matrix_.value_, ARstart,
                     ARindex, ARvalue);
dev_kkt_check::State state(
    model.num_col_, model.num_row_, model.a_matrix_.start_, Aend,
    model.a_matrix_.index_, model.a_matrix_.value_, ARstart, ARindex, ARvalue,
    model.col_cost_, model.col_lower_, model.col_upper_, model.row_lower_,
    model.row_upper_, flagCol, flagRow, sol.col_value, sol.col_dual,
    sol.row_value, sol.row_dual, basis.col_status, basis.row_status);
bool checkResult = dev_kkt_check::checkKkt(state, kktinfo);
if (checkResult && kktinfo.pass_bfs) {
  printf("kkt check of postsolved solution and basis passed\n");
  return;
}
size_t good = postsolve_stack.numReductions();
size_t bad = 0;
size_t reductionLim = (good + bad) / 2;

// good = 1734357, bad = 1734289;
// good = 1050606, bad = 1050605;
// good = 1811527, bad = 1811526;
// reductionLim = bad;
do {
  model = lp;
  model.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);

  {
    HPresolve presolve;
    presolve.setInput(model, options, options.presolve_reduction_limit);
    presolve.computeIntermediateMatrix(flagRow, flagCol, reductionLim);
  }
#if 1
  model = lp;
  model.integrality_.assign(lp.num_col_, HighsVarType::kContinuous);
  HPresolve presolve;
  presolve.setInput(model, options, options.presolve_reduction_limit);
  HighsPostsolveStack tmp;
  tmp.initializeIndexMaps(model.num_row_, model.num_col_);
  presolve.setReductionLimit(reductionLim);
  presolve.run(tmp);

  sol = reducedsol;
  basis = reducedbasis;
  postsolve_stack.undo(options, sol, basis, tmp.numReductions());

  HighsBasis temp_basis;
  HighsSolution temp_sol;
  temp_basis.col_status.resize(model.num_col_);
  temp_sol.col_dual.resize(model.num_col_);
  temp_sol.col_value.resize(model.num_col_);
  for (HighsInt i = 0; i != model.num_col_; ++i) {
    temp_sol.col_dual[i] = sol.col_dual[tmp.getOrigColIndex(i)];
    temp_sol.col_value[i] = sol.col_value[tmp.getOrigColIndex(i)];
    temp_basis.col_status[i] = basis.col_status[tmp.getOrigColIndex(i)];
  }

  temp_basis.row_status.resize(model.num_row_);
  temp_sol.row_dual.resize(model.num_row_);
  for (HighsInt i = 0; i != model.num_row_; ++i) {
    temp_sol.row_dual[i] = sol.row_dual[tmp.getOrigRowIndex(i)];
    temp_basis.row_status[i] = basis.row_status[tmp.getOrigRowIndex(i)];
  }
  temp_sol.row_value.resize(model.num_row_);
  calculateRowValuesQuad(model, sol);
  temp_basis.valid = true;
  temp_basis.useful = true;
  refineBasis(model, temp_sol, temp_basis);
  Highs highs;
  highs.passOptions(options);
  highs.passModel(model);
  temp_basis.debug_origin_name = "HPresolve::debug";
  highs.setBasis(temp_basis);
  // highs.writeModel("model.mps");
  // highs.writeBasis("bad.bas");
  highs.run();
  printf("simplex iterations with postsolved basis: %" HIGHSINT_FORMAT "\n",
         highs.getInfo().simplex_iteration_count);
  checkResult = highs.getInfo().simplex_iteration_count == 0;
#else

  if (reductionLim == good) break;

  Aend.assign(model.a_matrix_.start_.begin() + 1,
              model.a_matrix_.start_.end());
  highsSparseTranspose(model.num_row_, model.num_col_, model.a_matrix_.start_,
                       model.a_matrix_.index_, model.a_matrix_.value_,
                       ARstart, ARindex, ARvalue);
  sol = reducedsol;
  basis = reducedbasis;
  postsolve_stack.undo(options, sol, basis, reductionLim);

  calculateRowValuesQuad(model, sol);
  kktinfo = dev_kkt_check::initInfo();
  checkResult = dev_kkt_check::checkKkt(state, kktinfo);
  checkResult = checkResult && kktinfo.pass_bfs;
#endif
  if (bad == good - 1) break;

  if (checkResult) {
    good = reductionLim;
  } else {
    bad = reductionLim;
  }
  reductionLim = (bad + good) / 2;
  printf("binary search ongoing: good=%zu, bad=%zu\n", good, bad);
} while (true);

printf("binary search finished: good=%zu, bad=%zu\n", good, bad);
assert(false);
}

void HPresolve::computeIntermediateMatrix(std::vector<HighsInt>& flagRow,
                                        std::vector<HighsInt>& flagCol,
                                        size_t& numreductions) {
shrinkProblemEnabled = false;
HighsPostsolveStack stack;
stack.initializeIndexMaps(flagRow.size(), flagCol.size());
setReductionLimit(numreductions);
presolve(stack);
numreductions = stack.numReductions();

toCSC(model->a_matrix_.value_, model->a_matrix_.index_,
      model->a_matrix_.start_);

for (HighsInt i = 0; i != model->num_row_; ++i) flagRow[i] = !rowDeleted[i];
for (HighsInt i = 0; i != model->num_col_; ++i) flagCol[i] = !colDeleted[i];
}
*/
}  // namespace presolve
