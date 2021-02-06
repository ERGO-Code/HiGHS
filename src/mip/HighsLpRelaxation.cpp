#include "mip/HighsLpRelaxation.h"

#include <unordered_map>

#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsPseudocost.h"
#include "util/HighsCDouble.h"

void HighsLpRelaxation::LpRow::get(const HighsMipSolver& mipsolver, int& len,
                                   const int*& inds,
                                   const double*& vals) const {
  switch (origin) {
    case kCutPool:
      mipsolver.mipdata_->cutpool.getCut(index, len, inds, vals);
      break;
    case kModel:
      mipsolver.mipdata_->getRow(index, len, inds, vals);
  };
}

int HighsLpRelaxation::LpRow::getRowLen(const HighsMipSolver& mipsolver) const {
  switch (origin) {
    case kCutPool:
      return mipsolver.mipdata_->cutpool.getRowLength(index);
    case kModel:
      return mipsolver.mipdata_->ARstart_[index + 1] -
             mipsolver.mipdata_->ARstart_[index];
  };

  assert(false);
  return -1;
}

bool HighsLpRelaxation::LpRow::isIntegral(
    const HighsMipSolver& mipsolver) const {
  switch (origin) {
    case kCutPool:
      return mipsolver.mipdata_->cutpool.cutIsIntegral(index);
    case kModel:
      return mipsolver.mipdata_->rowintegral[index];
  };

  assert(false);
  return false;
}

double HighsLpRelaxation::LpRow::getMaxAbsVal(
    const HighsMipSolver& mipsolver) const {
  switch (origin) {
    case kCutPool:
      return mipsolver.mipdata_->cutpool.getMaxAbsCutCoef(index);
    case kModel:
      return mipsolver.mipdata_->maxAbsRowCoef[index];
  };

  assert(false);
  return 0.0;
}

double HighsLpRelaxation::slackLower(int row) const {
  double rowlower = rowLower(row);
  if (rowlower != -HIGHS_CONST_INF) return rowlower;

  switch (lprows[row].origin) {
    case LpRow::kCutPool:
      return mipsolver.mipdata_->domain.getMinCutActivity(
          mipsolver.mipdata_->cutpool, lprows[row].index);
    case LpRow::kModel:
      return mipsolver.mipdata_->domain.getMinActivity(lprows[row].index);
  };

  assert(false);
  return -HIGHS_CONST_INF;
}

double HighsLpRelaxation::slackUpper(int row) const {
  double rowupper = rowUpper(row);
  switch (lprows[row].origin) {
    case LpRow::kCutPool:
      return rowupper;
    case LpRow::kModel:
      if (rowupper != HIGHS_CONST_INF) return rowupper;
      return mipsolver.mipdata_->domain.getMaxActivity(lprows[row].index);
  };

  assert(false);
  return HIGHS_CONST_INF;
}

HighsLpRelaxation::HighsLpRelaxation(const HighsMipSolver& mipsolver)
    : mipsolver(mipsolver) {
  lpsolver.setHighsLogfile();
  lpsolver.setHighsOutput();
  lpsolver.setHighsOptionValue(
      "primal_feasibility_tolerance",
      mipsolver.options_mip_->mip_feasibility_tolerance);
  lpsolver.setHighsOptionValue(
      "dual_feasibility_tolerance",
      mipsolver.options_mip_->mip_feasibility_tolerance * 0.1);
  status = Status::NotSet;
  numlpiters = 0;
  epochs = 0;
  currentbasisstored = false;
}

HighsLpRelaxation::HighsLpRelaxation(const HighsLpRelaxation& other)
    : mipsolver(other.mipsolver),
      lprows(other.lprows),
      fractionalints(other.fractionalints),
      objective(other.objective),
      basischeckpoint(other.basischeckpoint),
      currentbasisstored(other.currentbasisstored) {
  lpsolver.setHighsLogfile();
  lpsolver.setHighsOutput();
  lpsolver.passHighsOptions(other.lpsolver.getHighsOptions());
  lpsolver.passModel(other.lpsolver.getLp());
  lpsolver.setBasis(other.lpsolver.getBasis());
  numlpiters = 0;
  epochs = 0;
}

void HighsLpRelaxation::loadModel() {
  HighsLp lpmodel = *mipsolver.model_;
  lpmodel.colLower_ = mipsolver.mipdata_->domain.colLower_;
  lpmodel.colUpper_ = mipsolver.mipdata_->domain.colUpper_;
  lprows.reserve(lpmodel.numRow_);
  for (int i = 0; i != lpmodel.numRow_; ++i) lprows.push_back(LpRow::model(i));
  lpmodel.integrality_.clear();
  lpsolver.passModel(std::move(lpmodel));
  mipsolver.mipdata_->domain.clearChangedCols();
}

double HighsLpRelaxation::computeBestEstimate(const HighsPseudocost& ps) const {
  HighsCDouble estimate = objective;

  for (const std::pair<int, double>& f : fractionalints)
    estimate += std::min(ps.getPseudocostUp(f.first, f.second),
                         ps.getPseudocostDown(f.first, f.second));

  return double(estimate);
}

void HighsLpRelaxation::addCuts(HighsCutSet& cutset) {
  int numcuts = cutset.numCuts();
  assert(lpsolver.getLp().numRow_ == (int)lpsolver.getLp().rowLower_.size());
  assert(lpsolver.getLp().numRow_ == (int)lprows.size());
  if (numcuts > 0) {
    status = Status::NotSet;
    currentbasisstored = false;
    basischeckpoint.reset();

    lprows.reserve(lprows.size() + numcuts);
    for (int i = 0; i != numcuts; ++i)
      lprows.push_back(LpRow::cut(cutset.cutindices[i]));

    lpsolver.addRows(numcuts, cutset.lower_.data(), cutset.upper_.data(),
                     cutset.ARvalue_.size(), cutset.ARstart_.data(),
                     cutset.ARindex_.data(), cutset.ARvalue_.data());
    assert(lpsolver.getLp().numRow_ == (int)lpsolver.getLp().rowLower_.size());
    cutset.clear();
  }
}

void HighsLpRelaxation::removeObsoleteRows() {
  int nlprows = numRows();
  int nummodelrows = getNumModelRows();
  std::vector<int> deletemask;

  int ndelcuts = 0;
  for (int i = nummodelrows; i != nlprows; ++i) {
    assert(lprows[i].origin == LpRow::Origin::kCutPool);
    if (lpsolver.getBasis().row_status[i] == HighsBasisStatus::BASIC) {
      if (ndelcuts == 0) deletemask.resize(nlprows);
      ++ndelcuts;
      deletemask[i] = 1;
      mipsolver.mipdata_->cutpool.lpCutRemoved(lprows[i].index);
    }
  }

  removeCuts(ndelcuts, deletemask);
}

void HighsLpRelaxation::removeCuts(int ndelcuts, std::vector<int>& deletemask) {
  assert(lpsolver.getLp().numRow_ == (int)lpsolver.getLp().rowLower_.size());
  if (ndelcuts > 0) {
    HighsBasis basis = lpsolver.getBasis();
    int nlprows = lpsolver.getNumRows();
    lpsolver.deleteRows(deletemask.data());
    for (int i = mipsolver.numRow(); i != nlprows; ++i) {
      if (deletemask[i] >= 0) {
        lprows[deletemask[i]] = lprows[i];
        basis.row_status[deletemask[i]] = basis.row_status[i];
      }
    }

    assert(lpsolver.getLp().numRow_ == (int)lpsolver.getLp().rowLower_.size());

    basis.row_status.resize(basis.row_status.size() - ndelcuts);
    lprows.resize(lprows.size() - ndelcuts);

    assert(lpsolver.getLp().numRow_ == (int)lprows.size());
    lpsolver.setBasis(basis);
    lpsolver.run();
  }
}

void HighsLpRelaxation::removeCuts() {
  assert(lpsolver.getLp().numRow_ == (int)lpsolver.getLp().rowLower_.size());
  int nlprows = lpsolver.getNumRows();
  int modelrows = mipsolver.numRow();

  lpsolver.deleteRows(modelrows, nlprows - 1);
  for (int i = modelrows; i != nlprows; ++i) {
    if (lprows[i].origin == LpRow::Origin::kCutPool)
      mipsolver.mipdata_->cutpool.lpCutRemoved(lprows[i].index);
  }
  lprows.resize(modelrows);
  assert(lpsolver.getLp().numRow_ == (int)lpsolver.getLp().rowLower_.size());
}

void HighsLpRelaxation::performAging() {
  assert(lpsolver.getLp().numRow_ == (int)lpsolver.getLp().rowLower_.size());

  size_t agelimit = mipsolver.options_mip_->mip_lp_age_limit;

  ++epochs;
  if (epochs % std::max(size_t(agelimit) / 2u, size_t(2)) != 0)
    agelimit = HIGHS_CONST_I_INF;
  else if (epochs < agelimit)
    agelimit = epochs;

  int nlprows = numRows();
  int nummodelrows = getNumModelRows();
  std::vector<int> deletemask;

  int ndelcuts = 0;
  for (int i = nummodelrows; i != nlprows; ++i) {
    assert(lprows[i].origin == LpRow::Origin::kCutPool);
    if (lpsolver.getBasis().row_status[i] == HighsBasisStatus::BASIC) {
      if (mipsolver.mipdata_->cutpool.ageLpCut(lprows[i].index, agelimit)) {
        if (ndelcuts == 0) deletemask.resize(nlprows);
        ++ndelcuts;
        deletemask[i] = 1;
      }
    } else {
      mipsolver.mipdata_->cutpool.resetAge(lprows[i].index);
    }
  }

  removeCuts(ndelcuts, deletemask);
}

void HighsLpRelaxation::flushDomain(HighsDomain& domain, bool continuous) {
  if (!domain.getChangedCols().empty()) {
    currentbasisstored = false;
    for (int col : domain.getChangedCols()) {
      if (!continuous &&
          mipsolver.variableType(col) == HighsVarType::CONTINUOUS)
        continue;
      lpsolver.changeColBounds(col, domain.colLower_[col],
                               domain.colUpper_[col]);
    }

    domain.clearChangedCols();
  }
}

bool HighsLpRelaxation::computeDualProof(const HighsDomain& globaldomain,
                                         double upperbound,
                                         std::vector<int>& inds,
                                         std::vector<double>& vals,
                                         double& rhs) const {
  std::vector<double> row_dual = lpsolver.getSolution().row_dual;

  const HighsLp& lp = lpsolver.getLp();

  assert(std::isfinite(upperbound));
  HighsCDouble upper = upperbound;

  for (int i = 0; i != lp.numRow_; ++i) {
    if (row_dual[i] < 0) {
      if (lp.rowLower_[i] != -HIGHS_CONST_INF)
        upper += row_dual[i] * lp.rowLower_[i];
      else
        row_dual[i] = 0;
    } else if (row_dual[i] > 0) {
      if (lp.rowUpper_[i] != HIGHS_CONST_INF)
        upper += row_dual[i] * lp.rowUpper_[i];
      else
        row_dual[i] = 0;
    }
  }

  inds.clear();
  vals.clear();
  for (int i = 0; i != lp.numCol_; ++i) {
    int start = lp.Astart_[i];
    int end = lp.Astart_[i + 1];

    HighsCDouble sum = lp.colCost_[i];

    for (int j = start; j != end; ++j) {
      if (row_dual[lp.Aindex_[j]] == 0) continue;
      sum += lp.Avalue_[j] * row_dual[lp.Aindex_[j]];
    }

    double val = double(sum);

    if (std::abs(val) <= mipsolver.options_mip_->small_matrix_value) continue;

    if (mipsolver.variableType(i) == HighsVarType::CONTINUOUS ||
        std::abs(val) <= mipsolver.mipdata_->feastol ||
        globaldomain.colLower_[i] == globaldomain.colUpper_[i]) {
      if (val < 0) {
        if (globaldomain.colUpper_[i] == HIGHS_CONST_INF) return false;
        upper -= val * globaldomain.colUpper_[i];
      } else {
        if (globaldomain.colLower_[i] == -HIGHS_CONST_INF) return false;

        upper -= val * globaldomain.colLower_[i];
      }

      continue;
    }

    vals.push_back(val);
    inds.push_back(i);
  }

  rhs = double(upper);
  assert(std::isfinite(rhs));
  globaldomain.tightenCoefficients(inds.data(), vals.data(), inds.size(), rhs);

  mipsolver.mipdata_->debugSolution.checkCut(inds.data(), vals.data(),
                                             inds.size(), rhs);

  return true;
}

void HighsLpRelaxation::storeDualInfProof() {
  assert(lpsolver.getModelStatus(true) == HighsModelStatus::PRIMAL_INFEASIBLE);

  int numrow = lpsolver.getNumRows();
  hasdualproof = false;
  lpsolver.getDualRay(hasdualproof);

  if (!hasdualproof) {
    HighsPrintMessage(mipsolver.options_mip_->output,
                      mipsolver.options_mip_->message_level, ML_VERBOSE,
                      "no dual ray stored\n");
    return;
  }

  dualproofinds.clear();
  dualproofvals.clear();
  dualproofrhs = HIGHS_CONST_INF;
  const HighsLp& lp = lpsolver.getLp();
  dualproofbuffer.resize(numrow);

  lpsolver.getDualRay(hasdualproof, dualproofbuffer.data());
  std::vector<double>& dualray = dualproofbuffer;

  HighsCDouble upper = 0.0;

  double maxval = 0;
  for (int i = 0; i != lp.numRow_; ++i)
    maxval = std::max(maxval, std::abs(dualray[i]));

  int expscal;
  std::frexp(maxval, &expscal);

  for (int i = 0; i != lp.numRow_; ++i) {
    dualray[i] = std::ldexp(dualray[i], -expscal);
    if (dualray[i] < 0) {
      if (lp.rowUpper_[i] == HIGHS_CONST_INF) dualray[i] = 0.0;
    } else if (dualray[i] > 0) {
      if (lp.rowLower_[i] == -HIGHS_CONST_INF) dualray[i] = 0.0;
    }
  }

  for (int i = 0; i != lp.numRow_; ++i) {
    if (dualray[i] < 0) {
      assert(lp.rowUpper_[i] != HIGHS_CONST_INF);
      upper -= dualray[i] * lp.rowUpper_[i];
    } else if (dualray[i] > 0) {
      assert(lp.rowLower_[i] != -HIGHS_CONST_INF);
      upper -= dualray[i] * lp.rowLower_[i];
    }
  }

  for (int i = 0; i != lp.numCol_; ++i) {
    int start = lp.Astart_[i];
    int end = lp.Astart_[i + 1];

    HighsCDouble sum = 0.0;

    for (int j = start; j != end; ++j) {
      if (dualray[lp.Aindex_[j]] == 0.0) continue;
      sum -= lp.Avalue_[j] * dualray[lp.Aindex_[j]];
    }

    double val = double(sum);

    if (std::abs(val) <= mipsolver.options_mip_->small_matrix_value) continue;

    if (mipsolver.variableType(i) == HighsVarType::CONTINUOUS ||
        std::abs(val) <= mipsolver.mipdata_->feastol ||
        mipsolver.mipdata_->domain.colLower_[i] ==
            mipsolver.mipdata_->domain.colUpper_[i]) {
      if (val < 0) {
        if (mipsolver.mipdata_->domain.colUpper_[i] == HIGHS_CONST_INF) return;
        upper -= val * mipsolver.mipdata_->domain.colUpper_[i];
      } else {
        if (mipsolver.mipdata_->domain.colLower_[i] == -HIGHS_CONST_INF) return;
        upper -= val * mipsolver.mipdata_->domain.colLower_[i];
      }

      continue;
    }

    dualproofvals.push_back(val);
    dualproofinds.push_back(i);
  }

  dualproofrhs = double(upper);
  mipsolver.mipdata_->domain.tightenCoefficients(
      dualproofinds.data(), dualproofvals.data(), dualproofinds.size(),
      dualproofrhs);

  mipsolver.mipdata_->debugSolution.checkCut(
      dualproofinds.data(), dualproofvals.data(), dualproofinds.size(),
      dualproofrhs);
}

void HighsLpRelaxation::storeDualUBProof() {
  dualproofinds.clear();
  dualproofvals.clear();
  dualproofrhs = HIGHS_CONST_INF;
  assert(lpsolver.getModelStatus(true) ==
         HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND);

  int numrow = lpsolver.getNumRows();
  bool hasdualray = false;
  lpsolver.getDualRay(hasdualray);

  if (!hasdualray) return;

  const HighsLp& lp = lpsolver.getLp();
  dualproofbuffer.resize(numrow);

  lpsolver.getDualRay(hasdualray, dualproofbuffer.data());
  std::vector<double>& dualray = dualproofbuffer;

  double scale = 0.0;

  for (int i = 0; i != lp.numRow_; ++i) {
    if (std::abs(dualray[i]) <= mipsolver.mipdata_->feastol) {
      dualray[i] = 0.0;
      continue;
    }

    if (scale * dualray[i] <= 0.0) {
      if (lp.rowUpper_[i] == HIGHS_CONST_INF) {
        if (scale == 0.0)
          scale = copysign(1.0, dualray[i]);
        else
          return;
      }
    }

    if (scale * dualray[i] >= 0.0) {
      if (lp.rowLower_[i] == -HIGHS_CONST_INF) {
        if (scale == 0.0)
          scale = -copysign(1.0, dualray[i]);
        else
          return;
      }
    }
  }

  if (scale == 0.0) scale = 1.0;

  assert(scale == 1.0);

  HighsCDouble upper =
      lpsolver.getHighsOptions().dual_objective_value_upper_bound;
  for (int i = 0; i != lp.numRow_; ++i) {
    if (dualray[i] == 0.0) continue;

    if (scale * dualray[i] < 0) {
      assert(lp.rowUpper_[i] != HIGHS_CONST_INF);
      upper -= scale * dualray[i] * lp.rowUpper_[i];
    } else {
      assert(lp.rowLower_[i] != -HIGHS_CONST_INF);
      upper -= scale * dualray[i] * lp.rowLower_[i];
    }
  }

  for (int i = 0; i != lp.numCol_; ++i) {
    int start = lp.Astart_[i];
    int end = lp.Astart_[i + 1];

    HighsCDouble sum = scale * mipsolver.colCost(i);

    for (int j = start; j != end; ++j) {
      if (dualray[lp.Aindex_[j]] == 0.0) continue;
      sum -= lp.Avalue_[j] * dualray[lp.Aindex_[j]];
    }

    double val = scale * double(sum);

    if (std::abs(val) <= 1e-12) continue;

    if (mipsolver.variableType(i) == HighsVarType::CONTINUOUS ||
        std::abs(val) < mipsolver.mipdata_->feastol ||
        mipsolver.mipdata_->domain.colLower_[i] ==
            mipsolver.mipdata_->domain.colUpper_[i]) {
      if (val < 0) {
        if (mipsolver.mipdata_->domain.colUpper_[i] == HIGHS_CONST_INF) return;
        upper -= val * mipsolver.mipdata_->domain.colUpper_[i];
      } else {
        if (mipsolver.mipdata_->domain.colLower_[i] == -HIGHS_CONST_INF) return;

        upper -= val * mipsolver.mipdata_->domain.colLower_[i];
      }

      continue;
    }

    dualproofvals.push_back(val);
    dualproofinds.push_back(i);
  }

  dualproofrhs = double(upper);
  mipsolver.mipdata_->domain.tightenCoefficients(
      dualproofinds.data(), dualproofvals.data(), dualproofinds.size(),
      dualproofrhs);

  mipsolver.mipdata_->debugSolution.checkCut(
      dualproofinds.data(), dualproofvals.data(), dualproofinds.size(),
      dualproofrhs);
}

bool HighsLpRelaxation::checkDualProof() const {
  if (!hasdualproof) return true;
  if (dualproofrhs == HIGHS_CONST_INF) return false;

  int len = dualproofinds.size();

  HighsCDouble viol = -dualproofrhs;

  const HighsLp& lp = lpsolver.getLp();

  for (int i = 0; i != len; ++i) {
    int col = dualproofinds[i];
    if (dualproofvals[i] > 0) {
      if (lp.colLower_[col] == -HIGHS_CONST_INF) return false;
      viol += dualproofvals[i] * lp.colLower_[col];
    } else {
      assert(dualproofvals[i] < 0);
      if (lp.colUpper_[col] == HIGHS_CONST_INF) return false;
      viol += dualproofvals[i] * lp.colUpper_[col];
    }
  }

  return viol > mipsolver.mipdata_->feastol;
}

bool HighsLpRelaxation::computeDualInfProof(const HighsDomain& globaldomain,
                                            std::vector<int>& inds,
                                            std::vector<double>& vals,
                                            double& rhs) {
  if (!hasdualproof) return false;

  assert(checkDualProof());

  inds = dualproofinds;
  vals = dualproofvals;
  rhs = dualproofrhs;
  return true;
}

void HighsLpRelaxation::recoverBasis() {
  if (basischeckpoint) lpsolver.setBasis(*basischeckpoint);
}

HighsLpRelaxation::Status HighsLpRelaxation::run(bool resolve_on_error) {
  HighsStatus callstatus;

  try {
    callstatus = lpsolver.run();
  } catch (const std::runtime_error&) {
    callstatus = HighsStatus::Error;
  }

  const HighsInfo& info = lpsolver.getHighsInfo();
  numlpiters += std::max(0, info.simplex_iteration_count);

  if (callstatus == HighsStatus::Error) {
    // first try to use the primal simplex solver starting from the last basis
    if (lpsolver.getHighsOptions().simplex_strategy == SIMPLEX_STRATEGY_DUAL) {
      lpsolver.setHighsOptionValue("simplex_strategy", SIMPLEX_STRATEGY_PRIMAL);
      recoverBasis();
      auto retval = run(false);
      lpsolver.setHighsOptionValue("simplex_strategy", SIMPLEX_STRATEGY_DUAL);

      return retval;
    }

    if (resolve_on_error) {
      // still an error: now try to solve with presolve from scratch
      lpsolver.clearSolver();
      lpsolver.setHighsOptionValue("simplex_strategy", SIMPLEX_STRATEGY_DUAL);
      lpsolver.setHighsOptionValue("presolve", "on");
      auto retval = run(false);
      lpsolver.setHighsOptionValue("presolve", "off");

      return retval;
    }

    recoverBasis();

    return Status::Error;
  }

  HighsModelStatus scaledmodelstatus = lpsolver.getModelStatus(true);
  switch (scaledmodelstatus) {
    case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      storeDualUBProof();
      if (checkDualProof()) return Status::Infeasible;

      return Status::Error;
    case HighsModelStatus::PRIMAL_DUAL_INFEASIBLE:
    case HighsModelStatus::PRIMAL_INFEASIBLE: {
      storeDualInfProof();
      if (checkDualProof()) {
        return Status::Infeasible;
      }
      hasdualproof = false;

      int scalestrategy = lpsolver.getHighsOptions().simplex_scale_strategy;
      if (scalestrategy != SIMPLEX_SCALE_STRATEGY_OFF) {
        lpsolver.setHighsOptionValue("simplex_scale_strategy",
                                     SIMPLEX_SCALE_STRATEGY_OFF);
        HighsBasis basis = lpsolver.getBasis();
        lpsolver.clearSolver();
        lpsolver.setBasis(basis);
        auto tmp = run(resolve_on_error);
        lpsolver.setHighsOptionValue("simplex_scale_strategy", scalestrategy);
        if (!scaledOptimal(tmp)) {
          lpsolver.clearSolver();
          lpsolver.setBasis(basis);
        }
        return tmp;
      }

      // trust the primal simplex result without scaling
      if (lpsolver.getModelStatus() == HighsModelStatus::PRIMAL_INFEASIBLE)
        return Status::Infeasible;

      // HighsLogMessage(mipsolver.options_mip_->logfile,
      //                 HighsMessageType::WARNING,
      //                 "LP failed to reliably determine infeasibility");

      // printf("error: unreliable infeasiblities, modelstatus = %d (scaled
      // %d)\n",
      //        (int)lpsolver.getModelStatus(),
      //        (int)lpsolver.getModelStatus(true));
      return Status::Error;
    }
    case HighsModelStatus::OPTIMAL:
      assert(info.max_primal_infeasibility >= 0);
      assert(info.max_dual_infeasibility >= 0);
      if (info.max_primal_infeasibility <= mipsolver.mipdata_->feastol &&
          info.max_dual_infeasibility <= mipsolver.mipdata_->feastol)
        return Status::Optimal;

      if (resolve_on_error) {
        // printf(
        //     "error: optimal with unscaled infeasibilities (primal:%g, "
        //     "dual:%g)\n",
        //     info.max_primal_infeasibility, info.max_dual_infeasibility);
        int scalestrategy = lpsolver.getHighsOptions().simplex_scale_strategy;
        if (scalestrategy != SIMPLEX_SCALE_STRATEGY_OFF) {
          lpsolver.setHighsOptionValue("simplex_scale_strategy",
                                       SIMPLEX_SCALE_STRATEGY_OFF);
          HighsBasis basis = lpsolver.getBasis();
          lpsolver.clearSolver();
          lpsolver.setBasis(basis);
          auto tmp = run(resolve_on_error);
          lpsolver.setHighsOptionValue("simplex_scale_strategy", scalestrategy);
          return tmp;
        }
      }

      if (info.max_primal_infeasibility <= mipsolver.mipdata_->feastol)
        return Status::UnscaledPrimalFeasible;

      if (info.max_dual_infeasibility <= mipsolver.mipdata_->feastol)
        return Status::UnscaledDualFeasible;

      return Status::UnscaledInfeasible;
    case HighsModelStatus::REACHED_ITERATION_LIMIT: {
      if (resolve_on_error) {
        // printf(
        //     "error: lpsolver reached iteration limit, resolving with basis "
        //     "from ipm\n");
        Highs ipm;
        ipm.passModel(lpsolver.getLp());
        ipm.setHighsOptionValue("solver", "ipm");
        ipm.setHighsOutput();
        ipm.setHighsLogfile();
        ipm.run();
        lpsolver.setBasis(ipm.getBasis());
        return run(false);
      }

      // printf("error: lpsolver reached iteration limit\n");
      return Status::Error;
    }
    // case HighsModelStatus::PRIMAL_DUAL_INFEASIBLE:
    // case HighsModelStatus::PRIMAL_INFEASIBLE:
    //  if (lpsolver.getModelStatus(false) == scaledmodelstatus)
    //    return Status::Infeasible;
    //  return Status::Error;
    default:
      // printf("error: lpsolver stopped with unexpected status %d\n",
      //        (int)scaledmodelstatus);
      HighsLogMessage(
          mipsolver.options_mip_->logfile, HighsMessageType::WARNING,
          "LP solved to unexpected status (%d)", (int)scaledmodelstatus);
      return Status::Error;
  }
}

HighsLpRelaxation::Status HighsLpRelaxation::resolveLp(HighsDomain* domain) {
  fractionalints.clear();

  bool solveagain;
  do {
    solveagain = false;
    if (domain) flushDomain(*domain);
    status = run();

    switch (status) {
      case Status::UnscaledInfeasible:
      case Status::UnscaledDualFeasible:
      case Status::UnscaledPrimalFeasible:
      case Status::Optimal: {
        std::unordered_map<int, std::pair<double, int>> fracints;
        const HighsSolution& sol = lpsolver.getSolution();

        HighsCDouble objsum = 0;
        bool roundable = true;

        for (int i : mipsolver.mipdata_->integral_cols) {
          // for the fractionality we assume that LP bounds are not violated
          // bounds that are violated by the unscaled LP are indicated by the
          // return status already
          double val = std::max(
              std::min(sol.col_value[i], lpsolver.getLp().colUpper_[i]),
              lpsolver.getLp().colLower_[i]);
          double intval = std::floor(val + 0.5);

          if (std::abs(val - intval) > mipsolver.mipdata_->feastol) {
            int col = i;
            if (roundable && mipsolver.mipdata_->uplocks[col] != 0 &&
                mipsolver.mipdata_->downlocks[col] != 0)
              roundable = false;

            const HighsCliqueTable::Substitution* subst =
                mipsolver.mipdata_->cliquetable.getSubstitution(col);
            while (subst != nullptr) {
              if (lpsolver.getLp().colLower_[subst->replace.col] ==
                  lpsolver.getLp().colUpper_[subst->replace.col]) {
                if (domain)
                  domain->fixCol(
                      col, subst->replace.weight(lpsolver.getLp().colLower_));
                else
                  break;
              }

              assert(lpsolver.getLp().colLower_[subst->replace.col] !=
                     lpsolver.getLp().colUpper_[subst->replace.col]);
              col = subst->replace.col;
              if (subst->replace.val == 0) val = 1.0 - val;

              subst = mipsolver.mipdata_->cliquetable.getSubstitution(col);
            }
            auto& pair = fracints[col];
            pair.first += val;
            pair.second += 1;
          }
        }

        if (domain && !domain->getChangedCols().empty()) {
          // printf("resolving due to fixings of substituted columns\n");
          solveagain = true;
          continue;
        }

        for (const auto& it : fracints) {
          fractionalints.emplace_back(
              it.first, it.second.first / (double)it.second.second);
        }

        if (roundable && !fractionalints.empty()) {
          std::vector<double> roundsol = sol.col_value;

          for (const std::pair<int, double>& fracint : fractionalints) {
            int col = fracint.first;

            if (mipsolver.mipdata_->uplocks[col] == 0 &&
                (mipsolver.colCost(col) < 0 ||
                 mipsolver.mipdata_->downlocks[col] != 0))
              roundsol[col] = std::ceil(fracint.second);
            else
              roundsol[col] = std::floor(fracint.second);
          }

          const auto& cliquesubst =
              mipsolver.mipdata_->cliquetable.getSubstitutions();
          for (int k = cliquesubst.size() - 1; k >= 0; --k) {
            if (cliquesubst[k].replace.val == 0)
              roundsol[cliquesubst[k].substcol] =
                  1 - roundsol[cliquesubst[k].replace.col];
            else
              roundsol[cliquesubst[k].substcol] =
                  roundsol[cliquesubst[k].replace.col];
          }

          for (int i = 0; i != mipsolver.numCol(); ++i)
            objsum += roundsol[i] * mipsolver.colCost(i);

          mipsolver.mipdata_->addIncumbent(roundsol, double(objsum), 'S');
          objsum = 0;
        }

        for (int i = 0; i != mipsolver.numCol(); ++i)
          objsum += sol.col_value[i] * mipsolver.colCost(i);

        objective = double(objsum);
        break;
      }
      case Status::Infeasible:
        objective = HIGHS_CONST_INF;
        break;
      default:
        break;
    }
  } while (solveagain);

  return status;
}