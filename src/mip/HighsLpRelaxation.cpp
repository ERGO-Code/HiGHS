#include "mip/HighsLpRelaxation.h"

#include "mip/HighsCutPool.h"
#include "mip/HighsDomain.h"
#include "mip/HighsMipSolver.h"
#include "mip/HighsMipSolverData.h"
#include "mip/HighsPseudocost.h"
#include "util/HighsCDouble.h"

HighsLpRelaxation::HighsLpRelaxation(const HighsMipSolver& mipsolver)
    : mipsolver(mipsolver) {
  lpsolver.passModel(*mipsolver.model_);
  lpsolver.setHighsOptionValue(
      "primal_feasibility_tolerance",
      mipsolver.options_mip_->mip_feasibility_tolerance);
  lpsolver.setHighsOptionValue(
      "dual_feasibility_tolerance",
      mipsolver.options_mip_->mip_feasibility_tolerance * 0.1);
  status = Status::NotSet;
  numlpiters = 0;
  currentbasisstored = false;
}

HighsLpRelaxation::HighsLpRelaxation(const HighsLpRelaxation& other)
    : mipsolver(other.mipsolver),
      lp2cutpoolindex(other.lp2cutpoolindex),
      fractionalints(other.fractionalints),
      objective(other.objective),
      basischeckpoint(other.basischeckpoint),
      currentbasisstored(other.currentbasisstored) {
  lpsolver.passModel(other.lpsolver.getLp());
  lpsolver.setBasis(other.lpsolver.getBasis());
  lpsolver.passHighsOptions(other.lpsolver.getHighsOptions());
  lpsolver.setHighsLogfile(NULL);
  numlpiters = 0;
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
  if (numcuts > 0) {
    status = Status::NotSet;
    currentbasisstored = false;
    lp2cutpoolindex.insert(lp2cutpoolindex.end(), cutset.cutindices.begin(),
                           cutset.cutindices.end());
    lpsolver.addRows(numcuts, cutset.lower_.data(), cutset.upper_.data(),
                     cutset.ARvalue_.size(), cutset.ARstart_.data(),
                     cutset.ARindex_.data(), cutset.ARvalue_.data());
    cutset.clear();
  }
}

void HighsLpRelaxation::removeCuts(int ndelcuts, std::vector<int>& deletemask) {
  if (ndelcuts > 0) {
    HighsBasis basis = lpsolver.getBasis();
    int nlprows = lpsolver.getNumRows();
    lpsolver.deleteRows(deletemask.data());
    for (int i = mipsolver.numRow(); i != nlprows; ++i) {
      if (deletemask[i] >= 0) {
        lp2cutpoolindex[deletemask[i] - mipsolver.numRow()] =
            lp2cutpoolindex[i - mipsolver.numRow()];
        basis.row_status[deletemask[i]] = basis.row_status[i];
      }
    }

    basis.row_status.resize(basis.row_status.size() - ndelcuts);
    lp2cutpoolindex.resize(lp2cutpoolindex.size() - ndelcuts);
    lpsolver.setBasis(basis);
    lpsolver.run();
  }
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
  const std::vector<double>& row_dual = lpsolver.getSolution().row_dual;

  const HighsLp& lp = lpsolver.getLp();

  assert(std::isfinite(upperbound));
  HighsCDouble upper = upperbound;

  for (int i = 0; i != lp.numRow_; ++i) {
    if (std::abs(row_dual[i]) <= mipsolver.mipdata_->feastol) continue;
    if (row_dual[i] > 0)
      upper += row_dual[i] * lp.rowUpper_[i];
    else
      upper += row_dual[i] * lp.rowLower_[i];
  }

  inds.clear();
  vals.clear();
  for (int i = 0; i != lp.numCol_; ++i) {
    int start = lp.Astart_[i];
    int end = lp.Astart_[i + 1];

    HighsCDouble sum = lp.colCost_[i];

    for (int j = start; j != end; ++j) {
      if (std::abs(row_dual[lp.Aindex_[j]]) <= mipsolver.mipdata_->feastol)
        continue;
      sum += lp.Avalue_[j] * row_dual[lp.Aindex_[j]];
    }

    double val = double(sum);

    if (std::abs(val) <= mipsolver.mipdata_->epsilon) continue;

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

#ifdef HIGHS_DEBUGSOL
  HighsCDouble debugactivity = 0;
  for (size_t i = 0; i != inds.size(); ++i)
    debugactivity += mip.debugSolution_[inds[i]] * vals[i];

  assert(debugactivity <= rhs + mipsolver.mipdata_->feastol);
#endif

  return true;
}

void HighsLpRelaxation::storeDualInfProof() {
  assert(lpsolver.getModelStatus(true) == HighsModelStatus::PRIMAL_INFEASIBLE);

  int numrow = lpsolver.getNumRows();
  bool hasdualray = false;
  lpsolver.getDualRay(hasdualray);

  if (!hasdualray) {
    printf("no dual ray stored\n");
    return;
  }

  dualproofinds.clear();
  dualproofvals.clear();
  dualproofrhs = HIGHS_CONST_INF;
  const HighsLp& lp = lpsolver.getLp();
  dualproofbuffer.resize(numrow);

  lpsolver.getDualRay(hasdualray, dualproofbuffer.data());
  std::vector<double>& dualray = dualproofbuffer;

  HighsCDouble upper = 0.0;
  double scale = 0.0;

  for (int i = 0; i != lp.numRow_; ++i) {
    if (std::abs(dualray[i]) <=
        lpsolver.getHighsOptions().dual_feasibility_tolerance) {
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

    HighsCDouble sum = 0.0;

    for (int j = start; j != end; ++j) {
      if (dualray[lp.Aindex_[j]] == 0.0) continue;
      sum -= lp.Avalue_[j] * dualray[lp.Aindex_[j]];
    }

    double val = scale * double(sum);

    if (std::abs(val) <= mipsolver.mipdata_->epsilon) continue;

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

#ifdef HIGHS_DEBUGSOL
  HighsCDouble debugactivity = 0;
  for (size_t i = 0; i != dualproofinds.size(); ++i)
    debugactivity += mip.debugSolution_[inds[i]] * dualproofvals[i];

  assert(debugactivity <= rhs + mipsolver.mipdata_->feastol);
#endif
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

#ifdef HIGHS_DEBUGSOL
  HighsCDouble debugactivity = 0;
  for (size_t i = 0; i != dualproofinds.size(); ++i)
    debugactivity += mip.debugSolution_[inds[i]] * dualproofvals[i];

  assert(debugactivity <= rhs + mipsolver.mipdata_->feastol);
#endif
}

bool HighsLpRelaxation::checkDualProof() const {
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
  assert(checkDualProof());

  inds = dualproofinds;
  vals = dualproofvals;
  rhs = dualproofrhs;
  return true;
}

void HighsLpRelaxation::recoverBasis() {
  if (basischeckpoint && basischeckpoint->valid_) {
    assert(isBasisRightSize(lpsolver.getLp(), *basischeckpoint));
    assert(isBasisConsistent(lpsolver.getLp(), *basischeckpoint));
    // basischeckpoint->row_status.resize(lpsolver.getNumRows(),
    //                                  HighsBasisStatus::BASIC);

    lpsolver.clearSolver();
    lpsolver.setBasis(*basischeckpoint);
  }
}

HighsLpRelaxation::Status HighsLpRelaxation::run(bool resolve_on_error) {
  HighsStatus callstatus;
  try {
    callstatus = lpsolver.run();
  } catch (const std::runtime_error&) {
    callstatus = HighsStatus::Error;
  }

  const HighsInfo& info = lpsolver.getHighsInfo();
  assert(info.max_primal_infeasibility >= 0);
  assert(info.max_dual_infeasibility >= 0);
  numlpiters += info.simplex_iteration_count;

  if (callstatus == HighsStatus::Error) {
    lpsolver.clearSolver();
    recoverBasis();
    if (resolve_on_error) return run(false);

    return Status::Error;
  }

  HighsModelStatus scaledmodelstatus = lpsolver.getModelStatus(true);

  switch (scaledmodelstatus) {
    case HighsModelStatus::REACHED_DUAL_OBJECTIVE_VALUE_UPPER_BOUND:
      storeDualUBProof();
      if (checkDualProof()) return Status::Infeasible;

      return Status::Error;
    case HighsModelStatus::PRIMAL_INFEASIBLE:
      storeDualInfProof();
      if (checkDualProof()) return Status::Infeasible;

      if (resolve_on_error) {
        int scalestrategy = lpsolver.getHighsOptions().simplex_scale_strategy;
        lpsolver.setHighsOptionValue("simplex_scale_strategy", 0);
        HighsBasis basis = lpsolver.getBasis();
        lpsolver.clearSolver();
        lpsolver.setBasis(basis);
        auto tmp = run(false);
        lpsolver.setHighsOptionValue("simplex_scale_strategy", scalestrategy);
        return tmp;
      }

      assert(false);
      return Status::Error;
    case HighsModelStatus::OPTIMAL:
      if (info.max_primal_infeasibility <= mipsolver.mipdata_->feastol &&
          info.max_dual_infeasibility <= mipsolver.mipdata_->feastol)
        return Status::Optimal;

      if (resolve_on_error) {
        int scalestrategy = lpsolver.getHighsOptions().simplex_scale_strategy;
        lpsolver.setHighsOptionValue("simplex_scale_strategy", 0);
        HighsBasis basis = lpsolver.getBasis();
        lpsolver.clearSolver();
        lpsolver.setBasis(basis);
        auto tmp = run(false);
        lpsolver.setHighsOptionValue("simplex_scale_strategy", scalestrategy);
        return tmp;
      }

      if (info.max_primal_infeasibility <= mipsolver.mipdata_->feastol) {
        return Status::UnscaledPrimalFeasible;
      }

      if (info.max_dual_infeasibility <= mipsolver.mipdata_->feastol) {
        return Status::UnscaledDualFeasible;
      }

      return Status::UnscaledInfeasible;
    case HighsModelStatus::PRIMAL_DUAL_INFEASIBLE:
    // case HighsModelStatus::PRIMAL_INFEASIBLE:
    //  if (lpsolver.getModelStatus(false) == scaledmodelstatus)
    //    return Status::Infeasible;
    //  return Status::Error;
    default:
      return Status::Error;
  }
}

HighsLpRelaxation::Status HighsLpRelaxation::resolveLp() {
  fractionalints.clear();

  status = run();

  switch (status) {
    case Status::UnscaledInfeasible:
    case Status::UnscaledDualFeasible:
    case Status::UnscaledPrimalFeasible:
    case Status::Optimal: {
      const HighsSolution& sol = lpsolver.getSolution();

      HighsCDouble objsum = 0;
      for (int i = 0; i != mipsolver.numCol(); ++i) {
        if (mipsolver.variableType(i) != HighsVarType::INTEGER) continue;

        // for the fractionality we assume that LP bounds are not violated
        // bounds that are violated by the unscaled LP are indicated by the
        // return status already
        double val =
            std::max(std::min(sol.col_value[i], lpsolver.getLp().colUpper_[i]),
                     lpsolver.getLp().colLower_[i]);
        double intval = std::floor(val + 0.5);

        if (std::abs(val - intval) > mipsolver.mipdata_->feastol)
          fractionalints.emplace_back(i, val);
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

  return status;
}