#include "mip/HighsDomain.h"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <queue>

#include "mip/HighsCutPool.h"
#include "mip/HighsMipSolverData.h"

static double activityContributionMin(double coef, const double& lb,
                                      const double& ub) {
  if (coef < 0) {
    if (ub == HIGHS_CONST_INF) return -HIGHS_CONST_INF;

    return coef * ub;
  } else {
    if (lb == -HIGHS_CONST_INF) return -HIGHS_CONST_INF;

    return coef * lb;
  }
}

static double activityContributionMax(double coef, const double& lb,
                                      const double& ub) {
  if (coef < 0) {
    if (lb == -HIGHS_CONST_INF) return HIGHS_CONST_INF;

    return coef * lb;
  } else {
    if (ub == HIGHS_CONST_INF) return HIGHS_CONST_INF;

    return coef * ub;
  }
}

HighsDomain::HighsDomain(HighsMipSolver& mipsolver) : mipsolver(&mipsolver) {
  colLower_ = mipsolver.model_->colLower_;
  colUpper_ = mipsolver.model_->colUpper_;
  changedcolsflags_.resize(mipsolver.numCol());
  changedcols_.reserve(mipsolver.numCol());
  infeasible_reason = Reason::unspecified();
  infeasible_ = false;
}

void HighsDomain::addCutpool(HighsCutPool& cutpool) {
  int cutpoolindex = cutpoolpropagation.size();
  cutpoolpropagation.emplace_back(cutpoolindex, this, cutpool);
}

HighsDomain::CutpoolPropagation::CutpoolPropagation(int cutpoolindex,
                                                    HighsDomain* domain,
                                                    HighsCutPool& cutpool)
    : cutpoolindex(cutpoolindex), domain(domain), cutpool(&cutpool) {
  cutpool.addPropagationDomain(this);
}

HighsDomain::CutpoolPropagation::CutpoolPropagation(
    const CutpoolPropagation& other)
    : cutpoolindex(other.cutpoolindex),
      domain(other.domain),
      cutpool(other.cutpool),
      activitycuts_(other.activitycuts_),
      activitycutsinf_(other.activitycutsinf_),
      activitycutversion_(other.activitycutversion_),
      propagatecutflags_(other.propagatecutflags_),
      propagatecutinds_(other.propagatecutinds_) {
  cutpool->addPropagationDomain(this);
}

HighsDomain::CutpoolPropagation::~CutpoolPropagation() {
  cutpool->removePropagationDomain(this);
}

void HighsDomain::CutpoolPropagation::cutAdded(int cut) {
  int start = cutpool->getMatrix().getRowStart(cut);
  int end = cutpool->getMatrix().getRowEnd(cut);
  const int* arindex = cutpool->getMatrix().getARindex();
  const double* arvalue = cutpool->getMatrix().getARvalue();

  if (int(activitycuts_.size()) <= cut) {
    activitycuts_.resize(cut + 1);
    activitycutsinf_.resize(cut + 1);
    propagatecutflags_.resize(cut + 1);
    activitycutversion_.resize(cut + 1);
  }

  activitycutversion_[cut] = cutpool->getModificationCount(cut);
  domain->computeMinActivity(start, end, arindex, arvalue,
                             activitycutsinf_[cut], activitycuts_[cut]);

  if (activitycutsinf_[cut] <= 1 && !propagatecutflags_[cut]) {
    markPropagateCut(cut);
    // propagatecutflags_[cut] = 1;
    // propagatecutinds_.push_back(cut);
  }
}

void HighsDomain::CutpoolPropagation::markPropagateCut(int cut) {
  if (!propagatecutflags_[cut] &&
      (activitycutsinf_[cut] == 1 ||
       (cutpool->getRhs()[cut] - activitycuts_[cut]) /
               cutpool->getMaxAbsCutCoef(cut) <
           1.0 - domain->mipsolver->mipdata_->feastol)) {
    propagatecutinds_.push_back(cut);
    propagatecutflags_[cut] = 1;
  }
}

void HighsDomain::CutpoolPropagation::updateActivityLbChange(int col,
                                                             double oldbound,
                                                             double newbound) {
  cutpool->getMatrix().forEachColumnEntry(col, [&](int row, double val) {
    if (val > 0) {
      double deltamin;

      assert(row < int(activitycutversion_.size()));

      if (activitycutversion_[row] != cutpool->getModificationCount(row)) {
        int start = cutpool->getMatrix().getRowStart(row);
        int end = cutpool->getMatrix().getRowEnd(row);
        const int* arindex = cutpool->getMatrix().getARindex();
        const double* arvalue = cutpool->getMatrix().getARvalue();

        domain->computeMinActivity(start, end, arindex, arvalue,
                                   activitycutsinf_[row], activitycuts_[row]);

        deltamin = HIGHS_CONST_INF;
      } else {
        if (oldbound == -HIGHS_CONST_INF) {
          --activitycutsinf_[row];
          deltamin = newbound * val;
        } else if (newbound == -HIGHS_CONST_INF) {
          ++activitycutsinf_[row];
          deltamin = -oldbound * val;
        } else {
          deltamin = (newbound - oldbound) * val;
        }
        activitycuts_[row] += deltamin;
      }

      if (activitycutsinf_[row] == 0 &&
          activitycuts_[row] - cutpool->getRhs()[row] >
              domain->mipsolver->mipdata_->feastol) {
        // todo, now that multiple cutpools are possible the index needs to be
        // encoded differently
        domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
        domain->infeasible_ = true;
        domain->infeasible_reason = Reason::cut(cutpoolindex, row);
      }

      if (deltamin > 0 && activitycutsinf_[row] <= 1 &&
          !propagatecutflags_[row]) {
        markPropagateCut(row);
        // propagatecutflags_[row] = 1;
        // propagatecutinds_.push_back(row);
      }
    }

    return true;
  });
}

void HighsDomain::CutpoolPropagation::updateActivityUbChange(int col,
                                                             double oldbound,
                                                             double newbound) {
  cutpool->getMatrix().forEachColumnEntry(col, [&](int row, double val) {
    if (val < 0) {
      double deltamin;

      assert(row < int(activitycutversion_.size()));

      if (activitycutversion_[row] != cutpool->getModificationCount(row)) {
        int start = cutpool->getMatrix().getRowStart(row);
        int end = cutpool->getMatrix().getRowEnd(row);
        const int* arindex = cutpool->getMatrix().getARindex();
        const double* arvalue = cutpool->getMatrix().getARvalue();

        domain->computeMinActivity(start, end, arindex, arvalue,
                                   activitycutsinf_[row], activitycuts_[row]);

        activitycutversion_[row] = cutpool->getModificationCount(row);

        deltamin = HIGHS_CONST_INF;
      } else {
        if (oldbound == HIGHS_CONST_INF) {
          --activitycutsinf_[row];
          deltamin = newbound * val;
        } else if (newbound == HIGHS_CONST_INF) {
          ++activitycutsinf_[row];
          deltamin = -oldbound * val;
        } else {
          deltamin = (newbound - oldbound) * val;
        }
        activitycuts_[row] += deltamin;
      }

      if (activitycutsinf_[row] == 0 &&
          activitycuts_[row] - cutpool->getRhs()[row] >
              domain->mipsolver->mipdata_->feastol) {
        domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
        domain->infeasible_ = true;
        domain->infeasible_reason = Reason::cut(cutpoolindex, row);
      }

      if (deltamin > 0 && activitycutsinf_[row] <= 1 &&
          !propagatecutflags_[row]) {
        markPropagateCut(row);
        // propagatecutflags_[row] = 1;
        // propagatecutinds_.push_back(row);
      }
    }

    return true;
  });
}

void HighsDomain::computeMinActivity(int start, int end, const int* ARindex,
                                     const double* ARvalue, int& ninfmin,
                                     HighsCDouble& activitymin) {
  activitymin = 0.0;
  ninfmin = 0;
  for (int j = start; j != end; ++j) {
    int col = ARindex[j];
    double val = ARvalue[j];

    assert(col < int(colLower_.size()));

    double contributionmin =
        activityContributionMin(val, colLower_[col], colUpper_[col]);

    if (contributionmin == -HIGHS_CONST_INF)
      ++ninfmin;
    else
      activitymin += contributionmin;
  }

  activitymin.renormalize();
}

void HighsDomain::computeMaxActivity(int start, int end, const int* ARindex,
                                     const double* ARvalue, int& ninfmax,
                                     HighsCDouble& activitymax) {
  activitymax = 0.0;
  ninfmax = 0;
  for (int j = start; j != end; ++j) {
    int col = ARindex[j];
    double val = ARvalue[j];

    assert(col < int(colLower_.size()));

    double contributionmin =
        activityContributionMax(val, colLower_[col], colUpper_[col]);

    if (contributionmin == HIGHS_CONST_INF)
      ++ninfmax;
    else
      activitymax += contributionmin;
  }

  activitymax.renormalize();
}

int HighsDomain::propagateRowUpper(const int* Rindex, const double* Rvalue,
                                   int Rlen, double Rupper,
                                   const HighsCDouble& minactivity, int ninfmin,
                                   HighsDomainChange* boundchgs) {
  assert(std::isfinite(double(minactivity)));
  if (ninfmin > 1) return 0;
  int numchgs = 0;
  for (int i = 0; i != Rlen; ++i) {
    HighsCDouble minresact;
    double actcontribution = activityContributionMin(
        Rvalue[i], colLower_[Rindex[i]], colUpper_[Rindex[i]]);
    if (ninfmin == 1) {
      if (actcontribution != -HIGHS_CONST_INF) continue;

      minresact = minactivity;
    } else {
      minresact = minactivity - actcontribution;
    }

    double bound = double((Rupper - minresact) / Rvalue[i]);

    if (Rvalue[i] > 0) {
      bool accept;

      if (mipsolver->variableType(Rindex[i]) != HighsVarType::CONTINUOUS) {
        bound = std::floor(bound + mipsolver->mipdata_->feastol);
        if (bound < colUpper_[Rindex[i]] &&
            colUpper_[Rindex[i]] - bound >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (colUpper_[Rindex[i]] == HIGHS_CONST_INF)
          accept = true;
        else if (bound + 1000.0 * mipsolver->mipdata_->feastol <
                     colUpper_[Rindex[i]] &&
                 (colUpper_[Rindex[i]] - bound) /
                         std::max(std::abs(colUpper_[Rindex[i]]),
                                  std::abs(bound)) >
                     0.3)
          accept = true;
        else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Upper, Rindex[i], bound};

    } else {
      bool accept;

      if (mipsolver->variableType(Rindex[i]) != HighsVarType::CONTINUOUS) {
        bound = std::ceil(bound - mipsolver->mipdata_->feastol);
        if (bound > colLower_[Rindex[i]] &&
            bound - colLower_[Rindex[i]] >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (colLower_[Rindex[i]] == -HIGHS_CONST_INF)
          accept = true;
        else if (bound - 1000.0 * mipsolver->mipdata_->feastol >
                     colLower_[Rindex[i]] &&
                 (bound - colLower_[Rindex[i]]) /
                         std::max(std::abs(colUpper_[Rindex[i]]),
                                  std::abs(bound)) >
                     0.3)
          accept = true;
        else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Lower, Rindex[i], bound};
    }
  }

  return numchgs;
}

int HighsDomain::propagateRowLower(const int* Rindex, const double* Rvalue,
                                   int Rlen, double Rlower,
                                   const HighsCDouble& maxactivity, int ninfmax,
                                   HighsDomainChange* boundchgs) {
  assert(std::isfinite(double(maxactivity)));
  if (ninfmax > 1) return 0;
  int numchgs = 0;
  for (int i = 0; i != Rlen; ++i) {
    HighsCDouble maxresact;
    double actcontribution = activityContributionMax(
        Rvalue[i], colLower_[Rindex[i]], colUpper_[Rindex[i]]);
    if (ninfmax == 1) {
      if (actcontribution != HIGHS_CONST_INF) continue;

      maxresact = maxactivity;
    } else {
      maxresact = maxactivity - actcontribution;
    }

    double bound = double((Rlower - maxresact) / Rvalue[i]);

    if (Rvalue[i] < 0) {
      bool accept;

      if (mipsolver->variableType(Rindex[i]) != HighsVarType::CONTINUOUS) {
        bound = std::floor(bound + mipsolver->mipdata_->feastol);
        if (bound < colUpper_[Rindex[i]] &&
            colUpper_[Rindex[i]] - bound >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (colUpper_[Rindex[i]] == HIGHS_CONST_INF)
          accept = true;
        else if (bound + 1000.0 * mipsolver->mipdata_->feastol <
                     colUpper_[Rindex[i]] &&
                 (colUpper_[Rindex[i]] - bound) /
                         std::max(std::abs(colUpper_[Rindex[i]]),
                                  std::abs(bound)) >
                     0.3)
          accept = true;
        else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Upper, Rindex[i], bound};
    } else {
      bool accept;

      if (mipsolver->variableType(Rindex[i]) != HighsVarType::CONTINUOUS) {
        bound = std::ceil(bound - mipsolver->mipdata_->feastol);
        if (bound > colLower_[Rindex[i]] &&
            bound - colLower_[Rindex[i]] >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (colLower_[Rindex[i]] == -HIGHS_CONST_INF)
          accept = true;
        else if (bound - 1000.0 * mipsolver->mipdata_->feastol >
                     colLower_[Rindex[i]] &&
                 (bound - colLower_[Rindex[i]]) /
                         std::max(std::abs(colUpper_[Rindex[i]]),
                                  std::abs(bound)) >
                     0.3)
          accept = true;
        else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Lower, Rindex[i], bound};
    }
  }

  return numchgs;
}

void HighsDomain::updateActivityLbChange(int col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  int start = mip->Astart_[col];
  int end = mip->Astart_[col + 1];

  for (int i = start; i != end; ++i) {
    if (mip->Avalue_[i] > 0) {
      double deltamin;
      if (oldbound == -HIGHS_CONST_INF) {
        --activitymininf_[mip->Aindex_[i]];
        deltamin = newbound * mip->Avalue_[i];
      } else if (newbound == -HIGHS_CONST_INF) {
        ++activitymininf_[mip->Aindex_[i]];
        deltamin = -oldbound * mip->Avalue_[i];
      } else {
        deltamin = (newbound - oldbound) * mip->Avalue_[i];
      }
      activitymin_[mip->Aindex_[i]] += deltamin;

      if (mip->rowUpper_[mip->Aindex_[i]] != HIGHS_CONST_INF &&
          activitymininf_[mip->Aindex_[i]] == 0 &&
          activitymin_[mip->Aindex_[i]] - mip->rowUpper_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_reason = Reason::modelRow(mip->Aindex_[i]);
      }

      if (deltamin > 0 && activitymininf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowUpper_[mip->Aindex_[i]] != HIGHS_CONST_INF) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    } else {
      double deltamax;
      if (oldbound == -HIGHS_CONST_INF) {
        --activitymaxinf_[mip->Aindex_[i]];
        deltamax = newbound * mip->Avalue_[i];
      } else if (newbound == -HIGHS_CONST_INF) {
        ++activitymaxinf_[mip->Aindex_[i]];
        deltamax = -oldbound * mip->Avalue_[i];
      } else {
        deltamax = (newbound - oldbound) * mip->Avalue_[i];
      }
      activitymax_[mip->Aindex_[i]] += deltamax;

      if (mip->rowLower_[mip->Aindex_[i]] != -HIGHS_CONST_INF &&
          activitymaxinf_[mip->Aindex_[i]] == 0 &&
          mip->rowLower_[mip->Aindex_[i]] - activitymax_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_reason = Reason::modelRow(mip->Aindex_[i]);
      }

      if (deltamax < 0 && activitymaxinf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowLower_[mip->Aindex_[i]] != -HIGHS_CONST_INF) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    }
  }

  for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
    cutpoolprop.updateActivityLbChange(col, oldbound, newbound);
}

void HighsDomain::updateActivityUbChange(int col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  int start = mip->Astart_[col];
  int end = mip->Astart_[col + 1];

  for (int i = start; i != end; ++i) {
    if (mip->Avalue_[i] > 0) {
      double deltamax;
      if (oldbound == HIGHS_CONST_INF) {
        --activitymaxinf_[mip->Aindex_[i]];
        deltamax = newbound * mip->Avalue_[i];
      } else if (newbound == HIGHS_CONST_INF) {
        ++activitymaxinf_[mip->Aindex_[i]];
        deltamax = -oldbound * mip->Avalue_[i];
      } else {
        deltamax = (newbound - oldbound) * mip->Avalue_[i];
      }
      activitymax_[mip->Aindex_[i]] += deltamax;

      if (mip->rowLower_[mip->Aindex_[i]] != -HIGHS_CONST_INF &&
          activitymaxinf_[mip->Aindex_[i]] == 0 &&
          mip->rowLower_[mip->Aindex_[i]] - activitymax_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = mip->Aindex_[i] + 1;
      }

      if (deltamax < 0 && activitymaxinf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowLower_[mip->Aindex_[i]] != -HIGHS_CONST_INF) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    } else {
      double deltamin;
      if (oldbound == HIGHS_CONST_INF) {
        --activitymininf_[mip->Aindex_[i]];
        deltamin = newbound * mip->Avalue_[i];
      } else if (newbound == HIGHS_CONST_INF) {
        ++activitymininf_[mip->Aindex_[i]];
        deltamin = -oldbound * mip->Avalue_[i];
      } else {
        deltamin = (newbound - oldbound) * mip->Avalue_[i];
      }

      activitymin_[mip->Aindex_[i]] += deltamin;

      if (mip->rowUpper_[mip->Aindex_[i]] != HIGHS_CONST_INF &&
          activitymininf_[mip->Aindex_[i]] == 0 &&
          activitymin_[mip->Aindex_[i]] - mip->rowUpper_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = mip->Aindex_[i] + 1;
      }

      if (deltamin > 0 && activitymininf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowUpper_[mip->Aindex_[i]] != HIGHS_CONST_INF) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    }
  }

  for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
    cutpoolprop.updateActivityUbChange(col, oldbound, newbound);
}

void HighsDomain::markPropagateCut(Reason reason) {
  switch (reason.type) {
    case Reason::kUnknown:
    case Reason::kBranching:
      break;
    case Reason::kModelRow:
      // markPropagate(reason.index);
      break;
    default:
      assert(reason.type >= 0 && reason.type < int(cutpoolpropagation.size()));
      cutpoolpropagation[reason.type].markPropagateCut(reason.index);
  }
}

void HighsDomain::markPropagate(int row) {
  // todo, check if std::min(maxactivity - lhs, rhs - minactivity) <  amax -
  // feastol and only mark in that case

  if (!propagateflags_[row]) {
    bool proplower = mipsolver->rowLower(row) != -HIGHS_CONST_INF &&
                     (activitymaxinf_[row] == 1 ||
                      (activitymax_[row] - mipsolver->rowLower(row)) /
                              mipsolver->mipdata_->maxAbsRowCoef[row] <
                          1.0 - mipsolver->mipdata_->feastol);
    bool propupper = mipsolver->rowUpper(row) != HIGHS_CONST_INF &&
                     (activitymininf_[row] == 1 ||
                      (mipsolver->rowUpper(row) - activitymin_[row]) /
                              mipsolver->mipdata_->maxAbsRowCoef[row] <
                          1.0 - mipsolver->mipdata_->feastol);

    if (proplower || propupper) {
      propagateinds_.push_back(row);
      propagateflags_[row] = 1;
    }
  }
}

void HighsDomain::computeRowActivities() {
  activitymin_.resize(mipsolver->numRow());
  activitymininf_.resize(mipsolver->numRow());
  activitymax_.resize(mipsolver->numRow());
  activitymaxinf_.resize(mipsolver->numRow());
  propagateflags_.resize(mipsolver->numRow());
  propagateinds_.reserve(mipsolver->numRow());

  for (int i = 0; i != mipsolver->numRow(); ++i) {
    int start = mipsolver->mipdata_->ARstart_[i];
    int end = mipsolver->mipdata_->ARstart_[i + 1];

    computeMinActivity(start, end, mipsolver->mipdata_->ARindex_.data(),
                       mipsolver->mipdata_->ARvalue_.data(), activitymininf_[i],
                       activitymin_[i]);
    computeMaxActivity(start, end, mipsolver->mipdata_->ARindex_.data(),
                       mipsolver->mipdata_->ARvalue_.data(), activitymaxinf_[i],
                       activitymax_[i]);

    if ((activitymininf_[i] <= 1 &&
         mipsolver->rowUpper(i) != HIGHS_CONST_INF) ||
        (activitymaxinf_[i] <= 1 &&
         mipsolver->rowLower(i) != -HIGHS_CONST_INF)) {
      markPropagate(i);
      // propagateflags_[i] = 1;
      // propagateinds_.push_back(i);
    }
  }
}

double HighsDomain::doChangeBound(const HighsDomainChange& boundchg) {
  double oldbound;

  if (boundchg.boundtype == HighsBoundType::Lower) {
    oldbound = colLower_[boundchg.column];
    colLower_[boundchg.column] = boundchg.boundval;
    updateActivityLbChange(boundchg.column, oldbound, boundchg.boundval);

    if (!changedcolsflags_[boundchg.column]) {
      changedcolsflags_[boundchg.column] = 1;
      changedcols_.push_back(boundchg.column);
    }
  } else {
    oldbound = colUpper_[boundchg.column];
    colUpper_[boundchg.column] = boundchg.boundval;
    updateActivityUbChange(boundchg.column, oldbound, boundchg.boundval);

    if (!changedcolsflags_[boundchg.column]) {
      changedcolsflags_[boundchg.column] = 1;
      changedcols_.push_back(boundchg.column);
    }
  }
  assert(oldbound != boundchg.boundval);

  return oldbound;
}

void HighsDomain::changeBound(HighsDomainChange boundchg, Reason reason) {
  assert(boundchg.column >= 0);
  assert(infeasible_ == 0);

  mipsolver->mipdata_->debugSolution.boundChangeAdded(
      *this, boundchg, reason.type == Reason::kBranching);

  if (boundchg.boundtype == HighsBoundType::Lower) {
    if (boundchg.boundval <= colLower_[boundchg.column]) return;
    if (boundchg.boundval > colUpper_[boundchg.column]) {
      if (boundchg.boundval - colUpper_[boundchg.column] >
          mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_reason = reason;
        return;
      }

      boundchg.boundval = colUpper_[boundchg.column];
      if (boundchg.boundval == colLower_[boundchg.column]) return;
    }
  } else {
    if (boundchg.boundval >= colUpper_[boundchg.column]) return;
    if (boundchg.boundval < colLower_[boundchg.column]) {
      if (colLower_[boundchg.column] - boundchg.boundval >
          mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_reason = reason;
        return;
      }

      boundchg.boundval = colLower_[boundchg.column];
      if (boundchg.boundval == colUpper_[boundchg.column]) return;
    }
  }

  double oldbound = doChangeBound(boundchg);

  prevboundval_.push_back(oldbound);
  domchgstack_.push_back(boundchg);
  domchgreason_.push_back(reason);
}

void HighsDomain::setDomainChangeStack(
    const std::vector<HighsDomainChange>& domchgstack) {
  infeasible_ = false;
  mipsolver->mipdata_->debugSolution.resetDomain(*this);

  prevboundval_.clear();
  domchgstack_.clear();
  domchgreason_.clear();
  int stacksize = domchgstack.size();
  for (int k = 0; k != stacksize; ++k) {
    if (domchgstack[k].boundtype == HighsBoundType::Lower &&
        domchgstack[k].boundval <= colLower_[domchgstack[k].column])
      continue;
    if (domchgstack[k].boundtype == HighsBoundType::Upper &&
        domchgstack[k].boundval >= colUpper_[domchgstack[k].column])
      continue;

    mipsolver->mipdata_->debugSolution.boundChangeAdded(*this, domchgstack[k],
                                                        true);

    changeBound(domchgstack[k], Reason::unspecified());

    if (infeasible_) break;
  }
}

HighsDomainChange HighsDomain::backtrack() {
  int k = domchgstack_.size() - 1;

  while (k >= 0) {
    double prevbound = prevboundval_[k];

    mipsolver->mipdata_->debugSolution.boundChangeRemoved(*this,
                                                          domchgstack_[k]);

    // change back to global bound
    doChangeBound(
        {domchgstack_[k].boundtype, domchgstack_[k].column, prevbound});

    if (domchgreason_[k].type == Reason::kBranching) break;

    --k;
  }

  if (infeasible_) {
    markPropagateCut(infeasible_reason);
    infeasible_reason = Reason::unspecified();
    infeasible_ = false;
  }

  int numreason = domchgreason_.size();
  for (int i = k + 1; i < numreason; ++i) markPropagateCut(domchgreason_[i]);

  if (k < 0) {
    domchgstack_.clear();
    prevboundval_.clear();
    domchgreason_.clear();
    return HighsDomainChange({HighsBoundType::Lower, -1, 0.0});
  }

  HighsDomainChange backtrackboundchg = domchgstack_[k];
  domchgstack_.erase(domchgstack_.begin() + k, domchgstack_.end());
  domchgreason_.resize(k);
  prevboundval_.resize(k);

  return backtrackboundchg;
}

void HighsDomain::propagate() {
  std::vector<int> propagateinds;

  auto havePropagationRows = [&]() {
    bool haverows = true;
    if (propagateinds_.empty()) {
      haverows = false;
      for (const auto& cutpoolprop : cutpoolpropagation) {
        if (!cutpoolprop.propagatecutinds_.empty()) {
          haverows = true;
          break;
        }
      }
    }
    return haverows;
  };

  if (!havePropagationRows()) return;

  size_t changedboundsize = 2 * mipsolver->mipdata_->ARvalue_.size();

  for (const auto& cutpoolprop : cutpoolpropagation)
    changedboundsize = std::max(
        changedboundsize, cutpoolprop.cutpool->getMatrix().nonzeroCapacity());

  std::unique_ptr<HighsDomainChange[]> changedbounds(
      new HighsDomainChange[changedboundsize]);

  while (havePropagationRows()) {
    if (!propagateinds_.empty()) {
      propagateinds.swap(propagateinds_);

      int propnnz = 0;
      int numproprows = propagateinds.size();
      for (int i = 0; i != numproprows; ++i) {
        int row = propagateinds[i];
        propagateflags_[row] = 0;
        propnnz += mipsolver->mipdata_->ARstart_[i + 1] -
                   mipsolver->mipdata_->ARstart_[i];
      }

      if (!infeasible_) {
        propRowNumChangedBounds_.assign(numproprows, 0);

        auto propagateIndex = [&](int k) {
          // for (int k = 0; k != numproprows; ++k) {
          int i = propagateinds[k];
          int start = mipsolver->mipdata_->ARstart_[i];
          int end = mipsolver->mipdata_->ARstart_[i + 1];
          int Rlen = end - start;
          const int* Rindex = &mipsolver->mipdata_->ARindex_[start];
          const double* Rvalue = &mipsolver->mipdata_->ARvalue_[start];
          int numchgs = 0;

          if (mipsolver->rowUpper(i) != HIGHS_CONST_INF) {
            // computeMinActivity(start, end, mipsolver->ARstart_.data(),
            // mipsolver->ARvalue_.data(), activitymininf_[i],
            //           activitymin_[i]);
            activitymin_[i].renormalize();
            numchgs = propagateRowUpper(
                Rindex, Rvalue, Rlen, mipsolver->rowUpper(i), activitymin_[i],
                activitymininf_[i], &changedbounds[2 * start]);
          }

          if (mipsolver->rowLower(i) != -HIGHS_CONST_INF) {
            // computeMaxActivity(start, end, mipsolver->ARstart_.data(),
            // mipsolver->ARvalue_.data(), activitymaxinf_[i],
            //           activitymax_[i]);
            activitymax_[i].renormalize();
            numchgs += propagateRowLower(
                Rindex, Rvalue, Rlen, mipsolver->rowLower(i), activitymax_[i],
                activitymaxinf_[i], &changedbounds[2 * start + numchgs]);
          }

          propRowNumChangedBounds_[k] = numchgs;
        };

        // printf("numproprows (model): %d\n", numproprows);

        for (int k = 0; k != numproprows; ++k) propagateIndex(k);

        for (int k = 0; k != numproprows; ++k) {
          if (propRowNumChangedBounds_[k] == 0) continue;
          int i = propagateinds[k];
          int start = 2 * mipsolver->mipdata_->ARstart_[i];
          int end = start + propRowNumChangedBounds_[k];
          for (int j = start; j != end && !infeasible_; ++j)
            changeBound(changedbounds[j], Reason::modelRow(i));

          if (infeasible_) break;
        }
      }

      propagateinds.clear();
    }

    const int numpools = cutpoolpropagation.size();
    for (int cutpool = 0; cutpool != numpools; ++cutpool) {
      auto& cutpoolprop = cutpoolpropagation[cutpool];
      if (!cutpoolprop.propagatecutinds_.empty()) {
        propagateinds.swap(cutpoolprop.propagatecutinds_);

        int propnnz = 0;
        int numproprows = propagateinds.size();

        for (int i = 0; i != numproprows; ++i) {
          int cut = propagateinds[i];
          cutpoolprop.propagatecutflags_[cut] = 0;
          propnnz += cutpoolprop.cutpool->getMatrix().getRowEnd(cut) -
                     cutpoolprop.cutpool->getMatrix().getRowStart(cut);
        }

        if (!infeasible_) {
          propRowNumChangedBounds_.assign(numproprows, 0);

          auto propagateIndex = [&](int k) {
            int i = propagateinds[k];

            int Rlen;
            const int* Rindex;
            const double* Rvalue;
            cutpoolprop.cutpool->getCut(i, Rlen, Rindex, Rvalue);

            if (cutpoolprop.activitycutversion_[i] !=
                cutpoolprop.cutpool->getModificationCount(i)) {
              cutpoolprop.activitycutversion_[i] =
                  cutpoolprop.cutpool->getModificationCount(i);
              int start = cutpoolprop.cutpool->getMatrix().getRowStart(i);
              if (start == -1) {
                cutpoolprop.activitycuts_[i] = 0;
                return;
              }
              int end = cutpoolprop.cutpool->getMatrix().getRowEnd(i);
              const int* arindex =
                  cutpoolprop.cutpool->getMatrix().getARindex();
              const double* arvalue =
                  cutpoolprop.cutpool->getMatrix().getARvalue();

              computeMinActivity(start, end, arindex, arvalue,
                                 cutpoolprop.activitycutsinf_[i],
                                 cutpoolprop.activitycuts_[i]);
            } else
              cutpoolprop.activitycuts_[i].renormalize();

            propRowNumChangedBounds_[k] = propagateRowUpper(
                Rindex, Rvalue, Rlen, cutpoolprop.cutpool->getRhs()[i],
                cutpoolprop.activitycuts_[i], cutpoolprop.activitycutsinf_[i],
                &changedbounds[cutpoolprop.cutpool->getMatrix().getRowStart(
                    i)]);
          };

          // printf("numproprows (cuts): %d\n", numproprows);

          for (int k = 0; k != numproprows; ++k) propagateIndex(k);

          for (int k = 0; k != numproprows; ++k) {
            if (propRowNumChangedBounds_[k] == 0) continue;
            int i = propagateinds[k];
            cutpoolprop.cutpool->resetAge(i);
            int start = cutpoolprop.cutpool->getMatrix().getRowStart(i);
            int end = start + propRowNumChangedBounds_[k];
            for (int j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::cut(cutpool, i));

            if (infeasible_) break;
          }
        }

        propagateinds.clear();
      }
    }
  }
}

void HighsDomain::tightenCoefficients(int* inds, double* vals, int len,
                                      double& rhs) const {
  HighsCDouble maxactivity = 0;

  for (int i = 0; i != len; ++i) {
    if (vals[i] > 0) {
      if (colUpper_[inds[i]] == HIGHS_CONST_INF) return;

      maxactivity += colUpper_[inds[i]] * vals[i];
    } else {
      if (colLower_[inds[i]] == -HIGHS_CONST_INF) return;

      maxactivity += colLower_[inds[i]] * vals[i];
    }
  }

  if (maxactivity - rhs > mipsolver->mipdata_->feastol) {
    HighsCDouble upper = rhs;
    HighsCDouble maxabscoef = double(maxactivity - rhs);
    int tightened = 0;
    for (int i = 0; i != len; ++i) {
      if (vals[i] > maxabscoef) {
        HighsCDouble delta = vals[i] - maxabscoef;
        upper -= delta * colUpper_[inds[i]];
        vals[i] = double(maxabscoef);
        ++tightened;
      } else if (vals[i] < -maxabscoef) {
        HighsCDouble delta = -vals[i] - maxabscoef;
        upper += delta * colLower_[inds[i]];
        vals[i] = -double(maxabscoef);
        ++tightened;
      }
    }

    if (tightened != 0) {
      // printf("tightened %d coefficients, rhs changed from %g to %g\n",
      //       tightened, rhs, double(upper));
      rhs = double(upper);
    }
  }
}

double HighsDomain::getMinCutActivity(const HighsCutPool& cutpool, int cut) {
  for (auto& cutpoolprop : cutpoolpropagation) {
    if (cutpoolprop.cutpool == &cutpool) {
      if (cutpool.getModificationCount(cut) !=
          cutpoolprop.activitycutversion_[cut]) {
        cutpoolprop.activitycutversion_[cut] =
            cutpoolprop.cutpool->getModificationCount(cut);
        int start = cutpoolprop.cutpool->getMatrix().getRowStart(cut);
        if (start == -1) {
          cutpoolprop.activitycuts_[cut] = 0;
          return -HIGHS_CONST_INF;
        }
        int end = cutpoolprop.cutpool->getMatrix().getRowEnd(cut);
        const int* arindex = cutpoolprop.cutpool->getMatrix().getARindex();
        const double* arvalue = cutpoolprop.cutpool->getMatrix().getARvalue();
        computeMinActivity(start, end, arindex, arvalue,
                           cutpoolprop.activitycutsinf_[cut],
                           cutpoolprop.activitycuts_[cut]);
      }
      return cutpoolprop.activitycutsinf_[cut] == 0
                 ? double(cutpoolprop.activitycuts_[cut])
                 : -HIGHS_CONST_INF;
    }
  }

  return -HIGHS_CONST_INF;
}

bool HighsDomain::isFixing(const HighsDomainChange& domchg) const {
  double otherbound = domchg.boundtype == HighsBoundType::Upper
                          ? colLower_[domchg.column]
                          : colUpper_[domchg.column];
  return std::abs(domchg.boundval - otherbound) <= mipsolver->mipdata_->epsilon;
}