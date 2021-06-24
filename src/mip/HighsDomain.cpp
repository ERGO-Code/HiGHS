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
    if (ub == kHighsInf) return -kHighsInf;

    return coef * ub;
  } else {
    if (lb == -kHighsInf) return -kHighsInf;

    return coef * lb;
  }
}

static double activityContributionMax(double coef, const double& lb,
                                      const double& ub) {
  if (coef < 0) {
    if (lb == -kHighsInf) return kHighsInf;

    return coef * lb;
  } else {
    if (ub == kHighsInf) return kHighsInf;

    return coef * ub;
  }
}

HighsDomain::HighsDomain(HighsMipSolver& mipsolver) : mipsolver(&mipsolver) {
  colLower_ = mipsolver.model_->colLower_;
  colUpper_ = mipsolver.model_->colUpper_;
  colLowerPos_.assign(mipsolver.numCol(), -1);
  colUpperPos_.assign(mipsolver.numCol(), -1);
  changedcolsflags_.resize(mipsolver.numCol());
  changedcols_.reserve(mipsolver.numCol());
  infeasible_reason = Reason::unspecified();
  infeasible_ = false;
}

void HighsDomain::addCutpool(HighsCutPool& cutpool) {
  HighsInt cutpoolindex = cutpoolpropagation.size();
  cutpoolpropagation.emplace_back(cutpoolindex, this, cutpool);
}

HighsDomain::CutpoolPropagation::CutpoolPropagation(HighsInt cutpoolindex,
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

void HighsDomain::CutpoolPropagation::cutAdded(HighsInt cut) {
  HighsInt start = cutpool->getMatrix().getRowStart(cut);
  HighsInt end = cutpool->getMatrix().getRowEnd(cut);
  const HighsInt* arindex = cutpool->getMatrix().getARindex();
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

void HighsDomain::CutpoolPropagation::markPropagateCut(HighsInt cut) {
  if (!propagatecutflags_[cut] &&
      (activitycutsinf_[cut] == 1 ||
       (cutpool->getRhs()[cut] - double(activitycuts_[cut]) <
        (1.0 - domain->mipsolver->mipdata_->feastol) *
            cutpool->getMaxAbsCutCoef(cut)))) {
    propagatecutinds_.push_back(cut);
    propagatecutflags_[cut] = 1;
  }
}

void HighsDomain::CutpoolPropagation::updateActivityLbChange(HighsInt col,
                                                             double oldbound,
                                                             double newbound) {
  assert(!domain->infeasible_);
  cutpool->getMatrix().forEachColumnEntry(col, [&](HighsInt row, double val) {
    if (val > 0) {
      double deltamin;

      assert(row < int(activitycutversion_.size()));

      if (activitycutversion_[row] != cutpool->getModificationCount(row)) {
        HighsInt start = cutpool->getMatrix().getRowStart(row);
        HighsInt end = cutpool->getMatrix().getRowEnd(row);
        const HighsInt* arindex = cutpool->getMatrix().getARindex();
        const double* arvalue = cutpool->getMatrix().getARvalue();

        domain->computeMinActivity(start, end, arindex, arvalue,
                                   activitycutsinf_[row], activitycuts_[row]);

        activitycutversion_[row] = cutpool->getModificationCount(row);

        deltamin = kHighsInf;
      } else {
        if (oldbound == -kHighsInf) {
          --activitycutsinf_[row];
          deltamin = newbound * val;
        } else if (newbound == -kHighsInf) {
          ++activitycutsinf_[row];
          deltamin = -oldbound * val;
        } else {
          deltamin = (newbound - oldbound) * val;
        }
        activitycuts_[row] += deltamin;
      }

      if (deltamin <= 0) return true;

      if (activitycutsinf_[row] == 0 &&
          activitycuts_[row] - cutpool->getRhs()[row] >
              domain->mipsolver->mipdata_->feastol) {
        // todo, now that multiple cutpools are possible the index needs to be
        // encoded differently
        domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
        domain->infeasible_ = true;
        domain->infeasible_pos = domain->domchgstack_.size();
        domain->infeasible_reason = Reason::cut(cutpoolindex, row);
        return false;
      }

      if (activitycutsinf_[row] <= 1 && !propagatecutflags_[row]) {
        markPropagateCut(row);
        // propagatecutflags_[row] = 1;
        // propagatecutinds_.push_back(row);
      }
    }

    return true;
  });

  if (domain->infeasible_) {
    assert(domain->infeasible_reason.type == cutpoolindex);
    assert(domain->infeasible_reason.index >= 0);
    std::swap(oldbound, newbound);
    cutpool->getMatrix().forEachColumnEntry(col, [&](HighsInt row, double val) {
      if (val > 0) {
        double deltamin;

        assert(row < int(activitycutversion_.size()));

        if (oldbound == -kHighsInf) {
          --activitycutsinf_[row];
          deltamin = newbound * val;
        } else if (newbound == -kHighsInf) {
          ++activitycutsinf_[row];
          deltamin = -oldbound * val;
        } else {
          deltamin = (newbound - oldbound) * val;
        }
        activitycuts_[row] += deltamin;

        if (domain->infeasible_reason.index == row) return false;
      }
      return true;
    });
  }
}

void HighsDomain::CutpoolPropagation::updateActivityUbChange(HighsInt col,
                                                             double oldbound,
                                                             double newbound) {
  assert(!domain->infeasible_);
  cutpool->getMatrix().forEachColumnEntry(col, [&](HighsInt row, double val) {
    if (val < 0) {
      double deltamin;

      assert(row < int(activitycutversion_.size()));

      if (activitycutversion_[row] != cutpool->getModificationCount(row)) {
        HighsInt start = cutpool->getMatrix().getRowStart(row);
        HighsInt end = cutpool->getMatrix().getRowEnd(row);
        const HighsInt* arindex = cutpool->getMatrix().getARindex();
        const double* arvalue = cutpool->getMatrix().getARvalue();

        domain->computeMinActivity(start, end, arindex, arvalue,
                                   activitycutsinf_[row], activitycuts_[row]);

        activitycutversion_[row] = cutpool->getModificationCount(row);

        deltamin = kHighsInf;
      } else {
        if (oldbound == kHighsInf) {
          --activitycutsinf_[row];
          deltamin = newbound * val;
        } else if (newbound == kHighsInf) {
          ++activitycutsinf_[row];
          deltamin = -oldbound * val;
        } else {
          deltamin = (newbound - oldbound) * val;
        }
        activitycuts_[row] += deltamin;
      }

      if (deltamin <= 0) return true;

      if (activitycutsinf_[row] == 0 &&
          activitycuts_[row] - cutpool->getRhs()[row] >
              domain->mipsolver->mipdata_->feastol) {
        domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
        domain->infeasible_ = true;
        domain->infeasible_pos = domain->domchgstack_.size();
        domain->infeasible_reason = Reason::cut(cutpoolindex, row);
        return false;
      }

      if (activitycutsinf_[row] <= 1 && !propagatecutflags_[row]) {
        markPropagateCut(row);
        // propagatecutflags_[row] = 1;
        // propagatecutinds_.push_back(row);
      }
    }

    return true;
  });

  if (domain->infeasible_) {
    assert(domain->infeasible_reason.type == cutpoolindex);
    assert(domain->infeasible_reason.index >= 0);
    std::swap(oldbound, newbound);
    cutpool->getMatrix().forEachColumnEntry(col, [&](HighsInt row, double val) {
      if (val < 0) {
        double deltamin;

        assert(row < int(activitycutversion_.size()));

        if (oldbound == kHighsInf) {
          --activitycutsinf_[row];
          deltamin = newbound * val;
        } else if (newbound == kHighsInf) {
          ++activitycutsinf_[row];
          deltamin = -oldbound * val;
        } else {
          deltamin = (newbound - oldbound) * val;
        }
        activitycuts_[row] += deltamin;

        if (domain->infeasible_reason.index == row) return false;
      }
      return true;
    });
  }
}

void HighsDomain::computeMinActivity(HighsInt start, HighsInt end,
                                     const HighsInt* ARindex,
                                     const double* ARvalue, HighsInt& ninfmin,
                                     HighsCDouble& activitymin) {
  activitymin = 0.0;
  ninfmin = 0;
  for (HighsInt j = start; j != end; ++j) {
    HighsInt col = ARindex[j];
    double val = ARvalue[j];

    assert(col < int(colLower_.size()));

    double contributionmin =
        activityContributionMin(val, colLower_[col], colUpper_[col]);

    if (contributionmin == -kHighsInf)
      ++ninfmin;
    else
      activitymin += contributionmin;
  }

  activitymin.renormalize();
}

void HighsDomain::computeMaxActivity(HighsInt start, HighsInt end,
                                     const HighsInt* ARindex,
                                     const double* ARvalue, HighsInt& ninfmax,
                                     HighsCDouble& activitymax) {
  activitymax = 0.0;
  ninfmax = 0;
  for (HighsInt j = start; j != end; ++j) {
    HighsInt col = ARindex[j];
    double val = ARvalue[j];

    assert(col < int(colLower_.size()));

    double contributionmin =
        activityContributionMax(val, colLower_[col], colUpper_[col]);

    if (contributionmin == kHighsInf)
      ++ninfmax;
    else
      activitymax += contributionmin;
  }

  activitymax.renormalize();
}

HighsInt HighsDomain::propagateRowUpper(const HighsInt* Rindex,
                                        const double* Rvalue, HighsInt Rlen,
                                        double Rupper,
                                        const HighsCDouble& minactivity,
                                        HighsInt ninfmin,
                                        HighsDomainChange* boundchgs) {
  assert(std::isfinite(double(minactivity)));
  if (ninfmin > 1) return 0;
  HighsInt numchgs = 0;
  for (HighsInt i = 0; i != Rlen; ++i) {
    HighsCDouble minresact;
    double actcontribution = activityContributionMin(
        Rvalue[i], colLower_[Rindex[i]], colUpper_[Rindex[i]]);
    if (ninfmin == 1) {
      if (actcontribution != -kHighsInf) continue;

      minresact = minactivity;
    } else {
      minresact = minactivity - actcontribution;
    }

    HighsCDouble boundVal = (Rupper - minresact) / Rvalue[i];
    if (std::abs(double(boundVal) * kHighsTiny) > mipsolver->mipdata_->feastol)
      continue;

    if (Rvalue[i] > 0) {
      bool accept;

      double bound;
      if (mipsolver->variableType(Rindex[i]) != HighsVarType::kContinuous) {
        bound = std::floor(double(boundVal + mipsolver->mipdata_->feastol));
        if (bound < colUpper_[Rindex[i]] &&
            colUpper_[Rindex[i]] - bound >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(double(boundVal) - colLower_[Rindex[i]]) <=
            mipsolver->mipdata_->epsilon)
          bound = colLower_[Rindex[i]];
        else
          bound = double(boundVal);
        if (colUpper_[Rindex[i]] == kHighsInf)
          accept = true;
        else if (bound + 1000.0 * mipsolver->mipdata_->feastol <
                 colUpper_[Rindex[i]]) {
          double relativeImprove = colUpper_[Rindex[i]] - bound;
          if (colLower_[Rindex[i]] != -kHighsInf)
            relativeImprove /= colUpper_[Rindex[i]] - colLower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(colUpper_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::kUpper, Rindex[i], bound};

    } else {
      bool accept;

      double bound;
      if (mipsolver->variableType(Rindex[i]) != HighsVarType::kContinuous) {
        bound = std::ceil(double(boundVal - mipsolver->mipdata_->feastol));
        if (bound > colLower_[Rindex[i]] &&
            bound - colLower_[Rindex[i]] >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(colUpper_[Rindex[i]] - double(boundVal)) <=
            mipsolver->mipdata_->epsilon)
          bound = colUpper_[Rindex[i]];
        else
          bound = double(boundVal);
        if (colLower_[Rindex[i]] == -kHighsInf)
          accept = true;
        else if (bound - 1000.0 * mipsolver->mipdata_->feastol >
                 colLower_[Rindex[i]]) {
          double relativeImprove = bound - colLower_[Rindex[i]];
          if (colUpper_[Rindex[i]] != kHighsInf)
            relativeImprove /= colUpper_[Rindex[i]] - colLower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(colLower_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::kLower, Rindex[i], bound};
    }
  }

  return numchgs;
}

HighsInt HighsDomain::propagateRowLower(const HighsInt* Rindex,
                                        const double* Rvalue, HighsInt Rlen,
                                        double Rlower,
                                        const HighsCDouble& maxactivity,
                                        HighsInt ninfmax,
                                        HighsDomainChange* boundchgs) {
  assert(std::isfinite(double(maxactivity)));
  if (ninfmax > 1) return 0;
  HighsInt numchgs = 0;
  for (HighsInt i = 0; i != Rlen; ++i) {
    HighsCDouble maxresact;
    double actcontribution = activityContributionMax(
        Rvalue[i], colLower_[Rindex[i]], colUpper_[Rindex[i]]);
    if (ninfmax == 1) {
      if (actcontribution != kHighsInf) continue;

      maxresact = maxactivity;
    } else {
      maxresact = maxactivity - actcontribution;
    }

    HighsCDouble boundVal = (Rlower - maxresact) / Rvalue[i];
    if (std::abs(double(boundVal) * kHighsTiny) > mipsolver->mipdata_->feastol)
      continue;

    if (Rvalue[i] < 0) {
      bool accept;

      double bound;
      if (mipsolver->variableType(Rindex[i]) != HighsVarType::kContinuous) {
        bound = std::floor(double(boundVal + mipsolver->mipdata_->feastol));
        if (bound < colUpper_[Rindex[i]] &&
            colUpper_[Rindex[i]] - bound >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(double(boundVal) - colLower_[Rindex[i]]) <=
            mipsolver->mipdata_->epsilon)
          bound = colLower_[Rindex[i]];
        else
          bound = double(boundVal);
        if (colUpper_[Rindex[i]] == kHighsInf)
          accept = true;
        else if (bound + 1000.0 * mipsolver->mipdata_->feastol <
                 colUpper_[Rindex[i]]) {
          double relativeImprove = colUpper_[Rindex[i]] - bound;
          if (colLower_[Rindex[i]] != -kHighsInf)
            relativeImprove /= colUpper_[Rindex[i]] - colLower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(colUpper_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::kUpper, Rindex[i], bound};
    } else {
      bool accept;

      double bound;
      if (mipsolver->variableType(Rindex[i]) != HighsVarType::kContinuous) {
        bound = std::ceil(double(boundVal - mipsolver->mipdata_->feastol));
        if (bound > colLower_[Rindex[i]] &&
            bound - colLower_[Rindex[i]] >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(colUpper_[Rindex[i]] - double(boundVal)) <=
            mipsolver->mipdata_->epsilon)
          bound = colUpper_[Rindex[i]];
        else
          bound = double(boundVal);
        if (colLower_[Rindex[i]] == -kHighsInf)
          accept = true;
        else if (bound - 1000.0 * mipsolver->mipdata_->feastol >
                 colLower_[Rindex[i]]) {
          double relativeImprove = bound - colLower_[Rindex[i]];
          if (colUpper_[Rindex[i]] != kHighsInf)
            relativeImprove /= colUpper_[Rindex[i]] - colLower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(colLower_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::kLower, Rindex[i], bound};
    }
  }

  return numchgs;
}

void HighsDomain::updateActivityLbChange(HighsInt col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  HighsInt start = mip->Astart_[col];
  HighsInt end = mip->Astart_[col + 1];

  assert(!infeasible_);

  for (HighsInt i = start; i != end; ++i) {
    if (mip->Avalue_[i] > 0) {
      double deltamin;
      if (oldbound == -kHighsInf) {
        --activitymininf_[mip->Aindex_[i]];
        deltamin = newbound * mip->Avalue_[i];
      } else if (newbound == -kHighsInf) {
        ++activitymininf_[mip->Aindex_[i]];
        deltamin = -oldbound * mip->Avalue_[i];
      } else {
        deltamin = (newbound - oldbound) * mip->Avalue_[i];
      }
      activitymin_[mip->Aindex_[i]] += deltamin;

      if (deltamin <= 0) continue;

      if (mip->rowUpper_[mip->Aindex_[i]] != kHighsInf &&
          activitymininf_[mip->Aindex_[i]] == 0 &&
          activitymin_[mip->Aindex_[i]] - mip->rowUpper_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowUpper(mip->Aindex_[i]);
        end = i + 1;
        break;
      }

      if (activitymininf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowUpper_[mip->Aindex_[i]] != kHighsInf) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    } else {
      double deltamax;
      if (oldbound == -kHighsInf) {
        --activitymaxinf_[mip->Aindex_[i]];
        deltamax = newbound * mip->Avalue_[i];
      } else if (newbound == -kHighsInf) {
        ++activitymaxinf_[mip->Aindex_[i]];
        deltamax = -oldbound * mip->Avalue_[i];
      } else {
        deltamax = (newbound - oldbound) * mip->Avalue_[i];
      }
      activitymax_[mip->Aindex_[i]] += deltamax;

      if (deltamax >= 0) continue;

      if (mip->rowLower_[mip->Aindex_[i]] != -kHighsInf &&
          activitymaxinf_[mip->Aindex_[i]] == 0 &&
          mip->rowLower_[mip->Aindex_[i]] - activitymax_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowLower(mip->Aindex_[i]);
        end = i + 1;
        break;
      }

      if (activitymaxinf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowLower_[mip->Aindex_[i]] != -kHighsInf) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    }
  }

  if (!infeasible_) {
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.updateActivityLbChange(col, oldbound, newbound);
  } else {
    assert(infeasible_reason.type == Reason::kModelRowLower ||
           infeasible_reason.type == Reason::kModelRowUpper);
    assert(infeasible_reason.index == mip->Aindex_[end - 1]);
  }

  if (infeasible_) {
    std::swap(oldbound, newbound);
    for (HighsInt i = start; i != end; ++i) {
      if (mip->Avalue_[i] > 0) {
        double deltamin;
        if (oldbound == -kHighsInf) {
          --activitymininf_[mip->Aindex_[i]];
          deltamin = newbound * mip->Avalue_[i];
        } else if (newbound == -kHighsInf) {
          ++activitymininf_[mip->Aindex_[i]];
          deltamin = -oldbound * mip->Avalue_[i];
        } else {
          deltamin = (newbound - oldbound) * mip->Avalue_[i];
        }
        activitymin_[mip->Aindex_[i]] += deltamin;
      } else {
        double deltamax;
        if (oldbound == -kHighsInf) {
          --activitymaxinf_[mip->Aindex_[i]];
          deltamax = newbound * mip->Avalue_[i];
        } else if (newbound == -kHighsInf) {
          ++activitymaxinf_[mip->Aindex_[i]];
          deltamax = -oldbound * mip->Avalue_[i];
        } else {
          deltamax = (newbound - oldbound) * mip->Avalue_[i];
        }
        activitymax_[mip->Aindex_[i]] += deltamax;
      }
    }

    return;
  }
}

void HighsDomain::updateActivityUbChange(HighsInt col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  HighsInt start = mip->Astart_[col];
  HighsInt end = mip->Astart_[col + 1];

  assert(!infeasible_);

  for (HighsInt i = start; i != end; ++i) {
    if (mip->Avalue_[i] > 0) {
      double deltamax;
      if (oldbound == kHighsInf) {
        --activitymaxinf_[mip->Aindex_[i]];
        deltamax = newbound * mip->Avalue_[i];
      } else if (newbound == kHighsInf) {
        ++activitymaxinf_[mip->Aindex_[i]];
        deltamax = -oldbound * mip->Avalue_[i];
      } else {
        deltamax = (newbound - oldbound) * mip->Avalue_[i];
      }
      activitymax_[mip->Aindex_[i]] += deltamax;

      if (deltamax >= 0) continue;

      if (mip->rowLower_[mip->Aindex_[i]] != -kHighsInf &&
          activitymaxinf_[mip->Aindex_[i]] == 0 &&
          mip->rowLower_[mip->Aindex_[i]] - activitymax_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowLower(mip->Aindex_[i]);
        end = i + 1;
        break;
      }

      if (activitymaxinf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowLower_[mip->Aindex_[i]] != -kHighsInf) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    } else {
      double deltamin;
      if (oldbound == kHighsInf) {
        --activitymininf_[mip->Aindex_[i]];
        deltamin = newbound * mip->Avalue_[i];
      } else if (newbound == kHighsInf) {
        ++activitymininf_[mip->Aindex_[i]];
        deltamin = -oldbound * mip->Avalue_[i];
      } else {
        deltamin = (newbound - oldbound) * mip->Avalue_[i];
      }

      activitymin_[mip->Aindex_[i]] += deltamin;

      if (deltamin <= 0) continue;

      if (mip->rowUpper_[mip->Aindex_[i]] != kHighsInf &&
          activitymininf_[mip->Aindex_[i]] == 0 &&
          activitymin_[mip->Aindex_[i]] - mip->rowUpper_[mip->Aindex_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowUpper(mip->Aindex_[i]);
        end = i + 1;
        break;
      }

      if (activitymininf_[mip->Aindex_[i]] <= 1 &&
          !propagateflags_[mip->Aindex_[i]] &&
          mip->rowUpper_[mip->Aindex_[i]] != kHighsInf) {
        markPropagate(mip->Aindex_[i]);
        // propagateflags_[mip->Aindex_[i]] = 1;
        // propagateinds_.push_back(mip->Aindex_[i]);
      }
    }
  }

  if (!infeasible_) {
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.updateActivityUbChange(col, oldbound, newbound);
  } else {
    assert(infeasible_reason.type == Reason::kModelRowLower ||
           infeasible_reason.type == Reason::kModelRowUpper);
    assert(infeasible_reason.index == mip->Aindex_[end - 1]);
  }

  if (infeasible_) {
    std::swap(oldbound, newbound);
    for (HighsInt i = start; i != end; ++i) {
      if (mip->Avalue_[i] > 0) {
        double deltamax;
        if (oldbound == kHighsInf) {
          --activitymaxinf_[mip->Aindex_[i]];
          deltamax = newbound * mip->Avalue_[i];
        } else if (newbound == kHighsInf) {
          ++activitymaxinf_[mip->Aindex_[i]];
          deltamax = -oldbound * mip->Avalue_[i];
        } else {
          deltamax = (newbound - oldbound) * mip->Avalue_[i];
        }
        activitymax_[mip->Aindex_[i]] += deltamax;
      } else {
        double deltamin;
        if (oldbound == kHighsInf) {
          --activitymininf_[mip->Aindex_[i]];
          deltamin = newbound * mip->Avalue_[i];
        } else if (newbound == kHighsInf) {
          ++activitymininf_[mip->Aindex_[i]];
          deltamin = -oldbound * mip->Avalue_[i];
        } else {
          deltamin = (newbound - oldbound) * mip->Avalue_[i];
        }

        activitymin_[mip->Aindex_[i]] += deltamin;
      }
    }

    return;
  }
}

void HighsDomain::markPropagateCut(Reason reason) {
  switch (reason.type) {
    case Reason::kUnknown:
    case Reason::kCliqueTable:
    case Reason::kBranching:
    case Reason::kModelRowLower:
    case Reason::kModelRowUpper:
      break;
    default:
      assert(reason.type >= 0 && reason.type < int(cutpoolpropagation.size()));
      cutpoolpropagation[reason.type].markPropagateCut(reason.index);
  }
}

void HighsDomain::markPropagate(HighsInt row) {
  // todo, check if std::min(maxactivity - lhs, rhs - minactivity) <  amax -
  // feastol and only mark in that case

  if (!propagateflags_[row]) {
    bool proplower = mipsolver->rowLower(row) != -kHighsInf &&
                     (activitymaxinf_[row] == 1 ||
                      (double(activitymax_[row]) - mipsolver->rowLower(row)) /
                              mipsolver->mipdata_->maxAbsRowCoef[row] <
                          1.0 - mipsolver->mipdata_->feastol);
    bool propupper = mipsolver->rowUpper(row) != kHighsInf &&
                     (activitymininf_[row] == 1 ||
                      (mipsolver->rowUpper(row) - double(activitymin_[row])) /
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

  for (HighsInt i = 0; i != mipsolver->numRow(); ++i) {
    HighsInt start = mipsolver->mipdata_->ARstart_[i];
    HighsInt end = mipsolver->mipdata_->ARstart_[i + 1];

    computeMinActivity(start, end, mipsolver->mipdata_->ARindex_.data(),
                       mipsolver->mipdata_->ARvalue_.data(), activitymininf_[i],
                       activitymin_[i]);
    computeMaxActivity(start, end, mipsolver->mipdata_->ARindex_.data(),
                       mipsolver->mipdata_->ARvalue_.data(), activitymaxinf_[i],
                       activitymax_[i]);

    if ((activitymininf_[i] <= 1 && mipsolver->rowUpper(i) != kHighsInf) ||
        (activitymaxinf_[i] <= 1 && mipsolver->rowLower(i) != -kHighsInf)) {
      markPropagate(i);
      // propagateflags_[i] = 1;
      // propagateinds_.push_back(i);
    }
  }
}

double HighsDomain::doChangeBound(const HighsDomainChange& boundchg) {
  double oldbound;

  if (boundchg.boundtype == HighsBoundType::kLower) {
    oldbound = colLower_[boundchg.column];
    colLower_[boundchg.column] = boundchg.boundval;
    if (!infeasible_)
      updateActivityLbChange(boundchg.column, oldbound, boundchg.boundval);

    if (!changedcolsflags_[boundchg.column]) {
      changedcolsflags_[boundchg.column] = 1;
      changedcols_.push_back(boundchg.column);
    }
  } else {
    oldbound = colUpper_[boundchg.column];
    colUpper_[boundchg.column] = boundchg.boundval;
    if (!infeasible_)
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
  HighsInt prevPos;
  if (boundchg.boundtype == HighsBoundType::kLower) {
    if (boundchg.boundval <= colLower_[boundchg.column]) return;
    if (boundchg.boundval > colUpper_[boundchg.column]) {
      if (boundchg.boundval - colUpper_[boundchg.column] >
          mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        if (!infeasible_) infeasible_pos = domchgstack_.size();
        infeasible_ = true;
        infeasible_reason = reason;
      } else {
        boundchg.boundval = colUpper_[boundchg.column];
        if (boundchg.boundval == colLower_[boundchg.column]) return;
      }
    }

    prevPos = colLowerPos_[boundchg.column];
    colLowerPos_[boundchg.column] = domchgstack_.size();
  } else {
    if (boundchg.boundval >= colUpper_[boundchg.column]) return;
    if (boundchg.boundval < colLower_[boundchg.column]) {
      if (colLower_[boundchg.column] - boundchg.boundval >
          mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        if (!infeasible_) infeasible_pos = domchgstack_.size();
        infeasible_ = true;
        infeasible_reason = reason;
      } else {
        boundchg.boundval = colLower_[boundchg.column];
        if (boundchg.boundval == colUpper_[boundchg.column]) return;
      }
    }

    prevPos = colUpperPos_[boundchg.column];
    colUpperPos_[boundchg.column] = domchgstack_.size();
  }

  if (reason.type == Reason::kBranching)
    branchPos_.push_back(domchgstack_.size());

  assert(prevPos < (HighsInt)domchgstack_.size());

  bool binary = isBinary(boundchg.column);

  double oldbound = doChangeBound(boundchg);

  prevboundval_.emplace_back(oldbound, prevPos);
  domchgstack_.push_back(boundchg);
  domchgreason_.push_back(reason);

  if (binary && !infeasible_)
    mipsolver->mipdata_->cliquetable.addImplications(*this, boundchg.column,
                                                     boundchg.boundval > 0.5);
}

void HighsDomain::setDomainChangeStack(
    const std::vector<HighsDomainChange>& domchgstack) {
  infeasible_ = false;
  mipsolver->mipdata_->debugSolution.resetDomain(*this);

  if (!domchgstack_.empty()) {
    for (const HighsDomainChange& domchg : domchgstack_) {
      if (domchg.boundtype == HighsBoundType::kLower)
        colLowerPos_[domchg.column] = -1;
      else
        colUpperPos_[domchg.column] = -1;
    }
  }

  prevboundval_.clear();
  domchgstack_.clear();
  domchgreason_.clear();
  branchPos_.clear();
  HighsInt stacksize = domchgstack.size();
  for (HighsInt k = 0; k != stacksize; ++k) {
    if (domchgstack[k].boundtype == HighsBoundType::kLower &&
        domchgstack[k].boundval <= colLower_[domchgstack[k].column])
      continue;
    if (domchgstack[k].boundtype == HighsBoundType::kUpper &&
        domchgstack[k].boundval >= colUpper_[domchgstack[k].column])
      continue;

    mipsolver->mipdata_->debugSolution.boundChangeAdded(*this, domchgstack[k],
                                                        true);

    changeBound(domchgstack[k], Reason::unspecified());

    if (infeasible_) break;
  }
}

HighsDomainChange HighsDomain::backtrack() {
  HighsInt k = HighsInt(domchgstack_.size()) - 1;
  bool old_infeasible = infeasible_;
  Reason old_reason = infeasible_reason;

  while (k >= 0) {
    double prevbound = prevboundval_[k].first;
    HighsInt prevpos = prevboundval_[k].second;
    assert(prevpos < k);

    mipsolver->mipdata_->debugSolution.boundChangeRemoved(*this,
                                                          domchgstack_[k]);

    if (domchgstack_[k].boundtype == HighsBoundType::kLower) {
      assert(colLowerPos_[domchgstack_[k].column] == k);
      colLowerPos_[domchgstack_[k].column] = prevpos;
    } else {
      assert(colUpperPos_[domchgstack_[k].column] == k);
      colUpperPos_[domchgstack_[k].column] = prevpos;
    }

    // change back to global bound
    doChangeBound(
        {domchgstack_[k].boundtype, domchgstack_[k].column, prevbound});

    if (infeasible_ && infeasible_pos == k) {
      assert(old_infeasible);
      assert(k == HighsInt(domchgstack_.size()) - 1);
      infeasible_ = false;
      infeasible_reason = Reason::unspecified();
    }

    if (domchgreason_[k].type == Reason::kBranching) {
      branchPos_.pop_back();
      break;
    }

    --k;
  }

  if (old_infeasible) {
    markPropagateCut(old_reason);
    infeasible_reason = Reason::unspecified();
    infeasible_ = false;
  }

  HighsInt numreason = domchgreason_.size();
  for (HighsInt i = k + 1; i < numreason; ++i)
    markPropagateCut(domchgreason_[i]);

  if (k < 0) {
    domchgstack_.clear();
    prevboundval_.clear();
    domchgreason_.clear();
    branchPos_.clear();
    return HighsDomainChange({HighsBoundType::kLower, -1, 0.0});
  }

  HighsDomainChange backtrackboundchg = domchgstack_[k];
  domchgstack_.erase(domchgstack_.begin() + k, domchgstack_.end());
  domchgreason_.resize(k);
  prevboundval_.resize(k);

  return backtrackboundchg;
}

bool HighsDomain::propagate() {
  std::vector<HighsInt> propagateinds;

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

  if (!havePropagationRows()) return false;

  size_t changedboundsize = 2 * mipsolver->mipdata_->ARvalue_.size();

  for (const auto& cutpoolprop : cutpoolpropagation)
    changedboundsize = std::max(
        changedboundsize, cutpoolprop.cutpool->getMatrix().nonzeroCapacity());

  std::unique_ptr<HighsDomainChange[]> changedbounds(
      new HighsDomainChange[changedboundsize]);

  while (havePropagationRows()) {
    if (!propagateinds_.empty()) {
      propagateinds.swap(propagateinds_);

      HighsInt propnnz = 0;
      HighsInt numproprows = propagateinds.size();
      for (HighsInt i = 0; i != numproprows; ++i) {
        HighsInt row = propagateinds[i];
        propagateflags_[row] = 0;
        propnnz += mipsolver->mipdata_->ARstart_[i + 1] -
                   mipsolver->mipdata_->ARstart_[i];
      }

      if (!infeasible_) {
        propRowNumChangedBounds_.assign(
            numproprows, std::make_pair(HighsInt{0}, HighsInt{0}));

        auto propagateIndex = [&](HighsInt k) {
          // for (HighsInt k = 0; k != numproprows; ++k) {
          HighsInt i = propagateinds[k];
          HighsInt start = mipsolver->mipdata_->ARstart_[i];
          HighsInt end = mipsolver->mipdata_->ARstart_[i + 1];
          HighsInt Rlen = end - start;
          const HighsInt* Rindex = mipsolver->mipdata_->ARindex_.data() + start;
          const double* Rvalue = mipsolver->mipdata_->ARvalue_.data() + start;

          if (mipsolver->rowUpper(i) != kHighsInf) {
            // computeMinActivity(start, end, mipsolver->ARstart_.data(),
            // mipsolver->ARvalue_.data(), activitymininf_[i],
            //           activitymin_[i]);
            activitymin_[i].renormalize();
            propRowNumChangedBounds_[k].first = propagateRowUpper(
                Rindex, Rvalue, Rlen, mipsolver->rowUpper(i), activitymin_[i],
                activitymininf_[i], &changedbounds[2 * start]);
          }

          if (mipsolver->rowLower(i) != -kHighsInf) {
            // computeMaxActivity(start, end, mipsolver->ARstart_.data(),
            // mipsolver->ARvalue_.data(), activitymaxinf_[i],
            //           activitymax_[i]);
            activitymax_[i].renormalize();
            propRowNumChangedBounds_[k].second = propagateRowLower(
                Rindex, Rvalue, Rlen, mipsolver->rowLower(i), activitymax_[i],
                activitymaxinf_[i],
                &changedbounds[2 * start + propRowNumChangedBounds_[k].first]);
          }
        };

        // printf("numproprows (model): %" HIGHSINT_FORMAT "\n", numproprows);

        for (HighsInt k = 0; k != numproprows; ++k) propagateIndex(k);

        for (HighsInt k = 0; k != numproprows; ++k) {
          if (propRowNumChangedBounds_[k].first != 0) {
            HighsInt i = propagateinds[k];
            HighsInt start = 2 * mipsolver->mipdata_->ARstart_[i];
            HighsInt end = start + propRowNumChangedBounds_[k].first;
            for (HighsInt j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::modelRowUpper(i));

            if (infeasible_) break;
          }
          if (propRowNumChangedBounds_[k].second != 0) {
            HighsInt i = propagateinds[k];
            HighsInt start = 2 * mipsolver->mipdata_->ARstart_[i] +
                             propRowNumChangedBounds_[k].first;
            HighsInt end = start + propRowNumChangedBounds_[k].second;
            for (HighsInt j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::modelRowLower(i));

            if (infeasible_) break;
          }
        }
      }

      propagateinds.clear();
    }

    const HighsInt numpools = cutpoolpropagation.size();
    for (HighsInt cutpool = 0; cutpool != numpools; ++cutpool) {
      auto& cutpoolprop = cutpoolpropagation[cutpool];
      if (!cutpoolprop.propagatecutinds_.empty()) {
        propagateinds.swap(cutpoolprop.propagatecutinds_);

        HighsInt propnnz = 0;
        HighsInt numproprows = propagateinds.size();

        for (HighsInt i = 0; i != numproprows; ++i) {
          HighsInt cut = propagateinds[i];
          cutpoolprop.propagatecutflags_[cut] = 0;
          propnnz += cutpoolprop.cutpool->getMatrix().getRowEnd(cut) -
                     cutpoolprop.cutpool->getMatrix().getRowStart(cut);
        }

        if (!infeasible_) {
          propRowNumChangedBounds_.assign(
              numproprows, std::make_pair(HighsInt{0}, HighsInt{0}));

          auto propagateIndex = [&](HighsInt k) {
            HighsInt i = propagateinds[k];

            HighsInt Rlen;
            const HighsInt* Rindex;
            const double* Rvalue;
            cutpoolprop.cutpool->getCut(i, Rlen, Rindex, Rvalue);

            if (cutpoolprop.activitycutversion_[i] !=
                cutpoolprop.cutpool->getModificationCount(i)) {
              cutpoolprop.activitycutversion_[i] =
                  cutpoolprop.cutpool->getModificationCount(i);
              HighsInt start = cutpoolprop.cutpool->getMatrix().getRowStart(i);
              if (start == -1) {
                cutpoolprop.activitycuts_[i] = 0;
                return;
              }
              HighsInt end = cutpoolprop.cutpool->getMatrix().getRowEnd(i);
              const HighsInt* arindex =
                  cutpoolprop.cutpool->getMatrix().getARindex();
              const double* arvalue =
                  cutpoolprop.cutpool->getMatrix().getARvalue();

              computeMinActivity(start, end, arindex, arvalue,
                                 cutpoolprop.activitycutsinf_[i],
                                 cutpoolprop.activitycuts_[i]);
            } else
              cutpoolprop.activitycuts_[i].renormalize();

            propRowNumChangedBounds_[k].first = propagateRowUpper(
                Rindex, Rvalue, Rlen, cutpoolprop.cutpool->getRhs()[i],
                cutpoolprop.activitycuts_[i], cutpoolprop.activitycutsinf_[i],
                &changedbounds[cutpoolprop.cutpool->getMatrix().getRowStart(
                    i)]);
          };

          // printf("numproprows (cuts): %" HIGHSINT_FORMAT "\n", numproprows);

          for (HighsInt k = 0; k != numproprows; ++k) propagateIndex(k);

          for (HighsInt k = 0; k != numproprows; ++k) {
            if (propRowNumChangedBounds_[k].first == 0) continue;
            HighsInt i = propagateinds[k];
            cutpoolprop.cutpool->resetAge(i);
            HighsInt start = cutpoolprop.cutpool->getMatrix().getRowStart(i);
            HighsInt end = start + propRowNumChangedBounds_[k].first;
            for (HighsInt j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::cut(cutpool, i));

            if (infeasible_) break;
          }
        }

        propagateinds.clear();
      }
    }
  }

  return true;
}

double HighsDomain::getColLowerPos(HighsInt col, HighsInt stackpos,
                                   HighsInt& pos) const {
  double lb = colLower_[col];
  pos = colLowerPos_[col];
  while (pos > stackpos) {
    lb = prevboundval_[pos].first;
    pos = prevboundval_[pos].second;
  }
  return lb;
}

double HighsDomain::getColUpperPos(HighsInt col, HighsInt stackpos,
                                   HighsInt& pos) const {
  double ub = colUpper_[col];
  pos = colUpperPos_[col];
  while (pos > stackpos) {
    ub = prevboundval_[pos].first;
    pos = prevboundval_[pos].second;
  }
  return ub;
}

bool HighsDomain::explainBoundChange(HighsInt pos,
                                     std::vector<HighsInt>& reasonDomChgs) {
  HighsDomain& globaldom = mipsolver->mipdata_->domain;
  switch (domchgreason_[pos].type) {
    case Reason::kUnknown:
    case Reason::kBranching:
      return false;
    case Reason::kCliqueTable: {
      HighsInt col = domchgreason_[pos].index >> 1;
      HighsInt val = domchgreason_[pos].index & 1;
      reasonDomChgs.clear();
      if (val) {
        assert(colLowerPos_[col] >= 0);
        assert(colLowerPos_[col] < domchgstack_.size());
        reasonDomChgs.push_back(colLowerPos_[col]);
      } else {
        assert(colUpperPos_[col] >= 0);
        assert(colUpperPos_[col] < domchgstack_.size());
        reasonDomChgs.push_back(colUpperPos_[col]);
      }

      break;
    }
    case Reason::kModelRowLower: {
      HighsInt rowIndex = domchgreason_[pos].index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      // get the coefficient value of the column for which we want to explain
      // the bound change
      double domchgVal = 0;
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];
        if (col == domchgstack_[pos].column) {
          domchgVal = vals[i];
          break;
        }
      }

      assert(domchgVal != 0);

      double maxAct = globaldom.getMaxActivity(rowIndex);
      // if the global min activity is infinite we do not filter the bounds
      // and just take all as it is rare and we would need to use the local
      // activity as base which is not available as easily as not all bound
      // changes on the stack are used. Using all local bounds of columns that
      // are present in the cut is however guaranteed to be sufficient even
      // though it might not be as strong.
      if (maxAct == kHighsInf) return false;

      // to explain the bound change we start from the bound constraint,
      // multiply it by the columns coefficient in the constraint. Then the
      // bound constraint is a0 * x0 >= b0 and the constraint
      // a0 * x0 + \sum ai * xi >= b. Let the max activity of \sum ai xi be M,
      // then the constraint yields the bound a0 * x0 >= b - M. If M is
      // sufficiently small the bound constraint is implied. Therefore we
      // decrease M by updating it with the stronger local bounds until
      // M <= b - b0 holds.
      double b0 = domchgstack_[pos].boundval;
      if (mipsolver->variableType(domchgstack_[pos].column) !=
          HighsVarType::kContinuous) {
        // in case of an integral variable the bound was rounded and can be
        // relaxed by 1-feastol
        if (domchgstack_[pos].boundtype == HighsBoundType::kLower)
          b0 -= (1.0 - mipsolver->mipdata_->feastol);
        else
          b0 += (1.0 - mipsolver->mipdata_->feastol);
      } else {
        // for a continuous variable we relax the bound by epsilon to
        // accomodate for tiny rounding errors
        if (domchgstack_[pos].boundtype == HighsBoundType::kLower)
          b0 -= mipsolver->mipdata_->epsilon;
        else
          b0 += mipsolver->mipdata_->epsilon;
      }

      // now multiply the bound constraint with the coefficient value in the
      // constraint to obtain b0
      b0 *= domchgVal;

      // compute the lower bound of M that is necessary
      double Mupper = mipsolver->rowLower(rowIndex) - b0;

      // M is the global residual activity initially
      double M = maxAct;
      if (domchgVal < 0)
        M -= domchgVal * globaldom.colLower_[domchgstack_[pos].column];
      else
        M -= domchgVal * globaldom.colUpper_[domchgstack_[pos].column];

      reasonDomChgs.clear();
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];
        if (col == domchgstack_[pos].column) continue;

        HighsInt boundpos;
        if (vals[i] > 0) {
          double ub = getColUpperPos(col, pos, boundpos);
          if (globaldom.colUpper_[col] <= ub) continue;
          M += vals[i] * (ub - globaldom.colUpper_[col]);
        } else {
          double lb = getColLowerPos(col, pos, boundpos);
          if (globaldom.colLower_[col] >= lb) continue;
          M += vals[i] * (lb - globaldom.colLower_[col]);
        }

        reasonDomChgs.push_back(boundpos);
        if (M < Mupper) break;
      }

      break;
    }
    case Reason::kModelRowUpper: {
      HighsInt rowIndex = domchgreason_[pos].index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      // get the coefficient value of the column for which we want to explain
      // the bound change
      double domchgVal = 0;
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];
        if (col == domchgstack_[pos].column) {
          domchgVal = vals[i];
          break;
        }
      }

      assert(domchgVal != 0);

      double minAct = globaldom.getMinActivity(rowIndex);
      // if the global min activity is infinite we do not filter the bounds
      // and just take all as it is rare and we would need to use the local
      // activity as base which is not available as easily as not all bound
      // changes on the stack are used. Using all local bounds of columns that
      // are present in the cut is however guaranteed to be sufficient even
      // though it might not be as strong.
      if (minAct == -kHighsInf) return false;

      // to explain the bound change we start from the bound constraint,
      // multiply it by the columns coefficient in the constraint. Then the
      // bound constraint is a0 * x0 <= b0 and the constraint
      // a0 * x0 + \sum ai * xi <= b. Let the min activity of \sum ai xi be M,
      // then the constraint yields the bound a0 * x0 <= b - M. If M is
      // sufficiently large the bound constraint is implied. Therefore we
      // increase M by updating it with the stronger local bounds until
      // M >= b - b0 holds.
      double b0 = domchgstack_[pos].boundval;
      if (mipsolver->variableType(domchgstack_[pos].column) !=
          HighsVarType::kContinuous) {
        // in case of an integral variable the bound was rounded and can be
        // relaxed by 1-feastol
        if (domchgstack_[pos].boundtype == HighsBoundType::kLower)
          b0 -= (1.0 - mipsolver->mipdata_->feastol);
        else
          b0 += (1.0 - mipsolver->mipdata_->feastol);
      } else {
        // for a continuous variable we relax the bound by epsilon to
        // accomodate for tiny rounding errors
        if (domchgstack_[pos].boundtype == HighsBoundType::kLower)
          b0 -= mipsolver->mipdata_->epsilon;
        else
          b0 += mipsolver->mipdata_->epsilon;
      }

      // now multiply the bound constraint with the coefficient value in the
      // constraint to obtain b0
      b0 *= domchgVal;

      // compute the lower bound of M that is necessary
      double Mlower = mipsolver->rowUpper(rowIndex) - b0;

      // M is the global residual activity initially
      double M = minAct;
      if (domchgVal < 0)
        M -= domchgVal * globaldom.colUpper_[domchgstack_[pos].column];
      else
        M -= domchgVal * globaldom.colLower_[domchgstack_[pos].column];

      reasonDomChgs.clear();
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];
        if (col == domchgstack_[pos].column) continue;

        HighsInt boundpos;
        if (vals[i] > 0) {
          double lb = getColLowerPos(col, pos, boundpos);
          if (globaldom.colLower_[col] >= lb) continue;
          M += vals[i] * (lb - globaldom.colLower_[col]);
        } else {
          double ub = getColUpperPos(col, pos, boundpos);
          if (globaldom.colUpper_[col] <= ub) continue;
          M += vals[i] * (ub - globaldom.colUpper_[col]);
        }

        reasonDomChgs.push_back(boundpos);
        if (M > Mlower) break;
      }

      break;
    }
    default: {
      assert(domchgreason_[pos].type >= 0);
      assert(domchgreason_[pos].type < cutpoolpropagation.size());
      HighsInt cutpoolIndex = domchgreason_[pos].type;
      HighsInt cutIndex = domchgreason_[pos].index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      cutpoolpropagation[cutpoolIndex].cutpool->getCut(cutIndex, len, inds,
                                                       vals);

      // get the coefficient value of the column for which we want to explain
      // the bound change
      double domchgVal = 0;
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];
        if (col == domchgstack_[pos].column) {
          domchgVal = vals[i];
          break;
        }
      }

      assert(domchgVal != 0);

      double minAct = globaldom.getMinCutActivity(
          *cutpoolpropagation[cutpoolIndex].cutpool, cutIndex);
      // if the global min activity is infinite we do not filter the bounds
      // and just take all as it is rare and we would need to use the local
      // activity as base which is not available as easily as not all bound
      // changes on the stack are used. Using all local bounds of columns that
      // are present in the cut is however guaranteed to be sufficient even
      // though it might not be as strong.
      if (minAct == -kHighsInf) return false;

      // to explain the bound change we start from the bound constraint,
      // multiply it by the columns coefficient in the constraint. Then the
      // bound constraint is a0 * x0 <= b0 and the constraint
      // a0 * x0 + \sum ai * xi <= b. Let the min activity of \sum ai xi be M,
      // then the constraint yields the bound a0 * x0 <= b - M. If M is
      // sufficiently large the bound constraint is implied. Therefore we
      // increase M by updating it with the stronger local bounds until
      // M >= b - b0 holds.
      double b0 = domchgstack_[pos].boundval;
      if (mipsolver->variableType(domchgstack_[pos].column) !=
          HighsVarType::kContinuous) {
        // in case of an integral variable the bound was rounded and can be
        // relaxed by 1-feastol
        if (domchgstack_[pos].boundtype == HighsBoundType::kLower)
          b0 -= (1.0 - mipsolver->mipdata_->feastol);
        else
          b0 += (1.0 - mipsolver->mipdata_->feastol);
      } else {
        // for a continuous variable we relax the bound by epsilon to
        // accomodate for tiny rounding errors
        if (domchgstack_[pos].boundtype == HighsBoundType::kLower)
          b0 -= mipsolver->mipdata_->epsilon;
        else
          b0 += mipsolver->mipdata_->epsilon;
      }

      // now multiply the bound constraint with the coefficient value in the
      // constraint to obtain b0
      b0 *= domchgVal;

      // compute the lower bound of M that is necessary
      double Mlower =
          cutpoolpropagation[cutpoolIndex].cutpool->getRhs()[cutIndex] - b0;

      // M is the global residual activity initially
      double M = minAct;
      if (domchgVal < 0)
        M -= domchgVal * globaldom.colUpper_[domchgstack_[pos].column];
      else
        M -= domchgVal * globaldom.colLower_[domchgstack_[pos].column];

      reasonDomChgs.clear();
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];
        if (col == domchgstack_[pos].column) continue;

        HighsInt boundpos;
        if (vals[i] > 0) {
          double lb = getColLowerPos(col, pos, boundpos);
          if (globaldom.colLower_[col] >= lb) continue;
          M += vals[i] * (lb - globaldom.colLower_[col]);
        } else {
          double ub = getColUpperPos(col, pos, boundpos);
          if (globaldom.colUpper_[col] <= ub) continue;
          M += vals[i] * (ub - globaldom.colUpper_[col]);
        }

        reasonDomChgs.push_back(boundpos);
        if (M > Mlower) break;
      }
    }
  }
  return true;
}

bool HighsDomain::explainInfeasibility(std::vector<HighsInt>& reasonDomChgs) {
  HighsDomain& globaldom = mipsolver->mipdata_->domain;
  switch (infeasible_reason.type) {
    case Reason::kUnknown:
    case Reason::kBranching:
      return false;
    case Reason::kCliqueTable: {
      HighsInt col = infeasible_reason.index >> 1;
      HighsInt val = infeasible_reason.index & 1;
      reasonDomChgs.clear();

      reasonDomChgs.push_back(infeasible_pos);
      if (val) {
        assert(colLowerPos_[col] >= 0);
        assert(colLowerPos_[col] < domchgstack_.size());
        assert(reasonDomChgs.back() != colLowerPos_[col]);
        reasonDomChgs.push_back(colLowerPos_[col]);
      } else {
        assert(colUpperPos_[col] >= 0);
        assert(colUpperPos_[col] < domchgstack_.size());
        assert(reasonDomChgs.back() != colUpperPos_[col]);
        reasonDomChgs.push_back(colUpperPos_[col]);
      }

      break;
    }
    case Reason::kModelRowLower: {
      HighsInt rowIndex = infeasible_reason.index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double maxAct = globaldom.getMaxActivity(rowIndex);
      // if the global min activity is infinite we do not filter the bounds
      // and just take all as it is rare and we would need to use the local
      // activity as base which is not available as easily as not all bound
      // changes on the stack are used. Using all local bounds of columns that
      // are present in the cut is however guaranteed to be sufficient even
      // though it might not be as strong.
      if (maxAct == kHighsInf) return false;

      // compute the lower bound of M that is necessary
      double Mupper = mipsolver->rowLower(rowIndex) -
                      std::max(1.0, std::abs(mipsolver->rowLower(rowIndex))) *
                          mipsolver->mipdata_->feastol;

      // M is the global residual activity initially
      double M = maxAct;

      reasonDomChgs.clear();
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];

        HighsInt boundpos;
        if (vals[i] > 0) {
          double ub = getColUpperPos(col, infeasible_pos, boundpos);
          if (globaldom.colUpper_[col] <= ub) continue;
          M += vals[i] * (ub - globaldom.colUpper_[col]);
        } else {
          double lb = getColLowerPos(col, infeasible_pos, boundpos);
          if (globaldom.colLower_[col] >= lb) continue;
          M += vals[i] * (lb - globaldom.colLower_[col]);
        }

        assert(boundpos >= 0);
        assert(boundpos < domchgstack_.size());

        reasonDomChgs.push_back(boundpos);
        if (M < Mupper) break;
      }

      break;
    }
    case Reason::kModelRowUpper: {
      HighsInt rowIndex = infeasible_reason.index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double minAct = globaldom.getMinActivity(rowIndex);
      // if the global min activity is infinite we do not filter the bounds
      // and just take all as it is rare and we would need to use the local
      // activity as base which is not available as easily as not all bound
      // changes on the stack are used. Using all local bounds of columns that
      // are present in the cut is however guaranteed to be sufficient even
      // though it might not be as strong.
      if (minAct == -kHighsInf) return false;

      // compute the lower bound of M that is necessary
      double Mlower = mipsolver->rowUpper(rowIndex) +
                      std::max(1.0, std::abs(mipsolver->rowUpper(rowIndex))) *
                          mipsolver->mipdata_->feastol;

      // M is the global residual activity initially
      double M = minAct;
      reasonDomChgs.clear();
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];

        HighsInt boundpos;
        if (vals[i] > 0) {
          double lb = getColLowerPos(col, infeasible_pos, boundpos);
          if (globaldom.colLower_[col] >= lb) continue;
          M += vals[i] * (lb - globaldom.colLower_[col]);
        } else {
          double ub = getColUpperPos(col, infeasible_pos, boundpos);
          if (globaldom.colUpper_[col] <= ub) continue;
          M += vals[i] * (ub - globaldom.colUpper_[col]);
        }

        assert(boundpos >= 0);
        assert(boundpos < domchgstack_.size());

        reasonDomChgs.push_back(boundpos);
        if (M > Mlower) break;
      }

      break;
    }
    default: {
      assert(infeasible_reason.type >= 0);
      assert(infeasible_reason.type < cutpoolpropagation.size());
      HighsInt cutpoolIndex = infeasible_reason.type;
      HighsInt cutIndex = infeasible_reason.index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      cutpoolpropagation[cutpoolIndex].cutpool->getCut(cutIndex, len, inds,
                                                       vals);

      double minAct = globaldom.getMinCutActivity(
          *cutpoolpropagation[cutpoolIndex].cutpool, cutIndex);
      // if the global min activity is infinite we do not filter the bounds
      // and just take all as it is rare and we would need to use the local
      // activity as base which is not available as easily as not all bound
      // changes on the stack are used. Using all local bounds of columns that
      // are present in the cut is however guaranteed to be sufficient even
      // though it might not be as strong.
      if (minAct == -kHighsInf) return false;

      // compute the lower bound of M that is necessary
      double Mlower =
          cutpoolpropagation[cutpoolIndex].cutpool->getRhs()[cutIndex] +
          std::max(1.0, std::abs(cutpoolpropagation[cutpoolIndex]
                                     .cutpool->getRhs()[cutIndex])) *
              mipsolver->mipdata_->feastol;

      // M is the global residual activity initially
      double M = minAct;
      reasonDomChgs.clear();
      for (HighsInt i = 0; i < len; ++i) {
        HighsInt col = inds[i];

        HighsInt boundpos;
        if (vals[i] > 0) {
          double lb = getColLowerPos(col, infeasible_pos, boundpos);
          if (globaldom.colLower_[col] >= lb) continue;
          M += vals[i] * (lb - globaldom.colLower_[col]);
        } else {
          double ub = getColUpperPos(col, infeasible_pos, boundpos);
          if (globaldom.colUpper_[col] <= ub) continue;
          M += vals[i] * (ub - globaldom.colUpper_[col]);
        }

        assert(boundpos >= 0);
        assert(boundpos < domchgstack_.size());

        reasonDomChgs.push_back(boundpos);
        if (M > Mlower) break;
      }
    }
  }
  return true;
}

void HighsDomain::conflictAnalysis() {
  std::set<HighsInt> reasonSideFrontier;
  std::priority_queue<HighsInt> resolveQueue;
  std::vector<HighsInt> reasonDomChgs;
  std::vector<HighsInt> inds;
  std::vector<double> vals;
  if (!infeasible_ || branchPos_.empty() ||
      infeasible_reason.type == Reason::kBranching ||
      infeasible_reason.type == Reason::kUnknown)
    return;

  reasonDomChgs.reserve(domchgstack_.size());
  if (!explainInfeasibility(reasonDomChgs)) return;

  HighsDomain& globaldom = mipsolver->mipdata_->domain;

  for (HighsInt domChgPos : reasonDomChgs) {
    if (reasonSideFrontier.insert(domChgPos).second &&
        !globaldom.isBinary(domchgstack_[domChgPos].column)) {
      resolveQueue.push(domChgPos);
    }
  }

  while (!resolveQueue.empty()) {
    HighsInt resolvePos = resolveQueue.top();
    resolveQueue.pop();
    if (!explainBoundChange(resolvePos, reasonDomChgs)) return;

    reasonSideFrontier.erase(resolvePos);
    for (HighsInt domChgPos : reasonDomChgs) {
      if (reasonSideFrontier.insert(domChgPos).second) {
        if (!globaldom.isBinary(domchgstack_[domChgPos].column))
          resolveQueue.push(domChgPos);
      }
    }
  }

  size_t maxSize = reasonSideFrontier.size();
  HighsInt currDepthEnd = kHighsIInf;

  for (HighsInt currDepth = branchPos_.size() - 1; currDepth >= 0;
       --currDepth) {
    for (auto it = reasonSideFrontier.lower_bound(branchPos_[currDepth] + 1);
         it != reasonSideFrontier.upper_bound(currDepthEnd); ++it)
      resolveQueue.push(*it);

    currDepthEnd = branchPos_[currDepth];

    if (resolveQueue.empty()) continue;

    size_t thisDepthSize = resolveQueue.size();

    size_t minQueueSize = 0;

    while (thisDepthSize > minQueueSize) {
      HighsInt resolvePos = resolveQueue.top();
      resolveQueue.pop();

      if (resolvePos > branchPos_[currDepth]) --thisDepthSize;

      if (explainBoundChange(resolvePos, reasonDomChgs)) {
        minQueueSize = 1;
        reasonSideFrontier.erase(resolvePos);
        for (HighsInt domChgPos : reasonDomChgs) {
          if (reasonSideFrontier.insert(domChgPos).second) {
            if (!globaldom.isBinary(domchgstack_[domChgPos].column) ||
                domChgPos > branchPos_[currDepth]) {
              if (domChgPos > branchPos_[currDepth]) ++thisDepthSize;
              resolveQueue.push(domChgPos);
            }
          }
        }
      } else if (!globaldom.isBinary(domchgstack_[resolvePos].column))
        return;
    }

    while (!resolveQueue.empty()) {
      HighsInt resolvePos = resolveQueue.top();
      resolveQueue.pop();
      if (globaldom.isBinary(domchgstack_[resolvePos].column)) continue;
      if (!explainBoundChange(resolvePos, reasonDomChgs)) return;

      reasonSideFrontier.erase(resolvePos);
      for (HighsInt domChgPos : reasonDomChgs) {
        if (reasonSideFrontier.insert(domChgPos).second) {
          if (!globaldom.isBinary(domchgstack_[domChgPos].column))
            resolveQueue.push(domChgPos);
        }
      }
    }

    if (minQueueSize != 0 && reasonSideFrontier.size() < maxSize) {
#if 0
      printf("UIP cut of depth %d:\n", currDepth);
      bool sep = false;
      for (HighsInt i : reasonSideFrontier) {
        char op =
            domchgstack_[i].boundtype == HighsBoundType::kLower ? '<' : '>';
        if (sep) {
          printf(" V (x%d %c %g)", domchgstack_[i].column, op,
                 domchgstack_[i].boundval);
        } else {
          printf("(x%d %c %g)", domchgstack_[i].column, op,
                 domchgstack_[i].boundval);
          sep = true;
        }
      }
      printf("\n");
#endif
      inds.clear();
      vals.clear();
      HighsInt rhs = reasonSideFrontier.size() - 1;
      for (HighsInt i : reasonSideFrontier) {
        HighsInt col = domchgstack_[i].column;
        assert(globaldom.isBinary(col));
        inds.push_back(col);
        if (domchgstack_[i].boundtype == HighsBoundType::kLower) {
          assert(domchgstack_[i].boundval == 1.0);
          vals.push_back(1.0);
        } else {
          vals.push_back(-1.0);
          rhs -= 1;
        }
      }

      mipsolver->mipdata_->cutpool.addCut(*mipsolver, inds.data(), vals.data(),
                                          inds.size(), rhs, true);
    }
  }
}

void HighsDomain::tightenCoefficients(HighsInt* inds, double* vals,
                                      HighsInt len, double& rhs) const {
  HighsCDouble maxactivity = 0;

  for (HighsInt i = 0; i != len; ++i) {
    if (vals[i] > 0) {
      if (colUpper_[inds[i]] == kHighsInf) return;

      maxactivity += colUpper_[inds[i]] * vals[i];
    } else {
      if (colLower_[inds[i]] == -kHighsInf) return;

      maxactivity += colLower_[inds[i]] * vals[i];
    }
  }

  HighsCDouble maxabscoef = maxactivity - rhs;
  if (maxabscoef > mipsolver->mipdata_->feastol) {
    HighsCDouble upper = rhs;
    HighsInt tightened = 0;
    for (HighsInt i = 0; i != len; ++i) {
      if (mipsolver->variableType(inds[i]) == HighsVarType::kContinuous)
        continue;
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
      // printf("tightened %" HIGHSINT_FORMAT " coefficients, rhs changed from
      // %g to %g\n",
      //       tightened, rhs, double(upper));
      rhs = double(upper);
    }
  }
}

double HighsDomain::getMinCutActivity(const HighsCutPool& cutpool,
                                      HighsInt cut) {
  for (auto& cutpoolprop : cutpoolpropagation) {
    if (cutpoolprop.cutpool == &cutpool) {
      if (cutpool.getModificationCount(cut) !=
          cutpoolprop.activitycutversion_[cut]) {
        cutpoolprop.activitycutversion_[cut] =
            cutpoolprop.cutpool->getModificationCount(cut);
        HighsInt start = cutpoolprop.cutpool->getMatrix().getRowStart(cut);
        if (start == -1) {
          cutpoolprop.activitycuts_[cut] = 0;
          return -kHighsInf;
        }
        HighsInt end = cutpoolprop.cutpool->getMatrix().getRowEnd(cut);
        const HighsInt* arindex = cutpoolprop.cutpool->getMatrix().getARindex();
        const double* arvalue = cutpoolprop.cutpool->getMatrix().getARvalue();
        computeMinActivity(start, end, arindex, arvalue,
                           cutpoolprop.activitycutsinf_[cut],
                           cutpoolprop.activitycuts_[cut]);
      }
      return cutpoolprop.activitycutsinf_[cut] == 0
                 ? double(cutpoolprop.activitycuts_[cut])
                 : -kHighsInf;
    }
  }

  return -kHighsInf;
}

bool HighsDomain::isFixing(const HighsDomainChange& domchg) const {
  double otherbound = domchg.boundtype == HighsBoundType::kUpper
                          ? colLower_[domchg.column]
                          : colUpper_[domchg.column];
  return std::abs(domchg.boundval - otherbound) <= mipsolver->mipdata_->epsilon;
}
