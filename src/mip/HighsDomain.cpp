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

#include "mip/HighsConflictPool.h"
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
  col_lower_ = mipsolver.model_->col_lower_;
  col_upper_ = mipsolver.model_->col_upper_;
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

void HighsDomain::addConflictPool(HighsConflictPool& conflictPool) {
  HighsInt conflictPoolIndex = conflictPoolPropagation.size();
  conflictPoolPropagation.emplace_back(conflictPoolIndex, this, conflictPool);
}

void HighsDomain::ConflictPoolPropagation::linkWatchedLiteral(
    HighsInt linkPos) {
  assert(watchedLiterals_[linkPos].domchg.column != -1);
  HighsInt& head =
      watchedLiterals_[linkPos].domchg.boundtype == HighsBoundType::kLower
          ? colLowerWatched_[watchedLiterals_[linkPos].domchg.column]
          : colUpperWatched_[watchedLiterals_[linkPos].domchg.column];

  watchedLiterals_[linkPos].prev = -1;
  watchedLiterals_[linkPos].next = head;
  if (head != -1) {
    watchedLiterals_[head].prev = linkPos;
    head = linkPos;
  }
}

void HighsDomain::ConflictPoolPropagation::unlinkWatchedLiteral(
    HighsInt linkPos) {
  if (watchedLiterals_[linkPos].domchg.column == -1) return;

  HighsInt& head =
      watchedLiterals_[linkPos].domchg.boundtype == HighsBoundType::kLower
          ? colLowerWatched_[watchedLiterals_[linkPos].domchg.column]
          : colUpperWatched_[watchedLiterals_[linkPos].domchg.column];
  watchedLiterals_[linkPos].domchg.column = -1;
  HighsInt prev = watchedLiterals_[linkPos].prev;
  HighsInt next = watchedLiterals_[linkPos].next;
  if (prev != -1)
    watchedLiterals_[prev].next = next;
  else
    head = next;

  if (next != -1) watchedLiterals_[next].prev = prev;
}

HighsDomain::ConflictPoolPropagation::ConflictPoolPropagation(
    HighsInt conflictpoolindex, HighsDomain* domain,
    HighsConflictPool& conflictpool_)
    : conflictpoolindex(conflictpoolindex),
      domain(domain),
      conflictpool_(&conflictpool_) {
  colLowerWatched_.resize(domain->mipsolver->numCol(), -1);
  colUpperWatched_.resize(domain->mipsolver->numCol(), -1);
  conflictpool_.addPropagationDomain(this);
}

HighsDomain::ConflictPoolPropagation::ConflictPoolPropagation(
    const ConflictPoolPropagation& other)
    : conflictpoolindex(other.conflictpoolindex),
      domain(other.domain),
      conflictpool_(other.conflictpool_),
      colLowerWatched_(other.colLowerWatched_),
      colUpperWatched_(other.colUpperWatched_),
      conflictFlag_(other.conflictFlag_),
      propagateConflictInds_(other.propagateConflictInds_),
      watchedLiterals_(other.watchedLiterals_) {
  conflictpool_->addPropagationDomain(this);
}

HighsDomain::ConflictPoolPropagation::~ConflictPoolPropagation() {
  conflictpool_->removePropagationDomain(this);
}

void HighsDomain::ConflictPoolPropagation::conflictDeleted(HighsInt conflict) {
  conflictFlag_[conflict] |= 8;
  unlinkWatchedLiteral(2 * conflict);
  unlinkWatchedLiteral(2 * conflict + 1);
}

void HighsDomain::ConflictPoolPropagation::conflictAdded(HighsInt conflict) {
  HighsInt start = conflictpool_->getConflictRanges()[conflict].first;
  HighsInt end = conflictpool_->getConflictRanges()[conflict].second;
  const std::vector<HighsDomainChange>& conflictEntries =
      conflictpool_->getConflictEntryVector();

  if (HighsInt(conflictFlag_.size()) <= conflict) {
    watchedLiterals_.resize(2 * conflict + 2);
    conflictFlag_.resize(conflict + 1);
  }

  HighsInt numWatched = 0;
  for (HighsInt i = start; i != end; ++i) {
    if (domain->isActive(conflictEntries[i])) continue;
    HighsInt col = conflictEntries[i].column;
    HighsInt watchPos = 2 * conflict + numWatched;
    watchedLiterals_[watchPos].domchg = conflictEntries[i];
    linkWatchedLiteral(watchPos);
    if (++numWatched == 2) break;
  }
  switch (numWatched) {
    case 0: {
      std::pair<HighsInt, HighsInt> latestActive[2];
      HighsInt numActive = 0;
      for (HighsInt i = start; i != end; ++i) {
        HighsInt pos = conflictEntries[i].boundtype == HighsBoundType::kLower
                           ? domain->colLowerPos_[conflictEntries[i].column]
                           : domain->colUpperPos_[conflictEntries[i].column];
        switch (numActive) {
          case 0:
            latestActive[0] = std::make_pair(pos, i);
            numActive = 1;
            break;
          case 1:
            latestActive[1] = std::make_pair(pos, i);
            numActive = 2;
            if (latestActive[0].first < latestActive[1].first)
              std::swap(latestActive[0], latestActive[1]);
            break;
          case 2:
            if (pos > latestActive[1].first) {
              latestActive[1] = std::make_pair(pos, i);
              if (latestActive[0].first < latestActive[1].first)
                std::swap(latestActive[0], latestActive[1]);
            }
        }
      }
      for (HighsInt i = 0; i < numActive; ++i) {
        HighsInt watchPos = 2 * conflict + i;
        watchedLiterals_[watchPos].domchg =
            conflictEntries[latestActive[i].second];
        linkWatchedLiteral(watchPos);
      }
      break;
    }
    case 1: {
      HighsInt latestActive = -1;
      HighsInt latestPos = -1;

      for (HighsInt i = start; i != end; ++i) {
        HighsInt pos = conflictEntries[i].boundtype == HighsBoundType::kLower
                           ? domain->colLowerPos_[conflictEntries[i].column]
                           : domain->colUpperPos_[conflictEntries[i].column];
        if (pos > latestPos) {
          latestActive = i;
          latestPos = pos;
        }
      }
      if (latestActive != -1) {
        HighsInt watchPos = 2 * conflict + 1;
        watchedLiterals_[watchPos].domchg = conflictEntries[latestActive];
        linkWatchedLiteral(watchPos);
      }
      break;
    }
    case 2:
      break;
  }

  conflictFlag_[conflict] = numWatched | (conflictFlag_[conflict] & 4);
  markPropagateConflict(conflict);
}

void HighsDomain::ConflictPoolPropagation::markPropagateConflict(
    HighsInt conflict) {
  if (conflictFlag_[conflict] < 2) {
    propagateConflictInds_.push_back(conflict);
    conflictFlag_[conflict] |= 4;
  }
}

void HighsDomain::ConflictPoolPropagation::updateActivityLbChange(
    HighsInt col, double oldbound, double newbound) {
  assert(!domain->infeasible_);

  const std::vector<HighsDomainChange>& conflictEntries =
      conflictpool_->getConflictEntryVector();

  for (HighsInt i = colLowerWatched_[col]; i != -1;
       i = watchedLiterals_[i].next) {
    HighsInt conflict = i >> 1;

    const HighsDomainChange& domchg = watchedLiterals_[i].domchg;
    HighsInt numInactiveDelta =
        (domchg.boundval > newbound) - (domchg.boundval > oldbound);
    if (numInactiveDelta != 0) {
      conflictFlag_[conflict] += numInactiveDelta;
      markPropagateConflict(conflict);
    }
  }
}

void HighsDomain::ConflictPoolPropagation::updateActivityUbChange(
    HighsInt col, double oldbound, double newbound) {
  assert(!domain->infeasible_);

  const std::vector<HighsDomainChange>& conflictEntries =
      conflictpool_->getConflictEntryVector();

  for (HighsInt i = colUpperWatched_[col]; i != -1;
       i = watchedLiterals_[i].next) {
    HighsInt conflict = i >> 1;

    const HighsDomainChange& domchg = watchedLiterals_[i].domchg;
    HighsInt numInactiveDelta =
        (domchg.boundval < newbound) - (domchg.boundval < oldbound);
    if (numInactiveDelta != 0) {
      conflictFlag_[conflict] += numInactiveDelta;
      markPropagateConflict(conflict);
    }
  }
}

void HighsDomain::ConflictPoolPropagation::propagateConflict(
    HighsInt conflict) {
  // remove propagate flag, but keep watched and deleted information
  conflictFlag_[conflict] &= (3 | 8);
  // if two inactive literals are watched or conflict has been deleted skip
  if (conflictFlag_[conflict] >= 2) return;

  if (domain->infeasible_) return;

  const std::vector<HighsDomainChange>& entries =
      conflictpool_->getConflictEntryVector();
  HighsInt start = conflictpool_->getConflictRanges()[conflict].first;
  if (start == -1) {
    unlinkWatchedLiteral(2 * conflict);
    unlinkWatchedLiteral(2 * conflict + 1);
    return;
  }
  HighsInt end = conflictpool_->getConflictRanges()[conflict].second;

  WatchedLiteral* watched = watchedLiterals_.data() + 2 * conflict;

  HighsInt inactive[2];
  HighsInt latestactive[2];
  HighsInt numInactive = 0;
  for (HighsInt i = start; i != end; ++i) {
    if (domain->isActive(entries[i])) continue;

    inactive[numInactive++] = i;
    if (numInactive == 2) break;
  }

  conflictFlag_[conflict] = numInactive;

  switch (numInactive) {
    case 0:
      assert(!domain->infeasible_);
      domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
      domain->infeasible_ = true;
      domain->infeasible_reason = Reason::cut(
          domain->cutpoolpropagation.size() + conflictpoolindex, conflict);
      domain->infeasible_pos = domain->domchgstack_.size();
      conflictpool_->resetAge(conflict);
      // printf("conflict propagation found infeasibility\n");
      break;
    case 1: {
      HighsDomainChange domchg = domain->flip(entries[inactive[0]]);
      if (!domain->isActive(domchg)) {
        domain->changeBound(
            domain->flip(entries[inactive[0]]),
            Reason::cut(domain->cutpoolpropagation.size() + conflictpoolindex,
                        conflict));
        conflictpool_->resetAge(conflict);
      }
      // printf("conflict propagation found bound change\n");
      break;
    }
    case 2: {
      if (watched[0].domchg != entries[inactive[0]]) {
        unlinkWatchedLiteral(2 * conflict);
        watched[0].domchg = entries[inactive[0]];
        linkWatchedLiteral(2 * conflict);
      }

      if (watched[1].domchg != entries[inactive[1]]) {
        unlinkWatchedLiteral(2 * conflict + 1);
        watched[1].domchg = entries[inactive[1]];
        linkWatchedLiteral(2 * conflict + 1);
      }

      return;
    }
  }
}

HighsDomain::CutpoolPropagation::CutpoolPropagation(HighsInt cutpoolindex,
                                                    HighsDomain* domain,
                                                    HighsCutPool& cutpool_)
    : cutpoolindex(cutpoolindex), domain(domain), cutpool(&cutpool_) {
  cutpool->addPropagationDomain(this);
}

HighsDomain::CutpoolPropagation::CutpoolPropagation(
    const CutpoolPropagation& other)
    : cutpoolindex(other.cutpoolindex),
      domain(other.domain),
      cutpool(other.cutpool),
      activitycuts_(other.activitycuts_),
      activitycutsinf_(other.activitycutsinf_),
      propagatecutflags_(other.propagatecutflags_),
      propagatecutinds_(other.propagatecutinds_),
      capacityThreshold_(other.capacityThreshold_) {
  cutpool->addPropagationDomain(this);
}

HighsDomain::CutpoolPropagation::~CutpoolPropagation() {
  cutpool->removePropagationDomain(this);
}

void HighsDomain::CutpoolPropagation::recomputeCapacityThreshold(HighsInt cut) {
  HighsInt start = cutpool->getMatrix().getRowStart(cut);
  HighsInt end = cutpool->getMatrix().getRowEnd(cut);
  const HighsInt* arindex = cutpool->getMatrix().getARindex();
  const double* arvalue = cutpool->getMatrix().getARvalue();
  capacityThreshold_[cut] = 0.0;
  for (HighsInt i = start; i < end; ++i) {
    if (domain->col_upper_[arindex[i]] == domain->col_lower_[arindex[i]])
      continue;

    double boundRange =
        domain->col_upper_[arindex[i]] - domain->col_lower_[arindex[i]];

    boundRange -= domain->variableType(arindex[i]) == HighsVarType::kContinuous
                      ? std::max(0.3 * boundRange,
                                 1000.0 * domain->mipsolver->mipdata_->feastol)
                      : domain->mipsolver->mipdata_->feastol;

    double threshold = std::abs(arvalue[i]) * boundRange;

    capacityThreshold_[cut] = std::max({capacityThreshold_[cut], threshold,
                                        domain->mipsolver->mipdata_->feastol});
  }
}

void HighsDomain::CutpoolPropagation::cutAdded(HighsInt cut, bool propagate) {
  if (!propagate) {
    if (domain != &domain->mipsolver->mipdata_->domain) return;
    HighsInt start = cutpool->getMatrix().getRowStart(cut);
    HighsInt end = cutpool->getMatrix().getRowEnd(cut);
    const HighsInt* arindex = cutpool->getMatrix().getARindex();
    const double* arvalue = cutpool->getMatrix().getARvalue();

    if (HighsInt(activitycuts_.size()) <= cut) {
      activitycuts_.resize(cut + 1);
      activitycutsinf_.resize(cut + 1);
      propagatecutflags_.resize(cut + 1, 2);
      capacityThreshold_.resize(cut + 1);
    }

    propagatecutflags_[cut] &= ~uint8_t{2};
    domain->computeMinActivity(start, end, arindex, arvalue,
                               activitycutsinf_[cut], activitycuts_[cut]);
  } else {
    HighsInt start = cutpool->getMatrix().getRowStart(cut);
    HighsInt end = cutpool->getMatrix().getRowEnd(cut);
    const HighsInt* arindex = cutpool->getMatrix().getARindex();
    const double* arvalue = cutpool->getMatrix().getARvalue();

    if (HighsInt(activitycuts_.size()) <= cut) {
      activitycuts_.resize(cut + 1);
      activitycutsinf_.resize(cut + 1);
      propagatecutflags_.resize(cut + 1, 2);
      capacityThreshold_.resize(cut + 1);
    }

    propagatecutflags_[cut] &= ~uint8_t{2};
    domain->computeMinActivity(start, end, arindex, arvalue,
                               activitycutsinf_[cut], activitycuts_[cut]);

    recomputeCapacityThreshold(cut);
    markPropagateCut(cut);
  }
}

void HighsDomain::CutpoolPropagation::cutDeleted(
    HighsInt cut, bool deletedOnlyForPropagation) {
  if (deletedOnlyForPropagation &&
      domain == &domain->mipsolver->mipdata_->domain) {
    assert(domain->branchPos_.empty());
    return;
  }

  if (cut < (HighsInt)propagatecutflags_.size()) propagatecutflags_[cut] |= 2;
}

void HighsDomain::CutpoolPropagation::markPropagateCut(HighsInt cut) {
  if (!propagatecutflags_[cut] &&
      (activitycutsinf_[cut] == 1 ||
       (cutpool->getRhs()[cut] - double(activitycuts_[cut]) <=
        capacityThreshold_[cut]))) {
    propagatecutinds_.push_back(cut);
    propagatecutflags_[cut] |= 1;
  }
}

void HighsDomain::CutpoolPropagation::updateActivityLbChange(HighsInt col,
                                                             double oldbound,
                                                             double newbound) {
  assert(!domain->infeasible_);

  if (newbound < oldbound) {
    cutpool->getMatrix().forEachNegativeColumnEntry(
        col, [&](HighsInt row, double val) {
          domain->updateThresholdLbChange(col, newbound, val,
                                          capacityThreshold_[row]);
          return true;
        });
  }

  cutpool->getMatrix().forEachPositiveColumnEntry(
      col, [&](HighsInt row, double val) {
        assert(val > 0);
        double deltamin;

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

        if (deltamin <= 0) {
          domain->updateThresholdLbChange(col, newbound, val,
                                          capacityThreshold_[row]);
          return true;
        }

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

        markPropagateCut(row);

        return true;
      });

  if (domain->infeasible_) {
    assert(domain->infeasible_reason.type == cutpoolindex);
    assert(domain->infeasible_reason.index >= 0);
    std::swap(oldbound, newbound);
    cutpool->getMatrix().forEachPositiveColumnEntry(
        col, [&](HighsInt row, double val) {
          assert(val > 0);
          double deltamin;

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

          return true;
        });
  }
}

void HighsDomain::CutpoolPropagation::updateActivityUbChange(HighsInt col,
                                                             double oldbound,
                                                             double newbound) {
  assert(!domain->infeasible_);

  if (newbound > oldbound) {
    cutpool->getMatrix().forEachNegativeColumnEntry(
        col, [&](HighsInt row, double val) {
          domain->updateThresholdUbChange(col, newbound, val,
                                          capacityThreshold_[row]);
          return true;
        });
  }

  cutpool->getMatrix().forEachNegativeColumnEntry(
      col, [&](HighsInt row, double val) {
        assert(val < 0);
        double deltamin;

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

        if (deltamin <= 0) {
          domain->updateThresholdUbChange(col, newbound, val,
                                          capacityThreshold_[row]);
          return true;
        }

        if (activitycutsinf_[row] == 0 &&
            activitycuts_[row] - cutpool->getRhs()[row] >
                domain->mipsolver->mipdata_->feastol) {
          domain->mipsolver->mipdata_->debugSolution.nodePruned(*domain);
          domain->infeasible_ = true;
          domain->infeasible_pos = domain->domchgstack_.size();
          domain->infeasible_reason = Reason::cut(cutpoolindex, row);
          return false;
        }

        markPropagateCut(row);

        return true;
      });

  if (domain->infeasible_) {
    assert(domain->infeasible_reason.type == cutpoolindex);
    assert(domain->infeasible_reason.index >= 0);
    std::swap(oldbound, newbound);
    cutpool->getMatrix().forEachNegativeColumnEntry(
        col, [&](HighsInt row, double val) {
          assert(val < 0);
          double deltamin;

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

          return true;
        });
  }
}

void HighsDomain::computeMinActivity(HighsInt start, HighsInt end,
                                     const HighsInt* ARindex,
                                     const double* ARvalue, HighsInt& ninfmin,
                                     HighsCDouble& activitymin) {
  if (infeasible_) {
    activitymin = 0.0;
    ninfmin = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < int(col_lower_.size()));

      HighsInt tmp;
      double lb = getColLowerPos(col, infeasible_pos - 1, tmp);
      double ub = getColUpperPos(col, infeasible_pos - 1, tmp);
      double contributionmin = activityContributionMin(val, lb, ub);

      if (contributionmin == -kHighsInf)
        ++ninfmin;
      else
        activitymin += contributionmin;
    }

    activitymin.renormalize();
  } else {
    activitymin = 0.0;
    ninfmin = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < int(col_lower_.size()));

      double contributionmin =
          activityContributionMin(val, col_lower_[col], col_upper_[col]);

      if (contributionmin == -kHighsInf)
        ++ninfmin;
      else
        activitymin += contributionmin;
    }

    activitymin.renormalize();
  }
}

void HighsDomain::computeMaxActivity(HighsInt start, HighsInt end,
                                     const HighsInt* ARindex,
                                     const double* ARvalue, HighsInt& ninfmax,
                                     HighsCDouble& activitymax) {
  if (infeasible_) {
    activitymax = 0.0;
    ninfmax = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < int(col_lower_.size()));

      HighsInt tmp;
      double lb = getColLowerPos(col, infeasible_pos - 1, tmp);
      double ub = getColUpperPos(col, infeasible_pos - 1, tmp);
      double contributionmin = activityContributionMax(val, lb, ub);

      if (contributionmin == kHighsInf)
        ++ninfmax;
      else
        activitymax += contributionmin;
    }

    activitymax.renormalize();
  } else {
    activitymax = 0.0;
    ninfmax = 0;
    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex[j];
      double val = ARvalue[j];

      assert(col < int(col_lower_.size()));

      double contributionmin =
          activityContributionMax(val, col_lower_[col], col_upper_[col]);

      if (contributionmin == kHighsInf)
        ++ninfmax;
      else
        activitymax += contributionmin;
    }

    activitymax.renormalize();
  }
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
        Rvalue[i], col_lower_[Rindex[i]], col_upper_[Rindex[i]]);
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
        if (bound < col_upper_[Rindex[i]] &&
            col_upper_[Rindex[i]] - bound >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(double(boundVal) - col_lower_[Rindex[i]]) <=
            mipsolver->mipdata_->epsilon)
          bound = col_lower_[Rindex[i]];
        else
          bound = double(boundVal);
        if (col_upper_[Rindex[i]] == kHighsInf)
          accept = true;
        else if (bound + 1000.0 * mipsolver->mipdata_->feastol <
                 col_upper_[Rindex[i]]) {
          double relativeImprove = col_upper_[Rindex[i]] - bound;
          if (col_lower_[Rindex[i]] != -kHighsInf)
            relativeImprove /= col_upper_[Rindex[i]] - col_lower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(col_upper_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kUpper};

    } else {
      bool accept;

      double bound;
      if (mipsolver->variableType(Rindex[i]) != HighsVarType::kContinuous) {
        bound = std::ceil(double(boundVal - mipsolver->mipdata_->feastol));
        if (bound > col_lower_[Rindex[i]] &&
            bound - col_lower_[Rindex[i]] >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(col_upper_[Rindex[i]] - double(boundVal)) <=
            mipsolver->mipdata_->epsilon)
          bound = col_upper_[Rindex[i]];
        else
          bound = double(boundVal);
        if (col_lower_[Rindex[i]] == -kHighsInf)
          accept = true;
        else if (bound - 1000.0 * mipsolver->mipdata_->feastol >
                 col_lower_[Rindex[i]]) {
          double relativeImprove = bound - col_lower_[Rindex[i]];
          if (col_upper_[Rindex[i]] != kHighsInf)
            relativeImprove /= col_upper_[Rindex[i]] - col_lower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(col_lower_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kLower};
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
        Rvalue[i], col_lower_[Rindex[i]], col_upper_[Rindex[i]]);
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
        if (bound < col_upper_[Rindex[i]] &&
            col_upper_[Rindex[i]] - bound >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(double(boundVal) - col_lower_[Rindex[i]]) <=
            mipsolver->mipdata_->epsilon)
          bound = col_lower_[Rindex[i]];
        else
          bound = double(boundVal);
        if (col_upper_[Rindex[i]] == kHighsInf)
          accept = true;
        else if (bound + 1000.0 * mipsolver->mipdata_->feastol <
                 col_upper_[Rindex[i]]) {
          double relativeImprove = col_upper_[Rindex[i]] - bound;
          if (col_lower_[Rindex[i]] != -kHighsInf)
            relativeImprove /= col_upper_[Rindex[i]] - col_lower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(col_upper_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kUpper};
    } else {
      bool accept;

      double bound;
      if (mipsolver->variableType(Rindex[i]) != HighsVarType::kContinuous) {
        bound = std::ceil(double(boundVal - mipsolver->mipdata_->feastol));
        if (bound > col_lower_[Rindex[i]] &&
            bound - col_lower_[Rindex[i]] >
                1000.0 * mipsolver->mipdata_->feastol * std::abs(bound))
          accept = true;
        else
          accept = false;
      } else {
        if (std::abs(col_upper_[Rindex[i]] - double(boundVal)) <=
            mipsolver->mipdata_->epsilon)
          bound = col_upper_[Rindex[i]];
        else
          bound = double(boundVal);
        if (col_lower_[Rindex[i]] == -kHighsInf)
          accept = true;
        else if (bound - 1000.0 * mipsolver->mipdata_->feastol >
                 col_lower_[Rindex[i]]) {
          double relativeImprove = bound - col_lower_[Rindex[i]];
          if (col_upper_[Rindex[i]] != kHighsInf)
            relativeImprove /= col_upper_[Rindex[i]] - col_lower_[Rindex[i]];
          else
            relativeImprove /=
                std::max(std::abs(col_lower_[Rindex[i]]), std::abs(bound));
          accept = relativeImprove >= 0.3;
        } else
          accept = false;
      }

      if (accept)
        boundchgs[numchgs++] = {bound, Rindex[i], HighsBoundType::kLower};
    }
  }

  return numchgs;
}

void HighsDomain::updateThresholdLbChange(HighsInt col, double newbound,
                                          double val, double& threshold) {
  if (newbound != col_upper_[col]) {
    double boundRange = (col_upper_[col] - newbound);

    boundRange -=
        variableType(col) == HighsVarType::kContinuous
            ? std::max(0.3 * boundRange, 1000.0 * mipsolver->mipdata_->feastol)
            : mipsolver->mipdata_->feastol;

    double thresholdNew = std::abs(val) * boundRange;

    // the new threshold is now the maximum of the new threshold and the current
    // one
    threshold =
        std::max({threshold, thresholdNew, mipsolver->mipdata_->feastol});
  }
}

void HighsDomain::updateThresholdUbChange(HighsInt col, double newbound,
                                          double val, double& threshold) {
  if (newbound != col_lower_[col]) {
    double boundRange = (newbound - col_lower_[col]);

    boundRange -=
        variableType(col) == HighsVarType::kContinuous
            ? std::max(0.3 * boundRange, 1000.0 * mipsolver->mipdata_->feastol)
            : mipsolver->mipdata_->feastol;

    double thresholdNew = std::abs(val) * boundRange;

    // the new threshold is now the maximum of the new threshold and the current
    // one
    threshold =
        std::max({threshold, thresholdNew, mipsolver->mipdata_->feastol});
  }
}

void HighsDomain::updateActivityLbChange(HighsInt col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  HighsInt start = mip->a_start_[col];
  HighsInt end = mip->a_start_[col + 1];

  assert(!infeasible_);

  for (HighsInt i = start; i != end; ++i) {
    if (mip->a_value_[i] > 0) {
      double deltamin;
      if (oldbound == -kHighsInf) {
        --activitymininf_[mip->a_index_[i]];
        deltamin = newbound * mip->a_value_[i];
      } else if (newbound == -kHighsInf) {
        ++activitymininf_[mip->a_index_[i]];
        deltamin = -oldbound * mip->a_value_[i];
      } else {
        deltamin = (newbound - oldbound) * mip->a_value_[i];
      }
      activitymin_[mip->a_index_[i]] += deltamin;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpminact;
        computeMinActivity(mipsolver->mipdata_->ARstart_[mip->a_index_[i]],
                           mipsolver->mipdata_->ARstart_[mip->a_index_[i] + 1],
                           mipsolver->mipdata_->ARindex_.data(),
                           mipsolver->mipdata_->ARvalue_.data(), tmpinf,
                           tmpminact);
        assert(std::abs(double(activitymin_[mip->a_index_[i]] - tmpminact)) <=
               mipsolver->mipdata_->feastol);
        assert(tmpinf == activitymininf_[mip->a_index_[i]]);
      }
#endif

      if (deltamin <= 0) {
        updateThresholdLbChange(col, newbound, mip->a_value_[i],
                                capacityThreshold_[mip->a_index_[i]]);
        continue;
      }

      if (mip->row_upper_[mip->a_index_[i]] != kHighsInf &&
          activitymininf_[mip->a_index_[i]] == 0 &&
          activitymin_[mip->a_index_[i]] - mip->row_upper_[mip->a_index_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowUpper(mip->a_index_[i]);
        end = i + 1;
        break;
      }

      if (activitymininf_[mip->a_index_[i]] <= 1 &&
          !propagateflags_[mip->a_index_[i]] &&
          mip->row_upper_[mip->a_index_[i]] != kHighsInf)
        markPropagate(mip->a_index_[i]);
    } else {
      double deltamax;
      if (oldbound == -kHighsInf) {
        --activitymaxinf_[mip->a_index_[i]];
        deltamax = newbound * mip->a_value_[i];
      } else if (newbound == -kHighsInf) {
        ++activitymaxinf_[mip->a_index_[i]];
        deltamax = -oldbound * mip->a_value_[i];
      } else {
        deltamax = (newbound - oldbound) * mip->a_value_[i];
      }
      activitymax_[mip->a_index_[i]] += deltamax;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpmaxact;
        computeMaxActivity(mipsolver->mipdata_->ARstart_[mip->a_index_[i]],
                           mipsolver->mipdata_->ARstart_[mip->a_index_[i] + 1],
                           mipsolver->mipdata_->ARindex_.data(),
                           mipsolver->mipdata_->ARvalue_.data(), tmpinf,
                           tmpmaxact);
        assert(std::abs(double(activitymax_[mip->a_index_[i]] - tmpmaxact)) <=
               mipsolver->mipdata_->feastol);
        assert(tmpinf == activitymaxinf_[mip->a_index_[i]]);
      }
#endif

      if (deltamax >= 0) continue;

      if (mip->row_lower_[mip->a_index_[i]] != -kHighsInf &&
          activitymaxinf_[mip->a_index_[i]] == 0 &&
          mip->row_lower_[mip->a_index_[i]] - activitymax_[mip->a_index_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowLower(mip->a_index_[i]);
        end = i + 1;
        break;
      }

      if (activitymaxinf_[mip->a_index_[i]] <= 1 &&
          !propagateflags_[mip->a_index_[i]] &&
          mip->row_lower_[mip->a_index_[i]] != -kHighsInf)
        markPropagate(mip->a_index_[i]);
    }
  }

  if (!infeasible_) {
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.updateActivityLbChange(col, oldbound, newbound);
  } else {
    assert(infeasible_reason.type == Reason::kModelRowLower ||
           infeasible_reason.type == Reason::kModelRowUpper);
    assert(infeasible_reason.index == mip->a_index_[end - 1]);
  }

  if (infeasible_) {
    std::swap(oldbound, newbound);
    for (HighsInt i = start; i != end; ++i) {
      if (mip->a_value_[i] > 0) {
        double deltamin;
        if (oldbound == -kHighsInf) {
          --activitymininf_[mip->a_index_[i]];
          deltamin = newbound * mip->a_value_[i];
        } else if (newbound == -kHighsInf) {
          ++activitymininf_[mip->a_index_[i]];
          deltamin = -oldbound * mip->a_value_[i];
        } else {
          deltamin = (newbound - oldbound) * mip->a_value_[i];
        }
        activitymin_[mip->a_index_[i]] += deltamin;
      } else {
        double deltamax;
        if (oldbound == -kHighsInf) {
          --activitymaxinf_[mip->a_index_[i]];
          deltamax = newbound * mip->a_value_[i];
        } else if (newbound == -kHighsInf) {
          ++activitymaxinf_[mip->a_index_[i]];
          deltamax = -oldbound * mip->a_value_[i];
        } else {
          deltamax = (newbound - oldbound) * mip->a_value_[i];
        }
        activitymax_[mip->a_index_[i]] += deltamax;
      }
    }

    return;
  } else {
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.updateActivityLbChange(col, oldbound, newbound);
  }
}

void HighsDomain::updateActivityUbChange(HighsInt col, double oldbound,
                                         double newbound) {
  auto mip = mipsolver->model_;
  HighsInt start = mip->a_start_[col];
  HighsInt end = mip->a_start_[col + 1];

  assert(!infeasible_);

  for (HighsInt i = start; i != end; ++i) {
    if (mip->a_value_[i] > 0) {
      double deltamax;
      if (oldbound == kHighsInf) {
        --activitymaxinf_[mip->a_index_[i]];
        deltamax = newbound * mip->a_value_[i];
      } else if (newbound == kHighsInf) {
        ++activitymaxinf_[mip->a_index_[i]];
        deltamax = -oldbound * mip->a_value_[i];
      } else {
        deltamax = (newbound - oldbound) * mip->a_value_[i];
      }
      activitymax_[mip->a_index_[i]] += deltamax;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpmaxact;
        computeMaxActivity(mipsolver->mipdata_->ARstart_[mip->a_index_[i]],
                           mipsolver->mipdata_->ARstart_[mip->a_index_[i] + 1],
                           mipsolver->mipdata_->ARindex_.data(),
                           mipsolver->mipdata_->ARvalue_.data(), tmpinf,
                           tmpmaxact);
        assert(std::abs(double(activitymax_[mip->a_index_[i]] - tmpmaxact)) <=
               mipsolver->mipdata_->feastol);
        assert(tmpinf == activitymaxinf_[mip->a_index_[i]]);
      }
#endif

      if (deltamax >= 0) {
        updateThresholdUbChange(col, newbound, mip->a_value_[i],
                                capacityThreshold_[mip->a_index_[i]]);
        continue;
      }

      if (mip->row_lower_[mip->a_index_[i]] != -kHighsInf &&
          activitymaxinf_[mip->a_index_[i]] == 0 &&
          mip->row_lower_[mip->a_index_[i]] - activitymax_[mip->a_index_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowLower(mip->a_index_[i]);
        end = i + 1;
        break;
      }

      if (activitymaxinf_[mip->a_index_[i]] <= 1 &&
          !propagateflags_[mip->a_index_[i]] &&
          mip->row_lower_[mip->a_index_[i]] != -kHighsInf) {
        markPropagate(mip->a_index_[i]);
        // propagateflags_[mip->a_index_[i]] = 1;
        // propagateinds_.push_back(mip->a_index_[i]);
      }
    } else {
      double deltamin;
      if (oldbound == kHighsInf) {
        --activitymininf_[mip->a_index_[i]];
        deltamin = newbound * mip->a_value_[i];
      } else if (newbound == kHighsInf) {
        ++activitymininf_[mip->a_index_[i]];
        deltamin = -oldbound * mip->a_value_[i];
      } else {
        deltamin = (newbound - oldbound) * mip->a_value_[i];
      }

      activitymin_[mip->a_index_[i]] += deltamin;

#ifndef NDEBUG
      {
        HighsInt tmpinf;
        HighsCDouble tmpminact;
        computeMinActivity(mipsolver->mipdata_->ARstart_[mip->a_index_[i]],
                           mipsolver->mipdata_->ARstart_[mip->a_index_[i] + 1],
                           mipsolver->mipdata_->ARindex_.data(),
                           mipsolver->mipdata_->ARvalue_.data(), tmpinf,
                           tmpminact);
        assert(std::abs(double(activitymin_[mip->a_index_[i]] - tmpminact)) <=
               mipsolver->mipdata_->feastol);
        assert(tmpinf == activitymininf_[mip->a_index_[i]]);
      }
#endif

      if (deltamin <= 0) continue;

      if (mip->row_upper_[mip->a_index_[i]] != kHighsInf &&
          activitymininf_[mip->a_index_[i]] == 0 &&
          activitymin_[mip->a_index_[i]] - mip->row_upper_[mip->a_index_[i]] >
              mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        infeasible_ = true;
        infeasible_pos = domchgstack_.size();
        infeasible_reason = Reason::modelRowUpper(mip->a_index_[i]);
        end = i + 1;
        break;
      }

      if (activitymininf_[mip->a_index_[i]] <= 1 &&
          !propagateflags_[mip->a_index_[i]] &&
          mip->row_upper_[mip->a_index_[i]] != kHighsInf) {
        markPropagate(mip->a_index_[i]);
        // propagateflags_[mip->a_index_[i]] = 1;
        // propagateinds_.push_back(mip->a_index_[i]);
      }
    }
  }

  if (!infeasible_) {
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.updateActivityUbChange(col, oldbound, newbound);
  } else {
    assert(infeasible_reason.type == Reason::kModelRowLower ||
           infeasible_reason.type == Reason::kModelRowUpper);
    assert(infeasible_reason.index == mip->a_index_[end - 1]);
  }

  if (infeasible_) {
    std::swap(oldbound, newbound);
    for (HighsInt i = start; i != end; ++i) {
      if (mip->a_value_[i] > 0) {
        double deltamax;
        if (oldbound == kHighsInf) {
          --activitymaxinf_[mip->a_index_[i]];
          deltamax = newbound * mip->a_value_[i];
        } else if (newbound == kHighsInf) {
          ++activitymaxinf_[mip->a_index_[i]];
          deltamax = -oldbound * mip->a_value_[i];
        } else {
          deltamax = (newbound - oldbound) * mip->a_value_[i];
        }
        activitymax_[mip->a_index_[i]] += deltamax;
      } else {
        double deltamin;
        if (oldbound == kHighsInf) {
          --activitymininf_[mip->a_index_[i]];
          deltamin = newbound * mip->a_value_[i];
        } else if (newbound == kHighsInf) {
          ++activitymininf_[mip->a_index_[i]];
          deltamin = -oldbound * mip->a_value_[i];
        } else {
          deltamin = (newbound - oldbound) * mip->a_value_[i];
        }

        activitymin_[mip->a_index_[i]] += deltamin;
      }
    }

    return;
  } else {
    for (ConflictPoolPropagation& conflictprop : conflictPoolPropagation)
      conflictprop.updateActivityUbChange(col, oldbound, newbound);
  }
}

void HighsDomain::recomputeCapacityThreshold(HighsInt row) {
  HighsInt start = mipsolver->mipdata_->ARstart_[row];
  HighsInt end = mipsolver->mipdata_->ARstart_[row + 1];

  capacityThreshold_[row] = 0.0;
  for (HighsInt i = start; i < end; ++i) {
    HighsInt col = mipsolver->mipdata_->ARindex_[i];

    if (col_upper_[col] == col_lower_[col]) continue;

    double boundRange = col_upper_[col] - col_lower_[col];

    boundRange -=
        variableType(col) == HighsVarType::kContinuous
            ? std::max(0.3 * boundRange, 1000.0 * mipsolver->mipdata_->feastol)
            : mipsolver->mipdata_->feastol;

    double threshold = std::abs(mipsolver->mipdata_->ARvalue_[i]) * boundRange;

    capacityThreshold_[row] = std::max(
        {capacityThreshold_[row], threshold, mipsolver->mipdata_->feastol});
  }
}

void HighsDomain::markPropagateCut(Reason reason) {
  switch (reason.type) {
    case Reason::kUnknown:
    case Reason::kCliqueTable:
    case Reason::kBranching:
    case Reason::kModelRowLower:
    case Reason::kModelRowUpper:
    case Reason::kConflictingBounds:
      break;
    default:
      assert(reason.type >= 0 &&
             reason.type < HighsInt(cutpoolpropagation.size() +
                                    conflictPoolPropagation.size()));
      if (reason.type < (HighsInt)cutpoolpropagation.size())
        cutpoolpropagation[reason.type].markPropagateCut(reason.index);
      else
        conflictPoolPropagation[reason.type - cutpoolpropagation.size()]
            .markPropagateConflict(reason.index);
  }
}

void HighsDomain::markPropagate(HighsInt row) {
  if (!propagateflags_[row]) {
    bool proplower = mipsolver->rowLower(row) != -kHighsInf &&
                     (activitymaxinf_[row] == 1 ||
                      (double(activitymax_[row]) - mipsolver->rowLower(row)) <=
                          capacityThreshold_[row]);
    bool propupper = mipsolver->rowUpper(row) != kHighsInf &&
                     (activitymininf_[row] == 1 ||
                      (mipsolver->rowUpper(row) - double(activitymin_[row])) <=
                          capacityThreshold_[row]);

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
  capacityThreshold_.resize(mipsolver->numRow());
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

    recomputeCapacityThreshold(i);

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
    oldbound = col_lower_[boundchg.column];
    col_lower_[boundchg.column] = boundchg.boundval;
    if (oldbound != boundchg.boundval) {
      if (!infeasible_)
        updateActivityLbChange(boundchg.column, oldbound, boundchg.boundval);

      if (!changedcolsflags_[boundchg.column]) {
        changedcolsflags_[boundchg.column] = 1;
        changedcols_.push_back(boundchg.column);
      }
    }
  } else {
    oldbound = col_upper_[boundchg.column];
    col_upper_[boundchg.column] = boundchg.boundval;
    if (oldbound != boundchg.boundval) {
      if (!infeasible_)
        updateActivityUbChange(boundchg.column, oldbound, boundchg.boundval);

      if (!changedcolsflags_[boundchg.column]) {
        changedcolsflags_[boundchg.column] = 1;
        changedcols_.push_back(boundchg.column);
      }
    }
  }

  return oldbound;
}

void HighsDomain::changeBound(HighsDomainChange boundchg, Reason reason) {
  assert(boundchg.column >= 0);
  assert(boundchg.column < (HighsInt)col_upper_.size());
  // assert(infeasible_ == 0);
  mipsolver->mipdata_->debugSolution.boundChangeAdded(
      *this, boundchg, reason.type == Reason::kBranching);
  HighsInt prevPos;
  if (boundchg.boundtype == HighsBoundType::kLower) {
    if (boundchg.boundval <= col_lower_[boundchg.column]) {
      if (reason.type != Reason::kBranching) return;
      boundchg.boundval = col_lower_[boundchg.column];
    }
    if (boundchg.boundval > col_upper_[boundchg.column]) {
      if (boundchg.boundval - col_upper_[boundchg.column] >
          mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        if (!infeasible_) {
          infeasible_pos = domchgstack_.size();
          infeasible_ = true;
          infeasible_reason = Reason::conflictingBounds(domchgstack_.size());
        }
      } else {
        boundchg.boundval = col_upper_[boundchg.column];
        if (boundchg.boundval == col_lower_[boundchg.column]) return;
      }
    }

    prevPos = colLowerPos_[boundchg.column];
    colLowerPos_[boundchg.column] = domchgstack_.size();
  } else {
    if (boundchg.boundval >= col_upper_[boundchg.column]) {
      if (reason.type != Reason::kBranching) return;
      boundchg.boundval = col_upper_[boundchg.column];
    }
    if (boundchg.boundval < col_lower_[boundchg.column]) {
      if (col_lower_[boundchg.column] - boundchg.boundval >
          mipsolver->mipdata_->feastol) {
        mipsolver->mipdata_->debugSolution.nodePruned(*this);
        if (!infeasible_) {
          infeasible_pos = domchgstack_.size();
          infeasible_ = true;
          infeasible_reason = Reason::conflictingBounds(domchgstack_.size());
        }
      } else {
        boundchg.boundval = col_lower_[boundchg.column];
        if (boundchg.boundval == col_upper_[boundchg.column]) return;
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

  if (binary && !infeasible_ && isFixed(boundchg.column))
    mipsolver->mipdata_->cliquetable.addImplications(
        *this, boundchg.column, col_lower_[boundchg.column] > 0.5);
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
        domchgstack[k].boundval <= col_lower_[domchgstack[k].column])
      continue;
    if (domchgstack[k].boundtype == HighsBoundType::kUpper &&
        domchgstack[k].boundval >= col_upper_[domchgstack[k].column])
      continue;

    changeBound(domchgstack[k], Reason::unspecified());

    if (infeasible_) break;
  }
}

void HighsDomain::setDomainChangeStack(
    const std::vector<HighsDomainChange>& domchgstack,
    const std::vector<HighsInt>& branchingPositions) {
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
  HighsInt nextBranchPos = -1;
  HighsInt k = 0;
  for (HighsInt branchPos : branchingPositions) {
    for (; k < branchPos; ++k) {
      if (domchgstack[k].boundtype == HighsBoundType::kLower &&
          domchgstack[k].boundval <= col_lower_[domchgstack[k].column])
        continue;
      if (domchgstack[k].boundtype == HighsBoundType::kUpper &&
          domchgstack[k].boundval >= col_upper_[domchgstack[k].column])
        continue;

      changeBound(domchgstack[k], Reason::unspecified());
      if (!infeasible_) propagate();
      if (infeasible_) return;
    }

    if (k == stacksize) return;

    // do not skip redundant branching changes, as we need to keep their status
    // as branching variables for computing correct stabilizers
    changeBound(domchgstack[k], Reason::branching());
    if (!infeasible_) propagate();
    if (infeasible_) return;
  }

  for (; k < stacksize; ++k) {
    if (domchgstack[k].boundtype == HighsBoundType::kLower &&
        domchgstack[k].boundval <= col_lower_[domchgstack[k].column])
      continue;
    if (domchgstack[k].boundtype == HighsBoundType::kUpper &&
        domchgstack[k].boundval >= col_upper_[domchgstack[k].column])
      continue;

    mipsolver->mipdata_->debugSolution.boundChangeAdded(*this, domchgstack[k],
                                                        true);

    changeBound(domchgstack[k], Reason::unspecified());
    if (!infeasible_) propagate();
    if (infeasible_) break;
  }
}

void HighsDomain::backtrackToGlobal() {
  HighsInt k = HighsInt(domchgstack_.size()) - 1;
  bool old_infeasible = infeasible_;
  Reason old_reason = infeasible_reason;

  if (infeasible_ && infeasible_pos == HighsInt(domchgstack_.size())) {
    assert(old_infeasible);
    assert(k == HighsInt(domchgstack_.size()) - 1);
    infeasible_ = false;
    infeasible_reason = Reason::unspecified();
  }

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

    if (prevbound != domchgstack_[k].boundval) {
      // change back to global bound
      doChangeBound(
          {prevbound, domchgstack_[k].column, domchgstack_[k].boundtype});
    }

    if (infeasible_ && infeasible_pos == k) {
      assert(old_infeasible);
      assert(k == HighsInt(domchgstack_.size()) - 1);
      infeasible_ = false;
      infeasible_reason = Reason::unspecified();
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

  domchgstack_.clear();
  prevboundval_.clear();
  domchgreason_.clear();
  branchPos_.clear();
}

HighsDomainChange HighsDomain::backtrack() {
  HighsInt k = HighsInt(domchgstack_.size()) - 1;
  bool old_infeasible = infeasible_;
  Reason old_reason = infeasible_reason;

  if (infeasible_ && infeasible_pos == HighsInt(domchgstack_.size())) {
    assert(old_infeasible);
    assert(k == HighsInt(domchgstack_.size()) - 1);
    infeasible_ = false;
    infeasible_reason = Reason::unspecified();
  }
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
        {prevbound, domchgstack_[k].column, domchgstack_[k].boundtype});

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
    return HighsDomainChange({0.0, -1, HighsBoundType::kLower});
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

      for (const auto& conflictprop : conflictPoolPropagation) {
        if (!conflictprop.propagateConflictInds_.empty()) {
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
    const HighsInt numConflictPools = conflictPoolPropagation.size();
    for (HighsInt conflictPool = 0; conflictPool < numConflictPools;
         ++conflictPool) {
      auto& conflictprop = conflictPoolPropagation[conflictPool];
      while (!conflictprop.propagateConflictInds_.empty()) {
        propagateinds.swap(conflictprop.propagateConflictInds_);

        for (HighsInt conflict : propagateinds)
          conflictprop.propagateConflict(conflict);

        propagateinds.clear();
      }
    }
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
          HighsInt i = propagateinds[k];

          if (propRowNumChangedBounds_[k].first != 0) {
            HighsInt start = 2 * mipsolver->mipdata_->ARstart_[i];
            HighsInt end = start + propRowNumChangedBounds_[k].first;
            for (HighsInt j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::modelRowUpper(i));

            if (infeasible_) break;
          }
          if (propRowNumChangedBounds_[k].second != 0) {
            HighsInt start = 2 * mipsolver->mipdata_->ARstart_[i] +
                             propRowNumChangedBounds_[k].first;
            HighsInt end = start + propRowNumChangedBounds_[k].second;
            for (HighsInt j = start; j != end && !infeasible_; ++j)
              changeBound(changedbounds[j], Reason::modelRowLower(i));

            if (infeasible_) break;
          }

          recomputeCapacityThreshold(i);
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
          cutpoolprop.propagatecutflags_[cut] &= 2;
          propnnz += cutpoolprop.cutpool->getMatrix().getRowEnd(cut) -
                     cutpoolprop.cutpool->getMatrix().getRowStart(cut);
        }

        if (!infeasible_) {
          propRowNumChangedBounds_.assign(
              numproprows, std::make_pair(HighsInt{0}, HighsInt{0}));

          auto propagateIndex = [&](HighsInt k) {
            // first check if cut is marked as deleted
            if (cutpoolprop.propagatecutflags_[k] & 2) return;
            HighsInt i = propagateinds[k];

            HighsInt Rlen;
            const HighsInt* Rindex;
            const double* Rvalue;
            cutpoolprop.cutpool->getCut(i, Rlen, Rindex, Rvalue);
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
            HighsInt i = propagateinds[k];
            if (propRowNumChangedBounds_[k].first != 0) {
              cutpoolprop.cutpool->resetAge(i);
              HighsInt start = cutpoolprop.cutpool->getMatrix().getRowStart(i);
              HighsInt end = start + propRowNumChangedBounds_[k].first;
              for (HighsInt j = start; j != end && !infeasible_; ++j)
                changeBound(changedbounds[j], Reason::cut(cutpool, i));
            }

            cutpoolprop.recomputeCapacityThreshold(i);

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
  double lb = col_lower_[col];
  pos = colLowerPos_[col];
  while (pos > stackpos || (pos != -1 && prevboundval_[pos].first == lb)) {
    lb = prevboundval_[pos].first;
    pos = prevboundval_[pos].second;
  }
  return lb;
}

double HighsDomain::getColUpperPos(HighsInt col, HighsInt stackpos,
                                   HighsInt& pos) const {
  double ub = col_upper_[col];
  pos = colUpperPos_[col];
  while (pos > stackpos || (pos != -1 && prevboundval_[pos].first == ub)) {
    ub = prevboundval_[pos].first;
    pos = prevboundval_[pos].second;
  }
  return ub;
}

void HighsDomain::conflictAnalysis(HighsConflictPool& conflictPool) {
  if (&mipsolver->mipdata_->domain == this) return;
  if (mipsolver->mipdata_->domain.infeasible() || !infeasible_) return;

  ConflictSet conflictSet(*this);
  // if (!conflictSet.resolvable(infeasible_pos)) return;

  conflictSet.conflictAnalysis(conflictPool);
}

void HighsDomain::conflictAnalysis(const HighsInt* proofinds,
                                   const double* proofvals, HighsInt prooflen,
                                   double proofrhs,
                                   HighsConflictPool& conflictPool) {
  if (&mipsolver->mipdata_->domain == this) return;
  if (mipsolver->mipdata_->domain.infeasible()) return;

  ConflictSet conflictSet(*this);
  conflictSet.conflictAnalysis(proofinds, proofvals, prooflen, proofrhs,
                               conflictPool);
}

void HighsDomain::conflictAnalyzeReconvergence(
    const HighsDomainChange& domchg, const HighsInt* proofinds,
    const double* proofvals, HighsInt prooflen, double proofrhs,
    HighsConflictPool& conflictPool) {
  if (&mipsolver->mipdata_->domain == this) return;
  ConflictSet conflictSet(*this);

  HighsInt ninfmin;
  HighsCDouble activitymin;
  mipsolver->mipdata_->domain.computeMinActivity(
      0, prooflen, proofinds, proofvals, ninfmin, activitymin);
  if (ninfmin != 0) return;

  if (!conflictSet.explainBoundChangeLeq(domchg, domchgstack_.size(), proofinds,
                                         proofvals, prooflen, proofrhs,
                                         double(activitymin)))
    return;

  if (conflictSet.resolvedDomainChanges.size() >
      0.3 * mipsolver->mipdata_->integral_cols.size())
    return;

  conflictSet.reconvergenceFrontier.insert(
      conflictSet.resolvedDomainChanges.begin(),
      conflictSet.resolvedDomainChanges.end());
  conflictSet.resolveDepth(conflictSet.reconvergenceFrontier, branchPos_.size(),
                           0);
  conflictPool.addReconvergenceCut(*this, conflictSet.reconvergenceFrontier,
                                   domchg);
}

void HighsDomain::tightenCoefficients(HighsInt* inds, double* vals,
                                      HighsInt len, double& rhs) const {
  HighsCDouble maxactivity = 0;

  for (HighsInt i = 0; i != len; ++i) {
    if (vals[i] > 0) {
      if (col_upper_[inds[i]] == kHighsInf) return;

      maxactivity += col_upper_[inds[i]] * vals[i];
    } else {
      if (col_lower_[inds[i]] == -kHighsInf) return;

      maxactivity += col_lower_[inds[i]] * vals[i];
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
        upper -= delta * col_upper_[inds[i]];
        vals[i] = double(maxabscoef);
        ++tightened;
      } else if (vals[i] < -maxabscoef) {
        HighsCDouble delta = -vals[i] - maxabscoef;
        upper += delta * col_lower_[inds[i]];
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
                                      HighsInt cut) const {
  for (auto& cutpoolprop : cutpoolpropagation) {
    if (cutpoolprop.cutpool == &cutpool) {
      // assert((cutpoolprop.propagatecutflags_[cut] & 2) == 0);

      return cut < (HighsInt)cutpoolprop.propagatecutflags_.size() &&
                     (cutpoolprop.propagatecutflags_[cut] & 2) == 0 &&
                     cutpoolprop.activitycutsinf_[cut] == 0
                 ? double(cutpoolprop.activitycuts_[cut])
                 : -kHighsInf;
    }
  }

  return -kHighsInf;
}

bool HighsDomain::isFixing(const HighsDomainChange& domchg) const {
  double otherbound = domchg.boundtype == HighsBoundType::kUpper
                          ? col_lower_[domchg.column]
                          : col_upper_[domchg.column];
  return std::abs(domchg.boundval - otherbound) <= mipsolver->mipdata_->epsilon;
}

HighsDomainChange HighsDomain::flip(const HighsDomainChange& domchg) const {
  if (domchg.boundtype == HighsBoundType::kLower) {
    HighsDomainChange flipped{domchg.boundval - mipsolver->mipdata_->feastol,
                              domchg.column, HighsBoundType::kUpper};
    if (mipsolver->variableType(domchg.column) != HighsVarType::kContinuous)
      flipped.boundval = std::floor(flipped.boundval);
    return flipped;
  } else {
    HighsDomainChange flipped{domchg.boundval + mipsolver->mipdata_->feastol,
                              domchg.column, HighsBoundType::kLower};
    if (mipsolver->variableType(domchg.column) != HighsVarType::kContinuous)
      flipped.boundval = std::ceil(flipped.boundval);
    return flipped;
  }
}

double HighsDomain::feastol() const { return mipsolver->mipdata_->feastol; }

HighsDomain::ConflictSet::ConflictSet(HighsDomain& localdom_)
    : localdom(localdom_),
      globaldom(localdom.mipsolver->mipdata_->domain),
      reasonSideFrontier(),
      reconvergenceFrontier(),
      resolveQueue(),
      resolvedDomainChanges() {}

bool HighsDomain::ConflictSet::explainBoundChangeGeq(
    const HighsDomainChange& domchg, HighsInt pos, const HighsInt* inds,
    const double* vals, HighsInt len, double rhs, double maxAct) {
  if (maxAct == kHighsInf) return false;

  // get the coefficient value of the column for which we want to explain
  // the bound change
  double domchgVal = 0;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    if (col == domchg.column) {
      domchgVal = vals[i];
      continue;
    }

    double delta;
    HighsInt numNodes;
    HighsInt boundpos;
    if (vals[i] > 0) {
      double ub = localdom.getColUpperPos(col, pos, boundpos);
      if (globaldom.col_upper_[col] <= ub) continue;
      delta = vals[i] * (ub - globaldom.col_upper_[col]);
      numNodes = nodequeue.numNodesDown(col);
    } else {
      double lb = localdom.getColLowerPos(col, pos, boundpos);
      if (globaldom.col_lower_[col] >= lb) continue;
      delta = vals[i] * (lb - globaldom.col_lower_[col]);
      numNodes = nodequeue.numNodesUp(col);
    }

    if (boundpos == -1) continue;

    resolveBuffer.emplace_back(delta, numNodes, boundpos);
  }

  if (domchgVal == 0) return false;

  std::sort(resolveBuffer.begin(), resolveBuffer.end(),
            [&](const std::tuple<double, HighsInt, HighsInt>& a,
                const std::tuple<double, HighsInt, HighsInt>& b) {
              double prioA = std::get<0>(a) * (std::get<1>(a) + 1);
              double prioB = std::get<0>(b) * (std::get<1>(b) + 1);
              return prioA > prioB;
            });

  // to explain the bound change we start from the bound constraint,
  // multiply it by the columns coefficient in the constraint. Then the
  // bound constraint is a0 * x0 >= b0 and the constraint
  // a0 * x0 + \sum ai * xi >= b. Let the max activity of \sum ai xi be M,
  // then the constraint yields the bound a0 * x0 >= b - M. If M is
  // sufficiently small the bound constraint is implied. Therefore we
  // decrease M by updating it with the stronger local bounds until
  // M <= b - b0 holds.
  double b0 = domchg.boundval;
  if (localdom.mipsolver->variableType(domchg.column) !=
      HighsVarType::kContinuous) {
    // in case of an integral variable the bound was rounded and can be
    // relaxed by 1-feastol. We use 1 - 10 * feastol for numerical safety.
    if (domchg.boundtype == HighsBoundType::kLower)
      b0 -= (1.0 - 10 * localdom.mipsolver->mipdata_->feastol);
    else
      b0 += (1.0 - 10 * localdom.mipsolver->mipdata_->feastol);
  } else {
    // for a continuous variable we relax the bound by epsilon to
    // accomodate for tiny rounding errors
    if (domchg.boundtype == HighsBoundType::kLower)
      b0 -= localdom.mipsolver->mipdata_->epsilon;
    else
      b0 += localdom.mipsolver->mipdata_->epsilon;
  }

  // now multiply the bound constraint with the coefficient value in the
  // constraint to obtain b0
  b0 *= domchgVal;

  // compute the lower bound of M that is necessary
  double Mupper = rhs - b0;

  // M is the global residual activity initially
  double M = maxAct;
  if (domchgVal < 0)
    M -= domchgVal * globaldom.col_lower_[domchg.column];
  else
    M -= domchgVal * globaldom.col_upper_[domchg.column];

  resolvedDomainChanges.clear();
  for (const std::tuple<double, HighsInt, HighsInt>& reasonDomchg :
       resolveBuffer) {
    M += std::get<0>(reasonDomchg);
    resolvedDomainChanges.push_back(std::get<2>(reasonDomchg));
    assert(resolvedDomainChanges.back() >= 0);
    assert(resolvedDomainChanges.back() < localdom.domchgstack_.size());
    if (M < Mupper) break;
  }

  if (M >= Mupper) {
    // printf("local bounds reach only value of %.12g, need at most
    // %.12g\n",
    //        M, Mupper);
    return false;
  }

  return true;
}

bool HighsDomain::ConflictSet::explainInfeasibility() {
  switch (localdom.infeasible_reason.type) {
    case Reason::kUnknown:
    case Reason::kBranching:
      return false;
    case Reason::kConflictingBounds: {
      resolvedDomainChanges.clear();
      HighsInt conflictingBoundPos = localdom.infeasible_reason.index;
      HighsInt col = localdom.domchgstack_[conflictingBoundPos].column;

      resolvedDomainChanges.push_back(conflictingBoundPos);

      HighsInt otherBoundPos;
      if (localdom.domchgstack_[conflictingBoundPos].boundtype ==
          HighsBoundType::kLower) {
        double ub =
            localdom.getColUpperPos(col, conflictingBoundPos, otherBoundPos);
        assert(localdom.domchgstack_[conflictingBoundPos].boundval - ub >
               +localdom.mipsolver->mipdata_->feastol);
      } else {
        double lb =
            localdom.getColLowerPos(col, conflictingBoundPos, otherBoundPos);
        assert(localdom.domchgstack_[conflictingBoundPos].boundval - lb <
               -localdom.mipsolver->mipdata_->feastol);
      }
      if (otherBoundPos != -1) resolvedDomainChanges.push_back(otherBoundPos);
      return true;
    }
    case Reason::kCliqueTable:
      assert(false);
      return false;
    case Reason::kModelRowLower: {
      HighsInt rowIndex = localdom.infeasible_reason.index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double maxAct = globaldom.getMaxActivity(rowIndex);

      return explainInfeasibilityGeq(
          inds, vals, len, localdom.mipsolver->rowLower(rowIndex), maxAct);
    }
    case Reason::kModelRowUpper: {
      HighsInt rowIndex = localdom.infeasible_reason.index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double minAct = globaldom.getMinActivity(rowIndex);

      return explainInfeasibilityLeq(
          inds, vals, len, localdom.mipsolver->rowUpper(rowIndex), minAct);
    }
    default:
      assert(localdom.infeasible_reason.type >= 0);
      assert(localdom.infeasible_reason.type <
             HighsInt(localdom.cutpoolpropagation.size() +
                      localdom.conflictPoolPropagation.size()));

      if (localdom.infeasible_reason.type <
          (HighsInt)localdom.cutpoolpropagation.size()) {
        HighsInt cutpoolIndex = localdom.infeasible_reason.type;
        HighsInt cutIndex = localdom.infeasible_reason.index;

        // retrieve the matrix values of the cut
        HighsInt len;
        const HighsInt* inds;
        const double* vals;
        localdom.cutpoolpropagation[cutpoolIndex].cutpool->getCut(cutIndex, len,
                                                                  inds, vals);

        double minAct = globaldom.getMinCutActivity(
            *localdom.cutpoolpropagation[cutpoolIndex].cutpool, cutIndex);

        return explainInfeasibilityLeq(inds, vals, len,
                                       localdom.cutpoolpropagation[cutpoolIndex]
                                           .cutpool->getRhs()[cutIndex],
                                       minAct);
      } else {
        HighsInt conflictPoolIndex = localdom.infeasible_reason.type -
                                     localdom.cutpoolpropagation.size();
        HighsInt conflictIndex = localdom.infeasible_reason.index;

        if (localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictFlag_[conflictIndex] &
            8)
          return false;

        // retrieve the conflict entries
        auto conflictRange =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictRanges()[conflictIndex];
        const HighsDomainChange* conflict =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictEntryVector()
                .data() +
            conflictRange.first;
        HighsInt len = conflictRange.second - conflictRange.first;

        return explainInfeasibilityConflict(conflict, len);
      }
  }

  return false;
}

bool HighsDomain::ConflictSet::explainInfeasibilityConflict(
    const HighsDomainChange* conflict, HighsInt len) {
  resolvedDomainChanges.clear();
  for (HighsInt i = 0; i < len; ++i) {
    if (globaldom.isActive(conflict[i])) continue;

    HighsInt pos;
    if (conflict[i].boundtype == HighsBoundType::kLower) {
      double lb = localdom.getColLowerPos(conflict[i].column,
                                          localdom.infeasible_pos, pos);
      if (pos == -1 || lb < conflict[i].boundval) return false;
    } else {
      double ub = localdom.getColUpperPos(conflict[i].column,
                                          localdom.infeasible_pos, pos);
      if (pos == -1 || ub > conflict[i].boundval) return false;
    }

    resolvedDomainChanges.push_back(pos);
  }

  return true;
}

bool HighsDomain::ConflictSet::explainInfeasibilityLeq(const HighsInt* inds,
                                                       const double* vals,
                                                       HighsInt len, double rhs,
                                                       double minAct) {
  if (minAct == -kHighsInf) return false;

  HighsInt infeasible_pos = kHighsIInf;
  if (localdom.infeasible_) infeasible_pos = localdom.infeasible_pos;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    double delta;
    HighsInt numNodes;
    HighsInt boundpos;
    if (vals[i] > 0) {
      double lb = localdom.getColLowerPos(col, infeasible_pos, boundpos);
      if (globaldom.col_lower_[col] >= lb) continue;
      delta = vals[i] * (lb - globaldom.col_lower_[col]);
      numNodes = nodequeue.numNodesUp(col);
    } else {
      double ub = localdom.getColUpperPos(col, infeasible_pos, boundpos);
      if (globaldom.col_upper_[col] <= ub) continue;
      delta = vals[i] * (ub - globaldom.col_upper_[col]);
      numNodes = nodequeue.numNodesDown(col);
    }

    if (boundpos == -1) continue;

    resolveBuffer.emplace_back(delta, numNodes, boundpos);
  }

  std::sort(resolveBuffer.begin(), resolveBuffer.end(),
            [&](const std::tuple<double, HighsInt, HighsInt>& a,
                const std::tuple<double, HighsInt, HighsInt>& b) {
              double prioA = std::get<0>(a) * (std::get<1>(a) + 1);
              double prioB = std::get<0>(b) * (std::get<1>(b) + 1);
              return prioA > prioB;
            });

  // compute the lower bound of M that is necessary
  double Mlower = rhs + std::max(10.0, std::abs(rhs)) *
                            localdom.mipsolver->mipdata_->feastol;
  // M is the global residual activity initially
  double M = minAct;
  resolvedDomainChanges.clear();
  for (const std::tuple<double, HighsInt, HighsInt>& reasonDomchg :
       resolveBuffer) {
    M += std::get<0>(reasonDomchg);
    resolvedDomainChanges.push_back(std::get<2>(reasonDomchg));
    assert(resolvedDomainChanges.back() >= 0);
    assert(resolvedDomainChanges.back() < localdom.domchgstack_.size());
    if (M > Mlower) break;
  }

  if (M <= Mlower) {
    // printf("local bounds reach only value of %.12g, need at least
    // %.12g\n",
    //        M, Mlower);
    return false;
  }

  return true;
}

bool HighsDomain::ConflictSet::explainInfeasibilityGeq(const HighsInt* inds,
                                                       const double* vals,
                                                       HighsInt len, double rhs,
                                                       double maxAct) {
  if (maxAct == kHighsInf) return false;

  HighsInt infeasible_pos = kHighsIInf;
  if (localdom.infeasible_) infeasible_pos = localdom.infeasible_pos;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    double delta;
    HighsInt numNodes;
    HighsInt boundpos;
    if (vals[i] > 0) {
      double ub = localdom.getColUpperPos(col, infeasible_pos, boundpos);
      if (globaldom.col_upper_[col] <= ub) continue;
      delta = vals[i] * (ub - globaldom.col_upper_[col]);
      numNodes = nodequeue.numNodesDown(col);
    } else {
      double lb = localdom.getColLowerPos(col, infeasible_pos, boundpos);
      if (globaldom.col_lower_[col] >= lb) continue;
      delta = vals[i] * (lb - globaldom.col_lower_[col]);
      numNodes = nodequeue.numNodesUp(col);
    }

    if (boundpos == -1) continue;

    resolveBuffer.emplace_back(delta, numNodes, boundpos);
  }

  std::sort(resolveBuffer.begin(), resolveBuffer.end(),
            [&](const std::tuple<double, HighsInt, HighsInt>& a,
                const std::tuple<double, HighsInt, HighsInt>& b) {
              double prioA = std::get<0>(a) * (std::get<1>(a) + 1);
              double prioB = std::get<0>(b) * (std::get<1>(b) + 1);
              return prioA > prioB;
            });

  // compute the lower bound of M that is necessary
  double Mupper = rhs - std::max(10.0, std::abs(rhs)) *
                            localdom.mipsolver->mipdata_->feastol;

  // M is the global residual activity initially
  double M = maxAct;

  resolvedDomainChanges.clear();
  for (const std::tuple<double, HighsInt, HighsInt>& reasonDomchg :
       resolveBuffer) {
    M += std::get<0>(reasonDomchg);
    resolvedDomainChanges.push_back(std::get<2>(reasonDomchg));
    assert(resolvedDomainChanges.back() >= 0);
    assert(resolvedDomainChanges.back() < localdom.domchgstack_.size());
    if (M < Mupper) break;
  }

  if (M >= Mupper) {
    // printf("local bounds reach only value of %.12g, need at most
    // %.12g\n",
    //        M, Mupper);
    return false;
  }

  return true;
}

bool HighsDomain::ConflictSet::explainBoundChange(HighsInt pos) {
  switch (localdom.domchgreason_[pos].type) {
    case Reason::kUnknown:
    case Reason::kBranching:
    case Reason::kConflictingBounds:
      return false;
    case Reason::kCliqueTable: {
      HighsInt col = localdom.domchgreason_[pos].index >> 1;
      HighsInt val = localdom.domchgreason_[pos].index & 1;
      resolvedDomainChanges.clear();
      HighsInt boundPos;
      if (val) {
        assert(localdom.colLowerPos_[col] >= 0);
        assert(localdom.colLowerPos_[col] < localdom.domchgstack_.size());

        localdom.getColLowerPos(col, pos, boundPos);
      } else {
        assert(localdom.colUpperPos_[col] >= 0);
        assert(localdom.colUpperPos_[col] < localdom.domchgstack_.size());

        localdom.getColUpperPos(col, pos, boundPos);
      }

      if (boundPos != -1) resolvedDomainChanges.push_back(boundPos);

      return true;
    }
    case Reason::kModelRowLower: {
      HighsInt rowIndex = localdom.domchgreason_[pos].index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double maxAct = globaldom.getMaxActivity(rowIndex);

      return explainBoundChangeGeq(localdom.domchgstack_[pos], pos, inds, vals,
                                   len, localdom.mipsolver->rowLower(rowIndex),
                                   maxAct);
    }
    case Reason::kModelRowUpper: {
      HighsInt rowIndex = localdom.domchgreason_[pos].index;

      // retrieve the matrix values of the cut
      HighsInt len;
      const HighsInt* inds;
      const double* vals;
      localdom.mipsolver->mipdata_->getRow(rowIndex, len, inds, vals);

      double minAct = globaldom.getMinActivity(rowIndex);

      return explainBoundChangeLeq(localdom.domchgstack_[pos], pos, inds, vals,
                                   len, localdom.mipsolver->rowUpper(rowIndex),
                                   minAct);
    }
    default:
      assert(localdom.domchgreason_[pos].type >= 0);
      assert(localdom.domchgreason_[pos].type <
             (HighsInt)(localdom.cutpoolpropagation.size() +
                        localdom.conflictPoolPropagation.size()));
      if (localdom.domchgreason_[pos].type <
          (HighsInt)localdom.cutpoolpropagation.size()) {
        HighsInt cutpoolIndex = localdom.domchgreason_[pos].type;
        HighsInt cutIndex = localdom.domchgreason_[pos].index;

        // retrieve the matrix values of the cut
        HighsInt len;
        const HighsInt* inds;
        const double* vals;
        localdom.cutpoolpropagation[cutpoolIndex].cutpool->getCut(cutIndex, len,
                                                                  inds, vals);
        double minAct = globaldom.getMinCutActivity(
            *localdom.cutpoolpropagation[cutpoolIndex].cutpool, cutIndex);

        return explainBoundChangeLeq(localdom.domchgstack_[pos], pos, inds,
                                     vals, len,
                                     localdom.cutpoolpropagation[cutpoolIndex]
                                         .cutpool->getRhs()[cutIndex],
                                     minAct);
      } else {
        HighsInt conflictPoolIndex = localdom.domchgreason_[pos].type -
                                     localdom.cutpoolpropagation.size();
        HighsInt conflictIndex = localdom.domchgreason_[pos].index;

        if (localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictFlag_[conflictIndex] &
            8)
          break;

        // retrieve the conflict entries
        auto conflictRange =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictRanges()[conflictIndex];
        const HighsDomainChange* conflict =
            localdom.conflictPoolPropagation[conflictPoolIndex]
                .conflictpool_->getConflictEntryVector()
                .data() +
            conflictRange.first;
        HighsInt len = conflictRange.second - conflictRange.first;

        return explainBoundChangeConflict(pos, conflict, len);
      }
  }

  return false;
}

bool HighsDomain::ConflictSet::explainBoundChangeConflict(
    HighsInt domchgPos, const HighsDomainChange* conflict, HighsInt len) {
  resolvedDomainChanges.clear();
  auto domchg = localdom.flip(localdom.domchgstack_[domchgPos]);
  bool foundDomchg = false;
  for (HighsInt i = 0; i < len; ++i) {
    if (!foundDomchg && conflict[i] == domchg) {
      foundDomchg = true;
      continue;
    }
    if (globaldom.isActive(conflict[i])) continue;

    HighsInt pos;
    if (conflict[i].boundtype == HighsBoundType::kLower)
      localdom.getColLowerPos(conflict[i].column, domchgPos - 1, pos);
    else
      localdom.getColUpperPos(conflict[i].column, domchgPos - 1, pos);

    if (pos == -1) continue;
    resolvedDomainChanges.push_back(pos);
  }

  return foundDomchg;
}

bool HighsDomain::ConflictSet::explainBoundChangeLeq(
    const HighsDomainChange& domchg, HighsInt pos, const HighsInt* inds,
    const double* vals, HighsInt len, double rhs, double minAct) {
  if (minAct == -kHighsInf) return false;
  // get the coefficient value of the column for which we want to explain
  // the bound change
  double domchgVal = 0;

  resolveBuffer.reserve(len);
  resolveBuffer.clear();
  const auto& nodequeue = localdom.mipsolver->mipdata_->nodequeue;
  for (HighsInt i = 0; i < len; ++i) {
    HighsInt col = inds[i];

    if (col == domchg.column) {
      domchgVal = vals[i];
      continue;
    }

    double delta;
    HighsInt numNodes;
    HighsInt boundpos;
    if (vals[i] > 0) {
      double lb = localdom.getColLowerPos(col, pos, boundpos);
      if (globaldom.col_lower_[col] >= lb) continue;
      delta = vals[i] * (lb - globaldom.col_lower_[col]);
      numNodes = nodequeue.numNodesUp(col);
    } else {
      double ub = localdom.getColUpperPos(col, pos, boundpos);
      if (globaldom.col_upper_[col] <= ub) continue;
      delta = vals[i] * (ub - globaldom.col_upper_[col]);
      numNodes = nodequeue.numNodesDown(col);
    }

    if (boundpos == -1) continue;

    resolveBuffer.emplace_back(delta, numNodes, boundpos);
  }

  if (domchgVal == 0) return false;

  std::sort(resolveBuffer.begin(), resolveBuffer.end(),
            [&](const std::tuple<double, HighsInt, HighsInt>& a,
                const std::tuple<double, HighsInt, HighsInt>& b) {
              double prioA = std::get<0>(a) * (std::get<1>(a) + 1);
              double prioB = std::get<0>(b) * (std::get<1>(b) + 1);
              return prioA > prioB;
            });

  assert(domchgVal != 0);

  // to explain the bound change we start from the bound constraint,
  // multiply it by the columns coefficient in the constraint. Then the
  // bound constraint is a0 * x0 <= b0 and the constraint
  // a0 * x0 + \sum ai * xi <= b. Let the min activity of \sum ai xi be M,
  // then the constraint yields the bound a0 * x0 <= b - M. If M is
  // sufficiently large the bound constraint is implied. Therefore we
  // increase M by updating it with the stronger local bounds until
  // M >= b - b0 holds.
  double b0 = domchg.boundval;
  if (localdom.mipsolver->variableType(domchg.column) !=
      HighsVarType::kContinuous) {
    // in case of an integral variable the bound was rounded and can be
    // relaxed by 1-feastol. We use 1 - 10 * feastol for numerical safety
    if (domchg.boundtype == HighsBoundType::kLower)
      b0 -= (1.0 - 10 * localdom.mipsolver->mipdata_->feastol);
    else
      b0 += (1.0 - 10 * localdom.mipsolver->mipdata_->feastol);
  } else {
    // for a continuous variable we relax the bound by epsilon to
    // accomodate for tiny rounding errors
    if (domchg.boundtype == HighsBoundType::kLower)
      b0 -= localdom.mipsolver->mipdata_->epsilon;
    else
      b0 += localdom.mipsolver->mipdata_->epsilon;
  }

  // now multiply the bound constraint with the coefficient value in the
  // constraint to obtain b0
  b0 *= domchgVal;

  // compute the lower bound of M that is necessary
  double Mlower = rhs - b0;

  // M is the global residual activity initially
  double M = minAct;
  if (domchgVal < 0)
    M -= domchgVal * globaldom.col_upper_[domchg.column];
  else
    M -= domchgVal * globaldom.col_lower_[domchg.column];

  resolvedDomainChanges.clear();
  for (const std::tuple<double, HighsInt, HighsInt>& reasonDomchg :
       resolveBuffer) {
    M += std::get<0>(reasonDomchg);
    resolvedDomainChanges.push_back(std::get<2>(reasonDomchg));
    assert(resolvedDomainChanges.back() >= 0);
    assert(resolvedDomainChanges.back() < localdom.domchgstack_.size());
    if (M > Mlower) break;
  }

  if (M <= Mlower) {
    // printf("local bounds reach only value of %.12g, need at least
    // %.12g\n",
    //        M, Mlower);
    return false;
  }

  return true;
}

void HighsDomain::ConflictSet::pushQueue(HighsInt domchgPos) {
  resolveQueue.push_back(domchgPos);
  std::push_heap(resolveQueue.begin(), resolveQueue.end());
}

HighsInt HighsDomain::ConflictSet::popQueue() {
  assert(!resolveQueue.empty());
  std::pop_heap(resolveQueue.begin(), resolveQueue.end());
  HighsInt elem = resolveQueue.back();
  resolveQueue.pop_back();
  return elem;
}

void HighsDomain::ConflictSet::clearQueue() { resolveQueue.clear(); }

HighsInt HighsDomain::ConflictSet::queueSize() { return resolveQueue.size(); }

bool HighsDomain::ConflictSet::resolvable(HighsInt domChgPos) {
  assert(domChgPos >= 0);
  assert(domChgPos < (HighsInt)localdom.domchgreason_.size());
  // printf("domchgPos: %d\n", domChgPos);
  // printf("stacksize: %ld\n", localdom.domchgreason_.size());
  switch (localdom.domchgreason_[domChgPos].type) {
    case Reason::kBranching:
    case Reason::kUnknown:
      return false;
    default:
      return true;
  }
}

HighsInt HighsDomain::ConflictSet::resolveDepth(std::set<HighsInt>& frontier,
                                                HighsInt depthLevel,
                                                HighsInt stopSize,
                                                HighsInt minResolve,
                                                bool increaseConflictScore) {
  clearQueue();
  HighsInt startPos =
      depthLevel == 0 ? 0 : localdom.branchPos_[depthLevel - 1] + 1;
  auto iterEnd = depthLevel == localdom.branchPos_.size()
                     ? frontier.end()
                     : frontier.upper_bound(localdom.branchPos_[depthLevel]);
  for (auto it = frontier.lower_bound(startPos); it != iterEnd; ++it) {
    assert(it != frontier.end());
    HighsInt domChgPos = *it;
    if (resolvable(domChgPos)) pushQueue(domChgPos);
  }

  HighsInt numResolved = 0;

  while (queueSize() > stopSize ||
         (queueSize() > 0 && numResolved < minResolve)) {
    HighsInt pos = popQueue();
    if (!explainBoundChange(pos)) continue;

    ++numResolved;
    frontier.erase(pos);
    for (HighsInt i : resolvedDomainChanges) {
      if (frontier.insert(i).second) {
        if (increaseConflictScore) {
          if (localdom.domchgstack_[i].boundtype == HighsBoundType::kLower)
            localdom.mipsolver->mipdata_->pseudocost.increaseConflictScoreUp(
                localdom.domchgstack_[i].column);
          else
            localdom.mipsolver->mipdata_->pseudocost.increaseConflictScoreDown(
                localdom.domchgstack_[i].column);
        }
        if (i >= startPos && resolvable(i)) pushQueue(i);
      }
    }
  }

  return numResolved;
}

HighsInt HighsDomain::ConflictSet::computeCuts(
    HighsInt depthLevel, HighsConflictPool& conflictPool) {
  HighsInt numResolved =
      resolveDepth(reasonSideFrontier, depthLevel, 1,
                   depthLevel == localdom.branchPos_.size() ? 1 : 0, true);

  HighsInt numConflicts = 0;
  if (numResolved != 0) {
    // add conflict cut
    localdom.mipsolver->mipdata_->debugSolution.checkConflictReasonFrontier(
        reasonSideFrontier, localdom.domchgstack_);

    conflictPool.addConflictCut(localdom, reasonSideFrontier);
    ++numConflicts;
  }

  // if the queue size is 1 then we have a resolvable UIP that is not the branch
  // vertex
  if (queueSize() == 1) {
    HighsInt uipPos = popQueue();
    clearQueue();

    // compute the UIP reconvergence cut
    reconvergenceFrontier.clear();
    reconvergenceFrontier.insert(uipPos);
    HighsInt numResolved = resolveDepth(reconvergenceFrontier, depthLevel, 0);

    if (numResolved != 0 && reconvergenceFrontier.count(uipPos) == 0) {
      localdom.mipsolver->mipdata_->debugSolution
          .checkConflictReconvergenceFrontier(reconvergenceFrontier, uipPos,
                                              localdom.domchgstack_);
      conflictPool.addReconvergenceCut(localdom, reconvergenceFrontier,
                                       localdom.domchgstack_[uipPos]);
      ++numConflicts;
    }
  }

  return numConflicts;
}

void HighsDomain::ConflictSet::conflictAnalysis(
    HighsConflictPool& conflictPool) {
  resolvedDomainChanges.reserve(localdom.domchgstack_.size());

  if (!explainInfeasibility()) return;

  localdom.mipsolver->mipdata_->pseudocost.increaseConflictWeight();
  for (HighsInt pos : resolvedDomainChanges) {
    if (localdom.domchgstack_[pos].boundtype == HighsBoundType::kLower)
      localdom.mipsolver->mipdata_->pseudocost.increaseConflictScoreUp(
          localdom.domchgstack_[pos].column);
    else
      localdom.mipsolver->mipdata_->pseudocost.increaseConflictScoreDown(
          localdom.domchgstack_[pos].column);
  }

  if (resolvedDomainChanges.size() >
      0.3 * localdom.mipsolver->mipdata_->integral_cols.size())
    return;

  reasonSideFrontier.insert(resolvedDomainChanges.begin(),
                            resolvedDomainChanges.end());

  localdom.mipsolver->mipdata_->debugSolution.checkConflictReasonFrontier(
      reasonSideFrontier, localdom.domchgstack_);

  HighsInt numConflicts = 0;
  for (HighsInt currDepth = localdom.branchPos_.size(); currDepth >= 0;
       --currDepth) {
    if (currDepth > 0) {
      // skip redundant branching changes which are just added for symmetry
      // handling
      HighsInt branchpos = localdom.branchPos_[currDepth - 1];
      if (localdom.domchgstack_[branchpos].boundval ==
          localdom.prevboundval_[branchpos].first)
        continue;
    }
    numConflicts += computeCuts(currDepth, conflictPool);

    if (numConflicts == 0) break;
  }
}

void HighsDomain::ConflictSet::conflictAnalysis(
    const HighsInt* proofinds, const double* proofvals, HighsInt prooflen,
    double proofrhs, HighsConflictPool& conflictPool) {
  resolvedDomainChanges.reserve(localdom.domchgstack_.size());

  HighsInt ninfmin;
  HighsCDouble activitymin;
  globaldom.computeMinActivity(0, prooflen, proofinds, proofvals, ninfmin,
                               activitymin);
  if (ninfmin != 0) return;

  if (!explainInfeasibilityLeq(proofinds, proofvals, prooflen, proofrhs,
                               double(activitymin)))
    return;

  localdom.mipsolver->mipdata_->pseudocost.increaseConflictWeight();
  for (HighsInt pos : resolvedDomainChanges) {
    if (localdom.domchgstack_[pos].boundtype == HighsBoundType::kLower)
      localdom.mipsolver->mipdata_->pseudocost.increaseConflictScoreUp(
          localdom.domchgstack_[pos].column);
    else
      localdom.mipsolver->mipdata_->pseudocost.increaseConflictScoreDown(
          localdom.domchgstack_[pos].column);
  }

  if (resolvedDomainChanges.size() >
      0.3 * localdom.mipsolver->mipdata_->integral_cols.size())
    return;

  reasonSideFrontier.insert(resolvedDomainChanges.begin(),
                            resolvedDomainChanges.end());

  assert(resolvedDomainChanges.size() == reasonSideFrontier.size());

  localdom.mipsolver->mipdata_->debugSolution.checkConflictReasonFrontier(
      reasonSideFrontier, localdom.domchgstack_);

  HighsInt numConflicts = 0;
  for (HighsInt currDepth = localdom.branchPos_.size(); currDepth >= 0;
       --currDepth) {
    if (currDepth > 0) {
      // skip redundant branching changes which are just added for symmetry
      // handling
      HighsInt branchpos = localdom.branchPos_[currDepth - 1];
      if (localdom.domchgstack_[branchpos].boundval ==
          localdom.prevboundval_[branchpos].first)
        continue;
    }
    numConflicts += computeCuts(currDepth, conflictPool);

    // at least in the highest depth level conflicts must be found, otherwise we
    // can immediately stop
    if (numConflicts == 0) break;
  }
}
