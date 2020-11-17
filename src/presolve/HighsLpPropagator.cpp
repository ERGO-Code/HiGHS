#include "presolve/HighsLpPropagator.h"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <queue>

#include "lp_data/HConst.h"

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

HighsLpPropagator::HighsLpPropagator(
    const std::vector<double>& colLower, const std::vector<double>& colUpper,
    const std::vector<double>& Avalue_, const std::vector<int>& Aindex_,
    const std::vector<int>& Astart_, const std::vector<int>& Aend_,
    const std::vector<double>& ARvalue_, const std::vector<int>& ARindex_,
    const std::vector<int>& ARstart_, const std::vector<int>& flagRow,
    const std::vector<int>& flagCol, const std::vector<double>& rowLower_,
    const std::vector<double>& rowUpper_)
    : Avalue_(Avalue_),
      Aindex_(Aindex_),
      Astart_(Astart_),
      Aend_(Aend_),
      ARvalue_(ARvalue_),
      ARindex_(ARindex_),
      ARstart_(ARstart_),
      flagRow(flagRow),
      flagCol(flagCol),
      rowLower_(rowLower_),
      rowUpper_(rowUpper_),
      colLower_(colLower),
      colUpper_(colUpper) {}

void HighsLpPropagator::computeMinActivity(int start, int end,
                                           const int* ARindex,
                                           const double* ARvalue, int& ninfmin,
                                           HighsCDouble& activitymin) {
  activitymin = 0.0;
  ninfmin = 0;
  for (int j = start; j != end; ++j) {
    int col = ARindex[j];
    if (!flagCol[col]) continue;
    double val = ARvalue[j];

    double contributionmin =
        activityContributionMin(val, colLower_[col], colUpper_[col]);

    if (contributionmin == -HIGHS_CONST_INF)
      ++ninfmin;
    else
      activitymin += contributionmin;
  }

  activitymin.renormalize();
}

void HighsLpPropagator::computeMaxActivity(int start, int end,
                                           const int* ARindex,
                                           const double* ARvalue, int& ninfmax,
                                           HighsCDouble& activitymax) {
  activitymax = 0.0;
  ninfmax = 0;
  for (int j = start; j != end; ++j) {
    int col = ARindex[j];
    if (!flagCol[col]) continue;
    double val = ARvalue[j];

    double contributionmin =
        activityContributionMax(val, colLower_[col], colUpper_[col]);

    if (contributionmin == HIGHS_CONST_INF)
      ++ninfmax;
    else
      activitymax += contributionmin;
  }

  activitymax.renormalize();
}

int HighsLpPropagator::propagateRowUpper(const int* Rindex,
                                         const double* Rvalue, int Rlen,
                                         double Rupper,
                                         const HighsCDouble& minactivity,
                                         int ninfmin,
                                         HighsDomainChange* boundchgs) {
  if (ninfmin > 1) return 0;
  int numchgs = 0;
  for (int i = 0; i != Rlen; ++i) {
    if (!flagCol[Rindex[i]]) continue;
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

      if (colUpper_[Rindex[i]] == HIGHS_CONST_INF)
        accept = true;
      else if (bound + 1e-3 < colUpper_[Rindex[i]] &&
               (colUpper_[Rindex[i]] - bound) /
                       std::max(std::abs(colUpper_[Rindex[i]]),
                                std::abs(bound)) >
                   0.3)
        accept = true;
      else
        accept = false;

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Upper, Rindex[i], bound};

    } else {
      bool accept;

      if (colLower_[Rindex[i]] == -HIGHS_CONST_INF)
        accept = true;
      else if (bound - 1e-3 > colLower_[Rindex[i]] &&
               (bound - colLower_[Rindex[i]]) /
                       std::max(std::abs(colUpper_[Rindex[i]]),
                                std::abs(bound)) >
                   0.3)
        accept = true;
      else
        accept = false;

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Lower, Rindex[i], bound};
    }
  }

  return numchgs;
}

int HighsLpPropagator::propagateRowLower(const int* Rindex,
                                         const double* Rvalue, int Rlen,
                                         double Rlower,
                                         const HighsCDouble& maxactivity,
                                         int ninfmax,
                                         HighsDomainChange* boundchgs) {
  if (ninfmax > 1) return 0;
  int numchgs = 0;
  for (int i = 0; i != Rlen; ++i) {
    if (!flagCol[Rindex[i]]) continue;
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

      if (colUpper_[Rindex[i]] == HIGHS_CONST_INF)
        accept = true;
      else if (bound + 1e-3 < colUpper_[Rindex[i]] &&
               (colUpper_[Rindex[i]] - bound) /
                       std::max(std::abs(colUpper_[Rindex[i]]),
                                std::abs(bound)) >
                   0.3)
        accept = true;
      else
        accept = false;

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Upper, Rindex[i], bound};
    } else {
      bool accept;

      if (colLower_[Rindex[i]] == -HIGHS_CONST_INF)
        accept = true;
      else if (bound - 1e-3 > colLower_[Rindex[i]] &&
               (bound - colLower_[Rindex[i]]) /
                       std::max(std::abs(colUpper_[Rindex[i]]),
                                std::abs(bound)) >
                   0.3)
        accept = true;
      else
        accept = false;

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Lower, Rindex[i], bound};
    }
  }

  return numchgs;
}

void HighsLpPropagator::updateActivityLbChange(int col, double oldbound,
                                               double newbound) {
  int start = Astart_[col];
  int end = Aend_[col];

  for (int i = start; i != end; ++i) {
    if (!flagRow[Aindex_[i]]) continue;
    if (Avalue_[i] > 0) {
      double deltamin;
      if (oldbound == -HIGHS_CONST_INF) {
        --activitymininf_[Aindex_[i]];
        deltamin = newbound * Avalue_[i];
      } else {
        deltamin = (newbound - oldbound) * Avalue_[i];
      }
      activitymin_[Aindex_[i]] += deltamin;

      if (rowUpper_[Aindex_[i]] != HIGHS_CONST_INF &&
          activitymininf_[Aindex_[i]] == 0 &&
          activitymin_[Aindex_[i]] - rowUpper_[Aindex_[i]] > 1e-6) {
        infeasible_ = Aindex_[i] + 1;
      }

      if (deltamin > 0 && activitymininf_[Aindex_[i]] <= 1 &&
          !propagateflags_[Aindex_[i]] &&
          rowUpper_[Aindex_[i]] != HIGHS_CONST_INF) {
        markPropagate(Aindex_[i]);
        // propagateflags_[Aindex_[i]] = 1;
        // propagateinds_.push_back(Aindex_[i]);
      }
    } else {
      double deltamax;
      if (oldbound == -HIGHS_CONST_INF) {
        --activitymaxinf_[Aindex_[i]];
        deltamax = newbound * Avalue_[i];
      } else {
        deltamax = (newbound - oldbound) * Avalue_[i];
      }
      activitymax_[Aindex_[i]] += deltamax;

      if (rowLower_[Aindex_[i]] != -HIGHS_CONST_INF &&
          activitymaxinf_[Aindex_[i]] == 0 &&
          rowLower_[Aindex_[i]] - activitymax_[Aindex_[i]] > 1e-6) {
        infeasible_ = Aindex_[i] + 1;
      }

      if (deltamax < 0 && activitymaxinf_[Aindex_[i]] <= 1 &&
          !propagateflags_[Aindex_[i]] &&
          rowLower_[Aindex_[i]] != -HIGHS_CONST_INF) {
        markPropagate(Aindex_[i]);
        // propagateflags_[Aindex_[i]] = 1;
        // propagateinds_.push_back(Aindex_[i]);
      }
    }
  }
}

void HighsLpPropagator::updateActivityUbChange(int col, double oldbound,
                                               double newbound) {
  int start = Astart_[col];
  int end = Aend_[col];

  for (int i = start; i != end; ++i) {
    if (!flagRow[Aindex_[i]]) continue;
    if (Avalue_[i] > 0) {
      double deltamax;
      if (oldbound == HIGHS_CONST_INF) {
        --activitymaxinf_[Aindex_[i]];
        deltamax = newbound * Avalue_[i];
      } else {
        deltamax = (newbound - oldbound) * Avalue_[i];
      }
      activitymax_[Aindex_[i]] += deltamax;

      if (rowLower_[Aindex_[i]] != -HIGHS_CONST_INF &&
          activitymaxinf_[Aindex_[i]] == 0 &&
          rowLower_[Aindex_[i]] - activitymax_[Aindex_[i]] > 1e-6) {
        infeasible_ = Aindex_[i] + 1;
      }

      if (deltamax < 0 && activitymaxinf_[Aindex_[i]] <= 1 &&
          !propagateflags_[Aindex_[i]] &&
          rowLower_[Aindex_[i]] != -HIGHS_CONST_INF) {
        markPropagate(Aindex_[i]);
        // propagateflags_[Aindex_[i]] = 1;
        // propagateinds_.push_back(Aindex_[i]);
      }
    } else {
      double deltamin;
      if (oldbound == HIGHS_CONST_INF) {
        --activitymininf_[Aindex_[i]];
        deltamin = newbound * Avalue_[i];
      } else {
        deltamin = (newbound - oldbound) * Avalue_[i];
      }

      activitymin_[Aindex_[i]] += deltamin;

      if (rowUpper_[Aindex_[i]] != HIGHS_CONST_INF &&
          activitymininf_[Aindex_[i]] == 0 &&
          activitymin_[Aindex_[i]] - rowUpper_[Aindex_[i]] > 1e-6) {
        infeasible_ = Aindex_[i] + 1;
      }

      if (deltamin > 0 && activitymininf_[Aindex_[i]] <= 1 &&
          !propagateflags_[Aindex_[i]] &&
          rowUpper_[Aindex_[i]] != HIGHS_CONST_INF) {
        markPropagate(Aindex_[i]);
        // propagateflags_[Aindex_[i]] = 1;
        // propagateinds_.push_back(Aindex_[i]);
      }
    }
  }
}

void HighsLpPropagator::markPropagate(int row) {
  // todo, check if std::min(maxactivity - lhs, rhs - minactivity) <  amax -
  // feastol and only mark in that case

  if (!propagateflags_[row] && flagRow[row]) {
    bool proplower = rowLower_[row] != -HIGHS_CONST_INF;
    bool propupper = rowUpper_[row] != HIGHS_CONST_INF;

    if (proplower || propupper) {
      propagateinds_.push_back(row);
      propagateflags_[row] = 1;
    }
  }
}

void HighsLpPropagator::computeRowActivities() {
  activitymin_.resize(rowLower_.size());
  activitymininf_.resize(rowLower_.size());
  activitymax_.resize(rowLower_.size());
  activitymaxinf_.resize(rowLower_.size());
  propagateflags_.resize(rowLower_.size());
  propagateinds_.reserve(rowLower_.size());

  for (int i = 0; i != rowLower_.size(); ++i) {
    if (!flagRow[i]) continue;
    int start = ARstart_[i];
    int end = ARstart_[i + 1];

    computeMinActivity(start, end, ARindex_.data(), ARvalue_.data(),
                       activitymininf_[i], activitymin_[i]);
    computeMaxActivity(start, end, ARindex_.data(), ARvalue_.data(),
                       activitymaxinf_[i], activitymax_[i]);

    if ((activitymininf_[i] <= 1 && rowUpper_[i] != HIGHS_CONST_INF) ||
        (activitymaxinf_[i] <= 1 && rowLower_[i] != -HIGHS_CONST_INF)) {
      markPropagate(i);
      // propagateflags_[i] = 1;
      // propagateinds_.push_back(i);
    }
  }
}

double HighsLpPropagator::doChangeBound(const HighsDomainChange& boundchg) {
  double oldbound;

  if (boundchg.boundtype == HighsBoundType::Lower) {
    oldbound = colLower_[boundchg.column];
    colLower_[boundchg.column] = boundchg.boundval;
    updateActivityLbChange(boundchg.column, oldbound, boundchg.boundval);
  } else {
    oldbound = colUpper_[boundchg.column];
    colUpper_[boundchg.column] = boundchg.boundval;
    updateActivityUbChange(boundchg.column, oldbound, boundchg.boundval);
  }

  return oldbound;
}

void HighsLpPropagator::changeBound(HighsDomainChange boundchg) {
  assert(boundchg.column >= 0);
  if (boundchg.boundtype == HighsBoundType::Lower) {
    if (boundchg.boundval > colUpper_[boundchg.column]) {
      if (boundchg.boundval - colUpper_[boundchg.column] > 1e-6) {
        infeasible_ = true;
        return;
      }

      boundchg.boundval = colUpper_[boundchg.column];
      if (boundchg.boundval == colLower_[boundchg.column]) return;
    }
  } else {
    if (boundchg.boundval < colLower_[boundchg.column]) {
      if (colLower_[boundchg.column] - boundchg.boundval > 1e-6) {
        infeasible_ = true;
        return;
      }

      boundchg.boundval = colLower_[boundchg.column];
      if (boundchg.boundval == colUpper_[boundchg.column]) return;
    }
  }

  doChangeBound(boundchg);
}

void HighsLpPropagator::propagate() {
  std::vector<int> propagateinds;

  if (propagateinds_.empty()) return;

  size_t changedboundsize = 2 * ARvalue_.size();
  std::unique_ptr<HighsDomainChange[]> changedbounds(
      new HighsDomainChange[changedboundsize]);

  while (!propagateinds_.empty()) {
    propagateinds.swap(propagateinds_);

    int propnnz = 0;
    int numproprows = propagateinds.size();
    for (int i = 0; i != numproprows; ++i) {
      int row = propagateinds[i];
      propagateflags_[row] = 0;
      propnnz += ARstart_[i + 1] - ARstart_[i];
    }

    if (!infeasible_) {
      std::vector<int> propRowNumChangedBounds_(numproprows);

      auto propagateIndex = [&](int k) {
        // for (int k = 0; k != numproprows; ++k) {
        int i = propagateinds[k];
        int start = ARstart_[i];
        int end = ARstart_[i + 1];
        int Rlen = end - start;
        const int* Rindex = &ARindex_[start];
        const double* Rvalue = &ARvalue_[start];
        int numchgs = 0;

        if (rowUpper_[i] != HIGHS_CONST_INF) {
          // computeMinActivity(start, end, mipsolver->ARstart_.data(),
          // mipsolver->ARvalue_.data(), activitymininf_[i],
          //           activitymin_[i]);
          activitymin_[i].renormalize();
          numchgs = propagateRowUpper(Rindex, Rvalue, Rlen, rowUpper_[i],
                                      activitymin_[i], activitymininf_[i],
                                      &changedbounds[2 * start]);
        }

        if (rowLower_[i] != -HIGHS_CONST_INF) {
          // computeMaxActivity(start, end, mipsolver->ARstart_.data(),
          // mipsolver->ARvalue_.data(), activitymaxinf_[i],
          //           activitymax_[i]);
          activitymax_[i].renormalize();
          numchgs += propagateRowLower(Rindex, Rvalue, Rlen, rowLower_[i],
                                       activitymax_[i], activitymaxinf_[i],
                                       &changedbounds[2 * start + numchgs]);
        }

        propRowNumChangedBounds_[k] = numchgs;
      };

      // printf("numproprows (model): %d\n", numproprows);

      for (int k = 0; k != numproprows; ++k) propagateIndex(k);

      for (int k = 0; k != numproprows; ++k) {
        if (propRowNumChangedBounds_[k] == 0) continue;
        int i = propagateinds[k];
        int start = 2 * ARstart_[i];
        int end = start + propRowNumChangedBounds_[k];
        for (int j = start; j != end; ++j) changeBound(changedbounds[j]);

        if (infeasible_) break;
      }
    }

    propagateinds.clear();
  }
}
