/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "presolve/HighsLpPropagator.h"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <queue>

#include "util/HighsUtils.h"

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
    const std::vector<HighsVarType>& integrality_, std::vector<double>& Avalue_,
    std::vector<HighsInt>& Aindex_, std::vector<HighsInt>& Astart_,
    std::vector<HighsInt>& Aend_, std::vector<double>& ARvalue_,
    std::vector<HighsInt>& ARindex_, std::vector<HighsInt>& ARstart_,
    const std::vector<HighsInt>& flagRow, const std::vector<HighsInt>& flagCol,
    std::vector<double>& rowLower_, std::vector<double>& rowUpper_)
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
      integrality_(integrality_),
      colLower_(colLower),
      colUpper_(colUpper) {}

void HighsLpPropagator::computeMinActivity(HighsInt start, HighsInt end,
                                           const HighsInt* ARindex,
                                           const double* ARvalue,
                                           HighsInt& ninfmin,
                                           HighsCDouble& activitymin) {
  activitymin = 0.0;
  ninfmin = 0;
  for (HighsInt j = start; j != end; ++j) {
    HighsInt col = ARindex[j];
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

void HighsLpPropagator::computeMaxActivity(HighsInt start, HighsInt end,
                                           const HighsInt* ARindex,
                                           const double* ARvalue,
                                           HighsInt& ninfmax,
                                           HighsCDouble& activitymax) {
  activitymax = 0.0;
  ninfmax = 0;
  for (HighsInt j = start; j != end; ++j) {
    HighsInt col = ARindex[j];
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

HighsInt HighsLpPropagator::propagateRowUpper(const HighsInt* Rindex,
                                              const double* Rvalue,
                                              HighsInt Rlen, double Rupper,
                                              const HighsCDouble& minactivity,
                                              HighsInt ninfmin,
                                              HighsDomainChange* boundchgs) {
  if (ninfmin > 1) return 0;
  HighsInt numchgs = 0;
  for (HighsInt i = 0; i != Rlen; ++i) {
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

      if (integrality_[Rindex[i]] != HighsVarType::CONTINUOUS) {
        bound = std::floor(bound + 1e-6);
        if (bound < colUpper_[Rindex[i]])
          accept = true;
        else
          accept = false;
      } else {
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
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Upper, Rindex[i], bound};

    } else {
      bool accept;

      if (integrality_[Rindex[i]] != HighsVarType::CONTINUOUS) {
        bound = std::ceil(bound - 1e-6);
        if (bound > colLower_[Rindex[i]])
          accept = true;
        else
          accept = false;
      } else {
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
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Lower, Rindex[i], bound};
    }
  }

  return numchgs;
}

HighsInt HighsLpPropagator::propagateRowLower(const HighsInt* Rindex,
                                              const double* Rvalue,
                                              HighsInt Rlen, double Rlower,
                                              const HighsCDouble& maxactivity,
                                              HighsInt ninfmax,
                                              HighsDomainChange* boundchgs) {
  if (ninfmax > 1) return 0;
  HighsInt numchgs = 0;
  for (HighsInt i = 0; i != Rlen; ++i) {
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

      if (integrality_[Rindex[i]] != HighsVarType::CONTINUOUS) {
        bound = std::floor(bound + 1e-6);
        if (bound < colUpper_[Rindex[i]])
          accept = true;
        else
          accept = false;
      } else {
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
      }

      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Upper, Rindex[i], bound};
    } else {
      bool accept;

      if (integrality_[Rindex[i]] != HighsVarType::CONTINUOUS) {
        bound = std::ceil(bound - 1e-6);
        if (bound > colLower_[Rindex[i]])
          accept = true;
        else
          accept = false;
      } else {
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
      }
      if (accept)
        boundchgs[numchgs++] = {HighsBoundType::Lower, Rindex[i], bound};
    }
  }

  return numchgs;
}

void HighsLpPropagator::updateActivityLbChange(HighsInt col, double oldbound,
                                               double newbound) {
  HighsInt start = Astart_[col];
  HighsInt end = Aend_[col];

  for (HighsInt i = start; i != end; ++i) {
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

void HighsLpPropagator::updateActivityUbChange(HighsInt col, double oldbound,
                                               double newbound) {
  HighsInt start = Astart_[col];
  HighsInt end = Aend_[col];

  for (HighsInt i = start; i != end; ++i) {
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

void HighsLpPropagator::markPropagate(HighsInt row) {
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

  assert(ARstart_.size() == rowLower_.size() + 1);
  const HighsInt numrow = rowLower_.size();
  for (HighsInt i = 0; i != numrow; ++i) {
    if (!flagRow[i]) continue;
    HighsInt start = ARstart_[i];
    HighsInt end = ARstart_[i + 1];

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
  ++numBoundChgs_;

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

HighsInt HighsLpPropagator::propagate() {
  std::vector<HighsInt> propagateinds;

  if (propagateinds_.empty()) return 0;

  HighsInt numchgs = -numBoundChgs_;
  size_t changedboundsize = 2 * ARvalue_.size();
  std::unique_ptr<HighsDomainChange[]> changedbounds(
      new HighsDomainChange[changedboundsize]);

  while (!propagateinds_.empty()) {
    propagateinds.swap(propagateinds_);

    HighsInt propnnz = 0;
    HighsInt numproprows = propagateinds.size();
    for (HighsInt i = 0; i != numproprows; ++i) {
      HighsInt row = propagateinds[i];
      propagateflags_[row] = 0;
      propnnz += ARstart_[i + 1] - ARstart_[i];
    }

    if (!infeasible_) {
      std::vector<HighsInt> propRowNumChangedBounds_(numproprows);

      auto propagateIndex = [&](HighsInt k) {
        // for (HighsInt k = 0; k != numproprows; ++k) {
        HighsInt i = propagateinds[k];
        HighsInt start = ARstart_[i];
        HighsInt end = ARstart_[i + 1];
        HighsInt Rlen = end - start;
        const HighsInt* Rindex = &ARindex_[start];
        const double* Rvalue = &ARvalue_[start];
        HighsInt numchgs = 0;

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

      // printf("numproprows (model): %" HIGHSINT_FORMAT "\n", numproprows);

      for (HighsInt k = 0; k != numproprows; ++k) propagateIndex(k);

      for (HighsInt k = 0; k != numproprows; ++k) {
        if (propRowNumChangedBounds_[k] == 0) continue;
        HighsInt i = propagateinds[k];
        HighsInt start = 2 * ARstart_[i];
        HighsInt end = start + propRowNumChangedBounds_[k];
        for (HighsInt j = start; j != end; ++j) changeBound(changedbounds[j]);

        if (infeasible_) break;
      }
    }

    propagateinds.clear();
  }

  numchgs += numBoundChgs_;
  return numchgs;
}

HighsInt HighsLpPropagator::tightenCoefficients() {
  HighsInt numrow = rowUpper_.size();
  HighsInt ntightenedtotal = 0;
  for (HighsInt i = 0; i != numrow; ++i) {
    if (!flagRow[i] ||
        (rowUpper_[i] != HIGHS_CONST_INF && rowLower_[i] != -HIGHS_CONST_INF))
      continue;

    HighsInt scale;

    if (rowUpper_[i] != HIGHS_CONST_INF) {
      if (activitymaxinf_[i] != 0) continue;

      if (activitymax_[i] - rowUpper_[i] <= 1e-6) continue;

      scale = 1;
    } else {
      if (activitymininf_[i] != 0) continue;

      if (rowLower_[i] - activitymin_[i] <= 1e-6) continue;

      scale = -1;
    }

    HighsCDouble maxactivity = scale == 1 ? activitymax_[i] : -activitymin_[i];
    HighsCDouble upper = scale == 1 ? rowUpper_[i] : -rowLower_[i];
    HighsCDouble maxabscoef = double(maxactivity - upper);
    HighsInt tightened = 0;
    const HighsInt start = ARstart_[i];
    const HighsInt end = ARstart_[i + 1];

    for (HighsInt j = start; j != end; ++j) {
      HighsInt col = ARindex_[j];
      if (!flagCol[col] || integrality_[col] == HighsVarType::CONTINUOUS)
        continue;

      double val = scale * ARvalue_[j];
      if (val > maxabscoef) {
        HighsCDouble delta = val - maxabscoef;
        upper -= delta * colUpper_[col];
        ARvalue_[j] = scale * double(maxabscoef);
        ++tightened;
      } else if (val < -maxabscoef) {
        HighsCDouble delta = -val - maxabscoef;
        upper += delta * colLower_[col];
        ARvalue_[j] = -scale * double(maxabscoef);
        ++tightened;
      }
    }

    if (tightened != 0) {
      if (scale == 1)
        rowUpper_[i] = double(upper);
      else
        rowLower_[i] = -double(upper);
      // printf("tightened %" HIGHSINT_FORMAT " coefficients, rhs changed from
      // %g to %g\n",
      //       tightened, rhs, double(upper));
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

      ntightenedtotal += tightened;
    }
  }

  if (ntightenedtotal != 0) {
    HighsInt transNcol = numrow;
    HighsInt transNrow = colLower_.size();

    highsSparseTranspose(transNrow, transNcol, ARstart_, ARindex_, ARvalue_,
                         Astart_, Aindex_, Avalue_);
    std::copy(Astart_.begin() + 1, Astart_.end(), Aend_.begin());
  }

  return ntightenedtotal;
}
