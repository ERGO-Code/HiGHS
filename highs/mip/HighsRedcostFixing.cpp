/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsRedcostFixing.h"

#include "mip/HighsMipSolverData.h"

std::vector<std::pair<double, HighsDomainChange>>
HighsRedcostFixing::getLurkingBounds(const HighsMipSolver& mipsolver,
                                     const HighsDomain& globaldom) const {
  std::vector<std::pair<double, HighsDomainChange>> domchgs;
  if (lurkingColLower.empty()) return domchgs;

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    for (const auto& lower : lurkingColLower[col]) {
      if (lower.second > globaldom.col_lower_[col])
        domchgs.emplace_back(
            lower.first, HighsDomainChange{static_cast<double>(lower.second),
                                           col, HighsBoundType::kLower});
    }

    for (const auto& upper : lurkingColUpper[col]) {
      if (upper.second < globaldom.col_upper_[col])
        domchgs.emplace_back(
            upper.first, HighsDomainChange{static_cast<double>(upper.second),
                                           col, HighsBoundType::kUpper});
    }
  }

  return domchgs;
}

void HighsRedcostFixing::propagateRootRedcost(const HighsMipSolver& mipsolver) {
  if (lurkingColLower.empty()) return;

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    lurkingColLower[col].erase(
        lurkingColLower[col].begin(),
        lurkingColLower[col].upper_bound(mipsolver.mipdata_->lower_bound));
    lurkingColUpper[col].erase(
        lurkingColUpper[col].begin(),
        lurkingColUpper[col].upper_bound(mipsolver.mipdata_->lower_bound));

    // By the staircase invariant (see addRootRedcost) the first active entry
    // carries the tightest bound; later active entries are implied by it.
    for (auto it =
             lurkingColLower[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColLower[col].end(); it = lurkingColLower[col].end()) {
      if (it->second > mipsolver.mipdata_->getDomain().col_lower_[col]) {
        mipsolver.mipdata_->getDomain().changeBound(
            HighsBoundType::kLower, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->getDomain().infeasible()) return;
      }
    }

    // By the staircase invariant (see addRootRedcost) the first active entry
    // carries the tightest bound; later active entries are implied by it.
    for (auto it =
             lurkingColUpper[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColUpper[col].end(); it = lurkingColUpper[col].end()) {
      if (it->second < mipsolver.mipdata_->getDomain().col_upper_[col]) {
        mipsolver.mipdata_->getDomain().changeBound(
            HighsBoundType::kUpper, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->getDomain().infeasible()) return;
      }
    }
  }

  mipsolver.mipdata_->getDomain().propagate();
}

void HighsRedcostFixing::propagateRedCost(const HighsMipSolver& mipsolver,
                                          HighsDomain& localdomain,
                                          HighsDomain& globaldom,
                                          const HighsLpRelaxation& lp,
                                          HighsConflictPool& conflictpool,
                                          HighsPseudocost& pseudocost,
                                          double upper_limit) {
  const std::vector<double>& lpredcost = lp.getSolution().col_dual;
  double lpobjective = lp.getObjective();
  HighsCDouble gap = static_cast<HighsCDouble>(upper_limit) - lpobjective;

  double tolerance = std::max(10 * mipsolver.mipdata_->feastol,
                              mipsolver.mipdata_->epsilon * double(gap));

  assert(!localdomain.infeasible());
  std::vector<HighsDomainChange> boundChanges;
  boundChanges.reserve(mipsolver.mipdata_->integral_cols.size());
  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    // lpobj + (col - bnd) * redcost <= cutoffbound
    // (col - bnd) * redcost <= gap
    // redcost * col <= gap + redcost * bnd
    //   case bnd is upper bound  => redcost < 0 :
    //      col >= gap/redcost + ub
    //      <=> redcost < gap / (lb - ub)
    //   case bnd is lower bound  => redcost > 0 :
    //      col <= gap/redcost + lb
    //      <=> redcost > gap / (ub - lb)
    if (localdomain.col_upper_[col] == localdomain.col_lower_[col]) continue;
    if (std::fabs(lpredcost[col]) <= tolerance) continue;

    double maxIncrease = lpredcost[col] * (localdomain.col_upper_[col] -
                                           localdomain.col_lower_[col]);
    if (maxIncrease > gap) {
      assert(localdomain.col_lower_[col] != -kHighsInf);
      assert(lpredcost[col] > tolerance);
      double newub =
          double(floor(gap / lpredcost[col] + localdomain.col_lower_[col] +
                       mipsolver.mipdata_->feastol));
      if (newub >= localdomain.col_upper_[col]) continue;
      assert(newub < localdomain.col_upper_[col]);

      if (globaldom.isBinary(col)) {
        boundChanges.emplace_back(
            HighsDomainChange{newub, col, HighsBoundType::kUpper});
      } else {
        localdomain.changeBound(HighsBoundType::kUpper, col, newub,
                                HighsDomain::Reason::unspecified());
        if (localdomain.infeasible()) return;
      }
    } else if (maxIncrease < -gap) {
      assert(localdomain.col_upper_[col] != kHighsInf);
      assert(lpredcost[col] < -tolerance);
      double newlb =
          double(ceil(gap / lpredcost[col] + localdomain.col_upper_[col] -
                      mipsolver.mipdata_->feastol));

      if (newlb <= localdomain.col_lower_[col]) continue;
      assert(newlb > localdomain.col_lower_[col]);

      if (globaldom.isBinary(col)) {
        boundChanges.emplace_back(
            HighsDomainChange{newlb, col, HighsBoundType::kLower});
      } else {
        localdomain.changeBound(HighsBoundType::kLower, col, newlb,
                                HighsDomain::Reason::unspecified());
        if (localdomain.infeasible()) return;
      }
    }
  }

  if (!boundChanges.empty()) {
    std::vector<HighsInt> inds;
    std::vector<double> vals;
    double rhs;

    if (boundChanges.size() <= 100 &&
        lp.computeDualProof(globaldom, upper_limit, inds, vals, rhs, false)) {
      bool addedConstraints = false;

      HighsInt oldNumConflicts = conflictpool.getNumConflicts();
      for (const HighsDomainChange& domchg : boundChanges) {
        if (localdomain.isActive(domchg)) continue;
        localdomain.conflictAnalyzeReconvergence(
            domchg, inds.data(), vals.data(),
            static_cast<HighsInt>(inds.size()), rhs, conflictpool, globaldom,
            pseudocost);
      }
      addedConstraints = conflictpool.getNumConflicts() != oldNumConflicts;

      if (addedConstraints) {
        localdomain.propagate();
        if (localdomain.infeasible()) return;

        boundChanges.erase(
            std::remove_if(boundChanges.begin(), boundChanges.end(),
                           [&](const HighsDomainChange& domchg) {
                             return localdomain.isActive(domchg);
                           }),
            boundChanges.end());
      }

      if (!boundChanges.empty()) {
        for (const HighsDomainChange& domchg : boundChanges) {
          localdomain.changeBound(domchg, HighsDomain::Reason::unspecified());
          if (localdomain.infeasible()) break;
        }

        if (!localdomain.infeasible()) localdomain.propagate();
      }
      // /printf("numConflicts: %d\n", numConflicts);
    } else {
      for (const HighsDomainChange& domchg : boundChanges) {
        localdomain.changeBound(domchg, HighsDomain::Reason::unspecified());
        if (localdomain.infeasible()) break;
      }

      if (!localdomain.infeasible()) localdomain.propagate();
    }
  }
}

void HighsRedcostFixing::addRootRedcost(const HighsMipSolver& mipsolver,
                                        const std::vector<double>& lpredcost,
                                        double lpobjective) {
  lurkingColLower.resize(mipsolver.numCol());
  lurkingColUpper.resize(mipsolver.numCol());

  // Provided domains won't be used (only used for dual proof)
  mipsolver.mipdata_->getLp().computeBasicDegenerateDuals(
      mipsolver.mipdata_->feastol, mipsolver.mipdata_->getDomain(),
      mipsolver.mipdata_->getDomain(), mipsolver.mipdata_->getConflictPool(),
      mipsolver.mipdata_->getPseudoCost(), false);

  // Compute maximum number of steps per column with large domain
  // max_steps = 2 ** k, k = max(5, min(10 ,round(log(|D| / 10)))),
  // D = {col : integral_cols | (ub - lb) >= 512}
  // This is to avoid doing 2**10 steps when there's many unbounded columns
  HighsInt numRedcostLargeDomainCols = 0;
  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    if ((mipsolver.mipdata_->getDomain().col_upper_[col] -
         mipsolver.mipdata_->getDomain().col_lower_[col]) >= 512 &&
        std::abs(lpredcost[col]) > mipsolver.mipdata_->feastol) {
      numRedcostLargeDomainCols++;
    }
  }
  HighsInt maxNumStepsExp = 10;
  int expshift = 0;
  std::frexp(numRedcostLargeDomainCols / 10, &expshift);
  if (expshift > 5) {
    expshift = std::min(expshift, static_cast<int>(maxNumStepsExp));
    maxNumStepsExp = maxNumStepsExp - expshift + 5;
  }
  HighsInt maxNumSteps = static_cast<HighsInt>(1ULL << maxNumStepsExp);

  // lambda for finding lurking bounds
  auto findLurkingBounds =
      [&](HighsInt col, HighsInt direction, HighsInt bound, HighsInt otherBound,
          bool isOtherBoundFinite, double lpObjective, double redCost,
          HighsInt maxNumSteps, HighsInt maxNumStepsExp,
          std::multimap<double, HighsInt>& lurkingBounds,
          std::multimap<double, HighsInt>& otherLurkingBounds) {
        if (direction * redCost == kHighsInf) {
          lurkingBounds.clear();
          otherLurkingBounds.clear();
          lurkingBounds.emplace(-kHighsInf, bound);
          return;
        }

        HighsInt lastBound;
        if (!isOtherBoundFinite)
          lastBound = bound + direction * maxNumSteps;
        else
          lastBound = otherBound - direction;

        HighsInt step = 1;
        HighsInt range = direction * (lastBound - bound);
        if (range > maxNumSteps)
          step = (range + maxNumSteps - 1) >> maxNumStepsExp;
        double shift = direction * (1 - 10 * mipsolver.mipdata_->feastol);
        step *= direction;

        // requiredCutoffBound grows linearly with the step index j, so jump
        // straight to the first j whose cutoff can exceed the lower bound
        // (the in-loop test still guards the boundary). Exact.
        HighsInt startBound = bound;
        {
          double denom = double(step) * redCost;
          double lbnd = mipsolver.mipdata_->lower_bound +
                        mipsolver.mipdata_->feastol;
          if (denom > 0 && lbnd > -kHighsInf) {
            double j = std::floor((lbnd - lpObjective - shift * redCost) / denom) - 1;
            if (j > 0) {
              double maxj = double(direction * (lastBound - bound)) /
                                double(direction * step) + 1;
              if (j > maxj) j = maxj;
              startBound = bound + HighsInt(j) * step;
            }
          }
        }
        for (HighsInt lurkingBound = startBound;
             direction * lurkingBound <= direction * lastBound;
             lurkingBound += step) {
          double fracBound = lurkingBound - bound + shift;
          double requiredCutoffBound = fracBound * redCost + lpObjective;
          if (requiredCutoffBound <
              mipsolver.mipdata_->lower_bound + mipsolver.mipdata_->feastol)
            continue;

          // Invariant (kept by this function and by the prefix erases in
          // propagateRootRedcost): along increasing key, direction*bound is
          // strictly increasing. Hence the entry at pos carries the smallest
          // bound among keys >= requiredCutoffBound, and the entries dominated
          // by the new one form a contiguous block right before pos. This
          // replaces the former O(n) scans per insertion (O(n^2) per column
          // with up to 1024 steps) by O(log n + #erased); result is identical.
          auto pos = lurkingBounds.lower_bound(requiredCutoffBound);
          if (pos != lurkingBounds.end() &&
              !(direction * pos->second >= direction * (lurkingBound + step)))
            continue;

          // store this lurking bound and drop the entries it dominates
          auto it = lurkingBounds.emplace_hint(pos, requiredCutoffBound,
                                               lurkingBound);
          while (it != lurkingBounds.begin()) {
            auto prev = std::prev(it);
            if (direction * prev->second >= direction * lurkingBound)
              lurkingBounds.erase(prev);
            else
              break;
          }
        }
      };

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    if (lpredcost[col] > mipsolver.mipdata_->feastol) {
      // col <= (cutoffbound - lpobj)/redcost + lb
      // so for lurkub = lb to ub - 1 we can compute the necessary cutoff
      // bound to reach this bound which is:
      //  lurkub = (cutoffbound - lpobj)/redcost + lb
      //  cutoffbound = (lurkub - lb) * redcost + lpobj
      findLurkingBounds(
          col, HighsInt{1},
          static_cast<HighsInt>(
              mipsolver.mipdata_->getDomain().col_lower_[col]),
          static_cast<HighsInt>(
              mipsolver.mipdata_->getDomain().col_upper_[col]),
          mipsolver.mipdata_->getDomain().col_upper_[col] != kHighsInf,
          lpobjective, lpredcost[col], maxNumSteps, maxNumStepsExp,
          lurkingColUpper[col], lurkingColLower[col]);
    } else if (lpredcost[col] < -mipsolver.mipdata_->feastol) {
      // col >= (cutoffbound - lpobj)/redcost + ub
      // so for lurklb = lb + 1 to ub we can compute the necessary cutoff
      // bound to reach this bound which is:
      //  lurklb = (cutoffbound - lpobj)/redcost + ub
      //  cutoffbound = (lurklb - ub) * redcost + lpobj
      findLurkingBounds(
          col, HighsInt{-1},
          static_cast<HighsInt>(
              mipsolver.mipdata_->getDomain().col_upper_[col]),
          static_cast<HighsInt>(
              mipsolver.mipdata_->getDomain().col_lower_[col]),
          mipsolver.mipdata_->getDomain().col_lower_[col] != -kHighsInf,
          lpobjective, lpredcost[col], maxNumSteps, maxNumStepsExp,
          lurkingColLower[col], lurkingColUpper[col]);
    }
  }
}
