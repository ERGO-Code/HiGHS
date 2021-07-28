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
#include "mip/HighsRedcostFixing.h"

#include "mip/HighsMipSolverData.h"

std::vector<std::pair<double, HighsDomainChange>>
HighsRedcostFixing::getLurkingBounds(const HighsMipSolver& mipsolver) const {
  std::vector<std::pair<double, HighsDomainChange>> domchgs;
  if (lurkingColLower.empty()) return domchgs;

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    for (auto it = lurkingColLower[col].begin();
         it != lurkingColLower[col].end(); ++it) {
      if (it->second > mipsolver.mipdata_->domain.col_lower_[col])
        domchgs.emplace_back(
            it->first,
            HighsDomainChange{(double)it->second, col, HighsBoundType::kLower});
    }

    for (auto it = lurkingColUpper[col].begin();
         it != lurkingColUpper[col].end(); ++it) {
      if (it->second < mipsolver.mipdata_->domain.col_upper_[col])
        domchgs.emplace_back(
            it->first,
            HighsDomainChange{(double)it->second, col, HighsBoundType::kUpper});
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

    for (auto it =
             lurkingColLower[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColLower[col].end(); ++it) {
      if (it->second > mipsolver.mipdata_->domain.col_lower_[col]) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::kLower, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }
    }

    for (auto it =
             lurkingColUpper[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColUpper[col].end(); ++it) {
      if (it->second < mipsolver.mipdata_->domain.col_upper_[col]) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::kUpper, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }
    }
  }

  mipsolver.mipdata_->domain.propagate();
}

void HighsRedcostFixing::propagateRedCost(const HighsMipSolver& mipsolver,
                                          HighsDomain& localdomain,
                                          const HighsLpRelaxation& lp) {
  const std::vector<double>& lpredcost = lp.getSolution().col_dual;
  double lpobjective = lp.getObjective();
  HighsCDouble gap =
      HighsCDouble(mipsolver.mipdata_->upper_limit) - lpobjective;

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
    //   case bnd is lower bound  => redcost > 0 :
    //      col <= gap/redcost + lb
    // Therefore for a tightening to be possible we need:
    //   redcost >= / <=  gap * (ub - lb)
    if (localdomain.col_upper_[col] == localdomain.col_lower_[col]) continue;

    double threshold = double((HighsCDouble(localdomain.col_upper_[col]) -
                               localdomain.col_lower_[col]) *
                                  gap +
                              tolerance);

    if ((localdomain.col_upper_[col] == kHighsInf &&
         lpredcost[col] > tolerance) ||
        lpredcost[col] > threshold) {
      assert(localdomain.col_lower_[col] != -kHighsInf);
      assert(lpredcost[col] > tolerance);
      double newub =
          double(floor(gap / lpredcost[col] + localdomain.col_lower_[col] +
                       mipsolver.mipdata_->feastol));
      if (newub >= localdomain.col_upper_[col]) continue;
      assert(newub < localdomain.col_upper_[col]);

      if (mipsolver.mipdata_->domain.isBinary(col)) {
        boundChanges.emplace_back(
            HighsDomainChange{newub, col, HighsBoundType::kUpper});
      } else {
        localdomain.changeBound(HighsBoundType::kUpper, col, newub,
                                HighsDomain::Reason::unspecified());
        if (localdomain.infeasible()) return;
      }
    } else if ((localdomain.col_lower_[col] == -kHighsInf &&
                lpredcost[col] < -tolerance) ||
               lpredcost[col] < -threshold) {
      assert(localdomain.col_upper_[col] != kHighsInf);
      assert(lpredcost[col] < -tolerance);
      double newlb =
          double(ceil(gap / lpredcost[col] + localdomain.col_upper_[col] -
                      mipsolver.mipdata_->feastol));

      if (newlb <= localdomain.col_lower_[col]) continue;
      assert(newlb > localdomain.col_lower_[col]);

      if (mipsolver.mipdata_->domain.isBinary(col)) {
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
        lp.computeDualProof(mipsolver.mipdata_->domain,
                            mipsolver.mipdata_->upper_limit, inds, vals, rhs,
                            false)) {
      bool addedConstraints = false;

      HighsInt oldNumConflicts =
          mipsolver.mipdata_->conflictPool.getNumConflicts();
      for (const HighsDomainChange& domchg : boundChanges) {
        if (localdomain.isActive(domchg)) continue;
        localdomain.conflictAnalyzeReconvergence(
            domchg, inds.data(), vals.data(), inds.size(), rhs,
            mipsolver.mipdata_->conflictPool);
      }
      addedConstraints =
          mipsolver.mipdata_->conflictPool.getNumConflicts() != oldNumConflicts;

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

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    if (lpredcost[col] > mipsolver.mipdata_->feastol) {
      // col <= (cutoffbound - lpobj)/redcost + lb
      // so for lurkub = lb to ub - 1 we can compute the necessary cutoff
      // bound to reach this bound which is:
      //  lurkub = (cutoffbound - lpobj)/redcost + lb
      //  cutoffbound = (lurkub - lb) * redcost + lpobj
      HighsInt lb = (HighsInt)mipsolver.mipdata_->domain.col_lower_[col];
      HighsInt maxub;
      if (mipsolver.mipdata_->domain.col_upper_[col] == kHighsInf)
        maxub = lb + 1024;
      else
        maxub = (HighsInt)std::floor(
            mipsolver.mipdata_->domain.col_upper_[col] - 0.5);

      HighsInt step = 1;
      if (maxub - lb > 1024) step = (maxub - lb + 1023) >> 10;

      for (HighsInt lurkub = lb; lurkub <= maxub; lurkub += step) {
        double fracbound = (lurkub - lb + 1) - 10 * mipsolver.mipdata_->feastol;
        double requiredcutoffbound = fracbound * lpredcost[col] + lpobjective;
        bool useful = true;

        // check if we already have a better lurking bound stored
        auto pos = lurkingColUpper[col].lower_bound(requiredcutoffbound);
        for (auto it = pos; it != lurkingColUpper[col].end(); ++it) {
          if (it->second < lurkub + step) {
            useful = false;
            break;
          }
        }

        if (!useful) continue;

        // we have no better lurking bound stored store this lurking bound and
        // check if it dominates one that is already stored
        auto it =
            lurkingColUpper[col].emplace_hint(pos, requiredcutoffbound, lurkub);

        auto i = lurkingColUpper[col].begin();
        while (i != it) {
          if (i->second >= lurkub) {
            auto del = i++;
            lurkingColUpper[col].erase(del);
          } else {
            ++i;
          }
        }
      }
    } else if (lpredcost[col] < -mipsolver.mipdata_->feastol) {
      // col >= (cutoffbound - lpobj)/redcost + ub
      // so for lurklb = lb + 1 to ub we can compute the necessary cutoff
      // bound to reach this bound which is:
      //  lurklb = (cutoffbound - lpobj)/redcost + ub
      //  cutoffbound = (lurklb - ub) * redcost + lpobj

      HighsInt ub = (HighsInt)mipsolver.mipdata_->domain.col_upper_[col];
      HighsInt minlb;
      if (mipsolver.mipdata_->domain.col_lower_[col] == -kHighsInf)
        minlb = ub - 1024;
      else
        minlb = (HighsInt)(mipsolver.mipdata_->domain.col_lower_[col] + 1.5);

      HighsInt step = 1;
      if (ub - minlb > 1024) step = (ub - minlb + 1023) >> 10;

      for (HighsInt lurklb = minlb; lurklb <= ub; lurklb += step) {
        double fracbound = (lurklb - ub - 1) + 10 * mipsolver.mipdata_->feastol;
        double requiredcutoffbound = fracbound * lpredcost[col] + lpobjective -
                                     mipsolver.mipdata_->feastol;
        bool useful = true;

        // check if we already have a better lurking bound stored
        auto pos = lurkingColLower[col].lower_bound(requiredcutoffbound);
        for (auto it = pos; it != lurkingColLower[col].end(); ++it) {
          if (it->second > lurklb - step) {
            useful = false;
            break;
          }
        }

        if (!useful) continue;

        // we have no better lurking bound stored store this lurking bound and
        // check if it dominates one that is already stored
        auto it =
            lurkingColLower[col].emplace_hint(pos, requiredcutoffbound, lurklb);

        auto i = lurkingColLower[col].begin();
        while (i != it) {
          if (i->second <= lurklb) {
            auto del = i++;
            lurkingColLower[col].erase(del);
          } else {
            ++i;
          }
        }
      }
    }
  }
}
