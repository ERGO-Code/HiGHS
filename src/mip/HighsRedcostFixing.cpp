/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include "mip/HighsRedcostFixing.h"

#include "mip/HighsMipSolverData.h"

void HighsRedcostFixing::propagateRootRedcost(const HighsMipSolver& mipsolver) {
  if (lurkingColLower.empty()) return;

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    for (auto it =
             lurkingColLower[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColLower[col].end(); ++it) {
      if (it->second > mipsolver.mipdata_->domain.colLower_[col]) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::Lower, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }
    }

    for (auto it =
             lurkingColUpper[col].lower_bound(mipsolver.mipdata_->upper_limit);
         it != lurkingColUpper[col].end(); ++it) {
      if (it->second < mipsolver.mipdata_->domain.colUpper_[col]) {
        mipsolver.mipdata_->domain.changeBound(
            HighsBoundType::Upper, col, (double)it->second,
            HighsDomain::Reason::unspecified());
        if (mipsolver.mipdata_->domain.infeasible()) return;
      }
    }
  }

  mipsolver.mipdata_->domain.propagate();
}

void HighsRedcostFixing::propagateRedCost(const HighsMipSolver& mipsolver,
                                          HighsDomain& localdomain,
                                          const std::vector<double>& lpredcost,
                                          double lpobjective) {
  HighsCDouble gap =
      HighsCDouble(mipsolver.mipdata_->upper_limit) - lpobjective;
  double tolerance = 10 * mipsolver.mipdata_->feastol;
  assert(!localdomain.infeasible());
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
    if (localdomain.colUpper_[col] == localdomain.colLower_[col]) continue;

    double threshold = double((HighsCDouble(localdomain.colUpper_[col]) -
                               localdomain.colLower_[col]) *
                                  gap +
                              tolerance);

    if ((localdomain.colUpper_[col] == HIGHS_CONST_INF &&
         lpredcost[col] > tolerance) ||
        lpredcost[col] > threshold) {
      assert(localdomain.colLower_[col] != -HIGHS_CONST_INF);
      assert(lpredcost[col] > tolerance);
      double newub =
          double(floor(gap / lpredcost[col] + localdomain.colLower_[col] +
                       mipsolver.mipdata_->feastol));
      if (newub >= localdomain.colUpper_[col]) continue;
      assert(newub < localdomain.colUpper_[col]);

      localdomain.changeBound(HighsBoundType::Upper, col, newub,
                              HighsDomain::Reason::unspecified());
      if (localdomain.infeasible()) return;
    } else if ((localdomain.colLower_[col] == -HIGHS_CONST_INF &&
                lpredcost[col] < -tolerance) ||
               lpredcost[col] < -threshold) {
      assert(localdomain.colUpper_[col] != HIGHS_CONST_INF);
      assert(lpredcost[col] < -tolerance);
      double newlb =
          double(ceil(gap / lpredcost[col] + localdomain.colUpper_[col] -
                      mipsolver.mipdata_->feastol));

      if (newlb <= localdomain.colLower_[col]) continue;
      assert(newlb > localdomain.colLower_[col]);

      localdomain.changeBound(HighsBoundType::Lower, col, newlb,
                              HighsDomain::Reason::unspecified());
      if (localdomain.infeasible()) return;
    }
  }

  localdomain.propagate();
}

void HighsRedcostFixing::addRootRedcost(const HighsMipSolver& mipsolver,
                                        const std::vector<double>& lpredcost,
                                        double lpobjective) {
  lurkingColLower.resize(mipsolver.numCol());
  lurkingColUpper.resize(mipsolver.numCol());

  for (HighsInt col : mipsolver.mipdata_->integral_cols) {
    if (lpredcost[col] > mipsolver.mipdata_->feastol) {
      // col <= (cutoffbound - lpobj)/redcost + lb
      // so for lurkub = lb to ub - 1 we can compute the necessary cutoff bound
      // to reach this bound which is:
      //  lurkub = (cutoffbound - lpobj)/redcost + lb
      //  cutoffbound = (lurkub - lb) * redcost + lpobj
      HighsInt lb = (HighsInt)mipsolver.mipdata_->domain.colLower_[col];
      HighsInt maxub;
      if (mipsolver.mipdata_->domain.colUpper_[col] == HIGHS_CONST_INF)
        maxub = lb + 10;
      else
        maxub = (HighsInt)(mipsolver.mipdata_->domain.colUpper_[col] - 0.5);

      HighsInt step = 1;
      if (maxub - lb > 1000) step = (maxub - lb + 999) / 1000;

      for (HighsInt lurkub = lb; lurkub <= maxub; lurkub += step) {
        double fracbound = (lurkub - lb + 1) - 10 * mipsolver.mipdata_->feastol;
        double requiredcutoffbound = fracbound * lpredcost[col] + lpobjective;
        bool useful = true;

        // check if we already have a better lurking bound stored
        auto pos = lurkingColUpper[col].lower_bound(requiredcutoffbound);
        for (auto it = pos; it != lurkingColUpper[col].end(); ++it) {
          if (it->second <= lurkub) {
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
      // so for lurklb = lb + 1 to ub we can compute the necessary cutoff bound
      // to reach this bound which is:
      //  lurklb = (cutoffbound - lpobj)/redcost + ub
      //  cutoffbound = (lurklb - ub) * redcost + lpobj

      HighsInt ub = (HighsInt)mipsolver.mipdata_->domain.colUpper_[col];
      HighsInt minlb;
      if (mipsolver.mipdata_->domain.colLower_[col] == -HIGHS_CONST_INF)
        minlb = ub - 10;
      else
        minlb = (HighsInt)(mipsolver.mipdata_->domain.colLower_[col] + 1.5);

      HighsInt step = 1;
      if (ub - minlb > 1000) step = (ub - minlb + 999) / 1000;

      for (HighsInt lurklb = minlb; lurklb <= ub; lurklb += step) {
        double fracbound = (lurklb - ub - 1) + 10 * mipsolver.mipdata_->feastol;
        double requiredcutoffbound = fracbound * lpredcost[col] + lpobjective -
                                     mipsolver.mipdata_->feastol;
        bool useful = true;

        // check if we already have a better lurking bound stored
        auto pos = lurkingColLower[col].lower_bound(requiredcutoffbound);
        for (auto it = pos; it != lurkingColLower[col].end(); ++it) {
          if (it->second >= lurklb) {
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
