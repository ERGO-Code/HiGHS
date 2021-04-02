/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_LP_PROPAGATOR_H_
#define HIGHS_LP_PROPAGATOR_H_

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsDomainChange.h"
#include "util/HighsCDouble.h"

/// propagates domains as part of LP presolve
/// final propagated bounds are relaxed by a wide enough margin
/// so that they cannot be used in any basic feasible solution
class HighsLpPropagator {
  std::vector<HighsCDouble> activitymin_;
  std::vector<HighsCDouble> activitymax_;
  std::vector<HighsInt> activitymininf_;
  std::vector<HighsInt> activitymaxinf_;
  std::vector<uint8_t> propagateflags_;
  std::vector<HighsInt> propagateinds_;

  std::vector<double>& Avalue_;
  std::vector<HighsInt>& Aindex_;
  std::vector<HighsInt>& Astart_;
  std::vector<HighsInt>& Aend_;

  std::vector<double>& ARvalue_;
  std::vector<HighsInt>& ARindex_;
  std::vector<HighsInt>& ARstart_;

  const std::vector<HighsInt>& flagRow;
  const std::vector<HighsInt>& flagCol;
  std::vector<double>& rowLower_;
  std::vector<double>& rowUpper_;
  const std::vector<HighsVarType>& integrality_;

  HighsInt infeasible_ = 0;
  HighsInt numBoundChgs_ = 0;

  void computeMinActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const double* ARvalue, HighsInt& ninfmin,
                          HighsCDouble& activitymin);

  void computeMaxActivity(HighsInt start, HighsInt end, const HighsInt* ARindex,
                          const double* ARvalue, HighsInt& ninfmax,
                          HighsCDouble& activitymax);

  HighsInt propagateRowUpper(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double Rupper,
                             const HighsCDouble& minactivity, HighsInt ninfmin,
                             HighsDomainChange* boundchgs);

  HighsInt propagateRowLower(const HighsInt* Rindex, const double* Rvalue,
                             HighsInt Rlen, double Rlower,
                             const HighsCDouble& maxactivity, HighsInt ninfmax,
                             HighsDomainChange* boundchgs);

  void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

  void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

  double doChangeBound(const HighsDomainChange& boundchg);

 public:
  std::vector<double> colLower_;
  std::vector<double> colUpper_;

  HighsLpPropagator(
      const std::vector<double>& colLower, const std::vector<double>& colUpper,
      const std::vector<HighsVarType>& integrality_,
      std::vector<double>& Avalue_, std::vector<HighsInt>& Aindex_,
      std::vector<HighsInt>& Astart_, std::vector<HighsInt>& Aend_,
      std::vector<double>& ARvalue_, std::vector<HighsInt>& ARindex_,
      std::vector<HighsInt>& ARstart_, const std::vector<HighsInt>& flagRow,
      const std::vector<HighsInt>& flagCol, std::vector<double>& rowLower_,
      std::vector<double>& rowUpper_);

  void markPropagate(HighsInt row);

  void computeRowActivities();

  bool infeasible() const { return infeasible_ != 0; }

  void changeBound(HighsDomainChange boundchg);

  HighsInt propagate();

  HighsInt tightenCoefficients();

  HighsInt getNumChangedBounds() const { return numBoundChgs_; }
};

#endif
