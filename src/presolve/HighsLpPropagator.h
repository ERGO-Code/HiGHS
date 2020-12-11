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
  std::vector<int> activitymininf_;
  std::vector<int> activitymaxinf_;
  std::vector<uint8_t> propagateflags_;
  std::vector<int> propagateinds_;

  std::vector<double>& Avalue_;
  std::vector<int>& Aindex_;
  std::vector<int>& Astart_;
  std::vector<int>& Aend_;

  std::vector<double>& ARvalue_;
  std::vector<int>& ARindex_;
  std::vector<int>& ARstart_;

  const std::vector<int>& flagRow;
  const std::vector<int>& flagCol;
  std::vector<double>& rowLower_;
  std::vector<double>& rowUpper_;
  const std::vector<HighsVarType>& integrality_;

  int infeasible_ = 0;
  int numBoundChgs_ = 0;

  void computeMinActivity(int start, int end, const int* ARindex,
                          const double* ARvalue, int& ninfmin,
                          HighsCDouble& activitymin);

  void computeMaxActivity(int start, int end, const int* ARindex,
                          const double* ARvalue, int& ninfmax,
                          HighsCDouble& activitymax);

  int propagateRowUpper(const int* Rindex, const double* Rvalue, int Rlen,
                        double Rupper, const HighsCDouble& minactivity,
                        int ninfmin, HighsDomainChange* boundchgs);

  int propagateRowLower(const int* Rindex, const double* Rvalue, int Rlen,
                        double Rlower, const HighsCDouble& maxactivity,
                        int ninfmax, HighsDomainChange* boundchgs);

  void updateActivityLbChange(int col, double oldbound, double newbound);

  void updateActivityUbChange(int col, double oldbound, double newbound);

  double doChangeBound(const HighsDomainChange& boundchg);

 public:
  std::vector<double> colLower_;
  std::vector<double> colUpper_;

  HighsLpPropagator(const std::vector<double>& colLower,
                    const std::vector<double>& colUpper,
                    const std::vector<HighsVarType>& integrality_,
                    std::vector<double>& Avalue_, std::vector<int>& Aindex_,
                    std::vector<int>& Astart_, std::vector<int>& Aend_,
                    std::vector<double>& ARvalue_, std::vector<int>& ARindex_,
                    std::vector<int>& ARstart_, const std::vector<int>& flagRow,
                    const std::vector<int>& flagCol,
                    std::vector<double>& rowLower_,
                    std::vector<double>& rowUpper_);

  void markPropagate(int row);

  void computeRowActivities();

  bool infeasible() const { return infeasible_ != 0; }

  void changeBound(HighsDomainChange boundchg);

  int propagate();

  int tightenCoefficients();

  int getNumChangedBounds() const { return numBoundChgs_; }
};

#endif