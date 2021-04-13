/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_CUTPOOL_H_
#define HIGHS_CUTPOOL_H_

#include <memory>
#include <unordered_map>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsDomain.h"
#include "mip/HighsDynamicRowMatrix.h"

class HighsLpRelaxation;

struct HighsCutSet {
  std::vector<HighsInt> cutindices;
  std::vector<HighsInt> ARstart_;
  std::vector<HighsInt> ARindex_;
  std::vector<double> ARvalue_;
  std::vector<double> lower_;
  std::vector<double> upper_;

  HighsInt numCuts() const { return cutindices.size(); }

  void resize(HighsInt nnz) {
    HighsInt ncuts = numCuts();
    lower_.resize(ncuts, -kHighsInf);
    upper_.resize(ncuts);
    ARstart_.resize(ncuts + 1);
    ARindex_.resize(nnz);
    ARvalue_.resize(nnz);
  }

  void clear() {
    cutindices.clear();
    upper_.clear();
    ARstart_.clear();
    ARindex_.clear();
    ARvalue_.clear();
  }

  bool empty() const { return cutindices.empty(); }
};

class HighsCutPool {
 private:
  HighsDynamicRowMatrix matrix_;
  std::vector<double> rhs_;
  std::vector<unsigned> modification_;
  std::vector<int16_t> ages_;
  std::vector<double> rownormalization_;
  std::vector<double> maxabscoef_;
  std::vector<uint8_t> rowintegral;
  std::unordered_multimap<uint32_t, int> supportmap;
  std::vector<HighsDomain::CutpoolPropagation*> propagationDomains;

  double bestObservedScore;
  double minScoreFactor;

  HighsInt agelim_;
  HighsInt softlimit_;
  HighsInt numLpCuts;
  std::vector<HighsInt> ageDistribution;

  bool isDuplicate(size_t hash, double norm, HighsInt* Rindex, double* Rvalue,
                   HighsInt Rlen, double rhs);

 public:
  HighsCutPool(HighsInt ncols, HighsInt agelim, HighsInt softlimit)
      : matrix_(ncols), agelim_(agelim), softlimit_(softlimit), numLpCuts(0) {
    ageDistribution.resize(agelim_ + 1);
    minScoreFactor = 0.9;
    bestObservedScore = 0.0;
  }
  const HighsDynamicRowMatrix& getMatrix() const { return matrix_; }

  const std::vector<double>& getRhs() const { return rhs_; }

  void resetAge(HighsInt cut) {
    if (ages_[cut] < 0)
      ages_[cut] = -1;
    else
      ages_[cut] = 0;
  }

  bool ageLpCut(HighsInt cut, HighsInt agelimit) {
    assert(ages_[cut] < 0);
    --ages_[cut];
    if (ages_[cut] < -agelimit) {
      ages_[cut] = 1;
      return true;
    }

    return false;
  }

  double getParallelism(HighsInt row1, HighsInt row2) const;

  void performAging();

  void lpCutRemoved(HighsInt cut);

  void addPropagationDomain(HighsDomain::CutpoolPropagation* domain) {
    propagationDomains.push_back(domain);
  }

  void removePropagationDomain(HighsDomain::CutpoolPropagation* domain) {
    for (HighsInt k = propagationDomains.size() - 1; k >= 0; --k) {
      if (propagationDomains[k] == domain) {
        propagationDomains.erase(propagationDomains.begin() + k);
        return;
      }
    }
  }

  void setAgeLimit(HighsInt agelim) {
    agelim_ = agelim;
    ageDistribution.resize(agelim_ + 1);
  }

  void separate(const std::vector<double>& sol, HighsDomain& domprop,
                HighsCutSet& cutset, double feastol);

  void separateLpCutsAfterRestart(HighsCutSet& cutset);

  bool cutIsIntegral(HighsInt cut) const { return rowintegral[cut]; }

  HighsInt getNumCuts() const {
    return matrix_.getNumRows() - matrix_.getNumDelRows();
  }

  double getMaxAbsCutCoef(HighsInt cut) const { return maxabscoef_[cut]; }

  HighsInt addCut(const HighsMipSolver& mipsolver, HighsInt* Rindex,
                  double* Rvalue, HighsInt Rlen, double rhs,
                  bool integral = false, bool extractCliques = true);

  HighsInt getRowLength(HighsInt row) const {
    return matrix_.getRowEnd(row) - matrix_.getRowStart(row);
  }

  unsigned getModificationCount(HighsInt cut) const {
    return modification_[cut];
  }

  void getCut(HighsInt cut, HighsInt& cutlen, const HighsInt*& cutinds,
              const double*& cutvals) const {
    HighsInt start = matrix_.getRowStart(cut);
    cutlen = matrix_.getRowEnd(cut) - start;
    cutinds = matrix_.getARindex() + start;
    cutvals = matrix_.getARvalue() + start;
  }
};

#endif
