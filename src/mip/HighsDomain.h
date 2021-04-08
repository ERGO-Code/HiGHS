/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_DOMAIN_H_
#define HIGHS_DOMAIN_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <vector>

#include "mip/HighsDomainChange.h"
#include "mip/HighsMipSolver.h"
#include "util/HighsCDouble.h"

class HighsCutPool;

class HighsDomain {
 public:
  struct Reason {
    HighsInt type;
    HighsInt index;

    enum {
      kBranching = -1,
      kUnknown = -2,
      kModelRow = -3,
      kCliqueTable = -4,
    };
    static Reason branching() { return Reason{kBranching, 0}; }
    static Reason unspecified() { return Reason{kUnknown, 0}; }
    static Reason cliqueTable() { return Reason{kCliqueTable, 0}; }
    static Reason modelRow(HighsInt row) { return Reason{kModelRow, row}; }
    static Reason cut(HighsInt cutpool, HighsInt cut) {
      return Reason{cutpool, cut};
    }
  };

  struct CutpoolPropagation {
    HighsInt cutpoolindex;
    HighsDomain* domain;
    HighsCutPool* cutpool;
    std::vector<HighsCDouble> activitycuts_;
    std::vector<HighsInt> activitycutsinf_;
    std::vector<unsigned> activitycutversion_;
    std::vector<uint8_t> propagatecutflags_;
    std::vector<HighsInt> propagatecutinds_;

    CutpoolPropagation(HighsInt cutpoolindex, HighsDomain* domain,
                       HighsCutPool& cutpool);

    CutpoolPropagation(const CutpoolPropagation& other);

    ~CutpoolPropagation();

    void cutAdded(HighsInt cut);

    void markPropagateCut(HighsInt cut);

    void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

    void updateActivityUbChange(HighsInt col, double oldbound, double newbound);
  };

 private:
  std::vector<uint8_t> changedcolsflags_;
  std::vector<HighsInt> changedcols_;

  std::vector<HighsInt> propRowNumChangedBounds_;

  std::vector<HighsDomainChange> domchgstack_;
  std::vector<Reason> domchgreason_;
  std::vector<double> prevboundval_;

  std::vector<HighsCDouble> activitymin_;
  std::vector<HighsCDouble> activitymax_;
  std::vector<HighsInt> activitymininf_;
  std::vector<HighsInt> activitymaxinf_;
  std::vector<uint8_t> propagateflags_;
  std::vector<HighsInt> propagateinds_;

  HighsMipSolver* mipsolver;

 private:
  std::deque<CutpoolPropagation> cutpoolpropagation;

  bool infeasible_ = 0;
  Reason infeasible_reason;
  HighsInt infeasible_pos;

  void updateActivityLbChange(HighsInt col, double oldbound, double newbound);

  void updateActivityUbChange(HighsInt col, double oldbound, double newbound);

  double doChangeBound(const HighsDomainChange& boundchg);

 public:
  std::vector<double> colLower_;
  std::vector<double> colUpper_;

  HighsDomain(HighsMipSolver& mipsolver);

  HighsDomain(const HighsDomain& other)
      : changedcolsflags_(other.changedcolsflags_),
        changedcols_(other.changedcols_),
        domchgstack_(other.domchgstack_),
        domchgreason_(other.domchgreason_),
        prevboundval_(other.prevboundval_),
        activitymin_(other.activitymin_),
        activitymax_(other.activitymax_),
        activitymininf_(other.activitymininf_),
        activitymaxinf_(other.activitymaxinf_),
        propagateflags_(other.propagateflags_),
        propagateinds_(other.propagateinds_),
        mipsolver(other.mipsolver),
        cutpoolpropagation(other.cutpoolpropagation),
        infeasible_(other.infeasible_),
        infeasible_reason(other.infeasible_reason),
        colLower_(other.colLower_),
        colUpper_(other.colUpper_) {
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.domain = this;
  }

  HighsDomain& operator=(const HighsDomain& other) {
    changedcolsflags_ = other.changedcolsflags_;
    changedcols_ = other.changedcols_;
    domchgstack_ = other.domchgstack_;
    domchgreason_ = other.domchgreason_;
    prevboundval_ = other.prevboundval_;
    activitymin_ = other.activitymin_;
    activitymax_ = other.activitymax_;
    activitymininf_ = other.activitymininf_;
    activitymaxinf_ = other.activitymaxinf_;
    propagateflags_ = other.propagateflags_;
    propagateinds_ = other.propagateinds_;
    mipsolver = other.mipsolver;
    cutpoolpropagation = other.cutpoolpropagation;
    infeasible_ = other.infeasible_;
    infeasible_reason = other.infeasible_reason;
    colLower_ = other.colLower_;
    colUpper_ = other.colUpper_;
    for (CutpoolPropagation& cutpoolprop : cutpoolpropagation)
      cutpoolprop.domain = this;
    return *this;
  }

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

  const std::vector<HighsInt>& getChangedCols() const { return changedcols_; }

  void addCutpool(HighsCutPool& cutpool);

  void clearChangedCols() {
    for (HighsInt i : changedcols_) changedcolsflags_[i] = 0;
    changedcols_.clear();
  }

  void clearChangedCols(HighsInt start) {
    HighsInt end = changedcols_.size();
    for (HighsInt i = start; i != end; ++i)
      changedcolsflags_[changedcols_[i]] = 0;

    changedcols_.resize(start);
  }

  void markPropagate(HighsInt row);

  void markPropagateCut(Reason reason);

  void computeRowActivities();

  bool infeasible() const { return infeasible_; }

  void changeBound(HighsDomainChange boundchg,
                   Reason reason = Reason::branching());

  void changeBound(HighsBoundType boundtype, HighsInt col, double boundval,
                   Reason reason = Reason::branching()) {
    changeBound({boundtype, col, boundval}, reason);
  }

  void fixCol(HighsInt col, double val, Reason reason = Reason::unspecified()) {
    assert(infeasible_ == 0);
    if (colLower_[col] < val)
      changeBound({HighsBoundType::Lower, col, val}, reason);

    if (infeasible_ == 0 && colUpper_[col] > val)
      changeBound({HighsBoundType::Upper, col, val}, reason);
  }

  HighsDomainChange backtrack();

  const std::vector<HighsDomainChange>& getDomainChangeStack() const {
    return domchgstack_;
  }

  const std::vector<Reason>& getDomainChangeReason() const {
    return domchgreason_;
  }

  std::vector<HighsDomainChange> getReducedDomainChangeStack() const {
    std::vector<HighsDomainChange> reducedstack;
    reducedstack.reserve(domchgstack_.size());
    for (const HighsDomainChange& domchg : domchgstack_) {
      // keep only the tightest bound change for each variable
      if ((domchg.boundtype == HighsBoundType::Lower &&
           colLower_[domchg.column] != domchg.boundval) ||
          (domchg.boundtype == HighsBoundType::Upper &&
           colUpper_[domchg.column] != domchg.boundval))
        continue;

      reducedstack.push_back(domchg);
    }
    reducedstack.shrink_to_fit();
    return reducedstack;
  }

  void setDomainChangeStack(const std::vector<HighsDomainChange>& domchgstack);

  bool propagate();

  void tightenCoefficients(HighsInt* inds, double* vals, HighsInt len,
                           double& rhs) const;

  double getMinActivity(HighsInt row) const {
    return activitymininf_[row] == 0 ? double(activitymin_[row])
                                     : -HIGHS_CONST_INF;
  }

  double getMaxActivity(HighsInt row) const {
    return activitymaxinf_[row] == 0 ? double(activitymax_[row])
                                     : HIGHS_CONST_INF;
  }

  double getMinCutActivity(const HighsCutPool& cutpool, HighsInt cut);

  bool isBinary(HighsInt col) const {
    return mipsolver->variableType(col) != HighsVarType::CONTINUOUS &&
           colLower_[col] == 0.0 && colUpper_[col] == 1.0;
  }

  HighsVarType variableType(HighsInt col) const {
    return mipsolver->variableType(col);
  }

  bool isFixed(HighsInt col) const { return colLower_[col] == colUpper_[col]; }

  bool isFixing(const HighsDomainChange& domchg) const;
};

#endif
