#ifndef HIGHS_DOMAIN_H_
#define HIGHS_DOMAIN_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <unordered_map>
#include <vector>

#include "mip/HighsDomainChange.h"
#include "mip/HighsMipSolver.h"
#include "util/HighsCDouble.h"

class HighsCutPool;

class HighsDomain {
 public:
  struct Reason {
    int type;
    int index;

    enum {
      kBranching = -1,
      kUnknown = -2,
      kModelRow = -3,
    };
    static Reason branching() { return Reason{kBranching, 0}; }
    static Reason unspecified() { return Reason{kUnknown, 0}; }
    static Reason modelRow(int row) { return Reason{kModelRow, row}; }
    static Reason cut(int cutpool, int cut) { return Reason{cutpool, cut}; }
  };

  struct CutpoolPropagation {
    int cutpoolindex;
    HighsDomain* domain;
    HighsCutPool* cutpool;
    std::vector<HighsCDouble> activitycuts_;
    std::vector<int> activitycutsinf_;
    std::vector<unsigned> activitycutversion_;
    std::vector<uint8_t> propagatecutflags_;
    std::vector<int> propagatecutinds_;

    CutpoolPropagation(int cutpoolindex, HighsDomain* domain,
                       HighsCutPool& cutpool);

    CutpoolPropagation(const CutpoolPropagation& other);

    ~CutpoolPropagation();

    void cutAdded(int cut);

    void markPropagateCut(int cut);

    void updateActivityLbChange(int col, double oldbound, double newbound);

    void updateActivityUbChange(int col, double oldbound, double newbound);
  };

 private:
  std::vector<uint8_t> changedcolsflags_;
  std::vector<int> changedcols_;

  std::vector<int> propRowNumChangedBounds_;

  std::vector<HighsDomainChange> domchgstack_;
  std::vector<Reason> domchgreason_;
  std::vector<double> prevboundval_;

  std::vector<HighsCDouble> activitymin_;
  std::vector<HighsCDouble> activitymax_;
  std::vector<int> activitymininf_;
  std::vector<int> activitymaxinf_;
  std::vector<uint8_t> propagateflags_;
  std::vector<int> propagateinds_;

  HighsMipSolver* mipsolver;

 private:
  std::deque<CutpoolPropagation> cutpoolpropagation;

  bool infeasible_ = 0;
  Reason infeasible_reason;

  void updateActivityLbChange(int col, double oldbound, double newbound);

  void updateActivityUbChange(int col, double oldbound, double newbound);

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

  const std::vector<int>& getChangedCols() const { return changedcols_; }

  void addCutpool(HighsCutPool& cutpool);

  void clearChangedCols() {
    for (int i : changedcols_) changedcolsflags_[i] = 0;
    changedcols_.clear();
  }

  void clearChangedCols(int start) {
    int end = changedcols_.size();
    for (int i = start; i != end; ++i) changedcolsflags_[changedcols_[i]] = 0;

    changedcols_.resize(start);
  }

  void markPropagate(int row);

  void markPropagateCut(Reason reason);

  void computeRowActivities();

  bool infeasible() const { return infeasible_; }

  void changeBound(HighsDomainChange boundchg,
                   Reason reason = Reason::branching());

  void changeBound(HighsBoundType boundtype, int col, double boundval,
                   Reason reason = Reason::branching()) {
    changeBound({boundtype, col, boundval}, reason);
  }

  void fixCol(int col, double val, Reason reason = Reason::unspecified()) {
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

  void propagate();

  void tightenCoefficients(int* inds, double* vals, int len, double& rhs) const;

  bool isBinary(int col) const {
    return mipsolver->variableType(col) != HighsVarType::CONTINUOUS &&
           colLower_[col] == 0.0 && colUpper_[col] == 1.0;
  }

  bool isFixed(int col) const { return colLower_[col] == colUpper_[col]; }

  bool isFixing(const HighsDomainChange& domchg) const;
};

#endif