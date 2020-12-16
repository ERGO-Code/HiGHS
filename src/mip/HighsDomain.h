#ifndef HIGHS_DOMAIN_H_
#define HIGHS_DOMAIN_H_

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "mip/HighsDomainChange.h"
#include "mip/HighsMipSolver.h"
#include "util/HighsCDouble.h"

class HighsCutPool;

class HighsDomain {
  std::vector<uint8_t> changedcolsflags_;
  std::vector<int> changedcols_;

  std::vector<int> propRowNumChangedBounds_;

  std::vector<HighsDomainChange> domchgstack_;
  std::vector<int> domchgreason_;
  std::vector<double> prevboundval_;

  std::vector<HighsCDouble> activitymin_;
  std::vector<HighsCDouble> activitymax_;
  std::vector<int> activitymininf_;
  std::vector<int> activitymaxinf_;
  std::vector<uint8_t> propagateflags_;
  std::vector<int> propagateinds_;

  std::vector<HighsCDouble> activitycuts_;
  std::vector<int> activitycutsinf_;
  std::vector<unsigned> activitycutversion_;
  std::vector<uint8_t> propagatecutflags_;
  std::vector<int> propagatecutinds_;

  HighsMipSolver* mipsolver;
  HighsCutPool* cutpool;
  HighsDomain* parentdomain;

  int infeasible_ = 0;

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
  HighsDomain(HighsMipSolver& mipsolver, HighsCutPool& cutpool);

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
        activitycuts_(other.activitycuts_),
        activitycutsinf_(other.activitycutsinf_),
        activitycutversion_(other.activitycutversion_),
        propagatecutflags_(other.propagatecutflags_),
        propagatecutinds_(other.propagatecutinds_),
        mipsolver(other.mipsolver),
        cutpool(other.cutpool),
        parentdomain(other.parentdomain),
        infeasible_(other.infeasible_),
        colLower_(other.colLower_),
        colUpper_(other.colUpper_) {}

  HighsDomain& operator=(const HighsDomain&) = default;

  HighsDomain createChildDomain() {
    HighsDomain childdomain(*this);
    childdomain.parentdomain = this;
    return childdomain;
  }

  const std::vector<int>& getChangedCols() const { return changedcols_; }

  void clearChangedCols() {
    for (int i : changedcols_) changedcolsflags_[i] = 0;
    changedcols_.clear();
  }

  void clearChangedCols(int start) {
    int end = changedcols_.size();
    for (int i = start; i != end; ++i) changedcolsflags_[changedcols_[i]] = 0;

    changedcols_.resize(start);
  }

  const std::vector<int>& getPropagateRows() const { return propagateinds_; }

  void clearPropagateRows(int start) {
    int end = propagateinds_.size();
    for (int i = start; i != end; ++i) propagateflags_[propagateinds_[i]] = 0;

    propagateinds_.resize(start);
  }

  const std::vector<int>& getPropagateCuts() const { return propagatecutinds_; }

  void clearPropagateCuts(int start) {
    int end = propagatecutinds_.size();
    for (int i = start; i != end; ++i)
      propagatecutflags_[propagatecutinds_[i]] = 0;

    propagatecutinds_.resize(start);
  }

  void setParentDomain(HighsDomain* parentdomain) {
    this->parentdomain = parentdomain;
  }

  void cutAdded(int cut);

  void markPropagateCut(int cut);

  void markPropagate(int row);

  void computeRowActivities();

  bool infeasible() const { return infeasible_ != 0; }

  void changeBound(HighsDomainChange boundchg, int reason = -1);

  void changeBound(HighsBoundType boundtype, int col, double boundval,
                   int reason = -1) {
    changeBound({boundtype, col, boundval}, reason);
  }

  void fixCol(int col, double val) {
    assert(infeasible_ == 0);
    if (colLower_[col] < val)
      changeBound({HighsBoundType::Lower, col, val}, -2);

    if (infeasible_ == 0 && colUpper_[col] > val)
      changeBound({HighsBoundType::Upper, col, val}, -2);
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