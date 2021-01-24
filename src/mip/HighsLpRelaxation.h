#ifndef HIGHS_LP_RELAXATION_H_
#define HIGHS_LP_RELAXATION_H_

#include <memory>

#include "Highs.h"
#include "mip/HighsMipSolver.h"

class HighsDomain;
class HighsCutSet;
class HighsPseudocost;

class HighsLpRelaxation {
 public:
  enum class Status {
    NotSet,
    Optimal,
    Infeasible,
    UnscaledDualFeasible,
    UnscaledPrimalFeasible,
    UnscaledInfeasible,
    Error,
  };

 private:
  struct LpRow {
    enum Origin {
      kModel,
      kCutPool,
    };

    Origin origin;
    int index;

    void get(const HighsMipSolver& mipsolver, int& len, const int*& inds,
             const double*& vals) const;

    static LpRow cut(int index) { return LpRow{kCutPool, index}; }
    static LpRow model(int index) { return LpRow{kModel, index}; }
  };

  const HighsMipSolver& mipsolver;
  Highs lpsolver;
  std::vector<int> lp2cutpoolindex;

  std::vector<LpRow> lprows;

  std::vector<std::pair<int, double>> fractionalints;
  std::vector<double> dualproofvals;
  std::vector<int> dualproofinds;
  std::vector<double> dualproofbuffer;
  double dualproofrhs;
  bool hasdualproof;
  double objective;
  std::shared_ptr<const HighsBasis> basischeckpoint;
  bool currentbasisstored;
  size_t numlpiters;
  Status status;

  void storeDualInfProof();

  void storeDualUBProof();

  bool checkDualProof() const;

 public:
  HighsLpRelaxation(const HighsMipSolver& mip);

  HighsLpRelaxation(const HighsLpRelaxation& other);

  int getCutIndex(int lpindex) const {
    return lp2cutpoolindex[lpindex - mipsolver.numRow()];
  }

  Status getStatus() const { return status; }

  size_t getNumLpIterations() const { return numlpiters; }

  bool integerFeasible() const {
    if ((status == Status::Optimal ||
         status == Status::UnscaledPrimalFeasible) &&
        fractionalints.empty())
      return true;

    return false;
  }

  double computeBestEstimate(const HighsPseudocost& ps) const;

  static bool scaledOptimal(Status status) {
    switch (status) {
      case Status::Optimal:
      case Status::UnscaledDualFeasible:
      case Status::UnscaledPrimalFeasible:
      case Status::UnscaledInfeasible:
        return true;
      default:
        return false;
    }
  }

  static bool unscaledPrimalFeasible(Status status) {
    switch (status) {
      case Status::Optimal:
      case Status::UnscaledPrimalFeasible:
        return true;
      default:
        return false;
    }
  }

  static bool unscaledDualFeasible(Status status) {
    switch (status) {
      case Status::Optimal:
      case Status::UnscaledDualFeasible:
        return true;
      default:
        return false;
    }
  }

  void recoverBasis();

  void setObjectiveLimit(double objlim = HIGHS_CONST_INF) {
    // lpsolver.setHighsOptionValue("dual_objective_value_upper_bound", objlim);
  }

  void storeBasis() {
    if (!currentbasisstored) {
      basischeckpoint = std::make_shared<HighsBasis>(lpsolver.getBasis());
      currentbasisstored = true;
    }
  }

  std::shared_ptr<const HighsBasis> getStoredBasis() const {
    return basischeckpoint;
  }

  void setStoredBasis(std::shared_ptr<const HighsBasis> basis) {
    basischeckpoint = basis;
    currentbasisstored = false;
  }

  const HighsMipSolver& getMipSolver() const { return mipsolver; }

  int getNumModelRows() const { return mipsolver.numRow(); }

  int getNumLpRows() const { return lpsolver.getNumRows(); }

  void addCuts(HighsCutSet& cutset);

  void removeCuts(int ndelcuts, std::vector<int>& deletemask);

  void removeCuts();

  void flushDomain(HighsDomain& domain, bool continuous = false);

  double rowLower(int row) const { return lpsolver.getLp().rowLower_[row]; }

  double rowUpper(int row) const { return lpsolver.getLp().rowUpper_[row]; }

  void getDualProof(const int*& inds, const double*& vals, double& rhs,
                    int& len) {
    inds = dualproofinds.data();
    vals = dualproofvals.data();
    rhs = dualproofrhs;
    len = dualproofinds.size();
  }

  bool computeDualProof(const HighsDomain& globaldomain, double upperbound,
                        std::vector<int>& inds, std::vector<double>& vals,
                        double& rhs) const;

  bool computeDualInfProof(const HighsDomain& globaldomain,
                           std::vector<int>& inds, std::vector<double>& vals,
                           double& rhs);

  Status resolveLp(HighsDomain* domain = nullptr);

  Status run(bool resolve_on_error = true);

  Highs& getLpSolver() { return lpsolver; }
  const Highs& getLpSolver() const { return lpsolver; }

  const std::vector<std::pair<int, double>>& getFractionalIntegers() const {
    return fractionalints;
  }

  std::vector<std::pair<int, double>>& getFractionalIntegers() {
    return fractionalints;
  }

  double getObjective() const { return objective; }

  void setIterationLimit(int limit = HIGHS_CONST_I_INF) {
    lpsolver.setHighsOptionValue("simplex_iteration_limit", limit);
  }
};

#endif