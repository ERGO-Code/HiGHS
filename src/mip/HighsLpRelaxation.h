/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_LP_RELAXATION_H_
#define HIGHS_LP_RELAXATION_H_

#include <cstdint>
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

    int getRowLen(const HighsMipSolver& mipsolver) const;

    bool isIntegral(const HighsMipSolver& mipsolver) const;

    double getMaxAbsVal(const HighsMipSolver& mipsolver) const;

    static LpRow cut(int index) { return LpRow{kCutPool, index}; }
    static LpRow model(int index) { return LpRow{kModel, index}; }
  };

  const HighsMipSolver& mipsolver;
  Highs lpsolver;

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
  int64_t numlpiters;
  double avgSolveIters;
  int64_t numSolved;
  size_t epochs;
  size_t maxNumFractional;
  Status status;

  void storeDualInfProof();

  void storeDualUBProof();

  bool checkDualProof() const;

 public:
  HighsLpRelaxation(const HighsMipSolver& mip);

  HighsLpRelaxation(const HighsLpRelaxation& other);

  void loadModel();

  void getRow(int row, int& len, const int*& inds, const double*& vals) const {
    if (row < mipsolver.numRow())
      assert(lprows[row].origin == LpRow::Origin::kModel);
    else
      assert(lprows[row].origin == LpRow::Origin::kCutPool);
    lprows[row].get(mipsolver, len, inds, vals);
  }

  bool isRowIntegral(int row) const {
    assert(row < (int)lprows.size());
    return lprows[row].isIntegral(mipsolver);
  }

  double getAvgSolveIters() { return avgSolveIters; }

  int getRowLen(int row) const { return lprows[row].getRowLen(mipsolver); }

  double getMaxAbsRowVal(int row) const {
    return lprows[row].getMaxAbsVal(mipsolver);
  }

  const HighsLp& getLp() const { return lpsolver.getLp(); }

  const HighsSolution& getSolution() const { return lpsolver.getSolution(); }

  double slackUpper(int row) const;

  double slackLower(int row) const;

  double rowLower(int row) const { return lpsolver.getLp().rowLower_[row]; }

  double rowUpper(int row) const { return lpsolver.getLp().rowUpper_[row]; }

  double colLower(int col) const {
    return col < lpsolver.getLp().numCol_
               ? lpsolver.getLp().colLower_[col]
               : slackLower(col - lpsolver.getLp().numCol_);
  }

  double colUpper(int col) const {
    return col < lpsolver.getLp().numCol_
               ? lpsolver.getLp().colUpper_[col]
               : slackUpper(col - lpsolver.getLp().numCol_);
  }

  bool isColIntegral(int col) const {
    return col < lpsolver.getLp().numCol_
               ? mipsolver.variableType(col) != HighsVarType::CONTINUOUS
               : isRowIntegral(col - lpsolver.getLp().numCol_);
  }

  double solutionValue(int col) const {
    return col < lpsolver.getLp().numCol_
               ? getSolution().col_value[col]
               : getSolution().row_value[col - lpsolver.getLp().numCol_];
  }

  Status getStatus() const { return status; }

  int64_t getNumLpIterations() const { return numlpiters; }

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
    if (!currentbasisstored && lpsolver.getBasis().valid_) {
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

  int numRows() const { return lpsolver.getNumRows(); }

  int numCols() const { return lpsolver.getNumCols(); }

  void addCuts(HighsCutSet& cutset);

  void performAging();

  void removeObsoleteRows(bool notifyPool = true);

  void removeCuts(int ndelcuts, std::vector<int>& deletemask);

  void removeCuts();

  void flushDomain(HighsDomain& domain, bool continuous = false);

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
