/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_IMPLICATIONS_H_
#define HIGHS_IMPLICATIONS_H_

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "mip/HighsDomain.h"
#include "mip/HighsDomainChange.h"
#include "util/HighsHashTree.h"

class HighsCliqueTable;
class HighsLpRelaxation;

class HighsImplications {
  HighsInt nextCleanupCall;

  struct Implics {
    std::vector<HighsDomainChange> implics;
    /* The "tentative" implications:
      A implication of type x_j \ge (\ell^1_j - \ell^0_j) x_k + \ell^0_j is called "tentative", if
      (1) c_j = 0
      (2) x_j is fixed by applying dual fixing in probing
      These implications can only be used to perform globally valid reductions.
      Therefore, special treatment is required.
    */
    std::vector<HighsDomainChange> implics_tentative;
    bool computed = false;
  };
  std::vector<Implics> implications;
  int64_t numImplications;
  int64_t numVarBounds;
  int64_t maxVarBounds;

  bool computeImplications(HighsInt col, bool val);

 public:
  struct VarBound {
    double coef;
    double constant;

    double minValue() const {
      return static_cast<double>(static_cast<HighsCDouble>(constant) +
                                 std::min(coef, 0.0));
    }
    double maxValue() const {
      return static_cast<double>(static_cast<HighsCDouble>(constant) +
                                 std::max(coef, 0.0));
    }
  };

 private:
  std::vector<HighsHashTree<HighsInt, VarBound>> vubs;
  std::vector<HighsHashTree<HighsInt, VarBound>> vlbs;

 public:
  const HighsMipSolver& mipsolver;
  std::vector<HighsSubstitution> substitutions;
  std::vector<HighsBool> colsubstituted;

  // vector used to derive global reductions from dfprobing
  std::vector<HighsInt> binaryInvolvedInds_;
  enum binaryFixType {
    kNoReduction          = 0b0000,
    kGlobalLower          = 0b1010,
    kGlobalUpper          = 0b0101,
    kSubstituteComplement = 0b1001,
    kSubstituteEqual      = 0b0110,
  };
  /*
      Possible values for binaryInvolvedFlags_
       0 (0000, kNoReduction): Not involved
       2 (0010): fixed to 0 in second side probing
       1 (0001): fixed to 1 in second side probing
       8 (1000): fixed to 0 in first side probing
       4 (0100): fixed to 1 in first side probing
       10(1010, kGlobalLower): fixed to 0 in both side probing (global fixing!)
       5 (0101, kGlobalUpper): fixed to 1 in both side probing (global fixing!)
       9 (1001, kSubstituteComplement): substitutation type 1 --- x1 + x2 = 1
       6 (0110, kSubstituteEqual): substitutation type 2 --- x1 = x2
  */
  std::vector<HighsInt> binaryInvolvedFlags_;

  HighsImplications(const HighsMipSolver& mipsolver) : mipsolver(mipsolver) {
    HighsInt numcol = mipsolver.numCol();
    implications.resize(2 * static_cast<size_t>(numcol));
    colsubstituted.resize(numcol);
    vubs.resize(numcol);
    vlbs.resize(numcol);
    nextCleanupCall = mipsolver.numNonzero();
    numImplications = 0;
    numVarBounds = 0;
    maxVarBounds = calcMaxVarBounds(numcol);

    binaryInvolvedInds_.reserve(numcol);
    binaryInvolvedFlags_.assign(numcol, 0b0000);
  }

  std::function<void(HighsInt, HighsInt, HighsInt, double)>
      storeLiftingOpportunity;

  void reset() {
    colsubstituted.clear();
    colsubstituted.shrink_to_fit();
    implications.clear();
    implications.shrink_to_fit();

    HighsInt numcol = mipsolver.numCol();
    implications.resize(2 * static_cast<size_t>(numcol));
    colsubstituted.resize(numcol);
    numImplications = 0;
    vubs.clear();
    vubs.shrink_to_fit();
    vubs.resize(numcol);
    vlbs.clear();
    vlbs.shrink_to_fit();
    vlbs.resize(numcol);
    numVarBounds = 0;
    maxVarBounds = calcMaxVarBounds(numcol);

    nextCleanupCall = mipsolver.numNonzero();
    binaryInvolvedInds_.reserve(numcol);
    binaryInvolvedFlags_.assign(numcol, 0b0000);

  }

  constexpr static int64_t calcMaxVarBounds(HighsInt numcol) {
    return int64_t{5000000} + 10 * static_cast<int64_t>(numcol);
  };

  HighsInt getNumImplications() const {
    return static_cast<HighsInt>(numImplications);
  }

  const std::vector<HighsDomainChange>& getImplications(HighsInt col, bool val,
                                                        bool& infeasible) {
    HighsInt loc = 2 * col + val;
    if (!implications[loc].computed)
      infeasible = computeImplications(col, val);
    else
      infeasible = false;

    assert(implications[loc].computed);

    return implications[loc].implics;
  }

  const std::vector<HighsDomainChange>& getImplications_tentative(HighsInt col, bool val) {
    HighsInt loc = 2 * col + val;
    return implications[loc].implics_tentative;
  }


  bool implicationsCached(HighsInt col, bool val) {
    HighsInt loc = 2 * col + val;
    return implications[loc].computed;
  }

  bool tooManyVarBounds() const { return numVarBounds >= maxVarBounds; }

  void strengthenVarBound(VarBound& vbnd, HighsInt multiplier) const;

  void addVUB(HighsInt col, HighsInt vubcol, double vubcoef,
              double vubconstant);

  void addVUB(HighsInt col, HighsInt vubcol, double vubcoef, double vubconstant,
              double colupperbound, bool colisinteger);

  void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef,
              double vlbconstant);

  void addVLB(HighsInt col, HighsInt vlbcol, double vlbcoef, double vlbconstant,
              double collowerbound, bool colisinteger);

  void columnTransformed(HighsInt col, double scale, double constant) {
    // Update variable bounds affected by transformation
    if (scale < 0) std::swap(vubs[col], vlbs[col]);

    auto transformVbd = [&](HighsInt, VarBound& vbd) {
      vbd.constant -= constant;
      vbd.constant /= scale;
      vbd.coef /= scale;
    };

    vlbs[col].for_each(transformVbd);
    vubs[col].for_each(transformVbd);

    // Update substitutions affected by transformation
    for (auto& substitution : substitutions) {
      if (substitution.substcol == col) {
        substitution.offset -= constant;
        substitution.offset /= scale;
        substitution.scale /= scale;
      }
    }
  }

  std::pair<HighsInt, VarBound> getBestVub(HighsInt col,
                                           const HighsSolution& lpSolution,
                                           double& bestUb,
                                           const HighsDomain& globaldom) const;

  std::pair<HighsInt, VarBound> getBestVlb(HighsInt col,
                                           const HighsSolution& lpSolution,
                                           double& bestLb,
                                           const HighsDomain& globaldom) const;

  bool runProbing(HighsInt col, HighsInt& numReductions);

  void rebuild(HighsInt ncols, const std::vector<HighsInt>& cIndex,
               const std::vector<HighsInt>& rIndex);

  void buildFrom(const HighsImplications& init);

  void separateImpliedBounds(const HighsLpRelaxation& lpRelaxation,
                             const std::vector<double>& sol,
                             HighsCutPool& cutpool, double feastol,
                             HighsDomain& globaldom, bool thread_safe);

  void cleanupVarbounds(HighsInt col);

  void cleanupVlb(HighsInt col, HighsInt vlbCol,
                  HighsImplications::VarBound& vlb, double lb, bool& redundant,
                  bool& infeasible, bool allowBoundChanges = true) const;

  void cleanupVub(HighsInt col, HighsInt vubCol,
                  HighsImplications::VarBound& vub, double ub, bool& redundant,
                  bool& infeasible, bool allowBoundChanges = true) const;

  void applyImplications(HighsDomain& domain, HighsInt col, HighsInt val);

  // collect tentative binary implications
  void recordTentativeCliques(bool val, const HighsDomainChange& bchg) {
    const int iCol = bchg.column;
    if (val == 0) { // probing x_k = 0
      if (bchg.boundtype == HighsBoundType::kLower) { // fixed to 1
        if (!isFixedTo1(val, iCol)) {
          if (binaryInvolvedFlags_[iCol] == 0)
            binaryInvolvedInds_.push_back(iCol);
          binaryInvolvedFlags_[iCol] += 0b0001; // 0001
        }
      }
      else { // fixed to 0
        if (!isFixedTo0(val, iCol)) {
          if (binaryInvolvedFlags_[iCol] == 0)
            binaryInvolvedInds_.push_back(iCol);
          binaryInvolvedFlags_[iCol] += 0b0010; // 0010
        }
      }
    }
    else { // probing x_k = 1
      if (bchg.boundtype == HighsBoundType::kLower) { // fixed to 1
        if (!isFixedTo1(val, iCol)) {
          if (binaryInvolvedFlags_[iCol] == 0)
            binaryInvolvedInds_.push_back(iCol);
          binaryInvolvedFlags_[iCol] += 0b0100; // 0100
        }
      }
      else { // fixed to 0
        if (!isFixedTo0(val, iCol)) {
          if (binaryInvolvedFlags_[iCol] == 0)
            binaryInvolvedInds_.push_back(iCol);
          binaryInvolvedFlags_[iCol] += 0b1000; // 1000
        }
      }
    }
  }
  // clear tentative binary implications
  void clearTentativeClique() {
    for (auto iCol : binaryInvolvedInds_)
      binaryInvolvedFlags_[iCol] = binaryFixType::kNoReduction;
    binaryInvolvedInds_.clear();
  }
  // tools for recordTentativeCliques
  bool isFixedTo0(bool val, HighsInt iCol) {
    if (binaryInvolvedFlags_[iCol] == 0)
      return false;

    uint8_t mask;
    if (val == 0) { // probing at x = 0, last two digits
      mask = 1 << (1);
      return (binaryInvolvedFlags_[iCol] & mask) != 0;
    }
    else { // probing at x = 1, first two digits
      mask = 1 << (3);
      return (binaryInvolvedFlags_[iCol] & mask) != 0;
    }
  }
  // tools for recordTentativeCliques
  bool isFixedTo1(bool val, HighsInt iCol) {
    if (binaryInvolvedFlags_[iCol] == 0)
      return false;

    uint8_t mask;
    if (val == 0) { // probing at x = 0, last two digits
      mask = 1;
      return (binaryInvolvedFlags_[iCol] & mask) != 0;
    }
    else { // probint at x = 1, first two digits
      mask = 1 << (2);
      return (binaryInvolvedFlags_[iCol] & mask) != 0;
    }
  }

};

#endif
