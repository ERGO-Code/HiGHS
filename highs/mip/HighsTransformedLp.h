/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsTransformedLp.h
 * @brief LP transformations useful for cutting plane separation. This includes
 * bound substitution with simple and variable bounds, handling of slack
 * variables, flipping the complementation of integers.
 */

#ifndef MIP_HIGHS_TRANSFORMED_LP_H_
#define MIP_HIGHS_TRANSFORMED_LP_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsImplications.h"
#include "util/HighsCDouble.h"
#include "util/HighsSparseVectorSum.h"
#include "util/HighsType.h"

class HighsLpRelaxation;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsTransformedLp {
 public:
  // sparse row
  struct Row {
    std::vector<HighsInt> inds;
    std::vector<double> vals;
    HighsCDouble rhs = 0.0;

    Row() = default;
    Row(std::vector<HighsInt> inds_, std::vector<double> vals_,
        HighsCDouble rhs_)
        : inds(std::move(inds_)), vals(std::move(vals_)), rhs(rhs_) {}

    // scales the right hand side and all coefficients
    void scale(double factor) {
      rhs *= factor;
      for (double& val : vals) val *= factor;
    }

    // removes zero coefficients
    void removeZeros() {
      // replace zeros by the last element to keep the removal cheap
      for (HighsInt i = static_cast<HighsInt>(vals.size()) - 1; i >= 0; --i) {
        if (vals[i] == 0.0) {
          inds[i] = inds.back();
          vals[i] = vals.back();
          inds.pop_back();
          vals.pop_back();
        }
      }
    }

    // returns the maximal absolute coefficient
    double maxAbsVal() const {
      double maxVal = 0.0;
      for (double val : vals) maxVal = std::max(std::abs(val), maxVal);
      return maxVal;
    }

    // returns the minimal absolute coefficient
    double minAbsVal() const {
      double minVal = kHighsInf;
      for (double val : vals) minVal = std::min(std::abs(val), minVal);
      return minVal;
    }

    // returns the violation of the row by the given solution
    HighsCDouble violation(const std::vector<double>& sol) const {
      HighsCDouble viol = -rhs;
      for (HighsInt i = 0; i != static_cast<HighsInt>(inds.size()); ++i)
        viol += static_cast<HighsCDouble>(sol[inds[i]]) * vals[i];
      return viol;
    }
  };

  // variable of a row in the transformed space
  struct TransformedVar {
    HighsInt col;
    double val;
    double upper;
    // solution value in the current complementation
    double solval;
    // solution value after flipping the complementation; swapping both
    // solution values restores them exactly when flipping back
    double flippedSolval;
    bool complemented;
    // whether the variable is integral; cut generation relaxes integral
    // variables with small coefficients
    bool integral;
  };

  // row in the transformed space
  struct TransformedRow {
    std::vector<TransformedVar> vars;
    HighsCDouble rhs;

    // scales the right hand side and all coefficients
    void scale(double factor) {
      rhs *= factor;
      for (TransformedVar& var : vars) var.val *= factor;
    }

    // removes variables with zero coefficients
    void removeZeros() {
      // replace zeros by the last element to keep the removal cheap
      for (HighsInt i = static_cast<HighsInt>(vars.size()) - 1; i >= 0; --i) {
        if (vars[i].val == 0.0) {
          vars[i] = vars.back();
          vars.pop_back();
        }
      }
    }

    // returns the maximal absolute coefficient
    double maxAbsVal() const {
      double maxVal = 0.0;
      for (const TransformedVar& var : vars)
        maxVal = std::max(std::abs(var.val), maxVal);
      return maxVal;
    }

    // flips the complementation of the given variable of this row
    void flipComplementation(TransformedVar& var) {
      // only variables with finite upper bounds can be complemented
      assert(var.upper != kHighsInf);

      // flip complementation
      var.complemented = !var.complemented;
      std::swap(var.solval, var.flippedSolval);
      rhs -= static_cast<HighsCDouble>(var.upper) * var.val;
      var.val = -var.val;
    }

    // removes the complementation of all variables
    void removeComplementation() {
      for (TransformedVar& var : vars) {
        if (var.complemented) flipComplementation(var);
      }
    }

    // copies the given indices, values, and right hand side
    void load(const std::vector<HighsInt>& inds,
              const std::vector<double>& vals, HighsCDouble rhs_) {
      const HighsInt numVars = inds.size();
      vars.resize(numVars);
      for (HighsInt i = 0; i != numVars; ++i) {
        vars[i].col = inds[i];
        vars[i].val = vals[i];
      }
      rhs = rhs_;
    }

    // copies the indices, values, and right hand side to the given row
    void store(Row& row) const {
      const HighsInt numVars = vars.size();
      row.inds.resize(numVars);
      row.vals.resize(numVars);
      for (HighsInt i = 0; i != numVars; ++i) {
        row.inds[i] = vars[i].col;
        row.vals[i] = vars[i].val;
      }
      row.rhs = rhs;
    }
  };

 private:
  const HighsLpRelaxation& lprelaxation;
  const HighsDomain& globaldom_;

  std::vector<std::pair<HighsInt, HighsImplications::VarBound>> bestVub;
  std::vector<std::pair<HighsInt, HighsImplications::VarBound>> bestVlb;
  std::vector<double> simpleLbDist;
  std::vector<double> simpleUbDist;
  std::vector<double> lbDist;
  std::vector<double> ubDist;
  std::vector<double> boundDist;
  enum class BoundType : uint8_t {
    kUnused,
    kSimpleUb,
    kSimpleLb,
    kVariableUb,
    kVariableLb,
  };
  std::vector<BoundType> boundTypes;
  HighsSparseVectorSum vectorsum;

 public:
  HighsTransformedLp(const HighsLpRelaxation& lprelaxation,
                     HighsImplications& implications,
                     const HighsDomain& globaldom);

  double boundDistance(HighsInt col) const { return boundDist[col]; }

  // transforms the given row into the transformed space
  bool transform(const std::vector<HighsInt>& inds,
                 const std::vector<double>& vals, HighsCDouble rhs,
                 TransformedRow& transformedRow, bool& integralPositive,
                 bool preferVbds = false, bool enforceSameBds = false);

  // transforms the given row back into the original space
  bool untransform(Row& row, bool integral = false);

  const HighsDomain& getGlobaldom() const { return globaldom_; }

  void initMultiRowTransform() {
    boundTypes.assign(boundTypes.size(), BoundType::kUnused);
  }
};

#endif
