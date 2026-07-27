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

#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsImplications.h"
#include "util/HighsCDouble.h"
#include "util/HighsInt.h"
#include "util/HighsSparseVectorSum.h"

class HighsLpRelaxation;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsTransformedLp {
 public:
  enum class BoundType : uint8_t {
    kSimpleUb,
    kSimpleLb,
    kVariableUb,
    kVariableLb,
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
  std::vector<BoundType> boundTypes;
  HighsSparseVectorSum vectorsum;

 public:
  HighsTransformedLp(const HighsLpRelaxation& lprelaxation,
                     HighsImplications& implications,
                     const HighsDomain& globaldom);

  double boundDistance(HighsInt col) const { return boundDist[col]; }

  bool transform(std::vector<double>& vals, std::vector<double>& upper,
                 std::vector<double>& solval, std::vector<HighsInt>& inds,
                 double& rhs, bool& integralPositive, bool preferVbds = false);

  bool untransform(std::vector<double>& vals, std::vector<HighsInt>& inds,
                   double& rhs, bool integral = false);

  // The bound substitution that the most recent call to transform() used for
  // the given column. A cut that is assembled from several transformed base
  // rows must be untransformed with the same substitutions that were used when
  // its coefficients were computed, so callers that keep a cut across further
  // transform() calls need to record these and restore them before calling
  // untransform().
  BoundType boundType(HighsInt col) const { return boundTypes[col]; }

  void setBoundType(HighsInt col, BoundType boundType) {
    boundTypes[col] = boundType;
  }

  const HighsDomain& getGlobaldom() const { return globaldom_; }
};

#endif
