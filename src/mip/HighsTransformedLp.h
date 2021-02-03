/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsTransformedLp.h
 * @brief LP transformations useful for cutting plane separation. This includes
 * bound substitution with simple and variable bounds, handling of slack
 * variables, flipping the complementation of integers.
 * @author Leona Gottwald
 */

#ifndef MIP_HIGHS_TRANSFORMED_LP_H_
#define MIP_HIGHS_TRANSFORMED_LP_H_

#include <cstdint>
#include <vector>

#include "lp_data/HConst.h"
#include "mip/HighsImplications.h"
#include "mip/HighsSparseVectorSum.h"
#include "util/HighsCDouble.h"

class HighsLpRelaxation;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsTransformedLp {
 private:
  const HighsLpRelaxation& lprelaxation;
  HighsImplications& implications;

  std::vector<const std::pair<const int, HighsImplications::VarBound>*> bestVub;
  std::vector<const std::pair<const int, HighsImplications::VarBound>*> bestVlb;
  std::vector<double> simpleLbDist;
  std::vector<double> simpleUbDist;
  std::vector<double> lbDist;
  std::vector<double> ubDist;
  std::vector<double> boundDist;
  enum class BoundType : uint8_t {
    kSimpleUb,
    kSimpleLb,
    kVariableUb,
    kVariableLb,
  };
  std::vector<BoundType> boundTypes;
  HighsSparseVectorSum vectorsum;

 public:
  HighsTransformedLp(const HighsLpRelaxation& lprelaxation,
                     HighsImplications& implications);

  double boundDistance(int col) const { return boundDist[col]; }

  bool transform(std::vector<double>& vals, std::vector<double>& upper,
                 std::vector<double>& solval, std::vector<int>& inds,
                 double& rhs, bool integralPositive, bool preferVbds = false);

  bool untransform(std::vector<double>& vals, std::vector<int>& inds,
                   double& rhs, bool integral = false);
};

#endif