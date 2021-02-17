/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file mip/HighsLpAggregator.h
 * @brief Class to aggregate rows of the LP
 *
 *
 * @author Leona Gottwald
 */

#ifndef MIP_HIGHS_LP_AGGREGATOR_H_
#define MIP_HIGHS_LP_AGGREGATOR_H_

#include <cstdint>
#include <vector>

#include "mip/HighsSparseVectorSum.h"

class HighsLpRelaxation;

/// Helper class to compute single-row relaxations from the current LP
/// relaxation by substituting bounds and aggregating rows
class HighsLpAggregator {
 private:
  const HighsLpRelaxation& lprelaxation;

  HighsSparseVectorSum vectorsum;

 public:
  HighsLpAggregator(const HighsLpRelaxation& lprelaxation);

  /// add an LP row to the aggregation using the given weight
  void addRow(int row, double weight);

  /// returns the current aggregation of LP rows. The aggregation includes slack
  /// variables so that it is always an equation with right hand side 0.
  void getCurrentAggregation(std::vector<int>& inds, std::vector<double>& vals,
                             bool negate);

  /// clear the current aggregation
  void clear();
};

#endif