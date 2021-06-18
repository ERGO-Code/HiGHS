/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/*    Authors: Julian Hall, Ivet Galabova, Qi Huangfu, Leona Gottwald    */
/*    and Michael Feldmeier                                              */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "mip/HighsLpAggregator.h"

#include "mip/HighsLpRelaxation.h"

HighsLpAggregator::HighsLpAggregator(const HighsLpRelaxation& lprelaxation)
    : lprelaxation(lprelaxation) {
  vectorsum.setDimension(lprelaxation.getLp().numRow_ +
                         lprelaxation.getLp().numCol_);
}

void HighsLpAggregator::addRow(HighsInt row, double weight) {
  HighsInt len;
  const double* vals;
  const HighsInt* inds;
  lprelaxation.getRow(row, len, inds, vals);

  for (HighsInt i = 0; i != len; ++i) vectorsum.add(inds[i], weight * vals[i]);

  vectorsum.add(lprelaxation.getLp().numCol_ + row, -weight);
}

void HighsLpAggregator::getCurrentAggregation(std::vector<HighsInt>& inds,
                                              std::vector<double>& vals,
                                              bool negate) {
  const double droptol =
      lprelaxation.getMipSolver().options_mip_->small_matrix_value;
  vectorsum.cleanup(
      [droptol](HighsInt col, double val) { return std::abs(val) <= droptol; });

  inds = vectorsum.getNonzeros();
  HighsInt len = inds.size();
  vals.resize(len);

  if (negate)
    for (HighsInt i = 0; i != len; ++i) vals[i] = -vectorsum.getValue(inds[i]);
  else
    for (HighsInt i = 0; i != len; ++i) vals[i] = vectorsum.getValue(inds[i]);
}

void HighsLpAggregator::clear() { vectorsum.clear(); }
