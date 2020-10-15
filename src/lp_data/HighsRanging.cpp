/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsRanging.cpp
 * @brief Compute LP ranging data for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsRanging.h"

#include <algorithm>
#include <cassert>

HighsStatus getHighsRanging(HighsRanging& ranging,
                            const HighsModelObject& highs_model_object) {
  if (highs_model_object.scaled_model_status_ != HighsModelStatus::OPTIMAL) {
    highsLogMessage(highs_model_object.options_.logfile, HighsMessageType::ERROR,
		    "Cannot get ranging without an optimal solution");
    return HighsStatus::Error;
  }
  if (!highs_model_object.simplex_lp_status_.valid) {
    highsLogMessage(highs_model_object.options_.logfile, HighsMessageType::ERROR,
		    "Cannot get ranging without a valid Simplex LP");
    return HighsStatus::Error;
  }
  const HighsSimplexInfo& simplex_info = highs_model_object.simplex_info_;
  const HighsSimplexBasis& simplex_basis = highs_model_object.simplex_basis_;
  const vector<int>& Bvalue_ = simplex_basis.BasicIndex_;
  vector<double> xi = Bvalue_;
  for (int i = 0; i < numRow; i++) {
    xi[i] = max(xi[i], Blower_[i]);
    xi[i] = min(xi[i], Bupper_[i]);
  }

  return HighsStatus::OK;
}
