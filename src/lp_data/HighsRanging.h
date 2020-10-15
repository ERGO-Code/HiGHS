/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsRanging.h
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_RANGING_H_
#define LP_DATA_HIGHS_RANGING_H_

#include <vector>

#include "lp_data/HighsModelObject.h"

struct HighsRangingRecord {
  std::vector<double> Value_;
  std::vector<double> Objective_;
  std::vector<int> InCol_;
  std::vector<int> OutCol_;
};

struct HighsRanging {
  HighsRangingRecord colCostUp;
  HighsRangingRecord colCostDown;
  HighsRangingRecord rowBoundUp;
  HighsRangingRecord rowBoundDown;
};

HighsStatus getHighsRanging(HighsRanging& ranging,
                            const HighsModelObject& highs_model_object);
#endif
