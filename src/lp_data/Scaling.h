/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/Scaling.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_SCALING_H_
#define LP_DATA_SCALING_H_

#include "HighsModelObject.h"
#include <cassert>
#include <vector>

void scaleHighsModel(HighsModelObject highs_model) {
  highs_model.lp_scaled_ = highs_model.lp_;
}
#endif
