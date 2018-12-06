/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2018 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HUtils.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHSUTILS_H_
#define LP_DATA_HIGHSUTILS_H_

#include "HConst.h"

class HighsUtils {
 public:

  // Logical check of double being +Infinity
  bool highs_isInfinity(double val);

};
#endif // LP_DATA_HIGHSUTILS_H_
