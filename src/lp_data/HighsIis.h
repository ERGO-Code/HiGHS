/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2024 by Julian Hall, Ivet Galabova,    */
/*    Leona Gottwald and Michael Feldmeier                               */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsIis.h
 * @brief Class-independent utilities for HiGHS
 */
#ifndef LP_DATA_HIGHSIIS_H_
#define LP_DATA_HIGHSIIS_H_

#include "lp_data/HighsLp.h"

struct HighsIis {
  bool valid = false;
  std::vector<HighsInt> col_index;
  std::vector<HighsInt> row_index;
  std::vector<HighsInt> col_bound;
  std::vector<HighsInt> row_bound;
  void invalidate();
};

HighsStatus getIisData(const HighsLp& lp,
                       const std::vector<double>& dual_ray_value,
                       HighsIis& iis);

#endif  // LP_DATA_HIGHSIIS_H_
