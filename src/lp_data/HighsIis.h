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

class HighsIis {
 public:
  HighsIis() {}

  void invalidate();
  void removeCddCol(const HighsInt cdd_col);
  void removeCddRow(const HighsInt cdd_row);
  HighsStatus getData(const HighsLp& lp, const HighsOptions& options,
		      const std::vector<double>& dual_ray_value);

  HighsStatus compute(const HighsLp& lp, const HighsOptions& options);

  bool inconsistentBounds(const HighsLp& lp, const HighsOptions& options);

  // Data members
  bool valid_ = false;
  HighsInt strategy_ = kIisStrategyMin;
  std::vector<HighsInt> col_index_;
  std::vector<HighsInt> row_index_;
  std::vector<HighsInt> col_bound_;
  std::vector<HighsInt> row_bound_;
};

#endif  // LP_DATA_HIGHSIIS_H_
