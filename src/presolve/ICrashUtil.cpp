/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file presolve/ICrashUtil.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */

#include "lp_data/HighsLp.h"
#include "util/HighsUtils.h"

void MuptiplyByTranspose(const HighsLp& lp, const std::vector<double>& v, std::vector<double>& result) {
  assert(result.size() == lp.numCol_);
  assert(v.size() == lp.numRow_);

  result.assign(lp.numCol_, 0);
  for (int col = 0; col < lp.numCol_; col++) {
    for (int k = lp.Astart_[col]; k < lp.Astart_[col + 1]; k++) {
      const int row = lp.Aindex_[k];
      result.at(col) += lp.Avalue_[k] * lp.rowUpper_[row];
    }
  }
}

void printMinorIterationDetails(const double iteration, const double col,
                                const double old_value, const double update,
                                const double ctx, const std::vector<double>& r,
                                const double quadratic_objective) {
  double rnorm = getNorm2(r);
  std::cout << "iter " << iteration;
  std::cout << ", col " << col;
  std::cout << ", update " << update;
  std::cout << ", old_value " << old_value;
  std::cout << ", new_value " << old_value + update;
  std::cout << ", ctx " << ctx;
  std::cout << ", r " << rnorm;
  std::cout << ", quadratic_objective " << quadratic_objective;
  std::cout << std::endl;
}
