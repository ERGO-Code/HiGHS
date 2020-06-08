/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsLp.h"

bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution) {
  if (solution.col_value.size() == (size_t)lp.numCol_ ||
      solution.col_dual.size() == (size_t)lp.numCol_ ||
      solution.row_value.size() == (size_t)lp.numRow_ ||
      solution.row_dual.size() == (size_t)lp.numRow_)
    return true;
  return false;
}
bool isBasisConsistent(const HighsLp& lp, const HighsBasis& basis) {
  if (basis.col_status.size() == (size_t)lp.numCol_ ||
      basis.row_status.size() == (size_t)lp.numRow_)
    return true;
  return false;
}

void clearSolutionUtil(HighsSolution& solution) {
  solution.col_dual.clear();
  solution.col_value.clear();
  solution.row_dual.clear();
  solution.row_value.clear();
}

void clearBasisUtil(HighsBasis& basis) {
  basis.row_status.clear();
  basis.col_status.clear();
  basis.valid_ = false;
}

void clearLp(HighsLp& lp) {
  lp.Astart_.clear();
  lp.Aindex_.clear();
  lp.Avalue_.clear();

  lp.col_names_.clear();
  lp.row_names_.clear();

  lp.sense_ = ObjSense::MINIMIZE;

  lp.colCost_.clear();
  lp.colLower_.clear();
  lp.colUpper_.clear();

  lp.integrality_.clear();
}

bool equalButForNames(const HighsLp& lp) {
  if (lp.numCol_ != lp.numCol_ || lp.numRow_ != lp.numRow_ ||
      lp.sense_ != lp.sense_ || lp.offset_ != lp.offset_ ||
      lp.model_name_ != lp.model_name_)
    return false;

  if (lp.colCost_ != lp.colCost_) return false;

  if (lp.colUpper_ != lp.colUpper_ || lp.colLower_ != lp.colLower_ ||
      lp.rowUpper_ != lp.rowUpper_ || lp.rowLower_ != lp.rowLower_)
    return false;

  if (lp.Astart_ != lp.Astart_ || lp.Aindex_ != lp.Aindex_ ||
      lp.Avalue_ != lp.Avalue_)
    return false;

  return true;
}

bool HighsLp::operator==(const HighsLp& lp) {
  if (!equalButForNames(lp)) return false;
  if (this->row_names_ != lp.row_names_ || this->col_names_ != lp.col_names_)
    return false;
  return true;
}