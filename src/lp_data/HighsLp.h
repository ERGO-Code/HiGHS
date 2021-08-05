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
/**@file lp_data/HighsLp.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_LP_H_
#define LP_DATA_HIGHS_LP_H_

#include <string>
#include <vector>

#include "lp_data/HConst.h"

class HighsLp;

class HighsLp {
 public:
  // Model data
  HighsInt num_col_ = 0;
  HighsInt num_row_ = 0;

  std::vector<HighsInt> a_start_;
  std::vector<HighsInt> a_index_;
  std::vector<double> a_value_;
  std::vector<double> col_cost_;
  std::vector<double> col_lower_;
  std::vector<double> col_upper_;
  std::vector<double> row_lower_;
  std::vector<double> row_upper_;

  MatrixFormat format_ = MatrixFormat::kNone;
  ObjSense sense_ = ObjSense::kMinimize;
  double offset_ = 0;

  std::string model_name_ = "";

  std::vector<std::string> col_names_;
  std::vector<std::string> row_names_;

  std::vector<HighsVarType> integrality_;

  bool operator==(const HighsLp& lp);
  bool equalButForNames(const HighsLp& lp) const;
  bool isMip() const;
  double objectiveValue(const std::vector<double>& solution) const;
  void clear();
};

#endif
