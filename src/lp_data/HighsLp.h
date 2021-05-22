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
  HighsInt numCol_ = 0;
  HighsInt numRow_ = 0;

  std::vector<HighsInt> Astart_;
  std::vector<HighsInt> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;

  MatrixOrientation orientation_ = MatrixOrientation::kNone;
  ObjSense sense_ = ObjSense::kMinimize;
  double offset_ = 0;

  std::string model_name_ = "";
  std::string lp_name_ = "";

  std::vector<std::string> col_names_;
  std::vector<std::string> row_names_;

  std::vector<HighsVarType> integrality_;

  bool operator==(const HighsLp& lp);
  bool equalButForNames(const HighsLp& lp);
  bool isMip();
  void clear();
};

#endif
