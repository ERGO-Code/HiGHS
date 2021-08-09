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

#include "lp_data/HStruct.h"
#include "lp_data/HighsSparseMatrix.h"
#include "simplex/SimplexStruct.h"  // for SimplexScale (temporary)

class HighsLp;

class HighsLp {
 public:
  HighsLp() { clear(); }
  // Model data
  HighsInt num_col_;
  HighsInt num_row_;

  std::vector<double> col_cost_;
  std::vector<double> col_lower_;
  std::vector<double> col_upper_;
  std::vector<double> row_lower_;
  std::vector<double> row_upper_;

  HighsSparseMatrix a_matrix_;

  ObjSense sense_;
  double offset_;

  std::string model_name_;

  std::vector<std::string> col_names_;
  std::vector<std::string> row_names_;

  std::vector<HighsVarType> integrality_;

  HighsScale scale_;
  bool is_scaled_;

  bool operator==(const HighsLp& lp);
  bool equalButForNames(const HighsLp& lp) const;
  bool isMip() const;
  double objectiveValue(const std::vector<double>& solution) const;
  void setMatrixDimensions();
  bool dimensionsOk(std::string message) const;
  bool equal(const SimplexScale& scale) const;
  void scaleClear();
  void applyScale(const SimplexScale& scale);
  void unapplyScale(const SimplexScale& scale);
  void applyScale();
  void unapplyScale();
  void clear();
};

#endif
