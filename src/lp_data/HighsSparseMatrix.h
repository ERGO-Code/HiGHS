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
/**@file lp_data/HighsSparseMatrix.h
 * @brief
 */
#ifndef LP_DATA_HIGHS_SPARSE_MATRIX_H_
#define LP_DATA_HIGHS_SPARSE_MATRIX_H_

#include <vector>

#include "lp_data/HConst.h"
#include "lp_data/HighsStatus.h"
#include "simplex/SimplexStruct.h"  //For SimplexScale until scaling is HighsScale
#include "util/HighsUtils.h"

class HighsSparseMatrix;

class HighsSparseMatrix {
 public:
  MatrixFormat format_ = MatrixFormat::kNone;
  HighsInt num_col_ = 0;
  HighsInt num_row_ = 0;
  std::vector<HighsInt> start_;
  std::vector<HighsInt> index_;
  std::vector<double> value_;

  bool operator==(const HighsSparseMatrix& matrix) const;
  void clear();
  void range(double& min_value, double& max_value) const;
  HighsStatus setFormat(
      const MatrixFormat desired_format = MatrixFormat::kColwise);
  void ensureColWise();
  void ensureRowWise();

  HighsStatus addCols(const HighsInt num_new_col, const HighsInt num_new_nz,
                      const HighsInt* new_matrix_start,
                      const HighsInt* new_matrix_index,
                      const double* new_matrix_value);
  HighsStatus addRows(const HighsInt num_new_row, const HighsInt num_new_nz,
                      const HighsInt* new_matrix_start,
                      const HighsInt* new_matrix_index,
                      const double* new_matrix_value);
  HighsStatus deleteCols(const HighsLogOptions& log_options,
                         const HighsIndexCollection& index_collection);
  HighsStatus deleteRows(const HighsLogOptions& log_options,
                         const HighsIndexCollection& index_collection);
  HighsStatus assessDimensions(const HighsLogOptions& log_options,
                               const std::string matrix_name);
  HighsStatus assess(const HighsLogOptions& log_options,
                     const std::string matrix_name,
                     const double small_matrix_value,
                     const double large_matrix_value);
  void scaleCol(const HighsInt col, const double colScale);
  void scaleRow(const HighsInt row, const double rowScale);
  void applyScale(const SimplexScale& scale);
  void unapplyScale(const SimplexScale& scale);
};

#endif
