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

//class HighsSparseMatrix;
//class HVector;

class HighsSparseMatrix {
 public:
  HighsSparseMatrix() { clear(); }
  MatrixFormat format_;
  HighsInt num_col_;
  HighsInt num_row_;
  std::vector<HighsInt> start_;
  std::vector<HighsInt> p_end_;
  std::vector<HighsInt> index_;
  std::vector<double> value_;

  bool operator==(const HighsSparseMatrix& matrix) const;
  void clear();
  bool isRowwise() const;
  bool isColwise() const;
  HighsInt num_nz() const;
  void range(double& min_value, double& max_value) const;
  HighsStatus setFormat(
      const MatrixFormat desired_format = MatrixFormat::kColwise);
  void ensureColWise();
  void ensureRowWise();

  HighsStatus addCols(const HighsInt num_new_col, const HighsInt num_new_nz,
                      const HighsInt* new_matrix_start,
                      const HighsInt* new_matrix_index,
                      const double* new_matrix_value,
                      const int8_t* in_partition = NULL);
  HighsStatus addRows(const HighsInt num_new_row, const HighsInt num_new_nz,
                      const HighsInt* new_matrix_start,
                      const HighsInt* new_matrix_index,
                      const double* new_matrix_value,
                      const int8_t* in_partition = NULL);
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
  // Methods for PRICE, including the creation and updating of the
  // partitioned row-wise matrix
  void createPartition(const HighsSparseMatrix& matrix,
                       const int8_t* in_partition = NULL);
  void priceByColumn(HVector& result, const HVector& vector) const;
  void priceByRow(HVector& result, const HVector& vector) const;
  void priceByRowWithSwitch(HVector& result, const HVector& vector,
                            const double expected_density,
                            const HighsInt from_row,
                            const double switch_density) const;
  void update(const HighsInt var_in, const HighsInt var_out);
  double computeDot(const HVector& vector, const HighsInt use_col) const;
  void collectAj(HVector& vector, const HighsInt use_col,
                 const double multiplier) const;

 private:
  void priceByRowDenseResult(HVector& result, const HVector& vector,
                             const HighsInt from_row);
};

#endif
