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
/**@file lp_data/HighsHessianUtils.cpp
 * @brief
 */
#include "model/HighsHessianUtils.h"

#include <algorithm>

#include "util/HighsSort.h"
//#include <cassert>

#include "lp_data/HighsModelUtils.h"

HighsStatus assessHessian(HighsHessian& hessian, const HighsOptions& options,
                          const ObjSense sense) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Assess the Hessian dimensions and vector sizes, returning on error
  call_status = assessMatrixDimensions(options.log_options, "Hessian",
                                       hessian.dim_, hessian.q_start_,
                                       hessian.q_index_, hessian.q_value_);
  return_status =
      interpretCallStatus(call_status, return_status, "assessMatrixDimensions");
  if (return_status == HighsStatus::kError) return return_status;

  // If the Hessian has no columns there is nothing left to test
  if (hessian.dim_ == 0) return HighsStatus::kOk;

  // Assess the Hessian matrix
  //
  // The start of column 0 must be zero.
  if (hessian.q_start_[0]) {
    highsLogUser(options.log_options, HighsLogType::kError,
                 "Hessian has nonzero value (%" HIGHSINT_FORMAT
                 ") for the start of column 0\n",
                 hessian.q_start_[0]);
    return HighsStatus::kError;
  }
  // Assess G, deferring the assessment of values (other than those
  // which are identically zero)
  call_status = assessMatrix(options.log_options, "Hessian", hessian.dim_,
                             hessian.dim_, hessian.q_start_, hessian.q_index_,
                             hessian.q_value_, 0, kHighsInf);
  return_status =
      interpretCallStatus(call_status, return_status, "assessMatrix");
  if (return_status == HighsStatus::kError) return return_status;

  if (hessian.format_ == HessianFormat::kSquare) {
    // Form Q = (G+G^T)/2
    call_status = normaliseHessian(options, hessian);
    return_status =
        interpretCallStatus(call_status, return_status, "normaliseHessian");
    if (return_status == HighsStatus::kError) return return_status;
  }
  if (kHessianFormatInternal == HessianFormat::kTriangular) {
    // Extract the triangular part of Q: lower triangle column-wise
    // or, equivalently, upper triangle row-wise, ensuring that the
    // diagonal entry comes first, unless it's zero
    call_status = extractTriangularHessian(options, hessian);
    return_status = interpretCallStatus(call_status, return_status,
                                        "extractTriangularHessian");
    if (return_status == HighsStatus::kError) return return_status;
  }

  // Assess Q
  call_status =
      assessMatrix(options.log_options, "Hessian", hessian.dim_, hessian.dim_,
                   hessian.q_start_, hessian.q_index_, hessian.q_value_,
                   options.small_matrix_value, options.large_matrix_value);
  return_status =
      interpretCallStatus(call_status, return_status, "assessMatrix");
  if (return_status == HighsStatus::kError) return return_status;

  HighsInt hessian_num_nz = hessian.q_start_[hessian.dim_];
  // If the Hessian has nonzeros, check its diagonal entries in the
  // context of the objective sense. It's OK to be identically zero,
  // since it will be ignored.
  if (hessian_num_nz)
    if (!okHessianDiagonal(options, hessian, sense))
      return_status = HighsStatus::kError;

  // If entries have been removed from the matrix, resize the index
  // and value vectors
  if ((HighsInt)hessian.q_index_.size() > hessian_num_nz)
    hessian.q_index_.resize(hessian_num_nz);
  if ((HighsInt)hessian.q_value_.size() > hessian_num_nz)
    hessian.q_value_.resize(hessian_num_nz);

  if (return_status != HighsStatus::kError) return_status = HighsStatus::kOk;
  if (return_status != HighsStatus::kOk)
    highsLogDev(options.log_options, HighsLogType::kInfo,
                "assessHessian returns HighsStatus = %s\n",
                HighsStatusToString(return_status).c_str());
  return return_status;
}

bool okHessianDiagonal(const HighsOptions& options, HighsHessian& hessian,
                       const ObjSense sense) {
  const double kSmallHessianDiagonalValue = options.small_matrix_value;
  double min_illegal_diagonal_value = kHighsInf;
  double max_illegal_diagonal_value = -kHighsInf;
  const HighsInt dim = hessian.dim_;
  const HighsInt sense_sign = (HighsInt)sense;
  HighsInt num_small_diagonal_value = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    double diagonal_value = 0;
    // Assumes that the diagonal entry is always first, unless it's
    // zero so doesn't appear
    HighsInt iEl = hessian.q_start_[iCol];
    if (hessian.q_index_[iEl] == iCol)
      diagonal_value = sense_sign * hessian.q_value_[iEl];
    if (diagonal_value <= kSmallHessianDiagonalValue) {
      min_illegal_diagonal_value =
          std::min(diagonal_value, min_illegal_diagonal_value);
      max_illegal_diagonal_value =
          std::max(diagonal_value, max_illegal_diagonal_value);
      num_small_diagonal_value++;
    }
  }

  if (num_small_diagonal_value) {
    if (sense == ObjSense::kMinimize) {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Hessian has %" HIGHSINT_FORMAT
                   " diagonal entries in [%g, %g] less than %g\n",
                   num_small_diagonal_value, min_illegal_diagonal_value,
                   max_illegal_diagonal_value, kSmallHessianDiagonalValue);
    } else {
      highsLogUser(options.log_options, HighsLogType::kError,
                   "Hessian has %" HIGHSINT_FORMAT
                   " diagonal entries in [%g, %g] greater than %g\n",
                   num_small_diagonal_value, -max_illegal_diagonal_value,
                   -min_illegal_diagonal_value, -kSmallHessianDiagonalValue);
    }
  }
  return num_small_diagonal_value == 0;
}

HighsStatus extractTriangularHessian(const HighsOptions& options,
                                     HighsHessian& hessian) {
  // Viewing the Hessian column-wise, remove any entries in the strict
  // upper triangle
  HighsStatus return_status = HighsStatus::kOk;
  const HighsInt dim = hessian.dim_;
  HighsInt nnz = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    double diagonal_value = 0;
    const HighsInt nnz0 = nnz;
    for (HighsInt iEl = hessian.q_start_[iCol];
         iEl < hessian.q_start_[iCol + 1]; iEl++) {
      HighsInt iRow = hessian.q_index_[iEl];
      if (iRow < iCol) continue;
      hessian.q_index_[nnz] = iRow;
      hessian.q_value_[nnz] = hessian.q_value_[iEl];
      if (iRow == iCol && nnz > nnz0) {
        // Diagonal entry is not first in column so swap it in
        hessian.q_index_[nnz] = hessian.q_index_[nnz0];
        hessian.q_value_[nnz] = hessian.q_value_[nnz0];
        hessian.q_index_[nnz0] = iRow;
        hessian.q_value_[nnz0] = hessian.q_value_[iEl];
      }
      nnz++;
    }
    hessian.q_start_[iCol] = nnz0;
  }
  const HighsInt num_ignored_nz = hessian.q_start_[dim] - nnz;
  assert(num_ignored_nz >= 0);
  if (hessian.format_ == HessianFormat::kTriangular && num_ignored_nz) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Ignored %" HIGHSINT_FORMAT
                 " entries of Hessian in opposite triangle\n",
                 num_ignored_nz);
    hessian.q_start_[dim] = nnz;
    return_status = HighsStatus::kWarning;
  }
  hessian.format_ = HessianFormat::kTriangular;
  return return_status;
}

void triangularToSquareHessian(const HighsHessian& hessian,
                               vector<HighsInt>& start, vector<HighsInt>& index,
                               vector<double>& value) {
  const HighsInt dim = hessian.dim_;
  if (dim <= 0) {
    start.assign(1, 0);
    return;
  }
  assert(hessian.format_ == HessianFormat::kTriangular);
  const HighsInt nnz = hessian.q_start_[dim];
  const HighsInt square_nnz = nnz + (nnz - dim);
  start.resize(dim + 1);
  index.resize(square_nnz);
  value.resize(square_nnz);
  vector<HighsInt> length;
  length.assign(dim, 0);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt iRow = hessian.q_index_[hessian.q_start_[iCol]];
    assert(iRow == iCol);
    length[iCol]++;
    for (HighsInt iEl = hessian.q_start_[iCol] + 1;
         iEl < hessian.q_start_[iCol + 1]; iEl++) {
      HighsInt iRow = hessian.q_index_[iEl];
      assert(iRow > iCol);
      length[iRow]++;
      length[iCol]++;
    }
  }
  start[0] = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    start[iCol + 1] = start[iCol] + length[iCol];
  assert(square_nnz == start[dim]);
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt iEl = hessian.q_start_[iCol];
    HighsInt iRow = hessian.q_index_[iEl];
    HighsInt toEl = start[iCol];
    index[toEl] = iRow;
    value[toEl] = hessian.q_value_[iEl];
    start[iCol]++;
    for (HighsInt iEl = hessian.q_start_[iCol] + 1;
         iEl < hessian.q_start_[iCol + 1]; iEl++) {
      HighsInt iRow = hessian.q_index_[iEl];
      HighsInt toEl = start[iRow];
      index[toEl] = iCol;
      value[toEl] = hessian.q_value_[iEl];
      start[iRow]++;
      toEl = start[iCol];
      index[toEl] = iRow;
      value[toEl] = hessian.q_value_[iEl];
      start[iCol]++;
    }
  }
  start[0] = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    start[iCol + 1] = start[iCol] + length[iCol];
}

HighsStatus normaliseHessian(const HighsOptions& options,
                             HighsHessian& hessian) {
  // Normalise the Hessian to be (Q + Q^T)/2, where Q is the matrix
  // supplied. This guarantees that what's used internally is
  // symmetric.
  //
  // So someone preferring to supply only the upper triangle would
  // have to double its values..
  HighsStatus return_status = HighsStatus::kOk;
  const HighsInt dim = hessian.dim_;
  const HighsInt hessian_num_nz = hessian.q_start_[dim];
  if (hessian_num_nz <= 0) return HighsStatus::kOk;
  bool warning_found = false;

  HighsHessian transpose;
  transpose.dim_ = dim;
  transpose.q_start_.resize(dim + 1);
  transpose.q_index_.resize(hessian_num_nz);
  transpose.q_value_.resize(hessian_num_nz);
  // Form transpose of Hessian
  vector<HighsInt> qr_length;
  qr_length.assign(dim, 0);
  for (HighsInt iEl = 0; iEl < hessian_num_nz; iEl++)
    qr_length[hessian.q_index_[iEl]]++;

  transpose.q_start_[0] = 0;
  for (HighsInt iRow = 0; iRow < dim; iRow++)
    transpose.q_start_[iRow + 1] = transpose.q_start_[iRow] + qr_length[iRow];
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    for (HighsInt iEl = hessian.q_start_[iCol];
         iEl < hessian.q_start_[iCol + 1]; iEl++) {
      HighsInt iRow = hessian.q_index_[iEl];
      HighsInt iRowEl = transpose.q_start_[iRow];
      transpose.q_index_[iRowEl] = iCol;
      transpose.q_value_[iRowEl] = hessian.q_value_[iEl];
      transpose.q_start_[iRow]++;
    }
  }

  transpose.q_start_[0] = 0;
  for (HighsInt iRow = 0; iRow < dim; iRow++)
    transpose.q_start_[iRow + 1] = transpose.q_start_[iRow] + qr_length[iRow];

  HighsHessian normalised;
  HighsInt normalised_num_nz = 0;
  HighsInt normalised_size = hessian_num_nz;
  normalised.dim_ = dim;
  normalised.q_start_.resize(dim + 1);
  normalised.q_index_.resize(normalised_size);
  normalised.q_value_.resize(normalised_size);
  vector<double> column_value;
  vector<HighsInt> column_index;
  column_index.resize(dim);
  column_value.assign(dim, 0.0);
  const double small_matrix_value = 0;
  HighsInt num_small_values = 0;
  double max_small_value = 0;
  double min_small_value = kHighsInf;
  normalised.q_start_[0] = 0;
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    HighsInt column_num_nz = 0;
    for (HighsInt iEl = hessian.q_start_[iCol];
         iEl < hessian.q_start_[iCol + 1]; iEl++) {
      HighsInt iRow = hessian.q_index_[iEl];
      column_value[iRow] = hessian.q_value_[iEl];
      column_index[column_num_nz] = iRow;
      column_num_nz++;
    }
    for (HighsInt iEl = transpose.q_start_[iCol];
         iEl < transpose.q_start_[iCol + 1]; iEl++) {
      HighsInt iRow = transpose.q_index_[iEl];
      if (column_value[iRow]) {
        column_value[iRow] += transpose.q_value_[iEl];
      } else {
        column_value[iRow] = transpose.q_value_[iEl];
        column_index[column_num_nz] = iRow;
        column_num_nz++;
      }
    }
    if (normalised_num_nz + column_num_nz > normalised_size) {
      normalised_size =
          std::max(normalised_num_nz + column_num_nz, 2 * normalised_size);
      normalised.q_index_.resize(normalised_size);
      normalised.q_value_.resize(normalised_size);
    }
    // Halve the values, zeroing and accounting for any small ones
    for (HighsInt ix = 0; ix < column_num_nz; ix++) {
      HighsInt iRow = column_index[ix];
      double value = 0.5 * column_value[iRow];
      double abs_value = std::fabs(value);
      bool ok_value = abs_value > small_matrix_value;
      if (!ok_value) {
        value = 0;
        if (max_small_value < abs_value) max_small_value = abs_value;
        if (min_small_value > abs_value) min_small_value = abs_value;
        num_small_values++;
      }
      column_value[iRow] = value;
    }
    // Decide whether to exploit sparsity in extracting the indices
    // and values of nonzeros
    const HighsInt kDimTolerance = 10;
    const double kDensityTolerance = 0.1;
    const double density = (1.0 * column_num_nz) / (1.0 * dim);
    HighsInt to_ix = dim;
    const bool exploit_sparsity =
        dim > kDimTolerance && density < kDensityTolerance;
    if (exploit_sparsity) {
      // Exploit sparsity
      to_ix = column_num_nz;
      sortSetData(column_num_nz, &column_index[0], NULL, NULL);
    } else {
      to_ix = dim;
    }
    for (HighsInt ix = 0; ix < to_ix; ix++) {
      HighsInt iRow;
      if (exploit_sparsity) {
        iRow = column_index[ix];
      } else {
        iRow = ix;
      }
      double value = column_value[iRow];
      if (value) {
        normalised.q_index_[normalised_num_nz] = iRow;
        normalised.q_value_[normalised_num_nz] = value;
        normalised_num_nz++;
        column_value[iRow] = 0;
      }
    }
    for (HighsInt iRow = 0; iRow < dim; iRow++) assert(column_value[iRow] == 0);
    normalised.q_start_[iCol + 1] = normalised_num_nz;
  }
  if (num_small_values) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Normalised Hessian contains %" HIGHSINT_FORMAT
                 " |values| in [%g, %g] "
                 "less than %g: ignored\n",
                 num_small_values, min_small_value, max_small_value,
                 small_matrix_value);
    warning_found = true;
  }
  // Replace the Hessian by the normalised form
  hessian = normalised;
  if (warning_found)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}

void reportHessian(const HighsLogOptions& log_options, const HighsInt dim,
                   const HighsInt num_nz, const HighsInt* start,
                   const HighsInt* index, const double* value) {
  if (dim <= 0) return;
  highsLogUser(log_options, HighsLogType::kInfo,
               "Hessian Index              Value\n");
  for (HighsInt col = 0; col < dim; col++) {
    highsLogUser(log_options, HighsLogType::kInfo,
                 "    %8" HIGHSINT_FORMAT " Start   %10" HIGHSINT_FORMAT "\n",
                 col, start[col]);
    HighsInt to_el = (col < dim - 1 ? start[col + 1] : num_nz);
    for (HighsInt el = start[col]; el < to_el; el++)
      highsLogUser(log_options, HighsLogType::kInfo,
                   "          %8" HIGHSINT_FORMAT " %12g\n", index[el],
                   value[el]);
  }
  highsLogUser(log_options, HighsLogType::kInfo,
               "             Start   %10" HIGHSINT_FORMAT "\n", num_nz);
}
