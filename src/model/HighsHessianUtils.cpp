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
//#include <cassert>

#include "lp_data/HighsModelUtils.h"

HighsStatus assessHessian(HighsHessian& hessian, const HighsOptions& options) {
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  // Assess the Hessian dimensions and vector sizes, returning on error
  call_status = assessMatrixDimensions(options.log_options, "Hessian",
				       hessian.dim_,
				       hessian.q_start_, hessian.q_index_, hessian.q_value_);
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
  call_status = assessMatrix(options.log_options, "Hessian",
			     hessian.dim_, hessian.dim_,
			     hessian.q_start_, hessian.q_index_, hessian.q_value_,
			     0, kHighsInf);
  return_status =
    interpretCallStatus(call_status, return_status, "assessMatrix");
  if (return_status == HighsStatus::kError) return return_status;

  call_status = normaliseHessian(options, hessian);
  return_status =
    interpretCallStatus(call_status, return_status, "normaliseHessian");
  if (return_status == HighsStatus::kError) return return_status;

  HighsInt hessian_num_nz = hessian.q_start_[hessian.dim_];
  // If entries have been removed from the matrix, resize the index
  // and value vectors
  if ((HighsInt)hessian.q_index_.size() > hessian_num_nz) hessian.q_index_.resize(hessian_num_nz);
  if ((HighsInt)hessian.q_value_.size() > hessian_num_nz) hessian.q_value_.resize(hessian_num_nz);

  if (return_status == HighsStatus::kError)
    return_status = HighsStatus::kError;
  else
    return_status = HighsStatus::kOk;
  if (return_status != HighsStatus::kOk)
    highsLogDev(options.log_options, HighsLogType::kInfo,
		"assessHessian returns HighsStatus = %s\n",
		HighsStatusToString(return_status).c_str());
  return return_status;
}

HighsStatus normaliseHessian(const HighsOptions& options, HighsHessian& hessian) {
  // Normalise the Hessian to be (Q + Q^T)/2, where Q is the matrix
  // supplied. This guarantees that what's used internally is
  // symmetric.
  //
  // So someone preferring to supply only the upper triangle would
  // have to double its values..
  HighsStatus return_status = HighsStatus::kOk;
  const HighsInt dim = hessian.dim_;
  const HighsInt hessian_num_nz = hessian.q_start_[dim];
  if (hessian_num_nz > 0 && dim <= 0) return HighsStatus::kError;
  if (hessian_num_nz <= 0) return HighsStatus::kOk;

  bool error_found = false;
  bool warning_found = false;

  HighsHessian transpose;
  transpose.q_start_.resize(dim+1);
  transpose.q_index_.resize(hessian_num_nz);
  transpose.q_value_.resize(hessian_num_nz);
  // Form transpose of Hessian
  vector<HighsInt> qr_length;
  qr_length.assign(dim, 0);
  for (HighsInt iEl = 0; iEl < hessian_num_nz; iEl++)
    qr_length[hessian.q_index_[iEl]]++;

  transpose.q_start_[0] = 0;
  for (HighsInt iRow = 0; iRow < dim; iRow++)
    transpose.q_start_[iRow+1] = transpose.q_start_[iRow] + qr_length[iRow];
  for (HighsInt iCol = 0; iCol < dim; iCol++) {
    for (HighsInt iEl = hessian.q_start_[iCol]; iEl < hessian.q_start_[iCol+1]; iEl++) {
      HighsInt iRow = hessian.q_index_[iEl];
      HighsInt iRowEl = transpose.q_start_[iRow];
      transpose.q_index_[iRowEl] = iCol;
      transpose.q_value_[iRowEl] = hessian.q_value_[iEl];
      transpose.q_start_[iRow]++;
    }
  }
  
  transpose.q_start_[0] = 0;
  for (HighsInt iRow = 0; iRow < dim; iRow++)
    transpose.q_start_[iRow+1] = transpose.q_start_[iRow] + qr_length[iRow];
							  
  HighsHessian normalised;
  HighsInt normalised_num_nz = 0;
  HighsInt normalised_size = hessian_num_nz;
  normalised.q_start_.resize(dim+1);
  normalised.q_index_.resize(normalised_size);
  normalised.q_value_.resize(normalised_size);
  HighsInt column_num_nz = 0;
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
    for (HighsInt iEl = hessian.q_start_[iCol]; iEl < hessian.q_start_[iCol+1]; iEl++) {
      HighsInt iRow = hessian.q_index_[iEl];
      column_value[iRow] = hessian.q_value_[iEl];
      column_index[column_num_nz] = iRow;
      column_num_nz++;
    }
    for (HighsInt iEl = transpose.q_start_[iCol]; iEl < transpose.q_start_[iCol+1]; iEl++) {
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
      normalised_size = std::max(normalised_num_nz + column_num_nz, 2*normalised_size);
      normalised.q_index_.resize(normalised_size);
      normalised.q_value_.resize(normalised_size);
    }
    for (HighsInt iEl = 0; iEl < column_num_nz; iEl++) {
      // Check the value
      HighsInt iRow = column_index[iEl];
      double abs_value = std::fabs(0.5 * column_value[iRow]);
      bool ok_value = abs_value > small_matrix_value;
      if (!ok_value) {
        if (max_small_value < abs_value) max_small_value = abs_value;
        if (min_small_value > abs_value) min_small_value = abs_value;
        num_small_values++;
      }
      if (ok_value) {
        normalised.q_index_[normalised_num_nz] = iRow;
        normalised.q_value_[normalised_num_nz] = column_value[iRow];
        normalised_num_nz++;
      }
      column_value[iRow] = 0;
    }
    normalised.q_start_[iCol+1] = normalised_num_nz;
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
  if (error_found)
    return_status = HighsStatus::kError;
  else if (warning_found)
    return_status = HighsStatus::kWarning;
  else
    return_status = HighsStatus::kOk;

  return return_status;
}
