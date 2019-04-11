/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsUtils.cpp
 * @brief Class-independent utilities for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsLpUtils.h"

#include <cassert>

#include "HConfig.h"
#include "io/HighsIO.h"
#include "lp_data/HighsModelUtils.h"
#include "util/HighsUtils.h"
#include "util/HighsSort.h"
#include "lp_data/HighsStatus.h"

HighsStatus checkLp(const HighsLp& lp) {
  // Check dimensions.
  if (lp.numCol_ < 0 || lp.numRow_ < 0)
    return HighsStatus::LpError;

  // Check vectors.
  if ((int)lp.colCost_.size() != lp.numCol_)
    return HighsStatus::LpError;

  if ((int)lp.colLower_.size() != lp.numCol_ ||
      (int)lp.colUpper_.size() != lp.numCol_)
    return HighsStatus::LpError;
  if ((int)lp.rowLower_.size() != lp.numRow_ ||
      (int)lp.rowUpper_.size() != lp.numRow_)
    return HighsStatus::LpError;

  for (int i = 0; i < lp.numRow_; i++)
    if (lp.rowLower_[i] < -HIGHS_CONST_INF || lp.rowUpper_[i] > HIGHS_CONST_INF)
      return HighsStatus::LpError;

  for (int j = 0; j < lp.numCol_; j++) {
    if (lp.colCost_[j] < -HIGHS_CONST_INF || lp.colCost_[j] > HIGHS_CONST_INF)
      return HighsStatus::LpError;

    if (lp.colLower_[j] < -HIGHS_CONST_INF || lp.colUpper_[j] > HIGHS_CONST_INF)
      return HighsStatus::LpError;
    if (lp.colLower_[j] > lp.colUpper_[j] + kBoundTolerance)
      return HighsStatus::LpError;
  }

  // Check matrix.
  if (lp.numCol_ == 0) return HighsStatus::OK;
  int lp_num_nz = lp.Astart_[lp.numCol_];
  if ((int)lp.Avalue_.size() < lp_num_nz) return HighsStatus::LpError;
  if (lp_num_nz < 0) return HighsStatus::LpError;
  if ((int)lp.Aindex_.size() < lp_num_nz) return HighsStatus::LpError;

  if (lp.numCol_ > 0) {
    if ((int)lp.Astart_.size() < lp.numCol_ + 1) return HighsStatus::LpError;
    // Was lp.Astart_[i] >= lp.nnz_ (below), but this is wrong when the
    // last column is empty. Need to check as follows, and also check
    // the entry lp.Astart_[lp.numCol_] > lp.nnz_
    for (int i = 0; i < lp.numCol_; i++) {
      if (lp.Astart_[i] > lp.Astart_[i + 1] || lp.Astart_[i] > lp_num_nz ||
	  lp.Astart_[i] < 0)
	return HighsStatus::LpError;
    }
    if (lp.Astart_[lp.numCol_] > lp_num_nz ||
	lp.Astart_[lp.numCol_] < 0)
      return HighsStatus::LpError;

    for (int k = 0; k < lp.nnz_; k++) {
      if (lp.Aindex_[k] < 0 || lp.Aindex_[k] >= lp.numRow_)
	return HighsStatus::LpError;
      if (lp.Avalue_[k] < -HIGHS_CONST_INF || lp.Avalue_[k] > HIGHS_CONST_INF)
	return HighsStatus::LpError;
    }
  }

  return HighsStatus::OK;
}

HighsStatus assessLp(HighsLp& lp, const HighsOptions& options, const bool normalise) {
  HighsStatus return_status = HighsStatus::NotSet;
  HighsStatus call_status;
  // Assess the LP dimensions and vector sizes, returning on error
  call_status = assessLpDimensions(lp);
  return_status = worseStatus(call_status, return_status);
  if (return_status == HighsStatus::Error) return HighsStatus::LpError;

  // If the LP has no columns there is nothing left to test
  // NB assessLpDimensions returns HighsStatus::Error if lp.numCol_ < 0
  if (lp.numCol_ == 0) return HighsStatus::OK;

  // From here, any LP has lp.numCol_ > 0 and lp.Astart_[lp.numCol_] exists (as the number of nonzeros)
  assert(lp.numCol_ > 0);

  // Assess the LP column costs
  call_status = assess_costs(0, lp.numCol_, true, 0, lp.numCol_, false, 0, NULL, false, NULL,
			     &lp.colCost_[0], options.infinite_cost);
  return_status = worseStatus(call_status, return_status);
  // Assess the LP column bounds
  call_status = assess_bounds("Col", 0, lp.numCol_, true, 0, lp.numCol_, false, 0, NULL, false, NULL,
			      &lp.colLower_[0], &lp.colUpper_[0], options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);
  // Assess the LP row bounds
  call_status = assess_bounds("Row", 0, lp.numRow_, true, 0, lp.numRow_, false, 0, NULL, false, NULL,
			      &lp.rowLower_[0], &lp.rowUpper_[0], options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);
  // Assess the LP matrix
  int lp_num_nz = lp.Astart_[lp.numCol_];
  call_status = assessMatrix(lp.numRow_, 0, lp.numCol_, lp.numCol_, lp_num_nz,
			     &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0],
			     options.small_matrix_value, options.large_matrix_value, normalise);
  lp.Astart_[lp.numCol_] = lp_num_nz;
  return_status = worseStatus(call_status, return_status);
  if (return_status == HighsStatus::Error) return_status = HighsStatus::LpError;
  else return_status = HighsStatus::OK;
#ifdef HiGHSDEV
  HighsLogMessage(HighsMessageType::INFO, "assess_lp returns HighsStatus = %s", HighsStatusToString(return_status).c_str());
#endif
  return return_status;
}

HighsStatus assessLpDimensions(const HighsLp& lp) {
  HighsStatus return_status = HighsStatus::NotSet;

  // Use error_found to track whether an error has been found in multiple tests
  bool error_found = false;

  // Don't expect the matrix_start_size to be legal if there are no columns
  bool check_matrix_start_size = lp.numCol_ > 0;

  // Assess column-related dimensions
  bool legal_num_col = lp.numCol_ >= 0;
  if (!legal_num_col) {
    HighsLogMessage(HighsMessageType::ERROR, "LP has illegal number of cols = %d\n", lp.numCol_);
    error_found = true;
  } else {
    // Check the size of the column vectors
    int col_cost_size = lp.colCost_.size();
    int col_lower_size = lp.colLower_.size();
    int col_upper_size = lp.colUpper_.size();
    int matrix_start_size = lp.Astart_.size();
    bool legal_col_cost_size = col_cost_size >= lp.numCol_;
    bool legal_col_lower_size = col_lower_size >= lp.numCol_;
    bool legal_col_upper_size = col_lower_size >= lp.numCol_;
    
    if (!legal_col_cost_size) {
      HighsLogMessage(HighsMessageType::ERROR, "LP has illegal colCost size = %d < %d\n", col_cost_size, lp.numCol_);
      error_found = true;
    }
    if (!legal_col_lower_size) {
      HighsLogMessage(HighsMessageType::ERROR, "LP has illegal colLower size = %d < %d\n", col_lower_size, lp.numCol_);
      error_found = true;
    }
    if (!legal_col_upper_size) {
      HighsLogMessage(HighsMessageType::ERROR, "LP has illegal colUpper size = %d < %d\n", col_upper_size, lp.numCol_);
      error_found = true;
    }
    if (check_matrix_start_size) {
      bool legal_matrix_start_size = matrix_start_size >= lp.numCol_+1;
      if (!legal_matrix_start_size) {
	HighsLogMessage(HighsMessageType::ERROR, "LP has illegal Astart size = %d < %d\n", matrix_start_size, lp.numCol_+1);
	error_found = true;
      }
    }
  }

  // Assess row-related dimensions
  bool legal_num_row = lp.numRow_ >= 0;
  if (!legal_num_row) {
    HighsLogMessage(HighsMessageType::ERROR, "LP has illegal number of rows = %d\n", lp.numRow_);
    error_found = true;
  } else {
    int row_lower_size = lp.rowLower_.size();
    int row_upper_size = lp.rowUpper_.size();
    bool legal_row_lower_size = row_lower_size >= lp.numRow_;
    bool legal_row_upper_size = row_lower_size >= lp.numRow_;
    if (!legal_row_lower_size) {
      HighsLogMessage(HighsMessageType::ERROR, "LP has illegal rowLower size = %d < %d\n", row_lower_size, lp.numRow_);
      error_found = true;
    }
    if (!legal_row_upper_size) {
      HighsLogMessage(HighsMessageType::ERROR, "LP has illegal rowUpper size = %d < %d\n", row_upper_size, lp.numRow_);
      error_found = true;
    }
  }

  // Assess matrix-related dimensions
  if (check_matrix_start_size) {
    int lp_num_nz = lp.Astart_[lp.numCol_];
    bool legal_num_nz = lp_num_nz >= 0;
    if (!legal_num_nz) {
      HighsLogMessage(HighsMessageType::ERROR, "LP has illegal number of nonzeros = %d\n", lp_num_nz);
      error_found = true;
    } else {
      int matrix_index_size = lp.Aindex_.size();
      int matrix_value_size = lp.Avalue_.size();
      bool legal_matrix_index_size = matrix_index_size >= lp_num_nz;
      bool legal_matrix_value_size = matrix_value_size >= lp_num_nz;
      if (!legal_matrix_index_size) {
	HighsLogMessage(HighsMessageType::ERROR, "LP has illegal Aindex size = %d < %d\n", matrix_index_size, lp_num_nz);
	error_found = true;
      }
      if (!legal_matrix_value_size) {
	HighsLogMessage(HighsMessageType::ERROR, "LP has illegal Avalue size = %d < %d\n", matrix_value_size, lp_num_nz);
	error_found = true;
      }
    }
  }
  if (error_found) return_status = HighsStatus::Error;
  else return_status = HighsStatus::OK;

  return return_status;  
}

HighsStatus assess_costs(const int col_ix_os,
			 const int col_dim,
			 const bool interval, const int from_col, const int to_col,
			 const bool set, const int num_set_entries, const int* col_set,
			 const bool mask, const int* col_mask,
			 const double* usr_col_cost,
			 const double infinite_cost) {
  // Check parameters for technique and, if OK set the loop limits - in iterator style
  int from_k;
  int to_k;
  HighsStatus return_status = assess_interval_set_mask(col_dim,
						       interval, from_col, to_col,
						       set, num_set_entries, col_set,
						       mask, col_mask,
						       from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  if (from_k >= to_k) return HighsStatus::OK;

  return_status = HighsStatus::NotSet;
  bool error_found = false;
  int usr_col;
  for (int k = from_k; k < to_k; k++) {
    if (interval || mask) {
      usr_col = k;
    } else {
      usr_col = col_set[k];
    }
    int col = col_ix_os + usr_col;
    if (mask && !col_mask[usr_col]) continue;
    double abs_cost = fabs(usr_col_cost[k]);
    bool legal_cost = abs_cost < infinite_cost;
    if (!legal_cost) {
      HighsLogMessage(HighsMessageType::ERROR, "Col  %12d has |cost| of %12g >= %12g",  col, abs_cost, infinite_cost);
      error_found = true;
    }
  }
  if (error_found) return_status = HighsStatus::Error;
  else return_status = HighsStatus::OK;
    
  return return_status;
}

HighsStatus assess_bounds(const char* type, const int ix_os,
			  const int ix_dim,
			  const bool interval, const int from_ix, const int to_ix,
			  const bool set, const int num_set_entries, const int* ix_set,
			  const bool mask, const int* ix_mask,
			  double* usr_lower, double* usr_upper,
			  const double infinite_bound, bool normalise) {
  // Check parameters for technique and, if OK set the loop limits - in iterator style
  int from_k;
  int to_k;
  HighsStatus return_status = assess_interval_set_mask(ix_dim,
						       interval, from_ix, to_ix,
						       set, num_set_entries, ix_set,
						       mask, ix_mask,
						       from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  if (from_k >= to_k) return HighsStatus::OK;

  return_status = HighsStatus::NotSet;
  bool error_found = false;
  bool warning_found = false;
  bool info_found = false;
  int num_infinite_lower_bound = 0;
  int num_infinite_upper_bound = 0;
  int usr_ix;
  for (int k = from_k; k < to_k; k++) {
    if (interval || mask) {
      usr_ix = k;
    } else {
      usr_ix = ix_set[k];
    }
    int ix = ix_os + usr_ix;
    if (mask && !ix_mask[usr_ix]) continue;

    if (!highs_isInfinity(-usr_lower[k])) {
      // Check whether a finite lower bound will be treated as -Infinity      
      bool infinite_lower_bound = usr_lower[k] <= -infinite_bound;
      if (infinite_lower_bound) {
	if (normalise) usr_lower[k] = -HIGHS_CONST_INF;
	num_infinite_lower_bound++;
      }
    }
    if (!highs_isInfinity(usr_upper[k])) {
      // Check whether a finite upper bound will be treated as Infinity      
      bool infinite_upper_bound = usr_upper[k] >= infinite_bound;
      if (infinite_upper_bound) {
	if (normalise) usr_upper[k] = HIGHS_CONST_INF;
	num_infinite_upper_bound++;
      }
    }
    // Check that the lower bound does not exceed the upper bound
    bool legalLowerUpperBound = usr_lower[k] <= usr_upper[k];
    if (!legalLowerUpperBound) {
      // Leave inconsistent bounds to be used to deduce infeasibility
      HighsLogMessage(HighsMessageType::WARNING, "%3s  %12d has inconsistent bounds [%12g, %12g]", type, ix, usr_lower[k], usr_upper[k]);
      warning_found = true;
    }
    // Check that the lower bound is not as much as +Infinity
    bool legalLowerBound = usr_lower[k] < infinite_bound;
    if (!legalLowerBound) {
      HighsLogMessage(HighsMessageType::ERROR, "%3s  %12d has lower bound of %12g >= %12g", type, ix, usr_lower[k], infinite_bound);
      error_found = true;
    }
    // Check that the upper bound is not as little as -Infinity
    bool legalUpperBound = usr_upper[k] > -infinite_bound;
    if (!legalUpperBound) {
      HighsLogMessage(HighsMessageType::ERROR, "%3s  %12d has upper bound of %12g <= %12g", type, ix, usr_upper[k], -infinite_bound);
      error_found = true;
    }
  }
  if (normalise) {
    if (num_infinite_lower_bound) {
      HighsLogMessage(HighsMessageType::INFO, "%3ss:%12d lower bounds exceeding %12g are treated as -Infinity", type, num_infinite_lower_bound, -infinite_bound);
      info_found = true;
    }
    if (num_infinite_upper_bound) {
      HighsLogMessage(HighsMessageType::INFO, "%3ss:%12d upper bounds exceeding %12g are treated as +Infinity", type, num_infinite_upper_bound, infinite_bound);
      info_found = true;
    }
  }

  if (error_found) return_status = HighsStatus::Error;
  else if (warning_found) return_status = HighsStatus::Warning;
  else if (info_found) return_status = HighsStatus::Info;
  else return_status = HighsStatus::OK;
    
  return return_status;
}

HighsStatus assessMatrix(const int vec_dim, const int from_ix, const int to_ix, const int num_vec,
			 int num_nz, int* Xstart, int* Xindex, double* Xvalue,
			 const double small_matrix_value, const double large_matrix_value, const bool normalise) {
  // Uses to_ix in iterator style
  if (from_ix < 0) return HighsStatus::Error;
  if (from_ix >= to_ix) return HighsStatus::OK;
  if (num_nz > 0 && vec_dim <= 0) return HighsStatus::Error;
  if (num_nz <= 0) return HighsStatus::OK;

  HighsStatus return_status = HighsStatus::NotSet;
  bool error_found = false;
  bool warning_found = false;

  // Warn the user if the first start is not zero
  int fromEl = Xstart[0];
  if (fromEl != 0) {
      HighsLogMessage(HighsMessageType::WARNING, "Matrix starts do not begin with 0");
      warning_found = true;
  }
  // Assess the starts
  // Set up previous_start for a fictitious previous empty packed vector
  int previous_start = std::max(0, Xstart[from_ix]);
  for (int ix = from_ix; ix < to_ix; ix++) {
    int this_start = Xstart[ix];
    bool this_start_too_small = this_start < previous_start;
    if (this_start_too_small) {
      HighsLogMessage(HighsMessageType::ERROR, "Matrix packed vector %d has illegal start of %d < %d = previous start", ix, this_start, previous_start);
      return HighsStatus::Error;
    }
    bool this_start_too_big = this_start > num_nz;
    if (this_start_too_big) {
      HighsLogMessage(HighsMessageType::ERROR, "Matrix packed vector %d has illegal start of %d > %d = number of nonzeros", ix, this_start, num_nz);
      return HighsStatus::Error;
    }
  }

  // Assess the indices and values
  // Count the number of acceptable indices/values
  int num_new_nz = Xstart[from_ix];
  int num_small_values = 0;
  // Set up a zeroed vector to detect duplicate indices
  vector<int> check_vector;
  if (vec_dim > 0) check_vector.assign(vec_dim, 0);
  for (int ix = from_ix; ix < to_ix; ix++) {
    int from_el = Xstart[ix];
    int to_el;
    if (ix < num_vec-1) {
      to_el = Xstart[ix+1]-1;
    } else {
      to_el = num_nz-1;
    }
    if (normalise) {
      // Account for any index-value pairs removed so far
      Xstart[ix] = num_new_nz;
    }
    for (int el = from_el; el <= to_el; el++) {
      int component = Xindex[el];
      // Check that the index is non-negative
      bool legal_component = component >= 0;
      if (!legal_component) {
	HighsLogMessage(HighsMessageType::ERROR, "Matrix packed vector %d, entry %d, is illegal index %d", ix, el, component);
	return HighsStatus::Error;
      }
      // Check that the index does not exceed the vector dimension
      legal_component = component < vec_dim;
      if (!legal_component) {
	HighsLogMessage(HighsMessageType::ERROR, "Matrix packed vector %d, entry %d, is illegal index %12d >= %d = vector dimension", ix, el, component, vec_dim);
	return HighsStatus::Error;
      }
      // Check that the index has not already ocurred
      legal_component = check_vector[component] == 0;
      if (!legal_component) {
	HighsLogMessage(HighsMessageType::ERROR, "Matrix packed vector %d, entry %d, is duplicate index %d", ix, el, component);	  
	return HighsStatus::Error;
      }
      // Check that the value is not too large
      double abs_value = fabs(Xvalue[el]);
      bool large_value = abs_value >= large_matrix_value;
      if (large_value) {
	HighsLogMessage(HighsMessageType::ERROR, "Matrix packed vector %d, entry %d, is large value |%g| >= %g", ix, el, abs_value, large_matrix_value);	  
	return HighsStatus::Error;
      }
      bool ok_value = abs_value > small_matrix_value;
      if (!ok_value) {
#ifdef HiGHSDEV
	HighsLogMessage(HighsMessageType::WARNING, "Matrix packed vector %d, entry %d, is small value |%g| <= %g", ix, el, abs_value, small_matrix_value);
#endif
	num_small_values++;
      }
      if (!ok_value && normalise) {
	Xindex[num_new_nz] = Xindex[el];
	Xvalue[num_new_nz] = Xvalue[el];
      } else {
	num_new_nz++;
      }
    }
    // Zero check_vector
    for (int el = Xstart[ix]; el <= num_new_nz-1; el++) check_vector[Xindex[el]] = 0;
#ifdef HiGHSDEV
    // Check zeroing of check vector
    for (int component = 0; component < vec_dim; component++) {
      if (check_vector[component]) error_found;
    }
    if (error_found) HighsLogMessage(HighsMessageType::ERROR, "assessMatrix: check_vector not zeroed");
#endif
  }
  if (num_small_values) {
    HighsLogMessage(HighsMessageType::WARNING, "Matrix packed vector contains %d values less than %g in magnitude: ignored", num_small_values, small_matrix_value);	  
    warning_found = true;
    if (normalise) {
      // Accommodate the loss of these values in any subsequent packed vectors
      for (int ix = to_ix; ix < num_vec; ix++) {
	int from_el = Xstart[ix];
	Xstart[ix] = num_new_nz;
	int to_el;
	if (ix < num_vec) {
	  to_el = Xstart[ix+1]-1;
	} else {
	  to_el = num_nz-1;
	}
	for (int el = Xstart[ix]; el <= to_el; el++) {
	  Xindex[num_new_nz] = Xindex[el];
	  Xvalue[num_new_nz] = Xvalue[el];
	  num_new_nz++;
	}
      }    
      num_nz = num_new_nz;
    }
  }
  if (error_found) return_status = HighsStatus::Error;
  else if (warning_found) return_status = HighsStatus::Warning;
  else return_status = HighsStatus::OK;

  return return_status;

}

HighsStatus add_lp_cols(HighsLp& lp,
			const int num_new_col, const double *XcolCost, const double *XcolLower,  const double *XcolUpper,
			const int num_new_nz, const int *XAstart, const int *XAindex, const double *XAvalue,
			const HighsOptions& options) {
  const bool valid_matrix = true;
  if (num_new_col < 0) return HighsStatus::Error;
  if (num_new_col == 0) return HighsStatus::OK;
  HighsStatus return_status = append_lp_cols(lp, num_new_col, XcolCost, XcolLower, XcolUpper,
					   num_new_nz, XAstart, XAindex, XAvalue,
					   options, valid_matrix);
  // Which of the following two??
  if (return_status == HighsStatus::Error) return HighsStatus::Error;
  //  if (return_status != HighsStatus::OK && return_status != HighsStatus::Warning) return return_status
  lp.numCol_ += num_new_col;
  return return_status;
}

HighsStatus append_lp_cols(HighsLp& lp,
			   const int num_new_col, const double *XcolCost, const double *XcolLower,  const double *XcolUpper,
			   const int num_new_nz, const int *XAstart, const int *XAindex, const double *XAvalue,
			   const HighsOptions& options, const bool valid_matrix) {
  if (num_new_col < 0) return HighsStatus::Error;
  if (num_new_col == 0) return HighsStatus::OK;
  HighsStatus return_status = HighsStatus::NotSet;
  int newNumCol = lp.numCol_ + num_new_col;
  // Assess the bounds and matrix indices, returning on error
  bool normalise = false;
  HighsStatus call_status;
  // Assess the column costs
  call_status = assess_costs(lp.numCol_, num_new_col, true, 0, num_new_col, false, 0, NULL, false, NULL,
			     (double*)XcolCost, options.infinite_cost);
  return_status = worseStatus(call_status, return_status);
  // Assess the column bounds
  call_status = assess_bounds("Col", lp.numCol_, num_new_col, true, 0, num_new_col, false, 0, NULL, false, NULL,
			     (double*)XcolLower, (double*)XcolUpper, options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);
  if (valid_matrix) {
    // Assess the matrix columns
    call_status = assessMatrix(lp.numRow_, 0, num_new_col, num_new_col, 
			       num_new_nz, (int*)XAstart, (int*)XAindex, (double*)XAvalue,
			       options.small_matrix_value, options.large_matrix_value, normalise);
    return_status = worseStatus(call_status, return_status);
  }
  if (return_status == HighsStatus::Error) return return_status;

  // Append the columns to the LP vectors and matrix
  call_status = append_cols_to_lp_vectors(lp, num_new_col, XcolCost, XcolLower, XcolUpper);
  return_status = worseStatus(call_status, return_status);
  if (return_status == HighsStatus::Error) return return_status;

  if (valid_matrix) {
    call_status = append_cols_to_lp_matrix(lp, num_new_col, num_new_nz, XAstart, XAindex, XAvalue);
    return_status = worseStatus(call_status, return_status);
    if (return_status == HighsStatus::Error) return return_status;
  }

  // Normalise the new LP column bounds
  normalise = true;
  call_status = assess_bounds("Col", lp.numCol_, num_new_col, true, 0, num_new_col, false, 0, NULL, false, NULL,
			     &lp.colLower_[0], &lp.colUpper_[0], options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);
  if (return_status == HighsStatus::Error) return return_status;
  if (valid_matrix) {
    // Normalise the new LP matrix columns
    int lp_num_nz = lp.Astart_[newNumCol];
    call_status = assessMatrix(lp.numRow_, 0, num_new_col, num_new_col,
			       lp_num_nz, &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0],
			       options.small_matrix_value, options.large_matrix_value, normalise);
    lp.Astart_[newNumCol] = lp_num_nz;
    return_status = worseStatus(call_status, return_status);
    if (return_status == HighsStatus::Error) return return_status;
  }
  return return_status;
}

HighsStatus append_cols_to_lp_vectors(HighsLp &lp, const int num_new_col,
				      const double*XcolCost, const double *XcolLower, const double *XcolUpper) {
  if (num_new_col < 0) return HighsStatus::Error;
  if (num_new_col == 0) return HighsStatus::OK;
  int newNumCol = lp.numCol_ + num_new_col;
  lp.colCost_.resize(newNumCol);
  lp.colLower_.resize(newNumCol);
  lp.colUpper_.resize(newNumCol);
  
  for (int col = 0; col < num_new_col; col++) {
    lp.colCost_[lp.numCol_ + col] = XcolCost[col];
    lp.colLower_[lp.numCol_ + col] = XcolLower[col];
    lp.colUpper_[lp.numCol_ + col] = XcolUpper[col];
  }
  return HighsStatus::OK;
}

HighsStatus add_lp_rows(HighsLp& lp,
			const int num_new_row, const double *XrowLower,  const double *XrowUpper,
			const int num_new_nz, const int *XARstart, const int *XARindex, const double *XARvalue,
			const HighsOptions& options) {
  const bool valid_matrix = true;
  if (num_new_row < 0) return HighsStatus::Error;
  if (num_new_row == 0) return HighsStatus::OK;
  HighsStatus return_status = append_lp_rows(lp, num_new_row, XrowLower, XrowUpper,
					     num_new_nz, XARstart, XARindex, XARvalue,
					     options, valid_matrix);
  // Which of the following two??
  if (return_status == HighsStatus::Error) return HighsStatus::Error;
  //  if (return_status != HighsStatus::OK && return_status != HighsStatus::Warning) return return_status
  lp.numRow_ += num_new_row;
  return return_status;
}

HighsStatus append_lp_rows(HighsLp& lp,
			   const int num_new_row, const double *XrowLower,  const double *XrowUpper,
			   const int num_new_nz, const int *XARstart, const int *XARindex, const double *XARvalue,
			   const HighsOptions& options, bool valid_matrix) {
  if (num_new_row < 0) return HighsStatus::Error;
  if (num_new_row == 0) return HighsStatus::OK;
  HighsStatus return_status = HighsStatus::NotSet;
  int new_num_row = lp.numRow_ + num_new_row;
  // Assess the bounds and matrix indices, returning on error
  bool normalise = false;
  HighsStatus call_status;
  // Assess the row bounds
  call_status = assess_bounds("Row", lp.numRow_, num_new_row, true, 0, num_new_row, false, 0, NULL, false, NULL,
			     (double*)XrowLower, (double*)XrowUpper, options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);
  if (valid_matrix) {
    // Assess the matrix columns
    call_status = assessMatrix(lp.numCol_, 0, num_new_row, num_new_row, 
			       num_new_nz, (int*)XARstart, (int*)XARindex, (double*)XARvalue,
			       options.small_matrix_value, options.large_matrix_value, normalise);
    return_status = worseStatus(call_status, return_status);
  }
  if (return_status == HighsStatus::Error) return return_status;

  // Append the rows to the LP vectors
  call_status = append_rows_to_lp_vectors(lp, num_new_row, XrowLower, XrowUpper);
  return_status = worseStatus(call_status, return_status);
  if (return_status == HighsStatus::Error) return return_status;

  // Normalise the new LP row bounds
  normalise = true;
  call_status = assess_bounds("Row", lp.numRow_, num_new_row, true, 0, num_new_row, false, 0, NULL, false, NULL,
			      &lp.rowLower_[0], &lp.rowUpper_[0], options.infinite_bound, normalise);
  return_status = worseStatus(call_status, return_status);

  if (valid_matrix) {
    // Copy the supplied row-wise matrix so it can be normalised before being appended
    int lc_num_new_nz = num_new_nz;
    int *lc_row_matrix_start = (int *)malloc(sizeof(int) * num_new_row);
    int *lc_row_matrix_index = (int *)malloc(sizeof(int) * lc_num_new_nz);
    double *lc_row_matrix_value = (double *)malloc(sizeof(double) * lc_num_new_nz);
    for (int row=0; row < num_new_row; row++) {
      lc_row_matrix_start[row] = XARstart[row];
    }
    for (int el=0; el < lc_num_new_nz; el++) {
      lc_row_matrix_index[el] = XARindex[el];
      lc_row_matrix_value[el] = XARvalue[el];
    }
    call_status = assessMatrix(lp.numCol_, 0, num_new_row, num_new_row, 
			       lc_num_new_nz, lc_row_matrix_start, lc_row_matrix_index, lc_row_matrix_value,
			       options.small_matrix_value, options.large_matrix_value, normalise);
    return_status = worseStatus(call_status, return_status);
    if (return_status == HighsStatus::Error) return return_status;
    // Append the matrix to the LP vectors
    call_status = append_rows_to_lp_matrix(lp, num_new_row, num_new_nz,
					   lc_row_matrix_start, lc_row_matrix_index, lc_row_matrix_value);
    return_status = worseStatus(call_status, return_status);
    if (return_status == HighsStatus::Error) return return_status;
  }
  return return_status;
}

HighsStatus append_rows_to_lp_vectors(HighsLp &lp, const int num_new_row,
				      const double *XrowLower, const double *XrowUpper) {
  if (num_new_row < 0) return HighsStatus::Error;
  if (num_new_row == 0) return HighsStatus::OK;
  int new_num_row = lp.numRow_ + num_new_row;
  lp.rowLower_.resize(new_num_row);
  lp.rowUpper_.resize(new_num_row);

  for (int row = 0; row < num_new_row; row++) {
    lp.rowLower_[lp.numRow_ + row] = XrowLower[row];
    lp.rowUpper_[lp.numRow_ + row] = XrowUpper[row];
  }
  return HighsStatus::OK;
}

HighsStatus append_cols_to_lp_matrix(HighsLp &lp, const int num_new_col,
				     const int num_new_nz, const int *XAstart, const int *XAindex, const double *XAvalue) {
  if (num_new_col < 0) return HighsStatus::Error;
  if (num_new_col == 0) return HighsStatus::OK;
  // Check that nonzeros aren't being appended to a matrix with no rows
  if (num_new_nz > 0 && lp.numRow_ <= 0) return HighsStatus::Error;
  // Determine the new number of columns in the matrix and resize the
  // starts accordingly. Even if no matrix entries are added, f they
  // are added later as rows it will be assumesd that the starts are
  // of the right size.
  int newNumCol = lp.numCol_ + num_new_col;
  lp.Astart_.resize(newNumCol + 1);
  // If adding columns to an empty LP then introduce the start for the fictitious column 0
  if (lp.numCol_ == 0) lp.Astart_[0] = 0;
  // If no nonzeros are bing added then there's nothing to do
  if (num_new_nz <= 0) return HighsStatus::OK;
  //
  // Adding a non-trivial matrix so determine the current number of nonzeros
  int current_num_nz = lp.Astart_[lp.numCol_];
  
  // Determine the new number of nonzeros and resize the column-wise matrix arrays accordingly
  int new_num_nz = current_num_nz + num_new_nz;
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);
  
  // Append the new columns
  for (int col = 0; col < num_new_col; col++) lp.Astart_[lp.numCol_ + col] = current_num_nz + XAstart[col];
  lp.Astart_[lp.numCol_ + num_new_col] = new_num_nz;
  
  for (int el = 0; el < num_new_nz; el++) {
    lp.Aindex_[current_num_nz + el] = XAindex[el];
    lp.Avalue_[current_num_nz + el] = XAvalue[el];
  }
  return HighsStatus::OK;
}

HighsStatus append_rows_to_lp_matrix(HighsLp &lp, const int num_new_row,
				     const int num_new_nz, const int *XARstart, const int *XARindex, const double *XARvalue) {
  if (num_new_row < 0) return HighsStatus::Error;
  if (num_new_row == 0) return HighsStatus::OK;
  // Check that nonzeros aren't being appended to a matrix with no columns
  if (num_new_nz > 0 && lp.numCol_ <= 0) return HighsStatus::Error;
  int new_num_row = lp.numRow_ + num_new_row;
  if (num_new_nz == 0) return HighsStatus::OK;
  int current_num_nz = lp.Astart_[lp.numCol_];
  vector<int> Alength;
  Alength.assign(lp.numCol_, 0);
  for (int el = 0; el < num_new_nz; el++)Alength[XARindex[el]]++;
  // Determine the new number of nonzeros and resize the column-wise matrix arrays
  int new_num_nz = current_num_nz + num_new_nz;
  lp.Aindex_.resize(new_num_nz);
  lp.Avalue_.resize(new_num_nz);

  // Append the new rows
  // Shift the existing columns to make space for the new entries
  int new_el = new_num_nz;
  for (int col = lp.numCol_ - 1; col >= 0; col--) {
    int start_col_plus_1 = new_el;
    new_el -= Alength[col];
    for (int el = lp.Astart_[col + 1] - 1; el >= lp.Astart_[col]; el--) {
      new_el--;
      lp.Aindex_[new_el] = lp.Aindex_[el];
      lp.Avalue_[new_el] = lp.Avalue_[el];
    }
    lp.Astart_[col + 1] = start_col_plus_1;
  }
  assert(new_el == 0);

  // Insert the new entries
  for (int row = 0; row < num_new_row; row++) {
    int first_el = XARstart[row];
    int last_el = (row < num_new_row - 1 ? XARstart[row + 1] : num_new_nz);
    for (int el = first_el; el < last_el; el++) {
      int col = XARindex[el];
      new_el = lp.Astart_[col + 1] - Alength[col];
      Alength[col]--;
      lp.Aindex_[new_el] = lp.numRow_ + row;
      lp.Avalue_[new_el] = XARvalue[el];
    }
  }
  return HighsStatus::OK;
}

HighsStatus delete_lp_cols(HighsLp &lp,
			   const bool interval, const int from_col, const int to_col,
			   const bool set, const int num_set_entries, const int* col_set,
			   const bool mask, int* col_mask,
			   const bool valid_matrix) {
  int new_num_col;
  HighsStatus call_status = delete_cols_from_lp_vectors(lp, new_num_col,
							interval, from_col, to_col,
							set, num_set_entries, col_set,
							mask, col_mask);
  if (call_status != HighsStatus::OK) return call_status;
  if (valid_matrix) {
    HighsStatus call_status = delete_cols_from_lp_matrix(lp,
							 interval, from_col, to_col,
							 set, num_set_entries, col_set,
							 mask, col_mask);
    if (call_status != HighsStatus::OK) return call_status;
  }
  lp.numCol_ = new_num_col;
  return HighsStatus::OK;
}

HighsStatus delete_cols_from_lp_vectors(HighsLp &lp, int &new_num_col,
					const bool interval, const int from_col, const int to_col,
					const bool set, const int num_set_entries, const int* col_set,
					const bool mask, const int* col_mask) {
  int from_k;
  int to_k;
  HighsStatus return_status = assess_interval_set_mask(lp.numCol_,
						     interval, from_col, to_col,
						     set, num_set_entries, col_set,
						     mask, col_mask,
						     from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  // Initialise new_num_col in case none is removed due to from_k >= to_k
  new_num_col = lp.numCol_;
  if (from_k >= to_k) return HighsStatus::OK;

  int delete_from_col;
  int delete_to_col;
  int keep_from_col;
  int keep_to_col = 0;
  int current_set_entry = 0;
  int col_dim = lp.numCol_;
  new_num_col = 0;
  for (int k = from_k; k < to_k; k++) {
    update_out_in_ix(col_dim,
		     interval, from_col, to_col,
		     set, num_set_entries, col_set,
		     mask, col_mask,
		     delete_from_col, delete_to_col,
		     keep_from_col, keep_to_col,
		     current_set_entry);
    if (delete_to_col == col_dim) break;
     assert(delete_to_col < col_dim);
     if (k == from_k) {
       // Account for the initial columns being kept
       new_num_col = delete_from_col;
     }
     for (int col = keep_from_col; col < keep_to_col; col++) {
       lp.colCost_[new_num_col] = lp.colCost_[col];
       lp.colLower_[new_num_col] = lp.colLower_[col];
       lp.colUpper_[new_num_col] = lp.colUpper_[col];
       new_num_col++;
     }
     if (keep_to_col == col_dim) break;
  }
  return HighsStatus::OK;
}

HighsStatus delete_cols_from_lp_matrix(HighsLp &lp,
				       const bool interval, const int from_col, const int to_col,
				       const bool set, const int num_set_entries, const int* col_set,
				       const bool mask, int* col_mask) {
  int from_k;
  int to_k;
  HighsStatus return_status = assess_interval_set_mask(lp.numCol_,
						       interval, from_col, to_col,
						       set, num_set_entries, col_set,
						       mask, col_mask,
						       from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  if (from_k >= to_k) return HighsStatus::OK;

  int delete_from_col;
  int delete_to_col;
  int keep_from_col;
  int keep_to_col = 0;
  int current_set_entry = 0;
  int col_dim = lp.numCol_;
  int new_num_col = 0;
  int new_num_nz = 0;
  for (int k = from_k; k < to_k; k++) {
    update_out_in_ix(col_dim,
		     interval, from_col, to_col,
		     set, num_set_entries, col_set,
		     mask, col_mask,
		     delete_from_col, delete_to_col,
		     keep_from_col, keep_to_col,
		     current_set_entry);
     if (k == from_k) {
       // Account for the initial columns being kept
       if (mask) {
	 for (int col = 0; col < delete_from_col; col++) {
	   col_mask[col] = new_num_col;
	   new_num_col++;
	 }
       } else {
	 new_num_col = delete_from_col;
       }
       new_num_nz = lp.Astart_[delete_from_col];
     }
     // Ensure that the starts of the deleted columns are zeroed to
     // avoid redundant start information for columns whose indices
     // are't used after the deletion takes place. In particular, if
     // all columns are deleted then something must be done to ensure
     // that the matrix isn't magially recreated by increasing the
     // number of columns from zero when there are no rows in the LP.
     for (int col = delete_from_col; col < delete_to_col; col++) lp.Astart_[col] = 0;
     for (int col = keep_from_col; col < keep_to_col; col++) {
       lp.Astart_[new_num_col] = new_num_nz + lp.Astart_[col] - lp.Astart_[keep_from_col];
       new_num_col++;
     }
     for (int el = lp.Astart_[keep_from_col]; el < lp.Astart_[keep_to_col]; el++) {
      lp.Aindex_[new_num_nz] = lp.Aindex_[el];
      lp.Avalue_[new_num_nz] = lp.Avalue_[el];
      new_num_nz++;
    }
    if (keep_to_col == col_dim) break;
  }
  // Ensure that the start of the spurious last column is zeroed so
  // that it doesn't give a positive number of matrix entries if the
  // number of columns in the LP is increased when there are no rows
  // in the LP.
  lp.Astart_[lp.numCol_] = 0;
  lp.Astart_[new_num_col] = new_num_nz;
  return HighsStatus::OK;
}

HighsStatus delete_lp_rows(HighsLp &lp,
			   const bool interval, const int from_row, const int to_row,
			   const bool set, const int num_set_entries, const int* row_set,
			   const bool mask, int* row_mask,
			   const bool valid_matrix) {
  int new_num_row;
  HighsStatus return_status = delete_rows_from_lp_vectors(lp, new_num_row,
							  interval, from_row, to_row,
							  set, num_set_entries, row_set,
							  mask, row_mask);
  if (return_status != HighsStatus::OK) return return_status;
  if (valid_matrix) {
    HighsStatus return_status = delete_rows_from_lp_matrix(lp,
							   interval, from_row, to_row,
							   set, num_set_entries, row_set,
							   mask, row_mask);
    if (return_status != HighsStatus::OK) return return_status;
  }
  lp.numRow_ = new_num_row;
  return HighsStatus::OK;
}

HighsStatus delete_rows_from_lp_vectors(HighsLp &lp, int &new_num_row,
					const bool interval, const int from_row, const int to_row,
					const bool set, const int num_set_entries, const int* row_set,
					const bool mask, const int* row_mask) {
  int from_k;
  int to_k;
  HighsStatus return_status = assess_interval_set_mask(lp.numRow_,
						     interval, from_row, to_row,
						     set, num_set_entries, row_set,
						     mask, row_mask,
						     from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  // Initialise new_num_row in case none is removed due to from_k >= to_k
  new_num_row = lp.numRow_;
  if (from_k >= to_k) return HighsStatus::OK;

  int delete_from_row;
  int delete_to_row;
  int keep_from_row;
  int keep_to_row = 0;
  int current_set_entry = 0;
  int row_dim = lp.numRow_;
  new_num_row = 0;
  for (int k = from_k; k < to_k; k++) {
    update_out_in_ix(row_dim,
		     interval, from_row, to_row,
		     set, num_set_entries, row_set,
		     mask, row_mask,
		     delete_from_row, delete_to_row,
		     keep_from_row, keep_to_row,
		     current_set_entry);
    if (delete_to_row == row_dim) break;
     assert(delete_to_row < row_dim);
     if (k == from_k) {
       // Account for the initial rows being kept
       new_num_row = delete_from_row;
     }
     for (int row = keep_from_row; row < keep_to_row; row++) {
       lp.rowLower_[new_num_row] = lp.rowLower_[row];
       lp.rowUpper_[new_num_row] = lp.rowUpper_[row];
       new_num_row++;
     }
     if (keep_to_row == row_dim) break;
  }
  return HighsStatus::OK;
}

HighsStatus delete_rows_from_lp_matrix(HighsLp &lp, 
				       const bool interval, const int from_row, const int to_row,
				       const bool set, const int num_set_entries, const int* row_set,
				       const bool mask, int* row_mask) {
  int from_k;
  int to_k;
  HighsStatus return_status = assess_interval_set_mask(lp.numRow_,
						     interval, from_row, to_row,
						     set, num_set_entries, row_set,
						     mask, row_mask,
						     from_k, to_k);
  if (return_status != HighsStatus::OK) return return_status;
  if (from_k >= to_k) return HighsStatus::OK;

  int delete_from_row;
  int delete_to_row;
  int keep_from_row;
  int row_dim = lp.numRow_;
  int keep_to_row;
  int current_set_entry;

  // Set up a row mask to indicate the new row index of kept rows and
  // -1 for deleted rows so that the kept entries in the column-wise
  // matrix can be identified and have their correct row index.
  int *new_index = (int *)malloc(sizeof(int) * lp.numRow_);
  int new_num_row = 0;
  if (!mask) {
    keep_to_row = 0;
    current_set_entry = 0;
    for (int k = from_k; k < to_k; k++) {
      update_out_in_ix(row_dim,
		       interval, from_row, to_row,
		       set, num_set_entries, row_set,
		       mask, row_mask,
		       delete_from_row, delete_to_row,
		       keep_from_row, keep_to_row,
		       current_set_entry);
      if (k == from_k) {
	// Account for any initial rows being kept
	for (int row = 0; row < delete_from_row; row++) {
	  new_index[row] = new_num_row;
	  new_num_row++;
	}
      }
      for (int row = delete_from_row; row < delete_to_row; row++) {
	new_index[row] = -1;
      }
      for (int row = keep_from_row; row < keep_to_row; row++) {
	new_index[row] = new_num_row;
	new_num_row++;
      }
      if (keep_to_row == row_dim) break;
    }
  } else {
    for (int row = 0; row < lp.numRow_; row++) {
      if (row_mask[row]) {
	new_index[row] = -1;
	row_mask[row] = new_index[row];
      } else {
	new_index[row] = new_num_row;
	row_mask[row] = new_index[row];
	new_num_row++;
      }
    }
  }
  int new_num_nz = 0;
  for (int col = 0; col < lp.numCol_; col++) {
    int from_el = lp.Astart_[col];
    lp.Astart_[col] = new_num_nz;
    for (int el = from_el; el < lp.Astart_[col+1]; el++) {
      int row = lp.Aindex_[el];
      int new_row = new_index[row];
      if (new_row >= 0) {
	lp.Aindex_[new_num_nz] = new_row;
	lp.Avalue_[new_num_nz] = lp.Avalue_[el];
	new_num_nz++;
      }
    }
  }
  lp.Astart_[lp.numCol_] = new_num_nz;
  return HighsStatus::OK;
}
    
HighsStatus change_lp_matrix_coefficient(HighsLp &lp, const int row, const int col, const double new_value) {
  if (row < 0 || row > lp.numRow_) return HighsStatus::Error;
  if (col < 0 || col > lp.numCol_) return HighsStatus::Error;
  int changeElement = -1;
  for (int el = lp.Astart_[col]; el < lp.Astart_[col + 1]; el++) {
    // printf("Column %d: Element %d is row %d. Is it %d?\n", col, el, lp.Aindex_[el], row);
    if (lp.Aindex_[el] == row) {
      changeElement = el;
      break;
    }
  }
  if (changeElement < 0) {
    //    printf("util_changeCoeff: Cannot find row %d in column %d\n", row, col);
    changeElement = lp.Astart_[col + 1];
    int new_num_nz = lp.Astart_[lp.numCol_] + 1;
    //    printf("model.util_changeCoeff: Increasing Nnonz from %d to %d\n",
    //    lp.Astart_[lp.numCol_], new_num_nz);
    lp.Aindex_.resize(new_num_nz);
    lp.Avalue_.resize(new_num_nz);
    for (int i = col + 1; i <= lp.numCol_; i++) lp.Astart_[i]++;
    for (int el = new_num_nz - 1; el > changeElement; el--) {
      lp.Aindex_[el] = lp.Aindex_[el - 1];
      lp.Avalue_[el] = lp.Avalue_[el - 1];
    }
  }
  lp.Aindex_[changeElement] = row;
  lp.Avalue_[changeElement] = new_value;
}

HighsStatus change_lp_costs(HighsLp &lp,
			    const bool interval, const int from_col, const int to_col,
			    const bool set, const int num_set_entries, const int* col_set,
			    const bool mask, const int* col_mask,
			    const double* usr_col_cost,
			    const double infinite_cost) {
  // Check parameters for technique and, if OK set the loop limits - in iterator style
  int from_k;
  int to_k;
  HighsStatus call_status = assess_interval_set_mask(lp.numCol_,
						     interval, from_col, to_col,
						     set, num_set_entries, col_set,
						     mask, col_mask,
						     from_k, to_k);
  HighsStatus return_status = HighsStatus::NotSet;
  if (call_status != HighsStatus::OK) {
    return_status = call_status;
    return return_status;
  }
  if (from_k >= to_k) return HighsStatus::OK;
  if (usr_col_cost == NULL) return HighsStatus::Error;

  // Assess the user costs and return on error
  call_status = assess_costs(0,
			     lp.numCol_,
			     interval, from_col, to_col,
			     set, num_set_entries, col_set,
			     mask, col_mask,
			     usr_col_cost, infinite_cost);
  if (call_status != HighsStatus::OK) {
    return_status = call_status;
    return return_status;
  }
  // Change the costs to the user-supplied costs, according to the technique
  int usr_col;
  for (int k = from_k; k < to_k; k++) {
    if (interval || mask) {
      usr_col = k;
    } else {
      usr_col = col_set[k];
    }
    int col = usr_col;
    if (mask && !col_mask[col]) continue;
    lp.colCost_[col] = usr_col_cost[k];
  }
  return HighsStatus::OK;
}

HighsStatus change_lp_col_bounds(
				 HighsLp &lp,
				 const bool interval, const int from_col, const int to_col,
				 const bool set, const int num_set_entries, const int* col_set,
				 const bool mask, const int* col_mask,
				 const double* usr_col_lower,
				 const double* usr_col_upper,
				 const double infinite_bound
				 ) {
  return change_bounds("col", &lp.colLower_[0], &lp.colUpper_[0], 
		       lp.numCol_,
		       interval, from_col, to_col,
		       set, num_set_entries, col_set,
		       mask, col_mask,
		       usr_col_lower, usr_col_upper,
		       infinite_bound);
}

HighsStatus change_lp_row_bounds(
				 HighsLp &lp,
				 const bool interval, const int from_row, const int to_row,
				 const bool set, const int num_set_entries, const int* row_set,
				 const bool mask, const int* row_mask,
				 const double* usr_row_lower,
				 const double* usr_row_upper,
				 const double infinite_bound
				 ) {
  return change_bounds("row", &lp.rowLower_[0], &lp.rowUpper_[0], 
		       lp.numRow_,
		       interval, from_row, to_row,
		       set, num_set_entries, row_set,
		       mask, row_mask,
		       usr_row_lower, usr_row_upper,
		       infinite_bound);
}

HighsStatus change_bounds(
			  const char* type,
			  double* lower,
			  double* upper,
			  const int ix_dim,
			  const bool interval, const int from_ix, const int to_ix,
			  const bool set, const int num_set_entries, const int* ix_set,
			  const bool mask, const int* ix_mask,
			  const double* usr_lower,
			  const double* usr_upper,
			  const double infinite_bound
			  ) {
  // Check parameters for technique and, if OK set the loop limits - in iterator style
  int from_k;
  int to_k;
  HighsStatus call_status = assess_interval_set_mask(ix_dim,
						       interval, from_ix, to_ix,
						       set, num_set_entries, ix_set,
						       mask, ix_mask,
						       from_k, to_k);
  HighsStatus return_status = HighsStatus::NotSet;
  if (call_status != HighsStatus::OK) {
    return_status = call_status;
    return return_status;
  }
  if (from_k >= to_k) return HighsStatus::OK;
  if (usr_lower == NULL) return HighsStatus::Error;
  if (usr_upper == NULL) return HighsStatus::Error;

  // Assess the user bounds and return on error
  bool normalise = false;
  call_status = assess_bounds(type, 0,
			      ix_dim,
			      interval, from_ix, to_ix,
			      set, num_set_entries, ix_set,
			      mask, ix_mask,
			      (double*)usr_lower, (double*)usr_upper, infinite_bound, normalise);
  if (call_status != HighsStatus::OK) {
    return_status = call_status;
    return return_status;
  }
  // Change the bounds to the user-supplied bounds, according to the technique
  int usr_ix;
  for (int k = from_k; k < to_k; k++) {
    if (interval || mask) {
      usr_ix = k;
    } else {
      usr_ix = ix_set[k];
    }
    int ix = usr_ix;
    if (mask && !ix_mask[ix]) continue;
    lower[ix] = usr_lower[k];
    upper[ix] = usr_upper[k];
  }
  normalise = true;
  call_status = assess_bounds(type, 0,
			      ix_dim,
			      interval, from_ix, to_ix,
			      set, num_set_entries, ix_set,
			      mask, ix_mask,
			      lower, upper, infinite_bound, normalise);
  if (call_status != HighsStatus::OK) {
    return_status = call_status;
    return return_status;
  }
  return HighsStatus::OK;
}

HighsStatus getLpCosts(const HighsLp& lp, const int from_col, const int to_col, double* XcolCost) {
  // Uses to_col in iterator style
  if (from_col < 0 || to_col > lp.numCol_) return HighsStatus::Error;
  if (from_col > to_col) return HighsStatus::OK;
  for (int col = from_col; col < to_col; col++) XcolCost[col - from_col] = lp.colCost_[col];
  return HighsStatus::OK;
}

HighsStatus getLpColBounds(const HighsLp& lp, const int from_col, const int to_col, double* XcolLower, double* XcolUpper) {
  // Uses to_col in iterator style
  if (from_col < 0 || to_col > lp.numCol_) return HighsStatus::Error;
  if (from_col > to_col) return HighsStatus::OK;
  for (int col = from_col; col < to_col; col++) {
    if (XcolLower != NULL) XcolLower[col - from_col] = lp.colLower_[col];
    if (XcolUpper != NULL) XcolUpper[col - from_col] = lp.colUpper_[col];
  }
  return HighsStatus::OK;
}

HighsStatus getLpRowBounds(const HighsLp& lp, const int from_row, const int to_row, double* XrowLower, double* XrowUpper) {
  // Uses to_row in iterator style
  if (from_row < 0 || to_row > lp.numRow_) return HighsStatus::Error;
  if (from_row > to_row) return HighsStatus::OK;
  for (int row = from_row; row < to_row; row++) {
    if (XrowLower != NULL) XrowLower[row - from_row] = lp.rowLower_[row];
    if (XrowUpper != NULL) XrowUpper[row - from_row] = lp.rowUpper_[row];
  }
  return HighsStatus::OK;
}

// Get a single coefficient from the matrix
HighsStatus getLpMatrixCoefficient(const HighsLp& lp, const int Xrow, const int Xcol, double *val) {
#ifdef HiGHSDEV
  printf("Called model.util_getCoeff(row=%d, col=%d)\n", Xrow, Xcol);
#endif
  if (Xrow < 0 || Xrow > lp.numRow_) return HighsStatus::Error;
  if (Xcol < 0 || Xcol > lp.numCol_) return HighsStatus::Error;

  int get_el = -1;
  for (int el = lp.Astart_[Xcol]; el < lp.Astart_[Xcol + 1]; el++) {
    if (lp.Aindex_[el] == Xrow) {
      get_el = el;
      break;
    }
  }
  if (get_el < 0) {
    *val = 0;
  } else {
    *val = lp.Avalue_[get_el];
  }
  return HighsStatus::OK;
}

// Methods for reporting an LP, including its row and column data and matrix
//
// Report the whole LP
void reportLp(const HighsLp &lp, const int report_level) {
  reportLpBrief(lp);
  if (report_level >= 1) {
    reportLpColVectors(lp);
    reportLpRowVectors(lp);
    if (report_level >= 2) reportLpColMatrix(lp);
  }
}

// Report the LP briefly
void reportLpBrief(const HighsLp &lp) {
  reportLpDimensions(lp);
  reportLpObjSense(lp);
}

// Report the LP dimensions
void reportLpDimensions(const HighsLp &lp) {
  int lp_num_nz;
  if (lp.numCol_ == 0) lp_num_nz = 0;
  else lp_num_nz = lp.Astart_[lp.numCol_];
  HighsPrintMessage(ML_MINIMAL,
                    "LP has %d columns, %d rows and %d nonzeros\n",
                    lp.numCol_, lp.numRow_, lp_num_nz);
}

// Report the LP objective sense
void reportLpObjSense(const HighsLp &lp) {
  if (lp.sense_ == OBJSENSE_MINIMIZE)
    HighsPrintMessage(ML_MINIMAL, "Objective sense is minimize\n");
  else if (lp.sense_ == OBJSENSE_MAXIMIZE)
    HighsPrintMessage(ML_MINIMAL, "Objective sense is maximize\n");
  else
    HighsPrintMessage(ML_MINIMAL,
                      "Objective sense is ill-defined as %d\n", lp.sense_);
}

// Report the vectors of LP column data
void reportLpColVectors(const HighsLp &lp) {
  if (lp.numCol_ <= 0) return;
  HighsPrintMessage(ML_VERBOSE,
                    "  Column        Lower        Upper         Cost\n");
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    HighsPrintMessage(ML_VERBOSE, "%8d %12g %12g %12g\n", iCol,
                      lp.colLower_[iCol], lp.colUpper_[iCol], lp.colCost_[iCol]);
  }
}

// Report the vectors of LP row data
void reportLpRowVectors(const HighsLp &lp) {
  if (lp.numRow_ <= 0) return;
  HighsPrintMessage(ML_VERBOSE,
                    "     Row        Lower        Upper\n");
  for (int iRow = 0; iRow < lp.numRow_; iRow++) {
    HighsPrintMessage(ML_VERBOSE, "%8d %12g %12g\n", iRow,
                      lp.rowLower_[iRow], lp.rowUpper_[iRow]);
  }
}

// Report the LP column-wise matrix
void reportLpColMatrix(const HighsLp &lp) {
  if (lp.numCol_ <= 0) return;
  reportMatrix("Column", lp.numCol_, lp.Astart_[lp.numCol_], &lp.Astart_[0], &lp.Aindex_[0], &lp.Avalue_[0]);
}

void reportMatrix(const char* message, const int num_col, const int num_nz, const int* start, const int* index, const double* value) {
  if (num_col <= 0) return;
  HighsPrintMessage(ML_VERBOSE, "%6s Index              Value\n", message);
  for (int col = 0; col < num_col; col++) {
    HighsPrintMessage(ML_VERBOSE, "    %8d Start   %10d\n", col, start[col]);
    int to_el = (col < num_col-1 ? start[col+1] : num_nz);
    for (int el = start[col]; el < to_el; el++) HighsPrintMessage(ML_VERBOSE, "          %8d %12g\n", index[el], value[el]);
  }
  HighsPrintMessage(ML_VERBOSE, "             Start   %10d\n", num_nz);
}

/*
void reportLpSolution(HighsModelObject &highs_model) {
  HighsLp lp = highs_model.simplex_lp_;
  reportLpBrief(lp);
  //  simplex_interface.report_simplex_solution_status();
  assert(lp.numCol_ > 0);
  assert(lp.numRow_ > 0);
  vector<double> colPrimal(lp.numCol_);
  vector<double> colDual(lp.numCol_);
  vector<int> colStatus(lp.numCol_);
  vector<double> rowPrimal(lp.numRow_);
  vector<double> rowDual(lp.numRow_);
  vector<int> rowStatus(lp.numRow_);
  //  util_getPrimalDualValues(colPrimal, colDual, rowPrimal, rowDual);
  //  if (util_convertWorkingToBaseStat(&colStatus[0], &rowStatus[0])) return;
  //  util_reportColVecSol(lp.numCol_, lp.colCost_, lp.colLower_, lp.colUpper_, colPrimal, colDual, colStatus);
  //  util_reportRowVecSol(lp.numRow_, lp.rowLower_, lp.rowUpper_, rowPrimal, rowDual, rowStatus);
}
*/



#ifdef HiGHSDEV
void util_analyseLp(const HighsLp &lp, const char *message) {
  printf("\n%s model data: Analysis\n", message);
  util_analyseVectorValues("Column costs", lp.numCol_, lp.colCost_, false);
  util_analyseVectorValues("Column lower bounds", lp.numCol_, lp.colLower_, false);
  util_analyseVectorValues("Column upper bounds", lp.numCol_, lp.colUpper_, false);
  util_analyseVectorValues("Row lower bounds", lp.numRow_, lp.rowLower_, false);
  util_analyseVectorValues("Row upper bounds", lp.numRow_, lp.rowUpper_, false);
  util_analyseVectorValues("Matrix sparsity", lp.Astart_[lp.numCol_], lp.Avalue_, true);
  util_analyseMatrixSparsity("Constraint matrix", lp.numCol_, lp.numRow_, lp.Astart_, lp.Aindex_);
  util_analyseModelBounds("Column", lp.numCol_, lp.colLower_, lp.colUpper_);
  util_analyseModelBounds("Row", lp.numRow_, lp.rowLower_, lp.rowUpper_);
}
#endif

HighsStatus convertBasis(const HighsLp& lp, const HighsBasis& basis,
                         HighsBasis_new& new_basis) {
  new_basis.col_status.clear();
  new_basis.row_status.clear();

  new_basis.col_status.resize(lp.numCol_);
  new_basis.row_status.resize(lp.numRow_);

  for (int col = 0; col < lp.numCol_; col++) {
    if (!basis.nonbasicFlag_[col]) {
      new_basis.col_status[col] = HighsBasisStatus::BASIC;
    } else if (basis.nonbasicMove_[col] == NONBASIC_MOVE_UP) {
        new_basis.col_status[col] = HighsBasisStatus::LOWER;
    } else if (basis.nonbasicMove_[col] == NONBASIC_MOVE_DN) {
        new_basis.col_status[col] =  HighsBasisStatus::UPPER;
    } else if (basis.nonbasicMove_[col] == NONBASIC_MOVE_ZE) {
      if (lp.colLower_[col] == lp.colUpper_[col]) {
          new_basis.col_status[col] =  HighsBasisStatus::LOWER;
      } else {
          new_basis.col_status[col] = HighsBasisStatus::ZERO;
      }
    } else {
      return HighsStatus::Error;
    }
  }


  for (int row = 0; row < lp.numRow_; row++) {
    int var = lp.numCol_ + row;
    if (!basis.nonbasicFlag_[var]) {
      new_basis.row_status[row] = HighsBasisStatus::BASIC;
    } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_DN) {
        new_basis.row_status[row] = HighsBasisStatus::LOWER;
    } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_UP) {
        new_basis.row_status[row] = HighsBasisStatus::UPPER;
    } else if (basis.nonbasicMove_[var] == NONBASIC_MOVE_ZE) {
      if (lp.rowLower_[row] == lp.rowUpper_[row]) {
          new_basis.row_status[row] = HighsBasisStatus::LOWER;
      } else {
          new_basis.row_status[row] = HighsBasisStatus::ZERO;
      }
    } else {
      return HighsStatus::Error;
    }
  }

  return HighsStatus::OK;
}

HighsBasis_new getHighsBasis(const HighsLp& lp, const HighsBasis& basis) {
  HighsBasis_new new_basis;
  HighsStatus result = convertBasis(lp, basis, new_basis);
  if (result != HighsStatus::OK)
    return HighsBasis_new();
  // Call Julian's code to translate basis once it's out of
  // SimplexInterface. Until it is out of SimplexInteface use code
  // I just added above which does the same but only returns an
  // error and not which basis index has an illegal value.
  return new_basis;
}

HighsStatus calculateColDuals(const HighsLp& lp, HighsSolution& solution) {
  assert(solution.row_dual.size() > 0);
  if (!isSolutionConsistent(lp, solution))
    return HighsStatus::Error;

  solution.col_dual.assign(lp.numCol_, 0);

  for (int col = 0; col < lp.numCol_; col++) {
    for (int i=lp.Astart_[col]; i<lp.Astart_[col+1]; i++) {
      const int row = lp.Aindex_[i];
      assert(row >= 0);
      assert(row < lp.numRow_);

      solution.col_dual[col] -= solution.row_dual[row] * lp.Avalue_[i];
    }
    solution.col_dual[col] += lp.colCost_[col];
  }

  return HighsStatus::OK;
}

HighsStatus calculateRowValues(const HighsLp& lp, HighsSolution& solution) {
  assert(solution.col_value.size() > 0);
  if (!isSolutionConsistent(lp, solution))
    return HighsStatus::Error;

  solution.row_value.clear();
  solution.row_value.assign(lp.numRow_, 0);

  for (int col = 0; col < lp.numCol_; col++) {
    for (int i=lp.Astart_[col]; i<lp.Astart_[col+1]; i++) {
      const int row = lp.Aindex_[i];
      assert(row >= 0);
      assert(row < lp.numRow_);

      solution.row_value[row] += solution.col_value[col] * lp.Avalue_[i];
    }
  }

  return HighsStatus::OK;
}

HighsStatus assess_interval_set_mask(const int ix_dim, 
				     const bool interval, const int from_ix, const int to_ix,
				     const bool set, int num_set_entries, const int* ix_set,
				     const bool mask, const int* ix_mask,
				     int &from_k, int &to_k) {
  // Check parameter for technique and, if OK, set the loop limits - in iterator style
  if (interval) {
    // Changing by interval: check the parameters and that check set and mask are false
    if (set) {
      HighsLogMessage(HighsMessageType::ERROR, "Index interval and set are both true");
      return HighsStatus::Error;
    }
    if (mask) {
      HighsLogMessage(HighsMessageType::ERROR, "Index interval and mask are both true");
      return HighsStatus::Error;
    }
    if (from_ix < 0) {
      HighsLogMessage(HighsMessageType::ERROR, "Index interval lower limit is %d < 0", from_ix);
      return HighsStatus::Error;
    }
    if (to_ix > ix_dim) {
      HighsLogMessage(HighsMessageType::ERROR, "Index interval upper limit is %d > %d", to_ix, ix_dim);
      return HighsStatus::Error;
    }
    from_k = from_ix;
    to_k = to_ix;
  } else if (set) {
    // Changing by set: check the parameters and check that interval and mask are false
    if (interval) {
      HighsLogMessage(HighsMessageType::ERROR, "Index set and interval are both true");
      return HighsStatus::Error;
    }
    if (mask) {
      HighsLogMessage(HighsMessageType::ERROR, "Index set and mask are both true");
      return HighsStatus::Error;
    }
    if (ix_set == NULL) {
      HighsLogMessage(HighsMessageType::ERROR, "Index set NULL");
      return HighsStatus::Error;
    }
    from_k = 0;
    to_k = num_set_entries;
      // Check that the values in the vector of integers are ascending
    int set_entry_upper = (int)ix_dim-1;
    bool ok = increasing_set_ok(ix_set, num_set_entries, 0, set_entry_upper);
    if (!ok) {
      HighsLogMessage(HighsMessageType::ERROR, "Index set is not ordered");
      return HighsStatus::Error;
    }
  } else if (mask) {
    // Changing by mask: check the parameters and check that set and interval are false
    if (ix_mask == NULL) {
      HighsLogMessage(HighsMessageType::ERROR, "Index mask is NULL");
      return HighsStatus::Error;
    }
    if (interval) {
      HighsLogMessage(HighsMessageType::ERROR, "Index mask and interval are both true");
      return HighsStatus::Error;
    }
    if (set) {
      HighsLogMessage(HighsMessageType::ERROR, "Index mask and set are both true");
      return HighsStatus::Error;
    }
    from_k = 0;
    to_k = ix_dim;
  } else {
    // No method defined
    HighsLogMessage(HighsMessageType::ERROR, "None of index interval, set or mask is true");
    return HighsStatus::Error;
  }
  return HighsStatus::OK;
}

void update_out_in_ix(const int ix_dim, 
		      const bool interval, const int from_ix, const int to_ix,
		      const bool set, const int num_set_entries, const int* ix_set,
		      const bool mask, const int* ix_mask,
		      int& out_from_ix, int& out_to_ix,
		      int& in_from_ix, int& in_to_ix,
		      int& current_set_entry) {
  
  if (interval) {
    out_from_ix = from_ix;
    out_to_ix = to_ix;
    in_from_ix = to_ix;
    in_to_ix = ix_dim;
  } else if (set) {
    out_from_ix = ix_set[current_set_entry];
    out_to_ix = out_from_ix+1;
    current_set_entry++;
    int current_set_entry0 = current_set_entry;
    for (int set_entry = current_set_entry0; set_entry < num_set_entries; set_entry++) {
      int ix = ix_set[set_entry];
      if (ix > out_to_ix) break;
      out_to_ix = ix_set[current_set_entry]+1;
      current_set_entry++;
    }
    in_from_ix = out_to_ix;
    if (current_set_entry < num_set_entries) {
      in_to_ix = ix_set[current_set_entry];
    } else {
      // Account for getting to the end of the set
      in_to_ix = ix_dim;
    }
  } else {
    out_from_ix = in_to_ix;
    out_to_ix = ix_dim;
    for (int ix = in_to_ix; ix < ix_dim; ix++) {
      if (!ix_mask[ix]) {
	out_to_ix = ix;
	break;	
      }
    }
    in_from_ix = out_to_ix;
    in_to_ix = ix_dim;
    for (int ix = out_to_ix; ix < ix_dim; ix++) {
      if (ix_mask[ix]) {
	in_to_ix = ix;
	break;	
      }
    }
  }
}

bool isColDataNull(const double *usr_col_cost, const double *usr_col_lower,  const double *usr_col_upper) {
  bool null_data = false;
  if (usr_col_cost == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied column costs are NULL");
    null_data = true;
  }
  if (usr_col_lower == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied column lower bounds are NULL");
    null_data = true;
  }
  if (usr_col_upper == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied column upper bounds are NULL");
    null_data = true;
  }
  return null_data;
}

bool isRowDataNull(const double *usr_row_lower,  const double *usr_row_upper) {
  bool null_data = false;
  if (usr_row_lower == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied row lower bounds are NULL");
    null_data = true;
  }
  if (usr_row_upper == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied row upper bounds are NULL");
    null_data = true;
  }
  return null_data;
}

bool isMatrixDataNull(const int *usr_matrix_start, const int *usr_matrix_index, const double *usr_matrix_value) {
  bool null_data = false;
  if (usr_matrix_start == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied matrix starts are NULL");
    null_data = true;
  }
  if (usr_matrix_index == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied matrix indices are NULL");
    null_data = true;
  }
  if (usr_matrix_value == NULL) {
    HighsLogMessage(HighsMessageType::ERROR, "User-supplied matrix values are NULL");
    null_data = true;
  }
  return null_data;
}

HighsLp transformIntoEqualityProblem(const HighsLp& lp) {
  assert(checkLp(lp) == HighsStatus::OK);

  // Copy lp.
  HighsLp equality_lp = lp;

  // Add slacks for each row with more than one bound.
  std::vector<double> rhs(lp.numRow_, 0);

  for (int row = 0; row < lp.numRow_; row++) {
    assert(equality_lp.Astart_[equality_lp.numCol_] == equality_lp.Avalue_.size());
    assert(equality_lp.Aindex_.size() == equality_lp.Avalue_.size());
    const int nnz = equality_lp.Astart_[equality_lp.numCol_];

    if (lp.rowLower_[row] == -HIGHS_CONST_INF &&
        lp.rowUpper_[row] == HIGHS_CONST_INF) {
      // free row
      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(-HIGHS_CONST_INF);
      equality_lp.colUpper_.push_back(HIGHS_CONST_INF);
    }
    else if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
             lp.rowUpper_[row] == HIGHS_CONST_INF) {
      // only lower bound
      rhs[row] = lp.rowLower_[row];

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(-1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(HIGHS_CONST_INF);
      equality_lp.colCost_.push_back(0);
    }
    else if (lp.rowLower_[row] == -HIGHS_CONST_INF &&
             lp.rowUpper_[row] < HIGHS_CONST_INF) {
      // only upper bound
      rhs[row] = lp.rowUpper_[row];

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(1.0);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(HIGHS_CONST_INF);
      equality_lp.colCost_.push_back(0);
    }
    else if (lp.rowLower_[row] > -HIGHS_CONST_INF &&
             lp.rowUpper_[row] < HIGHS_CONST_INF &&
             lp.rowLower_[row] != lp.rowUpper_[row]) {
      // both lower and upper bound that are different
      double rhs_value, coefficient;
      double difference = lp.rowUpper_[row] - lp.rowLower_[row];
      if (fabs(lp.rowLower_[row]) < fabs(lp.rowUpper_[row])) {
        rhs_value = lp.rowLower_[row];
        coefficient = -1;
      } else {
        rhs_value = lp.rowUpper_[row];
        coefficient = 1;
      }
      rhs[row] = rhs_value;

      equality_lp.Astart_.push_back(nnz + 1);
      equality_lp.Aindex_.push_back(row);
      equality_lp.Avalue_.push_back(coefficient);

      equality_lp.numCol_++;
      equality_lp.colLower_.push_back(0);
      equality_lp.colUpper_.push_back(difference);
      equality_lp.colCost_.push_back(0);
    }
    else if (lp.rowLower_[row] == lp.rowUpper_[row]) {
      // equality row
      rhs[row] = lp.rowLower_[row];
    } else {
      HighsPrintMessage(ML_ALWAYS,
                        "Unknown row type when adding slacks. \
                         Returning unmodified lp copy.");
      return lp;
    }
  }
  equality_lp.rowLower_ = rhs;
  equality_lp.rowUpper_ = rhs;

  return equality_lp;
}

// Given (P) returns (D) for the pair
// (P)
//    min c'x st Ax=b
//     st l <= x <= u
// (D)
//    max b'y + l'zl - u'zu
//     st A'y + zl - zu = c
//        y free, zl >=0, zu >= 0
HighsLp dualizeEqualityProblem(const HighsLp& lp) {
  assert(checkLp(lp) == HighsStatus::OK);
  assert(lp.sense_ == OBJSENSE_MINIMIZE);
  assert(lp.rowLower_ == lp.rowUpper_);

  HighsLp dual;
  const int ncols = lp.numRow_;
  const int nrows = lp.numCol_;

  dual.numRow_ = nrows;
  dual.rowLower_ = lp.colCost_;
  dual.rowUpper_ = lp.colCost_;

  // Add columns (y)
  dual.numCol_ = ncols;
  dual.colLower_.resize(ncols);
  dual.colUpper_.resize(ncols);
  dual.colCost_.resize(ncols);

  for (int col = 0; col < ncols; col++) {
    dual.colLower_[col] = -HIGHS_CONST_INF;
    dual.colUpper_[col] = HIGHS_CONST_INF;
    // cost b'y
    dual.colCost_[col] = lp.rowLower_[col];
  }

  // Get transpose of A
  int i, k;
  vector<int> iwork(lp.numRow_, 0);
  dual.Astart_.resize(lp.numRow_ + 1, 0);
  int AcountX = lp.Aindex_.size();
  dual.Aindex_.resize(AcountX);
  dual.Avalue_.resize(AcountX);
  for (int k = 0; k < AcountX; k++) iwork.at(lp.Aindex_.at(k))++;
  for (i = 1; i <= lp.numRow_; i++)
    dual.Astart_.at(i) = dual.Astart_.at(i - 1) + iwork.at(i - 1);
  for (i = 0; i < lp.numRow_; i++) iwork.at(i) = dual.Astart_.at(i);
  for (int iCol = 0; iCol < lp.numCol_; iCol++) {
    for (k = lp.Astart_.at(iCol); k < lp.Astart_.at(iCol + 1); k++) {
      int iRow = lp.Aindex_.at(k);
      int iPut = iwork.at(iRow)++;
      dual.Aindex_.at(iPut) = iCol;
      dual.Avalue_.at(iPut) = lp.Avalue_[k];
    }
  }

  // Add columns (zl)
  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colLower_[col] > -HIGHS_CONST_INF) {
      const int nnz = dual.Astart_[dual.numCol_];

      dual.colLower_.push_back(0);
      dual.colUpper_.push_back(HIGHS_CONST_INF);

      dual.colCost_.push_back(lp.colLower_[col]);

      // Add constaints
      dual.Astart_.push_back(nnz + 1);
      dual.Aindex_.push_back(col);
      dual.Avalue_.push_back(1.0);

      dual.numCol_++;
    }
  }

  // Add columns (zu)
  for (int col = 0; col < lp.numCol_; col++) {
    if (lp.colUpper_[col] < HIGHS_CONST_INF) {
      const int nnz = dual.Astart_[dual.numCol_];

      dual.colLower_.push_back(0);
      dual.colUpper_.push_back(HIGHS_CONST_INF);

      dual.colCost_.push_back(-lp.colUpper_[col]);

      // Add constaints
      dual.Astart_.push_back(nnz + 1);
      dual.Aindex_.push_back(col);
      dual.Avalue_.push_back(-1.0);

      dual.numCol_++;
    }
  }

  dual.offset_ = -lp.offset_;
  dual.sense_ = OBJSENSE_MAXIMIZE;

  HighsPrintMessage(ML_ALWAYS, "Dualized equality LP.\n");
  return dual;
}