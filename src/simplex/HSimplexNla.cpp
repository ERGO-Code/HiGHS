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
/**@file simplex/HSimplexNla.cpp
 *
 * @brief Interface to HFactor allowing non-HFactor updates, NLA-only
 * scaling and shifting of NLA analysis below simplex level.
 */
#include "simplex/HSimplexNla.h"

//#include <algorithm>
#include <stdio.h>

#ifdef OPENMP
#include "omp.h"
#endif

// using std::max;
// using std::min;
// using std::vector;

void HSimplexNla::setup(HighsInt num_col, HighsInt num_row,
                        const HighsInt* a_start, const HighsInt* a_index,
                        const double* a_value, HighsInt* base_index,
                        double factor_pivot_threshold, HighsOptions* options,
                        HighsTimer* timer, HighsSimplexAnalysis* analysis) {
  //  printf("In HSimplexNla::setup\n");
  num_col_ = num_col;
  num_row_ = num_row;
  a_start_ = a_start;
  a_index_ = a_index;
  a_value_ = a_value;
  scale_ = NULL;
  base_index_ = base_index;
  options_ = options;
  timer_ = timer;
  analysis_ = analysis;
  factor_.setup(num_col_, num_row_, a_start_, a_index_, a_value_, base_index_,
                factor_pivot_threshold, options_->factor_pivot_tolerance,
                options_->highs_debug_level, options_->output_flag,
                options_->log_file_stream, options_->log_to_console,
                options_->log_dev_level);
}

HighsInt HSimplexNla::invert() {
  //  printf("In HSimplexNla::invert\n");
  HighsTimerClock* factor_timer_clock_pointer = NULL;
  if (analysis_->analyse_factor_time) {
    HighsInt thread_id = 0;
#ifdef OPENMP
    thread_id = omp_get_thread_num();
#endif
    factor_timer_clock_pointer =
        analysis_->getThreadFactorTimerClockPtr(thread_id);
  }
  HighsInt rank_deficiency = factor_.build(factor_timer_clock_pointer);
  build_synthetic_tick_ = factor_.build_synthetic_tick;
  return rank_deficiency;
}

void HSimplexNla::btran(HVector& rhs, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  //  printf("In HSimplexNla::btran\n");
  if (scale_ != NULL) {
    // Scale the RHS
    printf("RHS:\n");
    const vector<double>& col_scale = scale_->col;
    const vector<double>& row_scale = scale_->row;
    HighsInt to_entry;
    const bool use_row_indices = sparseLoopStyle(rhs.count, num_row_, to_entry);
    for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
      HighsInt iCol;
      if (use_row_indices) {
	iCol = rhs.index[iEntry];
      } else {
	iCol = iEntry;
      }
      HighsInt iVar = base_index_[iCol];
      double prev_array_value = rhs.array[iCol];
      if (iVar < num_col_) {
	rhs.array[iCol] *= col_scale[iVar];
	printf("Entry %d - basic variable Col %d - scales by %g from %g to %g\n",
	       (int)iCol, (int)iVar, col_scale[iVar], prev_array_value, rhs.array[iCol]);
      } else {
	rhs.array[iCol] /= row_scale[iVar-num_col_];
	printf("Entry %d - basic variable Row %d - scales by %g from %g to %g\n",
	       (int)iCol, (int)(iVar-num_col_), row_scale[iVar-num_col_], prev_array_value, rhs.array[iCol]);
      }
    }
  }
  factor_.btranCall(rhs, expected_density, factor_timer_clock_pointer);
  if (scale_ != NULL) {
    // Scale the solution
    printf("Solution:\n");
    const vector<double>& col_scale = scale_->col;
    const vector<double>& row_scale = scale_->row;
    HighsInt to_entry;
    const bool use_row_indices = sparseLoopStyle(rhs.count, num_row_, to_entry);
    for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
      HighsInt iRow;
      if (use_row_indices) {
	iRow = rhs.index[iEntry];
      } else {
	iRow = iEntry;
      }
      double prev_array_value = rhs.array[iRow];
      rhs.array[iRow] *= row_scale[iRow];
      printf("Entry %d scales by %g from %g to %g\n",
	     (int)iRow, row_scale[iRow], prev_array_value, rhs.array[iRow]);
    }
  }
}

void HSimplexNla::ftran(HVector& rhs, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  //  printf("In HSimplexNla::ftran\n");
  factor_.ftranCall(rhs, expected_density, factor_timer_clock_pointer);
}

void HSimplexNla::update(HVector* aq, HVector* ep, HighsInt* iRow,
                         HighsInt* hint) {
  //  printf("In HSimplexNla::update\n");
  factor_.update(aq, ep, iRow, hint);
}

void HSimplexNla::setPivotThreshold(const double new_pivot_threshold) {
  factor_.setPivotThreshold(new_pivot_threshold);
}

void HSimplexNla::passScaleAndMatrixPointers(const HighsScale* scale,
					     const HighsInt* a_start, 
					     const HighsInt* a_index, 
					     const double* a_value) {
  scale_ = scale;
  a_start_ = a_start;
  a_index_ = a_index;
  a_value_ = a_value;
  factor_.setupMatrix(a_start, a_index, a_value);
}

bool HSimplexNla::sparseLoopStyle(const HighsInt count, const HighsInt dim,
				  HighsInt& to_entry) const {
  // Parameter to decide whether to use just the values in a HVector, or
  // use the indices of their nonzeros
  const double density_for_indexing = 0.4;
  const bool use_indices = count >= 0 && count < density_for_indexing * dim;
  if (use_indices) {
    to_entry = count;
  } else {
    to_entry = dim;
  }
  return use_indices;
}

