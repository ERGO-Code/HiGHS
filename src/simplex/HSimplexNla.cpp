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

void HSimplexNla::setup(const HighsLp* lp, HighsInt* base_index,
                        const double factor_pivot_threshold,
                        const HighsOptions* options, HighsTimer* timer,
                        HighsSimplexAnalysis* analysis) {
  lp_ = lp;
  scale_ = NULL;
  base_index_ = base_index;
  options_ = options;
  timer_ = timer;
  analysis_ = analysis;
  factor_.setup(lp_->num_col_, lp_->num_row_, &lp_->a_start_[0],
                &lp_->a_index_[0], &lp_->a_value_[0], base_index_,
                factor_pivot_threshold, options_->factor_pivot_tolerance,
                options_->highs_debug_level, options_->output_flag,
                options_->log_file_stream, options_->log_to_console,
                options_->log_dev_level);
}

HighsInt HSimplexNla::invert() {
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
  applyBasisMatrixColScale(rhs);
  factor_.btranCall(rhs, expected_density, factor_timer_clock_pointer);
  applyBasisMatrixRowScale(rhs);
}

void HSimplexNla::ftran(HVector& rhs, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  applyBasisMatrixRowScale(rhs);
  factor_.ftranCall(rhs, expected_density, factor_timer_clock_pointer);
  applyBasisMatrixColScale(rhs);
}

void HSimplexNla::update(HVector* aq, HVector* ep, HighsInt* iRow,
                         HighsInt* hint) {
  reportPackValue("  pack: aq Bf ", aq);
  reportPackValue("  pack: ep Bf ", ep);
  factor_.update(aq, ep, iRow, hint);
}

void HSimplexNla::transformForUpdate(HVector* aq, HVector* ep,
                                     const HighsInt variable_in,
                                     const HighsInt row_out) {
  if (scale_ == NULL) return;
  // For (\hat)aq, UPDATE needs packValue and array[row_out] to
  // correspond to \bar{B}^{-1}(R.aq.cq), but CB.\bar{B}^{-1}(R.aq)
  // has been computed.
  //
  // Hence packValue and array[row_out] need to be scaled by cq;
  //
  // array[row_out] has to be unscaled by the corresponding entry of
  // CB
  //
  reportPackValue("pack aq Bf ", aq);
  double scale_factor;
  if (variable_in < lp_->num_col_) {
    scale_factor = scale_->col[variable_in];
  } else {
    scale_factor = 1.0 / scale_->row[variable_in - lp_->num_col_];
  }
  for (HighsInt ix = 0; ix < aq->packCount; ix++)
    aq->packValue[ix] *= scale_factor;
  reportPackValue("pack aq Af ", aq);
  //
  // Now focus on the pivot value, aq->array[row_out]
  //
  // First scale by cq
  aq->array[row_out] *= scale_factor;
  //
  // Also have to unscale by cp
  HighsInt variable_out = base_index_[row_out];
  if (variable_out < lp_->num_col_) {
    scale_factor = scale_->col[variable_out];
  } else {
    scale_factor = 1.0 / scale_->row[variable_out - lp_->num_col_];
  }
  aq->array[row_out] /= scale_factor;
  // For (\hat)ep, UPDATE needs packValue to correspond to
  // \bar{B}^{-T}ep, but R.\bar{B}^{-T}(CB.ep) has been computed.
  //
  // Hence packValue needs to be unscaled by cp
  for (HighsInt ix = 0; ix < ep->packCount; ix++)
    ep->packValue[ix] /= scale_factor;
}

void HSimplexNla::setPivotThreshold(const double new_pivot_threshold) {
  factor_.setPivotThreshold(new_pivot_threshold);
}

void HSimplexNla::passScaleAndFactorMatrixPointers(
    const SimplexScale* scale, const HighsInt* factor_a_start,
    const HighsInt* factor_a_index, const double* factor_a_value) {
  scale_ = scale;
  factor_.setupMatrix(factor_a_start, factor_a_index, factor_a_value);
}

void HSimplexNla::applyBasisMatrixRowScale(HVector& rhs) const {
  if (scale_ == NULL) return;
  const vector<double>& col_scale = scale_->col;
  const vector<double>& row_scale = scale_->row;
  HighsInt to_entry;
  const bool use_row_indices =
      sparseLoopStyle(rhs.count, lp_->num_row_, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iRow;
    if (use_row_indices) {
      iRow = rhs.index[iEntry];
    } else {
      iRow = iEntry;
    }
    rhs.array[iRow] *= row_scale[iRow];
  }
}

void HSimplexNla::applyBasisMatrixColScale(HVector& rhs) const {
  if (scale_ == NULL) return;
  const vector<double>& col_scale = scale_->col;
  const vector<double>& row_scale = scale_->row;
  HighsInt to_entry;
  const bool use_row_indices =
      sparseLoopStyle(rhs.count, lp_->num_row_, to_entry);
  for (HighsInt iEntry = 0; iEntry < to_entry; iEntry++) {
    HighsInt iCol;
    if (use_row_indices) {
      iCol = rhs.index[iEntry];
    } else {
      iCol = iEntry;
    }
    HighsInt iVar = base_index_[iCol];
    if (iVar < lp_->num_col_) {
      rhs.array[iCol] *= col_scale[iVar];
    } else {
      rhs.array[iCol] /= row_scale[iVar - lp_->num_col_];
    }
  }
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

void HSimplexNla::reportArray(const std::string message,
                              const HVector* vector) const {
  if (!report) return;
  const HighsInt num_row = lp_->num_row_;
  if (num_row > 25) {
    reportArraySparse(message, vector);
  } else {
    printf("%s", message.c_str());
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow % 5 == 0) printf("\n");
      printf("%10.4g ", vector->array[iRow]);
    }
    printf("\n");
  }
}

void HSimplexNla::reportArraySparse(const std::string message,
                                    const HVector* vector) const {
  if (!report) return;
  const HighsInt num_row = lp_->num_row_;
  if (vector->count > 25) return;
  if (vector->count < num_row) {
    printf("%s", message.c_str());
    for (HighsInt en = 0; en < vector->count; en++) {
      HighsInt iRow = vector->index[en];
      if (en % 5 == 0) printf("\n");
      printf("(%4d %10.4g) ", (int)iRow, vector->array[iRow]);
    }
  } else {
    if (num_row > 25) return;
    printf("%s", message.c_str());
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow % 5 == 0) printf("\n");
      printf("%10.4g ", vector->array[iRow]);
    }
  }
  printf("\n");
}

void HSimplexNla::reportPackValue(const std::string message,
                                  const HVector* vector) const {
  if (!report) return;
  const HighsInt num_row = lp_->num_row_;
  if (vector->packCount > 25) return;
  printf("%s", message.c_str());
  for (HighsInt en = 0; en < vector->packCount; en++) {
    HighsInt iRow = vector->packIndex[en];
    if (en % 5 == 0) printf("\n");
    printf("(%4d %10.4g) ", (int)iRow, vector->packValue[en]);
  }
  printf("\n");
}
