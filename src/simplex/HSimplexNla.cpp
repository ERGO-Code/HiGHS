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

#include <stdio.h>

#include "simplex/HSimplex.h"

#ifdef OPENMP
#include "omp.h"
#endif

using std::vector;

void HSimplexNla::setup(const HighsLp* lp, HighsInt* base_index,
                        const HighsOptions* options, HighsTimer* timer,
                        HighsSimplexAnalysis* analysis,
                        const HighsSparseMatrix* factor_a_matrix,
                        const double factor_pivot_threshold) {
  this->setLpAndScalePointers(lp);
  this->base_index_ = base_index;
  this->options_ = options;
  this->timer_ = timer;
  this->analysis_ = analysis;
  this->report_ = false;
  this->factor_.setup(
      this->lp_->num_col_, this->lp_->num_row_, &factor_a_matrix->start_[0],
      &factor_a_matrix->index_[0], &factor_a_matrix->value_[0],
      this->base_index_, factor_pivot_threshold,
      this->options_->factor_pivot_tolerance, this->options_->highs_debug_level,
      this->options_->output_flag, this->options_->log_file_stream,
      this->options_->log_to_console, this->options_->log_dev_level);
  assert(debugCheckData("After HSimplexNla::setup") == HighsDebugStatus::kOk);
}

void HSimplexNla::setLpAndScalePointers(const HighsLp* for_lp) {
  this->lp_ = for_lp;
  this->scale_ = NULL;
  if (for_lp->scale_.has_scaling && !for_lp->is_scaled_)
    this->scale_ = &(for_lp->scale_);
}

void HSimplexNla::setBasicIndexPointers(HighsInt* basic_index) {
  this->base_index_ = basic_index;
  this->factor_.baseIndex = basic_index;
}

void HSimplexNla::setPointers(const HighsLp* for_lp,
                              const HighsSparseMatrix* factor_a_matrix,
                              HighsInt* base_index, const HighsOptions* options,
                              HighsTimer* timer,
                              HighsSimplexAnalysis* analysis) {
  this->setLpAndScalePointers(for_lp);
  if (factor_a_matrix) factor_.setupMatrix(factor_a_matrix);
  if (base_index) base_index_ = base_index;
  if (options) options_ = options;
  if (timer) timer_ = timer;
  if (analysis) analysis_ = analysis;
}

void HSimplexNla::clear() {
  lp_ = NULL;
  scale_ = NULL;
  base_index_ = NULL;
  options_ = NULL;
  timer_ = NULL;
  analysis_ = NULL;
  report_ = false;
  build_synthetic_tick_ = 0;
  this->frozenBasisClearAllData();
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
  // Clear any frozen basis updates
  frozenBasisClearAllUpdate();
  return rank_deficiency;
}

void HSimplexNla::btran(HVector& rhs, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  applyBasisMatrixColScale(rhs);
  frozenBtran(rhs);
  factor_.btranCall(rhs, expected_density, factor_timer_clock_pointer);
  applyBasisMatrixRowScale(rhs);
}

void HSimplexNla::ftran(HVector& rhs, const double expected_density,
                        HighsTimerClock* factor_timer_clock_pointer) const {
  applyBasisMatrixRowScale(rhs);
  factor_.ftranCall(rhs, expected_density, factor_timer_clock_pointer);
  frozenFtran(rhs);
  applyBasisMatrixColScale(rhs);
}

void HSimplexNla::frozenBtran(HVector& rhs) const {
  HighsInt frozen_basis_id = last_frozen_basis_id_;
  if (frozen_basis_id == kNoLink) return;
  // Apply any updates since the last frozen basis
  update_.btran(rhs);
  // Work through any updates associated with previously frozen basis
  frozen_basis_id = frozen_basis_[frozen_basis_id].prev_;
  if (frozen_basis_id == kNoLink) return;
  for (;;) {
    assert(frozen_basis_id != kNoLink);
    const FrozenBasis& frozen_basis = frozen_basis_[frozen_basis_id];
    frozen_basis.update_.btran(rhs);
    frozen_basis_id = frozen_basis.prev_;
    if (frozen_basis_id == kNoLink) break;
  }
}

void HSimplexNla::frozenFtran(HVector& rhs) const {
  // Work through any updates associated with previously frozen basis
  HighsInt frozen_basis_id = first_frozen_basis_id_;
  if (frozen_basis_id == kNoLink) return;
  for (;;) {
    assert(frozen_basis_id != kNoLink);
    if (frozen_basis_id == last_frozen_basis_id_) break;
    const FrozenBasis& frozen_basis = frozen_basis_[frozen_basis_id];
    frozen_basis.update_.ftran(rhs);
    frozen_basis_id = frozen_basis.next_;
  }
  // Now apply any updates since the last frozen basis
  update_.ftran(rhs);
}

void HSimplexNla::update(HVector* aq, HVector* ep, HighsInt* iRow,
                         HighsInt* hint) {
  reportPackValue("  pack: aq Bf ", aq);
  reportPackValue("  pack: ep Bf ", ep);
  factor_.refactor_info_.clear();
  if (last_frozen_basis_id_ == kNoLink) {
    factor_.update(aq, ep, iRow, hint);
  } else {
    *hint = update_.update(aq, iRow);
  }
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

void HSimplexNla::addCols(const HighsLp* updated_lp) {
  // Adding columns is easy, since they are nonbasic
  //
  // Set the pointers for the LP and scaling. The pointer to the
  // vector of basic variables isn't updated, since it hasn't been
  // resized. The HFactor matrix isn't needed until reinversion has to
  // be performed
  setLpAndScalePointers(updated_lp);
}

void HSimplexNla::addRows(const HighsLp* updated_lp, HighsInt* base_index,
                          const HighsSparseMatrix* scaled_ar_matrix) {
  // Adding rows is not so easy, since their slacks are basic
  //
  // Set the pointers for the LP, scaling and basic variables. The
  // HFactor matrix isn't needed until reinversion has to be performed
  setLpAndScalePointers(updated_lp);
  base_index_ = base_index;
  factor_.baseIndex = base_index;
  factor_.addRows(scaled_ar_matrix);
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

void HSimplexNla::reportArray(const std::string message, const HVector* vector,
                              const bool force) const {
  if (!report_ && !force) return;
  const HighsInt num_row = lp_->num_row_;
  if (num_row > 25) {
    reportArraySparse(message, vector, force);
  } else {
    printf("%s", message.c_str());
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow > 0 && iRow % 10 == 0)
        printf("\n                                 ");
      printf("%11.4g ", vector->array[iRow]);
    }
    printf("\n");
  }
}

void HSimplexNla::reportArraySparse(const std::string message,
                                    const HVector* vector,
                                    const bool force) const {
  if (!report_ && !force) return;
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
                                  const HVector* vector,
                                  const bool force) const {
  if (!report_ && !force) return;
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

HighsDebugStatus HSimplexNla::debugCheckData(const std::string message) const {
  std::string scale_status;
  if (scale_ == NULL) {
    scale_status = "NULL";
  } else {
    scale_status = "non-NULL";
  }
  //  if (options_->highs_debug_level < kHighsDebugLevelCheap) return
  //  HighsDebugStatus::kOk;
  HighsLp check_lp = *lp_;
  bool error0_found = false;
  bool error1_found = false;
  bool error2_found = false;
  bool error_found = false;
  const HighsInt* factor_Astart = factor_.getAstart();
  const HighsInt* factor_Aindex = factor_.getAindex();
  const double* factor_Avalue = factor_.getAvalue();

  if (scale_ == NULL) {
    if (factor_Astart != &(lp_->a_matrix_.start_[0])) error0_found = true;
    if (factor_Aindex != &(lp_->a_matrix_.index_[0])) error1_found = true;
    if (factor_Avalue != &(lp_->a_matrix_.value_[0])) error2_found = true;
    error_found = error0_found || error1_found || error2_found;
    if (error_found) {
      highsLogUser(options_->log_options, HighsLogType::kError,
                   "CheckNlaData: (%s) scale_ is %s lp_ - factor_ matrix "
                   "pointer errors\n",
                   message.c_str(), scale_status.c_str());
      if (error0_found)
        printf("a_matrix_.start_ pointer error: %p vs %p\n",
               (void*)factor_Astart, (void*)&(lp_->a_matrix_.start_[0]));
      if (error1_found) printf("a_matrix_.index pointer error\n");
      if (error2_found) printf("a_matrix_.value pointer error\n");
      assert(!error_found);
      return HighsDebugStatus::kLogicalError;
    }
  } else {
    check_lp.applyScale();
  }
  HighsInt error_col = -1;
  for (HighsInt iCol = 0; iCol < check_lp.num_col_ + 1; iCol++) {
    if (check_lp.a_matrix_.start_[iCol] != factor_Astart[iCol]) {
      error_col = iCol;
      break;
    }
  }
  error_found = error_col >= 0;
  if (error_found) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "CheckNlaData: (%s) scale_ is %s check_lp.a_matrix_.start_ != "
                 "factor_Astart for col %d (%d != %d)\n",
                 message.c_str(), scale_status.c_str(), (int)error_col,
                 (int)check_lp.a_matrix_.start_[error_col],
                 (int)factor_Astart[error_col]);
    assert(!error_found);
    return HighsDebugStatus::kLogicalError;
  }
  HighsInt nnz = check_lp.a_matrix_.numNz();
  HighsInt error_el = -1;
  for (HighsInt iEl = 0; iEl < nnz; iEl++) {
    if (check_lp.a_matrix_.index_[iEl] != factor_Aindex[iEl]) {
      error_el = iEl;
      break;
    }
  }
  error_found = error_el >= 0;
  if (error_found) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "CheckNlaData: (%s) scale_ is %s check_lp.a_matrix_.index_ != "
                 "factor_Aindex for el %d (%d != %d)\n",
                 message.c_str(), scale_status.c_str(), (int)error_el,
                 (int)check_lp.a_matrix_.index_[error_el],
                 (int)factor_Aindex[error_el]);
    assert(!error_found);
    return HighsDebugStatus::kLogicalError;
  }
  for (HighsInt iEl = 0; iEl < nnz; iEl++) {
    if (check_lp.a_matrix_.value_[iEl] != factor_Avalue[iEl]) {
      error_el = iEl;
      break;
    }
  }
  error_found = error_el >= 0;
  if (error_found) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "CheckNlaData: (%s) scale_ is %s check_lp.a_matrix_.value_ != "
                 "factor_Avalue for el %d (%g != %g)\n",
                 message.c_str(), scale_status.c_str(), (int)error_el,
                 check_lp.a_matrix_.value_[error_el], factor_Avalue[error_el]);
    assert(!error_found);
    return HighsDebugStatus::kLogicalError;
  }
  return HighsDebugStatus::kOk;
}
