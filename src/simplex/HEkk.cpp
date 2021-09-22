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
/**@file simplex/HEkk.cpp
 * @brief
 */
#include "simplex/HEkk.h"

#include "lp_data/HighsLpSolverObject.h"
#include "lp_data/HighsLpUtils.h"
#include "lp_data/HighsModelUtils.h"
#include "lp_data/HighsSolutionDebug.h"
#include "simplex/HEkkDual.h"
#include "simplex/HEkkPrimal.h"
#include "simplex/HSimplexDebug.h"
#include "simplex/HSimplexReport.h"
#include "simplex/SimplexTimer.h"

#ifdef OPENMP
#include "omp.h"
#endif

using std::fabs;
using std::max;
using std::min;

// using std::cout;
// using std::endl;

void HEkk::clear() {
  // Clears Ekk entirely. Clears all associated pointers, data scalars
  // and vectors, and the status values.
  this->clearEkkLp();
  this->clearEkkDualise();
  this->clearEkkData();
  this->clearEkkPointers();
  this->clearSimplexBasis(this->basis_);
  this->simplex_nla_.clear();
  this->clearEkkAllStatus();
}

void HEkk::clearEkkAllStatus() {
  // Clears the Ekk status entirely. Functionally junks all
  // information relating to the simplex solve, but doesn't clear the
  // associated data scalars and vectors.
  HighsSimplexStatus& status = this->status_;
  status.initialised_for_new_lp = false;
  status.initialised_for_solve = false;
  this->clearNlaStatus();
  this->clearEkkDataStatus();
}

void HEkk::clearEkkDataStatus() {
  // Just clears the Ekk status values associated with Ekk-specific
  // data: doesn't clear "initialised_for_new_lp", "initialised_for_solve" or
  // NLA status
  HighsSimplexStatus& status = this->status_;
  status.has_ar_matrix = false;
  status.has_dual_steepest_edge_weights = false;
  status.has_fresh_rebuild = false;
  status.has_dual_objective_value = false;
  status.has_primal_objective_value = false;
  status.has_dual_ray = false;
  status.has_primal_ray = false;
}

void HEkk::clearNlaStatus() {
  // Clears Ekk status values associated with NLA. Functionally junks
  // NLA, but doesn't clear the associated data scalars and vectors
  HighsSimplexStatus& status = this->status_;
  status.has_basis = false;
  status.has_nla = false;
  clearNlaInvertStatus();
}

void HEkk::clearNlaInvertStatus() {
  this->status_.has_invert = false;
  this->status_.has_fresh_invert = false;
}

void HEkk::clearEkkPointers() {
  this->options_ = NULL;
  this->timer_ = NULL;
}

void HEkk::clearEkkLp() {
  this->lp_.clear();
  lp_name_ = "";
}

void HEkk::clearEkkDualise() {
  this->original_col_cost_.clear();
  this->original_col_lower_.clear();
  this->original_col_upper_.clear();
  this->original_row_lower_.clear();
  this->original_row_upper_.clear();
  this->upper_bound_col_.clear();
  this->upper_bound_row_.clear();
}

void HEkk::clearEkkData() {
  // Clears Ekk-specific data scalars and vectors. Doesn't clear
  // status, as this is done elsewhere
  //
  // Doesn't clear the LP, simplex NLA or simplex basis as part of
  // clearing Ekk data, so that the simplex basis and HFactor instance
  // are maintained
  //
  // Does clear any frozen basis data
  if (this->status_.has_nla) this->simplex_nla_.frozenBasisClearAllData();

  // analysis_; No clear yet

  this->clearEkkDataInfo();
  model_status_ = HighsModelStatus::kNotset;
  // random_; Has no data
  this->workEdWt_ = NULL;
  this->workEdWtFull_ = NULL;

  this->simplex_in_scaled_space_ = false;
  this->ar_matrix_.clear();
  this->scaled_a_matrix_.clear();

  this->cost_scale_ = 1;
  this->iteration_count_ = 0;
  this->dual_simplex_cleanup_level_ = 0;
  this->dual_simplex_phase1_cleanup_level_ = 0;

  this->solve_bailout_ = false;
  this->called_return_from_solve_ = false;
  this->exit_algorithm_ = SimplexAlgorithm::kPrimal;
  this->return_primal_solution_status_ = 0;
  this->return_dual_solution_status_ = 0;

  this->build_synthetic_tick_ = 0.0;
  this->total_synthetic_tick_ = 0.0;
}

void HEkk::clearEkkDataInfo() {
  HighsSimplexInfo& info = this->info_;
  info.workCost_.clear();
  info.workDual_.clear();
  info.workShift_.clear();
  info.workLower_.clear();
  info.workUpper_.clear();
  info.workRange_.clear();
  info.workValue_.clear();
  info.workLowerShift_.clear();
  info.workUpperShift_.clear();
  info.baseLower_.clear();
  info.baseUpper_.clear();
  info.baseValue_.clear();
  info.numTotRandomValue_.clear();
  info.numTotPermutation_.clear();
  info.numColPermutation_.clear();
  info.devex_index_.clear();
  info.phase1_backtracking_test_done = false;
  info.phase2_backtracking_test_done = false;
  info.backtracking_ = false;
  info.valid_backtracking_basis_ = false;
  this->clearSimplexBasis(info.backtracking_basis_);
  info.backtracking_basis_costs_perturbed_ = 0;
  info.backtracking_basis_bounds_perturbed_ = 0;
  info.backtracking_basis_workShift_.clear();
  info.backtracking_basis_workLowerShift_.clear();
  info.backtracking_basis_workUpperShift_.clear();
  info.backtracking_basis_edge_weights_.clear();
  info.dual_ray_row_ = -1;
  info.dual_ray_sign_ = 0;
  info.primal_ray_col_ = -1;
  info.primal_ray_sign_ = 0;
  info.simplex_strategy = 0;
  info.dual_edge_weight_strategy = 0;
  info.primal_edge_weight_strategy = 0;
  info.price_strategy = 0;
  info.dual_simplex_cost_perturbation_multiplier = 1;
  info.primal_simplex_phase1_cost_perturbation_multiplier = 1;
  info.primal_simplex_bound_perturbation_multiplier = 1;
  info.allow_dual_steepest_edge_to_devex_switch = 0;
  info.dual_steepest_edge_weight_log_error_threshold = 0;
  info.run_quiet = false;
  info.store_squared_primal_infeasibility = false;
  info.report_simplex_inner_clock = false;
  info.report_simplex_outer_clock = false;
  info.report_simplex_phases_clock = false;
  info.report_HFactor_clock = false;
  info.analyse_lp = false;
  info.analyse_iterations = false;
  info.analyse_invert_form = false;
  info.allow_cost_perturbation = true;
  info.allow_bound_perturbation = true;
  info.costs_perturbed = false;
  info.bounds_perturbed = false;
  info.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  info.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  info.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  info.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
  info.dual_phase1_iteration_count = 0;
  info.dual_phase2_iteration_count = 0;
  info.primal_phase1_iteration_count = 0;
  info.primal_phase2_iteration_count = 0;
  info.primal_bound_swap = 0;
  info.min_threads = 1;
  info.num_threads = 1;
  info.max_threads = kHighsThreadLimit;
  info.multi_iteration = 0;
  info.update_count = 0;
  info.dual_objective_value = 0;
  info.primal_objective_value = 0;
  info.updated_dual_objective_value = 0;
  info.updated_primal_objective_value = 0;
  info.num_basic_logicals = 0;
}

void HEkk::clearEkkControlInfo() {
  HighsSimplexInfo& info = this->info_;
  info.control_iteration_count0 = 0;
  info.col_aq_density = 0.0;
  info.row_ep_density = 0.0;
  info.row_ap_density = 0.0;
  info.row_DSE_density = 0.0;
  info.col_basic_feasibility_change_density = 0.0;
  info.row_basic_feasibility_change_density = 0.0;
  info.col_BFRT_density = 0.0;
  info.primal_col_density = 0.0;
  info.dual_col_density = 0.0;
  info.costly_DSE_frequency = 0;
  info.num_costly_DSE_iteration = 0;
  info.costly_DSE_measure = 0;
  info.average_log_low_DSE_weight_error = 0;
  info.average_log_high_DSE_weight_error = 0;
}

void HEkk::clearEkkNlaInfo() {
  HighsSimplexInfo& info = this->info_;
  info.factor_pivot_threshold = 0;
  info.update_limit = 0;
}

void HEkk::clearSimplexBasis(SimplexBasis& simplex_basis) {
  // Clear a given simplex basis. Needs an argument because it's used
  // to clear the backtracking basis
  simplex_basis.basicIndex_.clear();
  simplex_basis.nonbasicFlag_.clear();
  simplex_basis.nonbasicMove_.clear();
}

void HotStart::clear() {
  this->valid = false;
  this->refactor_info.clear();
  this->nonbasicMove.clear();
}

void HEkk::clearHotStart() {
  this->hot_start_.clear();
  this->simplex_nla_.factor_.refactor_info_.clear();
}

void HEkk::invalidate() {
  this->status_.initialised_for_new_lp = false;
  assert(!this->status_.is_dualised);
  assert(!this->status_.is_permuted);
  this->status_.initialised_for_solve = false;
  this->invalidateBasisMatrix();
}

void HEkk::invalidateBasisMatrix() {
  // When the constraint matrix changes - dimensions or just (basic)
  // values, the simplex NLA becomes invalid, the simplex basis is
  // no longer valid, and
  this->status_.has_nla = false;
  invalidateBasis();
}

void HEkk::invalidateBasis() {
  // Invalidate the basis of the simplex LP, and all its other
  // basis-related properties
  this->status_.has_basis = false;
  this->invalidateBasisArtifacts();
}

void HEkk::invalidateBasisArtifacts() {
  // Invalidate the artifacts of the basis of the simplex LP
  this->status_.has_ar_matrix = false;
  this->status_.has_dual_steepest_edge_weights = false;
  this->status_.has_invert = false;
  this->status_.has_fresh_invert = false;
  this->status_.has_fresh_rebuild = false;
  this->status_.has_dual_objective_value = false;
  this->status_.has_primal_objective_value = false;
  this->status_.has_dual_ray = false;
  this->status_.has_primal_ray = false;
}

void HEkk::updateStatus(LpAction action) {
  assert(!this->status_.is_dualised);
  assert(!this->status_.is_permuted);
  switch (action) {
    case LpAction::kScale:
      this->invalidateBasisMatrix();
      this->clearHotStart();
      break;
    case LpAction::kNewCosts:
      this->status_.has_fresh_rebuild = false;
      this->status_.has_dual_objective_value = false;
      this->status_.has_primal_objective_value = false;
      break;
    case LpAction::kNewBounds:
      this->status_.has_fresh_rebuild = false;
      this->status_.has_dual_objective_value = false;
      this->status_.has_primal_objective_value = false;
      break;
    case LpAction::kNewBasis:
      this->invalidateBasis();
      this->clearHotStart();
      break;
    case LpAction::kNewCols:
      this->clear();
      this->clearHotStart();
      //    this->invalidateBasisArtifacts();
      break;
    case LpAction::kNewRows:
      if (kExtendInvertWhenAddingRows) {
        // Just clear Ekk data
        this->clearEkkData();
      } else {
        // Clear everything
        this->clear();
      }
      this->clearHotStart();
      //    this->invalidateBasisArtifacts();
      break;
    case LpAction::kDelCols:
      this->clear();
      this->clearHotStart();
      //    this->invalidateBasis();
      break;
    case LpAction::kDelNonbasicCols:
      this->clear();
      this->clearHotStart();
      //    this->invalidateBasis();
      break;
    case LpAction::kDelRows:
      this->clear();
      this->clearHotStart();
      //   this->invalidateBasis();
      break;
    case LpAction::kDelRowsBasisOk:
      assert(1 == 0);
      this->clearHotStart();
      //      info.lp_ = true;
      break;
    case LpAction::kScaledCol:
      this->invalidateBasisMatrix();
      this->clearHotStart();
      break;
    case LpAction::kScaledRow:
      this->invalidateBasisMatrix();
      this->clearHotStart();
      break;
    case LpAction::kHotStart:
      this->clearEkkData();  //
      this->clearNlaInvertStatus();
      break;
    case LpAction::kBacktracking:
      this->status_.has_ar_matrix = false;
      this->status_.has_fresh_rebuild = false;
      this->status_.has_dual_objective_value = false;
      this->status_.has_primal_objective_value = false;
      break;
    default:
      break;
  }
}

void HEkk::setNlaPointersForLpAndScale(const HighsLp& lp) {
  assert(status_.has_nla);
  simplex_nla_.setLpAndScalePointers(&lp);
}

void HEkk::setNlaPointersForTrans(const HighsLp& lp) {
  assert(status_.has_nla);
  assert(status_.has_basis);
  simplex_nla_.setLpAndScalePointers(&lp);
  simplex_nla_.base_index_ = &basis_.basicIndex_[0];
}

void HEkk::setNlaRefactorInfo() {
  simplex_nla_.factor_.refactor_info_ = this->hot_start_.refactor_info;
  simplex_nla_.factor_.refactor_info_.use = true;
}

void HEkk::btran(HVector& rhs, const double expected_density) {
  assert(status_.has_nla);
  simplex_nla_.btran(rhs, expected_density);
}

void HEkk::ftran(HVector& rhs, const double expected_density) {
  assert(status_.has_nla);
  simplex_nla_.ftran(rhs, expected_density);
}

void HEkk::moveLp(HighsLpSolverObject& solver_object) {
  // Move the incumbent LP to EKK
  HighsLp& incumbent_lp = solver_object.lp_;
  this->lp_ = std::move(incumbent_lp);
  incumbent_lp.is_moved_ = true;
  //
  // Invalidate the row-wise matrix
  this->status_.has_ar_matrix = false;
  //
  // The simplex algorithm runs in the same space as the LP that has
  // just been moved in. This is a scaled space if the LP is scaled.
  this->simplex_in_scaled_space_ = this->lp_.is_scaled_;
  //
  // Update other EKK pointers. Currently just pointers to the
  // HighsOptions and HighsTimer members of the Highs class that are
  // communicated by reference via the HighsLpSolverObject instance.
  this->setPointers(&solver_object.options_, &solver_object.timer_);
  // Initialise Ekk if this has not been done. Ekk isn't initialised
  // if moveLp hasn't been called for this instance of HiGHS, or if
  // the Ekk instance is junked due to removing rows from the LP
  this->initialiseEkk();
}

void HEkk::setPointers(HighsOptions* options, HighsTimer* timer) {
  this->options_ = options;
  this->timer_ = timer;
  this->analysis_.timer_ = this->timer_;
}

HighsSparseMatrix* HEkk::getScaledAMatrixPointer() {
  // Return a pointer to either the constraint matrix or a scaled copy
  // (that is a member of the HEkk class), with the latter returned if
  // the LP has scaling factors but is unscaled.
  HighsSparseMatrix* local_scaled_a_matrix = &(this->lp_.a_matrix_);
  if (this->lp_.scale_.has_scaling && !this->lp_.is_scaled_) {
    scaled_a_matrix_ = this->lp_.a_matrix_;
    scaled_a_matrix_.applyScale(this->lp_.scale_);
    local_scaled_a_matrix = &scaled_a_matrix_;
  }
  return local_scaled_a_matrix;
}

HighsStatus HEkk::dualise() {
  assert(lp_.a_matrix_.isColwise());
  original_num_col_ = lp_.num_col_;
  original_num_row_ = lp_.num_row_;
  original_num_nz_ = lp_.a_matrix_.numNz();
  original_offset_ = lp_.offset_;
  original_col_cost_ = lp_.col_cost_;
  original_col_lower_ = lp_.col_lower_;
  original_col_upper_ = lp_.col_upper_;
  original_row_lower_ = lp_.row_lower_;
  original_row_upper_ = lp_.row_upper_;
  // Reserve space for simple dual
  lp_.col_cost_.reserve(original_num_row_);
  lp_.col_lower_.reserve(original_num_row_);
  lp_.col_upper_.reserve(original_num_row_);
  lp_.row_lower_.reserve(original_num_col_);
  lp_.row_upper_.reserve(original_num_col_);
  // Invalidate the original data
  lp_.col_cost_.resize(0);
  lp_.col_lower_.resize(0);
  lp_.col_upper_.resize(0);
  lp_.row_lower_.resize(0);
  lp_.row_upper_.resize(0);
  // The bulk of the constraint matrix of the dual LP is the transpose
  // of the primal constraint matrix. This is obtained row-wise by
  // copying the matrix and flipping the dimensions
  HighsSparseMatrix dual_matrix = lp_.a_matrix_;
  dual_matrix.num_row_ = original_num_col_;
  dual_matrix.num_col_ = original_num_row_;
  dual_matrix.format_ = MatrixFormat::kRowwise;
  // The primal_bound_value vector accumulates the values of the
  // finite bounds on variables - or zero for a free variable - used
  // later to compute the offset and shift for the costs. Many of
  // these components will be zero - all for the case x>=0 - so
  // maintain a list of the nonzeros and corresponding indices. Don't
  // reserve space since they may not be needed
  vector<double> primal_bound_value;
  vector<HighsInt> primal_bound_index;
  const double inf = kHighsInf;
  for (HighsInt iCol = 0; iCol < original_num_col_; iCol++) {
    const double cost = original_col_cost_[iCol];
    const double lower = original_col_lower_[iCol];
    const double upper = original_col_upper_[iCol];
    double primal_bound = inf;
    double row_lower = inf;
    double row_upper = -inf;
    if (lower == upper) {
      // Fixed
      primal_bound = lower;
      // Dual activity a^Ty is free, implying dual for primal column
      // (slack for dual row) is free
      row_lower = -inf;
      row_upper = inf;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Treat as lower
        primal_bound = lower;
        // Dual activity a^Ty is bounded above by cost, implying dual
        // for primal column (slack for dual row) is non-negative
        row_lower = -inf;
        row_upper = cost;
        // Treat upper bound as additional constraint
        upper_bound_col_.push_back(iCol);
      } else {
        // Lower (since upper bound is infinite)
        primal_bound = lower;
        // Dual activity a^Ty is bounded above by cost, implying dual
        // for primal column (slack for dual row) is non-negative
        row_lower = -inf;
        row_upper = cost;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      primal_bound = upper;
      // Dual activity a^Ty is bounded below by cost, implying dual
      // for primal column (slack for dual row) is non-positive
      row_lower = cost;
      row_upper = inf;
    } else {
      // FREE
      //
      // Dual activity a^Ty is fixed by cost, implying dual for primal
      // column (slack for dual row) is fixed at zero
      primal_bound = 0;
      row_lower = cost;
      row_upper = cost;
    }
    assert(row_lower < inf);
    assert(row_upper > -inf);
    assert(primal_bound < inf);
    lp_.row_lower_.push_back(row_lower);
    lp_.row_upper_.push_back(row_upper);
    if (primal_bound) {
      primal_bound_value.push_back(primal_bound);
      primal_bound_index.push_back(iCol);
    }
  }
  for (HighsInt iRow = 0; iRow < original_num_row_; iRow++) {
    double lower = original_row_lower_[iRow];
    double upper = original_row_upper_[iRow];
    double col_cost = inf;
    double col_lower = inf;
    double col_upper = -inf;
    if (lower == upper) {
      // Equality constraint
      //
      // Dual variable has primal RHS as cost and is free
      col_cost = lower;
      col_lower = -inf;
      col_upper = inf;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Treat as lower
        col_cost = lower;
        col_lower = 0;
        col_upper = inf;
        // Treat upper bound as additional constraint
        upper_bound_row_.push_back(iRow);
      } else {
        // Lower (since upper bound is infinite)
        col_cost = lower;
        col_lower = 0;
        col_upper = inf;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      col_cost = upper;
      col_lower = -inf;
      col_upper = 0;
    } else {
      // FREE
      // Shouldn't get free rows, but handle them anyway
      col_cost = 0;
      col_lower = 0;
      col_upper = 0;
    }
    assert(col_lower < inf);
    assert(col_upper > -inf);
    assert(col_cost < inf);
    lp_.col_cost_.push_back(col_cost);
    lp_.col_lower_.push_back(col_lower);
    lp_.col_upper_.push_back(col_upper);
  }
  vector<HighsInt>& start = lp_.a_matrix_.start_;
  vector<HighsInt>& index = lp_.a_matrix_.index_;
  vector<double>& value = lp_.a_matrix_.value_;
  // Boxed variables and constraints yield extra columns in the dual LP
  HighsSparseMatrix extra_columns;
  extra_columns.ensureColwise();
  extra_columns.num_row_ = original_num_col_;
  HighsInt num_upper_bound_col = upper_bound_col_.size();
  HighsInt num_upper_bound_row = upper_bound_row_.size();
  HighsInt num_extra_col = 0;
  double one = 1;
  for (HighsInt iX = 0; iX < num_upper_bound_col; iX++) {
    HighsInt iCol = upper_bound_col_[iX];
    const double lower = original_col_lower_[iCol];
    const double upper = original_col_upper_[iCol];
    extra_columns.addVec(1, &iCol, &one);
    lp_.col_cost_.push_back(upper);
    lp_.col_lower_.push_back(-inf);
    lp_.col_upper_.push_back(0);
    num_extra_col++;
  }

  if (num_upper_bound_row) {
    // Need to identify the submatrix of constraint matrix rows
    // corresponding to those with a row index in
    // upper_bound_row_. When identifying numbers of entries in each
    // row of submatrix, use indirection to get corresponding row
    // index, with a dummy row for rows not in the submatrix.
    HighsInt dummy_row = num_upper_bound_row;
    vector<HighsInt> indirection;
    vector<HighsInt> count;
    indirection.assign(original_num_row_, dummy_row);
    count.assign(num_upper_bound_row + 1, 0);
    HighsInt extra_iRow = 0;
    for (HighsInt iX = 0; iX < num_upper_bound_row; iX++) {
      HighsInt iRow = upper_bound_row_[iX];
      indirection[iRow] = extra_iRow++;
      double upper = original_row_upper_[iRow];
      lp_.col_cost_.push_back(upper);
      lp_.col_lower_.push_back(-inf);
      lp_.col_upper_.push_back(0);
    }
    for (HighsInt iEl = 0; iEl < original_num_nz_; iEl++)
      count[indirection[index[iEl]]]++;
    extra_columns.start_.resize(num_upper_bound_col + num_upper_bound_row + 1);
    for (HighsInt iRow = 0; iRow < num_upper_bound_row; iRow++) {
      extra_columns.start_[num_upper_bound_col + iRow + 1] =
          extra_columns.start_[num_upper_bound_col + iRow] + count[iRow];
      count[iRow] = extra_columns.start_[num_upper_bound_col + iRow];
    }
    HighsInt extra_columns_num_nz =
        extra_columns.start_[num_upper_bound_col + num_upper_bound_row];
    extra_columns.index_.resize(extra_columns_num_nz);
    extra_columns.value_.resize(extra_columns_num_nz);
    for (HighsInt iCol = 0; iCol < original_num_col_; iCol++) {
      for (HighsInt iEl = start[iCol]; iEl < start[iCol + 1]; iEl++) {
        HighsInt iRow = indirection[index[iEl]];
        if (iRow < num_upper_bound_row) {
          HighsInt extra_columns_iEl = count[iRow];
          assert(extra_columns_iEl < extra_columns_num_nz);
          extra_columns.index_[extra_columns_iEl] = iCol;
          extra_columns.value_[extra_columns_iEl] = value[iEl];
          count[iRow]++;
        }
      }
    }
    extra_columns.num_col_ += num_upper_bound_row;
  }
  // Incorporate the cost shift by subtracting A*primal_bound from the
  // cost vector; compute the objective offset
  double delta_offset = 0;
  for (HighsInt iX = 0; iX < primal_bound_index.size(); iX++) {
    HighsInt iCol = primal_bound_index[iX];
    double multiplier = primal_bound_value[iX];
    delta_offset += multiplier * original_col_cost_[iCol];
    for (HighsInt iEl = start[iCol]; iEl < start[iCol + 1]; iEl++)
      lp_.col_cost_[index[iEl]] -= multiplier * value[iEl];
  }
  if (extra_columns.num_col_) {
    // Incorporate the cost shift by subtracting
    // extra_columns*primal_bound from the cost vector for the extra
    // dual variables
    //
    // Have to scatter the packed primal bound values into a
    // full-length vector
    //
    // ToDo Make this more efficient?
    vector<double> primal_bound;
    primal_bound.assign(original_num_col_, 0);
    for (HighsInt iX = 0; iX < primal_bound_index.size(); iX++)
      primal_bound[primal_bound_index[iX]] = primal_bound_value[iX];

    for (HighsInt iCol = 0; iCol < extra_columns.num_col_; iCol++) {
      double cost = lp_.col_cost_[original_num_row_ + iCol];
      for (HighsInt iEl = extra_columns.start_[iCol];
           iEl < extra_columns.start_[iCol + 1]; iEl++)
        cost -=
            primal_bound[extra_columns.index_[iEl]] * extra_columns.value_[iEl];
      lp_.col_cost_[original_num_row_ + iCol] = cost;
    }
  }
  lp_.offset_ += delta_offset;
  // Copy the row-wise dual LP constraint matrix and transpose it.
  // ToDo Make this more efficient
  lp_.a_matrix_ = dual_matrix;
  lp_.a_matrix_.ensureColwise();
  // Add the extra columns to the dual LP constraint matrix
  lp_.a_matrix_.addCols(extra_columns);

  HighsInt dual_num_col =
      original_num_row_ + num_upper_bound_col + num_upper_bound_row;
  HighsInt dual_num_row = original_num_col_;
  assert(dual_num_col == (int)lp_.col_cost_.size());
  assert(lp_.a_matrix_.num_col_ == dual_num_col);
  const bool ignore_scaling = true;
  if (!ignore_scaling) {
    // Flip any scale factors
    if (lp_.scale_.has_scaling) {
      std::vector<double> temp_scale = lp_.scale_.row;
      lp_.scale_.row = lp_.scale_.col;
      lp_.scale_.col = temp_scale;
      lp_.scale_.num_col = dual_num_col;
      lp_.scale_.num_row = dual_num_row;
    }
  }
  // Change optimzation sense
  if (lp_.sense_ == ObjSense::kMinimize) {
    lp_.sense_ = ObjSense::kMaximize;
  } else {
    lp_.sense_ = ObjSense::kMinimize;
  }
  // Flip LP dimensions
  lp_.num_col_ = dual_num_col;
  lp_.num_row_ = dual_num_row;
  status_.is_dualised = true;
  status_.has_basis = false;
  status_.has_ar_matrix = false;
  status_.has_nla = false;
  highsLogUser(options_->log_options, HighsLogType::kInfo,
               "Solving dual LP with %d columns", (int)dual_num_col);
  if (num_upper_bound_col + num_upper_bound_row) {
    highsLogUser(options_->log_options, HighsLogType::kInfo, " [%d extra from",
                 (int)dual_num_col - original_num_row_);
    if (num_upper_bound_col)
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   " %d boxed variable(s)", (int)num_upper_bound_col);
    if (num_upper_bound_col && num_upper_bound_row)
      highsLogUser(options_->log_options, HighsLogType::kInfo, " and");
    if (num_upper_bound_row)
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   " %d boxed constraint(s)", (int)num_upper_bound_row);
    highsLogUser(options_->log_options, HighsLogType::kInfo, "]");
  }
  highsLogUser(options_->log_options, HighsLogType::kInfo, " and %d rows\n",
               (int)dual_num_row);
  //  reportLp(options_->log_options, lp_, HighsLogType::kVerbose);
  return HighsStatus::kOk;
}

HighsStatus HEkk::undualise() {
  if (!this->status_.is_dualised) return HighsStatus::kOk;
  HighsInt dual_num_col = lp_.num_col_;
  HighsInt primal_num_tot = original_num_col_ + original_num_row_;
  // These two aren't used (yet)
  vector<double>& dual_work_dual = info_.workDual_;
  vector<double>& primal_work_value = info_.workValue_;
  // Take copies of the nonbasic information for the dual LP, since
  // its values will be over-written in constructing the corresponding
  // data for the primal problem
  vector<int8_t> dual_nonbasic_flag = basis_.nonbasicFlag_;
  vector<int8_t> dual_nonbasic_move = basis_.nonbasicMove_;
  vector<HighsInt>& primal_basic_index = basis_.basicIndex_;
  vector<int8_t>& primal_nonbasic_flag = basis_.nonbasicFlag_;
  vector<int8_t>& primal_nonbasic_move = basis_.nonbasicMove_;
  basis_.nonbasicFlag_.assign(primal_num_tot, kIllegalFlagValue);
  basis_.nonbasicMove_.assign(primal_num_tot, kIllegalMoveValue);
  basis_.basicIndex_.resize(0);
  const double inf = kHighsInf;
  // The number of dual rows is the number of primal columns, so all
  // dual basic variables are nonbasic in the primal problem.
  //
  // If there are extra dual columns due to upper bounds on boxed
  // primal variables/constraints, there will be an excess of dual
  // nonbasic variables for the required primal basic variables.
  //
  // For each pair of dual variables associated with a boxed primal
  // variable/constraint:
  //
  // * If one is basic then it yields a nonbasic primal variable, at
  //   a bound given by the basic dual
  //
  // * If both are nonbasic, then they yield a basic primal variable
  //
  // Keep track of the dual column added to handle upper bounds on
  // boxed variables/constraints
  HighsInt upper_bound_col = original_num_row_;
  for (HighsInt iCol = 0; iCol < original_num_col_; iCol++) {
    const double cost = original_col_cost_[iCol];
    const double lower = original_col_lower_[iCol];
    const double upper = original_col_upper_[iCol];
    HighsInt move = kIllegalMoveValue;
    HighsInt dual_variable = dual_num_col + iCol;
    bool dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
    if (lower == upper) {
      // Fixed
      if (dual_basic) move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (dual_basic) {
          // Primal variable is nonbasic at its lower bound
          move = kNonbasicMoveUp;
        } else {
          // Look at the corresponding dual variable for the upper bound
          dual_variable = upper_bound_col;
          dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
          if (dual_basic) {
            // Primal variable is nonbasic at its upper bound
            move = kNonbasicMoveDn;
          }
        }
        upper_bound_col++;
      } else {
        // Lower (since upper bound is infinite)
        if (dual_basic) move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      if (dual_basic) move = kNonbasicMoveDn;
    } else {
      // FREE
      //
      // Dual activity a^Ty is fixed by cost, implying dual for primal
      // column (slack for dual row) is fixed at zero
      if (dual_basic) {
        assert(4 == 0);
        move = kNonbasicMoveZe;
      }
    }
    if (dual_basic) {
      // Primal nonbasic column from basic dual row
      assert(move != kIllegalMoveValue);
      primal_nonbasic_flag[iCol] = kNonbasicFlagTrue;
      primal_nonbasic_move[iCol] = move;
    } else {
      // Primal basic column from nonbasic dual row
      primal_basic_index.push_back(iCol);
      primal_nonbasic_flag[iCol] = kNonbasicFlagFalse;
      primal_nonbasic_move[iCol] = 0;
    }
  }
  for (HighsInt iRow = 0; iRow < original_num_row_; iRow++) {
    double lower = original_row_lower_[iRow];
    double upper = original_row_upper_[iRow];
    HighsInt move = kIllegalMoveValue;
    HighsInt dual_variable = iRow;
    bool dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
    if (lower == upper) {
      // Equality constraint
      //
      // Dual variable has primal RHS as cost and is free
      if (dual_basic) {
        move = kNonbasicMoveZe;
      }
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (dual_basic) {
          // Primal variable is nonbasic at its lower bound
          move = kNonbasicMoveDn;
        } else {
          // Look at the corresponding dual variable for the upper bound
          dual_variable = upper_bound_col;
          dual_basic = dual_nonbasic_flag[dual_variable] == kNonbasicFlagFalse;
          if (dual_basic) {
            // Primal variable is nonbasic at its upper bound
            move = kNonbasicMoveUp;
          }
        }
        upper_bound_col++;
      } else {
        // Lower (since upper bound is infinite)
        if (dual_basic) {
          move = kNonbasicMoveDn;
        }
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      if (dual_basic) {
        move = kNonbasicMoveUp;
      }
    } else {
      // FREE
      if (dual_basic) {
        assert(14 == 0);
        move = kNonbasicMoveZe;
      }
      // Shouldn't get free rows, but handle them anyway
    }
    if (dual_basic) {
      // Primal nonbasic column from basic dual row
      assert(move != kIllegalMoveValue);
      primal_nonbasic_flag[original_num_col_ + iRow] = kNonbasicFlagTrue;
      primal_nonbasic_move[original_num_col_ + iRow] = move;
    } else {
      // Primal basic column from nonbasic dual row
      primal_basic_index.push_back(original_num_col_ + iRow);
      primal_nonbasic_flag[original_num_col_ + iRow] = kNonbasicFlagFalse;
      primal_nonbasic_move[original_num_col_ + iRow] = 0;
    }
  }
  const bool ignore_scaling = true;
  if (!ignore_scaling) {
    // Flip any scale factors
    if (lp_.scale_.has_scaling) {
      std::vector<double> temp_scale = lp_.scale_.row;
      lp_.scale_.row = lp_.scale_.col;
      lp_.scale_.col = temp_scale;
      lp_.scale_.col.resize(original_num_col_);
      lp_.scale_.row.resize(original_num_row_);
      lp_.scale_.num_col = original_num_col_;
      lp_.scale_.num_row = original_num_row_;
    }
  }
  // Change optimzation sense
  if (lp_.sense_ == ObjSense::kMinimize) {
    lp_.sense_ = ObjSense::kMaximize;
  } else {
    lp_.sense_ = ObjSense::kMinimize;
  }
  // Flip LP dimensions
  lp_.num_col_ = original_num_col_;
  lp_.num_row_ = original_num_row_;
  // Restore the original offset
  lp_.offset_ = original_offset_;
  // Copy back the costs and bounds
  lp_.col_cost_ = original_col_cost_;
  lp_.col_lower_ = original_col_lower_;
  lp_.col_upper_ = original_col_upper_;
  lp_.row_lower_ = original_row_lower_;
  lp_.row_upper_ = original_row_upper_;
  // The primal constraint matrix is available row-wise as the first
  // original_num_row_ vectors of the dual constratint matrix
  HighsSparseMatrix primal_matrix;
  primal_matrix.start_.resize(original_num_row_ + 1);
  primal_matrix.index_.resize(original_num_nz_);
  primal_matrix.value_.resize(original_num_nz_);

  for (HighsInt iCol = 0; iCol < original_num_row_ + 1; iCol++)
    primal_matrix.start_[iCol] = lp_.a_matrix_.start_[iCol];
  for (HighsInt iEl = 0; iEl < original_num_nz_; iEl++) {
    primal_matrix.index_[iEl] = lp_.a_matrix_.index_[iEl];
    primal_matrix.value_[iEl] = lp_.a_matrix_.value_[iEl];
  }
  primal_matrix.num_col_ = original_num_col_;
  primal_matrix.num_row_ = original_num_row_;
  primal_matrix.format_ = MatrixFormat::kRowwise;
  // Copy the row-wise primal LP constraint matrix and transpose it.
  // ToDo Make this more efficient
  lp_.a_matrix_ = primal_matrix;
  lp_.a_matrix_.ensureColwise();
  // Some sanity checks
  assert(lp_.num_col_ == original_num_col_);
  assert(lp_.num_row_ == original_num_row_);
  assert(lp_.a_matrix_.numNz() == original_num_nz_);
  HighsInt num_basic_variables = primal_basic_index.size();
  bool num_basic_variables_ok = num_basic_variables == original_num_row_;
  if (!num_basic_variables_ok)
    printf("HEkk::undualise: Have %d basic variables, not %d\n",
           (int)num_basic_variables, (int)original_num_row_);
  assert(num_basic_variables_ok);

  // Clear the data retained when solving dual LP
  clearEkkDualise();
  status_.is_dualised = false;
  // Now solve with this basis. Should just be a case of reinverting
  // and re-solving for optimal primal and dual values, but
  // numerically marginal LPs will need clean-up
  status_.has_basis = true;
  status_.has_ar_matrix = false;
  status_.has_nla = false;
  status_.has_invert = false;
  HighsInt primal_solve_iteration_count = -iteration_count_;
  printf(
      "HEkk::undualise() Solving the primal LP using the optimal basis of its "
      "dual\n");
  HighsStatus return_status = solve();
  primal_solve_iteration_count += iteration_count_;
  //  if (primal_solve_iteration_count)
  highsLogUser(options_->log_options, HighsLogType::kInfo,
               "Solving the primal LP (%s) using the optimal basis of its dual "
               "required %d simplex iterations\n",
               lp_.model_name_.c_str(), (int)primal_solve_iteration_count);
  return return_status;
}

HighsStatus HEkk::permute() {
  assert(1 == 0);
  return HighsStatus::kError;
}

HighsStatus HEkk::unpermute() {
  if (!this->status_.is_permuted) return HighsStatus::kOk;
  assert(1 == 0);
  return HighsStatus::kError;
}

HighsStatus HEkk::solve() {
  initialiseAnalysis();
  initialiseControl();

  if (analysis_.analyse_simplex_time)
    analysis_.simplexTimerStart(SimplexTotalClock);
  dual_simplex_cleanup_level_ = 0;
  dual_simplex_phase1_cleanup_level_ = 0;
  HighsStatus call_status = initialiseForSolve();
  if (call_status == HighsStatus::kError) return HighsStatus::kError;

  const HighsDebugStatus simplex_nla_status =
      simplex_nla_.debugCheckData("Before HEkk::solve()");
  const bool simplex_nla_ok = simplex_nla_status == HighsDebugStatus::kOk;
  if (!simplex_nla_ok) {
    highsLogUser(options_->log_options, HighsLogType::kError,
                 "Error in simplex NLA data\n");
    assert(simplex_nla_ok);
    return HighsStatus::kError;
  }

  assert(status_.has_basis);
  assert(status_.has_invert);
  assert(status_.initialised_for_solve);
  if (model_status_ == HighsModelStatus::kOptimal) return HighsStatus::kOk;

  HighsStatus return_status = HighsStatus::kOk;
  std::string algorithm_name;

  // Indicate that dual and primal rays are not known
  status_.has_dual_ray = false;
  status_.has_primal_ray = false;

  // Allow primal and dual perturbations in case a block on them is
  // hanging over from a previous call
  info_.allow_cost_perturbation = true;
  info_.allow_bound_perturbation = true;

  chooseSimplexStrategyThreads(*options_, info_);
  HighsInt& simplex_strategy = info_.simplex_strategy;

  // Initial solve according to strategy
  if (simplex_strategy == kSimplexStrategyPrimal) {
    algorithm_name = "primal";
    reportSimplexPhaseIterations(options_->log_options, iteration_count_, info_,
                                 true);
    highsLogUser(options_->log_options, HighsLogType::kInfo,
                 "Using EKK primal simplex solver\n");
    HEkkPrimal primal_solver(*this);
    workEdWt_ = NULL;
    workEdWtFull_ = NULL;
    call_status = primal_solver.solve();
    assert(called_return_from_solve_);
    return_status = interpretCallStatus(options_->log_options, call_status,
                                        return_status, "HEkkPrimal::solve");
  } else {
    algorithm_name = "dual";
    reportSimplexPhaseIterations(options_->log_options, iteration_count_, info_,
                                 true);
    // Solve, depending on the particular strategy
    if (simplex_strategy == kSimplexStrategyDualTasks) {
      highsLogUser(
          options_->log_options, HighsLogType::kInfo,
          "Using EKK parallel dual simplex solver - SIP with %" HIGHSINT_FORMAT
          " threads\n",
          info_.num_threads);
    } else if (simplex_strategy == kSimplexStrategyDualMulti) {
      highsLogUser(
          options_->log_options, HighsLogType::kInfo,
          "Using EKK parallel dual simplex solver - PAMI with %" HIGHSINT_FORMAT
          " threads\n",
          info_.num_threads);
    } else {
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   "Using EKK dual simplex solver - serial\n");
    }
    HEkkDual dual_solver(*this);
    workEdWt_ = dual_solver.getWorkEdWt();
    workEdWtFull_ = dual_solver.getWorkEdWtFull();
    call_status = dual_solver.solve();
    assert(called_return_from_solve_);
    return_status = interpretCallStatus(options_->log_options, call_status,
                                        return_status, "HEkkDual::solve");

    // Dual simplex solver may set model_status to be
    // kUnboundedOrInfeasible, and Highs::run() may not allow that to
    // be returned, so use primal simplex to distinguish
    if (model_status_ == HighsModelStatus::kUnboundedOrInfeasible &&
        !options_->allow_unbounded_or_infeasible) {
      HEkkPrimal primal_solver(*this);
      call_status = primal_solver.solve();
      assert(called_return_from_solve_);
      return_status = interpretCallStatus(options_->log_options, call_status,
                                          return_status, "HEkkPrimal::solve");
    }
  }
  reportSimplexPhaseIterations(options_->log_options, iteration_count_, info_);
  if (return_status == HighsStatus::kError) return return_status;
  highsLogDev(options_->log_options, HighsLogType::kInfo,
              "EKK %s simplex solver returns %" HIGHSINT_FORMAT
              " primal and %" HIGHSINT_FORMAT
              " dual infeasibilities: "
              "Status %s\n",
              algorithm_name.c_str(), info_.num_primal_infeasibilities,
              info_.num_dual_infeasibilities,
              utilModelStatusToString(model_status_).c_str());
  // Can model_status_ = HighsModelStatus::kNotset be returned?
  assert(model_status_ != HighsModelStatus::kNotset);

  if (analysis_.analyse_simplex_time) {
    analysis_.simplexTimerStop(SimplexTotalClock);
    analysis_.reportSimplexTimer();
  }
  if (analysis_.analyse_simplex_summary_data) analysis_.summaryReport();
  if (analysis_.analyse_factor_data) analysis_.reportInvertFormData();
  if (analysis_.analyse_factor_time) analysis_.reportFactorTimer();
  return return_status;
}

HighsStatus HEkk::cleanup() {
  // Clean up from a point with either primal or dual
  // infeasiblilities, but not both
  HighsStatus return_status = HighsStatus::kOk;
  HighsStatus call_status;
  if (info_.num_primal_infeasibilities) {
    // Primal infeasibilities, so should be just dual phase 2
    assert(!info_.num_dual_infeasibilities);
    // Use dual simplex (phase 2) with Devex pricing and no perturbation
    info_.simplex_strategy = kSimplexStrategyDual;
    info_.dual_simplex_cost_perturbation_multiplier = 0;
    info_.dual_edge_weight_strategy = kSimplexDualEdgeWeightStrategyDevex;
    HEkkDual dual_solver(*this);
    workEdWt_ = dual_solver.getWorkEdWt();
    workEdWtFull_ = dual_solver.getWorkEdWtFull();
    call_status = dual_solver.solve();
    assert(called_return_from_solve_);
    return_status =
        interpretCallStatus(this->options_->log_options, call_status,
                            return_status, "HEkkDual::solve");
    if (return_status == HighsStatus::kError) return return_status;
  } else {
    // Dual infeasibilities, so should be just primal phase 2
    assert(!info_.num_primal_infeasibilities);
    // Use primal simplex (phase 2) with no perturbation
    info_.simplex_strategy = kSimplexStrategyPrimal;
    info_.primal_simplex_bound_perturbation_multiplier = 0;
    HEkkPrimal primal_solver(*this);
    workEdWt_ = NULL;
    workEdWtFull_ = NULL;
    call_status = primal_solver.solve();
    assert(called_return_from_solve_);
    return_status =
        interpretCallStatus(this->options_->log_options, call_status,
                            return_status, "HEkkPrimal::solve");
    if (return_status == HighsStatus::kError) return return_status;
  }
  return return_status;
}

HighsStatus HEkk::setBasis() {
  // Set up nonbasicFlag and basicIndex for a logical basis
  const HighsInt num_col = lp_.num_col_;
  const HighsInt num_row = lp_.num_row_;
  const HighsInt num_tot = num_col + num_row;
  basis_.nonbasicFlag_.resize(num_tot);
  basis_.nonbasicMove_.resize(num_tot);
  basis_.basicIndex_.resize(num_row);
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    basis_.nonbasicFlag_[iCol] = kNonbasicFlagTrue;
    double lower = lp_.col_lower_[iCol];
    double upper = lp_.col_upper_[iCol];
    HighsInt move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed. Set to bound of LP that is closer to
        // zero
        if (move == kIllegalMoveValue) {
          if (fabs(lower) < fabs(upper)) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iCol] = move;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
    basis_.basicIndex_[iRow] = iVar;
  }
  info_.num_basic_logicals = num_row;
  status_.has_basis = true;
  return HighsStatus::kOk;
}

HighsStatus HEkk::setBasis(const HighsBasis& highs_basis) {
  // Shouldn't have to check the incoming basis since this is an
  // internal call, but it may be a basis that's set up internally
  // with errors :-) ...
  //
  HighsOptions& options = *options_;
  if (debugHighsBasisConsistent(options, lp_, highs_basis) ==
      HighsDebugStatus::kLogicalError) {
    highsLogDev(options_->log_options, HighsLogType::kError,
                "Supposed to be a Highs basis, but not valid\n");
    return HighsStatus::kError;
  }
  HighsInt num_col = lp_.num_col_;
  HighsInt num_row = lp_.num_row_;
  HighsInt num_tot = num_col + num_row;
  // Resize the basis in case none has yet been defined for this LP
  basis_.nonbasicFlag_.resize(num_tot);
  basis_.nonbasicMove_.resize(num_tot);
  basis_.basicIndex_.resize(num_row);

  HighsInt num_basic_variables = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    HighsInt iVar = iCol;

    const double lower = lp_.col_lower_[iCol];
    const double upper = lp_.col_upper_[iCol];
    if (highs_basis.col_status[iCol] == HighsBasisStatus::kBasic) {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
      basis_.nonbasicMove_[iVar] = 0;
      basis_.basicIndex_[num_basic_variables++] = iVar;
    } else {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagTrue;
      if (highs_basis.col_status[iCol] == HighsBasisStatus::kLower) {
        if (lower == upper) {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
        } else {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveUp;
        }
      } else if (highs_basis.col_status[iCol] == HighsBasisStatus::kUpper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveDn;
      } else {
        assert(highs_basis.col_status[iCol] == HighsBasisStatus::kZero);
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      }
    }
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    const double lower = lp_.row_lower_[iRow];
    const double upper = lp_.row_upper_[iRow];
    if (highs_basis.row_status[iRow] == HighsBasisStatus::kBasic) {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagFalse;
      basis_.nonbasicMove_[iVar] = 0;
      basis_.basicIndex_[num_basic_variables++] = iVar;
    } else {
      basis_.nonbasicFlag_[iVar] = kNonbasicFlagTrue;
      if (highs_basis.row_status[iRow] == HighsBasisStatus::kLower) {
        if (lower == upper) {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
        } else {
          basis_.nonbasicMove_[iVar] = kNonbasicMoveDn;
        }
      } else if (highs_basis.row_status[iRow] == HighsBasisStatus::kUpper) {
        basis_.nonbasicMove_[iVar] = kNonbasicMoveUp;
      } else {
        assert(highs_basis.row_status[iRow] == HighsBasisStatus::kZero);
        basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      }
    }
  }
  status_.has_basis = true;
  return HighsStatus::kOk;
}

/*
HighsStatus HEkk::setBasis(const SimplexBasis& basis) {
  // Shouldn't have to check the incoming basis since this is an
  // internal call, but it may be a basis that's set up internally
  // with errors :-) ...
  if (this->debugBasisConsistent() ==
      HighsDebugStatus::kLogicalError) {
    highsLogDev(this->options_->log_options, HighsLogType::kError,
                "Supposed to be a Highs basis, but not valid\n");
    return HighsStatus::kError;
  }
  this->basis_.nonbasicFlag_ = basis.nonbasicFlag_;
  this->basis_.nonbasicMove_ = basis.nonbasicMove_;
  this->basis_.basicIndex_ = basis.basicIndex_;
  this->status_.has_basis = true;
  return HighsStatus::kOk;
}
*/

void HEkk::addCols(const HighsLp& lp,
                   const HighsSparseMatrix& scaled_a_matrix) {
  // Should be extendSimplexLpRandomVectors
  //  if (valid_simplex_basis)
  //    appendBasicRowsToBasis(simplex_lp, simplex_basis, XnumNewRow);
  //  ekk_instance_.updateStatus(LpAction::kNewRows);
  //  if (valid_simplex_lp) {
  //    simplex_lp.num_row_ += XnumNewRow;
  //    ekk_instance_.initialiseSimplexLpRandomVectors();
  //  }
  //  if (valid_simplex_lp)
  //    assert(ekk_instance_.lp_.dimensionsOk("addCols - simplex"));
  if (this->status_.has_nla) this->simplex_nla_.addCols(&lp);
  this->updateStatus(LpAction::kNewCols);
}

void HEkk::addRows(const HighsLp& lp,
                   const HighsSparseMatrix& scaled_ar_matrix) {
  // Should be extendSimplexLpRandomVectors
  //  if (valid_simplex_basis)
  //    appendBasicRowsToBasis(simplex_lp, simplex_basis, XnumNewRow);
  //  ekk_instance_.updateStatus(LpAction::kNewRows);
  //  if (valid_simplex_lp) {
  //    simplex_lp.num_row_ += XnumNewRow;
  //    ekk_instance_.initialiseSimplexLpRandomVectors();
  //  }
  //  if (valid_simplex_lp)
  //    assert(ekk_instance_.lp_.dimensionsOk("addRows - simplex"));
  if (kExtendInvertWhenAddingRows && this->status_.has_nla) {
    this->simplex_nla_.addRows(&lp, &basis_.basicIndex_[0], &scaled_ar_matrix);
    setNlaPointersForTrans(lp);
    this->debugNlaCheckInvert("HEkk::addRows - on entry",
                              kHighsDebugLevelExpensive + 1);
  }
  // Update the number of rows in the simplex LP so that it's
  // consistent with simplex basis information
  this->lp_.num_row_ = lp.num_row_;
  this->updateStatus(LpAction::kNewRows);
}

void HEkk::deleteCols(const HighsIndexCollection& index_collection) {
  this->updateStatus(LpAction::kDelCols);
}
void HEkk::deleteRows(const HighsIndexCollection& index_collection) {
  this->updateStatus(LpAction::kDelRows);
}

void HEkk::unscaleSimplex(const HighsLp& incumbent_lp) {
  if (!this->simplex_in_scaled_space_) return;
  assert(incumbent_lp.scale_.has_scaling);
  const HighsInt num_col = incumbent_lp.num_col_;
  const HighsInt num_row = incumbent_lp.num_row_;
  const vector<double>& col_scale = incumbent_lp.scale_.col;
  const vector<double>& row_scale = incumbent_lp.scale_.row;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    const HighsInt iVar = iCol;
    const double factor = col_scale[iCol];
    this->info_.workCost_[iVar] /= factor;
    this->info_.workDual_[iVar] /= factor;
    this->info_.workShift_[iVar] /= factor;
    this->info_.workLower_[iVar] *= factor;
    this->info_.workUpper_[iVar] *= factor;
    this->info_.workRange_[iVar] *= factor;
    this->info_.workValue_[iVar] *= factor;
    this->info_.workLowerShift_[iVar] *= factor;
    this->info_.workUpperShift_[iVar] *= factor;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    const HighsInt iVar = num_col + iRow;
    const double factor = row_scale[iRow];
    this->info_.workCost_[iVar] *= factor;
    this->info_.workDual_[iVar] *= factor;
    this->info_.workShift_[iVar] *= factor;
    this->info_.workLower_[iVar] /= factor;
    this->info_.workUpper_[iVar] /= factor;
    this->info_.workRange_[iVar] /= factor;
    this->info_.workValue_[iVar] /= factor;
    this->info_.workLowerShift_[iVar] /= factor;
    this->info_.workUpperShift_[iVar] /= factor;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    double factor;
    const HighsInt iVar = this->basis_.basicIndex_[iRow];
    if (iVar < num_col) {
      factor = col_scale[iVar];
    } else {
      factor = 1.0 / row_scale[iVar - num_col];
    }
    this->info_.baseLower_[iRow] *= factor;
    this->info_.baseUpper_[iRow] *= factor;
    this->info_.baseValue_[iRow] *= factor;
  }
  this->simplex_in_scaled_space_ = false;
}
HighsSolution HEkk::getSolution() {
  HighsSolution solution;
  // Scatter the basic primal values
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++)
    info_.workValue_[basis_.basicIndex_[iRow]] = info_.baseValue_[iRow];
  // Zero the basic dual values
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++)
    info_.workDual_[basis_.basicIndex_[iRow]] = 0;

  // Now we can get the solution
  solution.col_value.resize(lp_.num_col_);
  solution.col_dual.resize(lp_.num_col_);
  solution.row_value.resize(lp_.num_row_);
  solution.row_dual.resize(lp_.num_row_);

  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    solution.col_value[iCol] = info_.workValue_[iCol];
    solution.col_dual[iCol] = (HighsInt)lp_.sense_ * info_.workDual_[iCol];
  }
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    solution.row_value[iRow] = -info_.workValue_[lp_.num_col_ + iRow];
    // @FlipRowDual negate RHS
    solution.row_dual[iRow] =
        -(HighsInt)lp_.sense_ * info_.workDual_[lp_.num_col_ + iRow];
  }
  solution.value_valid = true;
  solution.dual_valid = true;
  return solution;
}

HighsBasis HEkk::getHighsBasis(HighsLp& use_lp) const {
  HighsInt num_col = use_lp.num_col_;
  HighsInt num_row = use_lp.num_row_;
  HighsBasis highs_basis;
  highs_basis.col_status.resize(num_col);
  highs_basis.row_status.resize(num_row);
  assert(status_.has_basis);
  highs_basis.valid = false;
  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    HighsInt iVar = iCol;
    const double lower = use_lp.col_lower_[iCol];
    const double upper = use_lp.col_upper_[iCol];
    HighsBasisStatus basis_status = HighsBasisStatus::kNonbasic;
    if (!basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::kBasic;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveUp) {
      basis_status = HighsBasisStatus::kLower;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveDn) {
      basis_status = HighsBasisStatus::kUpper;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveZe) {
      if (lower == upper) {
        basis_status = HighsBasisStatus::kLower;
      } else {
        basis_status = HighsBasisStatus::kZero;
      }
    }
    highs_basis.col_status[iCol] = basis_status;
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = num_col + iRow;
    const double lower = use_lp.row_lower_[iRow];
    const double upper = use_lp.row_upper_[iRow];
    HighsBasisStatus basis_status = HighsBasisStatus::kNonbasic;
    if (!basis_.nonbasicFlag_[iVar]) {
      basis_status = HighsBasisStatus::kBasic;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveUp) {
      basis_status = HighsBasisStatus::kUpper;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveDn) {
      basis_status = HighsBasisStatus::kLower;
    } else if (basis_.nonbasicMove_[iVar] == kNonbasicMoveZe) {
      if (lower == upper) {
        basis_status = HighsBasisStatus::kLower;
      } else {
        basis_status = HighsBasisStatus::kZero;
      }
    }
    highs_basis.row_status[iRow] = basis_status;
  }
  highs_basis.valid = true;
  return highs_basis;
}

HighsInt HEkk::initialiseSimplexLpBasisAndFactor(
    const bool only_from_known_basis) {
  // If there's no basis, return error if the basis has to be known,
  // otherwise set a logical basis
  //
  // If the basis has to be known, then non-negative return value is
  // rank deficiency; negative return value is error. Otherwise, any
  // rank deficiency is handled and return is 0
  if (!status_.has_basis) {
    if (only_from_known_basis) {
      highsLogDev(options_->log_options, HighsLogType::kError,
                  "Simplex basis should be known but isn't\n");
      return -1;
    }
    setBasis();
  }

  //
  // The simplex NLA operates in the scaled space if the LP has
  // scaling factors. If they exist but haven't been applied, then the
  // simplex NLA needs a separate, scaled constraint matrix. Thus
  // getScaledAMatrixPointer() returns a pointer to either the
  // constraint matrix or a scaled copy (that is a member of the HEkk
  // class), with the latter returned if the LP has scaling factors
  // but is unscaled.
  //
  HighsSparseMatrix* local_scaled_a_matrix = getScaledAMatrixPointer();
  //
  // If simplex NLA is set up, pass the pointers that it uses. It
  // deduces any scaling factors that it must use by inspecting
  // whether the LP has scaling factors, and whether it is scaled.
  //
  // If simplex NLA is not set up, then it will be done if
  //
  if (this->status_.has_nla) {
    this->simplex_nla_.setPointers(&(this->lp_), local_scaled_a_matrix,
                                   &this->basis_.basicIndex_[0], this->options_,
                                   this->timer_, &(this->analysis_));
  } else {
    // todo @ Julian: this fails on glass4
    assert(info_.factor_pivot_threshold >= options_->factor_pivot_threshold);
    simplex_nla_.setup(&(this->lp_),                  //&lp_,
                       &this->basis_.basicIndex_[0],  //&basis_.basicIndex_[0],
                       this->options_,                // options_,
                       this->timer_,                  // timer_,
                       &(this->analysis_),            //&analysis_,
                       local_scaled_a_matrix,
                       this->info_.factor_pivot_threshold);
    status_.has_nla = true;
  }

  if (!status_.has_invert) {
    const HighsInt rank_deficiency = computeFactor();
    if (rank_deficiency) {
      // Basis is rank deficient
      if (only_from_known_basis) {
        // If only this basis should be used, then return error
        highsLogDev(options_->log_options, HighsLogType::kError,
                    "Supposed to be a full-rank basis, but incorrect\n");
        return rank_deficiency;
      }
      // Account for rank deficiency by correcing nonbasicFlag
      handleRankDeficiency();
      this->updateStatus(LpAction::kNewBasis);
      setNonbasicMove();
      status_.has_basis = true;
      status_.has_invert = true;
      status_.has_fresh_invert = true;
    }
    // Record the synthetic clock for INVERT, and zero it for UPDATE
    resetSyntheticClock();
  }
  assert(status_.has_invert);
  return 0;
}

void HEkk::handleRankDeficiency() {
  HFactor& factor = simplex_nla_.factor_;
  HighsInt rank_deficiency = factor.rank_deficiency;
  vector<HighsInt>& noPvC = factor.noPvC;
  vector<HighsInt>& noPvR = factor.noPvR;
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    HighsInt variable_in = lp_.num_col_ + noPvR[k];
    HighsInt variable_out = noPvC[k];
    basis_.nonbasicFlag_[variable_in] = kNonbasicFlagFalse;
    basis_.nonbasicFlag_[variable_out] = kNonbasicFlagTrue;
  }
  status_.has_ar_matrix = false;
}

// Private methods

void HEkk::initialiseEkk() {
  if (status_.initialised_for_new_lp) return;
  setSimplexOptions();
  initialiseControl();
  initialiseSimplexLpRandomVectors();
  simplex_nla_.clear();
  status_.initialised_for_new_lp = true;
}

bool HEkk::isUnconstrainedLp() {
  bool is_unconstrained_lp = lp_.num_row_ <= 0;
  if (is_unconstrained_lp)
    highsLogDev(
        options_->log_options, HighsLogType::kError,
        "HEkkDual::solve called for LP with non-positive (%" HIGHSINT_FORMAT
        ") number of constraints\n",
        lp_.num_row_);
  assert(!is_unconstrained_lp);
  return is_unconstrained_lp;
}

HighsStatus HEkk::initialiseForSolve() {
  const HighsInt error_return = initialiseSimplexLpBasisAndFactor();
  assert(!error_return);
  if (error_return) return HighsStatus::kError;
  assert(status_.has_basis);

  updateSimplexOptions();
  initialiseSimplexLpRandomVectors();
  initialisePartitionedRowwiseMatrix();  // Timed
  allocateWorkAndBaseArrays();
  initialiseCost(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
  initialiseBound(SimplexAlgorithm::kPrimal, kSolvePhaseUnknown, false);
  initialiseNonbasicValueAndMove();
  computePrimal();                // Timed
  computeDual();                  // Timed
  computeSimplexInfeasible();     // Timed
  computeDualObjectiveValue();    // Timed
  computePrimalObjectiveValue();  // Timed
  status_.initialised_for_solve = true;

  bool primal_feasible = info_.num_primal_infeasibilities == 0;
  bool dual_feasible = info_.num_dual_infeasibilities == 0;
  model_status_ = HighsModelStatus::kNotset;
  if (primal_feasible && dual_feasible)
    model_status_ = HighsModelStatus::kOptimal;
  return HighsStatus::kOk;
}

void HEkk::setSimplexOptions() {
  // Copy values of HighsOptions for the simplex solver
  // Currently most of these options are straight copies, but they
  // will become valuable when "choose" becomes a HiGHS strategy value
  // that will need converting into a specific simplex strategy value.
  //
  // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
  //
  info_.dual_edge_weight_strategy = options_->simplex_dual_edge_weight_strategy;
  info_.price_strategy = options_->simplex_price_strategy;
  info_.dual_simplex_cost_perturbation_multiplier =
      options_->dual_simplex_cost_perturbation_multiplier;
  info_.primal_simplex_bound_perturbation_multiplier =
      options_->primal_simplex_bound_perturbation_multiplier;
  info_.factor_pivot_threshold = options_->factor_pivot_threshold;
  info_.update_limit = options_->simplex_update_limit;
  random_.initialise(options_->random_seed);

  // Set values of internal options
  info_.store_squared_primal_infeasibility = true;
}

void HEkk::updateSimplexOptions() {
  // Update some simplex option values from HighsOptions when
  // (re-)solving an LP. Others aren't changed because better values
  // may have been learned due to solving this LP (possibly with some
  // modification) before.
  //
  // NB simplex_strategy is set by chooseSimplexStrategyThreads in each call
  //
  info_.dual_simplex_cost_perturbation_multiplier =
      options_->dual_simplex_cost_perturbation_multiplier;
  info_.primal_simplex_bound_perturbation_multiplier =
      options_->primal_simplex_bound_perturbation_multiplier;
}

void HEkk::initialiseSimplexLpRandomVectors() {
  const HighsInt num_col = lp_.num_col_;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  if (!num_tot) return;
  // Instantiate and (re-)initialise the random number generator
  //  HighsRandom random;
  HighsRandom& random = random_;
  //  random.initialise();

  if (num_col) {
    // Generate a random permutation of the column indices
    vector<HighsInt>& numColPermutation = info_.numColPermutation_;
    numColPermutation.resize(num_col);
    for (HighsInt i = 0; i < num_col; i++) numColPermutation[i] = i;
    random.shuffle(numColPermutation.data(), num_col);
  }

  // Re-initialise the random number generator and generate the
  // random vectors in the same order as hsol to maintain repeatable
  // performance
  // random.initialise();

  // Generate a random permutation of all the indices
  vector<HighsInt>& numTotPermutation = info_.numTotPermutation_;
  numTotPermutation.resize(num_tot);
  for (HighsInt i = 0; i < num_tot; i++) numTotPermutation[i] = i;
  random.shuffle(numTotPermutation.data(), num_tot);

  // Generate a vector of random reals
  info_.numTotRandomValue_.resize(num_tot);
  vector<double>& numTotRandomValue = info_.numTotRandomValue_;
  for (HighsInt i = 0; i < num_tot; i++) {
    numTotRandomValue[i] = random.fraction();
  }
}

void HEkk::chooseSimplexStrategyThreads(const HighsOptions& options,
                                        HighsSimplexInfo& info) {
  // Ensure that this is not called with an optimal basis
  assert(info.num_dual_infeasibilities > 0 ||
         info.num_primal_infeasibilities > 0);
  // Set the internal simplex strategy and number of threads for dual
  // simplex
  HighsInt& simplex_strategy = info.simplex_strategy;
  // By default, use the HighsOptions strategy. If this is
  // kSimplexStrategyChoose, then the strategy used will depend on
  // whether the current basis is primal feasible.
  simplex_strategy = options.simplex_strategy;
  if (simplex_strategy == kSimplexStrategyChoose) {
    // HiGHS is left to choose the simplex strategy
    if (info.num_primal_infeasibilities > 0) {
      // Not primal feasible, so use dual simplex
      simplex_strategy = kSimplexStrategyDual;
    } else {
      // Primal feasible. so use primal simplex
      simplex_strategy = kSimplexStrategyPrimal;
    }
  }
  // Set min/max_threads to correspond to serial code. They will be
  // set to other values if parallel options are used.
  info.min_threads = 1;
  info.max_threads = 1;
  // Record the min/max minimum number of HiGHS threads in the options
  const HighsInt highs_min_threads = options.highs_min_threads;
  const HighsInt highs_max_threads = options.highs_max_threads;
  HighsInt omp_max_threads = 0;
#ifdef OPENMP
  omp_max_threads = omp_get_max_threads();
#endif
  if (options.parallel == kHighsOnString &&
      simplex_strategy == kSimplexStrategyDual) {
    // The parallel strategy is on and the simplex strategy is dual so use
    // PAMI if there are enough OMP threads
    if (omp_max_threads >= kDualMultiMinThreads)
      simplex_strategy = kSimplexStrategyDualMulti;
  }
  //
  // If parallel stratgies are used, the minimum number of HiGHS threads used
  // will be set to be at least the minimum required for the strategy
  //
  // All this is independent of the number of OMP threads available,
  // since code with multiple HiGHS threads can be run in serial.
#ifdef OPENMP
  if (simplex_strategy == kSimplexStrategyDualTasks) {
    info.min_threads = max(kDualTasksMinThreads, highs_min_threads);
    info.max_threads = max(info.min_threads, highs_max_threads);
  } else if (simplex_strategy == kSimplexStrategyDualMulti) {
    info.min_threads = max(kDualMultiMinThreads, highs_min_threads);
    info.max_threads = max(info.min_threads, highs_max_threads);
  }
#endif
  // Set the number of HiGHS threads to be used to be the maximum
  // number to be used
  info.num_threads = info.max_threads;
  // Give a warning if the number of threads to be used is fewer than
  // the minimum number of HiGHS threads allowed
  if (info.num_threads < highs_min_threads) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Using %" HIGHSINT_FORMAT
                 " HiGHS threads for parallel strategy rather than "
                 "minimum number (%" HIGHSINT_FORMAT ") specified in options\n",
                 info.num_threads, highs_min_threads);
  }
  // Give a warning if the number of threads to be used is more than
  // the maximum number of HiGHS threads allowed
  if (info.num_threads > highs_max_threads) {
    highsLogUser(options.log_options, HighsLogType::kWarning,
                 "Using %" HIGHSINT_FORMAT
                 " HiGHS threads for parallel strategy rather than "
                 "maximum number (%" HIGHSINT_FORMAT ") specified in options\n",
                 info.num_threads, highs_max_threads);
  }
  // Give a warning if the number of threads to be used is fewer than
  // the number of OMP threads available
  if (info.num_threads > omp_max_threads) {
    highsLogUser(
        options.log_options, HighsLogType::kWarning,
        "Number of OMP threads available = %" HIGHSINT_FORMAT
        " < %" HIGHSINT_FORMAT
        " = Number of HiGHS threads "
        "to be used: Parallel performance will be less than anticipated\n",
        omp_max_threads, info.num_threads);
  }
}

bool HEkk::getNonsingularInverse(const HighsInt solve_phase) {
  assert(status_.has_basis);
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;
  // Take a copy of basicIndex from before INVERT to be used as the
  // saved ordering of basic variables - so reinvert will run
  // identically.
  const vector<HighsInt> basicIndex_before_compute_factor = basicIndex;
  // Save the number of updates performed in case it has to be used to determine
  // a limit
  const HighsInt simplex_update_count = info_.update_count;
  // Dual simplex edge weights are identified with rows, so must be
  // permuted according to INVERT. This must be done if workEdWt_ is
  // not NULL.
  const bool handle_edge_weights = workEdWt_ != NULL;
  // Scatter the edge weights so that, after INVERT, they can be
  // gathered according to the new permutation of basicIndex
  if (handle_edge_weights) {
    analysis_.simplexTimerStart(PermWtClock);
    for (HighsInt i = 0; i < lp_.num_row_; i++)
      workEdWtFull_[basicIndex[i]] = workEdWt_[i];
    analysis_.simplexTimerStop(PermWtClock);
  }

  // Call computeFactor to perform INVERT
  HighsInt rank_deficiency = computeFactor();
  //  if (rank_deficiency) {printf("Recover from rank_deficiency of %d\n",
  //  (int)rank_deficiency);}
  const bool artificial_rank_deficiency = false;  //  true;//
  if (artificial_rank_deficiency) {
    if (!info_.phase1_backtracking_test_done && solve_phase == kSolvePhase1) {
      // Claim rank deficiency to test backtracking
      printf("Phase1 (Iter %" HIGHSINT_FORMAT
             ") Claiming rank deficiency to test backtracking\n",
             iteration_count_);
      rank_deficiency = 1;
      info_.phase1_backtracking_test_done = true;
    } else if (!info_.phase2_backtracking_test_done &&
               solve_phase == kSolvePhase2) {
      // Claim rank deficiency to test backtracking
      printf("Phase2 (Iter %" HIGHSINT_FORMAT
             ") Claiming rank deficiency to test backtracking\n",
             iteration_count_);
      rank_deficiency = 1;
      info_.phase2_backtracking_test_done = true;
    }
  }
  if (rank_deficiency) {
    // Rank deficient basis, so backtrack to last full rank basis
    //
    // Get the last nonsingular basis - so long as there is one
    if (!getBacktrackingBasis(workEdWtFull_)) return false;
    // Record that backtracking is taking place
    info_.backtracking_ = true;
    this->updateStatus(LpAction::kBacktracking);
    HighsInt backtrack_rank_deficiency = computeFactor();
    // This basis has previously been inverted successfully, so it shouldn't be
    // singular
    if (backtrack_rank_deficiency) return false;
    // simplex update limit will be half of the number of updates
    // performed, so make sure that at least one update was performed
    if (simplex_update_count <= 1) return false;
    HighsInt use_simplex_update_limit = info_.update_limit;
    HighsInt new_simplex_update_limit = simplex_update_count / 2;
    info_.update_limit = new_simplex_update_limit;
    highsLogDev(options_->log_options, HighsLogType::kWarning,
                "Rank deficiency of %" HIGHSINT_FORMAT
                " after %" HIGHSINT_FORMAT
                " simplex updates, so "
                "backtracking: max updates reduced from %" HIGHSINT_FORMAT
                " to %" HIGHSINT_FORMAT "\n",
                rank_deficiency, simplex_update_count, use_simplex_update_limit,
                new_simplex_update_limit);
  } else {
    // Current basis is full rank so save it
    putBacktrackingBasis(basicIndex_before_compute_factor, workEdWtFull_);
    // Indicate that backtracking is not taking place
    info_.backtracking_ = false;
    // Reset the update limit in case this is the first successful
    // inversion after backtracking
    info_.update_limit = options_->simplex_update_limit;
  }
  if (handle_edge_weights) {
    // Gather the edge weights according to the permutation of
    // basicIndex after INVERT
    analysis_.simplexTimerStart(PermWtClock);
    for (HighsInt i = 0; i < lp_.num_row_; i++)
      workEdWt_[i] = workEdWtFull_[basicIndex[i]];
    analysis_.simplexTimerStop(PermWtClock);
  }
  return true;
}

bool HEkk::getBacktrackingBasis(double* scattered_edge_weights) {
  if (!info_.valid_backtracking_basis_) return false;
  basis_ = info_.backtracking_basis_;
  info_.costs_perturbed = info_.backtracking_basis_costs_perturbed_;
  info_.workShift_ = info_.backtracking_basis_workShift_;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  const bool handle_edge_weights = scattered_edge_weights != NULL;
  if (handle_edge_weights) {
    for (HighsInt iVar = 0; iVar < num_tot; iVar++)
      scattered_edge_weights[iVar] =
          info_.backtracking_basis_edge_weights_[iVar];
  }
  return true;
}

void HEkk::putBacktrackingBasis() {
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;
  const bool handle_edge_weights = workEdWt_ != NULL;
  if (handle_edge_weights) {
    analysis_.simplexTimerStart(PermWtClock);
    for (HighsInt i = 0; i < lp_.num_row_; i++)
      workEdWtFull_[basicIndex[i]] = workEdWt_[i];
    analysis_.simplexTimerStop(PermWtClock);
  }
  putBacktrackingBasis(basicIndex, workEdWtFull_);
}

void HEkk::putBacktrackingBasis(
    const vector<HighsInt>& basicIndex_before_compute_factor,
    double* scattered_edge_weights) {
  info_.valid_backtracking_basis_ = true;
  info_.backtracking_basis_ = basis_;
  info_.backtracking_basis_.basicIndex_ = basicIndex_before_compute_factor;
  info_.backtracking_basis_costs_perturbed_ = info_.costs_perturbed;
  info_.backtracking_basis_workShift_ = info_.workShift_;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  const bool handle_edge_weights = scattered_edge_weights != NULL;
  if (handle_edge_weights) {
    for (HighsInt iVar = 0; iVar < num_tot; iVar++)
      info_.backtracking_basis_edge_weights_[iVar] =
          scattered_edge_weights[iVar];
  }
}

void HEkk::computePrimalObjectiveValue() {
  analysis_.simplexTimerStart(ComputePrObjClock);
  info_.primal_objective_value = 0;
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iVar = basis_.basicIndex_[iRow];
    if (iVar < lp_.num_col_) {
      info_.primal_objective_value +=
          info_.baseValue_[iRow] * lp_.col_cost_[iVar];
    }
  }
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    if (basis_.nonbasicFlag_[iCol])
      info_.primal_objective_value +=
          info_.workValue_[iCol] * lp_.col_cost_[iCol];
  }
  info_.primal_objective_value *= cost_scale_;
  // Objective value calculation is done using primal values and
  // original costs so offset is vanilla
  info_.primal_objective_value += lp_.offset_;
  // Now have primal objective value
  status_.has_primal_objective_value = true;
  analysis_.simplexTimerStop(ComputePrObjClock);
}

void HEkk::computeDualObjectiveValue(const HighsInt phase) {
  analysis_.simplexTimerStart(ComputeDuObjClock);
  info_.dual_objective_value = 0;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (basis_.nonbasicFlag_[iCol]) {
      const double term = info_.workValue_[iCol] * info_.workDual_[iCol];
      if (term) {
        info_.dual_objective_value +=
            info_.workValue_[iCol] * info_.workDual_[iCol];
      }
    }
  }
  info_.dual_objective_value *= cost_scale_;
  if (phase != 1) {
    // In phase 1 the dual objective has no objective
    // shift. Otherwise, if minimizing the shift is added. If
    // maximizing, workCost (and hence workDual) are negated, so the
    // shift is subtracted. Hence the shift is added according to the
    // sign implied by sense_
    info_.dual_objective_value += ((HighsInt)lp_.sense_) * lp_.offset_;
  }
  // Now have dual objective value
  status_.has_dual_objective_value = true;
  analysis_.simplexTimerStop(ComputeDuObjClock);
}

bool HEkk::rebuildRefactor(HighsInt rebuild_reason) {
  bool refactor = info_.update_count > 0;
  double solution_error = 0;
  if (options_->no_unnecessary_rebuild_refactor) {
    // Consider whether not to refactor in rebuild
    //
    // Must have an INVERT just to consider this!
    assert(status_.has_invert);
    if (rebuild_reason == kRebuildReasonNo ||
        rebuild_reason == kRebuildReasonPossiblyOptimal ||
        rebuild_reason == kRebuildReasonPossiblyPhase1Feasible ||
        rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded ||
        rebuild_reason == kRebuildReasonPossiblyDualUnbounded ||
        rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex) {
      // By default, don't refactor!
      refactor = false;
      // Possibly revise the decision based on accuracy when solving a
      // test system
      const double error_tolerance =
          options_->rebuild_refactor_solution_error_tolerance;
      if (error_tolerance > 0) {
        solution_error = factorSolveError();
        refactor = solution_error > error_tolerance;
      }
    }
  }
  const bool report_refactorization = false;
  if (report_refactorization) {
    const std::string logic = refactor ? "   " : "no   ";
    if (info_.update_count &&
        rebuild_reason != kRebuildReasonSyntheticClockSaysInvert)
      printf(
          "%srefactorization after %4d updates and solution error = %11.4g for "
          "rebuild reason = %s\n",
          logic.c_str(), (int)info_.update_count, solution_error,
          rebuildReason(rebuild_reason).c_str());
  }
  return refactor;
}

HighsInt HEkk::computeFactor() {
  assert(status_.has_nla);
  if (status_.has_fresh_invert) return 0;
  //
  // Perform INVERT
  analysis_.simplexTimerStart(InvertClock);
  const HighsInt rank_deficiency = simplex_nla_.invert();
  analysis_.simplexTimerStop(InvertClock);
  //
  // Set up hot start information
  hot_start_.refactor_info = simplex_nla_.factor_.refactor_info_;
  hot_start_.nonbasicMove = basis_.nonbasicMove_;
  hot_start_.valid = true;

  if (analysis_.analyse_factor_data)
    analysis_.updateInvertFormData(simplex_nla_.factor_);

  HighsInt alt_debug_level = -1;
  if (rank_deficiency) alt_debug_level = kHighsDebugLevelCostly;
  debugNlaCheckInvert("HEkk::computeFactor - original", alt_debug_level);

  if (rank_deficiency) {
    // Have an invertible representation, but of B with column(s)
    // replacements due to singularity. So no (fresh) representation of
    // B^{-1}
    status_.has_invert = false;
    status_.has_fresh_invert = false;
  } else {
    // Now have a representation of B^{-1}, and it is fresh!
    status_.has_invert = true;
    status_.has_fresh_invert = true;
  }
  // Set the update count to zero since the corrected invertible
  // representation may be used for an initial basis. In any case the
  // number of updates shouldn't be positive
  info_.update_count = 0;

  return rank_deficiency;
}

void HEkk::resetSyntheticClock() {
  this->build_synthetic_tick_ = this->simplex_nla_.build_synthetic_tick_;
  this->total_synthetic_tick_ = 0;
}

void HEkk::initialisePartitionedRowwiseMatrix() {
  if (status_.has_ar_matrix) return;
  analysis_.simplexTimerStart(matrixSetupClock);
  ar_matrix_.createRowwisePartitioned(lp_.a_matrix_, &basis_.nonbasicFlag_[0]);
  assert(ar_matrix_.debugPartitionOk(&basis_.nonbasicFlag_[0]));
  analysis_.simplexTimerStop(matrixSetupClock);
  status_.has_ar_matrix = true;
}

void HEkk::setNonbasicMove() {
  const bool have_solution = false;
  // Don't have a simplex basis since nonbasicMove is not set up.

  // Assign nonbasicMove using as much information as is available
  double lower;
  double upper;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  basis_.nonbasicMove_.resize(num_tot);

  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      continue;
    }
    // Nonbasic variable
    if (iVar < lp_.num_col_) {
      lower = lp_.col_lower_[iVar];
      upper = lp_.col_upper_[iVar];
    } else {
      HighsInt iRow = iVar - lp_.num_col_;
      lower = -lp_.row_upper_[iRow];
      upper = -lp_.row_lower_[iRow];
    }
    HighsInt move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        //
        // Determine the bound to set the value to according to, in order of
        // priority
        //
        // 1. Any solution value
        if (have_solution) {
          double midpoint = 0.5 * (lower + upper);
          double value = info_.workValue_[iVar];
          if (value < midpoint) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
        // 2. Bound of original LP that is closer to zero
        if (move == kIllegalMoveValue) {
          if (fabs(lower) < fabs(upper)) {
            move = kNonbasicMoveUp;
          } else {
            move = kNonbasicMoveDn;
          }
        }
      } else {
        // Lower (since upper bound is infinite)
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      move = kNonbasicMoveDn;
    } else {
      // FREE
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iVar] = move;
  }
}

void HEkk::allocateWorkAndBaseArrays() {
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  info_.workCost_.resize(num_tot);
  info_.workDual_.resize(num_tot);
  info_.workShift_.resize(num_tot);

  info_.workLower_.resize(num_tot);
  info_.workUpper_.resize(num_tot);
  info_.workRange_.resize(num_tot);
  info_.workValue_.resize(num_tot);
  info_.workLowerShift_.resize(num_tot);
  info_.workUpperShift_.resize(num_tot);

  // Feel that it should be possible to resize this with in dual
  // solver, and only if Devex is being used, but a pointer to it
  // needs to be set up when constructing HDual
  info_.devex_index_.resize(num_tot);

  info_.baseLower_.resize(lp_.num_row_);
  info_.baseUpper_.resize(lp_.num_row_);
  info_.baseValue_.resize(lp_.num_row_);
}

void HEkk::initialiseLpColBound() {
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    info_.workLower_[iCol] = lp_.col_lower_[iCol];
    info_.workUpper_[iCol] = lp_.col_upper_[iCol];
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
    info_.workLowerShift_[iCol] = 0;
    info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowBound() {
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iCol = lp_.num_col_ + iRow;
    info_.workLower_[iCol] = -lp_.row_upper_[iRow];
    info_.workUpper_[iCol] = -lp_.row_lower_[iRow];
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
    info_.workLowerShift_[iCol] = 0;
    info_.workUpperShift_[iCol] = 0;
  }
}

void HEkk::initialiseCost(const SimplexAlgorithm algorithm,
                          const HighsInt solve_phase, const bool perturb) {
  // Copy the cost
  initialiseLpColCost();
  initialiseLpRowCost();
  info_.costs_perturbed = 0;
  // Primal simplex costs are either from the LP or set specially in phase 1
  if (algorithm == SimplexAlgorithm::kPrimal) return;
  // Dual simplex costs are either from the LP or perturbed
  if (!perturb || info_.dual_simplex_cost_perturbation_multiplier == 0) return;
  // Perturb the original costs, scale down if is too big
  const bool report_cost_perturbation =
      options_->output_flag && analysis_.analyse_simplex_summary_data;
  HighsInt num_original_nonzero_cost = 0;
  if (report_cost_perturbation)
    highsLogUser(options_->log_options, HighsLogType::kInfo,
                 "Cost perturbation for %s\n", lp_.model_name_.c_str());
  double bigc = 0;
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    const double abs_cost = fabs(info_.workCost_[i]);
    bigc = max(bigc, abs_cost);
    if (report_cost_perturbation && abs_cost) num_original_nonzero_cost++;
  }
  const HighsInt pct0 = (100 * num_original_nonzero_cost) / lp_.num_col_;
  double average_cost = 0;
  if (report_cost_perturbation) {
    if (num_original_nonzero_cost) {
      average_cost = bigc / num_original_nonzero_cost;
    } else {
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   "   STRANGE initial workCost has non nonzeros\n");
    }
    highsLogUser(options_->log_options, HighsLogType::kInfo,
                 "   Initially have %" HIGHSINT_FORMAT
                 " nonzero costs (%3" HIGHSINT_FORMAT
                 "%%) with bigc = "
                 "%g "
                 "and average = %g\n",
                 num_original_nonzero_cost, pct0, bigc, average_cost);
  }
  if (bigc > 100) {
    bigc = sqrt(sqrt(bigc));
    if (report_cost_perturbation)
      highsLogUser(options_->log_options, HighsLogType::kInfo,
                   "   Large so set bigc = sqrt(bigc) = %g\n", bigc);
  }

  // If there are few boxed variables, we will just use simple perturbation
  double boxedRate = 0;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt i = 0; i < num_tot; i++)
    boxedRate += (info_.workRange_[i] < 1e30);
  boxedRate /= num_tot;
  if (boxedRate < 0.01) {
    bigc = min(bigc, 1.0);
    if (report_cost_perturbation)
      highsLogUser(options_->log_options, HighsLogType::kInfo,

                   "   small boxedRate (%g) so set bigc = min(bigc, 1.0) = "
                   "%g\n",
                   boxedRate, bigc);
  }
  // Determine the perturbation base
  double base = 5e-7 * bigc;
  if (report_cost_perturbation)
    highsLogUser(options_->log_options, HighsLogType::kInfo,
                 "   Perturbation base = %g\n", base);

  // Now do the perturbation
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    double lower = lp_.col_lower_[i];
    double upper = lp_.col_upper_[i];
    double xpert = (fabs(info_.workCost_[i]) + 1) * base *
                   info_.dual_simplex_cost_perturbation_multiplier *
                   (1 + info_.numTotRandomValue_[i]);
    const double previous_cost = info_.workCost_[i];
    if (lower <= -kHighsInf && upper >= kHighsInf) {
      // Free - no perturb
    } else if (upper >= kHighsInf) {  // Lower
      info_.workCost_[i] += xpert;
    } else if (lower <= -kHighsInf) {  // Upper
      info_.workCost_[i] += -xpert;
    } else if (lower != upper) {  // Boxed
      info_.workCost_[i] += (info_.workCost_[i] >= 0) ? xpert : -xpert;
    } else {
      // Fixed - no perturb
    }
    if (report_cost_perturbation) {
      const double perturbation1 = fabs(info_.workCost_[i] - previous_cost);
      if (perturbation1)
        updateValueDistribution(perturbation1,
                                analysis_.cost_perturbation1_distribution);
    }
  }
  for (HighsInt i = lp_.num_col_; i < num_tot; i++) {
    double perturbation2 = (0.5 - info_.numTotRandomValue_[i]) *
                           info_.dual_simplex_cost_perturbation_multiplier *
                           1e-12;
    info_.workCost_[i] += perturbation2;
    if (report_cost_perturbation) {
      perturbation2 = fabs(perturbation2);
      updateValueDistribution(perturbation2,
                              analysis_.cost_perturbation2_distribution);
    }
  }
  info_.costs_perturbed = 1;
}

void HEkk::initialiseBound(const SimplexAlgorithm algorithm,
                           const HighsInt solve_phase, const bool perturb) {
  initialiseLpColBound();
  initialiseLpRowBound();
  info_.bounds_perturbed = 0;
  // Primal simplex bounds are either from the LP or perturbed
  if (algorithm == SimplexAlgorithm::kPrimal) {
    if (!perturb || info_.primal_simplex_bound_perturbation_multiplier == 0)
      return;
    // Perturb the bounds
    // Determine the smallest and largest finite lower/upper bounds
    HighsInt num_col = lp_.num_col_;
    HighsInt num_row = lp_.num_row_;
    HighsInt num_tot = num_col + num_row;
    double min_abs_lower = kHighsInf;
    double max_abs_lower = -1;
    double min_abs_upper = kHighsInf;
    double max_abs_upper = -1;
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      double abs_lower = fabs(info_.workLower_[iVar]);
      double abs_upper = fabs(info_.workUpper_[iVar]);
      if (abs_lower && abs_lower < kHighsInf) {
        min_abs_lower = min(abs_lower, min_abs_lower);
        max_abs_lower = max(abs_lower, max_abs_lower);
      }
      if (abs_upper && abs_upper < kHighsInf) {
        min_abs_upper = min(abs_upper, min_abs_upper);
        max_abs_upper = max(abs_upper, max_abs_upper);
      }
    }
    // printf(
    //     "Nonzero finite lower bounds in [%9.4g, %9.4g]; upper bounds in "
    //     "[%9.4g, %9.4g]\n",
    //     min_abs_lower, max_abs_lower, min_abs_upper, max_abs_upper);

    const double base =
        info_.primal_simplex_bound_perturbation_multiplier * 5e-7;
    for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
      double lower = info_.workLower_[iVar];
      double upper = info_.workUpper_[iVar];
      const bool fixed = lower == upper;
      // Don't perturb bounds of nonbasic fixed variables as they stay nonbasic
      if (basis_.nonbasicFlag_[iVar] == kNonbasicFlagTrue && fixed) continue;
      double random_value = info_.numTotRandomValue_[iVar];
      if (lower > -kHighsInf) {
        if (lower < -1) {
          lower -= random_value * base * (-lower);
        } else if (lower < 1) {
          lower -= random_value * base;
        } else {
          lower -= random_value * base * lower;
        }
        info_.workLower_[iVar] = lower;
      }
      if (upper < kHighsInf) {
        if (upper < -1) {
          upper += random_value * base * (-upper);
        } else if (upper < 1) {
          upper += random_value * base;
        } else {
          upper += random_value * base * upper;
        }
        info_.workUpper_[iVar] = upper;
      }
      info_.workRange_[iVar] = info_.workUpper_[iVar] - info_.workLower_[iVar];
      if (basis_.nonbasicFlag_[iVar] == kNonbasicFlagFalse) continue;
      // Set values of nonbasic variables
      if (basis_.nonbasicMove_[iVar] > 0) {
        info_.workValue_[iVar] = lower;
      } else if (basis_.nonbasicMove_[iVar] < 0) {
        info_.workValue_[iVar] = upper;
      }
    }
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      HighsInt iVar = basis_.basicIndex_[iRow];
      info_.baseLower_[iRow] = info_.workLower_[iVar];
      info_.baseUpper_[iRow] = info_.workUpper_[iVar];
    }
    info_.bounds_perturbed = 1;
    return;
  }
  // Dual simplex bounds are either from the LP or set to special values in
  // phase
  // 1
  assert(algorithm == SimplexAlgorithm::kDual);
  if (solve_phase == kSolvePhase2) return;

  // The dual objective is the sum of products of primal and dual
  // values for nonbasic variables. For dual simplex phase 1, the
  // primal bounds are set so that when the dual value is feasible, the
  // primal value is set to zero. Otherwise the value is +1/-1
  // according to the required sign of the dual, except for free
  // variables, where the bounds are [-1000, 1000]. Hence the dual
  // objective is the negation of the sum of infeasibilities, unless there are
  // free In Phase 1: change to dual phase 1 bound.
  const double inf = kHighsInf;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iCol = 0; iCol < num_tot; iCol++) {
    if (info_.workLower_[iCol] == -inf && info_.workUpper_[iCol] == inf) {
      // Don't change for row variables: they should never become
      // nonbasic when starting from a logical basis, and no crash
      // should make a free row nonbasic, but could an advanced basis
      // make a free row nonbasic.
      // But what it it happened?
      if (iCol >= lp_.num_col_) continue;
      info_.workLower_[iCol] = -1000,
      info_.workUpper_[iCol] = 1000;  // FREE
    } else if (info_.workLower_[iCol] == -inf) {
      info_.workLower_[iCol] = -1,
      info_.workUpper_[iCol] = 0;  // UPPER
    } else if (info_.workUpper_[iCol] == inf) {
      info_.workLower_[iCol] = 0,
      info_.workUpper_[iCol] = 1;  // LOWER
    } else {
      info_.workLower_[iCol] = 0,
      info_.workUpper_[iCol] = 0;  // BOXED or FIXED
    }
    info_.workRange_[iCol] = info_.workUpper_[iCol] - info_.workLower_[iCol];
  }
}

void HEkk::initialiseLpColCost() {
  double cost_scale_factor = pow(2.0, options_->cost_scale_factor);
  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    info_.workCost_[iCol] =
        (HighsInt)lp_.sense_ * cost_scale_factor * lp_.col_cost_[iCol];
    info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseLpRowCost() {
  for (HighsInt iCol = lp_.num_col_; iCol < lp_.num_col_ + lp_.num_row_;
       iCol++) {
    info_.workCost_[iCol] = 0;
    info_.workShift_[iCol] = 0;
  }
}

void HEkk::initialiseNonbasicValueAndMove() {
  // Initialise workValue and nonbasicMove from nonbasicFlag and
  // bounds, except for boxed variables when nonbasicMove is used to
  // set workValue=workLower/workUpper
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) {
      // Basic variable
      basis_.nonbasicMove_[iVar] = kNonbasicMoveZe;
      continue;
    }
    // Nonbasic variable
    const double lower = info_.workLower_[iVar];
    const double upper = info_.workUpper_[iVar];
    const HighsInt original_move = basis_.nonbasicMove_[iVar];
    double value;
    HighsInt move = kIllegalMoveValue;
    if (lower == upper) {
      // Fixed
      value = lower;
      move = kNonbasicMoveZe;
    } else if (!highs_isInfinity(-lower)) {
      // Finite lower bound so boxed or lower
      if (!highs_isInfinity(upper)) {
        // Finite upper bound so boxed
        if (original_move == kNonbasicMoveUp) {
          // Set at lower
          value = lower;
          move = kNonbasicMoveUp;
        } else if (original_move == kNonbasicMoveDn) {
          // Set at upper
          value = upper;
          move = kNonbasicMoveDn;
        } else {
          // Invalid nonbasicMove: correct and set value at lower
          value = lower;
          move = kNonbasicMoveUp;
        }
      } else {
        // Lower
        value = lower;
        move = kNonbasicMoveUp;
      }
    } else if (!highs_isInfinity(upper)) {
      // Upper
      value = upper;
      move = kNonbasicMoveDn;
    } else {
      // FREE
      value = 0;
      move = kNonbasicMoveZe;
    }
    assert(move != kIllegalMoveValue);
    basis_.nonbasicMove_[iVar] = move;
    info_.workValue_[iVar] = value;
  }
}

void HEkk::pivotColumnFtran(const HighsInt iCol, HVector& col_aq) {
  analysis_.simplexTimerStart(FtranClock);
  col_aq.clear();
  col_aq.packFlag = true;
  lp_.a_matrix_.collectAj(col_aq, iCol, 1);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordBefore(kSimplexNlaFtran, col_aq,
                                    info_.col_aq_density);
  simplex_nla_.ftran(col_aq, info_.col_aq_density,
                     analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaFtran, col_aq);
  HighsInt num_row = lp_.num_row_;
  const double local_col_aq_density = (double)col_aq.count / num_row;
  updateOperationResultDensity(local_col_aq_density, info_.col_aq_density);
  analysis_.simplexTimerStop(FtranClock);
}

void HEkk::unitBtran(const HighsInt iRow, HVector& row_ep) {
  analysis_.simplexTimerStart(BtranClock);
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = iRow;
  row_ep.array[iRow] = 1;
  row_ep.packFlag = true;
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordBefore(kSimplexNlaBtranEp, row_ep,
                                    info_.row_ep_density);
  simplex_nla_.btran(row_ep, info_.row_ep_density,
                     analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaBtranEp, row_ep);
  HighsInt num_row = lp_.num_row_;
  const double local_row_ep_density = (double)row_ep.count / num_row;
  updateOperationResultDensity(local_row_ep_density, info_.row_ep_density);
  analysis_.simplexTimerStop(BtranClock);
}

void HEkk::fullBtran(HVector& buffer) {
  // Performs BTRAN on the buffer supplied. Make sure that
  // buffer.count is large (>lp_.num_row_ to be sure) rather
  // than 0 if the indices of the RHS (and true value of buffer.count)
  // isn't known.
  analysis_.simplexTimerStart(BtranFullClock);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordBefore(kSimplexNlaBtranFull, buffer,
                                    info_.dual_col_density);
  simplex_nla_.btran(buffer, info_.dual_col_density,
                     analysis_.pointer_serial_factor_clocks);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaBtranFull, buffer);
  const double local_dual_col_density = (double)buffer.count / lp_.num_row_;
  updateOperationResultDensity(local_dual_col_density, info_.dual_col_density);
  analysis_.simplexTimerStop(BtranFullClock);
}

void HEkk::choosePriceTechnique(const HighsInt price_strategy,
                                const double row_ep_density,
                                bool& use_col_price,
                                bool& use_row_price_w_switch) {
  // By default switch to column PRICE when pi_p has at least this
  // density
  const double density_for_column_price_switch = 0.75;
  use_col_price = (price_strategy == kSimplexPriceStrategyCol) ||
                  (price_strategy == kSimplexPriceStrategyRowSwitchColSwitch &&
                   row_ep_density > density_for_column_price_switch);
  use_row_price_w_switch =
      price_strategy == kSimplexPriceStrategyRowSwitch ||
      price_strategy == kSimplexPriceStrategyRowSwitchColSwitch;
}

void HEkk::tableauRowPrice(const HVector& row_ep, HVector& row_ap) {
  analysis_.simplexTimerStart(PriceClock);
  const HighsInt solver_num_row = lp_.num_row_;
  const HighsInt solver_num_col = lp_.num_col_;
  const double local_density = 1.0 * row_ep.count / solver_num_row;
  bool use_col_price;
  bool use_row_price_w_switch;
  choosePriceTechnique(info_.price_strategy, local_density, use_col_price,
                       use_row_price_w_switch);
  if (analysis_.analyse_simplex_summary_data) {
    if (use_col_price) {
      const double expected_density = 1;
      analysis_.operationRecordBefore(kSimplexNlaPriceAp, row_ep,
                                      expected_density);
      analysis_.num_col_price++;
    } else if (use_row_price_w_switch) {
      analysis_.operationRecordBefore(kSimplexNlaPriceAp, row_ep,
                                      info_.row_ep_density);
      analysis_.num_row_price_with_switch++;
    } else {
      analysis_.operationRecordBefore(kSimplexNlaPriceAp, row_ep,
                                      info_.row_ep_density);
      analysis_.num_row_price++;
    }
  }
  row_ap.clear();
  if (use_col_price) {
    // Perform column-wise PRICE
    lp_.a_matrix_.priceByColumn(row_ap, row_ep);
  } else if (use_row_price_w_switch) {
    // Perform hyper-sparse row-wise PRICE, but switch if the density of row_ap
    // becomes extreme
    const double switch_density = kHyperPriceDensity;
    ar_matrix_.priceByRowWithSwitch(row_ap, row_ep, info_.row_ap_density, 0,
                                    switch_density);
  } else {
    // Perform hyper-sparse row-wise PRICE
    ar_matrix_.priceByRow(row_ap, row_ep);
  }
  if (use_col_price) {
    // Column-wise PRICE computes components corresponding to basic
    // variables, so zero these by exploiting the fact that, for basic
    // variables, nonbasicFlag[*]=0
    const int8_t* nonbasicFlag = &basis_.nonbasicFlag_[0];
    for (HighsInt iCol = 0; iCol < solver_num_col; iCol++)
      row_ap.array[iCol] *= nonbasicFlag[iCol];
  }
  // Update the record of average row_ap density
  const double local_row_ap_density = (double)row_ap.count / solver_num_col;
  updateOperationResultDensity(local_row_ap_density, info_.row_ap_density);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaPriceAp, row_ap);
  analysis_.simplexTimerStop(PriceClock);
}

void HEkk::fullPrice(const HVector& full_col, HVector& full_row) {
  analysis_.simplexTimerStart(PriceFullClock);
  full_row.clear();
  if (analysis_.analyse_simplex_summary_data) {
    const double expected_density = 1;
    analysis_.operationRecordBefore(kSimplexNlaPriceFull, full_col,
                                    expected_density);
  }
  lp_.a_matrix_.priceByColumn(full_row, full_col);
  if (analysis_.analyse_simplex_summary_data)
    analysis_.operationRecordAfter(kSimplexNlaPriceFull, full_row);
  analysis_.simplexTimerStop(PriceFullClock);
}

void HEkk::computePrimal() {
  analysis_.simplexTimerStart(ComputePrimalClock);
  const HighsInt num_row = lp_.num_row_;
  const HighsInt num_col = lp_.num_col_;
  // Setup a local buffer for the values of basic variables
  HVector primal_col;
  primal_col.setup(num_row);
  primal_col.clear();
  for (HighsInt i = 0; i < num_col + num_row; i++) {
    if (basis_.nonbasicFlag_[i] && info_.workValue_[i] != 0) {
      lp_.a_matrix_.collectAj(primal_col, i, info_.workValue_[i]);
    }
  }
  // It's possible that the buffer has no nonzeros, so performing
  // FTRAN is unnecessary. Not much of a saving, but the zero density
  // looks odd in the analysis!
  if (primal_col.count) {
    simplex_nla_.ftran(primal_col, info_.primal_col_density,
                       analysis_.pointer_serial_factor_clocks);
    const double local_primal_col_density = (double)primal_col.count / num_row;
    updateOperationResultDensity(local_primal_col_density,
                                 info_.primal_col_density);
  }
  for (HighsInt i = 0; i < num_row; i++) {
    HighsInt iCol = basis_.basicIndex_[i];
    info_.baseValue_[i] = -primal_col.array[i];
    info_.baseLower_[i] = info_.workLower_[iCol];
    info_.baseUpper_[i] = info_.workUpper_[iCol];
  }
  // Indicate that the primal infeasiblility information isn't known
  info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;

  analysis_.simplexTimerStop(ComputePrimalClock);
}

void HEkk::computeDual() {
  analysis_.simplexTimerStart(ComputeDualClock);
  // Create a local buffer for the pi vector
  HVector dual_col;
  dual_col.setup(lp_.num_row_);
  dual_col.clear();
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    const double value = info_.workCost_[basis_.basicIndex_[iRow]] +
                         info_.workShift_[basis_.basicIndex_[iRow]];
    if (value) {
      dual_col.index[dual_col.count++] = iRow;
      dual_col.array[iRow] = value;
    }
  }
  // Copy the costs in case the basic costs are all zero
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt i = 0; i < num_tot; i++)
    info_.workDual_[i] = info_.workCost_[i];

  if (dual_col.count) {
    fullBtran(dual_col);
    // Create a local buffer for the values of reduced costs
    HVector dual_row;
    dual_row.setup(lp_.num_col_);
    fullPrice(dual_col, dual_row);
    for (HighsInt i = 0; i < lp_.num_col_; i++)
      info_.workDual_[i] -= dual_row.array[i];
    for (HighsInt i = lp_.num_col_; i < num_tot; i++)
      info_.workDual_[i] -= dual_col.array[i - lp_.num_col_];
  }
  // Indicate that the dual infeasiblility information isn't known
  info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;

  analysis_.simplexTimerStop(ComputeDualClock);
}

void HEkk::computeDualInfeasibleWithFlips() {
  // Computes num/max/sum of dual infeasibliities according to
  // nonbasicMove, using the bounds only to identify free variables
  // and non-boxed. Fixed variables are assumed to have nonbasicMove=0
  // so that no dual infeasibility is counted for them. Indeed, when
  // called from cleanup() at the end of dual phase 1, nonbasicMove
  // relates to the phase 1 bounds, but workLower and workUpper will
  // have been set to phase 2 values!
  const double scaled_dual_feasibility_tolerance =
      options_->dual_feasibility_tolerance;
  // Possibly verify that nonbasicMove is correct for fixed variables
  //  debugFixedNonbasicMove(ekk_instance_);

  HighsInt num_dual_infeasibility = 0;
  double max_dual_infeasibility = 0;
  double sum_dual_infeasibility = 0;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;

  for (HighsInt iVar = 0; iVar < num_tot; iVar++) {
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double lower = info_.workLower_[iVar];
    const double upper = info_.workUpper_[iVar];
    const double dual = info_.workDual_[iVar];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else if (highs_isInfinity(-lower) || highs_isInfinity(upper)) {
      // Not free or boxed: any dual infeasibility is given by value
      // signed by nonbasicMove.
      //
      // For boxed variables, nonbasicMove may have the wrong sign for
      // dual, but nonbasicMove and the primal value can be flipped to
      // achieve dual feasiblility.
      dual_infeasibility = -basis_.nonbasicMove_[iVar] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  info_.num_dual_infeasibilities = num_dual_infeasibility;
  info_.max_dual_infeasibility = max_dual_infeasibility;
  info_.sum_dual_infeasibilities = sum_dual_infeasibility;
}

double HEkk::computeDualForTableauColumn(const HighsInt iVar,
                                         const HVector& tableau_column) {
  const vector<double>& workCost = info_.workCost_;
  const vector<HighsInt>& basicIndex = basis_.basicIndex_;

  double dual = info_.workCost_[iVar];
  for (HighsInt i = 0; i < tableau_column.count; i++) {
    HighsInt iRow = tableau_column.index[i];
    dual -= tableau_column.array[iRow] * workCost[basicIndex[iRow]];
  }
  return dual;
}

bool HEkk::correctDual(HighsInt* free_infeasibility_count) {
  const double tau_d = options_->dual_feasibility_tolerance;
  const double inf = kHighsInf;
  HighsInt workCount = 0;
  double flip_dual_objective_value_change = 0;
  double shift_dual_objective_value_change = 0;
  HighsInt num_flip = 0;
  HighsInt num_shift = 0;
  double sum_flip = 0;
  double sum_shift = 0;
  HighsInt num_shift_skipped = 0;
  const HighsInt num_tot = lp_.num_col_ + lp_.num_row_;
  for (HighsInt i = 0; i < num_tot; i++) {
    if (basis_.nonbasicFlag_[i]) {
      if (info_.workLower_[i] == -inf && info_.workUpper_[i] == inf) {
        // FREE variable
        workCount += (fabs(info_.workDual_[i]) >= tau_d);
      } else if (basis_.nonbasicMove_[i] * info_.workDual_[i] <= -tau_d) {
        if (info_.workLower_[i] != -inf && info_.workUpper_[i] != inf) {
          // Boxed variable = flip
          const HighsInt move = basis_.nonbasicMove_[i];
          flipBound(i);
          double flip = info_.workUpper_[i] - info_.workLower_[i];
          // Negative dual at lower bound (move=1): flip to upper
          // bound so objective contribution is change in value (flip)
          // times dual, being move*flip*dual
          //
          // Positive dual at upper bound (move=-1): flip to lower
          // bound so objective contribution is change in value
          // (-flip) times dual, being move*flip*dual
          double local_dual_objective_change = move * flip * info_.workDual_[i];
          local_dual_objective_change *= cost_scale_;
          flip_dual_objective_value_change += local_dual_objective_change;
          num_flip++;
          sum_flip += fabs(flip);
        } else {
          if (info_.allow_cost_perturbation) {
            // Other variable = shift
            info_.costs_perturbed = 1;
            std::string direction;
            double shift;
            if (basis_.nonbasicMove_[i] == 1) {
              direction = "  up";
              double dual = (1 + random_.fraction()) * tau_d;
              shift = dual - info_.workDual_[i];
              info_.workDual_[i] = dual;
              info_.workCost_[i] = info_.workCost_[i] + shift;
            } else {
              direction = "down";
              double dual = -(1 + random_.fraction()) * tau_d;
              shift = dual - info_.workDual_[i];
              info_.workDual_[i] = dual;
              info_.workCost_[i] = info_.workCost_[i] + shift;
            }
            double local_dual_objective_change = shift * info_.workValue_[i];
            local_dual_objective_change *= cost_scale_;
            shift_dual_objective_value_change += local_dual_objective_change;
            num_shift++;
            sum_shift += fabs(shift);
            highsLogDev(options_->log_options, HighsLogType::kVerbose,
                        "Move %s: cost shift = %g; objective change = %g\n",
                        direction.c_str(), shift, local_dual_objective_change);
          } else {
            // Shifting not permitted
            //
            // Before 07/01/20, these shifts were always done, but
            // doing it after cost perturbation has been removed can
            // lead to cycling when dual unboundedness (=> primal
            // infeasibility) has been detecteed in Phase 2, since the
            // shift removes dual infeasibilities, which are then
            // reinstated after the dual values are recomputed.
            //
            // ToDo: Not shifting leads to dual infeasibilities when an
            // LP is declared to be infeasible. Should go to
            // phase 1 primal simplex to "prove" infeasibility.
            //
            num_shift_skipped++;
          }
        }
      }
    }
  }
  if (num_shift_skipped) {
    highsLogDev(options_->log_options, HighsLogType::kError,
                "correctDual: Missed %d cost shifts\n", num_shift_skipped);
    // todo@Julian assert fails quickly on miplib2017 models momentum1, arki001,
    // glass4, and milo-v12-6-r2-40-1
    assert(!num_shift_skipped);
    return false;
  }
  if (num_flip)
    highsLogDev(options_->log_options, HighsLogType::kVerbose,
                "Performed %" HIGHSINT_FORMAT
                " flip(s): total = %g; objective change = %g\n",
                num_flip, sum_flip, flip_dual_objective_value_change);
  if (num_shift)
    highsLogDev(options_->log_options, HighsLogType::kDetailed,
                "Performed %" HIGHSINT_FORMAT
                " cost shift(s): total = %g; objective change = %g\n",
                num_shift, sum_shift, shift_dual_objective_value_change);
  *free_infeasibility_count = workCount;
  return true;
}

void HEkk::flipBound(const HighsInt iCol) {
  int8_t* nonbasicMove = &basis_.nonbasicMove_[0];
  const int8_t move = nonbasicMove[iCol] = -nonbasicMove[iCol];
  info_.workValue_[iCol] =
      move == 1 ? info_.workLower_[iCol] : info_.workUpper_[iCol];
}

bool HEkk::reinvertOnNumericalTrouble(
    const std::string method_name, double& numerical_trouble_measure,
    const double alpha_from_col, const double alpha_from_row,
    const double numerical_trouble_tolerance) {
  double abs_alpha_from_col = fabs(alpha_from_col);
  double abs_alpha_from_row = fabs(alpha_from_row);
  double min_abs_alpha = min(abs_alpha_from_col, abs_alpha_from_row);
  double abs_alpha_diff = fabs(abs_alpha_from_col - abs_alpha_from_row);
  numerical_trouble_measure = abs_alpha_diff / min_abs_alpha;
  const HighsInt update_count = info_.update_count;
  // Reinvert if the relative difference is large enough, and updates have been
  // performed
  const bool numerical_trouble =
      numerical_trouble_measure > numerical_trouble_tolerance;
  const bool reinvert = numerical_trouble && update_count > 0;
  debugReportReinvertOnNumericalTrouble(method_name, numerical_trouble_measure,
                                        alpha_from_col, alpha_from_row,
                                        numerical_trouble_tolerance, reinvert);
  if (reinvert) {
    // Consider increasing the Markowitz multiplier
    const double current_pivot_threshold = info_.factor_pivot_threshold;
    double new_pivot_threshold = 0;
    if (current_pivot_threshold < kDefaultPivotThreshold) {
      // Threshold is below default value, so increase it
      new_pivot_threshold =
          min(current_pivot_threshold * kPivotThresholdChangeFactor,
              kDefaultPivotThreshold);
    } else if (current_pivot_threshold < kMaxPivotThreshold) {
      // Threshold is below max value, so increase it if few updates have been
      // performed
      if (update_count < 10)
        new_pivot_threshold =
            min(current_pivot_threshold * kPivotThresholdChangeFactor,
                kMaxPivotThreshold);
    }
    if (new_pivot_threshold) {
      highsLogUser(options_->log_options, HighsLogType::kWarning,
                   "   Increasing Markowitz threshold to %g\n",
                   new_pivot_threshold);
      info_.factor_pivot_threshold = new_pivot_threshold;
      simplex_nla_.setPivotThreshold(new_pivot_threshold);
    }
  }
  return reinvert;
}

// The major model updates. Factor calls factor_.update; Matrix
// calls matrix_.update; updatePivots does everything---and is
// called from the likes of HDual::updatePivots
void HEkk::transformForUpdate(HVector* column, HVector* row_ep,
                              const HighsInt variable_in, HighsInt* row_out) {
  simplex_nla_.transformForUpdate(column, row_ep, variable_in, *row_out);
}

void HEkk::updateFactor(HVector* column, HVector* row_ep, HighsInt* iRow,
                        HighsInt* hint) {
  analysis_.simplexTimerStart(UpdateFactorClock);
  simplex_nla_.update(column, row_ep, iRow, hint);
  // Now have a representation of B^{-1}, but it is not fresh
  status_.has_invert = true;
  if (info_.update_count >= info_.update_limit)
    *hint = kRebuildReasonUpdateLimitReached;

  // Determine whether to reinvert based on the synthetic clock
  bool reinvert_syntheticClock =
      this->total_synthetic_tick_ >= this->build_synthetic_tick_;
  const bool performed_min_updates =
      info_.update_count >= kSyntheticTickReinversionMinUpdateCount;
  if (reinvert_syntheticClock && performed_min_updates)
    *hint = kRebuildReasonSyntheticClockSaysInvert;
  analysis_.simplexTimerStop(UpdateFactorClock);
  // Use the next level down for the debug level, since the cost of
  // checking the INVERT every iteration is an order more expensive
  // than checking after factorization.
  HighsInt alt_debug_level = options_->highs_debug_level - 1;
  // Forced expensive debug for development work
  // alt_debug_level = kHighsDebugLevelExpensive;
  HighsDebugStatus debug_status =
      debugNlaCheckInvert("HEkk::updateFactor", alt_debug_level);
  if (debug_status == HighsDebugStatus::kError) {
    *hint = kRebuildReasonPossiblySingularBasis;
  }
}

void HEkk::updatePivots(const HighsInt variable_in, const HighsInt row_out,
                        const HighsInt move_out) {
  analysis_.simplexTimerStart(UpdatePivotsClock);
  HighsInt variable_out = basis_.basicIndex_[row_out];

  // Incoming variable
  basis_.basicIndex_[row_out] = variable_in;
  basis_.nonbasicFlag_[variable_in] = 0;
  basis_.nonbasicMove_[variable_in] = 0;
  info_.baseLower_[row_out] = info_.workLower_[variable_in];
  info_.baseUpper_[row_out] = info_.workUpper_[variable_in];

  // Outgoing variable
  basis_.nonbasicFlag_[variable_out] = 1;
  if (info_.workLower_[variable_out] == info_.workUpper_[variable_out]) {
    info_.workValue_[variable_out] = info_.workLower_[variable_out];
    basis_.nonbasicMove_[variable_out] = 0;
  } else if (move_out == -1) {
    info_.workValue_[variable_out] = info_.workLower_[variable_out];
    basis_.nonbasicMove_[variable_out] = 1;
  } else {
    info_.workValue_[variable_out] = info_.workUpper_[variable_out];
    basis_.nonbasicMove_[variable_out] = -1;
  }
  // Update the dual objective value
  double nwValue = info_.workValue_[variable_out];
  double vrDual = info_.workDual_[variable_out];
  double dl_dual_objective_value = nwValue * vrDual;
  info_.updated_dual_objective_value += dl_dual_objective_value;
  info_.update_count++;
  // Update the number of basic logicals
  if (variable_out < lp_.num_col_) info_.num_basic_logicals++;
  if (variable_in < lp_.num_col_) info_.num_basic_logicals--;
  // No longer have a representation of B^{-1}, and certainly not
  // fresh!
  status_.has_invert = false;
  status_.has_fresh_invert = false;
  // Data are no longer fresh from rebuild
  status_.has_fresh_rebuild = false;
  analysis_.simplexTimerStop(UpdatePivotsClock);
}

void HEkk::updateMatrix(const HighsInt variable_in,
                        const HighsInt variable_out) {
  analysis_.simplexTimerStart(UpdateMatrixClock);
  ar_matrix_.update(variable_in, variable_out, lp_.a_matrix_);
  assert(ar_matrix_.debugPartitionOk(&basis_.nonbasicFlag_[0]));
  analysis_.simplexTimerStop(UpdateMatrixClock);
}

void HEkk::computeSimplexInfeasible() {
  computeSimplexPrimalInfeasible();
  computeSimplexDualInfeasible();
}

void HEkk::computeSimplexPrimalInfeasible() {
  // Computes num/max/sum of primal infeasibliities according to the
  // simplex bounds. This is used to determine optimality in dual
  // phase 1 and dual phase 2, albeit using different bounds in
  // workLower/Upper.
  analysis_.simplexTimerStart(ComputePrIfsClock);
  const double scaled_primal_feasibility_tolerance =
      options_->primal_feasibility_tolerance;
  HighsInt& num_primal_infeasibility = info_.num_primal_infeasibilities;
  double& max_primal_infeasibility = info_.max_primal_infeasibility;
  double& sum_primal_infeasibility = info_.sum_primal_infeasibilities;
  num_primal_infeasibility = 0;
  max_primal_infeasibility = 0;
  sum_primal_infeasibility = 0;

  for (HighsInt i = 0; i < lp_.num_col_ + lp_.num_row_; i++) {
    if (basis_.nonbasicFlag_[i]) {
      // Nonbasic column
      double value = info_.workValue_[i];
      double lower = info_.workLower_[i];
      double upper = info_.workUpper_[i];
      // @primal_infeasibility calculation
      double primal_infeasibility = 0;
      if (value < lower - scaled_primal_feasibility_tolerance) {
        primal_infeasibility = lower - value;
      } else if (value > upper + scaled_primal_feasibility_tolerance) {
        primal_infeasibility = value - upper;
      }
      if (primal_infeasibility > 0) {
        if (primal_infeasibility > scaled_primal_feasibility_tolerance)
          num_primal_infeasibility++;
        max_primal_infeasibility =
            std::max(primal_infeasibility, max_primal_infeasibility);
        sum_primal_infeasibility += primal_infeasibility;
      }
    }
  }
  for (HighsInt i = 0; i < lp_.num_row_; i++) {
    // Basic variable
    double value = info_.baseValue_[i];
    double lower = info_.baseLower_[i];
    double upper = info_.baseUpper_[i];
    // @primal_infeasibility calculation
    double primal_infeasibility = 0;
    if (value < lower - scaled_primal_feasibility_tolerance) {
      primal_infeasibility = lower - value;
    } else if (value > upper + scaled_primal_feasibility_tolerance) {
      primal_infeasibility = value - upper;
    }
    if (primal_infeasibility > 0) {
      if (primal_infeasibility > scaled_primal_feasibility_tolerance)
        num_primal_infeasibility++;
      max_primal_infeasibility =
          std::max(primal_infeasibility, max_primal_infeasibility);
      sum_primal_infeasibility += primal_infeasibility;
    }
  }
  analysis_.simplexTimerStop(ComputePrIfsClock);
}

void HEkk::computeSimplexDualInfeasible() {
  analysis_.simplexTimerStart(ComputeDuIfsClock);
  // Computes num/max/sum of dual infeasibilities in phase 1 and phase
  // 2 according to nonbasicMove. The bounds are only used to identify
  // free variables. Fixed variables are assumed to have
  // nonbasicMove=0 so that no dual infeasibility is counted for them.
  const double scaled_dual_feasibility_tolerance =
      options_->dual_feasibility_tolerance;
  HighsInt& num_dual_infeasibility = info_.num_dual_infeasibilities;
  double& max_dual_infeasibility = info_.max_dual_infeasibility;
  double& sum_dual_infeasibility = info_.sum_dual_infeasibilities;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (HighsInt iCol = 0; iCol < lp_.num_col_ + lp_.num_row_; iCol++) {
    if (!basis_.nonbasicFlag_[iCol]) continue;
    // Nonbasic column
    const double dual = info_.workDual_[iCol];
    const double lower = info_.workLower_[iCol];
    const double upper = info_.workUpper_[iCol];
    double dual_infeasibility = 0;
    if (highs_isInfinity(-lower) && highs_isInfinity(upper)) {
      // Free: any nonzero dual value is infeasible
      dual_infeasibility = fabs(dual);
    } else {
      // Not free: any dual infeasibility is given by the dual value
      // signed by nonbasicMove
      dual_infeasibility = -basis_.nonbasicMove_[iCol] * dual;
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  analysis_.simplexTimerStop(ComputeDuIfsClock);
}

void HEkk::computeSimplexLpDualInfeasible() {
  // Compute num/max/sum of dual infeasibliities according to the
  // bounds of the simplex LP. Assumes that boxed variables have
  // primal variable at the bound corresponding to the sign of the
  // dual so should only be used in dual phase 1 - where it's only
  // used for reporting after rebuilds.
  const double scaled_dual_feasibility_tolerance =
      options_->dual_feasibility_tolerance;
  HighsInt& num_dual_infeasibility =
      analysis_.num_dual_phase_1_lp_dual_infeasibility;
  double& max_dual_infeasibility =
      analysis_.max_dual_phase_1_lp_dual_infeasibility;
  double& sum_dual_infeasibility =
      analysis_.sum_dual_phase_1_lp_dual_infeasibility;
  num_dual_infeasibility = 0;
  max_dual_infeasibility = 0;
  sum_dual_infeasibility = 0;

  for (HighsInt iCol = 0; iCol < lp_.num_col_; iCol++) {
    HighsInt iVar = iCol;
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic column
    const double dual = info_.workDual_[iVar];
    const double lower = lp_.col_lower_[iCol];
    const double upper = lp_.col_upper_[iCol];
    double dual_infeasibility = 0;
    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free: any nonzero dual value is infeasible
        dual_infeasibility = fabs(dual);
      } else {
        // Only lower bounded: a negative dual is infeasible
        dual_infeasibility = -dual;
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded: a positive dual is infeasible
        dual_infeasibility = dual;
      } else {
        // Boxed or fixed: any dual value is feasible
        dual_infeasibility = 0;
      }
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
  for (HighsInt iRow = 0; iRow < lp_.num_row_; iRow++) {
    HighsInt iVar = lp_.num_col_ + iRow;
    if (!basis_.nonbasicFlag_[iVar]) continue;
    // Nonbasic row
    const double dual = -info_.workDual_[iVar];
    const double lower = lp_.row_lower_[iRow];
    const double upper = lp_.row_upper_[iRow];
    double dual_infeasibility = 0;
    if (highs_isInfinity(upper)) {
      if (highs_isInfinity(-lower)) {
        // Free: any nonzero dual value is infeasible
        dual_infeasibility = fabs(dual);
      } else {
        // Only lower bounded: a negative dual is infeasible
        dual_infeasibility = -dual;
      }
    } else {
      if (highs_isInfinity(-lower)) {
        // Only upper bounded: a positive dual is infeasible
        dual_infeasibility = dual;
      } else {
        // Boxed or fixed: any dual value is feasible
        dual_infeasibility = 0;
      }
    }
    if (dual_infeasibility > 0) {
      if (dual_infeasibility >= scaled_dual_feasibility_tolerance)
        num_dual_infeasibility++;
      max_dual_infeasibility =
          std::max(dual_infeasibility, max_dual_infeasibility);
      sum_dual_infeasibility += dual_infeasibility;
    }
  }
}

bool HEkk::sparseLoopStyle(const HighsInt count, const HighsInt dim,
                           HighsInt& to_entry) {
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

void HEkk::invalidatePrimalMaxSumInfeasibilityRecord() {
  info_.max_primal_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_primal_infeasibilities = kHighsIllegalInfeasibilityMeasure;
}

void HEkk::invalidatePrimalInfeasibilityRecord() {
  info_.num_primal_infeasibilities = kHighsIllegalInfeasibilityCount;
  invalidatePrimalMaxSumInfeasibilityRecord();
}

void HEkk::invalidateDualMaxSumInfeasibilityRecord() {
  info_.max_dual_infeasibility = kHighsIllegalInfeasibilityMeasure;
  info_.sum_dual_infeasibilities = kHighsIllegalInfeasibilityMeasure;
}

void HEkk::invalidateDualInfeasibilityRecord() {
  info_.num_dual_infeasibilities = kHighsIllegalInfeasibilityCount;
  invalidateDualMaxSumInfeasibilityRecord();
}

bool HEkk::bailoutOnTimeIterations() {
  if (solve_bailout_) {
    // Bailout has already been decided: check that it's for one of these
    // reasons
    assert(model_status_ == HighsModelStatus::kTimeLimit ||
           model_status_ == HighsModelStatus::kIterationLimit ||
           model_status_ == HighsModelStatus::kObjectiveBound ||
           model_status_ == HighsModelStatus::kObjectiveTarget);
  } else if (timer_->readRunHighsClock() > options_->time_limit) {
    solve_bailout_ = true;
    model_status_ = HighsModelStatus::kTimeLimit;
  } else if (iteration_count_ >= options_->simplex_iteration_limit) {
    solve_bailout_ = true;
    model_status_ = HighsModelStatus::kIterationLimit;
  }
  return solve_bailout_;
}

HighsStatus HEkk::returnFromSolve(const HighsStatus return_status) {
  // Always called before returning from HEkkPrimal/Dual::solve()
  if (solve_bailout_) {
    // If bailout has already been decided: check that it's for one of
    // these reasons
    assert(model_status_ == HighsModelStatus::kTimeLimit ||
           model_status_ == HighsModelStatus::kIterationLimit ||
           model_status_ == HighsModelStatus::kObjectiveBound ||
           model_status_ == HighsModelStatus::kObjectiveTarget);
  }
  // Check that returnFromSolve has not already been called: it should
  // be called exactly once per solve
  assert(!called_return_from_solve_);
  called_return_from_solve_ = true;
  info_.valid_backtracking_basis_ = false;

  // Initialise the status of the primal and dual solutions
  return_primal_solution_status_ = kSolutionStatusNone;
  return_dual_solution_status_ = kSolutionStatusNone;
  // Nothing more is known about the solve after an error return
  if (return_status == HighsStatus::kError) return return_status;

  // Check that an invert exists
  assert(status_.has_invert);

  // Determine a primal and dual solution, removing the effects of
  // perturbations and shifts
  //
  // Unless the solution is optimal, invalidate the infeasibility data
  if (model_status_ != HighsModelStatus::kOptimal) {
    invalidatePrimalInfeasibilityRecord();
    invalidateDualInfeasibilityRecord();
  }
  switch (model_status_) {
    case HighsModelStatus::kOptimal: {
      if (info_.num_primal_infeasibilities) {
        // Optimal - but not to desired primal feasibilit tolerance
        return_primal_solution_status_ = kSolutionStatusInfeasible;
      } else {
        return_primal_solution_status_ = kSolutionStatusFeasible;
      }
      if (info_.num_dual_infeasibilities) {
        // Optimal - but not to desired dual feasibilit tolerance
        return_dual_solution_status_ = kSolutionStatusInfeasible;
      } else {
        return_dual_solution_status_ = kSolutionStatusFeasible;
      }
      break;
    }
    case HighsModelStatus::kInfeasible: {
      // Primal simplex has identified primal infeasibility in phase 1, or
      // dual simplex has identified dual unboundedness in phase 2. In
      // both cases there should be no primal or dual perturbations
      assert(!info_.costs_perturbed && !info_.bounds_perturbed);
      if (exit_algorithm_ == SimplexAlgorithm::kPrimal) {
        // Reset the simplex costs and recompute duals after primal
        // phase 1
        initialiseCost(SimplexAlgorithm::kDual, kSolvePhase2);
        computeDual();
      }
      computeSimplexInfeasible();
      // Primal solution shouldn't be feasible
      assert(info_.num_primal_infeasibilities > 0);
      break;
    }
    case HighsModelStatus::kUnboundedOrInfeasible: {
      // Dual simplex has identified dual infeasibility in phase
      // 1. There should be no dual perturbations
      assert(exit_algorithm_ == SimplexAlgorithm::kDual);
      assert(!info_.costs_perturbed);
      // Reset the simplex bounds and recompute primals
      initialiseBound(SimplexAlgorithm::kDual, kSolvePhase2);
      computePrimal();
      computeSimplexInfeasible();
      // Dual solution shouldn't be feasible
      assert(info_.num_dual_infeasibilities > 0);
      break;
    }
    case HighsModelStatus::kUnbounded: {
      // Primal simplex has identified unboundedness in phase 2. There
      // should be no primal or dual perturbations
      assert(exit_algorithm_ == SimplexAlgorithm::kPrimal);
      assert(!info_.costs_perturbed && !info_.bounds_perturbed);
      computeSimplexInfeasible();
      // Primal solution should be feasible
      assert(info_.num_primal_infeasibilities == 0);
      break;
    }
    case HighsModelStatus::kObjectiveBound:
    case HighsModelStatus::kObjectiveTarget:
    case HighsModelStatus::kTimeLimit:
    case HighsModelStatus::kIterationLimit: {
      // Simplex has bailed out due to reaching the objecive cut-off,
      // time or iteration limit. Could happen anywhere (other than
      // the fist implying dual simplex)
      //
      // Reset the simplex bounds and recompute primals
      initialiseBound(SimplexAlgorithm::kDual, kSolvePhase2);
      initialiseNonbasicValueAndMove();
      computePrimal();
      // Reset the simplex costs and recompute duals
      initialiseCost(SimplexAlgorithm::kDual, kSolvePhase2);
      computeDual();
      computeSimplexInfeasible();
      break;
    }
    default: {
      std::string algorithm_name = "primal";
      if (exit_algorithm_ == SimplexAlgorithm::kDual) algorithm_name = "dual";
      highsLogDev(options_->log_options, HighsLogType::kError,
                  "EKK %s simplex solver returns status %s\n",
                  algorithm_name.c_str(),
                  utilModelStatusToString(model_status_).c_str());
      return HighsStatus::kError;
      break;
    }
  }
  assert(info_.num_primal_infeasibilities >= 0);
  assert(info_.num_dual_infeasibilities >= 0);
  if (info_.num_primal_infeasibilities == 0) {
    return_primal_solution_status_ = kSolutionStatusFeasible;
  } else {
    return_primal_solution_status_ = kSolutionStatusInfeasible;
  }
  if (info_.num_dual_infeasibilities == 0) {
    return_dual_solution_status_ = kSolutionStatusFeasible;
  } else {
    return_dual_solution_status_ = kSolutionStatusInfeasible;
  }
  computePrimalObjectiveValue();
  if (!options_->log_dev_level) {
    const bool force = true;
    analysis_.userInvertReport(force);
  }
  return return_status;
}

double HEkk::computeBasisCondition() {
  HighsInt solver_num_row = lp_.num_row_;
  HighsInt solver_num_col = lp_.num_col_;
  vector<double> bs_cond_x;
  vector<double> bs_cond_y;
  vector<double> bs_cond_z;
  vector<double> bs_cond_w;
  HVector row_ep;
  row_ep.setup(solver_num_row);

  const HighsInt* Astart = &lp_.a_matrix_.start_[0];
  const double* Avalue = &lp_.a_matrix_.value_[0];
  // Compute the Hager condition number estimate for the basis matrix
  const double expected_density = 1;
  bs_cond_x.resize(solver_num_row);
  bs_cond_y.resize(solver_num_row);
  bs_cond_z.resize(solver_num_row);
  bs_cond_w.resize(solver_num_row);
  // x = ones(n,1)/n;
  // y = A\x;
  double mu = 1.0 / solver_num_row;
  double norm_Binv = 0;
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = mu;
  row_ep.clear();
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
    double value = bs_cond_x[r_n];
    if (value) {
      row_ep.index[row_ep.count] = r_n;
      row_ep.array[r_n] = value;
      row_ep.count++;
    }
  }
  for (HighsInt ps_n = 1; ps_n <= 5; ps_n++) {
    row_ep.packFlag = false;
    simplex_nla_.ftran(row_ep, expected_density);

    // zeta = sign(y);
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_y[r_n] = row_ep.array[r_n];
      if (bs_cond_y[r_n] > 0)
        bs_cond_w[r_n] = 1.0;
      else if (bs_cond_y[r_n] < 0)
        bs_cond_w[r_n] = -1.0;
      else
        bs_cond_w[r_n] = 0.0;
    }
    // z=A'\zeta;
    row_ep.clear();
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
      double value = bs_cond_w[r_n];
      if (value) {
        row_ep.index[row_ep.count] = r_n;
        row_ep.array[r_n] = value;
        row_ep.count++;
      }
    }
    row_ep.packFlag = false;
    simplex_nla_.btran(row_ep, expected_density);
    double norm_z = 0.0;
    double ztx = 0.0;
    norm_Binv = 0.0;
    HighsInt argmax_z = -1;
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
      bs_cond_z[r_n] = row_ep.array[r_n];
      double abs_z_v = fabs(bs_cond_z[r_n]);
      if (abs_z_v > norm_z) {
        norm_z = abs_z_v;
        argmax_z = r_n;
      }
      ztx += bs_cond_z[r_n] * bs_cond_x[r_n];
      norm_Binv += fabs(bs_cond_y[r_n]);
    }
    if (norm_z <= ztx) break;
    // x = zeros(n,1);
    // x(fd_i) = 1;
    for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) bs_cond_x[r_n] = 0.0;
    row_ep.clear();
    row_ep.count = 1;
    row_ep.index[0] = argmax_z;
    row_ep.array[argmax_z] = 1.0;
    bs_cond_x[argmax_z] = 1.0;
  }
  double norm_B = 0.0;
  for (HighsInt r_n = 0; r_n < solver_num_row; r_n++) {
    HighsInt vr_n = basis_.basicIndex_[r_n];
    double c_norm = 0.0;
    if (vr_n < solver_num_col)
      for (HighsInt el_n = Astart[vr_n]; el_n < Astart[vr_n + 1]; el_n++)
        c_norm += fabs(Avalue[el_n]);
    else
      c_norm += 1.0;
    norm_B = max(c_norm, norm_B);
  }
  double cond_B = norm_Binv * norm_B;
  return cond_B;
}

void HEkk::initialiseAnalysis() {
  analysis_.setup(lp_name_, lp_, *options_, iteration_count_);
}

std::string HEkk::rebuildReason(const HighsInt rebuild_reason) {
  std::string rebuild_reason_string;
  if (rebuild_reason == kRebuildReasonNo) {
    rebuild_reason_string = "No reason";
  } else if (rebuild_reason == kRebuildReasonUpdateLimitReached) {
    rebuild_reason_string = "Update limit reached";
  } else if (rebuild_reason == kRebuildReasonSyntheticClockSaysInvert) {
    rebuild_reason_string = "Synthetic clock";
  } else if (rebuild_reason == kRebuildReasonPossiblyOptimal) {
    rebuild_reason_string = "Possibly optimal";
  } else if (rebuild_reason == kRebuildReasonPossiblyPhase1Feasible) {
    rebuild_reason_string = "Possibly phase 1 feasible";
  } else if (rebuild_reason == kRebuildReasonPossiblyPrimalUnbounded) {
    rebuild_reason_string = "Possibly primal unbounded";
  } else if (rebuild_reason == kRebuildReasonPossiblyDualUnbounded) {
    rebuild_reason_string = "Possibly dual unbounded";
  } else if (rebuild_reason == kRebuildReasonPossiblySingularBasis) {
    rebuild_reason_string = "Possibly singular basis";
  } else if (rebuild_reason == kRebuildReasonPrimalInfeasibleInPrimalSimplex) {
    rebuild_reason_string = "Primal infeasible in primal simplex";
  } else if (rebuild_reason == kRebuildReasonChooseColumnFail) {
    rebuild_reason_string = "Choose column failure";
  } else {
    rebuild_reason_string = "Unidentified";
    assert(1 == 0);
  }
  return rebuild_reason_string;
}

void HEkk::freezeBasis(HighsInt& frozen_basis_id) {
  assert(this->status_.has_invert);
  frozen_basis_id =
      this->simplex_nla_.freeze(this->basis_, info_.col_aq_density);
}

HighsStatus HEkk::unfreezeBasis(const HighsInt frozen_basis_id) {
  // Check that the ID passed is valid
  const bool valid_id = this->simplex_nla_.frozenBasisIdValid(frozen_basis_id);
  if (!valid_id) return HighsStatus::kError;
  const bool will_have_invert =
      this->simplex_nla_.frozenBasisHasInvert(frozen_basis_id);
  this->simplex_nla_.unfreeze(frozen_basis_id, basis_);
  // The pointers to simplex basis components have changed, so have to
  // tell simplex NLA to refresh the use of the pointer to the basic
  // indices
  this->simplex_nla_.setBasicIndexPointers(&basis_.basicIndex_[0]);
  updateStatus(LpAction::kNewBounds);
  // Indicate whether there is a valid factorization after unfreezing
  this->status_.has_invert = will_have_invert;
  // If there's no valid factorization, then there cannot be a fresh one
  if (!this->status_.has_invert) this->status_.has_fresh_invert = false;
  // Check for consistency
  if (!this->simplex_nla_.update_.valid_) assert(!this->status_.has_invert);
  //  printf("HEkk::unfreezeBasis: basis = ");
  //  for (HighsInt iRow = 0; iRow < (int)basis_.basicIndex_.size(); iRow++)
  //    printf(" %d", (int)basis_.basicIndex_[iRow]);
  //  printf("\n");
  return HighsStatus::kOk;
}

HighsStatus HEkk::frozenBasisAllDataClear() {
  return simplex_nla_.frozenBasisAllDataClear() ? HighsStatus::kOk
                                                : HighsStatus::kError;
}

double HEkk::factorSolveError() {
  // Cheap assessment of factor accuracy.
  //
  // Forms a random solution with at most 50 nonzeros, solves for
  // the corresponding RHS, and then checks the 50 solution values.
  const HighsInt num_col = this->lp_.num_col_;
  const HighsInt num_row = this->lp_.num_row_;
  const HighsSparseMatrix& a_matrix = this->lp_.a_matrix_;
  const vector<HighsInt>& base_index = this->basis_.basicIndex_;
  const HighsSparseMatrix& ar_matrix = this->ar_matrix_;
  HVector btran_rhs;
  HVector ftran_rhs;
  btran_rhs.setup(num_row);
  ftran_rhs.setup(num_row);

  // Solve for a random solution
  HighsRandom random(1);

  ftran_rhs.clear();
  const HighsInt ideal_solution_num_nz = 50;
  HighsInt solution_num_nz = min(ideal_solution_num_nz, (num_row + 1) / 2);
  assert(solution_num_nz > 0);
  vector<double> solution_value;
  vector<HighsInt> solution_index;
  vector<int8_t> solution_nonzero;
  solution_nonzero.assign(num_row, 0);
  for (;;) {
    HighsInt iRow = random.integer(num_row);
    assert(iRow < num_row);
    if (solution_nonzero[iRow]) continue;
    double value = random.fraction();
    solution_value.push_back(value);
    solution_index.push_back(iRow);
    solution_nonzero[iRow] = 1;
    HighsInt iCol = base_index[iRow];
    a_matrix.collectAj(ftran_rhs, iCol, value);
    if ((int)solution_value.size() == solution_num_nz) break;
  }

  btran_rhs.clear();
  vector<double> btran_solution;
  btran_solution.assign(num_row, 0);
  for (HighsInt iX = 0; iX < solution_value.size(); iX++)
    btran_solution[solution_index[iX]] = solution_value[iX];
  vector<double> btran_scattered_rhs;
  btran_scattered_rhs.assign(num_col + num_row, 0);
  for (HighsInt iX = 0; iX < solution_value.size(); iX++) {
    HighsInt iRow = solution_index[iX];
    for (HighsInt iEl = ar_matrix.p_end_[iRow];
         iEl < ar_matrix.start_[iRow + 1]; iEl++) {
      HighsInt iCol = ar_matrix.index_[iEl];
      btran_scattered_rhs[iCol] += ar_matrix.value_[iEl] * solution_value[iX];
    }
    HighsInt iCol = num_col + iRow;
    if (this->basis_.nonbasicFlag_[iCol] == 0)
      btran_scattered_rhs[iCol] = solution_value[iX];
  }
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iCol = base_index[iRow];
    if (btran_scattered_rhs[iCol] == 0) continue;
    btran_rhs.array[iRow] = btran_scattered_rhs[iCol];
    btran_rhs.index[btran_rhs.count++] = iRow;
  }

  const double expected_density = solution_num_nz * info_.col_aq_density;
  ftran(ftran_rhs, expected_density);
  btran(btran_rhs, expected_density);

  double ftran_solution_error = 0;
  for (HighsInt iX = 0; iX < solution_value.size(); iX++)
    ftran_solution_error =
        max(fabs(ftran_rhs.array[solution_index[iX]] - solution_value[iX]),
            ftran_solution_error);
  double btran_solution_error = 0;
  for (HighsInt iX = 0; iX < solution_value.size(); iX++)
    btran_solution_error =
        max(fabs(btran_rhs.array[solution_index[iX]] - solution_value[iX]),
            btran_solution_error);
  double solution_error = max(ftran_solution_error, btran_solution_error);
  return solution_error;
}
