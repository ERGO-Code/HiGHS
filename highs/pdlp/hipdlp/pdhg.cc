/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/pdhg.cc
 * @brief
 */
#include "pdhg.hpp"

#include <cmath>
#include <iostream>
#include <random>
#include <tuple>

#ifdef CUPDLP_GPU
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#endif

#include "defs.hpp"
#include "linalg.hpp"
#include "lp_data/HConst.h"
#include "pdhg_kernels.h"
#include "pdlp/cupdlp/cupdlp.h"  // For pdlpLogging
#include "pdlp_gpu_debug.hpp"
#include "restart.hpp"

#define PDHG_CHECK_INTERVAL 40
static constexpr double kDivergentMovement = 1e10;

using namespace std;

/*
void vecPrint(const std::vector<double>& vec, const char* name) {
  std::cout << name << ": [";
  for (size_t i = 0; i < vec.size(); ++i) {
    std::cout << vec[i];
    if (i < vec.size() - 1) std::cout << ", ";
  }
  std::cout << "]" << std::endl;
}
*/

void PDLPSolver::printConstraintInfo() {
  if (original_lp_ == nullptr) return;

  HighsInt nRows = original_lp_->num_row_;
  HighsInt nCols = original_lp_->num_col_;

  // Count constraint types BEFORE preprocessing
  HighsInt eq_count = 0, leq_count = 0, geq_count = 0, bound_count = 0,
           free_count = 0;

  for (HighsInt i = 0; i < nRows; ++i) {
    bool has_lower = original_lp_->row_lower_[i] > -kHighsInf;
    bool has_upper = original_lp_->row_upper_[i] < kHighsInf;

    if (has_lower && has_upper) {
      if (original_lp_->row_lower_[i] == original_lp_->row_upper_[i]) {
        eq_count++;
      } else {
        bound_count++;
      }
    } else if (has_lower) {
      geq_count++;
    } else if (has_upper) {
      leq_count++;
    } else {
      free_count++;
    }
  }

  // Count variable bound types
  HighsInt var_fixed = 0, var_lower_only = 0, var_upper_only = 0;
  HighsInt var_boxed = 0, var_free = 0;

  for (HighsInt i = 0; i < nCols; ++i) {
    bool has_lower = original_lp_->col_lower_[i] > -kHighsInf;
    bool has_upper = original_lp_->col_upper_[i] < kHighsInf;

    if (has_lower && has_upper) {
      if (original_lp_->col_lower_[i] == original_lp_->col_upper_[i]) {
        var_fixed++;
      } else {
        var_boxed++;
      }
    } else if (has_lower) {
      var_lower_only++;
    } else if (has_upper) {
      var_upper_only++;
    } else {
      var_free++;
    }
  }

  logger_.detailed("=== BEFORE PREPROCESSING ===");
  logger_.detailed("Rows: " + std::to_string(nRows) +
                   ", Cols: " + std::to_string(nCols));
  logger_.detailed("\nConstraint types:");
  logger_.detailed("  Equality constraints (=): " + std::to_string(eq_count));
  logger_.detailed("  One-sided inequality (>=): " + std::to_string(geq_count));
  logger_.detailed("  One-sided inequality (<=): " + std::to_string(leq_count));
  logger_.detailed("  Two-sided inequality: " + std::to_string(bound_count));
  logger_.detailed("  Free constraints: " + std::to_string(free_count));

  logger_.detailed("\nVariable bound types:");
  logger_.detailed("  Fixed variables (l = u): " + std::to_string(var_fixed));
  logger_.detailed("  Boxed variables (l <= x <= u): " +
                   std::to_string(var_boxed));
  logger_.detailed("  Lower bounded only (l <= x): " +
                   std::to_string(var_lower_only));
  logger_.detailed("  Upper bounded only (x <= u): " +
                   std::to_string(var_upper_only));
  logger_.detailed("  Free variables: " + std::to_string(var_free));

  logger_.detailed("\n=== AFTER PREPROCESSING ===");
  logger_.detailed("Rows: " + std::to_string(lp_.num_row_) +
                   ", Cols: " + std::to_string(lp_.num_col_));
  logger_.detailed("Equality rows (first " + std::to_string(num_eq_rows_) +
                   " rows)");
  logger_.detailed("Inequality rows (remaining " +
                   std::to_string(lp_.num_row_ - num_eq_rows_) + " rows)");
  logger_.detailed("Slack variables added: " +
                   std::to_string(lp_.num_col_ - nCols));

  // Count variable bounds in processed LP
  HighsInt proc_var_fixed = 0, proc_var_lower_only = 0, proc_var_upper_only = 0;
  HighsInt proc_var_boxed = 0, proc_var_free = 0;

  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    bool has_lower = lp_.col_lower_[i] > -kHighsInf;
    bool has_upper = lp_.col_upper_[i] < kHighsInf;

    if (has_lower && has_upper) {
      if (lp_.col_lower_[i] == lp_.col_upper_[i]) {
        proc_var_fixed++;
      } else {
        proc_var_boxed++;
      }
    } else if (has_lower) {
      proc_var_lower_only++;
    } else if (has_upper) {
      proc_var_upper_only++;
    } else {
      proc_var_free++;
    }
  }

  logger_.detailed("\nProcessed variable bound types:");
  logger_.detailed("  Fixed variables: " + std::to_string(proc_var_fixed));
  logger_.detailed("  Boxed variables: " + std::to_string(proc_var_boxed));
  logger_.detailed("  Lower bounded only: " +
                   std::to_string(proc_var_lower_only));
  logger_.detailed("  Upper bounded only: " +
                   std::to_string(proc_var_upper_only));
  logger_.detailed("  Free variables: " + std::to_string(proc_var_free));
}

void PDLPSolver::preprocessLp() {
#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockPreprocess);
#endif
  logger_.detailed(
      "Preprocessing LP using cupdlp formulation (slack variables for "
      "bounds)...");

  HighsInt nRows_orig = original_lp_->num_row_;
  HighsInt nCols_orig = original_lp_->num_col_;

  if (original_lp_->a_matrix_.isRowwise()) {
    logger_.info("Original LP matrix must be in column-wise format,");
  }

  lp_.offset_ = original_lp_->offset_;

  HighsInt num_new_cols = 0;
  HighsInt nEqs = 0;
  sense_origin_ = (original_lp_->sense_ == ObjSense::kMinimize) ? 1 : -1;
  constraint_types_.resize(nRows_orig);

  // 1. First pass: Classify constraints and count slack variables needed
  for (HighsInt i = 0; i < nRows_orig; ++i) {
    bool has_lower = original_lp_->row_lower_[i] > -kHighsInf;
    bool has_upper = original_lp_->row_upper_[i] < kHighsInf;

    if (has_lower && has_upper) {
      if (original_lp_->row_lower_[i] == original_lp_->row_upper_[i]) {
        constraint_types_[i] = EQ;
        nEqs++;
      } else {
        constraint_types_[i] = BOUND;
        num_new_cols++;
        nEqs++;
      }
    } else if (has_lower) {
      constraint_types_[i] = GEQ;
    } else if (has_upper) {
      constraint_types_[i] = LEQ;
    } else {
      constraint_types_[i] = FREE;
      num_new_cols++;
      nEqs++;
    }
  }

  // 2. Set new dimensions
  HighsLp& processed_lp = lp_;
  processed_lp.num_col_ = nCols_orig + num_new_cols;
  processed_lp.num_row_ = nRows_orig;
  original_num_col_ = nCols_orig;
  constraint_new_idx_.resize(nRows_orig);

  processed_lp.col_cost_.resize(processed_lp.num_col_);
  processed_lp.col_lower_.resize(processed_lp.num_col_);
  processed_lp.col_upper_.resize(processed_lp.num_col_);
  processed_lp.row_lower_.resize(processed_lp.num_row_);
  processed_lp.row_upper_.resize(processed_lp.num_row_);

  // 3. Determine row permutation: EQ/BOUND/FREE first, then LEQ/GEQ
  HighsInt eq_idx = 0;
  HighsInt ineq_idx = nEqs;
  for (HighsInt i = 0; i < nRows_orig; ++i) {
    if (constraint_types_[i] == EQ || constraint_types_[i] == BOUND ||
        constraint_types_[i] == FREE) {
      constraint_new_idx_[i] = eq_idx++;
    } else {
      constraint_new_idx_[i] = ineq_idx++;
    }
  }

  // 4. Create is_equality_row_ vector based on NEW ordering
  is_equality_row_.resize(nRows_orig, false);
  for (HighsInt i = 0; i < nRows_orig; ++i) {
    HighsInt new_row_idx = constraint_new_idx_[i];
    is_equality_row_[new_row_idx] =
        (constraint_types_[i] == EQ || constraint_types_[i] == BOUND ||
         constraint_types_[i] == FREE);
  }

  // 5. Populate costs and bounds for original and new slack variables
  for (HighsInt i = 0; i < nCols_orig; ++i) {
    processed_lp.col_cost_[i] = original_lp_->col_cost_[i];
    processed_lp.col_lower_[i] = original_lp_->col_lower_[i];
    processed_lp.col_upper_[i] = original_lp_->col_upper_[i];
  }

  HighsInt current_slack_col = nCols_orig;
  for (HighsInt i = 0; i < nRows_orig; ++i) {
    if (constraint_types_[i] == BOUND || constraint_types_[i] == FREE) {
      processed_lp.col_cost_[current_slack_col] = 0.0;
      processed_lp.col_lower_[current_slack_col] = original_lp_->row_lower_[i];
      processed_lp.col_upper_[current_slack_col] = original_lp_->row_upper_[i];
      current_slack_col++;
    }
  }

  // 6. Set the new RHS (row bounds) in PERMUTED order
  for (HighsInt i = 0; i < nRows_orig; ++i) {
    HighsInt new_row_idx = constraint_new_idx_[i];
    switch (constraint_types_[i]) {
      case EQ:
        processed_lp.row_lower_[new_row_idx] = original_lp_->row_lower_[i];
        processed_lp.row_upper_[new_row_idx] = original_lp_->row_upper_[i];
        break;
      case GEQ:
        processed_lp.row_lower_[new_row_idx] = original_lp_->row_lower_[i];
        processed_lp.row_upper_[new_row_idx] = kHighsInf;
        break;
      case LEQ:
        // Becomes -Ax >= -b
        processed_lp.row_lower_[new_row_idx] = -original_lp_->row_upper_[i];
        processed_lp.row_upper_[new_row_idx] = kHighsInf;
        break;
      case BOUND:
      case FREE:
        // Becomes Ax - z = 0
        processed_lp.row_lower_[new_row_idx] = 0.0;
        processed_lp.row_upper_[new_row_idx] = 0.0;
        break;
    }
  }

  // 7. Build the constraint matrix with proper row ordering
  // First, get the original matrix in column-wise format
  HighsSparseMatrix original_matrix = original_lp_->a_matrix_;
  original_matrix.ensureColwise();

  // Create new matrix
  HighsSparseMatrix& processed_matrix = processed_lp.a_matrix_;
  processed_matrix.clear();
  processed_matrix.format_ = MatrixFormat::kColwise;
  processed_matrix.num_col_ = processed_lp.num_col_;
  processed_matrix.num_row_ = processed_lp.num_row_;

  // Reserve space
  processed_matrix.start_.resize(processed_lp.num_col_ + 1);
  processed_matrix.index_.reserve(original_matrix.numNz() + num_new_cols);
  processed_matrix.value_.reserve(original_matrix.numNz() + num_new_cols);

  // Build column by column
  processed_matrix.start_[0] = 0;
  for (HighsInt col = 0; col < nCols_orig; ++col) {
    // For each column, add entries in the new row order
    std::vector<std::pair<HighsInt, double>> col_entries;

    for (HighsInt el = original_matrix.start_[col];
         el < original_matrix.start_[col + 1]; ++el) {
      HighsInt old_row = original_matrix.index_[el];
      HighsInt new_row = constraint_new_idx_[old_row];
      double val = original_matrix.value_[el];

      // Apply scaling for LEQ constraints
      if (constraint_types_[old_row] == LEQ) {
        val = -val;
      }

      col_entries.push_back({new_row, val});
    }

    // Sort by row index to maintain proper sparse matrix format
    std::sort(col_entries.begin(), col_entries.end());

    // Add to matrix
    for (const auto& entry : col_entries) {
      processed_matrix.index_.push_back(entry.first);
      processed_matrix.value_.push_back(entry.second);
    }

    processed_matrix.start_[col + 1] = processed_matrix.index_.size();
  }

  // Add slack variable columns
  current_slack_col = nCols_orig;
  for (HighsInt i = 0; i < nRows_orig; ++i) {
    if (constraint_types_[i] == BOUND || constraint_types_[i] == FREE) {
      HighsInt row_idx = constraint_new_idx_[i];
      processed_matrix.index_.push_back(row_idx);
      processed_matrix.value_.push_back(-1.0);
      processed_matrix.start_[current_slack_col + 1] =
          processed_matrix.index_.size();
      current_slack_col++;
    }
  }

  num_eq_rows_ = nEqs;
  scaling_.passLp(&processed_lp);
  unscaled_processed_lp_ = processed_lp;

  // 8. Compute and store norms of unscaled cost and rhs
  unscaled_c_norm_ = linalg::vector_norm(processed_lp.col_cost_);
  unscaled_rhs_norm_ = linalg::vector_norm(processed_lp.row_lower_);

  logger_.detailed("Preprocessing complete. New dimensions: " +
                   std::to_string(processed_lp.num_row_) + " rows, " +
                   std::to_string(processed_lp.num_col_) + " cols.");
  logger_.detailed(
      "Unscaled norms: ||c|| = " + std::to_string(unscaled_c_norm_) +
      ", ||b|| = " + std::to_string(unscaled_rhs_norm_));

  printConstraintInfo();
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockPreprocess);
#endif
}

PostSolveRetcode PDLPSolver::postprocess(HighsSolution& solution) {
  logger_.detailed("Post-solving the solution...");

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockPostprocess);
#endif
  if (original_lp_ == nullptr) {
    return PostSolveRetcode::DIMENSION_MISMATCH;
  }

  if (x_current_.empty() || y_current_.empty()) {
    return PostSolveRetcode::INVALID_SOLUTION;
  }

  if (x_current_.size() != static_cast<size_t>(lp_.num_col_) ||
      y_current_.size() != static_cast<size_t>(lp_.num_row_)) {
    logger_.info("Solution dimension mismatch: x_current size=" +
                 std::to_string(x_current_.size()) +
                 " vs expected=" + std::to_string(lp_.num_col_) +
                 ", y_current size=" + std::to_string(y_current_.size()) +
                 " vs expected=" + std::to_string(lp_.num_row_));
    return PostSolveRetcode::DIMENSION_MISMATCH;
  }

  // 2. Resize solution object to original dimensions
  solution.col_value.resize(original_num_col_);
  solution.row_value.resize(original_lp_->num_row_);
  solution.col_dual.resize(original_num_col_);
  solution.row_dual.resize(original_lp_->num_row_);

  // 3. Recover Primal Column Values (x)
  for (HighsInt i = 0; i < original_lp_->num_col_; ++i) {
    if (i >= (HighsInt)x_current_.size()) {
      logger_.info("Index " + std::to_string(i) +
                   " out of bounds for x_current_ of size " +
                   std::to_string(x_current_.size()));
      return PostSolveRetcode::DIMENSION_MISMATCH;
    }

    solution.col_value[i] = x_current_[i];

    if (!std::isfinite(solution.col_value[i])) {
      logger_.info("Non-finite primal variable value at index " +
                   std::to_string(i) + ": " +
                   std::to_string(solution.col_value[i]));
      return PostSolveRetcode::NUMERICAL_ERROR;
    }
  }

  double final_primal_objective = original_lp_->offset_;
  for (HighsInt i = 0; i < original_num_col_; ++i) {
    final_primal_objective +=
        original_lp_->col_cost_[i] * solution.col_value[i];
  }
  results_.primal_obj = final_primal_objective;

  // 4. Recover Dual Row Values (y)
  std::vector<double> y_reordered = y_current_;
  for (HighsInt orig_row = 0; orig_row < original_lp_->num_row_; ++orig_row) {
    HighsInt reordered_row = constraint_new_idx_[orig_row];

    // Get the dual value from the reordered position
    double dual_value = y_reordered[reordered_row];

    // Apply sign correction for LEQ constraints
    if (constraint_types_[orig_row] == LEQ) {
      solution.row_dual[orig_row] = -dual_value;
    } else {
      solution.row_dual[orig_row] = dual_value;
    }
  }

  // 5. Recover Primal Row Values (Ax)
  const std::vector<double>& row_scale = scaling_.GetRowScaling();
  const std::vector<double>& col_scale = scaling_.GetColScaling();
  if (scaling_.IsScaled()) {
    // Unscale matrix, costs, and rhs
    for (HighsInt iCol = 0; iCol < unscaled_processed_lp_.num_col_; ++iCol) {
      unscaled_processed_lp_.col_cost_[iCol] /= col_scale[iCol];
      for (HighsInt iEl = unscaled_processed_lp_.a_matrix_.start_[iCol];
           iEl < unscaled_processed_lp_.a_matrix_.start_[iCol + 1]; ++iEl) {
        HighsInt iRow = unscaled_processed_lp_.a_matrix_.index_[iEl];
        unscaled_processed_lp_.a_matrix_.value_[iEl] /=
            (row_scale[iRow] * col_scale[iCol]);
      }
    }
  }

  std::vector<double> ax_original(original_lp_->num_row_, 0.0);

  // Get original matrix in column-wise format
  HighsSparseMatrix orig_matrix = original_lp_->a_matrix_;
  orig_matrix.ensureColwise();

  // Compute Ax using only the original columns (not slack variables)
  for (HighsInt col = 0; col < original_num_col_; ++col) {
    double x_val = x_current_[col];  // Use unscaled x values

    for (HighsInt el = orig_matrix.start_[col];
         el < orig_matrix.start_[col + 1]; ++el) {
      HighsInt row = orig_matrix.index_[el];
      double a_val = orig_matrix.value_[el];
      ax_original[row] += a_val * x_val;
    }
  }

  // Now ax_original contains the correct row activity values
  for (HighsInt orig_row = 0; orig_row < original_lp_->num_row_; ++orig_row) {
    solution.row_value[orig_row] = ax_original[orig_row];
  }

  // 6. Recover Dual Column Values (Reduced Costs)
  // The duals on the variable bounds l <= x <= u are the reduced costs.
  // In the PDLP framework, these are given by dSlackPos - dSlackNeg.
  for (HighsInt i = 0; i < original_num_col_; ++i) {
    if (i >= (HighsInt)dSlackPos_.size() || i >= (HighsInt)dSlackNeg_.size()) {
      logger_.info("Index " + std::to_string(i) +
                   " out of bounds for dSlackPos/Neg of size " +
                   std::to_string(dSlackPos_.size()));
      return PostSolveRetcode::DIMENSION_MISMATCH;
    }
    solution.col_dual[i] = dSlackPos_[i] - dSlackNeg_[i];
    solution.col_dual[i] *= sense_origin_;
  }

  solution.value_valid = true;  // to do
  solution.dual_valid = true;
  logger_.detailed("Post-solve complete.");

#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockPostprocess);
#endif
  return PostSolveRetcode::OK;
}

void PDLPSolver::solve(std::vector<double>& x, std::vector<double>& y) {
#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockSolve);
#endif

  // 1. Initialization and Setup
#ifdef CUPDLP_GPU
  setupGpu();
#endif
  initializeStepSizes();
  initialize();  // Resets vectors, caches, and sets initial x_current,
                 // y_current
  best_primal_weight_ = primal_weight_;
  best_primal_dual_residual_gap_ = std::numeric_limits<double>::infinity();

  // Initialize internal solver state
  restart_scheme_.passParams(&params_);
  restart_scheme_.Initialize(results_);

  // Copy initial input x,y to internal state
  x_current_ = x;
  y_current_ = y;

  if (params_.use_halpern_restart) {
    x_anchor_ = x_current_;
    y_anchor_ = y_current_;
    halpern_iteration_ = 0;
#ifdef CUPDLP_GPU
    CUDA_CHECK(cudaMemcpy(d_x_anchor_, d_x_current_,
                          lp_.num_col_ * sizeof(double),
                          cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_y_anchor_, d_y_current_,
                          lp_.num_row_ * sizeof(double),
                          cudaMemcpyDeviceToDevice));
#endif
  }

  // Pre-calculate Ax and ATy for current iterate
#ifdef CUPDLP_GPU
  CUDA_CHECK(cudaMemset(d_x_sum_, 0, lp_.num_col_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_y_sum_, 0, lp_.num_row_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_x_avg_, 0, lp_.num_col_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_y_avg_, 0, lp_.num_row_ * sizeof(double)));
  sum_weights_gpu_ = 0.0;
  CUDA_CHECK(cudaMemcpy(d_x_current_, x.data(), lp_.num_col_ * sizeof(double),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_y_current_, y.data(), lp_.num_row_ * sizeof(double),
                        cudaMemcpyHostToDevice));
  linalgGpuAx(d_x_current_, d_ax_current_);
  linalgGpuATy(d_y_current_, d_aty_current_);
#else
  linalg::project_bounds(lp_, x_current_);
  linalg::Ax(lp_, x_current_, Ax_cache_);
  linalg::ATy(lp_, y_current_, ATy_cache_);
#endif

  logger_.print_iteration_header();

  // termination_status is not known. Setting it to NOTSET rather than
  // OPTIMAL means that if it's not set elsewhere then optimality is
  // not returned in error
  TerminationStatus termination_status = TerminationStatus::NOTSET;
  bool do_restart = false;
  double last_trial_fpe = std::numeric_limits<double>::infinity();

  final_iter_count_ = 0;

#ifdef CUPDLP_GPU
  bool graph_created = false;
  cudaGraphExec_t graphExec = nullptr;
#endif

  // 2. Main cuPDLPx-style Loop
  while (final_iter_count_ < params_.max_iterations) {
    // Check global time limit
    if (highs_timer_p->read() > params_.time_limit) {
      logger_.info("Time limit reached.");
      termination_status = TerminationStatus::TIMEOUT;
      break;
    }

    // -- Step 1 (Major, isolated for restart-FPE check) --
#ifdef CUPDLP_GPU
    CUDA_CHECK(cudaMemcpyAsync(d_primal_step_size_, &stepsize_.primal_step,
                               sizeof(double), cudaMemcpyHostToDevice,
                               gpu_stream_));
    CUDA_CHECK(cudaMemcpyAsync(d_dual_step_size_, &stepsize_.dual_step,
                               sizeof(double), cudaMemcpyHostToDevice,
                               gpu_stream_));
    CUDA_CHECK(cudaMemcpyAsync(d_halpern_iteration_, &halpern_iteration_,
                               sizeof(int), cudaMemcpyHostToDevice,
                               gpu_stream_));
    performHalpernPdhgStepGpu(true, 1);
#else
    performHalpernPdhgStep(true, 1);
#endif

    if (do_restart) {
#ifdef CUPDLP_GPU
      fpe_ = computeFixedPointErrorGpu();
#else
      fpe_ = computeFixedPointError();
#endif
      initial_fpe_ = fpe_;
      // if (DEBUG_MODE) std::cout << "[DEBUG] iter=" << final_iter_count_ << "
      // | Restarted with initial_fpe_ = " << std::setprecision(10) <<
      // initial_fpe_ << std::endl;
      do_restart = false;
    }

    // copy back x from gpu and print out
#ifdef CUPDLP_GPU
    if (final_iter_count_ == 9) {
      std::vector<double> x_next_host(lp_.num_col_);
      CUDA_CHECK(cudaMemcpy(x_next_host.data(), d_pdhg_primal_,
                            lp_.num_col_ * sizeof(double),
                            cudaMemcpyDeviceToHost));
      std::cout << "x_next (before 9): ";
      for (int i = 0; i < 3; i++)
        std::cout << std::fixed << std::setprecision(10) << x_next_host[i]
                  << " ";
      std::cout << "\n";
    }
#endif

    // -- Steps 2 to PDHG_CHECK_INTERVAL - 1 (Minor) --
#ifdef CUPDLP_GPU
    if (!graph_created) {
      CUDA_CHECK(
          cudaStreamBeginCapture(gpu_stream_, cudaStreamCaptureModeGlobal));

      for (int i = 2; i <= PDHG_CHECK_INTERVAL - 1; i++) {
        performHalpernPdhgStepGpu(false, i);
      }
      performHalpernPdhgStepGpu(true, PDHG_CHECK_INTERVAL);

      cudaGraph_t graph;
      CUDA_CHECK(cudaStreamEndCapture(gpu_stream_, &graph));
      CUDA_CHECK(cudaGraphInstantiate(&graphExec, graph, NULL, NULL, 0));
      CUDA_CHECK(cudaGraphDestroy(graph));
      graph_created = true;
    }

    CUDA_CHECK(cudaGraphLaunch(graphExec, gpu_stream_));
#else
    for (int i = 2; i <= PDHG_CHECK_INTERVAL - 1; i++) {
      performHalpernPdhgStep(false, i);
    }
    performHalpernPdhgStep(true, PDHG_CHECK_INTERVAL);
#endif

#ifdef CUPDLP_GPU
    if (final_iter_count_ == 9) {
      std::vector<double> x_next_host(lp_.num_col_);
      std::vector<double> x_current_host(lp_.num_col_);
      CUDA_CHECK(cudaMemcpy(x_current_host.data(), d_x_current_,
                            lp_.num_col_ * sizeof(double),
                            cudaMemcpyDeviceToHost));
      CUDA_CHECK(cudaMemcpy(x_next_host.data(), d_pdhg_primal_,
                            lp_.num_col_ * sizeof(double),
                            cudaMemcpyDeviceToHost));
      std::cout << "x_next (after 9): ";
      for (int i = 0; i < 3; i++)
        std::cout << std::fixed << std::setprecision(10) << x_next_host[i]
                  << " ";
      std::cout << "\n";
      std::cout << "x_current (after 9): ";
      for (int i = 0; i < 3; i++)
        std::cout << std::fixed << std::setprecision(10) << x_current_host[i]
                  << " ";
      std::cout << "\n";
    }
#endif

// Compute Error for Restart Check
#ifdef CUPDLP_GPU
    fpe_ = computeFixedPointErrorGpu();
#else
    fpe_ = computeFixedPointError();
#endif
    // if (DEBUG_MODE)    std::cout << "[DEBUG] iter=" << final_iter_count_ << "
    // | FPE at check: " << std::setprecision(10) << fpe_ << std::endl;

    halpern_iteration_ += PDHG_CHECK_INTERVAL;
    final_iter_count_ += PDHG_CHECK_INTERVAL;

    // -- Convergence Check --
    TerminationStatus check_status;
    bool should_terminate =
        runConvergenceCheckAndRestart(final_iter_count_, x, y, check_status);

    if (should_terminate) {
      termination_status = check_status;
      break;
    }

    // -- Check Adaptive Restart (Mimicking `should_do_adaptive_restart` in
    // solver.cu) --
    if (final_iter_count_ == PDHG_CHECK_INTERVAL) {
      do_restart = true;
    } else if (final_iter_count_ > PDHG_CHECK_INTERVAL) {
      if (fpe_ <= restart_scheme_.GetSufficientDecayFactor() * initial_fpe_) {
        do_restart = true;
      }
      if (fpe_ <= restart_scheme_.GetNecessaryDecayFactor() * initial_fpe_) {
        if (fpe_ > last_trial_fpe) {
          do_restart = true;
        }
      }
      if (halpern_iteration_ >=
          restart_scheme_.GetArtificialRestartThreshold() * final_iter_count_) {
        do_restart = true;
      }
    }
    last_trial_fpe = fpe_;

    // -- Perform Restart --
    if (do_restart) {
      if (params_.step_size_strategy == StepSizeStrategy::PID) {
        updatePrimalWeightAtRestart(results_);
        // if (DEBUG_MODE)        std::cout << "[DEBUG] iter=" <<
        // final_iter_count_ << " | Updated primal_weight_ to " <<
        // primal_weight_ << std::endl;
      }
#ifdef CUPDLP_GPU
      CUDA_CHECK(cudaMemcpy(d_x_anchor_, d_pdhg_primal_,
                            lp_.num_col_ * sizeof(double),
                            cudaMemcpyDeviceToDevice));
      CUDA_CHECK(cudaMemcpy(d_y_anchor_, d_pdhg_dual_,
                            lp_.num_row_ * sizeof(double),
                            cudaMemcpyDeviceToDevice));
      CUDA_CHECK(cudaMemcpy(d_x_current_, d_pdhg_primal_,
                            lp_.num_col_ * sizeof(double),
                            cudaMemcpyDeviceToDevice));
      CUDA_CHECK(cudaMemcpy(d_y_current_, d_pdhg_dual_,
                            lp_.num_row_ * sizeof(double),
                            cudaMemcpyDeviceToDevice));
#else
      x_anchor_ = x_next_;
      y_anchor_ = y_next_;
      x_current_ = x_next_;
      y_current_ = y_next_;
      linalg::Ax(lp_, x_current_, Ax_cache_);
      linalg::ATy(lp_, y_current_, ATy_cache_);
#endif

      halpern_iteration_ = 0;
      last_trial_fpe = std::numeric_limits<double>::infinity();
    }
  }

#ifdef CUPDLP_GPU
  if (graphExec) {
    CUDA_CHECK(cudaGraphExecDestroy(graphExec));
  }
#endif

  // 3. Loop Finished
  if (termination_status == TerminationStatus::NOTSET) {
    logger_.info("Iteration limit reached without convergence");
    termination_status = TerminationStatus::MAXITER;
  }

  solveReturn(termination_status);
}

double PDLPSolver::computeFixedPointError() {
  double primal_norm_sq = 0.0;
  double dual_norm_sq = 0.0;
  double cross_term = 0.0;

  std::vector<double> delta_x(lp_.num_col_);
  std::vector<double> delta_y(lp_.num_row_);

  for (int i = 0; i < lp_.num_col_; ++i) {
    delta_x[i] = x_next_[i] - reflected_x_[i];
    primal_norm_sq += delta_x[i] * delta_x[i];
  }

  for (int i = 0; i < lp_.num_row_; ++i) {
    delta_y[i] = y_next_[i] - reflected_y_[i];
    dual_norm_sq += delta_y[i] * delta_y[i];
  }

  std::vector<double> AT_delta_y(lp_.num_col_, 0.0);
  linalg::ATy(lp_, delta_y, AT_delta_y);

  for (int i = 0; i < lp_.num_col_; ++i) {
    cross_term += delta_x[i] * AT_delta_y[i];
  }

  double movement =
      primal_norm_sq * params_.omega + dual_norm_sq / params_.omega;
  double interaction = 2.0 * params_.eta * cross_term;

  return std::sqrt(std::max(0.0, movement + interaction));
}

#ifdef CUPDLP_GPU
double PDLPSolver::computeFixedPointErrorGpu() {
  double alpha_minus_one = -1.0;

  // 1. delta_x = x_next_ - reflected_x_
  // (Assuming d_pdhg_primal_ maps to x_next_ and d_x_next_ is used as
  // reflected_x_ in your minor/major steps)
  CUDA_CHECK(cudaMemcpy(d_delta_x_, d_pdhg_primal_,
                        a_num_cols_ * sizeof(double),
                        cudaMemcpyDeviceToDevice));
  CUBLAS_CHECK(cublasDaxpy(cublas_handle_, a_num_cols_, &alpha_minus_one,
                           d_x_next_, 1, d_delta_x_, 1));

  // 2. delta_y = y_next_ - reflected_y_
  CUDA_CHECK(cudaMemcpy(d_delta_y_, d_pdhg_dual_, a_num_rows_ * sizeof(double),
                        cudaMemcpyDeviceToDevice));
  CUBLAS_CHECK(cublasDaxpy(cublas_handle_, a_num_rows_, &alpha_minus_one,
                           d_y_next_, 1, d_delta_y_, 1));

  // 3. AT_delta_y = A^T * delta_y
  linalgGpuATy(d_delta_y_, d_AT_delta_y_);

  // 4. Compute norms and cross term
  double primal_norm = 0.0, dual_norm = 0.0, cross_term = 0.0;

  CUBLAS_CHECK(
      cublasDnrm2(cublas_handle_, a_num_cols_, d_delta_x_, 1, &primal_norm));
  CUBLAS_CHECK(
      cublasDnrm2(cublas_handle_, a_num_rows_, d_delta_y_, 1, &dual_norm));
  CUBLAS_CHECK(cublasDdot(cublas_handle_, a_num_cols_, d_delta_x_, 1,
                          d_AT_delta_y_, 1, &cross_term));

  double primal_norm_sq = primal_norm * primal_norm;
  double dual_norm_sq = dual_norm * dual_norm;

  double movement =
      primal_norm_sq * params_.omega + dual_norm_sq / params_.omega;
  double interaction = 2.0 * params_.eta * cross_term;

  return std::sqrt(std::max(0.0, movement + interaction));
}
#endif

bool PDLPSolver::runConvergenceCheckAndRestart(size_t iter,
                                               std::vector<double>& output_x,
                                               std::vector<double>& output_y,
                                               TerminationStatus& status) {
  // 1. Compute Average Iterate (GPU or CPU)
#ifdef CUPDLP_GPU
  computeAverageIterateGpu();
#else
  computeAverageIterate(Ax_avg_, ATy_avg_);
#endif

  SolverResults current_results;
  SolverResults average_results;
  bool current_converged = false;
  bool average_converged = false;

  // 2. Check Convergence (GPU or CPU)
#ifdef CUPDLP_GPU
  if (params_.use_halpern_restart && d_pdhg_primal_ != nullptr) {
    // For Halpern, "current iterate" convergence uses the PDHG projected
    // iterates Need Ax and ATy of pdhg iterates:
    linalgGpuAx(d_pdhg_primal_, d_ax_next_);  // reuse d_ax_next_ as scratch
    linalgGpuATy(d_pdhg_dual_, d_aty_next_);  // reuse d_aty_next_ as scratch
    current_converged = checkConvergenceGpu(
        iter, d_pdhg_primal_, d_pdhg_dual_, d_ax_next_, d_aty_next_,
        params_.tolerance, current_results, "[L-GPU]", d_dSlackPos_,
        d_dSlackNeg_,
        true);  // <-- TRUE: Use the Halpern slack saved in d_dSlackPos_
  } else {
    current_converged =
        checkConvergenceGpu(iter, d_x_current_, d_y_current_, d_ax_current_,
                            d_aty_current_, params_.tolerance, current_results,
                            "[L-GPU]", d_dSlackPos_, d_dSlackNeg_,
                            false);  // <-- FALSE: Recompute standard slack
    average_converged = checkConvergenceGpu(
        iter, d_x_avg_, d_y_avg_, d_ax_avg_, d_aty_avg_, params_.tolerance,
        average_results, "[A-GPU]", d_dSlackPosAvg_, d_dSlackNegAvg_,
        false);  // <-- FALSE: Recompute standard slack
  }

#else
  if (params_.use_halpern_restart && iter > 0) {
    std::vector<double> Ax_pdhg(lp_.num_row_);
    std::vector<double> ATy_pdhg(lp_.num_col_);
    linalg::Ax(lp_, x_next_, Ax_pdhg);
    linalg::ATy(lp_, y_next_, ATy_pdhg);
    current_converged = checkConvergence(
        iter, x_next_, y_next_, Ax_pdhg, ATy_pdhg, params_.tolerance,
        current_results, "[L]", dSlackPos_, dSlackNeg_);
  } else {
    current_converged = checkConvergence(
        iter, x_current_, y_current_, Ax_cache_, ATy_cache_, params_.tolerance,
        current_results, "[L]", dSlackPos_, dSlackNeg_);
  }
  average_converged = checkConvergence(iter, x_avg_, y_avg_, Ax_avg_, ATy_avg_,
                                       params_.tolerance, average_results,
                                       "[A]", dSlackPosAvg_, dSlackNegAvg_);
#endif

  // Determine whether to log iterations
  bool iteration_log = logger_.getConsoleLevel() >= LogLevel::kDetailed;
  double time_now = highs_timer_p->read();
  iteration_log = time_now > last_logger_time + kHipdlpLoggerFrequency;
  if (iteration_log) {
    logger_.print_iteration_stats(iter, average_results, current_eta_,
                                  time_now);
    last_logger_time = time_now;
  }
  // 3. Handle Convergence Success
  if (current_converged || average_converged) {
    bool prefer_avg = average_converged;

    // Copy result to output vectors
#ifdef CUPDLP_GPU
    double* src_x = prefer_avg ? d_x_avg_
                               : (params_.use_halpern_restart ? d_pdhg_primal_
                                                              : d_x_current_);
    double* src_y = prefer_avg ? d_y_avg_
                               : (params_.use_halpern_restart ? d_pdhg_dual_
                                                              : d_y_current_);
    double* src_sp = prefer_avg ? d_dSlackPosAvg_ : d_dSlackPos_;
    double* src_sn = prefer_avg ? d_dSlackNegAvg_ : d_dSlackNeg_;

    CUDA_CHECK(cudaMemcpy(output_x.data(), src_x, lp_.num_col_ * sizeof(double),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(output_y.data(), src_y, lp_.num_row_ * sizeof(double),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dSlackPos_.data(), src_sp,
                          lp_.num_col_ * sizeof(double),
                          cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(dSlackNeg_.data(), src_sn,
                          lp_.num_col_ * sizeof(double),
                          cudaMemcpyDeviceToHost));
#else
    if (prefer_avg) {
      output_x = x_avg_;
      output_y = y_avg_;
      dSlackPos_ = dSlackPosAvg_;
      dSlackNeg_ = dSlackNegAvg_;
    } else if (params_.use_halpern_restart) {
      output_x = x_next_;
      output_y = y_next_;
    } else {
      output_x = x_current_;
      output_y = y_current_;
    }
#endif

    final_iter_count_ = iter;
    results_ = prefer_avg ? average_results : current_results;
    logger_.info((prefer_avg ? "Average" : "Current") +
                 std::string(" solution converged"));

    status = TerminationStatus::OPTIMAL;
    return true;  // Stop
  }

  results_ = current_results;
  return false;
}

// ----------------------------------------------------------------------------
// Helper: Perform PDHG Step
// ----------------------------------------------------------------------------
void PDLPSolver::performPdhgStep() {
#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockIterateUpdate);
#endif

  // CPU/GPU agnostic logic for step strategy
  switch (params_.step_size_strategy) {
    case StepSizeStrategy::FIXED:
      updateIteratesFixed();
      break;
    case StepSizeStrategy::ADAPTIVE:
      updateIteratesAdaptive();
      break;
    case StepSizeStrategy::MALITSKY_POCK:
      // Note: Error handling for MP failure simplified for brevity
      updateIteratesMalitskyPock(false);
      break;
    case StepSizeStrategy::PID:
      assert(111 == 999);
  }

#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockIterateUpdate);
#endif
}

// ----------------------------------------------------------------------------
// Helper: Perform Halpern Step
// ----------------------------------------------------------------------------
void PDLPSolver::performHalpernPdhgStep(bool is_major, int k_offset) {
  double primal_step = stepsize_.primal_step;
  double dual_step = stepsize_.dual_step;
  double rho = params_.halpern_gamma;

  int current_k = halpern_iteration_ + k_offset;
  double w = static_cast<double>(current_k) / (current_k + 1.0);

  // --- STEP 1: Primal projection + reflection ---
  if (is_major) {
    if (halpern_dual_slack_next_.size() != static_cast<size_t>(lp_.num_col_)) {
      halpern_dual_slack_next_.resize(lp_.num_col_, 0.0);
    }
    halpern_dual_slack_next_valid_ = true;
  }
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    double temp =
        x_current_[i] - primal_step * (lp_.col_cost_[i] - ATy_cache_[i]);
    double temp_proj =
        linalg::project_box(temp, lp_.col_lower_[i], lp_.col_upper_[i]);
    if (is_major) {
      x_next_[i] = temp_proj;  // pdhg_primal (for convergence check)
      halpern_dual_slack_next_[i] = (temp_proj - temp) / primal_step;
    }
    reflected_x_[i] = 2.0 * temp_proj - x_current_[i];
  }

  // --- STEP 2: SpMV on reflected primal ---
  linalg::Ax(lp_, reflected_x_, Ax_next_);

  // --- STEP 3: Dual projection + reflection ---
  for (HighsInt j = 0; j < lp_.num_row_; j++) {
    double temp = y_current_[j] / dual_step - Ax_next_[j];
    double lower_proj = -lp_.row_upper_[j];
    double upper_proj = -lp_.row_lower_[j];
    double temp_proj = std::max(lower_proj, std::min(temp, upper_proj));

    double pdhg_dual_j = (temp - temp_proj) * dual_step;

    if (is_major) {
      y_next_[j] = pdhg_dual_j;
    }
    reflected_y_[j] = 2.0 * pdhg_dual_j - y_current_[j];
  }

  // --- STEP 4: Halpern blend ---
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    double blended = rho * reflected_x_[i] + (1.0 - rho) * x_current_[i];
    x_current_[i] = w * blended + (1.0 - w) * x_anchor_[i];
  }
  for (HighsInt j = 0; j < lp_.num_row_; j++) {
    double blended = rho * reflected_y_[j] + (1.0 - rho) * y_current_[j];
    y_current_[j] = w * blended + (1.0 - w) * y_anchor_[j];
  }

  // --- STEP 5: Recompute ATy ---
  linalg::ATy(lp_, y_current_, ATy_cache_);
}

// ----------------------------------------------------------------------------
// Helper: Accumulate Averages
// ----------------------------------------------------------------------------
void PDLPSolver::accumulateAverages(size_t iter) {
#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockAverageIterate);
#endif

  HighsInt inner_iter = iter - restart_scheme_.GetLastRestartIter();
  double dMeanStepSize = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);

#ifdef CUPDLP_GPU
  if (params_.use_halpern_restart) {
    // If Halpern, we average the 'current' blended iterate
    launchKernelUpdateAverages_wrapper(d_x_sum_, d_y_sum_, d_x_current_,
                                       d_y_current_, dMeanStepSize,
                                       lp_.num_col_, lp_.num_row_, gpu_stream_);
    sum_weights_gpu_ += dMeanStepSize;
  } else {
    // If standard, we average the 'next' iterate computed by PDHG
    updateAverageIteratesGpu(inner_iter);
  }
#else
  if (params_.use_halpern_restart) {
    updateAverageIterates(x_current_, y_current_, inner_iter);
  } else {
    updateAverageIterates(x_next_, y_next_, inner_iter);
  }
#endif

#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockAverageIterate);
#endif
}

// ----------------------------------------------------------------------------
// Helper: Prepare for Next Iteration
// ----------------------------------------------------------------------------
void PDLPSolver::prepareNextIteration() {
  if (params_.use_halpern_restart) {
    // Halpern update already wrote into x_current_, no swap needed
    // Just ensure caches are correct if not done in PerformHalpernStep
    return;
  } else {
    // Standard PDHG: x_current becomes x_next
#ifdef CUPDLP_GPU
    CUDA_CHECK(cudaMemcpy(d_x_current_, d_x_next_,
                          lp_.num_col_ * sizeof(double),
                          cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_y_current_, d_y_next_,
                          lp_.num_row_ * sizeof(double),
                          cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_ax_current_, d_ax_next_,
                          lp_.num_row_ * sizeof(double),
                          cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_aty_current_, d_aty_next_,
                          lp_.num_col_ * sizeof(double),
                          cudaMemcpyDeviceToDevice));
#else
    x_current_ = x_next_;
    y_current_ = y_next_;
    Ax_cache_ = Ax_next_;
    ATy_cache_ = ATy_next_;
#endif
  }
}
#ifdef CUPDLP_GPU
void PDLPSolver::performHalpernPdhgStepGpu(bool is_major, int k_offset) {
  linalgGpuATy(d_y_current_, d_aty_current_);

  if (is_major) {
    launchKernelHalpernPrimalMajor_wrapper(
        d_x_current_, d_pdhg_primal_, d_x_next_ /*reflected*/,
        d_aty_current_ /*dual_product*/, d_col_cost_, d_col_lower_,
        d_col_upper_, d_primal_step_size_, a_num_cols_,
        d_dSlackPos_ /*dual_slack*/, gpu_stream_);
  } else {
    launchKernelHalpernPrimalMinor_wrapper(
        d_x_current_, d_x_next_ /*reflected*/, d_aty_current_ /*dual_product*/,
        d_col_cost_, d_col_lower_, d_col_upper_, d_primal_step_size_,
        a_num_cols_, gpu_stream_);
  }

  linalgGpuAx(d_x_next_, d_ax_next_);

  if (is_major) {
    launchKernelHalpernDualMajor_wrapper(
        d_y_current_, d_pdhg_dual_, d_y_next_ /*reflected*/,
        d_ax_next_ /*primal_product*/, d_row_lower_, d_is_equality_row_,
        d_dual_step_size_, a_num_rows_, gpu_stream_);
  } else {
    launchKernelHalpernDualMinor_wrapper(
        d_y_current_, d_y_next_ /*reflected*/, d_ax_next_ /*primal_product*/,
        d_row_lower_, d_is_equality_row_, d_dual_step_size_, a_num_rows_,
        gpu_stream_);
  }

  double rho = params_.halpern_gamma;

  launchKernelHalpernBlend_wrapper(d_x_current_, d_x_next_ /*reflected*/,
                                   d_x_anchor_, d_halpern_iteration_, k_offset,
                                   rho, a_num_cols_, gpu_stream_);
  launchKernelHalpernBlend_wrapper(d_y_current_, d_y_next_ /*reflected*/,
                                   d_y_anchor_, d_halpern_iteration_, k_offset,
                                   rho, a_num_rows_, gpu_stream_);
}
#endif

void PDLPSolver::solveReturn(const TerminationStatus term_code) {
#ifdef CUPDLP_GPU
  cleanupGpu();
#endif
  results_.term_code = term_code;
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockSolve);
#endif
}

void PDLPSolver::initialize() {
  // Initialize x and y based on the LP problem
  x_current_.resize(lp_.num_col_, 0.0);
  y_current_.resize(lp_.num_row_, 0.0);
  x_next_.resize(lp_.num_col_, 0.0);
  y_next_.resize(lp_.num_row_, 0.0);

  x_at_last_restart_ = x_current_;
  y_at_last_restart_ = y_current_;
  x_sum_.resize(lp_.num_col_, 0.0);
  y_sum_.resize(lp_.num_row_, 0.0);

  sum_weights_ = 0.0;

  x_avg_.resize(lp_.num_col_, 0.0);
  y_avg_.resize(lp_.num_row_, 0.0);

  reflected_x_.resize(lp_.num_col_, 0.0);
  reflected_y_.resize(lp_.num_row_, 0.0);
  x_anchor_.resize(lp_.num_col_, 0.0);
  y_anchor_.resize(lp_.num_row_, 0.0);
  halpern_iteration_ = 0;

  Ax_cache_.resize(lp_.num_row_, 0.0);
  ATy_cache_.resize(lp_.num_col_, 0.0);
  Ax_next_.resize(lp_.num_row_, 0.0);
  ATy_next_.resize(lp_.num_col_, 0.0);

  // Member variables added back in HPP:
  Ax_avg_.resize(lp_.num_row_, 0.0);
  ATy_avg_.resize(lp_.num_col_, 0.0);
  K_times_x_diff_.resize(lp_.num_row_, 0.0);

  dSlackPos_.resize(lp_.num_col_, 0.0);
  dSlackNeg_.resize(lp_.num_col_, 0.0);
  dSlackPosAvg_.resize(lp_.num_col_, 0.0);
  dSlackNegAvg_.resize(lp_.num_col_, 0.0);
  halpern_dual_slack_next_.resize(lp_.num_col_, 0.0);
  halpern_dual_slack_next_valid_ = false;
}

// Update primal weight
void PDLPSolver::computeStepSizeRatio(PrimalDualParams& working_params) {
  // 1. Calculate the L2 norm of the difference between current and last-restart
  // iterates.
  double primal_diff_norm = linalg::diffTwoNorm(x_at_last_restart_, x_current_);
  double dual_diff_norm = linalg::diffTwoNorm(y_at_last_restart_, y_current_);

  double dMeanStepSize = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);

  // 2. Update the primal weight (beta = omega^2) if movements are significant.
  if (std::min(primal_diff_norm, dual_diff_norm) > 1e-10) {
    double beta_update_ratio = dual_diff_norm / primal_diff_norm;
    double old_beta = stepsize_.beta;

    double dLogBetaUpdate =
        0.5 * std::log(beta_update_ratio) + 0.5 * std::log(std::sqrt(old_beta));
    stepsize_.beta = std::exp(2.0 * dLogBetaUpdate);
  }

  stepsize_.primal_step = dMeanStepSize / std::sqrt(stepsize_.beta);
  stepsize_.dual_step = stepsize_.primal_step * stepsize_.beta;
  working_params.eta = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);
  working_params.omega = std::sqrt(stepsize_.beta);
  restart_scheme_.UpdateBeta(stepsize_.beta);
}

void PDLPSolver::updateAverageIterates(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       HighsInt inner_iter) {
  double dMeanStepSize = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockAverageIterateUpdateX);
#endif
  for (size_t i = 0; i < x.size(); ++i) x_sum_[i] += x[i] * dMeanStepSize;
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockAverageIterateUpdateX);
#endif

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockAverageIterateUpdateY);
#endif
  for (size_t i = 0; i < y.size(); ++i) y_sum_[i] += y[i] * dMeanStepSize;
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockAverageIterateUpdateY);
#endif

  sum_weights_ += dMeanStepSize;
}

void PDLPSolver::computeAverageIterate(std::vector<double>& ax_avg,
                                       std::vector<double>& aty_avg) {
  double dPrimalScale = sum_weights_ > 1e-10 ? 1.0 / sum_weights_ : 1.0;
  double dDualScale = sum_weights_ > 1e-10 ? 1.0 / sum_weights_ : 1.0;

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockAverageIterateComputeX);
#endif
  for (size_t i = 0; i < x_avg_.size(); ++i)
    x_avg_[i] = x_sum_[i] * dPrimalScale;
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockAverageIterateComputeX);
#endif

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockAverageIterateComputeY);
#endif
  for (size_t i = 0; i < y_avg_.size(); ++i) y_avg_[i] = y_sum_[i] * dDualScale;
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockAverageIterateComputeY);
#endif

#if PDLP_DEBUG_LOG
  debug_pdlp_data_.x_average_norm = linalg::vector_norm_squared(x_avg_);
#endif

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockAverageIterateMatrixMultiply);
#endif
  linalg::Ax(lp_, x_avg_, ax_avg);
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockAverageIterateMatrixMultiply);
#endif

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockAverageIterateMatrixTransposeMultiply);
#endif
  linalg::ATy(lp_, y_avg_, aty_avg);
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockAverageIterateMatrixTransposeMultiply);
#endif

#if PDLP_DEBUG_LOG
  debug_pdlp_data_.ax_average_norm = linalg::vector_norm_squared(ax_avg);
  debug_pdlp_data_.aty_average_norm = linalg::vector_norm_squared(aty_avg);
#endif
}

// lambda = c - proj_{\Lambda}(c - K^T y)
std::vector<double> PDLPSolver::computeLambda(
    const std::vector<double>& y, const std::vector<double>& ATy_vector) {
  std::vector<double> lambda(lp_.num_col_, 0.0);
  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    double residual = lp_.col_cost_[i] - ATy_vector[i];

    if (lp_.col_lower_[i] <= -kHighsInf && lp_.col_upper_[i] >= kHighsInf) {
      lambda[i] = 0;
    } else if (lp_.col_lower_[i] <= -kHighsInf) {
      lambda[i] = std::min(0.0, residual);
    } else if (lp_.col_upper_[i] >= kHighsInf) {
      lambda[i] = std::max(0.0, residual);
    } else {
      lambda[i] = residual;
    }
  }
  return lambda;
}

double PDLPSolver::computePrimalFeasibility(
    const std::vector<double>& Ax_vector) {
  std::vector<double> primal_residual(lp_.num_row_, 0.0);

  // Compute Ax - rhs
  for (HighsInt i = 0; i < lp_.num_row_; ++i) {
    primal_residual[i] = Ax_vector[i] - lp_.row_lower_[i];

    // For inequality constraints, project to negative part
    if (!is_equality_row_[i]) {
      primal_residual[i] = std::min(0.0, primal_residual[i]);
    }
  }

  // Apply scaling if needed
  if (scaling_.IsScaled()) {
    const auto& row_scale = scaling_.GetRowScaling();
    for (HighsInt i = 0; i < lp_.num_row_; ++i) {
      primal_residual[i] *= row_scale[i];
    }
  }

  return linalg::vector_norm(primal_residual);
}

void PDLPSolver::computeDualSlacks(const std::vector<double>& dualResidual,
                                   std::vector<double>& dSlackPos,
                                   std::vector<double>& dSlackNeg) {
  // Ensure vectors are correctly sized
  if (dSlackPos.size() != static_cast<size_t>(lp_.num_col_))
    dSlackPos.resize(lp_.num_col_);
  if (dSlackNeg.size() != static_cast<size_t>(lp_.num_col_))
    dSlackNeg.resize(lp_.num_col_);

  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    if (params_.use_halpern_restart) {
      // For Halpern current iterate checks, use major-step dual slack:
      // dual_slack = (proj(temp) - temp) / primal_step.
      // For other cases (e.g. average iterate), fall back to sign/bounds
      // projection.
      const bool use_cached_halpern_slack = (&dSlackPos == &dSlackPos_) &&
                                            (&dSlackNeg == &dSlackNeg_) &&
                                            halpern_dual_slack_next_valid_ &&
                                            (halpern_dual_slack_next_.size() ==
                                             static_cast<size_t>(lp_.num_col_));

      double dual_slack = 0.0;
      if (use_cached_halpern_slack) {
        dual_slack = halpern_dual_slack_next_[i];
      } else {
        const bool has_lower = lp_.col_lower_[i] > -kHighsInf;
        const bool has_upper = lp_.col_upper_[i] < kHighsInf;
        if (has_lower && has_upper) {
          dual_slack = dualResidual[i];
        } else if (has_lower) {
          dual_slack = std::max(0.0, dualResidual[i]);
        } else if (has_upper) {
          dual_slack = std::min(0.0, dualResidual[i]);
        }
      }

      dSlackPos[i] = std::max(0.0, dual_slack);
      dSlackNeg[i] = std::max(0.0, -dual_slack);
    } else {
      // Compute positive slack (for lower bounds)
      // CUPDLP: max(dual_residual, 0) * hasLower
      if (lp_.col_lower_[i] > -kHighsInf) {
        dSlackPos[i] = std::max(0.0, dualResidual[i]);
      } else {
        dSlackPos[i] = 0.0;
      }

      // Compute negative slack (for upper bounds)
      // CUPDLP: -min(dual_residual, 0) * hasUpper
      if (lp_.col_upper_[i] < kHighsInf) {
        dSlackNeg[i] = std::max(0.0, -dualResidual[i]);
      } else {
        dSlackNeg[i] = 0.0;
      }
    }
  }
}

double PDLPSolver::computeDualFeasibility(const std::vector<double>& ATy_vector,
                                          std::vector<double>& dSlackPos,
                                          std::vector<double>& dSlackNeg) {
  std::vector<double> dualResidual(lp_.num_col_, 0.0);
  // dualResidual = c-A'y
  dualResidual = linalg::vector_subtrac(lp_.col_cost_, ATy_vector);

  // Call the refactored function to populate dSlackPos and dSlackNeg
  computeDualSlacks(dualResidual, dSlackPos, dSlackNeg);

  std::vector<double> dual_residual(lp_.num_col_);

  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    // Matching CUPDLP: c - A'y - dSlackPos + dSlackNeg
    dual_residual[i] = dualResidual[i] - dSlackPos[i] + dSlackNeg[i];
  }

  // Apply scaling if needed
  if (scaling_.IsScaled()) {
    const auto& col_scale = scaling_.GetColScaling();
    for (HighsInt i = 0; i < lp_.num_col_; ++i) {
      dual_residual[i] *= col_scale[i];
    }
  }

  double dual_feasibility = linalg::vector_norm(dual_residual);

  return dual_feasibility;
}

std::tuple<double, double, double, double, double>
PDLPSolver::computeDualityGap(const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& lambda) {
  double qTy = 0.0;
  for (HighsInt i = 0; i < lp_.num_row_; ++i) {
    if (lp_.row_lower_[i] > -kHighsInf) {
      qTy += lp_.row_lower_[i] * y[i];
    } else {
      qTy += lp_.row_upper_[i] * y[i];
    }
  }

  double lT_lambda_plus = 0.0;
  double uT_lambda_minus = 0.0;

  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    if (lp_.col_lower_[i] > -kHighsInf) {
      lT_lambda_plus += lp_.col_lower_[i] * std::max(0.0, lambda[i]);
    }
    if (lp_.col_upper_[i] < kHighsInf) {
      uT_lambda_minus += lp_.col_upper_[i] * std::max(0.0, -lambda[i]);
    }
  }

  double dual_objective = qTy + lT_lambda_plus - uT_lambda_minus;

  double cTx = 0.0;
  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    cTx += lp_.col_cost_[i] * x[i];
  }

  double duality_gap = abs(dual_objective - cTx);
  return std::make_tuple(duality_gap, qTy, lT_lambda_plus, uT_lambda_minus,
                         cTx);
}

double PDLPSolver::computeDualObjective(const std::vector<double>& y,
                                        const std::vector<double>& dSlackPos,
                                        const std::vector<double>& dSlackNeg) {
  double dual_obj = lp_.offset_;

  // Compute b'y (or rhs'y in cuPDLP notation)
  for (HighsInt i = 0; i < lp_.num_row_; ++i) {
    dual_obj += lp_.row_lower_[i] * y[i];
  }

  // Add contribution from lower bounds: l'*slackPos
  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    if (lp_.col_lower_[i] > -kHighsInf) {
      dual_obj += lp_.col_lower_[i] * dSlackPos[i];
    }
  }

  // Subtract contribution from upper bounds: u'*slackNeg
  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    if (lp_.col_upper_[i] < kHighsInf) {
      dual_obj -= lp_.col_upper_[i] * dSlackNeg[i];
    }
  }

  return dual_obj;
}

bool PDLPSolver::checkConvergence(
    const HighsInt iter, const std::vector<double>& x,
    const std::vector<double>& y, const std::vector<double>& ax_vector,
    const std::vector<double>& aty_vector, double epsilon,
    SolverResults& results, const char* type,
    // Add slack vectors as non-const references
    std::vector<double>& dSlackPos, std::vector<double>& dSlackNeg) {
  // computeDualSlacks is now called inside computeDualFeasibility

  // Compute primal feasibility
  double primal_feasibility = computePrimalFeasibility(ax_vector);
  results.primal_feasibility = primal_feasibility;

  // Compute dual feasibility
  // This will populate dSlackPos and dSlackNeg
  double dual_feasibility =
      computeDualFeasibility(aty_vector, dSlackPos, dSlackNeg);
  results.dual_feasibility = dual_feasibility;

  // Compute objectives
  double primal_obj = lp_.offset_;
  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    primal_obj += lp_.col_cost_[i] * x[i];
  }
  results.primal_obj = primal_obj;

  // Pass the now-populated slack vectors to computeDualObjective
  double dual_obj = computeDualObjective(y, dSlackPos, dSlackNeg);
  results.dual_obj = dual_obj;

  // Compute duality gap
  double duality_gap = primal_obj - dual_obj;
  results.duality_gap = std::abs(duality_gap);

  // Compute relative gap (matching cuPDLP formula)
  double relative_obj_gap =
      std::abs(duality_gap) / (1.0 + std::abs(primal_obj) + std::abs(dual_obj));
  results.relative_obj_gap = relative_obj_gap;

#if PDLP_DEBUG_LOG
  debugPdlpFeasOptLog(debug_pdlp_log_file_, iter, primal_obj, dual_obj,
                      relative_obj_gap,
                      primal_feasibility / (1.0 + unscaled_rhs_norm_),
                      dual_feasibility / (1.0 + unscaled_c_norm_), type);
#endif

  // Check convergence criteria (matching cuPDLP)
  bool primal_feasible =
      primal_feasibility < epsilon * (1.0 + unscaled_rhs_norm_);
  bool dual_feasible = dual_feasibility < epsilon * (1.0 + unscaled_c_norm_);
  bool gap_small = relative_obj_gap < epsilon;

  return primal_feasible && dual_feasible && gap_small;
}

double PDLPSolver::PowerMethod() {
  const HighsLp& lp = lp_;
  double op_norm_sq = 1.0;
  if (lp.num_col_ == 0 || lp.num_row_ == 0) return op_norm_sq;

  // Parameters for the power method
  const HighsInt max_iter = 20;
  const double tol = 1e-6;

  // Initialize a random vector x
  std::vector<double> x_vec(lp.num_col_);
  std::random_device rd;
  std::mt19937 engine_fixed_seed(12345);  // gen(rd());
  std::uniform_real_distribution<> dis(-1.0, 1.0);
  for (HighsInt i = 0; i < lp.num_col_; ++i) {
    x_vec[i] = dis(engine_fixed_seed);
  }
  linalg::normalize(x_vec);  // Assumes a normalize function in linalg

  const HighsInt kYanyuPowerMethod = 0;
  const HighsInt kYanyuPowerMethodDev = 1;
  const HighsInt kCuPdlpAATPowerMethod = 2;

  const HighsInt power_method =
      // kYanyuPowerMethod;
      // kYanyuPowerMethodDev;
      kCuPdlpAATPowerMethod;
  highsLogDev(params_.log_options_, HighsLogType::kInfo, "Power method: %s\n",
              power_method == kYanyuPowerMethod      ? "Yanyu"
              : power_method == kYanyuPowerMethodDev ? "Yanyu dev"
                                                     : "CuPdlp-C");
  // Dev version of Yanyu power method (based on A'A) has
  //
  // * First iterate as vector of ones
  //
  // * No "convergence" check

  // cuPDLP-C version of power method is based on AA'
  //
  // Allocate memory for matrix-vector products
  std::vector<double> y_vec;
  std::vector<double> z_vec;
  if (power_method < kCuPdlpAATPowerMethod) {
    y_vec.resize(lp.num_row_);
    z_vec.resize(lp.num_col_);
  } else {
    y_vec.resize(lp.num_col_);
    z_vec.resize(lp.num_row_);
  }

  double op_norm_sq_old = 0.0;
  LogLevel log_level = logger_.getLogLevel();
  HighsInt log_iters =
      log_level == LogLevel::kVerbose || log_level == LogLevel::kDebug;

  if (log_iters)
    highsLogUser(params_.log_options_, HighsLogType::kInfo,
                 "It       lambda   dl_lambda\n");

  if (power_method == kYanyuPowerMethodDev) {
    // Start from a vector
    x_vec.assign(lp.num_col_, 1);
  } else if (power_method == kCuPdlpAATPowerMethod) {
    x_vec.assign(lp.num_row_, 1);
  }
  const HighsSparseMatrix& matrix = lp.a_matrix_;
  double lambda = 0.0;
  double previous_lambda = lambda;
  for (HighsInt iter = 0; iter < max_iter; ++iter) {
    if (power_method == kYanyuPowerMethod) {
      // Original Yanyu power method
      //
      // Compute A'Ax = A'(Ax)
      linalg::Ax(lp, x_vec, y_vec);
      linalg::ATy(lp, y_vec, z_vec);  // Note: ATy computes A' * vector

      // Estimate the squared operator norm (largest eigenvalue of A'A)
      op_norm_sq = linalg::dot(
          x_vec, z_vec);  // Assumes a dot product function in linalg

      // Check for convergence
      if (std::abs(op_norm_sq - op_norm_sq_old) < tol * op_norm_sq) {
        return op_norm_sq;
      }
      double dl_op_norm_sq = std::fabs(op_norm_sq - op_norm_sq_old);
      op_norm_sq_old = op_norm_sq;

      // Prepare for the next iteration
      linalg::normalize(z_vec);  // Normalize the result
      x_vec = z_vec;
      if (log_iters)
        highsLogUser(params_.log_options_, HighsLogType::kInfo,
                     "%2d %12.6g %11.4g\n", HighsInt(iter), op_norm_sq,
                     dl_op_norm_sq);
    } else {
      if (power_method == kYanyuPowerMethodDev) {
        // Yanyu power method without "convergence" check
        //
        // Compute z = A'Ax = A'(Ax)
        matrix.product(y_vec, x_vec);
        matrix.productTranspose(z_vec, y_vec);
        // Compute the Rayleigh quotient of x which, since x'x=1, is
        //
        // lambda = x'A'Ax = x'(A'Ax) = x'z
        lambda = linalg::dot(x_vec, z_vec);
        // x = z / norm(z)
        double z_norm = std::sqrt(linalg::dot(z_vec, z_vec));
        for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
          x_vec[iCol] = z_vec[iCol] / z_norm;

      } else {
        // cuPDLP-C power method
        //
        // Compute z = AA'x = A(A'x)
        matrix.productTranspose(y_vec, x_vec);
        matrix.product(z_vec, y_vec);
        // q = z / norm(z)
        double z_norm = std::sqrt(linalg::dot(z_vec, z_vec));
        for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
          z_vec[iRow] /= z_norm;
        // Compute the Rayleigh quotient of q which, since q'q=1, is
        //
        // lambda = q'AA'q = w^Tw = ||w||_2
        //
        // where w = A'q
        //
        // y_vec is no longer needed, so w is stored in it
        matrix.productTranspose(y_vec, z_vec);
        lambda = linalg::dot(y_vec, y_vec);
        x_vec = z_vec;
      }
      double dl_lambda = std::fabs(lambda - previous_lambda);
      previous_lambda = lambda;
      if (log_iters)
        highsLogUser(params_.log_options_, HighsLogType::kInfo,
                     "%2d %12.6g %11.4g\n", HighsInt(iter), lambda, dl_lambda);
    }
  }
  if (power_method != kYanyuPowerMethod) op_norm_sq = lambda;
  // If the method did not converge within max_iter
  return op_norm_sq;
}

void PDLPSolver::setup(const HighsOptions& options, HighsTimer& timer) {
  logger_.setLevel(options.log_dev_level);
  logger_.passHighsLogOptions(options.log_options);
  logger_.print_header();
  highs_timer_p = &timer;
  highsLogUser(options.log_options, HighsLogType::kInfo,
               "Using HiPDLP first order PDLP solver on a %s\n",
#ifdef CUPDLP_GPU
               "GPU"
#else
                "CPU: performance may be disappointing!"
#endif
  );
#ifdef CUPDLP_GPU
  HighsInt n_devices = 0;
  cudaGetDeviceCount(&n_devices);
  if (n_devices != 1)
    highsLogUser(
        options.log_options, HighsLogType::kInfo,
        "Number of CUDA-enabled devices is %d: device 0 will be used\n",
        n_devices);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  highsLogUser(options.log_options, HighsLogType::kInfo, "Cuda device: %s\n",
               prop.name);
  highsLogUser(options.log_options, HighsLogType::kInfo,
               "Global memory available on device (GB): %f\n",
               (float(prop.totalGlobalMem)) / 1e9);
#endif

#if PDLP_PROFILE
  hipdlp_clocks_.timer_pointer_ = &timer;
  hipdlp_timer_.initialiseHipdlpClocks(hipdlp_clocks_);
#endif

#if PDLP_DEBUG_LOG
  debug_pdlp_log_file_ = fopen("HiPDLP.log", "w");
#endif

  params_.initialise();
  //  params.eta = 0; Not set in parse_options_file
  //  params.omega = 0; Not set in parse_options_file
  params_.tolerance = options.pdlp_optimality_tolerance;
  if (options.kkt_tolerance != kDefaultKktTolerance)
    params_.tolerance = options.kkt_tolerance;
  params_.max_iterations = options.pdlp_iteration_limit;
  params_.device_type = Device::CPU;
  params_.time_limit = options.time_limit;

  params_.use_ruiz_scaling = false;
  params_.use_pc_scaling = false;
  params_.use_l2_scaling = false;
  if ((options.pdlp_features_off & kPdlpScalingOff) == 0) {
    // Use scaling: now see which
    params_.use_ruiz_scaling = options.pdlp_scaling_mode & kPdlpScalingRuiz;
    params_.use_pc_scaling = options.pdlp_scaling_mode & kPdlpScalingPC;
    params_.use_l2_scaling = options.pdlp_scaling_mode & kPdlpScalingL2;
  }
  params_.ruiz_iterations = options.pdlp_ruiz_iterations;
  //  params_.ruiz_norm = INFINITY; Not set in parse_options_file
  //  params_.pc_alpha = 1.0; Not set in parse_options_file

  // Restart strategy maps 0/1/2 to RestartStrategy
  params_.restart_strategy = RestartStrategy::NO_RESTART;
  params_.use_halpern_restart = false;
  if ((options.pdlp_features_off & kPdlpRestartOff) == 0) {
    // Use restart: now see which
    if (options.pdlp_restart_strategy == kPdlpRestartStrategyFixed) {
      params_.restart_strategy = RestartStrategy::FIXED_RESTART;
    } else if (options.pdlp_restart_strategy == kPdlpRestartStrategyAdaptive) {
      params_.restart_strategy = RestartStrategy::ADAPTIVE_RESTART;
    } else if (options.pdlp_restart_strategy == kPdlpRestartStrategyHalpern) {
      params_.restart_strategy = RestartStrategy::ADAPTIVE_RESTART;
      params_.use_halpern_restart = true;
    }
  }
  //  params_.fixed_restart_interval = 0; Not set in parse_options_file

  params_.step_size_strategy = StepSizeStrategy::FIXED;
  if ((options.pdlp_features_off & kPdlpAdaptiveStepSizeOff) == 0) {
    // Use adaptive step size: now see which
    if (options.pdlp_step_size_strategy == kPdlpStepSizeStrategyAdaptive) {
      params_.step_size_strategy = StepSizeStrategy::ADAPTIVE;
    } else if (options.pdlp_step_size_strategy ==
               kPdlpStepSizeStrategyMalitskyPock) {
      params_.step_size_strategy = StepSizeStrategy::MALITSKY_POCK;
    } else if (options.pdlp_step_size_strategy == kPdlpStepSizeStrategyPid) {
      params_.step_size_strategy = StepSizeStrategy::PID;
    }
  }
  //  params_.malitsky_pock_params.initialise(); Not set in parse_options_file
  //  params_.adaptive_linesearch_params.initialise(); Not set in
  //  parse_options_file
  scaling_.passParams(&params_);
  restart_scheme_.passParams(&params_);
  // Copy what's needed to use HiGHS logging
  params_.log_options_ = options.log_options;
  // log the options
  logger_.print_params(params_);
}

void PDLPSolver::scaleProblem() {
  scaling_.passLp(&lp_);
  scaling_.passParams(&params_);
  scaling_.Initialize(lp_);
  scaling_.scaleProblem();
}

void PDLPSolver::unscaleSolution(std::vector<double>& x,
                                 std::vector<double>& y) {
  scaling_.unscaleSolution(x, y);

  x_current_ = x;
  y_current_ = y;

  const std::vector<double>& col_scale = scaling_.GetColScaling();
  if (!dSlackPos_.empty() && col_scale.size() == dSlackPos_.size()) {
    for (size_t i = 0; i < dSlackPos_.size(); ++i) {
      dSlackPos_[i] *= col_scale[i];
      dSlackNeg_[i] *= col_scale[i];
    }
  }
}

void PDLPSolver::logSummary() {
  logger_.print_summary(results_, final_iter_count_, highs_timer_p->read());
}

void PrimalDualParams::initialise() {
  this->eta = 0;
  this->omega = 0;
  this->tolerance = 0;
  this->max_iterations = 0;
  this->device_type = Device::CPU;
  this->time_limit = kHighsInf;
  this->restart_strategy = RestartStrategy::NO_RESTART;
  this->fixed_restart_interval = 0;
  this->use_halpern_restart = false;
  this->halpern_gamma = 1.0;
  this->use_ruiz_scaling = false;
  this->use_pc_scaling = false;
  this->use_l2_scaling = false;
  this->ruiz_iterations = 10;
  this->ruiz_norm = INFINITY;
  this->pc_alpha = 1.0;
  this->k_p = 0.99;
  this->k_i = 0.01;
  this->k_d = 0.0;
  this->step_size_strategy = StepSizeStrategy::FIXED;
  this->malitsky_pock_params.initialise();
  this->adaptive_linesearch_params.initialise();
  this->log_options_.clear();
}

void MalitskyPockParams::initialise() {
  this->step_size_interpolation = 0.5;  // Between 0 and 1
  this->step_size_downscaling_factor = 0.7;
  this->linesearch_contraction_factor = 0.99;
}

void AdaptiveLinesearchParams::initialise() {
  this->step_size_reduction_exponent = 0.3;
  this->step_size_growth_exponent = 0.6;
}

// =============================================================================
//  SECTION 4: Step Update Methods (from step.cc)
// =============================================================================

void PDLPSolver::initializeStepSizes() {
  double cost_norm_sq = linalg::vector_norm_squared(lp_.col_cost_);
  double rhs_norm_sq = linalg::vector_norm_squared(lp_.row_lower_);

  if (std::min(cost_norm_sq, rhs_norm_sq) > 1e-6) {
    stepsize_.beta = cost_norm_sq / rhs_norm_sq;
  } else {
    stepsize_.beta = 1.0;
  }

  // Match initial primal weight calculation from cuPDLPx
  params_.omega = (unscaled_c_norm_ + 1.0) / (unscaled_rhs_norm_ + 1.0);
  primal_weight_ = params_.omega;

  // Initialize step sizes based on strategy
  if (params_.step_size_strategy == StepSizeStrategy::FIXED ||
      params_.step_size_strategy == StepSizeStrategy::PID) {
    // Use power method for fixed step size
    // const double op_norm_sq = PowerMethod();
    double op_norm_sq = PowerMethod();

    stepsize_.power_method_lambda = op_norm_sq;

    const double safety_factor = 0.998;
    double base_step = safety_factor / std::sqrt(op_norm_sq);

    params_.eta = base_step;
    stepsize_.primal_step = base_step / params_.omega;
    stepsize_.dual_step = base_step * params_.omega;

    highsLogDev(
        params_.log_options_, HighsLogType::kInfo,
        "Initial step sizes from power method lambda = %g: primal step= "
        "%g; dual step = %g, eta = %g, omega = %g\n",
        op_norm_sq, stepsize_.primal_step, stepsize_.dual_step, params_.eta,
        params_.omega);
  } else {
    // Use matrix infinity norm for adaptive step size
    // Compute infinity norm of matrix elements
    double mat_elem_norm_inf = 0.0;
    const HighsSparseMatrix& matrix = lp_.a_matrix_;
    for (HighsInt i = 0; i < matrix.numNz(); ++i) {
      mat_elem_norm_inf =
          std::max(mat_elem_norm_inf, std::abs(matrix.value_[i]));
    }

    if (mat_elem_norm_inf < 1e-10) {
      mat_elem_norm_inf = 1.0;  // Avoid division by zero
    }

    // Initialize step sizes using infinity norm
    stepsize_.primal_step =
        (1.0 / mat_elem_norm_inf) / std::sqrt(stepsize_.beta);
    stepsize_.dual_step = stepsize_.primal_step * stepsize_.beta;

    highsLogDev(params_.log_options_, HighsLogType::kInfo,
                "Initial step sizes from matrix inf-norm = %g: primal = %g; "
                "dual = %g\n",
                mat_elem_norm_inf, stepsize_.primal_step, stepsize_.dual_step);
  }
}

void PDLPSolver::updatePrimalWeightAtRestart(const SolverResults& results) {
  // Compute delta = pdhg_iterate - anchor
  double primal_dist = 0.0;
  double dual_dist = 0.0;

#ifdef CUPDLP_GPU
  primal_dist = computeDiffNormCuBLAS(d_pdhg_primal_, d_x_anchor_, a_num_cols_);
  dual_dist = computeDiffNormCuBLAS(d_pdhg_dual_, d_y_anchor_, a_num_rows_);
#else
  // CPU version
  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    double diff = x_next_[i] - x_anchor_[i];
    primal_dist += diff * diff;
  }
  for (HighsInt j = 0; j < lp_.num_row_; ++j) {
    double diff = y_next_[j] - y_anchor_[j];
    dual_dist += diff * diff;
  }
  primal_dist = std::sqrt(primal_dist);
  dual_dist = std::sqrt(dual_dist);
#endif

  double rel_primal = results.primal_feasibility / (1.0 + unscaled_rhs_norm_);
  double rel_dual = results.dual_feasibility / (1.0 + unscaled_c_norm_);
  double ratio_infeas = (rel_primal > 0.0) ? (rel_dual / rel_primal) : 1e300;

  const bool silent = true;
  // cuPDLPx-style weight update (PID control)
  if (primal_dist > 1e-16 && dual_dist > 1e-16 && primal_dist < 1e12 &&
      dual_dist < 1e12 && ratio_infeas > 1e-8 && ratio_infeas < 1e8) {
    double error =
        std::log(dual_dist) - std::log(primal_dist) - std::log(primal_weight_);

    primal_weight_error_sum_ =
        params_.i_smooth * primal_weight_error_sum_ + error;
    double delta_error = error - primal_weight_last_error_;

    primal_weight_ *=
        std::exp(params_.k_p * error + params_.k_i * primal_weight_error_sum_ +
                 params_.k_d * delta_error);

    primal_weight_last_error_ = error;
    if (!silent)
      logger_.info("Primal weight updated: " + std::to_string(primal_weight_));
  } else {
    // Revert to best known weight
    primal_weight_ = best_primal_weight_;
    primal_weight_error_sum_ = 0.0;
    primal_weight_last_error_ = 0.0;
    if (!silent)
      logger_.info(
          "Weight update failed (bad norms/ratio), reverted to best: " +
          std::to_string(primal_weight_));
  }

  // Track best weight
  double gap = (rel_primal > 0.0 && rel_dual > 0.0)
                   ? std::abs(std::log10(rel_dual / rel_primal))
                   : best_primal_dual_residual_gap_;
  if (gap < best_primal_dual_residual_gap_) {
    best_primal_dual_residual_gap_ = gap;
    best_primal_weight_ = primal_weight_;
  }

  double eta =
      std::sqrt(stepsize_.primal_step * stepsize_.dual_step);  // invariant
  stepsize_.beta = primal_weight_ * primal_weight_;
  stepsize_.primal_step = eta / primal_weight_;
  stepsize_.dual_step = eta * primal_weight_;
  params_.omega = primal_weight_;
  restart_scheme_.UpdateBeta(stepsize_.beta);
}

std::vector<double> PDLPSolver::updateX(const std::vector<double>& x,
                                        const std::vector<double>& aty,
                                        double primal_step) {
  std::vector<double> x_new(lp_.num_col_);
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    double gradient = lp_.col_cost_[i] - aty[i];
    x_new[i] = linalg::project_box(x[i] - primal_step * gradient,
                                   lp_.col_lower_[i], lp_.col_upper_[i]);
  }
  return x_new;
}

std::vector<double> PDLPSolver::updateY(const std::vector<double>& y,
                                        const std::vector<double>& ax,
                                        const std::vector<double>& ax_next,
                                        double dual_step) {
  std::vector<double> y_new(lp_.num_row_);
  for (HighsInt j = 0; j < lp_.num_row_; j++) {
    double extr_ax = 2 * ax_next[j] - ax[j];
    bool is_equality = (lp_.row_lower_[j] == lp_.row_upper_[j]);
    double q = lp_.row_lower_[j];
    double dual_update = y[j] + dual_step * (q - extr_ax);
    y_new[j] =
        is_equality ? dual_update : linalg::project_non_negative(dual_update);
  }
  return y_new;
}

void PDLPSolver::updateIteratesFixed() {
#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockProjectX);
#endif
  x_next_ = updateX(x_current_, ATy_cache_, stepsize_.primal_step);
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockProjectX);
#endif

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockMatrixMultiply);
#endif
  linalg::Ax(lp_, x_next_, Ax_next_);
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockMatrixMultiply);
#endif

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockProjectY);
#endif
  y_next_ = updateY(y_current_, Ax_cache_, Ax_next_, stepsize_.dual_step);
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockProjectY);
#endif

#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockMatrixTransposeMultiply);
#endif
  linalg::ATy(lp_, y_next_, ATy_next_);
#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockMatrixTransposeMultiply);
#endif

#ifdef CUPDLP_GPU
  // Add this check before the memcpy
  if (d_x_next_ == nullptr) {
    std::cerr << "Error1: d_x_next_ is null!" << std::endl;
    return;
  }
  launchKernelUpdateX(stepsize_.primal_step);
  CUDA_CHECK(cudaDeviceSynchronize());
  linalgGpuAx(d_x_next_, d_ax_next_);
  CUDA_CHECK(cudaDeviceSynchronize());
  launchKernelUpdateY(stepsize_.dual_step);
  CUDA_CHECK(cudaDeviceSynchronize());
  linalgGpuATy(d_y_next_, d_aty_next_);
  CUDA_CHECK(cudaDeviceSynchronize());

  // Add this check before the memcpy
  if (d_x_next_ == nullptr) {
    std::cerr << "Error2: d_x_next_ is null!" << std::endl;
    return;
  }
#endif
}

void PDLPSolver::updateIteratesAdaptive() {
#if PDLP_PROFILE
  hipdlpTimerStart(kHipdlpClockStepSizeAdjustment);
#endif

  const double MIN_ETA = 1e-6;
  const double MAX_ETA = 1.0;

  bool accepted_step = false;
  HighsInt inner_iterations = 0;
  HighsInt num_rejected_steps = 0;

  double dStepSizeUpdate =
      std::sqrt(stepsize_.primal_step * stepsize_.dual_step);

  // Compute candidate solution
  std::vector<double> x_candidate = x_current_;    // Start from current x
  std::vector<double> y_candidate = y_current_;    // Start from current y
  std::vector<double> ax_candidate = Ax_cache_;    // Start from current Ax
  std::vector<double> aty_candidate = ATy_cache_;  // Start from current ATy
  std::vector<double> xupdate = x_next_;
  std::vector<double> yupdate = y_next_;
  std::vector<double> axupdate = Ax_next_;
  std::vector<double> atyupdate = ATy_next_;
  while (!accepted_step) {
    stepsize_.step_size_iter++;  // nStepSizeIter
    inner_iterations++;
    /* cupdlp does not have a max iteration limit here
    if (inner_iterations >= 60) {
      std::cerr << "Warning: Adaptive line search exceeded 60 iterations."
                << std::endl;
      // Force accept the last candidate
      break;
    }
    */
    // Calculate step sizes for this iteration
    double primal_step_update = dStepSizeUpdate / std::sqrt(stepsize_.beta);
    double dual_step_update = dStepSizeUpdate * std::sqrt(stepsize_.beta);

#ifdef CUPDLP_GPU
    launchKernelUpdateX_wrapper(d_x_next_,       // Output (Trial)
                                d_x_current_,    // Input (Base)
                                d_aty_current_,  // Input (Base ATy)
                                d_col_cost_, d_col_lower_, d_col_upper_,
                                primal_step_update, a_num_cols_, gpu_stream_);
    linalgGpuAx(d_x_next_, d_ax_next_);
    launchKernelUpdateY_wrapper(d_y_next_,      // Output (Trial)
                                d_y_current_,   // Input (Base)
                                d_ax_current_,  // Input (Base Ax)
                                d_ax_next_,     // Input (Trial Ax)
                                d_row_lower_, d_is_equality_row_,
                                dual_step_update, a_num_rows_, gpu_stream_);
    linalgGpuATy(d_y_next_, d_aty_next_);
    double movement =
        computeMovementGpu(d_x_next_, d_x_current_, d_y_next_, d_y_current_);
    double nonlinearity = computeNonlinearityGpu(d_x_next_, d_x_current_,
                                                 d_aty_next_, d_aty_current_);
#else
    // Primal update
#if PDLP_PROFILE
    hipdlpTimerStart(kHipdlpClockProjectX);
#endif
    xupdate = updateX(x_candidate, aty_candidate, primal_step_update);
#if PDLP_PROFILE
    hipdlpTimerStop(kHipdlpClockProjectX);
#endif

#if PDLP_PROFILE
    hipdlpTimerStart(kHipdlpClockMatrixMultiply);
#endif
    linalg::Ax(lp_, xupdate, axupdate);
#if PDLP_PROFILE
    hipdlpTimerStop(kHipdlpClockMatrixMultiply);
#endif

    // Dual update with timing
#if PDLP_PROFILE
    hipdlpTimerStart(kHipdlpClockProjectY);
#endif
    yupdate = updateY(y_candidate, ax_candidate, axupdate, dual_step_update);
#if PDLP_PROFILE
    hipdlpTimerStop(kHipdlpClockProjectY);
#endif

#if PDLP_PROFILE
    hipdlpTimerStart(kHipdlpClockMatrixTransposeMultiply);
#endif
    linalg::ATy(lp_, yupdate, atyupdate);
#if PDLP_PROFILE
    hipdlpTimerStop(kHipdlpClockMatrixTransposeMultiply);
#endif

    // Compute deltas
    std::vector<double> delta_x(lp_.num_col_);
    std::vector<double> delta_y(lp_.num_row_);
    std::vector<double> delta_aty(lp_.num_col_);

    for (size_t i = 0; i < x_candidate.size(); ++i) {
      delta_x[i] = xupdate[i] - x_candidate[i];
    }
    for (size_t i = 0; i < y_candidate.size(); ++i) {
      delta_y[i] = yupdate[i] - y_candidate[i];
    }
    for (size_t i = 0; i < aty_candidate.size(); ++i) {
      delta_aty[i] = atyupdate[i] - aty_candidate[i];
    }

    // Check numerical stability
    /*
    if (!CheckNumericalStability(delta_x, delta_y)) {
      std::cerr << "Numerical instability detected" << std::endl;
      dStepSizeUpdate *= 0.5;  // Drastically reduce step size
      continue;
    }
    */

    // Compute movement and nonlinearity
    double movement = computeMovement(delta_x, delta_y);
    double nonlinearity = computeNonlinearity(delta_x, delta_aty);
#endif
    // Compute step size limit
    double step_size_limit = (nonlinearity != 0.0)
                                 ? (movement / std::fabs(nonlinearity))
                                 : std::numeric_limits<double>::infinity();

    /*
    highsLogDev(*log_options_, HighsLogType::kInfo,
                 "Iteration %d: eta = %g, movement = %g nonlinearity = %g, limit
    = %g\n", HighsInt(inner_iterations), current_eta, movement, nonlinearity.
    step_size_limit);
    */
    if (dStepSizeUpdate <= step_size_limit) {
      accepted_step = true;
    } else {
      num_rejected_steps++;
    }

    // Compute new step size
    double first_term =
        (std::isinf(step_size_limit))
            ? step_size_limit
            : (1.0 - std::pow(stepsize_.step_size_iter + 1.0,
                              -params_.adaptive_linesearch_params
                                   .step_size_reduction_exponent)) *
                  step_size_limit;

    double second_term =
        (1.0 +
         std::pow(
             stepsize_.step_size_iter + 1.0,
             -params_.adaptive_linesearch_params.step_size_growth_exponent)) *
        dStepSizeUpdate;

    dStepSizeUpdate = std::min(first_term, second_term);
    // current_eta_ = std::max(MIN_ETA, std::min(MAX_ETA, current_eta_));
  }

  x_next_ = xupdate;
  y_next_ = yupdate;
  Ax_next_ = axupdate;
  ATy_next_ = atyupdate;

  current_eta_ = dStepSizeUpdate;
  stepsize_.primal_step = dStepSizeUpdate / std::sqrt(stepsize_.beta);
  stepsize_.dual_step = dStepSizeUpdate * std::sqrt(stepsize_.beta);

#if PDLP_PROFILE
  hipdlpTimerStop(kHipdlpClockStepSizeAdjustment);
#endif
}

bool PDLPSolver::updateIteratesMalitskyPock(bool first_iteration) {
  std::vector<double> x_candidate(lp_.num_col_);
  return true;
}

// =============================================================================
//  SECTION 5: Step Helper Methods (from step.cc)
// =============================================================================

bool PDLPSolver::CheckNumericalStability(const std::vector<double>& delta_x,
                                         const std::vector<double>& delta_y) {
  double movement = computeMovement(delta_x, delta_y);

  if (movement == 0.0) {
    highsLogUser(log_options_, HighsLogType::kInfo,
                 "Warning: Zero movement detected - numerical termination\n");
    return false;
  }

  if (movement > kDivergentMovement) {
    highsLogUser(log_options_, HighsLogType::kInfo,
                 "Warning: Divergent movement detected: %g\n", movement);
    return false;
  }
  return true;  // Placeholder
}

double PDLPSolver::computeMovement(const std::vector<double>& delta_primal,
                                   const std::vector<double>& delta_dual) {
  double primal_squared_norm = 0.0;
  for (const auto& val : delta_primal) {
    primal_squared_norm += val * val;
  }

  double dual_squared_norm = 0.0;
  for (const auto& val : delta_dual) {
    dual_squared_norm += val * val;
  }

  double primal_weight = std::sqrt(stepsize_.beta);
  return (0.5 * primal_weight * primal_squared_norm) +
         (0.5 / primal_weight) * dual_squared_norm;
}

double PDLPSolver::computeNonlinearity(const std::vector<double>& delta_primal,
                                       const std::vector<double>& delta_aty) {
  // Nonlinearity = |Dx' * D(A'y)|
  double nonlinearity = 0.0;
  for (size_t i = 0; i < delta_primal.size(); ++i) {
    nonlinearity += delta_primal[i] * delta_aty[i];
  }
  return nonlinearity;  // cupdlp does not take absolute value
}

#if PDLP_PROFILE
void PDLPSolver::reportHipdlpTimer() {
  hipdlp_timer_.reportHipdlpCoreClock(hipdlp_clocks_);
  hipdlp_timer_.reportHipdlpSolveClock(hipdlp_clocks_);
  hipdlp_timer_.reportHipdlpIterateUpdateClock(hipdlp_clocks_);
  hipdlp_timer_.reportHipdlpAverageIterateClock(hipdlp_clocks_);
  hipdlp_timer_.reportHipdlpMatrixMultiplyClock(hipdlp_clocks_);
}
#endif

#if PDLP_PROFILE
void PDLPSolver::hipdlpTimerStart(const HighsInt hipdlp_clock) {
  HighsInt highs_timer_clock = hipdlp_clocks_.clock_[hipdlp_clock];
  hipdlp_clocks_.timer_pointer_->start(highs_timer_clock);
}
#endif

#if PDLP_PROFILE
void PDLPSolver::hipdlpTimerStop(const HighsInt hipdlp_clock) {
  HighsInt highs_timer_clock = hipdlp_clocks_.clock_[hipdlp_clock];
  hipdlp_clocks_.timer_pointer_->stop(highs_timer_clock);
}
#endif

#if PDLP_DEBUG_LOG
void PDLPSolver::closeDebugLog() {
  if (debug_pdlp_log_file_) fclose(debug_pdlp_log_file_);
}
#endif

// =============================================================================
//  SECTION 5: GPU Part
// =============================================================================
#ifdef CUPDLP_GPU
void PDLPSolver::setupGpu() {
  CUDA_CHECK(cudaStreamCreate(&gpu_stream_));

  // 1. Initialize cuSPARSE
  CUSPARSE_CHECK(cusparseCreate(&cusparse_handle_));
  CUBLAS_CHECK(cublasCreate(&cublas_handle_));
  CUSPARSE_CHECK(cusparseSetStream(cusparse_handle_, gpu_stream_));
  CUBLAS_CHECK(cublasSetStream(cublas_handle_, gpu_stream_));
  // 2. Get matrix data from lp_ (CSC)
  a_num_rows_ = lp_.num_row_;
  a_num_cols_ = lp_.num_col_;
  a_nnz_ = lp_.a_matrix_.numNz();

  HighsSparseMatrix lp_csr = lp_.a_matrix_;
  lp_csr.ensureRowwise();

  const std::vector<HighsInt>& h_a_row_ptr = lp_csr.start_;
  const std::vector<HighsInt>& h_a_col_ind = lp_csr.index_;
  const std::vector<double>& h_a_val = lp_csr.value_;

  // 3. Allocate and copy A's CSR data to GPU
  CUDA_CHECK(
      cudaMalloc((void**)&d_a_row_ptr_, (a_num_rows_ + 1) * sizeof(HighsInt)));
  CUDA_CHECK(cudaMalloc((void**)&d_a_col_ind_, a_nnz_ * sizeof(HighsInt)));
  CUDA_CHECK(cudaMalloc((void**)&d_a_val_, a_nnz_ * sizeof(double)));
  CUDA_CHECK(cudaMemcpy(d_a_row_ptr_, h_a_row_ptr.data(),
                        (a_num_rows_ + 1) * sizeof(HighsInt),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_a_col_ind_, h_a_col_ind.data(),
                        a_nnz_ * sizeof(HighsInt), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_a_val_, h_a_val.data(), a_nnz_ * sizeof(double),
                        cudaMemcpyHostToDevice));

  CUSPARSE_CHECK(cusparseCreateCsr(&mat_a_csr_, a_num_rows_, a_num_cols_,
                                   a_nnz_, d_a_row_ptr_, d_a_col_ind_, d_a_val_,
                                   CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                   CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

  // 4. Create matrix AT in CSR format = A in CSC
  const std::vector<HighsInt>& h_at_row_ptr = lp_.a_matrix_.start_;
  const std::vector<HighsInt>& h_at_col_ind = lp_.a_matrix_.index_;
  const std::vector<double>& h_at_val = lp_.a_matrix_.value_;

  CUDA_CHECK(cudaMalloc((void**)&d_at_row_ptr_,
                        (a_num_cols_ + 1) * sizeof(HighsInt)));  // Fixed!
  CUDA_CHECK(cudaMalloc((void**)&d_at_col_ind_, a_nnz_ * sizeof(HighsInt)));
  CUDA_CHECK(cudaMalloc((void**)&d_at_val_, a_nnz_ * sizeof(double)));

  CUDA_CHECK(cudaMemcpy(d_at_row_ptr_, h_at_row_ptr.data(),
                        (a_num_cols_ + 1) * sizeof(HighsInt),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_at_col_ind_, h_at_col_ind.data(),
                        a_nnz_ * sizeof(HighsInt), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_at_val_, h_at_val.data(), a_nnz_ * sizeof(double),
                        cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMalloc(&d_halpern_iteration_, sizeof(int)));
  CUDA_CHECK(cudaMemset(d_halpern_iteration_, 0, sizeof(int)));
  CUDA_CHECK(cudaMalloc(&d_primal_step_size_, sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_dual_step_size_, sizeof(double)));

  // Create AT descriptor with swapped dimensions
  CUSPARSE_CHECK(cusparseCreateCsr(
      &mat_a_T_csr_, a_num_cols_, a_num_rows_, a_nnz_,  // Dimensions swapped!
      d_at_row_ptr_, d_at_col_ind_, d_at_val_, CUSPARSE_INDEX_32I,
      CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

  CUDA_CHECK(cudaMalloc(&d_col_cost_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_col_lower_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_col_upper_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_row_lower_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_is_equality_row_, a_num_rows_ * sizeof(bool)));
  CUDA_CHECK(cudaMalloc(&d_x_current_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_y_current_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x_avg_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_y_avg_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x_next_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_y_next_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x_at_last_restart_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_y_at_last_restart_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x_anchor_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_y_anchor_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(
      cudaMalloc(&d_x_temp_diff_norm_result_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(
      cudaMalloc(&d_y_temp_diff_norm_result_, a_num_rows_ * sizeof(double)));
  if (params_.use_halpern_restart) {
    CUDA_CHECK(cudaMalloc(&d_pdhg_primal_, a_num_cols_ * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_pdhg_dual_, a_num_rows_ * sizeof(double)));
    CUDA_CHECK(cudaMemset(d_pdhg_primal_, 0, a_num_cols_ * sizeof(double)));
    CUDA_CHECK(cudaMemset(d_pdhg_dual_, 0, a_num_rows_ * sizeof(double)));
  }
  CUDA_CHECK(cudaMalloc(&d_delta_x_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_delta_y_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_AT_delta_y_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_ax_current_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_aty_current_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_ax_next_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_aty_next_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_ax_avg_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_aty_avg_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_x_sum_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_y_sum_, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_convergence_results_, 4 * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_dSlackPos_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_dSlackNeg_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_dSlackPosAvg_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_dSlackNegAvg_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_col_scale_, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_row_scale_, a_num_rows_ * sizeof(double)));

  CUSPARSE_CHECK(
      cusparseCreateDnVec(&vec_x_desc_, a_num_cols_, d_x_current_, CUDA_R_64F));
  CUSPARSE_CHECK(
      cusparseCreateDnVec(&vec_y_desc_, a_num_rows_, d_y_current_, CUDA_R_64F));
  CUSPARSE_CHECK(cusparseCreateDnVec(&vec_ax_desc_, a_num_rows_, d_ax_current_,
                                     CUDA_R_64F));
  CUSPARSE_CHECK(cusparseCreateDnVec(&vec_aty_desc_, a_num_cols_,
                                     d_aty_current_, CUDA_R_64F));

  CUDA_CHECK(cudaMemcpy(d_col_cost_, lp_.col_cost_.data(),
                        a_num_cols_ * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_col_lower_, lp_.col_lower_.data(),
                        a_num_cols_ * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_col_upper_, lp_.col_upper_.data(),
                        a_num_cols_ * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_row_lower_, lp_.row_lower_.data(),
                        a_num_rows_ * sizeof(double), cudaMemcpyHostToDevice));
  std::vector<uint8_t> temp_equality(a_num_rows_);
  for (HighsInt i = 0; i < a_num_rows_; ++i) {
    temp_equality[i] = is_equality_row_[i] ? 1 : 0;
  }

  // Copy to device
  CUDA_CHECK(cudaMemcpy(d_is_equality_row_, temp_equality.data(),
                        a_num_rows_ * sizeof(uint8_t), cudaMemcpyHostToDevice));

  // 6. Preallocate buffer for cuSPARSE SpMV
  // Buffer for Ax
  double alpha = 1.0;
  double beta = 0.0;
  cusparseDnVecDescr_t vec_x, vec_ax;
  CUSPARSE_CHECK(
      cusparseCreateDnVec(&vec_x, a_num_cols_, d_x_current_, CUDA_R_64F));
  CUSPARSE_CHECK(
      cusparseCreateDnVec(&vec_ax, a_num_rows_, d_ax_current_, CUDA_R_64F));

  CUSPARSE_CHECK(cusparseSpMV_bufferSize(
      cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_a_csr_,
      vec_x, &beta, vec_ax, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2,
      &spmv_buffer_size_ax_));
  CUDA_CHECK(cudaMalloc(&d_spmv_buffer_ax_, spmv_buffer_size_ax_));

  CUSPARSE_CHECK(cusparseSpMV_preprocess(
      cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_a_csr_,
      vec_x_desc_, &beta, vec_ax_desc_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2,
      d_spmv_buffer_ax_));

  CUSPARSE_CHECK(cusparseDestroyDnVec(vec_x));
  CUSPARSE_CHECK(cusparseDestroyDnVec(vec_ax));

  // Buffer for ATy
  cusparseDnVecDescr_t vec_y, vec_aty;
  CUSPARSE_CHECK(
      cusparseCreateDnVec(&vec_y, a_num_rows_, d_y_current_, CUDA_R_64F));
  CUSPARSE_CHECK(
      cusparseCreateDnVec(&vec_aty, a_num_cols_, d_aty_current_, CUDA_R_64F));
  CUSPARSE_CHECK(cusparseSpMV_bufferSize(
      cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_a_T_csr_,
      vec_y, &beta, vec_aty, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2,
      &spmv_buffer_size_aty_));
  CUDA_CHECK(cudaMalloc(&d_spmv_buffer_aty_, spmv_buffer_size_aty_));

  CUSPARSE_CHECK(cusparseSpMV_preprocess(
      cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_a_T_csr_,
      vec_y_desc_, &beta, vec_aty_desc_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2,
      d_spmv_buffer_aty_));

  CUSPARSE_CHECK(cusparseDestroyDnVec(vec_y));
  CUSPARSE_CHECK(cusparseDestroyDnVec(vec_aty));

  CUDA_CHECK(cudaMemset(d_x_current_, 0, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_y_current_, 0, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_x_avg_, 0, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_y_avg_, 0, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_x_next_, 0, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_y_next_, 0, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_ax_current_, 0, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_ax_next_, 0, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_ax_avg_, 0, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_aty_avg_, 0, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_x_sum_, 0, a_num_cols_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_y_sum_, 0, a_num_rows_ * sizeof(double)));
  CUDA_CHECK(cudaMemset(d_aty_current_, 0, a_num_cols_ * sizeof(double)));
  sum_weights_gpu_ = 0.0;

  if (scaling_.IsScaled()) {
    CUDA_CHECK(cudaMemcpy(d_col_scale_, scaling_.GetColScaling().data(),
                          a_num_cols_ * sizeof(double),
                          cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_row_scale_, scaling_.GetRowScaling().data(),
                          a_num_rows_ * sizeof(double),
                          cudaMemcpyHostToDevice));
  } else {
    cudaFree(d_col_scale_);
    d_col_scale_ = nullptr;
    cudaFree(d_row_scale_);
    d_row_scale_ = nullptr;
  }

  size_t max_size = std::max(a_num_cols_, a_num_rows_);
  CUDA_CHECK(cudaMalloc(&d_buffer_, max_size * sizeof(double)));
  CUDA_CHECK(cudaMalloc(&d_buffer2_, max_size * sizeof(double)));

  highsLogDev(params_.log_options_, HighsLogType::kInfo,
              "GPU setup complete. Matrix A (CSR) and A^T (CSR) transferred "
              "to device.\n");
}

void PDLPSolver::cleanupGpu() {
  if (gpu_stream_) CUDA_CHECK(cudaStreamDestroy(gpu_stream_));
  if (cusparse_handle_) CUSPARSE_CHECK(cusparseDestroy(cusparse_handle_));
  if (cublas_handle_) CUBLAS_CHECK(cublasDestroy(cublas_handle_));
  if (mat_a_csr_) CUSPARSE_CHECK(cusparseDestroySpMat(mat_a_csr_));
  if (mat_a_T_csr_) CUSPARSE_CHECK(cusparseDestroySpMat(mat_a_T_csr_));
  CUDA_CHECK(cudaFree(d_a_row_ptr_));
  CUDA_CHECK(cudaFree(d_a_col_ind_));
  CUDA_CHECK(cudaFree(d_a_val_));
  CUDA_CHECK(cudaFree(d_at_row_ptr_));
  CUDA_CHECK(cudaFree(d_at_col_ind_));
  CUDA_CHECK(cudaFree(d_at_val_));
  CUDA_CHECK(cudaFree(d_col_cost_));
  CUDA_CHECK(cudaFree(d_col_lower_));
  CUDA_CHECK(cudaFree(d_col_upper_));
  CUDA_CHECK(cudaFree(d_row_lower_));
  CUDA_CHECK(cudaFree(d_is_equality_row_));
  CUDA_CHECK(cudaFree(d_x_at_last_restart_));
  CUDA_CHECK(cudaFree(d_y_at_last_restart_));
  CUDA_CHECK(cudaFree(d_halpern_iteration_));
  CUDA_CHECK(cudaFree(d_primal_step_size_));
  CUDA_CHECK(cudaFree(d_dual_step_size_));
  if (d_x_anchor_) CUDA_CHECK(cudaFree(d_x_anchor_));
  if (d_y_anchor_) CUDA_CHECK(cudaFree(d_y_anchor_));
  if (d_pdhg_primal_) CUDA_CHECK(cudaFree(d_pdhg_primal_));
  if (d_pdhg_dual_) CUDA_CHECK(cudaFree(d_pdhg_dual_));
  CUDA_CHECK(cudaFree(d_delta_x_));
  CUDA_CHECK(cudaFree(d_delta_y_));
  CUDA_CHECK(cudaFree(d_AT_delta_y_));
  CUDA_CHECK(cudaFree(d_x_temp_diff_norm_result_));
  CUDA_CHECK(cudaFree(d_y_temp_diff_norm_result_));
  CUDA_CHECK(cudaFree(d_x_current_));
  CUDA_CHECK(cudaFree(d_y_current_));
  CUDA_CHECK(cudaFree(d_x_next_));
  CUDA_CHECK(cudaFree(d_y_next_));
  CUDA_CHECK(cudaFree(d_ax_current_));
  CUDA_CHECK(cudaFree(d_aty_current_));
  CUDA_CHECK(cudaFree(d_ax_avg_));
  CUDA_CHECK(cudaFree(d_aty_avg_));
  CUDA_CHECK(cudaFree(d_ax_next_));
  CUDA_CHECK(cudaFree(d_aty_next_));
  CUDA_CHECK(cudaFree(d_x_sum_));
  CUDA_CHECK(cudaFree(d_y_sum_));
  CUDA_CHECK(cudaFree(d_spmv_buffer_ax_));
  CUDA_CHECK(cudaFree(d_spmv_buffer_aty_));
  CUDA_CHECK(cudaFree(d_convergence_results_));
  CUDA_CHECK(cudaFree(d_dSlackPos_));
  CUDA_CHECK(cudaFree(d_dSlackNeg_));
  CUDA_CHECK(cudaFree(d_dSlackPosAvg_));
  CUDA_CHECK(cudaFree(d_dSlackNegAvg_));
  if (d_col_scale_) CUDA_CHECK(cudaFree(d_col_scale_));
  if (d_row_scale_) CUDA_CHECK(cudaFree(d_row_scale_));
  CUDA_CHECK(cudaFree(d_buffer_));
  CUDA_CHECK(cudaFree(d_buffer2_));
  if (vec_x_desc_) CUSPARSE_CHECK(cusparseDestroyDnVec(vec_x_desc_));
  if (vec_y_desc_) CUSPARSE_CHECK(cusparseDestroyDnVec(vec_y_desc_));
  if (vec_ax_desc_) CUSPARSE_CHECK(cusparseDestroyDnVec(vec_ax_desc_));
  if (vec_aty_desc_) CUSPARSE_CHECK(cusparseDestroyDnVec(vec_aty_desc_));
}

void PDLPSolver::linalgGpuAx(const double* d_x_in, double* d_ax_out) {
  // Ax = 1.0 * A * x + 0.0 * ax
  double alpha = 1.0;
  double beta = 0.0;

  CUSPARSE_CHECK(cusparseDnVecSetValues(vec_x_desc_, (void*)d_x_in));
  CUSPARSE_CHECK(cusparseDnVecSetValues(vec_ax_desc_, (void*)d_ax_out));
  if (spmv_buffer_size_ax_ == 0) {
    CUSPARSE_CHECK(cusparseSpMV_bufferSize(
        cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, mat_a_csr_,
        vec_x_desc_, &beta, vec_ax_desc_, CUDA_R_64F, CUSPARSE_SPMV_CSR_ALG2,
        &spmv_buffer_size_ax_));
    CUDA_CHECK(cudaMalloc(&d_spmv_buffer_ax_, spmv_buffer_size_ax_));
  }

  CUSPARSE_CHECK(
      cusparseSpMV(cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                   mat_a_csr_, vec_x_desc_, &beta, vec_ax_desc_, CUDA_R_64F,
                   CUSPARSE_SPMV_CSR_ALG2, d_spmv_buffer_ax_));
}

void PDLPSolver::linalgGpuATy(const double* d_y_in, double* d_aty_out) {
  // ATy = 1.0 * A^T * y + 0.0 * aty
  double alpha = 1.0;
  double beta = 0.0;

  CUSPARSE_CHECK(cusparseDnVecSetValues(vec_y_desc_, (void*)d_y_in));
  CUSPARSE_CHECK(cusparseDnVecSetValues(vec_aty_desc_, (void*)d_aty_out));
  if (spmv_buffer_size_aty_ == 0) {
    CUSPARSE_CHECK(cusparseSpMV_bufferSize(
        cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
        mat_a_T_csr_, vec_y_desc_, &beta, vec_aty_desc_, CUDA_R_64F,
        CUSPARSE_SPMV_CSR_ALG2, &spmv_buffer_size_aty_));
    CUDA_CHECK(cudaMalloc(&d_spmv_buffer_aty_, spmv_buffer_size_aty_));
  }
  CUSPARSE_CHECK(
      cusparseSpMV(cusparse_handle_, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha,
                   mat_a_T_csr_, vec_y_desc_, &beta, vec_aty_desc_, CUDA_R_64F,
                   CUSPARSE_SPMV_CSR_ALG2, d_spmv_buffer_aty_));
}

void PDLPSolver::launchKernelUpdateX(double primal_step) {
  launchKernelUpdateX_wrapper(d_x_next_, d_x_current_, d_aty_current_,
                              d_col_cost_, d_col_lower_, d_col_upper_,
                              primal_step, a_num_cols_, gpu_stream_);
  CUDA_CHECK(cudaGetLastError());
}

void PDLPSolver::launchKernelUpdateY(double dual_step) {
  launchKernelUpdateY_wrapper(d_y_next_, d_y_current_, d_ax_current_,
                              d_ax_next_, d_row_lower_, d_is_equality_row_,
                              dual_step, a_num_rows_, gpu_stream_);
  CUDA_CHECK(cudaGetLastError());
}

void PDLPSolver::launchKernelUpdateAverages(double weight) {
  launchKernelUpdateAverages_wrapper(d_x_sum_, d_y_sum_, d_x_next_, d_y_next_,
                                     weight, a_num_cols_, a_num_rows_,
                                     gpu_stream_);
  CUDA_CHECK(cudaGetLastError());
}

void PDLPSolver::launchKernelScaleVector(double* d_out, const double* d_in,
                                         double scale, int n) {
  launchKernelScaleVector_wrapper(d_out, d_in, scale, n, gpu_stream_);
  CUDA_CHECK(cudaGetLastError());
}

bool PDLPSolver::checkConvergenceGpu(const HighsInt iter, const double* d_x,
                                     const double* d_y, const double* d_ax,
                                     const double* d_aty, double epsilon,
                                     SolverResults& results, const char* type,
                                     double* d_slackPos_out,
                                     double* d_slackNeg_out,
                                     bool use_halpern_slack) {
  launchCheckConvergenceKernels_wrapper(
      d_convergence_results_, d_slackPos_out, d_slackNeg_out, d_x, d_y, d_ax,
      d_aty, d_col_cost_, d_row_lower_, d_col_lower_, d_col_upper_,
      d_is_equality_row_, d_col_scale_, d_row_scale_, lp_.num_col_,
      lp_.num_row_, use_halpern_slack, gpu_stream_);

  // copy 4 doubles back to CPU

  double h_results[4];
  CUDA_CHECK(cudaMemcpy(h_results, d_convergence_results_, 4 * sizeof(double),
                        cudaMemcpyDeviceToHost));

  double primal_feas_sq = h_results[0];
  double dual_feas_sq = h_results[1];
  double primal_obj = h_results[2] + lp_.offset_;
  double dual_obj = h_results[3] + lp_.offset_;

  results.primal_feasibility = std::sqrt(primal_feas_sq);
  results.dual_feasibility = std::sqrt(dual_feas_sq);
  results.primal_obj = primal_obj;
  results.dual_obj = dual_obj;

  double duality_gap = primal_obj - dual_obj;
  results.duality_gap = std::abs(duality_gap);
  results.relative_obj_gap =
      std::abs(duality_gap) / (1.0 + std::abs(primal_obj) + std::abs(dual_obj));

#if PDLP_DEBUG_LOG
  debugPdlpFeasOptLog(debug_pdlp_log_file_, iter, primal_obj, dual_obj,
                      results.relative_obj_gap,
                      results.primal_feasibility / (1.0 + unscaled_rhs_norm_),
                      results.dual_feasibility / (1.0 + unscaled_c_norm_),
                      type);
#endif

  bool primal_feasible =
      results.primal_feasibility < epsilon * (1.0 + unscaled_rhs_norm_);
  bool dual_feasible =
      results.dual_feasibility < epsilon * (1.0 + unscaled_c_norm_);
  bool gap_small = results.relative_obj_gap < epsilon;

  return primal_feasible && dual_feasible && gap_small;
}

void PDLPSolver::computeStepSizeRatioGpu(PrimalDualParams& working_params) {
  // 1. Compute ||x_last - x_current||^2 using cuBLAS
  double primal_diff_norm =
      computeDiffNormCuBLAS(d_x_at_last_restart_, d_x_current_, a_num_cols_);

  // 2. Compute ||y_last - y_current||^2 using cuBLAS
  double dual_diff_norm =
      computeDiffNormCuBLAS(d_y_at_last_restart_, d_y_current_, a_num_rows_);

  double dMeanStepSize = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);

  // 3. Update beta (same CPU logic)
  if (std::min(primal_diff_norm, dual_diff_norm) > 1e-10) {
    double beta_update_ratio = dual_diff_norm / primal_diff_norm;
    double old_beta = stepsize_.beta;
    double dLogBetaUpdate =
        0.5 * std::log(beta_update_ratio) + 0.5 * std::log(std::sqrt(old_beta));
    stepsize_.beta = std::exp(2.0 * dLogBetaUpdate);
  }

  // Update step sizes
  stepsize_.primal_step = dMeanStepSize / std::sqrt(stepsize_.beta);
  stepsize_.dual_step = stepsize_.primal_step * stepsize_.beta;
  working_params.eta = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);
  working_params.omega = std::sqrt(stepsize_.beta);
  restart_scheme_.UpdateBeta(stepsize_.beta);
}

void PDLPSolver::updateAverageIteratesGpu(HighsInt inner_iter) {
  double dMeanStepSize = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);

  launchKernelUpdateAverages(dMeanStepSize);

  sum_weights_gpu_ += dMeanStepSize;
}

void PDLPSolver::computeAverageIterateGpu() {
  double dScale = sum_weights_gpu_ > 1e-10 ? 1.0 / sum_weights_gpu_ : 1.0;

  // x_avg = x_sum * scale
  launchKernelScaleVector(d_x_avg_, d_x_sum_, dScale, a_num_cols_);

  // y_avg = y_sum * scale
  launchKernelScaleVector(d_y_avg_, d_y_sum_, dScale, a_num_rows_);

  // Recompute Ax_avg and ATy_avg on GPU
  linalgGpuAx(d_x_avg_, d_ax_avg_);
  linalgGpuATy(d_y_avg_, d_aty_avg_);

#if PDLP_DEBUG_LOG
  // copy x_avg to host
  CUDA_CHECK(cudaMemcpy(x_avg_.data(), d_x_avg_, a_num_cols_ * sizeof(double),
                        cudaMemcpyDeviceToHost));
  debug_pdlp_data_.x_average_norm = linalg::vector_norm_squared(x_avg_);
#endif
}

double PDLPSolver::computeMovementGpu(const double* d_x_new,
                                      const double* d_x_old,
                                      const double* d_y_new,
                                      const double* d_y_old) {
  // 1. Compute ||x_new - x_old|| using cuBLAS
  double primal_diff_norm =
      computeDiffNormCuBLAS(d_x_new, d_x_old, a_num_cols_);

  // 2. Compute ||y_new - y_old|| using cuBLAS
  double dual_diff_norm = computeDiffNormCuBLAS(d_y_new, d_y_old, a_num_rows_);

  // 3. Combine on CPU
  double primal_weight = std::sqrt(stepsize_.beta);
  double primal_diff_sq = primal_diff_norm * primal_diff_norm;
  double dual_diff_sq = dual_diff_norm * dual_diff_norm;

  return (0.5 * primal_weight * primal_diff_sq) +
         (0.5 / primal_weight * dual_diff_sq);
}

double PDLPSolver::computeNonlinearityGpu(const double* d_x_new,
                                          const double* d_x_old,
                                          const double* d_aty_new,
                                          const double* d_aty_old) {
  // 1. Compute delta_x = x_new - x_old
  CUDA_CHECK(cudaMemcpy(d_buffer_, d_x_new, a_num_cols_ * sizeof(double),
                        cudaMemcpyDeviceToDevice));
  double alpha = -1.0;
  CUBLAS_CHECK(cublasDaxpy(cublas_handle_, a_num_cols_, &alpha, d_x_old, 1,
                           d_buffer_, 1));

  // 2. Compute delta_aty = aty_new - aty_old
  CUDA_CHECK(cudaMemcpy(d_buffer2_, d_aty_new, a_num_cols_ * sizeof(double),
                        cudaMemcpyDeviceToDevice));
  CUBLAS_CHECK(cublasDaxpy(cublas_handle_, a_num_cols_, &alpha, d_aty_old, 1,
                           d_buffer2_, 1));

  // 3. Compute dot product: delta_x' * delta_aty
  double result;
  CUBLAS_CHECK(cublasDdot(cublas_handle_, a_num_cols_, d_buffer_, 1, d_buffer2_,
                          1, &result));

  return result;
}

double PDLPSolver::computeDiffNormCuBLAS(const double* d_a, const double* d_b,
                                         HighsInt n) {
  // 1. Copy a to buffer: buffer = a
  CUDA_CHECK(
      cudaMemcpy(d_buffer_, d_a, n * sizeof(double), cudaMemcpyDeviceToDevice));

  // 2. buffer = buffer - b  (using cuBLAS axpy)
  double alpha = -1.0;
  CUBLAS_CHECK(cublasDaxpy(cublas_handle_, n, &alpha, d_b, 1, d_buffer_, 1));

  // 3. result = ||buffer||_2  (using cuBLAS nrm2)
  double norm;
  CUBLAS_CHECK(cublasDnrm2(cublas_handle_, n, d_buffer_, 1, &norm));

  return norm;
}
#endif
