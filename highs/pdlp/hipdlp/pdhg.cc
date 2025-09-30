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

#include "HConst.h"
#include "linalg.hpp"
#include "pdlp/cupdlp/cupdlp.h"  // For pdlpLogging
#include "restart.hpp"
#include "defs.hpp"

#define PDHG_CHECK_INTERVAL 40
static constexpr double kDivergentMovement = 1e10;

using namespace std;

void PDLPSolver::preprocessLp() {
  logger_.info(
      "Preprocessing LP using cupdlp formulation (slack variables for "
      "bounds)...");

  int nRows_orig = original_lp_->num_row_;
  int nCols_orig = original_lp_->num_col_;

  if (original_lp_->a_matrix_.isRowwise()) {
    logger_.info("Original LP matrix must be in column-wise format,");
  }

  int num_new_cols = 0;
  int nEqs = 0;
  constraint_types_.resize(nRows_orig);

  // 1. First pass: Classify constraints and count slack variables needed
  for (int i = 0; i < nRows_orig; ++i) {
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
  int eq_idx = 0;
  int ineq_idx = nEqs;
  for (int i = 0; i < nRows_orig; ++i) {
    if (constraint_types_[i] == EQ || constraint_types_[i] == BOUND ||
        constraint_types_[i] == FREE) {
      constraint_new_idx_[i] = eq_idx++;
    } else {
      constraint_new_idx_[i] = ineq_idx++;
    }
  }

  // 4. Create is_equality_row_ vector based on NEW ordering
  is_equality_row_.resize(nRows_orig, false);
  for (int i = 0; i < nRows_orig; ++i) {
    int new_row_idx = constraint_new_idx_[i];
    is_equality_row_[new_row_idx] = (constraint_types_[i] == EQ || 
                                      constraint_types_[i] == BOUND || 
                                      constraint_types_[i] == FREE);
  }

  // 5. Populate costs and bounds for original and new slack variables
  for (int i = 0; i < nCols_orig; ++i) {
    processed_lp.col_cost_[i] = original_lp_->col_cost_[i];
    processed_lp.col_lower_[i] = original_lp_->col_lower_[i];
    processed_lp.col_upper_[i] = original_lp_->col_upper_[i];
  }

  int current_slack_col = nCols_orig;
  for (int i = 0; i < nRows_orig; ++i) {
    if (constraint_types_[i] == BOUND || constraint_types_[i] == FREE) {
      processed_lp.col_cost_[current_slack_col] = 0.0;
      processed_lp.col_lower_[current_slack_col] = original_lp_->row_lower_[i];
      processed_lp.col_upper_[current_slack_col] = original_lp_->row_upper_[i];
      current_slack_col++;
    }
  }

  // 6. Set the new RHS (row bounds) in PERMUTED order
  for (int i = 0; i < nRows_orig; ++i) {
    int new_row_idx = constraint_new_idx_[i];
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
  for (int col = 0; col < nCols_orig; ++col) {
    // For each column, add entries in the new row order
    std::vector<std::pair<int, double>> col_entries;
    
    for (int el = original_matrix.start_[col]; el < original_matrix.start_[col + 1]; ++el) {
      int old_row = original_matrix.index_[el];
      int new_row = constraint_new_idx_[old_row];
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
  for (int i = 0; i < nRows_orig; ++i) {
    if (constraint_types_[i] == BOUND || constraint_types_[i] == FREE) {
      int row_idx = constraint_new_idx_[i];
      processed_matrix.index_.push_back(row_idx);
      processed_matrix.value_.push_back(-1.0);
      processed_matrix.start_[current_slack_col + 1] = processed_matrix.index_.size();
      current_slack_col++;
    }
  }

  num_eq_rows_ = nEqs;
  scaling_.passLp(&processed_lp);
  unscaled_processed_lp_ = processed_lp;

  // 8. Compute and store norms of unscaled cost and rhs
  unscaled_c_norm_ = linalg::vector_norm(processed_lp.col_cost_);
  unscaled_rhs_norm_ = linalg::vector_norm(processed_lp.row_lower_);

  logger_.info("Preprocessing complete. New dimensions: " +
               std::to_string(processed_lp.num_row_) + " rows, " +
               std::to_string(processed_lp.num_col_) + " cols.");
  logger_.info("Unscaled norms: ||c|| = " + std::to_string(unscaled_c_norm_) +
               ", ||b|| = " + std::to_string(unscaled_rhs_norm_));
}

PostSolveRetcode PDLPSolver::postprocess(HighsSolution& solution) {
  logger_.info("Post-solving the solution...");

  if (original_lp_ == nullptr) {
    return PostSolveRetcode::DIMENSION_MISMATCH;
  }

  if (x_current_.empty() || y_current_.empty()) {
    return PostSolveRetcode::INVALID_SOLUTION;
  }

  if (x_current_.size() != lp_.num_col_ || y_current_.size() != lp_.num_row_) {
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
  for (int i = 0; i < original_lp_->num_col_; ++i) {
    if (i >= (int)x_current_.size()) {
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
  for (int i = 0; i < original_num_col_; ++i) {
    final_primal_objective +=
        original_lp_->col_cost_[i] * solution.col_value[i];
  }
  results_.primal_obj = final_primal_objective;

  // 4. Recover Dual Row Values (y)
  for (int i = 0; i < original_lp_->num_row_; ++i) {
    if (constraint_types_[i] == LEQ) {
      solution.row_dual[i] = -y_current_[i];
    } else {
      solution.row_dual[i] = y_current_[i];
    }
  }

  // 5. Recover Primal Row Values (Ax)
  const std::vector<double>& row_scale = scaling_.GetRowScaling();
  const std::vector<double>& col_scale = scaling_.GetColScaling();
  if (scaling_.IsScaled()) {
    // Unscale matrix, costs, and rhs
    for (int iCol = 0; iCol < unscaled_processed_lp_.num_col_; ++iCol) {
      unscaled_processed_lp_.col_cost_[iCol] /= col_scale[iCol];
      for (int iEl = unscaled_processed_lp_.a_matrix_.start_[iCol];
           iEl < unscaled_processed_lp_.a_matrix_.start_[iCol + 1]; ++iEl) {
        int iRow = unscaled_processed_lp_.a_matrix_.index_[iEl];
        unscaled_processed_lp_.a_matrix_.value_[iEl] /=
            (row_scale[iRow] * col_scale[iCol]);
      }
    }
  }

  std::vector<double> ax_current_(unscaled_processed_lp_.num_row_);
  linalg::Ax(unscaled_processed_lp_, x_current_, ax_current_);

  int slack_variable_idx = original_num_col_;
  for (int i = 0; i < original_lp_->num_row_; ++i) {
    if (constraint_types_[i] == BOUND || constraint_types_[i] == FREE) {
      solution.row_value[i] = x_current_[slack_variable_idx++];
    } else if (constraint_types_[i] == LEQ) {
      // We transformed Ax <= b to -Ax >= -b. The original row value is Ax.
      // The calculated ax_current_ is for -Ax, so we flip the sign back.
      solution.row_value[i] = -ax_current_[i];
    } else {  // EQ, GEQ
      solution.row_value[i] = ax_current_[i];
    }
  }

  // 6. Recover Dual Column Values (Reduced Costs)
  // The duals on the variable bounds l <= x <= u are the reduced costs.
  // In the PDLP framework, these are given by dSlackPos - dSlackNeg.
  for (int i = 0; i < original_num_col_; ++i) {
    if (i >= (int)dSlackPos_.size() || i >= (int)dSlackNeg_.size()) {
      logger_.info("Index " + std::to_string(i) +
                   " out of bounds for dSlackPos/Neg of size " +
                   std::to_string(dSlackPos_.size()));
      return PostSolveRetcode::DIMENSION_MISMATCH;
    }
    solution.col_dual[i] = dSlackPos_[i] - dSlackNeg_[i];
  }

  solution.value_valid = true;  // to do
  solution.dual_valid = true;
  logger_.info("Post-solve complete.");

  return PostSolveRetcode::OK;
}

void PDLPSolver::solve(std::vector<double>& x, std::vector<double>& y) {
  Timer solver_timer;

  const HighsLp& lp = lp_;

  debug_pdlp_log_file_ = fopen("HiPDLP.log", "w");
  assert(debug_pdlp_log_file_);

  // --- 0. Using PowerMethod to estimate the largest eigenvalue ---
  const double op_norm_sq = PowerMethod();
  stepsize_ =
      InitializeStepSizesPowerMethod(op_norm_sq);
  const double fixed_eta = 0.99 / sqrt(op_norm_sq);
  PrimalDualParams working_params = params_;

  working_params.omega = std::sqrt(stepsize_.dual_step / stepsize_.primal_step);
  working_params.eta = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);
  current_eta_ = working_params.eta;

  highsLogUser(params_.log_options_, HighsLogType::kInfo,
               "Using power method step sizes: eta = %g, omega = %g\n",
               working_params.eta, working_params.omega);

  highsLogUser(params_.log_options_, HighsLogType::kInfo,
               "Initial step sizes from power method lambda = %g: primal = %g; "
               "dual = %g\n",
               stepsize_.power_method_lambda, stepsize_.primal_step,
               stepsize_.dual_step);

  // --- 1. Initialization ---
  restart_scheme_.passLogOptions(&working_params.log_options_);
  restart_scheme_.passDebugPdlpLogFile(debug_pdlp_log_file_);
  restart_scheme_.passDebugPdlpData(&debug_pdlp_data_);
  Initialize();  // Sets initial x, y and results_
  restart_scheme_.passParams(&working_params);
  restart_scheme_.Initialize(results_);

  // Use member variables for iterates for cleaner state management
  x_current_ = x;
  y_current_ = y;
  num_rejected_steps_ = 0;
  bool first_malitsky_iteration = true;
  ratio_last_two_step_sizes_ = 1.0;
  bool using_malitsky_averaging =
      (params_.step_size_strategy == StepSizeStrategy::MALITSKY_POCK);
  bool primal_average_initialized = false;

  // Initialize vectors for matrix-vector products
  std::vector<double> Ax_avg(lp.num_row_, 0.0);
  std::vector<double> ATy_avg(lp.num_col_, 0.0);

  // Store iterates at last restart for primal weight update
  x_at_last_restart_ = x_current_;
  y_at_last_restart_ = y_current_;

  logger_.print_iteration_header();

  // --- 2. Main PDHG Loop ---
  debugPdlpIterHeaderLog(debug_pdlp_log_file_);
  debugPdlpDataInitialise(&debug_pdlp_data_);
  debug_pdlp_data_.ax_average_norm = 0.0;
  debug_pdlp_data_.aty_average_norm = 0.0;
  debug_pdlp_data_.x_average_norm = 0.0;
  debug_pdlp_data_.ax_norm = 0.0;

  for (int iter = 0; iter < params_.max_iterations; ++iter) {
    debugPdlpIterLog(debug_pdlp_log_file_, iter, &debug_pdlp_data_,
                     restart_scheme_.getBeta(), stepsize_.primal_step,
                     stepsize_.dual_step);

    // Check time limit
    if (solver_timer.read() > params_.time_limit) {
      logger_.info("Time limit reached.");
      final_iter_count_ = iter;
      results_.term_code = TerminationStatus::TIMEOUT;
      return;
    }

    // --- 3. Convergence and Restart Check (BEFORE iterate update) ---
    bool bool_checking = (iter < 10) ||
                         (iter == (params_.max_iterations - 1)) ||
                         (iter > 0 && (iter % PDHG_CHECK_INTERVAL == 0));

    if (bool_checking) {
      ComputeAverageIterate(Ax_avg, ATy_avg);
      // Compute matrix-vector products for current and average iterates
      linalg::Ax(lp, x_current_, Ax_cache_);
      linalg::ATy(lp, y_current_, ATy_cache_);

      // Reset the average iterate accumulation
      int inner_iter = iter - restart_scheme_.GetLastRestartIter();

      // Compute residuals and convergence metrics
      SolverResults current_results;
      SolverResults average_results;

      // Compute residuals for current iterate
      bool current_converged = CheckConvergence(
          iter, x_current_, y_current_, Ax_cache_, ATy_cache_,
          params_.tolerance, current_results, "[L]");

      // Compute residuals for average iterate
      bool average_converged =
          CheckConvergence(iter, x_avg_, y_avg_, Ax_avg, ATy_avg,
                           params_.tolerance, average_results, "[A]");

      debugPdlpIterHeaderLog(debug_pdlp_log_file_);

      // Print iteration statistics
      logger_.print_iteration_stats(iter, average_results, current_eta_);

      // Check for convergence
      if (current_converged) {
        logger_.info("Current solution converged in " + std::to_string(iter) +
                     " iterations.");
        final_iter_count_ = iter;
        x = x_current_;
        y = y_current_;
        results_ = current_results;
        results_.term_code = TerminationStatus::OPTIMAL;
        return;
      }

      if (average_converged) {
        logger_.info("Average solution converged in " + std::to_string(iter) +
                     " iterations.");
        final_iter_count_ = iter;
        x = x_avg_;
        y = y_avg_;
        results_ = average_results;
        results_.term_code = TerminationStatus::OPTIMAL;
        return;
      }

      // --- 4. Restart Check (using computed results) ---
      RestartInfo restart_info =
          restart_scheme_.Check(iter, current_results, average_results);

      if (restart_info.should_restart) {
        std::vector<double> restart_x, restart_y;
        if (restart_info.restart_to_average) {
          // Restart from the average iterate
          restart_x = x_avg_;
          restart_y = y_avg_;
        } else {
          // Restart from the current iterate
          restart_x = x_current_;
          restart_y = y_current_;
        }

        // Perform the primal weight update using z^{n,0} and z^{n-1,0}
        PDHG_Compute_Step_Size_Ratio(working_params, restart_x, restart_y,
                                     x_at_last_restart_, y_at_last_restart_);
        current_eta_ = working_params.eta;
        restart_scheme_.passParams(&working_params);
        restart_scheme_.UpdateBeta(working_params.omega * working_params.omega);

        x_at_last_restart_ = restart_x;  // Current becomes the new last
        y_at_last_restart_ = restart_y;

        std::fill(x_sum_.begin(), x_sum_.end(), 0.0);
        std::fill(y_sum_.begin(), y_sum_.end(), 0.0);
        sum_weights_ = 0.0;

        // Set x_current_ and y_current_ for the next iteration's starting point
        x_current_ = restart_x;
        y_current_ = restart_y;

        // Recompute Ax and ATy for the restarted iterates
        linalg::Ax(lp, x_current_, Ax_cache_);
        debug_pdlp_data_.ax_norm = linalg::vector_norm(Ax_cache_);
        linalg::ATy(lp, y_current_, ATy_cache_);
        debug_pdlp_data_.aty_norm = linalg::vector_norm(ATy_cache_);
        restart_scheme_.SetLastRestartIter(iter);
      }
    } else {
      // If not checking, still need Ax and ATy for the update
      linalg::Ax(lp, x_current_, Ax_cache_);
      debug_pdlp_data_.ax_norm = linalg::vector_norm(Ax_cache_);
      linalg::ATy(lp, y_current_, ATy_cache_);
      debug_pdlp_data_.aty_norm = linalg::vector_norm(ATy_cache_);
    }

    // --- 5. Core PDHG Update Step ---
    bool step_success = true;

    // Store current iterates before update (for next iteration's x_current_,
    // y_current_)
    x_next_ = x_current_;
    y_next_ = y_current_;

    switch (params_.step_size_strategy) {
      case StepSizeStrategy::FIXED:
        UpdateIteratesFixed();
        break;

      case StepSizeStrategy::ADAPTIVE:
        UpdateIteratesAdaptive(iter);
        break;

      case StepSizeStrategy::MALITSKY_POCK:
        step_success = UpdateIteratesMalitskyPock(
            first_malitsky_iteration);

        if (!step_success) {
          std::cerr << "Malitsky-Pock step failed at iteration " << iter
                    << std::endl;
          // Reset to average and terminate
          x = x_avg_;
          y = y_avg_;
          return;
        }
    }

    // Compute ATy for the new iterate
    linalg::ATy(lp, y_next_, ATy_cache_);

    // --- 6. Update Average Iterates ---
    // The number of iterations since the last restart
    int inner_iter = iter - restart_scheme_.GetLastRestartIter();
    UpdateAverageIterates(x_next_, y_next_, working_params, inner_iter);;

    // --- 7. Prepare for next iteration ---
    x_current_ = x_next_;
    y_current_ = y_next_;
    // iteration
  }

  // --- 8. Handle Max Iterations Reached ---
  logger_.info("Max iterations reached without convergence.");
  final_iter_count_ = params_.max_iterations;

  // Return the average solution
  x = x_avg_;
  y = y_avg_;

  results_.term_code = TerminationStatus::TIMEOUT;
  return;
}

void PDLPSolver::solveReturn() {
  if (debug_pdlp_log_file_) fclose(debug_pdlp_log_file_);
}

void PDLPSolver::Initialize() {
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

  Ax_cache_.resize(lp_.num_row_, 0.0);
  ATy_cache_.resize(lp_.num_col_, 0.0);
  Ax_next_.resize(lp_.num_row_, 0.0);
  ATy_next_.resize(lp_.num_col_, 0.0);
  K_times_x_diff_.resize(lp_.num_row_, 0.0);
  dSlackPos_.resize(lp_.num_col_, 0.0);
  dSlackNeg_.resize(lp_.num_col_, 0.0);
}

// Update primal weight
void PDLPSolver::PDHG_Compute_Step_Size_Ratio(
    PrimalDualParams& working_params, const std::vector<double>& x_n_0,
    const std::vector<double>& y_n_0, const std::vector<double>& x_n_minus_1_0,
    const std::vector<double>& y_n_minus_1_0) {
  // 1. Calculate the L2 norm of the difference between current and last-restart
  // iterates.
  //double primal_diff_norm = linalg::diffTwoNorm(x_n_0, x_n_minus_1_0);
  double primal_diff_norm = linalg::diffTwoNorm(x_at_last_restart_, x_current_);
  double norm_xLastRestart = linalg::vector_norm(x_at_last_restart_);
  double norm_xCurrent = linalg::vector_norm(x_current_);

  //double dual_diff_norm = linalg::diffTwoNorm(y_n_0, y_n_minus_1_0);
  double dual_diff_norm = linalg::diffTwoNorm(y_at_last_restart_, y_current_);

  double dMeanStepSize = std::sqrt(stepsize_.primal_step *
                                       stepsize_.dual_step);


  // 2. Update the primal weight (beta = omega^2) if movements are significant.
  if (std::min(primal_diff_norm, dual_diff_norm) > 1e-10) {
    double beta_update_ratio = dual_diff_norm / primal_diff_norm;
    double old_beta = working_params.omega * working_params.omega;

    double dLogBetaUpdate = 0.5 * std::log(beta_update_ratio) +
                            0.5 * std::log(std::sqrt(old_beta));
    stepsize_.beta = std::exp(2.0 * dLogBetaUpdate);
  }

  stepsize_.primal_step = dMeanStepSize / std::sqrt(stepsize_.beta);
  stepsize_.dual_step = stepsize_.primal_step * stepsize_.beta;
  working_params.eta = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);
  working_params.omega = std::sqrt(stepsize_.beta);
}

void PDLPSolver::UpdateAverageIterates(const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const PrimalDualParams& params,
                                       int inner_iter) {

  double dMeanStepSize = std::sqrt(stepsize_.primal_step * stepsize_.dual_step);

  // Debug: print what we're adding
    double x_norm_before = linalg::vector_norm(x_sum_);

  for (size_t i = 0; i < x.size(); ++i) {
    x_sum_[i] += x[i] * dMeanStepSize;
  }
  for (size_t i = 0; i < y.size(); ++i) {
    y_sum_[i] += y[i] * dMeanStepSize;
  }

  sum_weights_ += dMeanStepSize;
}

void PDLPSolver::ComputeAverageIterate(std::vector<double>& ax_avg,
                                       std::vector<double>& aty_avg) {
  double dPrimalScale = sum_weights_ > 1e-10 ? 1.0 / sum_weights_ : 1.0;
  double dDualScale = sum_weights_ > 1e-10 ? 1.0 / sum_weights_ : 1.0;

  for (size_t i = 0; i < x_avg_.size(); ++i) {
    x_avg_[i] = x_sum_[i] * dPrimalScale;
  }
  for (size_t i = 0; i < y_avg_.size(); ++i) {
    y_avg_[i] = y_sum_[i] * dDualScale; 
  }

  debug_pdlp_data_.x_average_norm = linalg::vector_norm_squared(x_avg_);
  linalg::Ax(lp_, x_avg_, ax_avg);
  linalg::ATy(lp_, y_avg_, aty_avg);
  debug_pdlp_data_.ax_average_norm = linalg::vector_norm_squared(ax_avg);
  debug_pdlp_data_.aty_average_norm = linalg::vector_norm_squared(aty_avg);
}

// lambda = c - proj_{\Lambda}(c - K^T y)
std::vector<double> PDLPSolver::ComputeLambda(
    const std::vector<double>& y, const std::vector<double>& ATy_vector) {
  std::vector<double> lambda(lp_.num_col_, 0.0);
  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    double residual = lp_.col_cost_[i] - ATy_vector[i];

    if (lp_.col_lower_[i] <= -kHighsInf && lp_.col_upper_[i] >= kHighsInf) {
      // Case 1: Unbounded variable (l_i = -inf, u_i = +inf) -> lambda_i = 0
      lambda[i] = 0;
    } else if (lp_.col_lower_[i] <= -kHighsInf) {
      // Case 2: Only upper bound (l_i = -inf, u_i is finite) -> lambda_i in R^-
      lambda[i] = std::min(0.0, residual);
    } else if (lp_.col_upper_[i] >= kHighsInf) {
      // Case 3: Only lower bound (l_i is finite, u_i = +inf) -> lambda_i in R^+
      lambda[i] = std::max(0.0, residual);
    } else {
      // Case 4: Boxed variable (l_i and u_i are finite) -> lambda_i in R
      lambda[i] = residual;
    }
  }
  return lambda;
}

double PDLPSolver::ComputePrimalFeasibility(
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

void PDLPSolver::ComputeDualSlacks(const std::vector<double>& ATy_vector) {
  // Ensure vectors are correctly sized
  if (dSlackPos_.size() != lp_.num_col_) dSlackPos_.resize(lp_.num_col_);
  if (dSlackNeg_.size() != lp_.num_col_) dSlackNeg_.resize(lp_.num_col_);

  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    double reduced_cost = lp_.col_cost_[i] - ATy_vector[i];

    // Slack for lower bound l_i <= x_i (dual variable is non-negative)
    dSlackPos_[i] =
        (lp_.col_lower_[i] > -kHighsInf) ? std::max(0.0, reduced_cost) : 0.0;

    // Slack for upper bound x_i <= u_i (dual variable is non-positive, so slack
    // is non-negative)
    dSlackNeg_[i] =
        (lp_.col_upper_[i] < kHighsInf) ? std::max(0.0, -reduced_cost) : 0.0;
  }
}

double PDLPSolver::ComputeDualFeasibility(
    const std::vector<double>& ATy_vector) {
  ComputeDualSlacks(ATy_vector);  // This updates dSlackPos_ and dSlackNeg_

  std::vector<double> dual_residual(lp_.num_col_);

  for (HighsInt i = 0; i < lp_.num_col_; ++i) {
    // Matching CUPDLP: c - A'y - dSlackPos + dSlackNeg
    dual_residual[i] = lp_.col_cost_[i] - ATy_vector[i] 
                      - dSlackPos_[i] + dSlackNeg_[i];
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
PDLPSolver::ComputeDualityGap(const std::vector<double>& x,
                              const std::vector<double>& y,
                              const std::vector<double>& lambda) {
  double qTy = 0.0;
  for (int i = 0; i < lp_.num_row_; ++i) {
    // Assumes Gx >= h and Ax = b are combined into a single matrix constraint
    // where the lower bounds represent h and b. This needs to be consistent
    // with how the problem is transformed.
    if (lp_.row_lower_[i] > -kHighsInf) {
      qTy += lp_.row_lower_[i] * y[i];
    } else {
      // For equality constraints Ax=b, often represented as Ax>=b and -Ax>=-b.
      // If b is in row_lower_, this is fine. If it's split, logic may need
      // adjustment.
      qTy += lp_.row_upper_[i] * y[i];
    }
  }

  double lT_lambda_plus = 0.0;
  double uT_lambda_minus = 0.0;

  for (int i = 0; i < lp_.num_col_; ++i) {
    if (lp_.col_lower_[i] > -kHighsInf) {
      lT_lambda_plus +=
          lp_.col_lower_[i] * std::max(0.0, lambda[i]);  // l^T * lambda^+
    }
    if (lp_.col_upper_[i] < kHighsInf) {
      uT_lambda_minus +=
          lp_.col_upper_[i] * std::max(0.0, -lambda[i]);  // u^T * lambda^-
    }
  }

  double dual_objective = qTy + lT_lambda_plus - uT_lambda_minus;

  double cTx = 0.0;
  for (int i = 0; i < lp_.num_col_; ++i) {
    cTx += lp_.col_cost_[i] * x[i];
  }

  double duality_gap = abs(dual_objective - cTx);

  // For the termination criteria, you need the components of the dual objective
  return std::make_tuple(duality_gap, qTy, lT_lambda_plus, uT_lambda_minus,
                         cTx);
}

double PDLPSolver::ComputeDualObjective(const std::vector<double>& y) {
  double dual_obj = 0.0;

  // Compute b'y (or rhs'y in cuPDLP notation)
  for (int i = 0; i < lp_.num_row_; ++i) {
    dual_obj += lp_.row_lower_[i] * y[i];
  }

  // Add contribution from lower bounds: l'*slackPos
  for (int i = 0; i < lp_.num_col_; ++i) {
    if (lp_.col_lower_[i] > -kHighsInf) {
      dual_obj += lp_.col_lower_[i] * dSlackPos_[i];
    }
  }

  // Subtract contribution from upper bounds: u'*slackNeg
  for (int i = 0; i < lp_.num_col_; ++i) {
    if (lp_.col_upper_[i] < kHighsInf) {
      dual_obj -= lp_.col_upper_[i] * dSlackNeg_[i];
    }
  }

  return dual_obj;
}

bool PDLPSolver::CheckConvergence(const int iter, const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& ax_vector,
                                  const std::vector<double>& aty_vector,
                                  double epsilon, SolverResults& results,
                                  const char* type) {
  // Compute dual slacks first
  ComputeDualSlacks(aty_vector);

  // Compute primal feasibility
  double primal_feasibility = 
      ComputePrimalFeasibility(ax_vector);
  results.primal_feasibility = primal_feasibility;

  // Compute dual feasibility
  double dual_feasibility = ComputeDualFeasibility(aty_vector);
  results.dual_feasibility = dual_feasibility;

  // Compute objectives
  double primal_obj = 0.0;
  for (int i = 0; i < lp_.num_col_; ++i) {
    primal_obj += lp_.col_cost_[i] * x[i];
  }
  results.primal_obj = primal_obj;

  double dual_obj = ComputeDualObjective(y);
  results.dual_obj = dual_obj;

  // Compute duality gap
  double duality_gap = primal_obj - dual_obj;
  results.duality_gap = std::abs(duality_gap);

  // Compute relative gap (matching cuPDLP formula)
  double relative_obj_gap =
      std::abs(duality_gap) / (1.0 + std::abs(primal_obj) + std::abs(dual_obj));
  results.relative_obj_gap = relative_obj_gap;

  debugPdlpFeasOptLog(debug_pdlp_log_file_, iter, primal_obj, dual_obj,
                      relative_obj_gap, primal_feasibility / (1.0 + unscaled_rhs_norm_),
                      dual_feasibility / (1.0 + unscaled_c_norm_), type);

  // Check convergence criteria (matching cuPDLP)
  bool primal_feasible = primal_feasibility < epsilon * (1.0 + unscaled_rhs_norm_);
  bool dual_feasible = dual_feasibility < epsilon * (1.0 + unscaled_c_norm_);
  bool gap_small = relative_obj_gap < epsilon;

  return primal_feasible && dual_feasible && gap_small;
}

double PDLPSolver::PowerMethod() {
  const HighsLp& lp = lp_;
  double op_norm_sq = 1.0;
  if (lp.num_col_ == 0 || lp.num_row_ == 0) return op_norm_sq;

  // Parameters for the power method
  const int max_iter = 20;
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
  highsLogUser(params_.log_options_, HighsLogType::kInfo, "Power method: %s\n",
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
  int log_iters =
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
  for (int iter = 0; iter < max_iter; ++iter) {
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
                     "%2d %12.6g %11.4g\n", int(iter), op_norm_sq,
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
                     "%2d %12.6g %11.4g\n", int(iter), lambda, dl_lambda);
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

  params_.initialise();
  //  params.eta = 0; Not set in parse_options_file
  //  params.omega = 0; Not set in parse_options_file
  params_.tolerance = options.pdlp_optimality_tolerance;
  if (options.kkt_tolerance != kDefaultKktTolerance) {
    params_.tolerance = options.kkt_tolerance;
  }
  params_.max_iterations = options.pdlp_iteration_limit;
  params_.device_type = Device::CPU;
  // HiPDLP has its own timer, so set its time limit according to
  // the time remaining with respect to the HiGHS time limit (if
  // finite)s
  double time_limit = options.time_limit;
  if (time_limit < kHighsInf) {
    time_limit -= timer.read();
    time_limit = std::max(0.0, time_limit);
  }
  params_.time_limit = time_limit;

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
  if ((options.pdlp_features_off & kPdlpRestartOff) == 0) {
    // Use restart: now see which
    if (options.pdlp_restart_strategy == kPdlpRestartStrategyFixed) {
      params_.restart_strategy = RestartStrategy::FIXED_RESTART;
    } else if (options.pdlp_restart_strategy == kPdlpRestartStrategyAdaptive) {
      params_.restart_strategy = RestartStrategy::ADAPTIVE_RESTART;
    }
  }
  //  params_.fixed_restart_interval = 0; Not set in parse_options_file
  //  params_.use_halpern_restart = false; Not set in parse_options_file

  params_.step_size_strategy = StepSizeStrategy::FIXED;
  if ((options.pdlp_features_off & kPdlpAdaptiveStepSizeOff) == 0) {
    // Use adaptive step size: now see which
    if (options.pdlp_step_size_strategy == kPdlpStepSizeStrategyAdaptive) {
      params_.step_size_strategy = StepSizeStrategy::ADAPTIVE;
    } else if (options.pdlp_step_size_strategy ==
               kPdlpStepSizeStrategyMalitskyPock) {
      params_.step_size_strategy = StepSizeStrategy::MALITSKY_POCK;
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
      dSlackPos_[i] /= col_scale[i];
      dSlackNeg_[i] /= col_scale[i];
    }
  }
}

void PDLPSolver::logSummary() {
  logger_.print_summary(results_, final_iter_count_, total_timer.read());
}

void PrimalDualParams::initialise() {
  this->eta = 0;
  this->omega = 0;
  this->tolerance = 0;
  this->max_iterations = 0;
  this->device_type = Device::CPU;
  this->time_limit = 3600.0;
  this->restart_strategy = RestartStrategy::NO_RESTART;
  this->fixed_restart_interval = 0;
  this->use_halpern_restart = false;
  this->use_ruiz_scaling = false;
  this->use_pc_scaling = false;
  this->use_l2_scaling = false;
  this->ruiz_iterations = 10;
  this->ruiz_norm = INFINITY;
  this->pc_alpha = 1.0;
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

StepSizeConfig PDLPSolver::InitializeStepSizesPowerMethod(double op_norm_sq) {
  StepSizeConfig config;
  config.power_method_lambda = op_norm_sq;
  double cost_norm = linalg::compute_cost_norm(lp_, 2.0);
  double rhs_norm = linalg::compute_rhs_norm(lp_, 2.0);

  highsLogUser(params_.log_options_, HighsLogType::kInfo,
               "Cost norm: %g, RHS norm: %g\n", cost_norm, rhs_norm);
  config.beta = cost_norm * cost_norm / (rhs_norm * rhs_norm + 1e-10);

  const double safety_factor = 0.8;
  double base_step = safety_factor / std::sqrt(op_norm_sq);

  config.primal_step = base_step / std::sqrt(config.beta);
  config.dual_step = base_step * std::sqrt(config.beta);
  return config;
}

std::vector<double> PDLPSolver::UpdateX() {
  std::vector<double> x_new(lp_.num_col_);
  linalg::ATy(lp_, y_current_, ATy_cache_);
  debug_pdlp_data_.aty_norm = linalg::vector_norm(ATy_cache_);
  for (HighsInt i = 0; i < lp_.num_col_; i++) {
    double gradient = lp_.col_cost_[i] - ATy_cache_[i];
    x_new[i] = linalg::project_box(x_current_[i] - stepsize_.primal_step * gradient,
                                   lp_.col_lower_[i], lp_.col_upper_[i]);
  }
  return x_new;
}

std::vector<double> PDLPSolver::UpdateY() {
  std::vector<double> y_new(lp_.num_row_);
  for (HighsInt j = 0; j < lp_.num_row_; j++) {
    double extr_ax = 2 * Ax_next_[j] - Ax_cache_[j];
    bool is_equality = (lp_.row_lower_[j] == lp_.row_upper_[j]);
    double q = lp_.row_lower_[j];
    double dual_update = y_current_[j] + stepsize_.dual_step * (q - extr_ax);
    y_new[j] = is_equality ? dual_update : linalg::project_non_negative(dual_update);
  }
  return y_new;
}

void PDLPSolver::UpdateIteratesFixed() {
  x_next_ = UpdateX();
  linalg::Ax(lp_, x_next_, Ax_next_);
  y_next_ = UpdateY();
}

void PDLPSolver::UpdateIteratesAdaptive(int& step_size_iter_count) {
  const double MIN_ETA = 1e-6;
  const double MAX_ETA = 1.0;

  bool accepted_step = false;
  int inner_iterations = 0;
  int num_rejected_steps = 0;

while (!accepted_step) {
    inner_iterations++;

    if (inner_iterations >= 60) {
      std::cerr << "Warning: Adaptive line search exceeded 60 iterations."
                << std::endl;
      // Force accept the last candidate
      break;
    }

    // Compute candidate solution
    std::vector<double> x_candidate(lp_.num_col_);
    std::vector<double> y_candidate(lp_.num_row_);
    std::vector<double> ax_candidate(lp_.num_row_);
    std::vector<double> aty_candidate(lp_.num_col_);

    // Primal update
    x_candidate = UpdateX();
    linalg::Ax(lp_, x_candidate, ax_candidate);

    // Dual update
    y_candidate = UpdateY();
    linalg::ATy(lp_, y_candidate, aty_candidate);

    // Compute deltas
    std::vector<double> delta_x(lp_.num_col_);
    std::vector<double> delta_y(lp_.num_row_);
    std::vector<double> delta_aty(lp_.num_col_);

    for (size_t i = 0; i < x_candidate.size(); ++i) {
      delta_x[i] = x_candidate[i] - x_current_[i];
    }
    for (size_t i = 0; i < y_candidate.size(); ++i) {
      delta_y[i] = y_candidate[i] - y_current_[i];
    }
    for (size_t i = 0; i < aty_candidate.size(); ++i) {
      delta_aty[i] = aty_candidate[i] - ATy_cache_[i];
    }

    // Check numerical stability
    if (!CheckNumericalStability(delta_x, delta_y)) {
      std::cerr << "Numerical instability detected" << std::endl;
      current_eta_ *= 0.5;  // Drastically reduce step size
      continue;
    }

    // Compute movement and nonlinearity
    double movement = ComputeMovement(delta_x, delta_y);
    double nonlinearity =
        ComputeNonlinearity(delta_x, delta_aty);

    // Compute step size limit
    double step_size_limit = (nonlinearity > 1e-12)
                                 ? (movement / (2.0 * nonlinearity))
                                 :  // in cupdlp-c, the factor is 1
                                 std::numeric_limits<double>::infinity();

    /*
    highsLogUser(*log_options_, HighsLogType::kInfo,
                 "Iteration %d: eta = %g, movement = %g nonlinearity = %g, limit
    = %g\n", int(inner_iterations), current_eta, movement, nonlinearity.
    step_size_limit);
    */
    if (current_eta_ <= step_size_limit) {
      // Accept the step
      x_next_ = x_candidate;
      y_next_ = y_candidate;
      Ax_next_ = ax_candidate;
      ATy_cache_ = aty_candidate;
      accepted_step = true;
    } else {
      num_rejected_steps++;
    }

    // Compute new step size
    double first_term =
        (std::isinf(step_size_limit))
            ? step_size_limit
            : (1.0 - std::pow(step_size_iter_count + 1.0,
                              -params_.adaptive_linesearch_params
                                   .step_size_reduction_exponent)) *
                  step_size_limit;

    double second_term =
        (1.0 +
         std::pow(
             step_size_iter_count + 1.0,
             -params_.adaptive_linesearch_params.step_size_growth_exponent)) *
        current_eta_;

    current_eta_ = std::min(first_term, second_term);
    current_eta_ = std::max(MIN_ETA, std::min(MAX_ETA, current_eta_));
  }
}

bool PDLPSolver::UpdateIteratesMalitskyPock(
    bool first_iteration) {
    std::vector<double> x_candidate(lp_.num_col_);
    return true; 
}

// =============================================================================
//  SECTION 5: Step Helper Methods (from step.cc)
// =============================================================================

bool PDLPSolver::CheckNumericalStability(const std::vector<double>& delta_x,
                                         const std::vector<double>& delta_y) {
  double movement = ComputeMovement(delta_x, delta_y);

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
    return true; // Placeholder
}

double PDLPSolver::ComputeMovement(const std::vector<double>& delta_primal,
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

double PDLPSolver::ComputeNonlinearity(const std::vector<double>& delta_primal,
                                       const std::vector<double>& delta_aty) {
  // Nonlinearity = |Δx' * Δ(A'y)|
  double nonlinearity = 0.0;
  for (size_t i = 0; i < delta_primal.size(); ++i) {
    nonlinearity += delta_primal[i] * delta_aty[i];
  }
  return std::abs(nonlinearity);
}
