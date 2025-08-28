#include "pdhg.hpp"
#include "linalg.hpp"
#include "restart.hpp"
#include "step.hpp"
//#include "Highs.h"
//#include "highs_interface.hpp"
#include <iostream>
#include <tuple>
#include <cmath>
#include <random>

#define PDHG_CHECK_INTERVAL 40


using namespace std;

int PDLPSolver::GetIterationCount() const {
    return final_iter_count_;
}

void PDLPSolver::PreprocessLp(const HighsLp& original_lp, HighsLp& processed_lp) {
    logger_.info("Preprocessing LP using cupdlp formulation (slack variables for bounds)...");

    int nRows_orig = original_lp.num_row_;
    int nCols_orig = original_lp.num_col_;

    //std::vector<ConstraintType> constraint_types(nRows_orig);
    int num_new_cols = 0;
    int nEqs = 0;
    constraint_types_.resize(nRows_orig);
    constraint_new_idx_.resize(nRows_orig);

    // 1. First pass: Classify constraints and count slack variables needed
    for (int i = 0; i < nRows_orig; ++i) {
        bool has_lower = original_lp.row_lower_[i] > -kHighsInf;
        bool has_upper = original_lp.row_upper_[i] < kHighsInf;

        if (has_lower && has_upper) {
            if (original_lp.row_lower_[i] == original_lp.row_upper_[i]) {
                constraint_types_[i] = EQ;
            } else {
                constraint_types_[i] = BOUND;
                num_new_cols++; // Need one slack variable for each ranged constraint
            }
        } else if (has_lower) {
            constraint_types_[i] = GEQ;
        } else if (has_upper) {
            constraint_types_[i] = LEQ;
        } else {
            constraint_types_[i] = FREE; // Free rows become bounded equalities Ax=0, z=[-inf,inf]
            num_new_cols++;
        }
    }

    // 2. Set new dimensions
    processed_lp.num_col_ = nCols_orig + num_new_cols;
    processed_lp.num_row_ = nRows_orig;
    original_num_col_ = nCols_orig; // store for postsolve

    // 3. Resize all vectors in the processed_lp
    processed_lp.col_cost_.resize(processed_lp.num_col_);
    processed_lp.col_lower_.resize(processed_lp.num_col_);
    processed_lp.col_upper_.resize(processed_lp.num_col_);
    processed_lp.row_lower_.resize(processed_lp.num_row_);
    processed_lp.row_upper_.resize(processed_lp.num_row_);
    
    // 4. Populate costs and bounds for original and new slack variables
    for (int i = 0; i < nCols_orig; ++i) {
        processed_lp.col_cost_[i] = original_lp.col_cost_[i];
        processed_lp.col_lower_[i] = original_lp.col_lower_[i];
        processed_lp.col_upper_[i] = original_lp.col_upper_[i];
    }
    
    int current_slack_col = nCols_orig;
    for (int i = 0; i < nRows_orig; ++i) {
        if (constraint_types_[i] == BOUND || constraint_types_[i] == FREE) {
            processed_lp.col_cost_[current_slack_col] = 0.0;
            processed_lp.col_lower_[current_slack_col] = original_lp.row_lower_[i];
            processed_lp.col_upper_[current_slack_col] = original_lp.row_upper_[i];
            current_slack_col++;
        }
    }

    // 5. Formulate the new constraints matrix (A') and RHS (b')
    //
    // Take a copy of the original LP's constraint matrix since we
    // need it rowwise and it's const
    HighsSparseMatrix original_matrix = original_lp.a_matrix_;
    original_matrix.ensureRowwise();
    // Set up the processed constraint matrix as an empty row-wise
    // matrix that can have nCols_orig columns
    HighsSparseMatrix& processed_matrix = processed_lp.a_matrix_;
    processed_matrix.clear();
    processed_matrix.num_col_ = nCols_orig;
    processed_matrix.format_ = MatrixFormat::kRowwise;
    HighsSparseMatrix row;
    const double negative_one = -1.0;
    const double* negative_one_p = &negative_one;
    for (int i = 0; i < nRows_orig; ++i) {
      // Get the row from the original constraint matrix
      original_matrix.getRow(i, row);
      // Scale the row by -1 if it's a LEQ
      if (constraint_types_[i] == LEQ)
	row.scaleRows(negative_one_p);
      // Add the row to the processed constraint matrix
      processed_matrix.addRows(row);
    }
    // Convert the processed constraint matrix to column-wise orientation
    processed_matrix.ensureColwise();

    // Add slack variable entries (-1) for BOUND and FREE constraints
    current_slack_col = nCols_orig;
    // Set up a negated identity column as a column-wise matrix with
    // nRows_orig rows and a single column containing -1 in row 0 (for
    // now)
    HighsSparseMatrix col;
    assert(col.isColwise());
    col.num_col_ = 1;
    col.num_row_ = nRows_orig;
    col.start_.push_back(1);
    col.index_.push_back(0);
    col.value_.push_back(-1);
    for (int iRow = 0; iRow < nRows_orig; ++iRow) {
        if (constraint_types_[iRow] == BOUND || constraint_types_[iRow] == FREE) {
	  // Put the -1 in row iRow, and add the column to the matrix
	  col.index_[0] = iRow;
	  processed_matrix.addCols(col);
	  current_slack_col++;
        }
    }

    // 6. Set the new RHS (row bounds)
    for (int i = 0; i < nRows_orig; ++i) {
        switch (constraint_types_[i]) {
            case EQ:
                processed_lp.row_lower_[i] = original_lp.row_lower_[i];
                processed_lp.row_upper_[i] = original_lp.row_upper_[i];
                break;
            case GEQ:
                processed_lp.row_lower_[i] = original_lp.row_lower_[i];
                processed_lp.row_upper_[i] = kHighsInf;
                break;
            case LEQ:
                // Becomes -Ax >= -b
                processed_lp.row_lower_[i] = -original_lp.row_upper_[i];
                processed_lp.row_upper_[i] = kHighsInf;
                break;
            case BOUND:
            case FREE:
                // Becomes Ax - z = 0
                processed_lp.row_lower_[i] = 0.0;
                processed_lp.row_upper_[i] = 0.0;
                break;
        }
    }

    // 7. Convert COO matrix to CSC for the processed_lp
    //
    // Already achieved by construction
    logger_.info("Preprocessing complete. New dimensions: " + std::to_string(processed_lp.num_row_) + " rows, " + std::to_string(processed_lp.num_col_) + " cols.");
}

void PDLPSolver::Postsolve(const HighsLp& original_lp, HighsLp& processed_lp,
                           const std::vector<double>& x_processed, 
                           const std::vector<double>& y_processed,
                           HighsSolution& solution) {

    logger_.info("Post-solving the solution...");

    std::vector<double> x_unscaled = x_processed;
    std::vector<double> y_unscaled = y_processed;

    std::vector<double> dSlackPos_unscaled = dSlackPos_;
    std::vector<double> dSlackNeg_unscaled = dSlackNeg_;
    
    // 1. Unscale the solution vectors first
    scaling_.UnscaleSolution(x_unscaled, y_unscaled);
    const auto& row_scale = scaling_.GetRowScaling(); // Assumes getter exists
    const auto& col_scale = scaling_.GetColScaling(); // Assumes getter exists
    if (col_scale.size() == dSlackPos_unscaled.size()){
        for (size_t i = 0; i < col_scale.size(); ++i) {
            dSlackPos_unscaled[i] /= col_scale[i];
            dSlackNeg_unscaled[i] /= col_scale[i];
        }
    }

    // 2. Resize solution object to original dimensions
    solution.col_value.resize(original_num_col_);
    solution.row_value.resize(original_lp.num_row_);
    solution.col_dual.resize(original_num_col_);
    solution.row_dual.resize(original_lp.num_row_);

    // 3. Recover Primal Column Values (x)
    // This is the easy part: just take the first 'original_num_col_' elements.
    for (int i = 0; i < original_num_col_; ++i) {
        solution.col_value[i] = x_unscaled[i];
    }

    double final_primal_objective = original_lp.offset_;
    for (int i = 0; i < original_num_col_; ++i) {
        final_primal_objective += original_lp.col_cost_[i] * solution.col_value[i];
    }
    results_.primal_obj = final_primal_objective;

    // 4. Recover Dual Row Values (y)
    // This requires reversing the sign flip for LEQ constraints.
    for (int i = 0; i < original_lp.num_row_; ++i) {
        if (constraint_types_[i] == LEQ) {
            solution.row_dual[i] = -y_unscaled[i];
        } else {
            solution.row_dual[i] = y_unscaled[i];
        }
    }
    
    // 5. Recover Primal Row Values (Ax)
    // This requires re-calculating Ax with the unscaled solution and using the slack values.
    HighsLp unscaled_processed_lp = processed_lp;

    if(scaling_.IsScaled()) {
        // Unscale matrix, costs, and rhs
        for (int iCol = 0; iCol < unscaled_processed_lp.num_col_; ++iCol) {
            unscaled_processed_lp.col_cost_[iCol] /= col_scale[iCol];
            for (int iEl = unscaled_processed_lp.a_matrix_.start_[iCol]; iEl < unscaled_processed_lp.a_matrix_.start_[iCol + 1]; ++iEl) {
                int iRow = unscaled_processed_lp.a_matrix_.index_[iEl];
                unscaled_processed_lp.a_matrix_.value_[iEl] /= (row_scale[iRow] * col_scale[iCol]);
            }
        }
    }

    std::vector<double> ax_unscaled(unscaled_processed_lp.num_row_);
    linalg::Ax(unscaled_processed_lp, x_unscaled, ax_unscaled);

    int slack_variable_idx = original_num_col_;
    for (int i = 0; i < original_lp.num_row_; ++i) {
        if (constraint_types_[i] == BOUND || constraint_types_[i] == FREE) {
            solution.row_value[i] = x_unscaled[slack_variable_idx++];
        } else if (constraint_types_[i] == LEQ) {
            // We transformed Ax <= b to -Ax >= -b. The original row value is Ax.
            // The calculated ax_unscaled is for -Ax, so we flip the sign back.
            solution.row_value[i] = -ax_unscaled[i];
        } else { // EQ, GEQ
            solution.row_value[i] = ax_unscaled[i];
        }
    }

    // 6. Recover Dual Column Values (Reduced Costs)
    // The duals on the variable bounds l <= x <= u are the reduced costs.
    // In the PDLP framework, these are given by dSlackPos - dSlackNeg.
    for (int i = 0; i < original_num_col_; ++i) {
        solution.col_dual[i] = dSlackPos_unscaled[i] - dSlackNeg_unscaled[i];
    }

    solution.value_valid = true;
    solution.dual_valid = true; // We can now set this to true
    logger_.info("Post-solve complete.");
}

PDLPSolver::PDLPSolver(Logger& logger) : logger_(logger) {}

void PDLPSolver::Solve(HighsLp & original_lp, const PrimalDualParams& params, std::vector<double>& x, std::vector<double>& y) {
    Timer solver_timer;
    if (original_lp.num_row_ == 0 || original_lp.a_matrix_.start_.size() < 2) {
        std::cout << "Empty LP problem." << std::endl;
        return;
    }

    HighsSolution solution;

    // 1. Preprocess original LP to handle ranged constraints
    logger_.info("Preprocessing LP to handle ranged constraints...");
    HighsLp lp;
    PreprocessLp(original_lp, lp);

    // 2. Pass the PREPROCESSED LP to HiGHS and run presolve
    /*
    logger_.info("Running HiGHS presolve on the preprocessed problem...");
    HighsLp lp;
    if (!highs_interface_.passModelAndPresolve(preprocessed_lp, lp)) {
        logger_.info("Failed during HiGHS passModel or presolve on preprocessed LP.");
        results_.term_code = TerminationStatus::ERROR;
        return;
    }

    HighsModelStatus status1 = highs_interface_.getModelStatus();
    if (status1 == HighsModelStatus::kOptimal) {
        // If solved, postsolve to get the solution for the preprocessed problem
        logger_.info("Problem solved to optimality by HiGHS presolve.");

        return;
    }
*/
    // --- SCALING ---
    ScalingParams scaling_params;
    scaling_params.method = params.scaling_method;
    scaling_params.use_ruiz = params.use_ruiz_scaling;
    scaling_params.use_pc = params.use_pc_scaling;
    scaling_params.use_l2 = params.use_l2_scaling;
    scaling_params.ruiz_iterations = params.ruiz_iterations;
    scaling_params.ruiz_norm = params.ruiz_norm;
    scaling_params.pc_alpha = params.pc_alpha;
    scaling_.Initialize(lp);
    scaling_.ScaleProblem(lp, scaling_params);

    // --- 0.Using PowerMethod to estimate the largest eigenvalue ---
    double op_norm_sq = 1.0; // Default value
    HighsStatus status = PowerMethod(lp, op_norm_sq);
    if (status == HighsStatus::kWarning) {
        std::cout << "Warning: Power method for step size did not converge. Using default." << std::endl;
    }

    // Set step sizes based on the operator norm to ensure convergence
    // A safe choice satisfying eta * omega * ||A||^2 < 1
    const double fixed_eta = 0.99 / sqrt(op_norm_sq);
    current_eta_ = 1.0; // Initial step size for adaptive strategy
    const double omega = 0.99 / sqrt(op_norm_sq);
    PrimalDualParams working_params = params;
    working_params.omega = omega;
    working_params.eta = fixed_eta;
    std::cout << "Using power method step sizes: eta = " << fixed_eta << ", omega = " << omega << std::endl;
    printf("Initial step sizes from power method lambda = %g: primal = %g; dual = %g\n",
	   op_norm_sq, fixed_eta/omega, fixed_eta*omega);
    // --- 1. Initialization ---
    Initialize(lp, x, y); // Sets initial x, y and results_
    restart_scheme_.Initialize(params, results_);

    // Use member variables for iterates for cleaner state management
    x_current_ = x;
    y_current_ = y;
    num_rejected_steps_ = 0; 
    bool first_malitsky_iteration = true;
    ratio_last_two_step_sizes_ = 1.0; 
    bool using_malitsky_averaging = (params.step_size_strategy == StepSizeStrategy::MALITSKY_POCK);
    bool primal_average_initialized = false;

    // We need Ax for the y-update and for computing residuals
    std::vector<double> Ax_current(lp.num_row_, 0.0);
    std::vector<double> ATy_current(lp.num_col_, 0.0);
    linalg::Ax(lp, x_current_, Ax_current);
    linalg::ATy(lp, y_current_, ATy_current);
    std::vector<double> Ax_new(lp.num_row_, 0.0);
    std::vector<double> ATy_new(lp.num_col_, 0.0);

    logger_.print_iteration_header();

    // --- 2. Main PDHG Loop ---
    // A single loop handles max iterations, convergence, and restarts.
    for (int iter = 0; iter < params.max_iterations; ++iter) {
        if (solver_timer.read() > params.time_limit) {
            logger_.info("Time limit reached.");
            final_iter_count_ = iter;
            
            // Set termination status
            results_.term_code = TerminationStatus::TIMEOUT;
            
            // Ensure the final summary is accurate by using the last computed results
            // (The `results_` member should already hold the last iteration's data)
            
            // Unscale solution before returning
            scaling_.UnscaleSolution(x, y);

            // Update objective with final unscaled solution
            double primal_objective = 0.0;
            for (HighsInt i = 0; i < original_lp.num_col_; ++i) {
                primal_objective += original_lp.col_cost_[i] * x[i];
            }
            results_.primal_obj = primal_objective;

            return; // Exit the function
        }

        // --- 3. Core PDHG Update Step ---
        bool step_success = true;
        
        switch (params.step_size_strategy) {
            case StepSizeStrategy::FIXED:
                step::UpdateIteratesFixed(lp, working_params, fixed_eta, x, y,
                                         Ax_new, x_current_, y_current_, Ax_current);
                break;
                
            case StepSizeStrategy::ADAPTIVE:
                step::UpdateIteratesAdaptive(lp, working_params, x, y, Ax_new, 
                                           x_current_, y_current_, Ax_current, 
                                           ATy_current, current_eta_, iter);
                break;
                
            case StepSizeStrategy::MALITSKY_POCK:
                step_success = step::UpdateIteratesMalitskyPock(
                    lp, working_params, x, y, Ax_new,
                    x_current_, y_current_, Ax_current, ATy_current,
                    current_eta_, ratio_last_two_step_sizes_, num_rejected_steps_,
                    first_malitsky_iteration);
                
                if (!step_success) {
                    std::cerr << "Malitsky-Pock step failed at iteration " << iter << std::endl;
                    // Reset to average and terminate
                    x = x_avg_;
                    y = y_avg_;
                    //scaling_.UnscaleSolution(x, y);
                    return;
                }
        }

        linalg::ATy(lp, y, ATy_new);
        
        // --- 4. Update Average Iterates ---
        // The number of iterations since the last restart
        int inner_iter = iter - restart_scheme_.GetLastRestartIter();
        UpdateAverageIterates(x, y, working_params, inner_iter);

        // --- 5. Convergence Check ---
        bool bool_checking = (iter < 10) || (iter == (params.max_iterations -1)) || (iter > 0 && (iter % PDHG_CHECK_INTERVAL == 0));
        if (bool_checking) { 
            // To check convergence on the average, you need A*x_avg and A^T*y_avg
            SolverResults current_results;
            SolverResults average_results;
            std::vector<double> Ax_avg(lp.num_row_, 0.0), ATy_avg(lp.num_col_, 0.0);
            linalg::Ax(lp, x_avg_, Ax_avg);
            linalg::ATy(lp, y_avg_, ATy_avg);

            bool average_converged = CheckConvergence(lp, x_avg_, y_avg_, Ax_avg, ATy_avg, 
                                                  params.tolerance, average_results);
        
            logger_.print_iteration_stats(iter, average_results, current_eta_);
           
            if (average_converged) {
                logger_.info("Average solution converged in " + std::to_string(iter) + " iterations.");
                final_iter_count_ = iter;
                results_ = average_results;

                // Use the average solution
                x = x_avg_;
                y = y_avg_;
                
                // Unscale solution
                //scaling_.UnscaleSolution(x, y);
                
                // Calculate primal objective on original problem
                double primal_objective = 0.0;
                for (HighsInt i = 0; i < original_lp.num_col_; ++i) {
                    primal_objective += original_lp.col_cost_[i] * x[i];
                }
                
                // Update results with the converged metrics
                results_ = average_results;
                results_.term_code = TerminationStatus::OPTIMAL;
                Postsolve(original_lp, lp, x, y, solution);
                return;
            }

            bool current_converged = CheckConvergence(lp, x, y, Ax_new, ATy_new, 
                                                    params.tolerance, current_results);
            
            // Log iteration stats (using current results for display)
            logger_.print_iteration_stats(iter, current_results, current_eta_);

            if (current_converged) {
                logger_.info("Current solution converged in " + std::to_string(iter) + " iterations.");
                final_iter_count_ = iter;
                results_ = current_results; 
                
                // Unscale solution
                //scaling_.UnscaleSolution(x, y);
                
                // Calculate primal objective on original problem
                double primal_objective = 0.0;
                for (HighsInt i = 0; i < original_lp.num_col_; ++i) {
                    primal_objective += original_lp.col_cost_[i] * x[i];
                }
                
                // Update results with the converged metrics
                results_ = current_results;
                results_.term_code = TerminationStatus::OPTIMAL;
                Postsolve(original_lp, lp, x, y, solution);
                return;
            }
        

            // --- 6. Restart Check ---
            RestartInfo restart_info = restart_scheme_.Check(iter, current_results, average_results);

            if (restart_info.should_restart) {
                logger_.info("Restarting at iteration " + std::to_string(iter));
                std::vector<double> restart_x, restart_y;
                if (restart_info.restart_to_average) {
                    // Restart from the average iterate
                    restart_x = x_avg_;
                    restart_y = y_avg_;
                } else {
                    // Restart from the current iterate (for adaptive scheme)
                    restart_x = x;
                    restart_y = y;
                }

                // Perform the primal weight update using z^{n,0} and z^{n-1,0}
                //PDHG_Compute_Step_Size_Ratio(working_params,
                //                         restart_x, restart_y,
                //                         x_at_last_restart_, y_at_last_restart_);

                x_at_last_restart_ = restart_x; // Current becomes the new last
                y_at_last_restart_ = restart_y;

                // Set x_current_ and y_current_ for the next iteration's starting point
                x_current_ = restart_x;
                y_current_ = restart_y;

                // Reset the average iterate accumulation by passing -1 as a signal
                UpdateAverageIterates(x_current_, y_current_, working_params, -1);
            } else {
                // Standard update: next iterate's starting point is the current iterate
                x_current_ = x;
                y_current_ = y;
            }

            Ax_current = Ax_new; 
            ATy_current = ATy_new;
        }
    }

    // --- 7. Handle Max Iterations Reached ---
    logger_.info("Max iterations reached without convergence.");
    final_iter_count_ = params.max_iterations;
    //scaling_.UnscaleSolution(x, y);

    // Calculate primal objective on original problem
    double primal_objective = 0.0;
    for (HighsInt i = 0; i < original_lp.num_col_; ++i) {
        primal_objective += original_lp.col_cost_[i] * x[i];
    }

    // Update the final pieces of the results
    results_.primal_obj = primal_objective;
    results_.term_code = TerminationStatus::TIMEOUT;
    Postsolve(original_lp, lp, x_avg_, y_avg_, solution);
    return;
}

void PDLPSolver::Initialize(const HighsLp & lp, std::vector<double>& x, std::vector<double>& y) {
    // Initialize x and y based on the LP problem
    x.resize(lp.num_col_,0.0);
    y.resize(lp.num_row_, 0.0);

    x_at_last_restart_ = x;
    y_at_last_restart_ = y;
    x_sum_.resize(lp.num_col_, 0.0);
    y_sum_.resize(lp.num_row_, 0.0);

    sum_weights_ = 0.0;

    x_avg_.resize(lp.num_col_, 0.0);
    y_avg_.resize(lp.num_row_, 0.0);


    Ax_cache_.resize(lp.num_row_);
    ATy_cache_.resize(lp.num_col_);
    K_times_x_diff_.resize(lp.num_row_);
    dSlackPos_.resize(lp.num_col_, 0.0);
    dSlackNeg_.resize(lp.num_col_, 0.0);
}

// Update primal weight
void PDLPSolver::PDHG_Compute_Step_Size_Ratio(
    PrimalDualParams& working_params,
    const std::vector<double>& x_n_0,
    const std::vector<double>& y_n_0,
    const std::vector<double>& x_n_minus_1_0,
    const std::vector<double>& y_n_minus_1_0
) {
    // 1. Calculate the L2 norm of the difference between current and last-restart iterates.
    double primal_diff_norm = linalg::diffTwoNorm(x_n_0, x_n_minus_1_0);
    double dual_diff_norm = linalg::diffTwoNorm(y_n_0, y_n_minus_1_0);

    // 2. Update the primal weight (beta = omega^2) if movements are significant.
    if (std::min(primal_diff_norm, dual_diff_norm) > 1e-10) {
        double beta_update_ratio = dual_diff_norm / primal_diff_norm;
        double current_beta = working_params.omega * working_params.omega;

        double dLogBetaUpdate = 0.5 * std::log(beta_update_ratio) + 0.5 * std::log(std::sqrt(current_beta));
        double new_beta = std::exp(2.0 * dLogBetaUpdate);
        
        working_params.omega = std::sqrt(new_beta);
    }
}


void PDLPSolver::UpdateAverageIterates(const std::vector<double>& x, const std::vector<double>& y, const PrimalDualParams& params, int inner_iter) {
    if (inner_iter == 0 || inner_iter == -1) {
        x_avg_ = x;
        y_avg_ = y;

        x_sum_ = std::vector<double>(x.size(), 0.0);
        y_sum_ = std::vector<double>(y.size(), 0.0);

        sum_weights_ = 0.0;
    } else {
        for (size_t i = 0; i < x.size(); ++i) {
            x_sum_[i]  += x[i] * params.eta;
        }
        for (size_t i = 0; i < y.size(); ++i) {
            y_sum_[i]  += y[i] * params.eta;
        }
        sum_weights_ += params.eta;
        for (size_t i = 0; i < x.size(); ++i) {
            x_avg_[i] = x_sum_[i] / sum_weights_;
        }
        for (size_t i = 0; i < y.size(); ++i) {
            y_avg_[i] = y_sum_[i] / sum_weights_;
        }
    }
}

// lambda = c - proj_{\Lambda}(c - K^T y)
std::vector<double> PDLPSolver::ComputeLambda(const HighsLp& lp, const std::vector<double>& y, const std::vector<double>& ATy_vector) {
    std::vector<double> lambda(lp.num_col_, 0.0);
    for (HighsInt i = 0; i < lp.num_col_; ++i) {
        double residual = lp.col_cost_[i] - ATy_vector[i];

        if (lp.col_lower_[i] <= -kHighsInf && lp.col_upper_[i] >= kHighsInf) {
            // Case 1: Unbounded variable (l_i = -inf, u_i = +inf) -> lambda_i = 0
            lambda[i] = 0;
        } else if (lp.col_lower_[i] <= -kHighsInf) {
            // Case 2: Only upper bound (l_i = -inf, u_i is finite) -> lambda_i in R^-
            lambda[i] = std::min(0.0, residual);
        } else if (lp.col_upper_[i] >= kHighsInf) {
            // Case 3: Only lower bound (l_i is finite, u_i = +inf) -> lambda_i in R^+
            lambda[i] = std::max(0.0, residual);
        } else {
            // Case 4: Boxed variable (l_i and u_i are finite) -> lambda_i in R
            lambda[i] = residual;
        }
    }
    return lambda;
}

std::pair<double,double> PDLPSolver::ComputePrimalFeasibility(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& Ax_vector) {
    double primal_feasibility_squared = 0.0;
    std::vector<double> Ax_minus_b(lp.num_row_, 0.0);
    for (HighsInt i = 0; i < lp.num_row_; ++i) {
        double Ax = Ax_vector[i];
        if (lp.row_lower_[i] == lp.row_upper_[i]){ 
            // Equality constraint
            Ax_minus_b[i] = Ax - lp.row_lower_[i];
        } else { // Inequality constraint
            Ax_minus_b[i] = std::max(0.0, Ax - lp.row_upper_[i]);
        }
        primal_feasibility_squared += Ax_minus_b[i] * Ax_minus_b[i];
    }

    double primal_feasibility = sqrt(primal_feasibility_squared);

    // Compute norm of q
    double q_norm = 0.0;
    for (HighsInt i = 0; i < lp.num_row_; ++i) {
        q_norm += lp.row_lower_[i] * lp.row_lower_[i];
    }
    q_norm = sqrt(q_norm);

    return std::make_pair(primal_feasibility, q_norm);
}

void PDLPSolver::ComputeDualSlacks(const HighsLp& lp, const std::vector<double>& ATy_vector) {
    // Ensure vectors are correctly sized
    if (dSlackPos_.size() != lp.num_col_) dSlackPos_.resize(lp.num_col_);
    if (dSlackNeg_.size() != lp.num_col_) dSlackNeg_.resize(lp.num_col_);

    for (HighsInt i = 0; i < lp.num_col_; ++i) {
        double reduced_cost = lp.col_cost_[i] - ATy_vector[i];
        
        // Slack for lower bound l_i <= x_i (dual variable is non-negative)
        dSlackPos_[i] = (lp.col_lower_[i] > -kHighsInf) ? std::max(0.0, reduced_cost) : 0.0;
        
        // Slack for upper bound x_i <= u_i (dual variable is non-positive, so slack is non-negative)
        dSlackNeg_[i] = (lp.col_upper_[i] < kHighsInf) ? std::max(0.0, -reduced_cost) : 0.0;
    }
}

std::pair<double, double> PDLPSolver::ComputeDualFeasibility(const HighsLp& lp, const std::vector<double>& ATy_vector) {
    double dual_feasibility_squared = 0.0;
    std::vector<double> dual_residual(lp.num_col_);

    for (HighsInt i = 0; i < lp.num_col_; ++i) {
        // The dual feasibility residual is c - A'y - (dSlackPos - dSlackNeg)
        // where (dSlackPos - dSlackNeg) is the dual variable for the box constraints.
        dual_residual[i] = lp.col_cost_[i] - ATy_vector[i] - (dSlackPos_[i] - dSlackNeg_[i]);
        dual_feasibility_squared += dual_residual[i] * dual_residual[i];
    }

    double dual_feasibility = sqrt(dual_feasibility_squared);
    double c_norm = linalg::vector_norm(lp.col_cost_);
    return std::make_pair(dual_feasibility, c_norm);
}

std::tuple<double, double, double, double, double> PDLPSolver::ComputeDualityGap(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& lambda) {
     double qTy = 0.0;
    for (int i = 0; i < lp.num_row_; ++i) {
        // Assumes Gx >= h and Ax = b are combined into a single matrix constraint
        // where the lower bounds represent h and b. This needs to be consistent
        // with how the problem is transformed.
        if (lp.row_lower_[i] > -kHighsInf) {
            qTy += lp.row_lower_[i] * y[i];
        } else {
            // For equality constraints Ax=b, often represented as Ax>=b and -Ax>=-b.
            // If b is in row_lower_, this is fine. If it's split, logic may need adjustment.
            qTy += lp.row_upper_[i] * y[i];
        }
    }

    double lT_lambda_plus = 0.0;
    double uT_lambda_minus = 0.0;

    for (int i = 0; i < lp.num_col_; ++i) {
        if (lp.col_lower_[i] > -kHighsInf) {
            lT_lambda_plus += lp.col_lower_[i] * std::max(0.0, lambda[i]); // l^T * lambda^+
        }
        if (lp.col_upper_[i] < kHighsInf) {
            uT_lambda_minus += lp.col_upper_[i] * std::max(0.0, -lambda[i]); // u^T * lambda^-
        }
    }

    double dual_objective = qTy + lT_lambda_plus - uT_lambda_minus;

    double cTx = 0.0;
    for (int i = 0; i < lp.num_col_; ++i) {
        cTx += lp.col_cost_[i] * x[i];
    }

    double duality_gap = abs(dual_objective - cTx);
    
    // For the termination criteria, you need the components of the dual objective
    return std::make_tuple(duality_gap, qTy, lT_lambda_plus, uT_lambda_minus, cTx);
}

static double ComputeWeightedNorm(const std::vector<double>& x1, const std::vector<double>& y1,
                            const std::vector<double>& x2, const std::vector<double>& y2, double omega) {
    double norm = 0.0;
    for (size_t i = 0; i < x1.size(); ++i) {
        double diff_x = x1[i] - x2[i];
        norm += (1.0 / omega) * diff_x * diff_x;
    }
    for (size_t i = 0; i < y1.size(); ++i) {
        double diff_y = y1[i] - y2[i];
        norm += omega * diff_y * diff_y;
    }
    return std::sqrt(norm);
}

double ComputeNormalizedDualityGap(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& x_ref, const std::vector<double>& y_ref, double omega, const std::vector<double>& ax, const std::vector<double>& aty){
    // 1. Calculate radius r = ||z - z_ref||_omega
    double r = ComputeWeightedNorm(x, y, x_ref, y_ref, omega);
    if (r < 1e-9) return 0.0;

    // 2. Calculate Lagrangian gradient components
    // g_x = c - A^T * y
    std::vector<double> g_x(lp.num_col_, 0.0);

    //to do 
    return 0.0; 
}
 
bool PDLPSolver::CheckConvergence(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& ax_vector, const std::vector<double>& aty_vector, double epsilon, SolverResults& results) {
    ComputeDualSlacks(lp, aty_vector);
    std::vector<double> lambda = ComputeLambda(lp, y, aty_vector);

    auto [primal_feasibility, q_norm] = ComputePrimalFeasibility(lp, x, ax_vector);
    results.primal_feasibility = primal_feasibility;

    auto [dual_feasibility, c_norm] = ComputeDualFeasibility(lp, aty_vector);
    results.dual_feasibility = dual_feasibility;
    auto [duality_gap, qTy, lTlambda_plus, uTlambda_minus, cTx] = ComputeDualityGap(lp, x, y, lambda);
    results.duality_gap = duality_gap;

    bool primal_feasible = primal_feasibility <= epsilon * (1 + q_norm);
    bool dual_feasible = dual_feasibility <= epsilon * (1 + c_norm);

    double dual_objective_value = qTy + lTlambda_plus - uTlambda_minus;
    double denominator = 1 + abs(dual_objective_value);
    if (denominator == 0) {
        logger_.info("Denominator is zero, cannot check duality gap.");
    }
    bool duality_gap_small = duality_gap <= epsilon * denominator;


    return primal_feasible && dual_feasible && duality_gap_small;
}

HighsStatus PDLPSolver::PowerMethod(HighsLp &lp, double& op_norm_sq) {
    if (lp.num_col_ == 0 || lp.num_row_ == 0) {
        op_norm_sq = 1.0;
        return HighsStatus::kOk;
    }

    // Parameters for the power method
    const int max_iter = 20;
    const double tol = 1e-6;

    // Initialize a random vector x
    std::vector<double> x_vec(lp.num_col_);
    std::random_device rd;
    std::mt19937 engine_fixed_seed(12345); //gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    for (HighsInt i = 0; i < lp.num_col_; ++i) {
        x_vec[i] = dis(engine_fixed_seed);
    }
    linalg::normalize(x_vec); // Assumes a normalize function in linalg

    const HighsInt kYanyuPowerMethod = 0;
    const HighsInt kYanyuPowerMethodDev = 1;
    const HighsInt kCuPdlpAATPowerMethod = 2;
    
    const HighsInt power_method =
       kYanyuPowerMethod;
      // kYanyuPowerMethodDev;
      // kCuPdlpAATPowerMethod;
    printf("Power method: %s\n", power_method == kYanyuPowerMethod ? "Yanyu" : power_method == kYanyuPowerMethodDev ? "Yanyu dev" : "CuPdlp-C");
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
    int log_iters = log_level == LogLevel::kVerbose || log_level == LogLevel::kDebug;
     log_iters = 1;

    if (log_iters) printf("It       lambda   dl_lambda\n");

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
        linalg::ATy(lp, y_vec, z_vec); // Note: ATy computes A' * vector

        // Estimate the squared operator norm (largest eigenvalue of A'A)
        op_norm_sq = linalg::dot(x_vec, z_vec); // Assumes a dot product function in linalg
        
        // Check for convergence
        if (std::abs(op_norm_sq - op_norm_sq_old) < tol * op_norm_sq) {
            return HighsStatus::kOk;
        }
	double dl_op_norm_sq = std::fabs(op_norm_sq - op_norm_sq_old);
        op_norm_sq_old = op_norm_sq;

        // Prepare for the next iteration
        linalg::normalize(z_vec); // Normalize the result
        x_vec = z_vec;
	if (log_iters) printf("%2d %12.6g %11.4g\n", iter, op_norm_sq, dl_op_norm_sq);
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
	if (log_iters) printf("%2d %12.6g %11.4g\n", iter, lambda, dl_lambda);
      }
    }
    if (power_method != kYanyuPowerMethod) op_norm_sq = lambda;
    // If the method did not converge within max_iter
    return HighsStatus::kWarning;
}
