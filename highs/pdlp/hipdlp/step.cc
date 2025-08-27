#include "step.hpp"
#include "linalg.hpp"
#include <cmath>
#include <iostream>
#include <limits>

namespace step {

static constexpr double kDivergentMovement = 1e10;

bool CheckNumericalStability(const std::vector<double>& delta_x, 
                             const std::vector<double>& delta_y,
                             double omega) {
    // Using omega as primal_weight for consistency
    double movement = ComputeMovement(delta_x, delta_y, omega);
    
    if (movement == 0.0) {
        std::cout << "Warning: Zero movement detected - numerical termination" << std::endl;
        return false;
    }
    
    if (movement > kDivergentMovement) {
        std::cout << "Warning: Divergent movement detected: " << movement << std::endl;
        return false;
    }
    
    return true;
}

double ComputeMovement(const std::vector<double>& delta_primal,
                      const std::vector<double>& delta_dual,
                      double primal_weight) {
    double primal_squared_norm = 0.0;
    for (const auto& val : delta_primal) {
        primal_squared_norm += val * val;
    }
    
    double dual_squared_norm = 0.0;
    for (const auto& val : delta_dual) {
        dual_squared_norm += val * val;
    }
    
    return (0.5 * primal_weight * primal_squared_norm) +
           (0.5 / primal_weight) * dual_squared_norm;
}

double ComputeNonlinearity(const std::vector<double>& delta_primal,
                           const std::vector<double>& delta_aty,
                           const std::vector<double>& aty_current,
                           const std::vector<double>& aty_new) {
    // Nonlinearity = |Δx' * Δ(A'y)|
    double nonlinearity = 0.0;
    for (size_t i = 0; i < delta_primal.size(); ++i) {
        nonlinearity += delta_primal[i] * delta_aty[i];
    }
    return std::abs(nonlinearity);
}

// --- Implementation of functions from step.hpp ---

void UpdateX(std::vector<double>& x_new, const std::vector<double>& x_current, 
             const HighsLp& lp, const std::vector<double>& y_current,
             double eta, double omega) 
{
    std::vector<double> ATy_cache(lp.num_col_);
    linalg::ATy(lp, y_current, ATy_cache);
    
    for (HighsInt i = 0; i < lp.num_col_; i++) {
        double gradient = lp.col_cost_[i] - ATy_cache[i];
        x_new[i] = linalg::project_box(x_current[i] - (eta / omega) * gradient, lp.col_lower_[i], lp.col_upper_[i]);
    }
}

void UpdateY(std::vector<double>& y_new, const std::vector<double>& y_current, 
             const HighsLp& lp, const std::vector<double>& ax_new,
             const std::vector<double>& ax_current, double eta, double omega) 
{
    for (HighsInt j = 0; j < lp.num_row_; j++) {
        double extr_ax = 2 * ax_new[j] - ax_current[j];
        
        //Row corresponds to an Ax = b constraint
        bool is_equality_constraint = (lp.row_lower_[j] == lp.row_upper_[j]);
        //Row corresponds to a Gx >= h constraint
        bool is_inequality_constraint = (lp.row_lower_[j] > -kHighsInf && lp.row_upper_[j] > -kHighsInf);

        double q;
        if (is_equality_constraint) {
            q = lp.row_lower_[j];
        } else if (is_inequality_constraint) {
            q = lp.row_lower_[j];
        } else {
            std::cerr << "lp.row_lower_[j]: " << lp.row_lower_[j] << std::endl;
            std::cerr << "lp.row_upper_[j]: " << lp.row_upper_[j] << std::endl;
            std::cerr << "Error: Invalid constraint type for row " << j << std::endl;
            continue; // Skip invalid constraints
        }

        double dual_update = y_current[j] + eta * omega * (q - extr_ax);

        if (is_equality_constraint) {
            y_new[j] = dual_update; // No projection needed for equality constraints
        } else {
            y_new[j] = linalg::project_non_negative(dual_update);
        }
    }
}

void UpdateIteratesFixed(const HighsLp& lp, const PrimalDualParams& params, double fixed_eta,
                        std::vector<double>& x_new, std::vector<double>& y_new,
                        std::vector<double>& ax_new, 
                        const std::vector<double>& x_current, 
                        const std::vector<double>& y_current,
                        const std::vector<double>& ax_current)
{
    UpdateX(x_new, x_current, lp, y_current, fixed_eta, params.omega);
    linalg::Ax(lp, x_new, ax_new);
    UpdateY(y_new, y_current, lp, ax_new, ax_current, fixed_eta, params.omega);
}

void UpdateIteratesAdaptive(const HighsLp& lp, const PrimalDualParams& params,
                           std::vector<double>& x_new, std::vector<double>& y_new,
                           std::vector<double>& ax_new, 
                           const std::vector<double>& x_current,
                           const std::vector<double>& y_current,
                           const std::vector<double>& ax_current,
                           const std::vector<double>& aty_current, 
                           double& current_eta,
                           int& step_size_iter_count)
{
    const double MIN_ETA = 1e-6;
    const double MAX_ETA = 1.0;
    
    bool accepted_step = false;
    int inner_iterations = 0;
    int num_rejected_steps = 0;
    
    while (!accepted_step) {
        inner_iterations++;
        
        if (inner_iterations >= 60) {
            std::cerr << "Warning: Adaptive line search exceeded 60 iterations." << std::endl;
            // Force accept the last candidate
            break;
        }
        
        // Compute candidate solution
        std::vector<double> x_candidate(lp.num_col_);
        std::vector<double> y_candidate(lp.num_row_);
        std::vector<double> ax_candidate(lp.num_row_);
        std::vector<double> aty_candidate(lp.num_col_);
        
        // Primal update
        UpdateX(x_candidate, x_current, lp, y_current, current_eta, params.omega);
        linalg::Ax(lp, x_candidate, ax_candidate);
        
        // Dual update 
        UpdateY(y_candidate, y_current, lp, ax_candidate, ax_current, current_eta, params.omega);
        linalg::ATy(lp, y_candidate, aty_candidate);
        
        // Compute deltas
        std::vector<double> delta_x(lp.num_col_);
        std::vector<double> delta_y(lp.num_row_);
        std::vector<double> delta_aty(lp.num_col_);
        
        for (size_t i = 0; i < x_candidate.size(); ++i) {
            delta_x[i] = x_candidate[i] - x_current[i];
        }
        for (size_t i = 0; i < y_candidate.size(); ++i) {
            delta_y[i] = y_candidate[i] - y_current[i];
        }
        for (size_t i = 0; i < aty_candidate.size(); ++i) {
            delta_aty[i] = aty_candidate[i] - aty_current[i];
        }
        
        // Check numerical stability
        if (!CheckNumericalStability(delta_x, delta_y, params.omega)) {
            std::cerr << "Numerical instability detected" << std::endl;
            current_eta *= 0.5;  // Drastically reduce step size
            continue;
        }
        
        // Compute movement and nonlinearity
        double movement = ComputeMovement(delta_x, delta_y, params.omega);
        double nonlinearity = ComputeNonlinearity(delta_x, delta_aty, aty_current, aty_candidate);
        
        // Compute step size limit
        double step_size_limit = (nonlinearity > 1e-12) ? 
                                 (movement / (2.0 * nonlinearity)) :  // in cupdlp-c, the factor is 1
                                 std::numeric_limits<double>::infinity();
        
        //std::cout << "Iteration " << inner_iterations << ": eta = " << current_eta 
        //          << ", movement = " << movement << ", nonlinearity = " << nonlinearity 
        //          << ", limit = " << step_size_limit << std::endl;

        
        if (current_eta <= step_size_limit) {
            // Accept the step
            x_new = x_candidate;
            y_new = y_candidate;
            ax_new = ax_candidate;
            accepted_step = true;
        } else {
            num_rejected_steps++;
        }
        
        // Compute new step size
        double first_term = (std::isinf(step_size_limit)) ? step_size_limit :
            (1.0 - std::pow(step_size_iter_count + 1.0, -params.adaptive_linesearch_params.step_size_reduction_exponent)) * step_size_limit;
        
        double second_term = (1.0 + std::pow(step_size_iter_count+ 1.0, 
            -params.adaptive_linesearch_params.step_size_growth_exponent)) * current_eta;
        
        current_eta = std::min(first_term, second_term);
        current_eta = std::max(MIN_ETA, std::min(MAX_ETA, current_eta));
    }
}

bool UpdateIteratesMalitskyPock(
    const HighsLp& lp, const PrimalDualParams& params,
    std::vector<double>& x_new, std::vector<double>& y_new,
    std::vector<double>& ax_new,
    const std::vector<double>& x_current,
    const std::vector<double>& y_current,
    const std::vector<double>& ax_current,
    const std::vector<double>& aty_current,
    double& current_eta,
    double& ratio_last_two_step_sizes,
    int& num_rejected_steps,
    bool first_iteration)  // Add this parameter to track first iteration
{
    // Step 1: Compute primal update first (without extrapolation)
    const double primal_step_size = current_eta / params.omega;
    
    std::vector<double> x_candidate(lp.num_col_);
    UpdateX(x_candidate, x_current, lp, y_current, primal_step_size, params.omega);
    
    // Compute Ax for the primal candidate
    std::vector<double> ax_candidate(lp.num_row_);
    linalg::Ax(lp, x_candidate, ax_candidate);
    
    // Step 2: Compute the dilating coefficient for new step size
    double dilating_coeff = 1.0 + (params.malitsky_pock_params.step_size_interpolation * 
                                   (std::sqrt(1.0 + ratio_last_two_step_sizes) - 1.0));
    
    double new_primal_step_size = primal_step_size * dilating_coeff;
    
    // Step 3: Line search loop
    bool accepted_step = false;
    int inner_iterations = 0;
    
    while (!accepted_step && inner_iterations < 60) {
        inner_iterations++;
        
        const double new_ratio = new_primal_step_size / primal_step_size;
        
        // CRITICAL: The dual weight should be primal_weight^2, and dual step = dual_weight * primal_step
        // In your case, primal_weight = omega, so dual_weight = omega^2
        const double dual_step_size = params.omega * params.omega * new_primal_step_size;
        
        // Step 4: Compute dual update with CORRECT extrapolation
        std::vector<double> y_candidate(lp.num_row_);
        
        // The extrapolation for dual uses the new ratio (theta in the paper)
        // Formula: y_new = y_current + dual_step * (b - A * x_extrapolated)
        // where x_extrapolated = x_candidate + new_ratio * (x_candidate - x_current)
        
        for (HighsInt j = 0; j < lp.num_row_; j++) {
            // Compute extrapolated Ax value
            // This is equivalent to: (1 + new_ratio) * ax_candidate - new_ratio * ax_current
            double ax_extrapolated = ax_candidate[j] + new_ratio * (ax_candidate[j] - ax_current[j]);
            
            bool is_equality = (lp.row_lower_[j] == lp.row_upper_[j]);
            double b = lp.row_lower_[j];  // Your 'q' is their 'b'
            
            // Dual gradient step
            double dual_update = y_current[j] + dual_step_size * (b - ax_extrapolated);
            
            if (is_equality) {
                y_candidate[j] = dual_update;
            } else {
                y_candidate[j] = linalg::project_non_negative(dual_update);
            }
        }
        
        // Step 5: Compute A^T * y_candidate for line search
        std::vector<double> aty_candidate(lp.num_col_);
        linalg::ATy(lp, y_candidate, aty_candidate);
        
        // Step 6: Compute norms for line search condition
        // Delta dual = y_candidate - y_current
        double delta_dual_norm_sq = 0.0;
        for (size_t i = 0; i < y_current.size(); ++i) {
            double diff = y_candidate[i] - y_current[i];
            delta_dual_norm_sq += diff * diff;
        }
        double delta_dual_norm = std::sqrt(delta_dual_norm_sq);
        
        // Delta dual product = A^T * (y_candidate - y_current)
        double delta_dual_prod_norm_sq = 0.0;
        for (size_t i = 0; i < aty_current.size(); ++i) {
            double diff = aty_candidate[i] - aty_current[i];
            delta_dual_prod_norm_sq += diff * diff;
        }
        double delta_dual_prod_norm = std::sqrt(delta_dual_prod_norm_sq);
        
        // Step 7: Line search condition (CRITICAL - must match reference exactly)
        if (params.omega * new_primal_step_size * delta_dual_prod_norm <= 
            params.malitsky_pock_params.linesearch_contraction_factor * delta_dual_norm) {
            
            // Accept the step
            current_eta = new_primal_step_size * params.omega;  // Store the combined step size
            ratio_last_two_step_sizes = new_ratio;
            
            x_new = x_candidate;
            y_new = y_candidate;
            ax_new = ax_candidate;
            
            // Check numerical stability
            std::vector<double> delta_x(lp.num_col_);
            std::vector<double> delta_y(lp.num_row_);
            
            for (size_t i = 0; i < x_candidate.size(); ++i) {
                delta_x[i] = x_candidate[i] - x_current[i];
            }
            for (size_t i = 0; i < y_candidate.size(); ++i) {
                delta_y[i] = y_candidate[i] - y_current[i];
            }
            
            if (!CheckNumericalStability(delta_x, delta_y, params.omega)) {
                std::cerr << "Numerical instability in Malitsky-Pock step" << std::endl;
                return false;
            }
            
            accepted_step = true;
            
            if (inner_iterations > 1) {
                std::cout << "Malitsky-Pock: accepted after " << inner_iterations 
                          << " line search iterations" << std::endl;
            }
            
        } else {
            // Reduce step size and try again
            new_primal_step_size *= params.malitsky_pock_params.step_size_downscaling_factor;
            
            if (inner_iterations % 10 == 0) {
                std::cout << "Malitsky-Pock line search: iteration " << inner_iterations 
                          << ", reducing step to " << new_primal_step_size << std::endl;
            }
        }
    }
    
    if (!accepted_step) {
        std::cerr << "Malitsky-Pock: Failed to find acceptable step after 60 iterations" << std::endl;
        return false;
    }
    
    num_rejected_steps += (inner_iterations - 1);
    
    return true;
}

} // namespace step
