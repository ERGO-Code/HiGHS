/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-07-09 14:54:26
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-02 18:32:14
 * @FilePath: /cupdlp-CPP/src/restart.cc
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
// restart.cc
#include "restart.hpp"
#include <algorithm>
#include <cmath>

// Initializes the restart scheme with parameters and initial results
void RestartScheme::Initialize(const PrimalDualParams& params, const SolverResults& results) {
    strategy_ = params.restart_strategy;
    fixed_restart_interval_ = params.fixed_restart_interval;
    last_restart_iter_ = 0;

      // Stepsize ratio,
  //  β = dBeta = dDualStep / dPrimalStep,
  //    in the paper, primal weight is the ω:
  //    ω = √β
    beta_ = std::sqrt(params.omega);

    last_restart_iter_ = 0;
    last_restart_score_ = std::numeric_limits<double>::infinity();
    last_candidate_score_ = std::numeric_limits<double>::infinity();
}

// Computes a weighted score to evaluate solution quality
double RestartScheme::ComputeRestartScore(const SolverResults& results) {
    double weight_squared = beta_;
    double primal_feas_sq = results.primal_feasibility * results.primal_feasibility;
    double dual_feas_sq = results.dual_feasibility * results.dual_feasibility;
    double gap_sq = results.duality_gap * results.duality_gap;

    return std::sqrt(weight_squared * primal_feas_sq + dual_feas_sq / weight_squared + gap_sq);
}

// Main logic to check if a restart is needed
RestartInfo RestartScheme::Check(int current_iter, const SolverResults& current_results, const SolverResults& average_results) {
    RestartInfo info;
    if (current_iter == last_restart_iter_) {
        double current_score = ComputeRestartScore(current_results);
        last_restart_score_ = current_score;
        last_candidate_score_ = current_score;
        return info; // Not enough progress since last restart
    }

    switch (strategy_) {
        case RestartStrategy::NO_RESTART:
            break; // Do nothing

        case RestartStrategy::FIXED_RESTART:
            if (current_iter - last_restart_iter_ >= fixed_restart_interval_) {
                info.should_restart = true;
                info.restart_to_average = true; // Fixed restarts typically use the average
            }
            break;

        case RestartStrategy::ADAPTIVE_RESTART: {
            double current_score = ComputeRestartScore(current_results);
            double average_score = ComputeRestartScore(average_results);
            
            // Choose the best candidate (current vs. average) based on the score
            double candidate_score = std::min(current_score, average_score);
            info.restart_to_average = (average_score < current_score);


            // 1. Artificial Restart Check 
            bool artificial_restart = false;
            if (current_iter > 64){
                artificial_restart = (current_iter - last_restart_iter_) >= (0.36 * current_iter);
            }
             
            // Adaptive restart conditions
            bool sufficient_decay = (candidate_score < sufficient_decay_factor_ * last_restart_score_);
            bool necessary_decay = (candidate_score < necessary_decay_factor_ * last_restart_score_) && (candidate_score > last_candidate_score_);

            if (artificial_restart) {
                std::cout << "Artificial restart triggered at iteration " << current_iter << std::endl;
                info.should_restart = true;
            } else if (sufficient_decay) {
                std::cout << "Sufficient decay triggered at iteration " << current_iter << std::endl;
                info.should_restart = true;
            } else if (necessary_decay) {
                std::cout << "Necessary decay triggered at iteration " << current_iter << std::endl;
                info.should_restart = true;
            } else {
                info.should_restart = false;
                last_candidate_score_ = std::min(last_candidate_score_, candidate_score);
            }
            
            break;
        }
    }

    if (info.should_restart) {
        last_restart_iter_ = current_iter;
        double candidate_score = std::min(ComputeRestartScore(current_results), ComputeRestartScore(average_results));
        last_restart_score_ = candidate_score;
        last_candidate_score_ = last_restart_score_; 
    }

    return info;
}