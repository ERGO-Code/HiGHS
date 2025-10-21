/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/restart.cc
 * @brief
 */
// restart.cc
#include "restart.hpp"

#include <algorithm>
#include <cmath>

#include "io/HighsIO.h"          // For pdlpLogging
#include "pdlp/cupdlp/cupdlp.h"  // For pdlpLogging

// Helper function to compute the score
static double compute_score(double beta, double primal_feas, double dual_feas, double duality_gap) {
    return std::sqrt(beta * primal_feas * primal_feas +
                     dual_feas * dual_feas / beta +
                     duality_gap * duality_gap);
}

// Initializes the restart scheme with parameters and initial results
void RestartScheme::Initialize(const SolverResults& results) {
  strategy_ = params_->restart_strategy;
  fixed_restart_interval_ = params_->fixed_restart_interval;
  last_restart_iter_ = 0;

  // Stepsize ratio,
  //  β = dBeta = dDualStep / dPrimalStep,
  //    in the paper, primal weight is the ω:
  //    ω = √β
  beta_ = params_->omega * params_->omega;
}

/*
// Main logic to check if a restart is needed
RestartInfo RestartScheme::Check(int current_iter,
                                 const SolverResults& current_results,
                                 const SolverResults& average_results) {
  RestartInfo info;

  if (current_iter <= last_restart_iter_) {
    double current_score = ComputeRestartScore(current_results);
    last_restart_score_ = current_score;
    last_candidate_score_ = current_score;
    return info;  // Not enough progress since last restart
  }

  switch (strategy_) {
    case RestartStrategy::NO_RESTART:
      break;  // Do nothing

    case RestartStrategy::FIXED_RESTART:
      if (current_iter - last_restart_iter_ >= fixed_restart_interval_) {
        info.should_restart = true;
        info.restart_to_average =
            true;  // Fixed restarts typically use the average
      }
      break;

    case RestartStrategy::ADAPTIVE_RESTART: {
      double current_score = ComputeRestartScore(current_results);
      double average_score = ComputeRestartScore(average_results);
      //debugPdlpRestartLog(debug_pdlp_log_file_, current_iter, current_score,
      //                    average_score);

      // Choose the best candidate (current vs. average) based on the score
      double candidate_score = std::min(current_score, average_score);
      info.restart_to_average = (average_score < current_score);

      // 1. Artificial Restart Check
      bool artificial_restart = false;
      //      if (current_iter > 64) {
      artificial_restart =
          (current_iter - last_restart_iter_) >= (0.36 * current_iter);
      //      }

      // Adaptive restart conditions
      bool sufficient_decay =
          (candidate_score < sufficient_decay_factor_ * last_restart_score_);
      bool necessary_decay =
          (candidate_score < necessary_decay_factor_ * last_restart_score_) &&
          (candidate_score > last_candidate_score_);

      if (artificial_restart) {
        highsLogUser(*log_options_, HighsLogType::kInfo,
                     "Artificial restart triggered at iteration %d\n",
                     int(current_iter));
        info.should_restart = true;
      } else if (sufficient_decay) {
        highsLogUser(*log_options_, HighsLogType::kInfo,
                     "Sufficient decay triggered at iteration %d\n",
                     int(current_iter));
        info.should_restart = true;
      } else if (necessary_decay) {
        highsLogUser(*log_options_, HighsLogType::kInfo,
                     "Necessary decay triggered at iteration %d\n",
                     int(current_iter));
        info.should_restart = true;
      } else {
        info.should_restart = false;
        last_candidate_score_ =
            std::min(last_candidate_score_, candidate_score);
      }

      break;
    }
  }

  if (info.should_restart) {
    last_restart_iter_ = current_iter;
    double candidate_score = std::min(ComputeRestartScore(current_results),
                                      ComputeRestartScore(average_results));
    last_restart_score_ = candidate_score;
    last_candidate_score_ = last_restart_score_;
  }

  return info;
}
*/

RestartInfo RestartScheme::Check(int current_iter,
                                 const SolverResults& current_results,
                                 const SolverResults& average_results) {
  if (strategy_ != RestartStrategy::ADAPTIVE_RESTART) {
    if (strategy_ == RestartStrategy::FIXED_RESTART && 
        current_iter - last_restart_iter_ >= fixed_restart_interval_) {
        return RestartInfo(true, true);
    }
    return RestartInfo(false, false); // ✨ Corrected return
  }

  // Base Case: First check after a restart. Initialize state and exit.
  if (current_iter == last_restart_iter_) {
    primal_feas_last_restart_ = current_results.primal_feasibility;
    dual_feas_last_restart_ = current_results.dual_feasibility;
    duality_gap_last_restart_ = current_results.duality_gap;

    primal_feas_last_candidate_ = current_results.primal_feasibility;
    dual_feas_last_candidate_ = current_results.dual_feasibility;
    duality_gap_last_candidate_ = current_results.duality_gap;
    
    return {false, false};
  }

  // 1. Calculate scores and determine the best candidate
  double mu_current = compute_score(beta_, current_results.primal_feasibility,
                                    current_results.dual_feasibility, current_results.duality_gap);
  double mu_average = compute_score(beta_, average_results.primal_feasibility,
                                    average_results.dual_feasibility, average_results.duality_gap);
  debugPdlpRestartLog(debug_pdlp_log_file_, current_iter, mu_current,
                          mu_average);
  double candidate_score;
  bool restart_to_average = false;
  if (mu_current < mu_average) {
    candidate_score = mu_current;
    restart_to_average = false;
  } else {
    candidate_score = mu_average;
    restart_to_average = true;
  }

  // 2. Check the three restart conditions in order
  bool should_restart = false;
  PDHG_restart_choice restart_choice_for_logic = PDHG_NO_RESTART;
 
  if ((current_iter - last_restart_iter_) >= 0.36 * current_iter) {
    // Condition 1: Artificial Restart
    should_restart = true;
  } else {

    double mu_last_restart = compute_score(beta_, primal_feas_last_restart_,
                                           dual_feas_last_restart_, duality_gap_last_restart_);
    
    if (candidate_score < sufficient_decay_factor_ * mu_last_restart) {
      // Condition 2: Sufficient Decay
      should_restart = true;
    } else {
      double mu_last_candidate = compute_score(beta_, primal_feas_last_candidate_,
                                                 dual_feas_last_candidate_, duality_gap_last_candidate_);
      
      if (candidate_score < necessary_decay_factor_ * mu_last_restart && candidate_score > mu_last_candidate) {
        // Condition 3: Necessary Decay
        should_restart = true;
      }
    }
  }
  
  // 3. ALWAYS update the "last candidate" metrics for the next check
  if (restart_to_average) {
    primal_feas_last_candidate_ = average_results.primal_feasibility;
    dual_feas_last_candidate_ = average_results.dual_feasibility;
    duality_gap_last_candidate_ = average_results.duality_gap;
  } else {
    primal_feas_last_candidate_ = current_results.primal_feasibility;
    dual_feas_last_candidate_ = current_results.dual_feasibility;
    duality_gap_last_candidate_ = current_results.duality_gap;
  }
/*
  if (should_restart ) {
    if (restart_to_average){
      std::cout << "Last restart was iter " << last_restart_iter_ << ": average\n";
    } else {
      std::cout << "Last restart was iter " << last_restart_iter_ << ": current\n";
    }
  }
*/  
  return RestartInfo(should_restart, restart_to_average);
}