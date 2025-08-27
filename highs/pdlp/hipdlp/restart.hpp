/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-07-09 14:54:26
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-11 15:19:54
 * @FilePath: /cupdlp-CPP/include/restart.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
// restart.hpp
#ifndef RESTART_HPP
#define RESTART_HPP

#include <vector>
#include <cmath>
#include <iostream>
#include "solver_results.hpp"
#include "scaling.hpp"

enum class Device {
    CPU,
    GPU 
};

enum class RestartStrategy {
    NO_RESTART,
    FIXED_RESTART,
    ADAPTIVE_RESTART
};

enum class StepSizeStrategy {
    FIXED,
    ADAPTIVE,
    MALITSKY_POCK
};

struct MalitskyPockParams {
    double step_size_interpolation = 0.5;  // Between 0 and 1
    double step_size_downscaling_factor = 0.7;
    double linesearch_contraction_factor = 0.99;
  void initialise();
};

struct AdaptiveLinesearchParams {
    double step_size_reduction_exponent = 0.3;
    double step_size_growth_exponent = 0.6;
  void initialise();
};

struct  PrimalDualParams {
    double eta;
    double omega;
    double tolerance;
    size_t max_iterations;
    Device device_type;
    double time_limit = 3600.0;

    // Restart parameters
    RestartStrategy restart_strategy;
    int fixed_restart_interval;

    bool use_halpern_restart = false;

    // Scaling parameters
    ScalingMethod scaling_method = ScalingMethod::NONE;
    bool use_ruiz_scaling = false;
    bool use_pc_scaling = false;
    bool use_l2_scaling = false;
    
    // Ruiz scaling parameters
    int ruiz_iterations = 10;
    double ruiz_norm = INFINITY;
    
    // Pock-Chambolle scaling parameters
    double pc_alpha = 1.0;

    // Step sizes strategy
    StepSizeStrategy step_size_strategy = StepSizeStrategy::FIXED;

    MalitskyPockParams malitsky_pock_params;
    AdaptiveLinesearchParams adaptive_linesearch_params;
  void initialise();
};

// Struct to communicate restart decisions
struct RestartInfo {
    bool should_restart = false;
    bool restart_to_average = false; // If true, restart to average; otherwise, to current
};

class RestartScheme {
public:
    RestartScheme() = default;

    void Initialize(const PrimalDualParams& params, const SolverResults& results);

    // Checks if a restart should be performed based on the chosen strategy
    RestartInfo Check(int current_iter, const SolverResults& current_results, const SolverResults& average_results);

    int GetLastRestartIter() const { return last_restart_iter_; }
private:
    // Computes a merit score for a given set of residuals
    double ComputeRestartScore(const SolverResults& results);

    RestartStrategy strategy_ = RestartStrategy::NO_RESTART;
    int fixed_restart_interval_ = 100;
    int last_restart_iter_ = 0;
    double beta_;

    // State for adaptive restart
    double last_restart_score_ = 1;
    double last_candidate_score_ = 1;
    double sufficient_decay_factor_ = 0.2;
    double necessary_decay_factor_ = 0.8;
};

#endif // RESTART_HPP
