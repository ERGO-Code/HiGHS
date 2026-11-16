#pragma once
#include <cuda_runtime.h>

struct GpuSolverMetrics {
    double primal_feas;
    double dual_feas;
    double primal_obj;
    double dual_obj;
    double duality_gap;
};

struct GpuRestartConfig {
    double sufficient_decay_factor; // e.g. 0.2
    double necessary_decay_factor;  // e.g. 0.8
    int fixed_restart_interval;
    int restart_strategy; // 0=Fixed, 1=Adaptive
};

struct GpuRestartState {
    int last_restart_iter;
    double beta; 
    
    // Last Restart Snapshot
    double primal_feas_last_restart;
    double dual_feas_last_restart;
    double duality_gap_last_restart;

    // Last Candidate Snapshot
    double primal_feas_last_candidate;
    double dual_feas_last_candidate;
    double duality_gap_last_candidate;
};

struct GpuRestartInfo {
    int should_restart;      // 0 or 1
    int restart_to_average;  // 0 or 1
};