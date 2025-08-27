/*
 * @Author: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @Date: 2025-07-11 22:24:59
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-07-14 12:04:26
 * @FilePath: /cupdlp-CPP/include/solver_results.hpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef SOLVER_RESULTS_HPP
#define SOLVER_RESULTS_HPP

#include <vector>

enum class TerminationStatus {
    OPTIMAL,
    INFEASIBLE,
    UNBOUNDED,
    TIMEOUT,
    FEASIBLE,
    WARNING,
    ERROR
};

enum class TerminationIterate {
    LAST_ITERATE,
    AVERAGE_ITERATE
};

struct SolverResults {
    TerminationStatus term_code = TerminationStatus::TIMEOUT;
    TerminationIterate term_iterate = TerminationIterate::LAST_ITERATE;

    double primal_obj = 0.0;
    double dual_obj = 0.0;
    double duality_gap = 0.0;
    double complementarity = 0.0;
    double primal_feasibility = 0.0;
    double dual_feasibility = 0.0;
    double relative_obj_gap = 0.0;

    // --- Averaged Metrics ---
    double primal_obj_average = 0.0;
    double dual_obj_average = 0.0;
    double duality_gap_average = 0.0;
    double primal_feasibility_average = 0.0;
    double dual_feasibility_average = 0.0;

    // --- Residual Vectors ---
    std::vector<double> primal_residual;
    std::vector<double> dual_residual;
    std::vector<double> primal_residual_average;
    std::vector<double> dual_residual_average;

    // --- Infeasibility Detection Metrics ---
    TerminationStatus primal_term_code = TerminationStatus::FEASIBLE;
    TerminationStatus dual_term_code = TerminationStatus::FEASIBLE;
    double primal_infeasibility_obj = 0.0;
    double dual_infeasibility_obj = 0.0;
    double primal_infeasibility_res = 0.0;
    double dual_infeasibility_res = 0.0;

    // --- Restart Statistics ---
    double primal_feasibility_last_restart = 0.0;
    double dual_feasibility_last_restart = 0.0;
    double duality_gap_last_restart = 0.0;
};

#endif // SOLVER_RESULTS_HPP