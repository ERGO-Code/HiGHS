/*
 * @Author: Yanyu000 earthazyy@hotmail.com
 * @Date: 2025-05-12 12:43:12
 * @LastEditors: Zhou Yanyu（周妍妤） 47125824+Yanyu000@users.noreply.github.com
 * @LastEditTime: 2025-08-21 16:27:35
 * @FilePath: /cupdlp-CPP/include/pdhg.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef PDHG_HPP
#define PDHG_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <tuple>
#include <memory>

#include "Highs.h"

#include "linalg.hpp"
#include "solver_results.hpp"
#include "restart.hpp"
#include "scaling.hpp"
#include "step.hpp"
#include "logger.hpp"
//#include "highs_interface.hpp"

// --- Define Macros ---
// Enable or disable GPU usage
#define USE_GPU 1
// Debug mode
#define DEBUG_MODE 1
enum ConstraintType { EQ, GEQ, LEQ, BOUND, FREE };

// --- Classes ---
class PDLPSolver {
    public:
        PDLPSolver(Logger& logger);
        
        // Main solve method
        void Solve(HighsLp & lp, const PrimalDualParams& params, std::vector<double>& x, std::vector<double>& y);

        const SolverResults& GetResults() const {
            return results_;
        }

        void SetScalingParams(const ScalingParams& scaling_params) {
            scaling_params_ = scaling_params;
        }

        int GetIterationCount() const;

    private:
        //Problem data
        HighsLp lp_;
        PrimalDualParams params_;
  //        HighsInterface highs_interface_;
        Logger& logger_;
        int final_iter_count_ = 0;
        int original_num_col_;
        std::vector<int> constraint_new_idx_; // stores permutation for rows
        std::vector<ConstraintType> constraint_types_; // stores type of each original row

        //Iterates and state
        std::vector<double>* x_ = nullptr;
        std::vector<double>* y_ = nullptr;
        std::vector<double> x_current_;
        std::vector<double> y_current_;
        // Restart State
        std::vector<double> x_avg_, y_avg_;
        std::vector<double> x_sum_, y_sum_;
        double sum_weights_ = 0.0;
        std::vector<double> x_outer_start_, y_outer_start_;
        std::vector<double> x_prev_outer_start_, y_prev_outer_start_;
        RestartScheme restart_scheme_;

        //  Halpern restart
        std::vector<double> x_initial_;
        std::vector<double> y_initial_;

        std::vector<double> dSlackPos_;
        std::vector<double> dSlackNeg_;

        // Scaling
        ScalingParams scaling_params_;
        Scaling scaling_;

        //HighsStatus TransformGxLeqToGeq(HighsLp& lp);
        void PreprocessLp(const HighsLp& original_lp, HighsLp& processed_lp);
        void Postsolve(const HighsLp& original_lp, HighsLp& processed_lp,
                           const std::vector<double>& x_processed, 
                           const std::vector<double>& y_processed,
                           HighsSolution& solution);
        
        // Helper functions
        void Initialize(const HighsLp & lp, std::vector<double>& x, std::vector<double>& y);
        bool RunPresolve(const HighsLp& original_lp, Highs& highs, HighsLp& presolved_lp);
        SolverResults results_;
        TerminationStatus SolvePresolvedProblem(HighsLp& presolved_lp, const PrimalDualParams& params, std::vector<double>& x, std::vector<double>& y);
        void PostsolveAndFinalize(const std::vector<double>& presolved_x, const std::vector<double>& presolved_y, std::vector<double>& final_x, std::vector<double>& final_y);

        // Check convergence 
        double ComputeWeightedNorm(const std::vector<double>& x1, const std::vector<double>& y1,
                                   const std::vector<double>& x2, const std::vector<double>& y2, double omega);
        std::vector<double> ComputeLambda(const HighsLp& lp, const std::vector<double>& y, const std::vector<double>& ATy_vector);
        std::pair<double,double> ComputePrimalFeasibility(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& Ax_vector);
        void ComputeDualSlacks(const HighsLp& lp, const std::vector<double>& ATy_vector);
        std::pair<double, double> ComputeDualFeasibility(const HighsLp& lp, const std::vector<double>& ATy_vector);
        std::tuple<double, double, double, double, double> ComputeDualityGap(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& lambda);
        double ComputeKKTError(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& lambda, double omega);
        double ComputeNormalizedDualityGap(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& x_ref, const std::vector<double>& y_ref, double omega, const std::vector<double>& ax, const std::vector<double>& aty);
        bool CheckConvergence(const HighsLp& lp, const std::vector<double>& x, const std::vector<double>& y,const std::vector<double>& ax_vector, const std::vector<double>& aty_vector, double epsilon, SolverResults& results);

        // GPU specific functions

        //step size 
        double current_eta_;
        //void UpdateIteratesAdaptive(HighsLp& lp, const PrimalDualParams& params, std::vector<double>& x, std::vector<double>& y);
        HighsStatus PowerMethod(HighsLp & lp, double& lambda);
        double ratio_last_two_step_sizes_ = 1.0; // state for Malitsky-Pock adaptive step size
        int num_rejected_steps_ = 0; // state for adaptive linesearch
        static constexpr double kDivergentMovement = 1e10;

        // Primal weight update
        std::vector<double> x_at_last_restart_;
        std::vector<double> y_at_last_restart_;

        void PDHG_Compute_Step_Size_Ratio(
                PrimalDualParams& working_params,
                const std::vector<double>& x_n_0, // Corresponds to z^{n,0}
                const std::vector<double>& y_n_0,
                const std::vector<double>& x_n_minus_1_0, // Corresponds to z^{n-1,0}
                const std::vector<double>& y_n_minus_1_0
            );

        //restart
        bool CheckRestartCondition(const HighsLp& lp, const PrimalDualParams& params, int inner_iter, std::vector<double>& x_cand, std::vector<double>& y_cand);
        void PerformRestart(const std::vector<double>& x_restart, const std::vector<double>& y_restart, int inner_iter, const PrimalDualParams& params, const HighsLp& lp);
        void UpdateAverageIterates(const std::vector<double>& x, const std::vector<double>& y, const PrimalDualParams& params, int inner_iter);

        // Cache Matrix-Vector products
        std::vector<double> Ax_cache_;
        std::vector<double> ATy_cache_;
        std::vector<double> K_times_x_diff_;

        double mu_candidate_ = 0.0;
        double mu_prev_candidate_ = 0.0;
        double mu_outer_start_ = 0.0;

        int outer_iter_ = 0;
        int last_outer_loop_iter_count_ = 0;
};

#endif 
