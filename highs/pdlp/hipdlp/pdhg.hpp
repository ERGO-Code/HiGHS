/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file pdlp/hipdlp/pdhg.hpp
 * @brief
 */
#ifndef PDLP_HIPDLP_PDHG_HPP
#define PDLP_HIPDLP_PDHG_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Highs.h"
#include "linalg.hpp"
#include "logger.hpp"
#include "restart.hpp"
#include "scaling.hpp"
#include "solver_results.hpp"
#include "step.hpp"

// --- Define Macros ---
// Enable or disable GPU usage
#define USE_GPU 1
// Debug mode
#define DEBUG_MODE 1
// --- Classes ---
class PDLPSolver {
 public:
  void setup(const HighsOptions& options, HighsTimer& timer);
  void preprocessLp();
  void scaleProblem();
  void solve(std::vector<double>& x, std::vector<double>& y);
  void unscaleSolution(std::vector<double>& x, std::vector<double>& y);
  PostSolveRetcode postprocess(HighsSolution& solution);
  void setSolution(const std::vector<double>& col_value,
                   const std::vector<double>& row_dual);
  void getSolution(std::vector<double>& col_value,
                   std::vector<double>& row_dual);

  void passLp(const HighsLp* lp) { original_lp_ = lp; }
  TerminationStatus getTerminationCode() const { return results_.term_code; }
  int getIterationCount() const { return final_iter_count_; }
  void logSummary();

  void solveReturn();
  FILE* debug_pdlp_log_file_ = nullptr;
  DebugPdlpData debug_pdlp_data_;

 private:
  // Problem data
  HighsLp lp_;                  // The problem to solve
  const HighsLp* original_lp_;  // The original problem (for postsolve)
  HighsLp unscaled_processed_lp_;
  PrimalDualParams params_;
  Logger logger_;
  PdlpStep step_;

  int final_iter_count_ = 0;
  int original_num_col_;
  int original_num_row_;
  int num_eq_rows_;  // Number of equality-like rows (EQ, BOUND, FREE)
  std::vector<int> constraint_orig_idx_;
  std::vector<int> constraint_new_idx_;  // stores permutation for rows
  std::vector<ConstraintType>
      constraint_types_;  // stores type of each original row

  // Solver state using PdlpIterate
  PdlpIterate current_;        // z^k
  PdlpIterate next_;           // z^{k+1}
  PdlpIterate average_;        // Running average
  PdlpIterate restart_point_;  // Point to restart from

  // Iterates and state
  std::vector<double>* x_ = nullptr;
  std::vector<double>* y_ = nullptr;
  std::vector<double> x_current_;
  std::vector<double> y_current_;
  // Scaling
  Scaling scaling_;
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

  Timer total_timer;

  // Helper functions
  void Initialize(const HighsLp& lp, std::vector<double>& x,
                  std::vector<double>& y);
  SolverResults results_;

  void PostsolveAndFinalize(const std::vector<double>& presolved_x,
                            const std::vector<double>& presolved_y,
                            std::vector<double>& final_x,
                            std::vector<double>& final_y);

  // Check convergence
  double ComputeWeightedNorm(const std::vector<double>& x1,
                             const std::vector<double>& y1,
                             const std::vector<double>& x2,
                             const std::vector<double>& y2, double omega);
  std::vector<double> ComputeLambda(const std::vector<double>& y,
                                    const std::vector<double>& ATy_vector);
  double ComputeDualObjective(const std::vector<double>& y);
  std::pair<double, double> ComputePrimalFeasibility(
      const std::vector<double>& x, const std::vector<double>& Ax_vector);
  void ComputeDualSlacks(const std::vector<double>& ATy_vector);
  std::pair<double, double> ComputeDualFeasibility(
      const std::vector<double>& ATy_vector);
  std::tuple<double, double, double, double, double> ComputeDualityGap(
      const std::vector<double>& x, const std::vector<double>& y,
      const std::vector<double>& lambda);
  double ComputeKKTError(const std::vector<double>& x,
                         const std::vector<double>& y,
                         const std::vector<double>& lambda, double omega);
  bool CheckConvergence(const int iter, const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& ax_vector,
                        const std::vector<double>& aty_vector, double epsilon,
                        SolverResults& results, const char* type);

  // GPU specific functions

  // step size
  double current_eta_;
  // void UpdateIteratesAdaptive(HighsLp& lp, const PrimalDualParams& params,
  // std::vector<double>& x, std::vector<double>& y);
  double PowerMethod();

  double ratio_last_two_step_sizes_ =
      1.0;                      // state for Malitsky-Pock adaptive step size
  int num_rejected_steps_ = 0;  // state for adaptive linesearch
  static constexpr double kDivergentMovement = 1e10;

  // Primal weight update
  std::vector<double> x_at_last_restart_;
  std::vector<double> y_at_last_restart_;

  void PDHG_Compute_Step_Size_Ratio(
      PrimalDualParams& working_params,
      const std::vector<double>& x_n_0,  // Corresponds to z^{n,0}
      const std::vector<double>& y_n_0,
      const std::vector<double>& x_n_minus_1_0,  // Corresponds to z^{n-1,0}
      const std::vector<double>& y_n_minus_1_0);

  // restart
  bool CheckRestartCondition(const PrimalDualParams& params, int inner_iter,
                             std::vector<double>& x_cand,
                             std::vector<double>& y_cand);
  void PerformRestart(const std::vector<double>& x_restart,
                      const std::vector<double>& y_restart, int inner_iter,
                      const PrimalDualParams& params, const HighsLp& lp);
  void UpdateAverageIterates(const std::vector<double>& x,
                             const std::vector<double>& y,
                             const PrimalDualParams& params, int inner_iter);

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
