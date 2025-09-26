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
#ifndef PDHG_HPP
#define PDHG_HPP

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

// Forward declaration for a struct defined in the .cc file
struct StepSizeConfig;

class PDLPSolver {
 public:
  // --- Public API ---
  void setup(const HighsOptions& options, HighsTimer& timer);
  void passLp(const HighsLp* lp) { original_lp_ = lp; }
  void preprocessLp();
  void scaleProblem();
  void solve(std::vector<double>& x, std::vector<double>& y);
  void unscaleSolution(std::vector<double>& x, std::vector<double>& y);
  PostSolveRetcode postprocess(HighsSolution& solution);
  void logSummary();
  void solveReturn();

  // --- Getters ---
  TerminationStatus getTerminationCode() const { return results_.term_code; }
  int getIterationCount() const { return final_iter_count_; }
  int getnCol() const { return lp_.num_col_; }
  int getnRow() const { return lp_.num_row_; }

  // --- Debugging ---
  FILE* debug_pdlp_log_file_ = nullptr;
  DebugPdlpData debug_pdlp_data_;

 private:
  // --- Core Algorithm Logic ---
  void Initialize();
  bool CheckConvergence(const int iter, const std::vector<double>& x,
                        const std::vector<double>& y,
                        const std::vector<double>& ax_vector,
                        const std::vector<double>& aty_vector, double epsilon,
                        SolverResults& results, const char* type);
  void UpdateAverageIterates(const std::vector<double>& x,
                             const std::vector<double>& y,
                             const PrimalDualParams& params, int inner_iter);
  void ComputeAverageIterate(std::vector<double>& ax_avg,
                                       std::vector<double>& aty_avg);
  double PowerMethod();

  // --- Step Update Methods (previously in Step) ---
  StepSizeConfig InitializeStepSizesPowerMethod(double op_norm_sq);
  std::vector<double> UpdateX();
  std::vector<double> UpdateY();
  void UpdateIteratesFixed();
  void UpdateIteratesAdaptive(int& step_size_iter_count);
  bool UpdateIteratesMalitskyPock(bool first_malitsky_iteration);

  // --- Step Size Helper Methods (previously in PdlpStep) ---
  bool CheckNumericalStability(const std::vector<double>& delta_x,
                               const std::vector<double>& delta_y);
  double ComputeMovement(const std::vector<double>& delta_primal,
                         const std::vector<double>& delta_dual);
  double ComputeNonlinearity(const std::vector<double>& delta_primal,
                             const std::vector<double>& delta_aty);

  // --- Feasibility, Duality, and KKT Checks ---
  std::vector<double> ComputeLambda(const std::vector<double>& y,
                                    const std::vector<double>& ATy_vector);
  double ComputeDualObjective(const std::vector<double>& y);
  double ComputePrimalFeasibility(
      const std::vector<double>& Ax_vector);
  void ComputeDualSlacks(const std::vector<double>& ATy_vector);
  double ComputeDualFeasibility(
      const std::vector<double>& ATy_vector);
  std::tuple<double, double, double, double, double> ComputeDualityGap(
      const std::vector<double>& x, const std::vector<double>& y,
      const std::vector<double>& lambda);
  void PDHG_Compute_Step_Size_Ratio(
      PrimalDualParams& working_params, const std::vector<double>& x_n_0,
      const std::vector<double>& y_n_0, const std::vector<double>& x_n_minus_1_0,
      const std::vector<double>& y_n_minus_1_0);

  // --- Problem Data and Parameters ---
  HighsLp lp_;
  const HighsLp* original_lp_;
  HighsLp unscaled_processed_lp_;
  PrimalDualParams params_;
  StepSizeConfig stepsize_;
  Logger logger_;
  HighsLogOptions log_options_;
  SolverResults results_;
  int original_num_col_;
  int num_eq_rows_;
  std::vector<bool> is_equality_row_; 
  std::vector<int> constraint_new_idx_;
  std::vector<ConstraintType> constraint_types_;

  // --- Solver State ---
  int final_iter_count_ = 0;
  std::vector<double> x_current_, y_current_;
  std::vector<double> x_next_, y_next_;
  std::vector<double> x_avg_, y_avg_;
  std::vector<double> x_sum_, y_sum_;
  double sum_weights_ = 0.0;
  double current_eta_ = 0.0;
  double ratio_last_two_step_sizes_ = 1.0;
  int num_rejected_steps_ = 0;
  std::vector<double> dSlackPos_;
  std::vector<double> dSlackNeg_;
  Timer total_timer;

  // --- Scaling ---
  Scaling scaling_;
  double unscaled_rhs_norm_ = 0.0;
  double unscaled_c_norm_ = 0.0;

  // --- Restarting ---
  RestartScheme restart_scheme_;
  std::vector<double> x_at_last_restart_;
  std::vector<double> y_at_last_restart_;

  // --- Caching for Matrix-Vector Products ---
  std::vector<double> Ax_cache_;
  std::vector<double> ATy_cache_;
  std::vector<double> Ax_next_, ATy_next_;
  std::vector<double> K_times_x_diff_;
};

#endif