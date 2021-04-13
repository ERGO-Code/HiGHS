/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2021 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/SimplexStruct.h
 * @brief Structs for HiGHS simplex solvers
 */
#ifndef SIMPLEX_SIMPLEXSTRUCT_H_
#define SIMPLEX_SIMPLEXSTRUCT_H_

#include "HConfig.h"
#include "simplex/SimplexConst.h"

struct SimplexBasis {
  // The basis for the simplex method consists of basicIndex,
  // nonbasicFlag and nonbasicMove. If HighsSimplexLpStatus has_basis
  // is true then it is assumed that basicIndex_ and nonbasicFlag_ are
  // self-consistent and correpond to the dimensions of an associated
  // HighsLp, but the basis matrix B is not necessarily nonsingular.
  std::vector<HighsInt> basicIndex_;
  std::vector<int8_t> nonbasicFlag_;
  std::vector<int8_t> nonbasicMove_;
};

struct HighsSimplexLpStatus {
  // Status of LP solved by the simplex method and its data
  bool initialised = false;
  bool valid = false;
  bool scaling_tried = false;
  bool has_basis = false;   // The simplex LP has a valid simplex basis
  bool has_matrix = false;  // The HMatrix matrices are valid
  bool has_factor_arrays =
      false;  // Has the arrays for the representation of B^{-1}
  bool has_dual_steepest_edge_weights = false;  // The DSE weights are known
  bool has_nonbasic_dual_values = false;  // The nonbasic dual values are known
  bool has_basic_primal_values = false;   // The basic primal values are known
  bool has_invert =
      false;  // The representation of B^{-1} corresponds to the current basis
  bool has_fresh_invert = false;  // The representation of B^{-1} corresponds to
                                  // the current basis and is fresh
  bool has_fresh_rebuild = false;  // The data are fresh from rebuild
  bool has_dual_objective_value =
      false;  // The dual objective function value is known
  bool has_primal_objective_value =
      false;                    // The dual objective function value is known
  bool has_dual_ray = false;    // A dual unbounded ray is known
  bool has_primal_ray = false;  // A primal unbounded ray is known
};

struct HighsSimplexInfo {
  // Simplex information regarding primal solution, dual solution and
  // objective for this Highs Model Object. This is information which
  // should be retained from one run to the next in order to provide
  // hot starts.
  //
  // Part of working model which are assigned and populated as much as
  // possible when a model is being defined

  // workCost: Originally just costs from the model but, in solve(), may
  // be perturbed or set to alternative values in Phase I??
  //
  // workDual: Values of the dual variables corresponding to
  // workCost. Latter not known until solve() is called since B^{-1}
  // is required to compute them. Knowledge of them is indicated by
  // has_nonbasic_dual_values
  //
  // workShift: Values added to workCost in order that workDual
  // remains feasible, thereby remaining dual feasible in phase 2
  //
  std::vector<double> workCost_;
  std::vector<double> workDual_;
  std::vector<double> workShift_;

  // workLower/workUpper: Originally just lower (upper) bounds from
  // the model but, in solve(), may be perturbed or set to
  // alternative values in Phase I??
  //
  // workRange: Distance between lower and upper bounds
  //
  // workValue: Values of the nonbasic variables corresponding to
  // workLower/workUpper and the basis. Always known.
  //
  std::vector<double> workLower_;
  std::vector<double> workUpper_;
  std::vector<double> workRange_;
  std::vector<double> workValue_;
  std::vector<double> workLowerShift_;
  std::vector<double> workUpperShift_;
  //
  // baseLower/baseUpper/baseValue: Lower and upper bounds on the
  // basic variables and their values. Latter not known until solve()
  // is called since B^{-1} is required to compute them. Knowledge of
  // them is indicated by has_basic_primal_values
  //
  std::vector<double> baseLower_;
  std::vector<double> baseUpper_;
  std::vector<double> baseValue_;
  //
  // Vectors of random reals for column cost perturbation, a random
  // permutation of all indices for CHUZR and a random permutation of
  // column indices for permuting the columns
  std::vector<double> numTotRandomValue_;
  std::vector<HighsInt> numTotPermutation_;
  std::vector<HighsInt> numColPermutation_;

  std::vector<HighsInt> devex_index_;

  // Data for backtracking in the event of a singular basis
  HighsInt phase1_backtracking_test_done = false;
  HighsInt phase2_backtracking_test_done = false;
  bool backtracking_ = false;
  bool valid_backtracking_basis_ = false;
  SimplexBasis backtracking_basis_;
  HighsInt backtracking_basis_costs_perturbed_;
  HighsInt backtracking_basis_bounds_perturbed_;
  std::vector<double> backtracking_basis_workShift_;
  std::vector<double> backtracking_basis_workLowerShift_;
  std::vector<double> backtracking_basis_workUpperShift_;
  std::vector<double> backtracking_basis_edge_weights_;

  // Dual and primal ray vectors
  HighsInt dual_ray_row_;
  HighsInt dual_ray_sign_;
  HighsInt primal_ray_col_;
  HighsInt primal_ray_sign_;

  // Options from HighsOptions for the simplex solver
  HighsInt simplex_strategy;
  HighsInt dual_edge_weight_strategy;
  HighsInt primal_edge_weight_strategy;
  HighsInt price_strategy;

  double dual_simplex_cost_perturbation_multiplier;
  double primal_simplex_phase1_cost_perturbation_multiplier = 1;
  double primal_simplex_bound_perturbation_multiplier;
  double factor_pivot_threshold;
  HighsInt update_limit;

  // Simplex control parameters from HSA
  HighsInt control_iteration_count0;
  double col_aq_density;
  double row_ep_density;
  double row_ap_density;
  double row_DSE_density;
  double col_basic_feasibility_change_density;
  double row_basic_feasibility_change_density;
  double col_BFRT_density;
  double primal_col_density;
  double dual_col_density;
  // For control of switch from DSE to Devex in dual simplex
  bool allow_dual_steepest_edge_to_devex_switch;
  double dual_steepest_edge_weight_log_error_threshold;
  double costly_DSE_frequency;
  HighsInt num_costly_DSE_iteration;
  double average_log_low_DSE_weight_error;
  double average_log_high_DSE_weight_error;
  // Needed globally??

  // Internal options - can't be changed externally
  bool run_quiet = false;
  bool store_squared_primal_infeasibility = false;

  //  bool analyse_lp_solution = true;
  // Options for reporting timing
  bool report_simplex_inner_clock = false;
  bool report_simplex_outer_clock = false;
  bool report_simplex_phases_clock = false;
  bool report_HFactor_clock = false;
  // Option for analysing the LP simplex iterations, INVERT time and rebuild
  // time
  bool analyse_lp = false;
  bool analyse_iterations = false;
  bool analyse_invert_form = false;
  //  bool analyse_invert_condition = false;
  //  bool analyse_invert_time = false;
  //  bool analyse_rebuild_time = false;

  // Simplex runtime information
  bool allow_cost_perturbation = true;
  bool costs_perturbed = false;
  bool bounds_perturbed = false;

  HighsInt num_primal_infeasibility = -1;
  double max_primal_infeasibility;
  double sum_primal_infeasibility;
  HighsInt num_dual_infeasibility = -1;
  double max_dual_infeasibility;
  double sum_dual_infeasibility;

  // Records of cumulative iteration counts - updated at the end of a phase
  HighsInt dual_phase1_iteration_count = 0;
  HighsInt dual_phase2_iteration_count = 0;
  HighsInt primal_phase1_iteration_count = 0;
  HighsInt primal_phase2_iteration_count = 0;
  HighsInt primal_bound_swap = 0;

  HighsInt min_threads = 1;
  HighsInt num_threads = 1;
  HighsInt max_threads = kHighsThreadLimit;

  // Cutoff for PAMI
  double pami_cutoff = 0.95;

  // Info on PAMI iterations
  HighsInt multi_iteration = 0;

  // Number of UPDATE operations performed - should be zeroed when INVERT is
  // performed
  HighsInt update_count;
  // Value of dual objective - only set when computed from scratch in dual
  // rebuild()
  double dual_objective_value;
  // Value of primal objective - only set when computed from scratch in primal
  // rebuild()
  double primal_objective_value;

  // Value of dual objective that is updated in dual simplex solver
  double updated_dual_objective_value;
  // Value of primal objective that is updated in primal simplex solver
  double updated_primal_objective_value;
  // Number of logical variables in the basis
  HighsInt num_basic_logicals;
};

#endif /* SIMPLEX_SIMPLEXSTRUCT_H_ */
