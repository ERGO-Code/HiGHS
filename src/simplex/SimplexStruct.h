/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/SimplexStruct.h
 * @brief Structs for HiGHS simplex solvers
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
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
  std::vector<int> basicIndex_;
  std::vector<int> nonbasicFlag_;
  std::vector<int> nonbasicMove_;
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
  SimplexSolutionStatus solution_status =
      SimplexSolutionStatus::UNSET;  // The solution status is UNSET
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
  std::vector<int> numTotPermutation_;
  std::vector<int> numColPermutation_;

  std::vector<int> devex_index_;

  // Data for backtracking in the event of a singular basis
  int phase1_backtracking_test_done = false;
  int phase2_backtracking_test_done = false;
  bool backtracking_ = false;
  bool valid_backtracking_basis_ = false;
  SimplexBasis backtracking_basis_;
  int backtracking_basis_costs_perturbed_;
  int backtracking_basis_bounds_perturbed_;
  std::vector<double> backtracking_basis_workShift_;
  std::vector<double> backtracking_basis_workLowerShift_;
  std::vector<double> backtracking_basis_workUpperShift_;
  std::vector<double> backtracking_basis_edge_weights_;

  // Dual and primal ray vectors
  int dual_ray_row_;
  int dual_ray_sign_;
  int primal_ray_col_;
  int primal_ray_sign_;

  // Options from HighsOptions for the simplex solver
  int simplex_strategy;
  int dual_edge_weight_strategy;
  int primal_edge_weight_strategy;
  int price_strategy;

  double dual_simplex_cost_perturbation_multiplier;
  double primal_simplex_phase1_cost_perturbation_multiplier = 1;
  double primal_simplex_bound_perturbation_multiplier;
  double factor_pivot_threshold;
  int update_limit;

  // Simplex control parameters from HSA

  // Internal options - can't be changed externally
  bool run_quiet = false;
  bool store_squared_primal_infeasibility = false;
#ifndef HiGHSDEV
  bool analyse_lp_solution = false;  // true;//
#else
  bool analyse_lp_solution = true;
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
  bool analyse_invert_condition = false;
  bool analyse_invert_time = false;
  bool analyse_rebuild_time = false;
#endif
  // Simplex runtime information
  bool allow_cost_perturbation = true;
  bool costs_perturbed = false;
  bool bounds_perturbed = false;

  int num_primal_infeasibilities = -1;
  double max_primal_infeasibility;
  double sum_primal_infeasibilities;
  int num_dual_infeasibilities = -1;
  double max_dual_infeasibility;
  double sum_dual_infeasibilities;

  // Records of cumulative iteration counts - updated at the end of a phase
  int dual_phase1_iteration_count = 0;
  int dual_phase2_iteration_count = 0;
  int primal_phase1_iteration_count = 0;
  int primal_phase2_iteration_count = 0;
  int primal_bound_swap = 0;

  int min_threads = 1;
  int num_threads = 1;
  int max_threads = HIGHS_THREAD_LIMIT;

  // Cutoff for PAMI
  double pami_cutoff = 0.95;

  // Info on PAMI iterations
  int multi_iteration = 0;

  // Number of UPDATE operations performed - should be zeroed when INVERT is
  // performed
  int update_count;
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
  int num_basic_logicals;
};

#endif /* SIMPLEX_SIMPLEXSTRUCT_H_ */
