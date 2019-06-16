/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsLp.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_LP_H_
#define LP_DATA_HIGHS_LP_H_

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "HConfig.h"
#include "lp_data/HConst.h" // For HiGHS strategy options
#include "simplex/SimplexConst.h" // For simplex strategy options

enum class LpAction {
    DUALISE = 0,
    PERMUTE,
    SCALE,
    NEW_COSTS,
    NEW_BOUNDS,
    NEW_BASIS,
    NEW_COLS,
    NEW_ROWS,
    DEL_COLS,
    DEL_ROWS,
    DEL_ROWS_BASIS_OK
    };

class HighsLp {
 public:
  // Model data
  int numCol_ = 0;
  int numRow_ = 0;
  int numInt_ = 0;
  int nnz_ = 0;

  std::vector<int> Astart_;
  std::vector<int> Aindex_;
  std::vector<double> Avalue_;
  std::vector<double> colCost_;
  std::vector<double> colLower_;
  std::vector<double> colUpper_;
  std::vector<double> rowLower_;
  std::vector<double> rowUpper_;

  // sense 1 = minimize, -1 = maximize
  int sense_ = 1;
  double offset_ = 0;

  std::string model_name_ = "";

  std::vector<std::string> row_names_;
  std::vector<std::string> col_names_;

  std::vector<int> integrality_;

  bool operator==(const HighsLp& lp) {
    if (numCol_ != lp.numCol_ || numRow_ != lp.numRow_ || nnz_ != lp.nnz_ ||
        sense_ != lp.sense_ || offset_ != lp.offset_ ||
        model_name_ != lp.model_name_)
      return false;

    if (row_names_ != lp.row_names_ || col_names_ != lp.col_names_)
      return false;

    if (colCost_ != lp.colCost_)
      return false;

    if (colUpper_ != lp.colUpper_ || colLower_ != lp.colLower_ ||
        rowUpper_ != lp.rowUpper_ || rowLower_ != lp.rowLower_)
      return false;

    if (Astart_ != lp.Astart_ || Aindex_ != lp.Aindex_ ||
        Avalue_ != lp.Avalue_)
      return false;

    return true;
  }
};

// Cost, column and row scaling factors
struct HighsScale {
  double cost_;
  std::vector<double> col_;
  std::vector<double> row_;
};

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
  bool valid = false;
  bool is_dualised = false;
  bool is_permuted = false;
  bool is_scaled = false;
  bool has_basis = false; // The LP has a valid simplex basis
  bool has_matrix_col_wise = false; // The LP has a column-wise constraint matrix
  bool has_matrix_row_wise = false; // The LP has a row-wise constraint matrix
  bool has_factor_arrays = false; // Has the arrays for the representation of B^{-1}
  bool has_dual_steepest_edge_weights = false; // The DSE weights are known
  bool has_nonbasic_dual_values = false; // The nonbasic dual values are known
  bool has_basic_primal_values = false;// The basic primal values are known
  bool has_invert = false; // The representation of B^{-1} corresponds to the current basis
  bool has_fresh_invert = false; // The representation of B^{-1} corresponds to the current basis and is fresh
  bool has_fresh_rebuild = false; // The data are fresh from rebuild
  bool has_dual_objective_value = false; // The dual objective function value is known
  bool has_primal_objective_value = false; // The dual objective function value is known
  SimplexSolutionStatus solution_status = SimplexSolutionStatus::UNSET; // The solution status is UNSET
};

struct HighsSimplexInfo {
  // Simplex information regarding primal and dual solution, objective
  // and iteration counts for this Highs Model Object. This is
  // information which should be retained from one run to the next in
  // order to provide hot starts.
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
  // workShift: WTF
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

  // Values of iClock for simplex timing clocks
  std::vector<int> clock_;
  //
  // Options from HighsOptions for the simplex solver
  double highs_run_time_limit;
  SimplexStrategy simplex_strategy;
  SimplexDualEdgeWeightStrategy dual_edge_weight_strategy;
  SimplexPrimalEdgeWeightStrategy primal_edge_weight_strategy;
  SimplexPriceStrategy price_strategy;

  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  bool perturb_costs;
  int update_limit;
  int iteration_limit;
  double dual_objective_value_upper_bound;
  
  // Internal options - can't be changed externally

  bool analyseLpSolution;
#ifdef HiGHSDEV
  // Options for reporting timing
  bool report_simplex_inner_clock;
  bool report_simplex_outer_clock;
  bool report_simplex_phases_clock;
  // Option for analysing the LP simplex iterations, INVERT time and rebuild time
  bool analyseLp;
  bool analyseSimplexIterations;
  bool analyse_invert_time;
  bool analyseRebuildTime;
#endif
  // Simplex runtime information
  int costs_perturbed = 0;
  // Cumulative iteration count - updated in simplex solvers
  int iteration_count = 0;
  // Records of cumulative iteration counts - updated at the end of a phase
  int dual_phase1_iteration_count = 0;
  int dual_phase2_iteration_count = 0;
  int primal_phase1_iteration_count = 0;
  int primal_phase2_iteration_count = 0;

  // Cutoff for PAMI
  double pami_cutoff = 0.95;

  // Info on PAMI iterations
  int multi_iteration = 0;

  // Number of UPDATE operations performed - should be zeroed when INVERT is performed
  int update_count;
  // Value of dual objective - only set when computed from scratch in dual rebuild()
  double dual_objective_value;
  // Value of primal objective - only set when computed from scratch in primal rebuild()
  double primal_objective_value;


  // Value of dual objective that is updated in dual simplex solver
  double updated_dual_objective_value;
  // Value of primal objective that is updated in primal simplex solver
  double updated_primal_objective_value;
  // Number of logical variables in the basis 
  int num_basic_logicals;
  // Number/max/sum of primal and dual infeasibilities
  int num_primal_infeasibilities;
  double max_primal_infeasibility;
  double sum_primal_infeasibilities;
  int num_dual_infeasibilities;
  double max_dual_infeasibility;
  double sum_dual_infeasibilities;

#ifdef HiGHSDEV
  // Analysis of INVERT
  int total_inverts;
  double total_invert_time;
#endif

  /*
#ifdef HiGHSDEV
  // Move this to Simplex class once it's created
  vector<int> historyColumnIn;
  vector<int> historyColumnOut;
  vector<double> historyAlpha;
#endif
  */

};

struct HighsSolution {
  std::vector<double> col_value;
  std::vector<double> col_dual;
  std::vector<double> row_value;
  std::vector<double> row_dual;
};

// To be the basis representation given back to the user. Values of
// HighsBasisStatus are defined in HConst.h
struct HighsBasis {
  bool valid_ = false;
  std::vector<HighsBasisStatus> col_status;
  std::vector<HighsBasisStatus> row_status;
};

struct HighsRanging {
  std::vector<double> colCostRangeUpValue_;
  std::vector<double> colCostRangeUpObjective_;
  std::vector<int>    colCostRangeUpInCol_;
  std::vector<int>    colCostRangeUpOutCol_;
  std::vector<double> colCostRangeDnValue_;
  std::vector<double> colCostRangeDnObjective_;
  std::vector<int>    colCostRangeDnInCol_;
  std::vector<int>    colCostRangeDnOutCol_;
  std::vector<double> rowBoundRangeUpValue_;
  std::vector<double> rowBoundRangeUpObjective_;
  std::vector<int>    rowBoundRangeUpInCol_;
  std::vector<int>    rowBoundRangeUpOutCol_;
  std::vector<double> rowBoundRangeDnValue_;
  std::vector<double> rowBoundRangeDnObjective_;
  std::vector<int>    rowBoundRangeDnInCol_;
  std::vector<int>    rowBoundRangeDnOutCol_;
};

// Make sure the dimensions of solution are the same as numRow_ and numCol_.
bool isSolutionConsistent(const HighsLp& lp, const HighsSolution& solution);

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
//void checkStatus(HighsStatus status);


#endif
