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

#include "HConfig.h"
#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "HConst.h" // For HiGHS strategy options
#include "SimplexConst.h" // For simplex strategy options

// The free parser also reads fixed format MPS files but the fixed
// parser does not read free mps files.
enum class HighsMpsParserType { free, fixed };

/** SCIP/HiGHS Objective sense */
enum objSense
{
  OBJSENSE_MINIMIZE = 1,
  OBJSENSE_MAXIMIZE = -1
};

// For now, but later change so HiGHS properties are string based so that new
// options (for debug and testing too) can be added easily. The options below
// are just what has been used to parse options from argv.
// todo: when creating the new options don't forget underscores for class
// variables but no underscores for struct
struct HighsOptions {
  std::string filenames = "";

  // Options passed through the command line

  ParallelOption parallel_option = ParallelOption::DEFAULT;
  PresolveOption presolve_option = PresolveOption::DEFAULT;
  CrashOption crash_option = CrashOption::DEFAULT;
  SimplexOption simplex_option = SimplexOption::DEFAULT;
  double highs_run_time_limit = HIGHS_RUN_TIME_LIMIT_DEFAULT;
  double infinite_bound = INFINITE_BOUND_DEFAULT;
  double small_matrix_value = SMALL_MATRIX_VALUE_DEFAULT;

  bool pami = 0;
  bool sip = 0;
  bool scip = 0;
  SimplexStrategy simplex_strategy = SimplexStrategy::DEFAULT;
  SimplexCrashStrategy simplex_crash_strategy = SimplexCrashStrategy::DEFAULT;
  HighsMpsParserType parser_type = HighsMpsParserType::free;

  SimplexDualEdgeWeightStrategy simplex_dual_edge_weight_strategy = SimplexDualEdgeWeightStrategy::DEFAULT;
  SimplexPriceStrategy simplex_price_strategy = SimplexPriceStrategy::DEFAULT;

  // Options not passed through the command line

  // Options for HighsPrintMessage and HighsLogMessage
  FILE* logfile = stdout;
  FILE* output = stdout;
  unsigned int messageLevel = 0;

  // Declare HighsOptions for an LP model, any solver and simplex solver, setting the default value
  //
  // For an LP model
  //
  // Try to solve the dual of the LP
  bool transpose_solver_lp = false;
  // Perform LP scaling
  bool scale_solver_lp = true;
  // Permute the columns of the LP randomly to aid load distribution in block parallelism
  bool permute_solver_lp = false;
  // Perform LP bound tightening
  bool tighten_solver_lp = false;
  //
  // For any solver
  //
 // primal feasibility (dual optimality) tolerance
  double primal_feasibility_tolerance = PRIMAL_FEASIBILITY_TOLERANCE_DEFAULT;
  // dual feasibility (primal optimality) tolerance
  double dual_feasibility_tolerance = DUAL_FEASIBILITY_TOLERANCE_DEFAULT;


  // Upper bound on dual objective value
  double dual_objective_value_upper_bound = DUAL_OBJECTIVE_VALUE_UPPER_BOUND_DEFAULT;
  //
  // For the simplex solver
  //
  bool simplex_perturb_costs = true;
  // Maximum number of simplex iterations
  int simplex_iteration_limit = SIMPLEX_ITERATION_LIMIT_DEFAULT;
  int simplex_update_limit = SIMPLEX_UPDATE_LIMIT_DEFAULT;

  bool clean_up = false;
};

class HighsLp {
 public:
  // Model data
  int numCol_ = 0;
  int numRow_ = 0;
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

};

// HiGHS status
enum class HighsStatus {
  OK,
  Init,
  LpError,
  OptionsError,
  PresolveError,
  SolutionError,
  PostsolveError,
  NotImplemented,
  ReachedDualObjectiveUpperBound,
  Unbounded,
  Infeasible,
  Feasible,
  Optimal,
  Timeout
};

enum class HighsInputStatus {
  OK,
  FileNotFound,
  ErrorMatrixDimensions,
  ErrorMatrixIndices,
  ErrorMatrixStart,
  ErrorMatrixValue,
  ErrorColBounds,
  ErrorRowBounds,
  ErrorObjective
};

// Cost, column and row scaling factors
struct HighsScale {
  double cost_;
  std::vector<double> col_;
  std::vector<double> row_;
};

struct HighsBasis {
  std::vector<int> basicIndex_;
  std::vector<int> nonbasicFlag_;
  std::vector<int> nonbasicMove_;
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
  // mlFg_haveNonbasicDuals.
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
  // them is indicated by mlFg_haveBasicPrimals.
  //
  std::vector<double> baseLower_;
  std::vector<double> baseUpper_;
  std::vector<double> baseValue_;
  //
  // Vectors of random reals for column cost perturbation, a random
  // permutation of all indices for CHUZR and a random permutation of
  // column indices for shuffling the columns
  std::vector<double> numTotRandomValue_;
  std::vector<int> numTotPermutation_;
  std::vector<int> numColPermutation_;

  // Values of iClock for simplex timing clocks
  std::vector<int> clock_;
  //
  // Options from HighsOptions for the simplex solver
  double highs_run_time_limit;
  SimplexStrategy simplex_strategy;
  SimplexCrashStrategy crash_strategy;
  SimplexDualEdgeWeightStrategy dual_edge_weight_strategy;
  SimplexPriceStrategy price_strategy;

  double primal_feasibility_tolerance;
  double dual_feasibility_tolerance;
  bool perturb_costs;
  int update_limit;
  int iteration_limit;
  double dual_objective_value_upper_bound;
  
  // Options for the LP to be solved
  bool transpose_solver_lp;
  bool scale_solver_lp;
  bool permute_solver_lp;
  bool tighten_solver_lp;
  // Internal options - can't be changed externally

  // Options for reporting timing
  bool reportSimplexInnerClock;
  bool reportSimplexOuterClock;
  bool reportSimplexPhasesClock;
#ifdef HiGHSDEV
  // Option for analysing simplex iterations, INVERT time and rebuild time
  bool analyseLp;
  bool analyseSimplexIterations;
  bool analyseLpSolution;
  bool analyseInvertTime;
  bool analyseRebuildTime;
#endif
  // Solved LP status
  bool transposed_solver_lp = false;
  bool scaled_solver_lp = false;
  bool permuted_solver_lp = false;
  bool tightened_solver_lp = false;

  // Simplex status

  // Simplex runtime information
  SimplexSolutionStatus solution_status = SimplexSolutionStatus::UNSET;
  int costs_perturbed = 0;
  // Cumulative iteration count - updated in simplex solvers
  int iteration_count = 0;
  // Records of cumulative iteration counts - updated at the end of a phase
  int dual_phase1_iteration_count = 0;
  int dual_phase2_iteration_count = 0;
  int primal_phase1_iteration_count = 0;
  int primal_phase2_iteration_count = 0;

  // Number of UPDATE operations performed - should be zeroed when INVERT is performed
  int update_count;
  // Value of dual objective - only set when computed from scratch in rebuild()
  double dualObjectiveValue;


  // Value of dual objective that is updated in dual simplex solver
  double updatedDualObjectiveValue;
  // Number of logical variables in the basis 
  int num_basic_logicals;

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
  std::vector<double> colValue_;
  std::vector<double> colDual_;
  std::vector<double> rowValue_;
  std::vector<double> rowDual_;
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

// Return a string representation of HighsStatus.
std::string HighsStatusToString(HighsStatus status);

// Return a string representation of ParseStatus.
std::string HighsInputStatusToString(HighsInputStatus status);

// If debug this method terminates the program when the status is not OK. If
// standard build it only prints a message.
void checkStatus(HighsStatus status);

HighsInputStatus checkLp(const HighsLp& lp);

#endif
