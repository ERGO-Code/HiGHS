/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsOptions.h
 * @brief 
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef LP_DATA_HIGHS_OPTIONS_H_
#define LP_DATA_HIGHS_OPTIONS_H_

#include <cstring>

#include "HConst.h"
#include "HighsLp.h"
#include "Presolve.h"
#include "HighsModelObject.h"
#include "SimplexConst.h"
#include "cxxopts.hpp"

// The free parser also reads fixed format MPS files but the fixed
// parser does not read free mps files.
enum class HighsMpsParserType
{
  free,
  fixed
};

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
struct HighsOptions
{
  std::string filename = "";
  std::string options_file = "";

  // Options passed through the command line

  ParallelOption parallel_option = ParallelOption::DEFAULT;
  PresolveOption presolve_option = PresolveOption::DEFAULT;
  CrashOption crash_option = CrashOption::DEFAULT;
  SimplexOption simplex_option = SimplexOption::DEFAULT;
  bool ipx = false;
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
  FILE *logfile = stdout;
  FILE *output = stdout;
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

// Used only for options allowed for the user. For other options see setOptionValue.
bool setUserOptionValue(HighsOptions& options, const std::string& option, const std::string& value);

// Used for extended options read from file or set internally.
bool setOptionValue(HighsOptions& options, const std::string& option, const std::string& value);

// Called before sovle. This would check whether tolerances are set to correct values and
// all options are consistent.
bool checkOptionsValue(HighsOptions& options);

#endif