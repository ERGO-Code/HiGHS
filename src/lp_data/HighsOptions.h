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

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "io/HighsIO.h"
#include "presolve/Presolve.h"
#include "simplex/SimplexConst.h"

enum class OptionStatus
{
  OK = 0,
  NO_FILE,
  UNKNOWN_OPTION,
  ILLEGAL_VALUE
};

const string on_string = "on";
const string off_string = "off";

const string fixed_string = "fixed";
const string free_string = "free";

// Strings for command line options
const string file_string = "file";
const string presolve_string = "presolve";
const string crash_string = "crash";
const string parallel_string = "parallel";
const string simplex_string = "simplex";
const string ipm_string = "ipm";
const string highs_run_time_limit_string = "highs_run_time_limit";
const string simplex_iteration_limit_string = "simplex_iteration_limit";
const string options_file_string = "options_file";
const string parser_type_string = "parser_type";
const string mip_string = "mip";
const string find_feasibility_string = "find_feasibility";
const string find_feasibility_strategy_string = "feasibility_strategy";
const string find_feasibility_dualize_string = "feasibility_dualize"; 

// Strings for file options
const string infinite_cost_string = "infinite_cost";
const string infinite_bound_string = "infinite_bound";
const string small_matrix_value_string = "small_matrix_value";
const string large_matrix_value_string = "large_matrix_value";
const string primal_feasibility_tolerance_string = "primal_feasibility";
const string dual_feasibility_tolerance_string = "dual_feasibility";
const string dual_objective_value_upper_bound_string = "dual_objective_value_upper_bound";

const string simplex_strategy_string = "simplex_strategy";
const string simplex_crash_strategy_string = "simplex_crash_strategy";
const string simplex_dual_edge_weight_strategy_string = "simplex_dual_edge_weight_strategy";
const string simplex_price_strategy_string = "simplex_price_strategy";

const string message_level_string = "message_level";

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
  PresolveOption presolve_option = PresolveOption::DEFAULT;
  CrashOption crash_option = CrashOption::DEFAULT;
  ParallelOption parallel_option = ParallelOption::DEFAULT;
  SimplexOption simplex_option = SimplexOption::DEFAULT;
  bool ipx = false;
  double highs_run_time_limit = HIGHS_RUN_TIME_LIMIT_DEFAULT;
  int simplex_iteration_limit = SIMPLEX_ITERATION_LIMIT_DEFAULT;
  HighsMpsParserType parser_type = HighsMpsParserType::fixed;

  // Options not passed through the command line
  double infinite_cost = INFINITE_COST_DEFAULT;
  double infinite_bound = INFINITE_BOUND_DEFAULT;
  double small_matrix_value = SMALL_MATRIX_VALUE_DEFAULT;
  double large_matrix_value = LARGE_MATRIX_VALUE_DEFAULT;
  double primal_feasibility_tolerance = PRIMAL_FEASIBILITY_TOLERANCE_DEFAULT;
  double dual_feasibility_tolerance = DUAL_FEASIBILITY_TOLERANCE_DEFAULT;
  double dual_objective_value_upper_bound = DUAL_OBJECTIVE_VALUE_UPPER_BOUND_DEFAULT;
  SimplexStrategy simplex_strategy = SimplexStrategy::DEFAULT;
  SimplexCrashStrategy simplex_crash_strategy = SimplexCrashStrategy::DEFAULT;
  SimplexDualEdgeWeightStrategy simplex_dual_edge_weight_strategy = SimplexDualEdgeWeightStrategy::DEFAULT;
  SimplexPriceStrategy simplex_price_strategy = SimplexPriceStrategy::DEFAULT;

  int allow_superbasic = false;

  bool pami = 0;
  bool sip = 0;
  bool scip = 0;

  // Options for HighsPrintMessage and HighsLogMessage
  FILE *logfile = stdout;
  FILE *output = stdout;
  unsigned int messageLevel = ML_MINIMAL;

  void (*printmsgcb) (unsigned int level, const char* msg, void* msgcb_data) = NULL;
  void (*logmsgcb) (HighsMessageType type, const char* msg, void* msgcb_data) = NULL;
  void* msgcb_data = NULL;

  // Declare HighsOptions for an LP model, any solver and simplex solver, setting the default value
  //
  // For an LP model
  //
  // Try to solve the dual of the LP
  bool transpose_simplex_lp = false;
  // Perform LP scaling
  bool scale_simplex_lp = true;
  // Permute the columns of the LP randomly to aid load distribution in block parallelism
  bool permute_simplex_lp = false;
  // Perform LP bound tightening
  bool tighten_simplex_lp = false;
  //
  // For any solver
  //

  //
  // For the simplex solver
  //
  bool simplex_perturb_costs = true;
  // Maximum number of simplex updates
  int simplex_update_limit = SIMPLEX_UPDATE_LIMIT_DEFAULT;

  bool clean_up = false;
  bool find_feasibility = false;
  FeasibilityStrategy feasibility_strategy = FeasibilityStrategy::kApproxComponentWise;
  bool feasibility_strategy_dualize = false;

  bool mip = false;
};

OptionStatus setOptionValue(HighsOptions& options, const std::string& option, const std::string& value);

// Called before solve. This would check whether tolerances are set to correct values and
// all options are consistent.
OptionStatus checkOptionsValue(HighsOptions& options);

OptionStatus setPresolveValue(HighsOptions& options, const std::string& value);
OptionStatus setCrashValue(HighsOptions& options, const std::string& value);
OptionStatus setParallelValue(HighsOptions& options, const std::string& value);
OptionStatus setSimplexValue(HighsOptions& options, const std::string& value);
OptionStatus setIpmValue(HighsOptions& options, const std::string& value);
OptionStatus setHighsRunTimeLimitValue(HighsOptions& options, const double& value);
OptionStatus setSimplexIterationLimitValue(HighsOptions& options, const int& value);
OptionStatus setParserTypeValue(HighsOptions& options, const std::string& value);
OptionStatus setMipValue(HighsOptions& options, const std::string& value);
OptionStatus setFindFeasibilityValue(HighsOptions& options, const std::string& value);
OptionStatus setFindFeasibilityStrategyValue(HighsOptions& options, const std::string& value);
OptionStatus setFindFeasibilityDualizeValue(HighsOptions& options, const std::string& value);
 

OptionStatus setInfiniteCostValue(HighsOptions& options, const double& value);
OptionStatus setInfiniteBoundValue(HighsOptions& options, const double& value);
OptionStatus setSmallMatrixValueValue(HighsOptions& options, const double& value);
OptionStatus setLargeMatrixValueValue(HighsOptions& options, const double& value);
OptionStatus setPrimalFeasibilityToleranceValue(HighsOptions& options, const double& value);
OptionStatus setDualFeasibilityToleranceValue(HighsOptions& options, const double& value);
OptionStatus setDualObjectiveValueUpperBoundValue(HighsOptions& options, const double& value);
OptionStatus setSimplexStrategyValue(HighsOptions& options, const int& value);
OptionStatus setSimplexCrashStrategyValue(HighsOptions& options, const int& value);
OptionStatus setSimplexDualEdgeWeightStrategyValue(HighsOptions& options, const int& value);
OptionStatus setSimplexPriceStrategyValue(HighsOptions& options, const int& value);

OptionStatus setMessageLevelValue(HighsOptions& options, const int& value);

SimplexStrategy intToSimplexStrategy(const int& value);
SimplexCrashStrategy intToSimplexCrashStrategy(const int& value);
SimplexDualEdgeWeightStrategy intToSimplexDualEdgeWeightStrategy(const int& value);
SimplexPriceStrategy intToSimplexPriceStrategy(const int& value);

#endif

