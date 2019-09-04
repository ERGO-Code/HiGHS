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

#include "io/HighsIO.h"
#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"
#include "presolve/Presolve.h"
#include "simplex/SimplexConst.h"

enum class OptionStatus { OK = 0, NO_FILE, UNKNOWN_OPTION, ILLEGAL_VALUE };

class OptionRecord {
 public:
  HighsOptionType type;
  std::string name;
  std::string description;
  bool advanced;
  
  OptionRecord(HighsOptionType Xtype, std::string Xname, std::string Xdescription, bool Xadvanced) {
    this->type = Xtype;
    this->name = Xname;
    this->description = Xdescription;
    this->advanced = Xadvanced;
  }
  
  ~OptionRecord() {}
};

class OptionRecordBool : public OptionRecord {
 public:
  bool* value;
  bool default_value;
 OptionRecordBool(
		 std::string Xname,
		 std::string Xdescription,
		 bool Xadvanced,
		 bool* Xvalue_pointer,
		 bool Xdefault_value) : OptionRecord(
						    HighsOptionType::BOOL,
						    Xname,
						    Xdescription,
						    Xadvanced)  {
    advanced = Xadvanced;
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }
  
  void assignvalue(bool Xvalue) {
    *value = Xvalue;
  }
  
  ~OptionRecordBool() {}
};

class OptionRecordInt : public OptionRecord {
 public:
  int* value;
  int lower_bound;
  int upper_bound;
  int default_value;
 OptionRecordInt(
		std::string Xname,
		std::string Xdescription,
		bool Xadvanced,
		int* Xvalue_pointer,
		int Xlower_bound,
		int Xupper_bound,
		int Xdefault_value) : OptionRecord(
						  HighsOptionType::INT,
						  Xname,
						  Xdescription,
						  Xadvanced) {
    value = Xvalue_pointer;
    lower_bound = Xlower_bound;
    upper_bound = Xupper_bound;
    default_value = Xdefault_value;
    *value = default_value;
  }

void assignvalue(int Xvalue) {
  *value = Xvalue;
}

~OptionRecordInt() {}
};

class OptionRecordDouble : public OptionRecord {
 public:
  double* value;
  double lower_bound;
  double upper_bound;
  double default_value;
  OptionRecordDouble(std::string Xname,
		    std::string Xdescription,
		    bool Xadvanced,
		    double* Xvalue_pointer,
		    double Xlower_bound,
		    double Xupper_bound,
		    double Xdefault_value) : OptionRecord(
							 HighsOptionType::DOUBLE,
							 Xname,
							 Xdescription,
							 Xadvanced)  {
    value = Xvalue_pointer;
    lower_bound = Xlower_bound;
    upper_bound = Xupper_bound;
    default_value = Xdefault_value;
    *value = default_value;
  }

void assignvalue(double Xvalue) {
  *value = Xvalue;
}

~OptionRecordDouble() {}
};

class OptionRecordString : public OptionRecord {
 public:
  std::string* value;
  std::string default_value;
  OptionRecordString(
		    std::string Xname,
		    std::string Xdescription,
		    bool Xadvanced,
		    std::string* Xvalue_pointer,
		    std::string Xdefault_value) : OptionRecord(
							      HighsOptionType::STRING,
							      Xname,
							      Xdescription,
							      Xadvanced)  {
    value = Xvalue_pointer;
    default_value = Xdefault_value;
    *value = default_value;
  }

void assignvalue(std::string Xvalue) {
  *value = Xvalue;
}

~OptionRecordString() {}
};

inline const char* bool2string(bool b);

bool commandLineOffChooseOnOk(const string& value);
bool commandLineSolverOk(const string& value);

bool boolFromString(const std::string value, bool& bool_value);

OptionStatus getOptionIndex(const std::string& name, const std::vector<OptionRecord*>& option_records, int& index);

OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const bool value);
OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const int value);
OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const double value);
OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const std::string value);

OptionStatus setOptionValue(OptionRecordBool& option, const bool value);
OptionStatus setOptionValue(OptionRecordInt& option, const int value);
OptionStatus setOptionValue(OptionRecordDouble& option, const double value);
OptionStatus setOptionValue(OptionRecordString& option, std::string const value);

OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, bool& value);
OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, int& value);
OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, double& value);
OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, std::string& value);

void reportOptions(FILE* file, const std::vector<OptionRecord*>& option_records, const bool force_report=false);
void reportOption(FILE* file, const OptionRecordBool& option, const bool force_report=false);
void reportOption(FILE* file, const OptionRecordInt& option, const bool force_report=false);
void reportOption(FILE* file, const OptionRecordDouble& option, const bool force_report=false);
void reportOption(FILE* file, const OptionRecordString& option, const bool force_report=false);

//======================================


const string simplex_string = "simplex";
const string ipm_string = "ipm";
const int KEEP_N_ROWS_DELETE_ROWS = -1;
const int KEEP_N_ROWS_DELETE_ENTRIES = 0;
const int KEEP_N_ROWS_KEEP_ROWS = 1;

// Strings for command line options
const string model_file_string = "model_file";
const string presolve_string = "presolve";
const string solver_string = "solver";
const string parallel_string = "parallel";
const string time_limit_string = "time_limit";
const string options_file_string = "options_file";

/** SCIP/HiGHS Objective sense */
enum objSense { OBJSENSE_MINIMIZE = 1, OBJSENSE_MAXIMIZE = -1 };

// For now, but later change so HiGHS properties are string based so that new
// options (for debug and testing too) can be added easily. The options below
// are just what has been used to parse options from argv.
// todo: when creating the new options don't forget underscores for class
// variables but no underscores for struct
class HighsOptions {
 public:
  HighsOptions() {
    OptionRecordBool* record_bool;
    OptionRecordInt* record_int;
    OptionRecordDouble* record_double;
    OptionRecordString* record_string;
    // Options read from the command line
    record_string = new OptionRecordString(model_file_string,
					   "Model file",
					   false, &model_file,
					   FILENAME_DEFAULT);
    records.push_back(record_string);
    record_string = new OptionRecordString(presolve_string,
					   "Presolve option: \"off\", \"choose\" or \"on\"",
					   false, &presolve,
					   choose_string);
    records.push_back(record_string);
    record_string = new OptionRecordString(solver_string,
					   "Solver option: \"simplex\", \"choose\" or \"ipm\"",
					   false, &solver,
					   choose_string);
    records.push_back(record_string);
    record_string = new OptionRecordString(parallel_string,
					   "Parallel option: \"off\", \"choose\" or \"on\"",
					   false, &parallel,
					   choose_string);
    records.push_back(record_string);
    record_double = new OptionRecordDouble(time_limit_string,
					   "Time limit",
					   false, &time_limit,
					   0, HIGHS_CONST_INF, HIGHS_CONST_INF);
    records.push_back(record_double);
    record_string = new OptionRecordString(options_file_string,
					   "Options file",
					   false, &options_file,
					   FILENAME_DEFAULT);
    records.push_back(record_string);
    // Options read from the file
    record_int = new OptionRecordInt("simplex_iteration_limit",
				     "Iteration limit for simplex solver",
				     false, &simplex_iteration_limit,
				     0, HIGHS_CONST_I_INF, HIGHS_CONST_I_INF);
    records.push_back(record_int);

    // Advanced options
    record_bool = new OptionRecordBool("run_as_hsol",
				       "Run HiGHS simplex solver as if it were hsol",
				       true, &run_as_hsol,
				       false);
    records.push_back(record_bool);
    record_bool = new OptionRecordBool("mps_parser_type_free",
				       "Use the free format MPS file reader",
				       true, &mps_parser_type_free,
				       true);
    records.push_back(record_bool);
    record_int = new OptionRecordInt("keep_n_rows",
				     "For multiple N-rows in MPS files: delete rows / delete entries / keep rows -1/0/1",
				     true, &keep_n_rows,
				     KEEP_N_ROWS_DELETE_ROWS, KEEP_N_ROWS_KEEP_ROWS, KEEP_N_ROWS_DELETE_ROWS);
    records.push_back(record_int);

  }
  std::vector<OptionRecord*> records;

  // Options read from the command line
  std::string model_file;
  std::string presolve;
  std::string solver;
  std::string parallel;
  double time_limit;
  std::string options_file;
  
  // Options read from the file
  int simplex_iteration_limit;

  // Advanced options
  bool run_as_hsol;
  bool mps_parser_type_free;
  int keep_n_rows;

  double infinite_cost = INFINITE_COST_DEFAULT;
  double infinite_bound = INFINITE_BOUND_DEFAULT;
  double small_matrix_value = SMALL_MATRIX_VALUE_DEFAULT;
  double large_matrix_value = LARGE_MATRIX_VALUE_DEFAULT;
  int allowed_simplex_scale_factor = ALLOWED_SIMPLEX_SCALE_FACTOR_DEFAULT;
  double primal_feasibility_tolerance = PRIMAL_FEASIBILITY_TOLERANCE_DEFAULT;
  double dual_feasibility_tolerance = DUAL_FEASIBILITY_TOLERANCE_DEFAULT;
  double dual_objective_value_upper_bound =
      DUAL_OBJECTIVE_VALUE_UPPER_BOUND_DEFAULT;
  SimplexStrategy simplex_strategy = SimplexStrategy::DEFAULT;
  SimplexDualiseStrategy simplex_dualise_strategy =
      SimplexDualiseStrategy::DEFAULT;
  SimplexPermuteStrategy simplex_permute_strategy =
      SimplexPermuteStrategy::DEFAULT;
  SimplexScaleStrategy simplex_scale_strategy = SimplexScaleStrategy::DEFAULT;
  SimplexCrashStrategy simplex_crash_strategy = SimplexCrashStrategy::DEFAULT;
  SimplexDualEdgeWeightStrategy simplex_dual_edge_weight_strategy =
      SimplexDualEdgeWeightStrategy::DEFAULT;
  SimplexPrimalEdgeWeightStrategy simplex_primal_edge_weight_strategy =
      SimplexPrimalEdgeWeightStrategy::DEFAULT;
  SimplexPriceStrategy simplex_price_strategy = SimplexPriceStrategy::DEFAULT;

  bool simplex_initial_condition_check = true;
  double simplex_initial_condition_tolerance =
      SIMPLEX_INITIAL_CONDITION_TOLERANCE_DEFAULT;
  double dual_steepest_edge_weight_log_error_threshhold =
    DUAL_STEEPEST_EDGE_WEIGHT_LOG_ERROR_THRESHHOLD_DEFAULT;
  int allow_superbasic = false;

  bool pami = 0;
  bool sip = 0;
  bool scip = 0;

  // Options for HighsPrintMessage and HighsLogMessage
  FILE* logfile = stdout;
  FILE* output = stdout;
  unsigned int messageLevel = ML_MINIMAL;

  void (*printmsgcb)(unsigned int level, const char* msg,
                     void* msgcb_data) = NULL;
  void (*logmsgcb)(HighsMessageType type, const char* msg,
                   void* msgcb_data) = NULL;
  void* msgcb_data = NULL;

  // Declare HighsOptions for an LP model, any solver and simplex solver,
  // setting the default value
  //
  // For the simplex solver
  //
  bool simplex_perturb_costs = true;
  // Maximum number of simplex updates
  int simplex_update_limit = SIMPLEX_UPDATE_LIMIT_DEFAULT;

  bool find_feasibility = false;
  FeasibilityStrategy feasibility_strategy =
      FeasibilityStrategy::kApproxComponentWise;
  bool feasibility_strategy_dualize = false;

  bool mip = false;

};


// Called before solve. This would check whether tolerances are set to correct
// values and all options are consistent.
OptionStatus checkOptionsValue(HighsOptions& options);
void setHsolOptions(HighsOptions& options);

SimplexStrategy intToSimplexStrategy(const int& value);
SimplexDualiseStrategy intToSimplexDualiseStrategy(const int& value);
SimplexPermuteStrategy intToSimplexPermuteStrategy(const int& value);
SimplexScaleStrategy intToSimplexScaleStrategy(const int& value);
SimplexCrashStrategy intToSimplexCrashStrategy(const int& value);
SimplexDualEdgeWeightStrategy intToSimplexDualEdgeWeightStrategy(
    const int& value);
SimplexPrimalEdgeWeightStrategy intToSimplexPrimalEdgeWeightStrategy(
    const int& value);
SimplexPriceStrategy intToSimplexPriceStrategy(const int& value);

#endif
