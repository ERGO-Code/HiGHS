/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file lp_data/HighsOptions.cpp
 * @brief
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "lp_data/HighsOptions.h"
#include "io/HighsIO.h"

inline const char* bool2string(bool b) { return b ? "true" : "false"; }

std::string optionEntryType2string(const HighsOptionType type) {
  if (type == HighsOptionType::BOOL) {
    return "bool";
  } else if (type == HighsOptionType::INT) {
    return "int";
  } else if (type == HighsOptionType::DOUBLE) {
    return "double";
  } else
    return "string";
}

bool commandLineOffChooseOnOk(const string& value) {
  if (value == off_string || value == choose_string || value == on_string) return true;
  HighsLogMessage(HighsMessageType::ERROR, "Command line option value \"%s\" is not one of \"%s\", \"%s\" or \"%s\"\n",
		  value.c_str(), off_string.c_str(), choose_string.c_str(), on_string.c_str());
  return false;
}

bool commandLineSolverOk(const string& value) {
  if (value == simplex_string || value == choose_string || value == ipm_string) return true;
  HighsLogMessage(HighsMessageType::ERROR, "Command line option value \"%s\" is not one of \"%s\", \"%s\" or \"%s\"\n",
		  value.c_str(), simplex_string.c_str(), choose_string.c_str(), ipm_string.c_str());
  return false;
}

bool boolFromString(const std::string value, bool& bool_value) {
  if (
      value == "t" ||
      value == "true" ||
      value == "T" ||
      value == "True" ||
      value == "TRUE") {
    bool_value = true;
  } else if (
      value == "f" ||
      value == "false" ||
      value == "F" ||
      value == "False" ||
      value == "FALSE") {
    bool_value = false;
  } else {
    return false;
  }
  return true;
}

OptionStatus getOptionIndex(const std::string& name, const std::vector<OptionRecord*>& option_records, int& index) {
  int num_options = option_records.size();
  for (index=0; index<num_options; index++) if (option_records[index]->name == name) return OptionStatus::OK;
  HighsLogMessage(HighsMessageType::ERROR, "getOptionIndex: Option \"%s\" is unknown", name.c_str());
  return OptionStatus::UNKNOWN_OPTION;
}

OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const bool value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::BOOL) {
    HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Option \"%s\" cannot be assigned a bool", name.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return setOptionValue(((OptionRecordBool*)option_records[index])[0], value);
}

OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const int value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::INT) {
    HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Option \"%s\" cannot be assigned an int", name.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return setOptionValue(((OptionRecordInt*)option_records[index])[0], value);
}

OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const double value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::DOUBLE) {
    HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Option \"%s\" cannot be assigned a double", name.c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  return setOptionValue(((OptionRecordDouble*)option_records[index])[0], value);
}

OptionStatus setOptionValue(const std::string& name, std::vector<OptionRecord*>& option_records, const std::string value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type == HighsOptionType::BOOL) {
    bool bool_value;
    bool return_status = boolFromString(value, bool_value);
    if (!return_status) {
      HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Value \"%s\" cannot be interpreted as a bool", value.c_str());
      return OptionStatus::ILLEGAL_VALUE;
    }
    return setOptionValue(((OptionRecordBool*)option_records[index])[0], bool_value);
  } else if (type == HighsOptionType::INT) {
    return setOptionValue(((OptionRecordInt*)option_records[index])[0], atoi(value.c_str()));
  } else if (type == HighsOptionType::DOUBLE) {
    return setOptionValue(((OptionRecordDouble*)option_records[index])[0], atof(value.c_str()));
  } else {
    return setOptionValue(((OptionRecordString*)option_records[index])[0], value);
  }
}

OptionStatus setOptionValue(OptionRecordBool& option, const bool value) {
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus setOptionValue(OptionRecordInt& option, const int value) {
  if (value < option.lower_bound) {
    HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Trying to set option \"%s\" to value %d below lower bound of %d",
	   option.name.c_str(), value, option.lower_bound);
    return OptionStatus::ILLEGAL_VALUE;
  } else if (value > option.upper_bound) {
    HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Trying to set option \"%s\" to value %d above upper bound of %d",
	   option.name.c_str(), value, option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus setOptionValue(OptionRecordDouble& option, const double value) {
  if (value < option.lower_bound) {
    HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Trying to set option \"%s\" to value %g below lower bound of %g",
	   option.name.c_str(), value, option.lower_bound);
    return OptionStatus::ILLEGAL_VALUE;
  } else if (value > option.upper_bound) {
    HighsLogMessage(HighsMessageType::ERROR, "setOptionValue: Trying to set option \"%s\" to value %g above upper bound of %g",
	   option.name.c_str(), value, option.upper_bound);
    return OptionStatus::ILLEGAL_VALUE;
  }
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus setOptionValue(OptionRecordString& option, const std::string value) {
  // Setting a string option: check that value is OK
  if (option.name == presolve_string) {
    if (!commandLineOffChooseOnOk(value)) return OptionStatus::ILLEGAL_VALUE;
  } else if (option.name == solver_string) {
    if (!commandLineSolverOk(value)) return OptionStatus::ILLEGAL_VALUE;
  } else if (option.name == parallel_string) {
    if (!commandLineOffChooseOnOk(value)) return OptionStatus::ILLEGAL_VALUE;
  }
  option.assignvalue(value);
  return OptionStatus::OK;
}

OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, bool& value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::BOOL) {
    HighsLogMessage(HighsMessageType::ERROR, "getOptionValue: Option \"%s\" requires value of type %s, not bool", name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordBool option = ((OptionRecordBool*)option_records[index])[0];
  value = *(option.value);
  return OptionStatus::OK;
}

OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, int& value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::INT) {
    HighsLogMessage(HighsMessageType::ERROR, "getOptionValue: Option \"%s\" requires value of type %s, not int", name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordInt option = ((OptionRecordInt*)option_records[index])[0];
  value = *(option.value);
  return OptionStatus::OK;
}

OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, double& value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::DOUBLE) {
    HighsLogMessage(HighsMessageType::ERROR, "getOptionValue: Option \"%s\" requires value of type %s, not double", name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordDouble option = ((OptionRecordDouble*)option_records[index])[0];
  value = *(option.value);
  return OptionStatus::OK;
}

OptionStatus getOptionValue(const std::string& name, const std::vector<OptionRecord*>& option_records, std::string& value) {
  int index;
  OptionStatus status = getOptionIndex(name, option_records, index);
  if (status != OptionStatus::OK) return status;
  HighsOptionType type = option_records[index]->type;
  if (type != HighsOptionType::STRING) {
    HighsLogMessage(HighsMessageType::ERROR, "getOptionValue: Option \"%s\" requires value of type %s, not string", name.c_str(), optionEntryType2string(type).c_str());
    return OptionStatus::ILLEGAL_VALUE;
  }
  OptionRecordString option = ((OptionRecordString*)option_records[index])[0];
  value = *(option.value);
  return OptionStatus::OK;
}

void reportOptions(FILE* file, const std::vector<OptionRecord*>& option_records, const bool force_report) {
  int num_options = option_records.size();
  for (int i=0; i<num_options; i++) {
    HighsOptionType type = option_records[i]->type;
    //    fprintf(file, "\n# Option %1d\n", i);
    if (type == HighsOptionType::BOOL) {
      reportOption(file, ((OptionRecordBool*)option_records[i])[0], force_report);
    } else if (type == HighsOptionType::INT) {
      reportOption(file, ((OptionRecordInt*)option_records[i])[0], force_report);
    } else if (type == HighsOptionType::DOUBLE) {
      reportOption(file, ((OptionRecordDouble*)option_records[i])[0], force_report);
    } else {
      reportOption(file, ((OptionRecordString*)option_records[i])[0], force_report);
    } 
  }
}

void reportOption(FILE* file, const OptionRecordBool& option, const bool force_report) {
  if (force_report || option.default_value != *(option.value)) {
    fprintf(file, "\n# %s\n", option.description.c_str());
    fprintf(file, "# [type: bool, advanced: %s, range: {false, true}, default: %s]\n",
	   bool2string(option.advanced),
	   bool2string(option.default_value));
    fprintf(file, "%s = %s\n", option.name.c_str(), bool2string(*(option.value)));
  }
}

void reportOption(FILE* file, const OptionRecordInt& option, const bool force_report) {
  if (force_report || option.default_value != *(option.value)) {
    fprintf(file, "\n# %s\n", option.description.c_str());
    fprintf(file, "# [type: int, advanced: %s, range: {%d, %d}, default: %d]\n",
	   bool2string(option.advanced),
	   option.lower_bound,
	   option.upper_bound,
	   option.default_value);
    fprintf(file, "%s = %d\n", option.name.c_str(), *(option.value));
  }
}

void reportOption(FILE* file, const OptionRecordDouble& option, const bool force_report) {
  if (force_report || option.default_value != *(option.value)) {
    fprintf(file, "\n# %s\n", option.description.c_str());
    fprintf(file, "# [type: double, advanced: %s, range: [%g, %g], default: %g]\n",
	   bool2string(option.advanced),
	   option.lower_bound,
	   option.upper_bound,
	   option.default_value);
    fprintf(file, "%s = %g\n", option.name.c_str(), *(option.value));
  }
}

void reportOption(FILE* file, const OptionRecordString& option, const bool force_report) {
  // Don't report for the options file if writing to an options file
  if (
      //file != stdout &&
      option.name == options_file_string) return;
  if (force_report || option.default_value != *(option.value)) {
    fprintf(file, "\n# %s\n", option.description.c_str());
    fprintf(file, "# [type: string, advanced: %s, default: \"%s\"]\n",
	   bool2string(option.advanced),
	   option.default_value.c_str());
    fprintf(file, "%s = %s\n", option.name.c_str(), (*(option.value)).c_str());
  }
}

//==========================================

OptionStatus checkOptionsValue(HighsOptions& options) {
  return OptionStatus::OK;
}

// Set values of options so that HiGHS runs as Hsol
void setHsolOptions(HighsOptions& options) {
  // Set command line options to their hsol values
  options.presolve = OPTION_OFF;
  options.solver = SOLVER_OPTION_SIMPLEX;
  options.parallel = OPTION_OFF;
  options.time_limit = HIGHS_CONST_INF;
  
  options.simplex_iteration_limit = HIGHS_CONST_I_INF;
  options.mps_parser_type_free = false;
  options.keep_n_rows = KEEP_N_ROWS_KEEP_ROWS;
  options.infinite_cost = HIGHS_CONST_INF;
  options.infinite_bound = HIGHS_CONST_INF;
  options.small_matrix_value = 0;
  options.large_matrix_value = HIGHS_CONST_INF;
  options.allowed_simplex_scale_factor = HIGHS_CONST_I_INF;
  options.primal_feasibility_tolerance = 1e-7;
  options.dual_feasibility_tolerance = 1e-7;
  options.dual_objective_value_upper_bound = HIGHS_CONST_INF;
  options.simplex_strategy = SimplexStrategy::DUAL_PLAIN;
  options.simplex_dualise_strategy = SimplexDualiseStrategy::OFF;
  options.simplex_permute_strategy = SimplexPermuteStrategy::OFF;
  options.simplex_scale_strategy = SimplexScaleStrategy::HSOL;
  options.simplex_crash_strategy = SimplexCrashStrategy::OFF;
  options.simplex_dual_edge_weight_strategy = SimplexDualEdgeWeightStrategy::STEEPEST_EDGE;
  options.simplex_primal_edge_weight_strategy = SimplexPrimalEdgeWeightStrategy::DANTZIG;
  options.simplex_price_strategy = SimplexPriceStrategy::ROW;
}

OptionStatus setMessageLevelValue(HighsOptions& options, const int& value) {
  HighsSetMessagelevel(value);
  options.messageLevel = value;
  return OptionStatus::OK;
}

SimplexStrategy intToSimplexStrategy(const int& value) {
  if (value == (int)SimplexStrategy::CHOOSE) return SimplexStrategy::CHOOSE;
  if (value == (int)SimplexStrategy::DUAL_PLAIN)
    return SimplexStrategy::DUAL_PLAIN;
  if (value == (int)SimplexStrategy::DUAL_MULTI)
    return SimplexStrategy::DUAL_MULTI;
  if (value == (int)SimplexStrategy::PRIMAL) return SimplexStrategy::PRIMAL;
  return SimplexStrategy::DEFAULT;
}
SimplexDualiseStrategy intToSimplexDualiseStrategy(const int& value) {
  if (value == (int)SimplexDualiseStrategy::OFF)
    return SimplexDualiseStrategy::OFF;
  if (value == (int)SimplexDualiseStrategy::CHOOSE)
    return SimplexDualiseStrategy::CHOOSE;
  if (value == (int)SimplexDualiseStrategy::ON)
    return SimplexDualiseStrategy::ON;
  return SimplexDualiseStrategy::DEFAULT;
}
SimplexPermuteStrategy intToSimplexPermuteStrategy(const int& value) {
  if (value == (int)SimplexPermuteStrategy::OFF)
    return SimplexPermuteStrategy::OFF;
  if (value == (int)SimplexPermuteStrategy::CHOOSE)
    return SimplexPermuteStrategy::CHOOSE;
  if (value == (int)SimplexPermuteStrategy::ON)
    return SimplexPermuteStrategy::ON;
  return SimplexPermuteStrategy::DEFAULT;
}
SimplexScaleStrategy intToSimplexScaleStrategy(const int& value) {
  if (value == (int)SimplexScaleStrategy::OFF) return SimplexScaleStrategy::OFF;
  if (value == (int)SimplexScaleStrategy::CHOOSE)
    return SimplexScaleStrategy::CHOOSE;
  if (value == (int)SimplexScaleStrategy::HSOL)
    return SimplexScaleStrategy::HSOL;
  if (value == (int)SimplexScaleStrategy::HIGHS)
    return SimplexScaleStrategy::HIGHS;
  return SimplexScaleStrategy::DEFAULT;
}
SimplexCrashStrategy intToSimplexCrashStrategy(const int& value) {
  if (value == (int)SimplexCrashStrategy::LTSSF_K)
    return SimplexCrashStrategy::LTSSF_K;
  if (value == (int)SimplexCrashStrategy::LTSSF_PRI)
    return SimplexCrashStrategy::LTSSF_PRI;
  if (value == (int)SimplexCrashStrategy::BIXBY)
    return SimplexCrashStrategy::BIXBY;
  printf("Crash strategy option %d cannot be parsed yet!\n", value);
  return SimplexCrashStrategy::DEFAULT;
}
SimplexDualEdgeWeightStrategy intToSimplexDualEdgeWeightStrategy(
    const int& value) {
  if (value == (int)SimplexDualEdgeWeightStrategy::DANTZIG)
    return SimplexDualEdgeWeightStrategy::DANTZIG;
  if (value == (int)SimplexDualEdgeWeightStrategy::DEVEX)
    return SimplexDualEdgeWeightStrategy::DEVEX;
  if (value == (int)SimplexDualEdgeWeightStrategy::STEEPEST_EDGE)
    return SimplexDualEdgeWeightStrategy::STEEPEST_EDGE;
  if (value ==
      (int)SimplexDualEdgeWeightStrategy::STEEPEST_EDGE_TO_DEVEX_SWITCH)
    return SimplexDualEdgeWeightStrategy::STEEPEST_EDGE_TO_DEVEX_SWITCH;
  return SimplexDualEdgeWeightStrategy::DEFAULT;
}
SimplexPrimalEdgeWeightStrategy intToSimplexPrimalEdgeWeightStrategy(
    const int& value) {
  if (value == (int)SimplexPrimalEdgeWeightStrategy::DANTZIG)
    return SimplexPrimalEdgeWeightStrategy::DANTZIG;
  if (value == (int)SimplexPrimalEdgeWeightStrategy::DEVEX)
    return SimplexPrimalEdgeWeightStrategy::DEVEX;
  return SimplexPrimalEdgeWeightStrategy::DEFAULT;
}
SimplexPriceStrategy intToSimplexPriceStrategy(const int& value) {
  if (value == (int)SimplexPriceStrategy::COL) return SimplexPriceStrategy::COL;
  if (value == (int)SimplexPriceStrategy::ROW) return SimplexPriceStrategy::ROW;
  if (value == (int)SimplexPriceStrategy::ROW_SWITCH)
    return SimplexPriceStrategy::ROW_SWITCH;
  if (value == (int)SimplexPriceStrategy::ROW_SWITCH_COL_SWITCH)
    return SimplexPriceStrategy::ROW_SWITCH_COL_SWITCH;
  return SimplexPriceStrategy::DEFAULT;
}
