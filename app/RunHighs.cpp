/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2019 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file ../app/RunHighs.cpp
 * @brief HiGHS main
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#include "Highs.h"
#include "HighsIO.h"
#include "HighsOptions.h"
#include "HighsRuntimeOptions.h"
#include "HighsTimer.h"
#include "LoadProblem.h"

void HiGHSRun(FILE* output, const int message_level, const char* message = nullptr) {
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "Running HiGHS %d.%d.%d [date: %s, git hash: %s]\n",
		    HIGHS_VERSION_MAJOR,
		    HIGHS_VERSION_MINOR,
		    HIGHS_VERSION_PATCH,
		    HIGHS_COMPILATION_DATE,
		    HIGHS_GITHASH);
  HighsPrintMessage(output, message_level, ML_ALWAYS,
                    "Copyright (c) 2019 ERGO-Code under MIT licence terms\n\n");
#ifdef HiGHSDEV
  // Report on preprocessing macros
  if (message != nullptr) {
    HighsPrintMessage(output, message_level, ML_ALWAYS,
		      "In %s\n", message);
  }
#ifdef OPENMP
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "OPENMP           is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "OPENMP           is not defined\n");
#endif

#ifdef SCIP_DEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "SCIP_DEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "SCIP_DEV         is not defined\n");
#endif

#ifdef HiGHSDEV
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "HiGHSDEV         is     defined\n");
#else
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "HiGHSDEV         is not defined\n");
#endif
  HighsPrintMessage(output, message_level, ML_ALWAYS,
		    "Built with CMAKE_BUILD_TYPE=%s\n", CMAKE_BUILD_TYPE);
#endif
}

int main(int argc, char** argv) {
  std::stringstream ss;
  HiGHSRun(stdout, ML_ALWAYS);
  HighsStatus return_status;

  Highs highs;
  //  highs.writeHighsOptions("HiGHS.set");

  // Load user options.
  HighsOptions options;
  bool options_ok = loadOptions(argc, argv, options);
  if (!options_ok) return 0;

  FILE* output;
  int message_level;
  output = options.output;
  message_level = options.message_level;

  bool run_quiet = false;//true; // 
  if (run_quiet) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, 
		      "In main: running highs.run() quietly\n");
  }
  bool force_options_file = false;  // true;//
  if (force_options_file) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, 
		      "In main: set options.options_file = Options.set so vscode can be "
		      "used to debug\n");
    options.options_file = "Options.set";
    if (!loadOptionsFromFile(options)) {
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"In main: fail return from loadOptionsFromFile\n");
      return (int)HighsStatus::Error;
    }
  }

  output = options.output;
  message_level = options.message_level;

  return_status = highs.passHighsOptions(options);
  if (return_status != HighsStatus::OK) {
    if (return_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"HighsStatus::Warning return from passHighsOptions\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS, 
			"In main: fail return from passHighsOptions\n");
      return (int)return_status;
    }
  }

  if (run_quiet) {
    highs.setHighsLogfile(NULL);
    highs.setHighsOutput(NULL);
  }

  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  if (read_status == HighsStatus::Error) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "Error loading file\n");
    return (int)HighsStatus::Error;
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS,
		      "LP       : %s\n", lp.model_name_.c_str());
    HighsPrintMessage(output, message_level, ML_ALWAYS,
		      "Rows     : %d\n", lp.numRow_);
    HighsPrintMessage(output, message_level, ML_ALWAYS,
		      "Cols     : %d\n", lp.numCol_);
    HighsPrintMessage(output, message_level, ML_ALWAYS,
		      "Nonzeros : %d\n", lp.Avalue_.size());
    if (lp.numInt_)
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"Integer  : %d\n", lp.numInt_);
    HighsPrintMessage(output, message_level, ML_ALWAYS, "\n");
  }

  HighsStatus init_status = highs.passModel(lp);
  if (init_status != HighsStatus::OK) {
    if (init_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"HighsStatus::Warning return setting HighsLp\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS, "Error setting HighsLp\n");
      return (int)HighsStatus::Error;
    }
  }

  /*
  HighsStatus write_status;
  write_status = highs.writeModel("write.mps"); 
  if (write_status != HighsStatus::OK) {
    if (write_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      HighsPrintMessage(output, message_level, ML_ALWAYS,
                        "HighsStatus::Warning return from highs.writeModel\n");
#endif
    } else {
      HighsPrintMessage(output, message_level, ML_ALWAYS, 
                        "Error return from highs.writeModel\n");
    }
  }
  */

  // Write all the options to an options file
  //  highs.writeHighsOptions("Highs.set", false);
  // Write all the options as HTML
  //  highs.writeHighsOptions("Highs.html", false);
  // Possibly report options settings
  highs.writeHighsOptions("");  //, false);

  if (run_quiet) HighsPrintMessage(output, message_level, ML_ALWAYS, "Before calling highs.run()\n");
  HighsStatus run_status = highs.run();
  if (run_quiet) HighsPrintMessage(output, message_level, ML_ALWAYS, "After calling highs.run()\n");
  std::string statusname = HighsStatusToString(run_status);

  if (run_status == HighsStatus::Error) {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "HiGHS status: %s\n", statusname.c_str());
  } else {
    HighsPrintMessage(output, message_level, ML_ALWAYS, "\n");
    HighsModelStatus model_status = highs.getModelStatus();
    HighsModelStatus scaled_model_status = highs.getModelStatus(true);
    HighsInfo highs_info = highs.getHighsInfo();
    if (model_status != scaled_model_status) {
      if (scaled_model_status == HighsModelStatus::OPTIMAL) {
	// The scaled model has been solved to optimality, but not the
	// unscaled model, flag this up, but report the scaled model
	// status
	HighsPrintMessage(output, message_level, ML_ALWAYS,
			  "Primal infeasibility: %10.3e (%d)\n",
			  highs_info.max_primal_infeasibility,
			  highs_info.num_primal_infeasibilities);
	HighsPrintMessage(output, message_level, ML_ALWAYS,
			  "Dual   infeasibility: %10.3e (%d)\n",
			  highs_info.max_dual_infeasibility,
			  highs_info.num_dual_infeasibilities);
	model_status = scaled_model_status;
      }
    }
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"Model   status      : %s\n",
			highs.highsModelStatusToString(model_status).c_str());
      /*
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"Primal  status      : %s\n",
			highs.highsPrimalDualStatusToString(highs_info.primal_status).c_str());
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"Dual    status      : %s\n",
			highs.highsPrimalDualStatusToString(highs_info.dual_status).c_str());
      */
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"Simplex   iterations: %d\n",
			highs_info.simplex_iteration_count);
    if (highs_info.ipm_iteration_count)
        HighsPrintMessage(output, message_level, ML_ALWAYS,
			  "IPM       iterations: %d\n",
			  highs_info.ipm_iteration_count);
    if (highs_info.crossover_iteration_count)
        HighsPrintMessage(output, message_level, ML_ALWAYS,
			  "Crossover iterations: %d\n",
			  highs_info.crossover_iteration_count);
    if (model_status == HighsModelStatus::OPTIMAL) {
      double objective_function_value;
      highs.getHighsInfoValue("objective_function_value", objective_function_value);
      HighsPrintMessage(output, message_level, ML_ALWAYS,
			"Objective value     : %13.6e\n",
			objective_function_value);
    }

  // Possibly write the solution to a file
  if (options.write_solution_to_file)
    highs.writeSolution(options.solution_file, options.write_solution_pretty);
  }
  /*
  highs.writeSolution("", true);
  highs.writeSolution("", false);
  highs.writeHighsInfo("");
  highs.writeHighsInfo("HighsInfo.html");
  */
  return (int)run_status;
}
