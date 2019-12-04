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

void HiGHSRun(const char* message = nullptr) {
  std::stringstream ss;
  ss << "Running HiGHS " << HIGHS_VERSION_MAJOR << "." << HIGHS_VERSION_MINOR
     << "." << HIGHS_VERSION_PATCH << " [date: " << HIGHS_COMPILATION_DATE
     << ", git hash: " << HIGHS_GITHASH << "]" << std::endl;

  HighsPrintMessage(ML_ALWAYS, ss.str().c_str());
  HighsPrintMessage(ML_ALWAYS,
                    "Copyright (c) 2019 ERGO-Code under MIT licence terms\n\n");

  if (message != nullptr) {
  }

#ifdef HiGHSDEV
  // Report on preprocessing macros
  if (message != nullptr) {
    std::cout << "In " << message << std::endl;
  }
#ifdef OPENMP
  std::cout << "OPENMP           is     defined" << std::endl;
#else
  std::cout << "OPENMP           is not defined" << std::endl;
#endif

#ifdef SCIP_DEV
  std::cout << "SCIP_DEV         is     defined" << std::endl;
#else
  std::cout << "SCIP_DEV         is not defined" << std::endl;
#endif

#ifdef HiGHSDEV
  std::cout << "HiGHSDEV         is     defined" << std::endl;
#else
  std::cout << "HiGHSDEV         is not defined" << std::endl;
#endif
  std::cout << "Built with CMAKE_BUILD_TYPE=" << CMAKE_BUILD_TYPE << std::endl;

#endif
}

int main(int argc, char** argv) {
  HiGHSRun();
  HighsStatus return_status;

  Highs highs;
  //  highs.writeHighsOptions("HiGHS.set");

  // Load user options.
  HighsOptions options;
  bool options_ok = loadOptions(argc, argv, options);
  if (!options_ok) return 0;

  bool force_options_file = false;  // true;//
  if (force_options_file) {
    printf(
        "In main: set options.options_file = Options.set so vscode can be "
        "used to debug\n");
    options.options_file = "Options.set";
    if (!loadOptionsFromFile(options)) {
      printf("In main: fail return from loadOptionsFromFile\n");
      return (int)HighsStatus::Error;
    }
  }

  return_status = highs.passHighsOptions(options);
  if (return_status != HighsStatus::OK) {
    if (return_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      printf("HighsStatus::Warning return from passHighsOptions\n");
#endif
    } else {
      printf("In main: fail return from passHighsOptions\n");
      return (int)return_status;
    }
  }

  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  if (read_status == HighsStatus::Error) {
    std::cout << "Error loading file" << std::endl;
    return (int)HighsStatus::Error;
  } else {
    std::stringstream message;
    message << "LP       : " << lp.model_name_.c_str() << std::endl;
    message << "Rows     : " << lp.numRow_ << std::endl;
    message << "Cols     : " << lp.numCol_ << std::endl;
    message << "Nonzeros : " << lp.Avalue_.size() << std::endl;
    if (lp.numInt_)
      message << "Integer  : " << lp.numInt_ << std::endl;
    message  << std::endl;
    std::cout << message.str();
  }

  HighsStatus init_status = highs.passModel(lp);
  if (init_status != HighsStatus::OK) {
    if (init_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      printf("HighsStatus::Warning return setting HighsLp\n");
#endif
    } else {
      HighsPrintMessage(ML_ALWAYS, "Error setting HighsLp\n");
      return (int)HighsStatus::Error;
    }
  }

  /*
  HighsStatus write_status;
  write_status = highs.writeModel("write.mps"); 
  if (write_status != HighsStatus::OK) {
    if (write_status == HighsStatus::Warning) {
#ifdef HiGHSDEV
      printf("HighsStatus::Warning return from highs.writeModel\n");
#endif
    } else {
      printf("Error return from highs.writeModel\n");
    }
  }
  */

  //  highs.options_ = options;
  HighsStatus run_status = highs.run();
  std::string statusname = HighsStatusToString(run_status);

  if (run_status == HighsStatus::Error) {
    HighsPrintMessage(ML_ALWAYS, "HiGHS status: %s\n", statusname.c_str());
  } else {
    std::stringstream message;
    message << std::endl;
    HighsModelStatus model_status = highs.getModelStatus();
    HighsModelStatus scaled_model_status = highs.getModelStatus(true);
    HighsInfo highs_info = highs.getHighsInfo();
    if (model_status != scaled_model_status) {
      if (scaled_model_status == HighsModelStatus::OPTIMAL) {
	// The scaled model has been solved to optimality, but not the
	// unscaled model, flag this up, but report the scaled model
	// status
	double max_primal_infeasibility = highs_info.max_primal_infeasibility;
	double max_dual_infeasibility = highs_info.max_dual_infeasibility;
	//	message << std::setprecision(9);
	message << "Primal infeasibility: " << max_primal_infeasibility << std::endl;
	message << "Dual   infeasibility: " << max_dual_infeasibility << std::endl;
	model_status = scaled_model_status;
      }
    }
    message << "Model   status      : " << highs.highsModelStatusToString(model_status) << std::endl;
    message << "Primal  status      : " << highs.highsPrimalDualStatusToString(highs_info.primal_status) << std::endl;
    message << "Dual    status      : " << highs.highsPrimalDualStatusToString(highs_info.dual_status) << std::endl;
    message << "Simplex   iterations: " << highs_info.simplex_iteration_count << std::endl;
    if (highs_info.ipm_iteration_count)
      message << "IPM       iterations: " << highs_info.ipm_iteration_count << std::endl;
    if (highs_info.crossover_iteration_count)
      message << "Crossover iterations: " << highs_info.crossover_iteration_count << std::endl;
    if (model_status == HighsModelStatus::OPTIMAL) {
      double objective_function_value;
      highs.getHighsInfoValue("objective_function_value", objective_function_value);
      message << "Objective value     : " << std::scientific << objective_function_value << std::endl;
    }
    message << std::endl;
    std::cout << message.str();

  // Possibly write the solution to a file
  if (options.write_solution_to_file)
    highs.writeSolution(options.solution_file, options.write_solution_pretty);
  }
  /*
  highs.writeHighsInfo("HighsInfo.dat");
  highs.writeHighsInfo("HighsInfo.html");
  */
  return (int)run_status;
}
