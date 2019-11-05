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

  // Load user options.
  HighsOptions options;
  bool options_ok = loadOptions(argc, argv, options);
  if (!options_ok) return 0;

  bool force_options_file = false;  // true;// 
  if (force_options_file && options.options_file.size() > 0) {
    printf(
        "In main: set options.options_file = options_file so vscode can be "
        "used to debug\n");
    options.options_file = "options_file";
    if (!loadOptionsFromFile(options)) {
      printf("In main: fail return from loadOptionsFromFile\n");
      return (int)HighsStatus::Error;
    }
  }
  if (options.run_as_hsol) setHsolOptions(options);

  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  if (read_status != HighsStatus::OK) {
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

  Highs highs;
  //  highs.writeHighsOptions("HiGHS.set");

  HighsStatus init_status = highs.passModel(lp);
  if (init_status != HighsStatus::OK) {
    HighsPrintMessage(ML_ALWAYS, "Error setting HighsLp.\n");
    return (int)HighsStatus::Error;
  }

  HighsStatus run_status;
  /*
  run_status = highs.writeModel("write.mps"); 
  if (run_status != HighsStatus::OK) printf("Error return from highs.writeModel\n");
  */

  highs.options_ = options;
  run_status = highs.run();
  std::string statusname = HighsStatusToString(run_status);

  if (run_status != HighsStatus::OK) {
    HighsPrintMessage(ML_ALWAYS, "HiGHS status: %s\n", statusname.c_str());
  } else {
    HighsModelStatus model_status = highs.getModelStatus();
    int iteration_count = highs.getIterationCount();
    std::stringstream message;
    message << std::endl;
    message << "Run status : " << highs.highsModelStatusToString(model_status) << std::endl;
    message << "Iterations : " << iteration_count << std::endl;

    if (model_status == HighsModelStatus::OPTIMAL) {
      double dual_objective_value = highs.getObjectiveValue();
      message << "Objective  : " << std::scientific << dual_objective_value << std::endl;
    }
    message << std::endl;
    std::cout << message.str();
  }
  return (int)run_status;
}
