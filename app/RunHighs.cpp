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
#include "HighsTimer.h"
#include "HighsIO.h"
#include "HighsOptions.h"
#include "HighsRuntimeOptions.h"
#include "LoadProblem.h"

void HiGHSRun(const char *message = nullptr) {
  std::stringstream ss;
  ss << "Running HiGHS " << HIGHS_VERSION_MAJOR << "." << HIGHS_VERSION_MINOR
     << "." << HIGHS_VERSION_PATCH << " [date: " << HIGHS_COMPILATION_DATE
     << ", git hash: " << HIGHS_GITHASH << "]" << std::endl;

  HighsPrintMessage(ML_ALWAYS, ss.str().c_str());
  HighsPrintMessage(ML_ALWAYS, "Copyright (c) 2019 ERGO-Code under MIT licence terms\n\n");

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
};

int main(int argc, char **argv) {
  HiGHSRun();

  // Load user options.
  HighsOptions options;
  bool options_ok = loadOptions(argc, argv, options);

  if (!options_ok) return 0;

  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  if (read_status != HighsStatus::OK) {
    HighsPrintMessage(ML_ALWAYS, "Error loading file.\n");
    return (int)HighsStatus::LpError;
  } else {
    HighsPrintMessage(ML_MINIMAL, "LP       : %s\n",
                      lp.model_name_.c_str());
    HighsPrintMessage(ML_MINIMAL,
                      "Rows     : %d\nCols     : %d\nNonzeros : %d\n\n",
                      lp.numRow_, lp.numCol_, lp.Avalue_.size());
  }

  Highs highs;
  
  HighsStatus init_status = highs.initializeLp(lp);
  if (init_status != HighsStatus::OK) {
    HighsPrintMessage(ML_ALWAYS, "Error setting HighsLp.\n");
    return (int)HighsStatus::LpError;
  }

  HighsStatus run_status = highs.run();
  std::string statusname = HighsStatusToString(run_status);
  if (run_status != HighsStatus::OK &&
      run_status != HighsStatus::Optimal)
    HighsPrintMessage(ML_ALWAYS, "Highs status: %s\n", statusname.c_str());

  return 0;
}
