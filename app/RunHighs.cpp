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
#include "HighsSetup.h"
#include "HighsTimer.h"
#include "LoadProblem.h"

int main(int argc, char **argv) {
  // Initialise timer
  HighsTimer timer;
  double start_time = timer.getWallTime();

  HiGHSRun();

  // Load user options.
  HighsOptions options;
  HighsStatus init_status = loadOptions(argc, argv, options);

  if (init_status != HighsStatus::OK) return 0;

  HighsLp lp;
  HighsInputStatus read_status = loadLpFromFile(options, lp);
  if (read_status != HighsInputStatus::OK) {
    std::string message = "Error loading file: status = " +
                          HighsInputStatusToString(read_status) + "\n";
    HighsLogMessage(HighsMessageType::INFO, message.c_str());
    return (int)HighsStatus::LpError;
  }

  Highs highs(options);
  HighsSolution solution;

  HighsStatus run_status = highs.run(lp, solution);

  double end_time = timer.getWallTime();
  HighsLogMessage(HighsMessageType::INFO, "HiGHS run ended after %12g seconds",
                  end_time - start_time);
  return 0;
}
