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
    HighsLogMessage(HighsMessageType::INFO, "Error when parsing file\n");
    return (int)HighsStatus::LpError;
  }

  Highs highs(options);
  HighsSolution solution;

  HighsStatus run_status = highs.run(lp, solution);

  double end_time = timer.getWallTime();
  HighsLogMessage(HighsMessageType::INFO, "HiGHS run ended after %12g seconds\n", end_time-start_time);
  
  checkStatus(run_status);

  
  return 0;
}
