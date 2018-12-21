#include "HighsSetup.h"
#include "HighsTimer.h"
#include "LoadProblem.h"

int main(int argc, char **argv) {
  // Initialise timer
  HighsTimer timer;
  int loadClock = timer.clockDef("Load", " Ld");
  int runClock = timer.clockDef("Run", "Run");
  //  timer.reset();
  HiGHSRun();

  // Load user options.
  HighsOptions options;
  HighsStatus init_status = loadOptions(argc, argv, options);

  if (init_status != HighsStatus::OK) return 0;

  HighsLp lp;
  timer.start(loadClock);
  HighsInputStatus read_status = loadLpFromFile(options, lp);
  timer.stop(loadClock);
  if (read_status != HighsInputStatus::OK) {
    HighsLogMessage(HighsMessageType::INFO, "Error when parsing file\n");
    return (int)HighsStatus::LpError;
  }

  Highs highs(options);
  HighsSolution solution;

  timer.start(runClock);
  HighsStatus run_status = highs.run(lp, solution);
  timer.stop(runClock);

  checkStatus(run_status);

  // Report times
  std::vector<int> clockList{loadClock, runClock};
  timer.report(clockList);
  return 0;
}
