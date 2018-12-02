#include "HighsSetup.h"
#include "LoadProblem.h"

int main(int argc, char **argv) {
  // Load user options.
  HighsOptions options;
 // HighsStatus init_status = loadOptions(argc, argv, options);

  // Use to replace old HighsOptions.
  HighsStringOptions options_;
  loadOptions(argc, argv, options_);

  HighsLp lp;
  HighsInputStatus read_status = loadLpFromFile(options, lp);
  if (read_status != HighsInputStatus::OK) {
    return (int) HighsStatus::LpError;
  }

  Highs highs(options);
  HighsSolution solution;

  HighsStatus run_status = highs.run(lp, solution);
  checkStatus(run_status);

  return 0;
}