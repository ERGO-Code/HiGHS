#include "HighsSetup.h"

int main(int argc, char **argv) {
  // Load user options.
  HighsOptions options;
  HighsStatus init_status = loadOptions(argc, argv, options);

  if (init_status != HighsStatus::OK) {
    printHelp(argv[0]);
    return 0;
  }

  // todo:
  // Read Lp from a file.
  HighsLp lp;

  /*
    Status read_status = loadLpFromFile(options, lp);
    checkStatus(read_status);
  */

  Highs highs(options);
  HighsSolution solution;

  HighsStatus run_status = highs.run(lp, solution);
  checkStatus(run_status);

  return 0;
}