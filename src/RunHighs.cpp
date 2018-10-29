#include "HighsSetup.h"

int main(int argc, char **argv) {
  // Load user options.
  Options options;
  Status init_status = loadOptions(argc, argv, options);

  if (init_status == Status::InputError) {
    printHelp(argv[0]);
    return 0;
  }
  checkStatus(init_status);

  // todo:
  // Read LpData from a file.
  LpData lp;
/*
  Status read_status = loadLpFromFile(options, lp);
  checkStatus(read_status);
*/

  Highs highs(options);
  Solution solution;

  Status run_status = highs.run(lp, solution);
  checkStatus(run_status);

  return 0;
}