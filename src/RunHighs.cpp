#include "HApp.cpp"

int main_(int argc, char **argv) {
  // Load user options.
  Options options;
  Status init_status = loadOptions(options);

  if (init_status == Status::InputError) {
    printHelp(argv[0]);
    return 0;
  }
  checkStatus(init_status);
  
  // Read LpData from a file.
  LpData lp;
  Status read_status = loadLpDataFromFile(options_.filename, lp);
  checkStatus(read_status);

  Highs highs(options);
  Solution solution(lp.getDimensions());

  Status run_status = highs.run();
  checkStatus(run_status);


  
  return 0;
}