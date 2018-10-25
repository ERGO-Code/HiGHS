#include "HApp.cpp"

int main_(int argc, char **argv) {
  // Load user options.
  Options options;
  Status status = loadOptions(options);
  checkStatus(status);
  
  // Read LpData from a file.
  LpData lp;
  status = loadLpDataFromFile(options_.filename, lp);
  checkStatus(status);

  Highs highs(options);
  highs.run();


 

  return 0;
}