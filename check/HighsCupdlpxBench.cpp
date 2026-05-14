#include <iostream>
#include <string>

#include "ExperimentalHighsCupdlpx.h"

namespace {

void printUsage(const char* executable) {
  std::cerr << "Usage: " << executable << " <model_file>\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    printUsage(argv[0]);
    return 1;
  }

  ExperimentalHighsCupdlpxRun result;
  std::string error_message;
  const HighsStatus status = runExperimentalHighsCupdlpxBenchmark(
      argv[1], result, error_message);
  if (status == HighsStatus::kError) {
    std::cerr << error_message << "\n";
    return 1;
  }

  std::cout << experimentalHighsCupdlpxCsvLine(result) << "\n";
  return 0;
}
