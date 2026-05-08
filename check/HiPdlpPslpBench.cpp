#include <iostream>
#include <string>

#include "ExperimentalPslpHiPdlp.h"

namespace {

void printUsage(const char* executable) {
  std::cerr << "Usage: " << executable
            << " <model_file> --mode=plain|pslp\n";
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    printUsage(argv[0]);
    return 1;
  }

  const std::string model_file = argv[1];
  const std::string mode_flag = argv[2];
  const std::string mode_prefix = "--mode=";
  if (mode_flag.compare(0, mode_prefix.size(), mode_prefix) != 0) {
    printUsage(argv[0]);
    return 1;
  }

  ExperimentalHiPdlpPslpRun result;
  std::string error_message;
  const HighsStatus status = runExperimentalHiPdlpPslpBenchmark(
      model_file, mode_flag.substr(mode_prefix.size()), result, error_message);
  if (status == HighsStatus::kError) {
    std::cerr << error_message << "\n";
    return 1;
  }

  std::cout << experimentalHiPdlpPslpCsvLine(result) << "\n";
  return 0;
}
