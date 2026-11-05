// Main file for unit testing with catch. We manually create the main function
// and initialize the HiGHS external API. Do not define CATCH_CONFIG_MAIN
// anywhere else. See correct-print-options in TestSetup.cpp.
#define CATCH_CONFIG_RUNNER
#include "Highs.h"
#include "catch.hpp"

int main(int argc, char* argv[]) {
  Highs::initializeExtras();  // ignore if not available, will throw
                              // errors in tests if needed
  return Catch::Session().run(argc, argv);
}