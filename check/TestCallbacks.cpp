#include "Highs.h"
#include "catch.hpp"
const bool dev_run = true;

// Callback that provides user logging

// const int highs_callback_type
//  , const HighsCallbackDataOut& callback_data_out,
//  HighsCallbackDataIn& callback_data_in
static void userHighsCallback(HighsLogType type, const char* message,
                              void* deprecated) {
  if (dev_run)
    printf("userHighsCallback:%s", message);  // %d", highs_callback_type);
}

TEST_CASE("highs-run-callback", "[callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.setHighsCallback(userHighsCallback);
  highs.readModel(filename);
  highs.run();
}
