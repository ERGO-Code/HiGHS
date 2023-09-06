#include "Highs.h"
#include "catch.hpp"
const bool dev_run = true;
const HighsInt kLogUserCallbackNoData = -1;
const HighsInt kLogUserCallbackData = 99;

static void userHighsCallback(const int highs_callback_type,
                              const char* message, void* user_callback_data,
                              const HighsCallbackDataOut& callback_data_out,
                              HighsCallbackDataIn& callback_data_in) {
  // Extract local_callback_data from user_callback_data unless it
  // is nullptr
  const int local_callback_data =
      user_callback_data
          ? static_cast<int>(reinterpret_cast<intptr_t>(user_callback_data))
          : kLogUserCallbackNoData;
  if (user_callback_data) {
    REQUIRE(local_callback_data == kLogUserCallbackData);
  } else {
    REQUIRE(local_callback_data == kLogUserCallbackNoData);
  }
  if (dev_run) {
    if (highs_callback_type == kHighsCallbackLogging) {
      printf("userHighsCallback(type %2d; data %2d): %s", highs_callback_type,
             local_callback_data, message);
    } else if (highs_callback_type == kHighsCallbackInterrupt) {
      printf(
          "userHighsCallback(type %2d; data %2d): %s with iteration count = "
          "%d\n",
          highs_callback_type, local_callback_data, message,
          callback_data_out.simplex_iteration_count);
      callback_data_in.user_interrupt =
          callback_data_out.simplex_iteration_count > 30;
    }
  }
}

TEST_CASE("highs-callback-logging", "[highs-callback]") {
  // Uses userHighsCallback to start logging lines with
  // "userHighsCallback(kLogUserCallbackData): " since
  // Highs::setHighsCallback has second argument
  // p_highs_user_callback_data
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  int highs_user_callback_data = kLogUserCallbackData;
  void* p_highs_user_callback_data =
      reinterpret_cast<void*>(static_cast<intptr_t>(highs_user_callback_data));
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setHighsCallback(userHighsCallback, p_highs_user_callback_data);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-interrupt", "[highs-callback]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setHighsCallback(userHighsCallback);
  highs.readModel(filename);
  highs.run();
}
