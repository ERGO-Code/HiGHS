#include "Highs.h"
#include "catch.hpp"
const bool dev_run = true;

static void userHighsCallback(const int highs_callback_type,
                              const HighsCallbackDataOut& callback_data_out,
                              HighsCallbackDataIn& callback_data_in) {
  if (dev_run)
    printf("userHighsCallback: %d with simplex iteration count = %d\n",
	   highs_callback_type,
	   callback_data_out.simplex_iteration_count);
  callback_data_in.user_interrupt = callback_data_out.simplex_iteration_count > 30;
}

TEST_CASE("highs-callback-interrupt", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.setHighsCallback(userHighsCallback);
  highs.readModel(filename);
  highs.run();
}
