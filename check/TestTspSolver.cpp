#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
#include "lp_data/HighsCallback.h"

const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

HighsCallbackFunctionType userDefineLazyConstraints =
    [](int callback_type, const std::string& message,
       const HighsCallbackDataOut* data_out, HighsCallbackDataIn* data_in,
       void* user_callback_data) {
      if (dev_run) {
        printf("userDefineLazyConstraints:\n");
      }
    };

TEST_CASE("tsp-p01", "[highs_test_tsp_solver]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/p01.mps";
  const double optimal_obective_value = 263;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  //  MipData user_callback_data;
  highs.setCallback(userDefineLazyConstraints);  //, p_user_callback_data);
  printf("Calling highs.setCallback\n");
  highs.startCallback(kCallbackMipDefineLazyConstraints);
  highs.run();
  REQUIRE(highs.getObjectiveValue() == optimal_obective_value);
}
