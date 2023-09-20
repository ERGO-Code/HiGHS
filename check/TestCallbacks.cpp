#include <cstdio>
#include <cstring>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

const double egout_optimal_objective = 568.1007;
const double egout_objective_target = 610;
const HighsInt adlittle_ipm_iteration_limit = 5;
const HighsInt adlittle_simplex_iteration_limit = 30;

const HighsInt kLogBufferSize = kIoBufferSize;
const HighsInt kUserCallbackNoData = -1;
const HighsInt kUserCallbackData = 99;

char printed_log[kLogBufferSize];

using std::memset;
using std::strcmp;
using std::strcpy;
using std::strlen;
using std::strncmp;
using std::strstr;

// Callback that saves message for comparison
static void myLogCallback(const int callback_type, const char* message,
                          const HighsCallbackDataOut* data_out,
                          HighsCallbackDataIn* data_in,
                          void* user_callback_data) {
  strcpy(printed_log, message);
}

static void userInterruptCallback(const int callback_type, const char* message,
                                  const HighsCallbackDataOut* data_out,
                                  HighsCallbackDataIn* data_in,
                                  void* user_callback_data) {
  // Extract local_callback_data from user_callback_data unless it
  // is nullptr
  if (callback_type == kCallbackMipImprovingSolution) {
    // Use local_callback_data to maintain the objective value from
    // the previous callback
    assert(user_callback_data);
    // Extract the double value pointed to from void* user_callback_data
    const double local_callback_data = *(double*)user_callback_data;
    if (dev_run)
      printf(
          "userCallback(type %2d; data %11.4g): %s with objective %g and "
          "solution[0] = %g\n",
          callback_type, local_callback_data, message,
          data_out->objective_function_value, data_out->mip_solution[0]);
    REQUIRE(local_callback_data >= data_out->objective_function_value);
    // Update the double value pointed to from void* user_callback_data
    *(double*)user_callback_data = data_out->objective_function_value;
  } else {
    const int local_callback_data =
        user_callback_data
            ? static_cast<int>(reinterpret_cast<intptr_t>(user_callback_data))
            : kUserCallbackNoData;
    if (user_callback_data) {
      REQUIRE(local_callback_data == kUserCallbackData);
    } else {
      REQUIRE(local_callback_data == kUserCallbackNoData);
    }
    if (callback_type == kCallbackLogging) {
      if (dev_run)
        printf("userInterruptCallback(type %2d; data %2d): %s", callback_type,
               local_callback_data, message);
    } else if (callback_type == kCallbackSimplexInterrupt) {
      if (dev_run)
        printf(
            "userInterruptCallback(type %2d; data %2d): %s with iteration "
            "count = "
            "%d\n",
            callback_type, local_callback_data, message,
            int(data_out->simplex_iteration_count));
      data_in->user_interrupt =
          data_out->simplex_iteration_count > adlittle_simplex_iteration_limit;
    } else if (callback_type == kCallbackIpmInterrupt) {
      if (dev_run)
        printf(
            "userInterruptCallback(type %2d; data %2d): %s with iteration "
            "count = "
            "%d\n",
            callback_type, local_callback_data, message,
            int(data_out->ipm_iteration_count));
      data_in->user_interrupt =
          data_out->ipm_iteration_count > adlittle_ipm_iteration_limit;
      data_in->user_interrupt = true;
    } else if (callback_type == kCallbackMipInterrupt) {
      if (dev_run)
        printf(
            "userInterruptCallback(type %2d; data %2d): %s with Bounds "
            "(%11.4g, %11.4g); Gap = %11.4g; Objective = "
            "%g\n",
            callback_type, local_callback_data, message,
            data_out->mip_dual_bound, data_out->mip_primal_bound,
            data_out->mip_gap, data_out->objective_function_value);
      data_in->user_interrupt =
          data_out->objective_function_value < egout_objective_target;
    }
  }
}

static void userDataCallback(const int callback_type, const char* message,
                             const HighsCallbackDataOut* data_out,
                             HighsCallbackDataIn* data_in,
                             void* user_callback_data) {
  assert(callback_type == kCallbackMipInterrupt ||
         callback_type == kCallbackMipLogging ||
         callback_type == kCallbackMipImprovingSolution);
  if (dev_run)
    printf("userDataCallback: Node count = %" PRId64
           "; Time = %6.2f; "
           "Bounds (%11.4g, %11.4g); Gap = %11.4g; Objective = %11.4g: %s\n",
           data_out->mip_node_count, data_out->running_time,
           data_out->mip_dual_bound, data_out->mip_primal_bound,
           data_out->mip_gap, data_out->objective_function_value, message);
}

TEST_CASE("my-callback-logging", "[highs-callback]") {
  bool output_flag = true;  // Still runs quietly
  bool log_to_console = false;
  HighsInt log_dev_level = kHighsLogDevLevelInfo;
  HighsLogOptions log_options;
  log_options.clear();
  log_options.log_stream = stdout;
  log_options.output_flag = &output_flag;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  log_options.user_callback = myLogCallback;
  log_options.user_callback_active = true;

  highsLogDev(log_options, HighsLogType::kInfo, "Hi %s!", "HiGHS");
  if (dev_run) printf("Log callback yields \"%s\"\n", printed_log);
  REQUIRE(strcmp(printed_log, "Hi HiGHS!") == 0);

  // Check that nothing is printed if the type is VERBOSE when
  // log_dev_level is kHighsLogDevLevelInfo;
  *printed_log = '\0';
  highsLogDev(log_options, HighsLogType::kVerbose, "Hi %s!", "HiGHS");
  REQUIRE(*printed_log == '\0');

  {
    char long_message[sizeof(printed_log)];
    memset(long_message, 'H', sizeof(long_message));
    long_message[sizeof(long_message) - 2] = '\0';
    long_message[sizeof(long_message) - 1] = '\n';
    highsLogDev(log_options, HighsLogType::kInfo, long_message);
    if (dev_run) printf("Log callback yields \"%s\"\n", printed_log);
    REQUIRE(strncmp(printed_log, "HHHH", 4) == 0);
    REQUIRE(strlen(printed_log) <= sizeof(printed_log));
  }

  highsLogUser(log_options, HighsLogType::kInfo, "Hello %s!\n", "HiGHS");
  REQUIRE(strlen(printed_log) > 9);
  REQUIRE(strcmp(printed_log, "Hello HiGHS!\n") == 0);

  {
    char long_message[sizeof(printed_log)];
    memset(long_message, 'H', sizeof(long_message));
    long_message[sizeof(long_message) - 2] = '\0';
    long_message[sizeof(long_message) - 1] = '\n';
    highsLogUser(log_options, HighsLogType::kWarning, long_message);
    if (dev_run) printf("Log callback yields \"%s\"\n", printed_log);
    REQUIRE(strstr(printed_log, "HHHH") != nullptr);
    REQUIRE(strlen(printed_log) <= sizeof(printed_log));
  }
}

TEST_CASE("highs-callback-logging", "[highs-callback]") {
  // Uses userInterruptCallback to start logging lines with
  // "userInterruptCallback(kUserCallbackData): " since
  // Highs::setCallback has second argument p_user_callback_data
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  int user_callback_data = kUserCallbackData;
  void* p_user_callback_data =
      reinterpret_cast<void*>(static_cast<intptr_t>(user_callback_data));
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userInterruptCallback, p_user_callback_data);
  highs.startCallback(kCallbackLogging);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-simplex-interrupt", "[highs-callback]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userInterruptCallback);
  highs.startCallback(kCallbackSimplexInterrupt);
  highs.readModel(filename);
  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInterrupt);
  REQUIRE(highs.getInfo().simplex_iteration_count > adlittle_simplex_iteration_limit);
}

TEST_CASE("highs-callback-ipm-interrupt", "[highs-callback]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userInterruptCallback);
  highs.startCallback(kCallbackIpmInterrupt);
  highs.readModel(filename);
  highs.setOptionValue("solver", kIpmString);
  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInterrupt);
  REQUIRE(highs.getInfo().ipm_iteration_count > adlittle_ipm_iteration_limit);
}

TEST_CASE("highs-callback-mip-interrupt", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setCallback(userInterruptCallback);
  highs.startCallback(kCallbackMipInterrupt);
  highs.readModel(filename);
  HighsStatus status = highs.run();
  REQUIRE(status == HighsStatus::kWarning);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInterrupt);
  REQUIRE(highs.getInfo().objective_function_value > egout_optimal_objective);
}

TEST_CASE("highs-callback-mip-improving", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  double user_callback_data = kHighsInf;
  void* p_user_callback_data = (void*)(&user_callback_data);
  highs.setCallback(userInterruptCallback, p_user_callback_data);
  highs.startCallback(kCallbackMipImprovingSolution);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-mip-data", "[highs-callback]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setCallback(userDataCallback);
  highs.startCallback(kCallbackMipImprovingSolution);
  highs.startCallback(kCallbackMipLogging);
  highs.readModel(filename);
  highs.run();
}
