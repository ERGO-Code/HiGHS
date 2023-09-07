#include <cstdio>
#include <cstring>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = true;

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

static void userCallback(const int callback_type, const char* message,
                         const HighsCallbackDataOut* data_out,
                         HighsCallbackDataIn* data_in,
                         void* user_callback_data) {
  // Extract local_callback_data from user_callback_data unless it
  // is nullptr
  const int local_callback_data =
      user_callback_data
          ? static_cast<int>(reinterpret_cast<intptr_t>(user_callback_data))
          : kUserCallbackNoData;
  if (user_callback_data) {
    REQUIRE(local_callback_data == kUserCallbackData);
  } else {
    REQUIRE(local_callback_data == kUserCallbackNoData);
  }
  if (dev_run) {
    if (callback_type == kHighsCallbackLogging) {
      printf("userCallback(type %2d; data %2d): %s", callback_type,
             local_callback_data, message);
    } else if (callback_type == kHighsCallbackInterrupt) {
      printf("userCallback(type %2d; data %2d): %s with iteration count = "
	     "%d\n",
	     callback_type, local_callback_data, message,
	     data_out->simplex_iteration_count);
      data_in->user_interrupt =	data_out->simplex_iteration_count > 30;
    } else if (callback_type == kHighsCallbackMipImprovingSolution) {
      printf("userCallback(type %2d; data %2d): %s with objective %g and solution[0] = %g\n",
	     callback_type, local_callback_data, message,
	     data_out->objective, data_out->col_value[0]);
    }
  }
}

TEST_CASE("my-callback-logging", "[highs-callback]") {
  bool output_flag = true;
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
  // Uses userCallback to start logging lines with
  // "userCallback(kUserCallbackData): " since
  // Highs::setCallback has second argument p_user_callback_data
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  int user_callback_data = kUserCallbackData;
  void* p_user_callback_data =
      reinterpret_cast<void*>(static_cast<intptr_t>(user_callback_data));
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userCallback, p_user_callback_data);
  highs.startCallback(kHighsCallbackLogging);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-interrupt", "[highs-callback]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setCallback(userCallback);
  highs.startCallback(kHighsCallbackInterrupt);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("highs-callback-mip-improving", "[highs-callback]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/egout.mps";
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.setCallback(userCallback);
  highs.startCallback(kHighsCallbackMipImprovingSolution);
  highs.readModel(filename);
  highs.run();
}
