#include <cstdio>
#include <cstring>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

const HighsInt kLogBufferSize = kIoBufferSize;
const HighsInt kLogUserCallbackNoData = -1;
const HighsInt kLogUserCallbackData = 99;

char alt_printed_log[kLogBufferSize];

using std::memset;
using std::strcmp;
using std::strcpy;
using std::strlen;
using std::strncmp;
using std::strstr;

// Callback that saves message for comparison
static void myLogCallback(HighsLogType type, const char* message,
                          void* user_log_callback_data) {
  strcpy(alt_printed_log, message);
}

TEST_CASE("run-callback", "[highs_io]") {
  // Uses userLogCallback to start logging lines with
  // "userLogCallback(kLogUserCallbackNoData): " since
  // Highs::setLogCallback has no second argument so
  // user_log_callback_data will be nullptr
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  // highs.setLogCallback(userLogCallback);
  highs.readModel(filename);
  highs.run();

  highs.resetGlobalScheduler(true);
}

TEST_CASE("log-callback", "[highs_io]") {
  bool output_flag = true;
  bool log_to_console = false;
  HighsInt log_dev_level = kHighsLogDevLevelInfo;
  HighsLogOptions log_options;
  log_options.log_stream = stdout;
  log_options.output_flag = &output_flag;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  log_options.user_log_callback = myLogCallback;

  highsLogDev(log_options, HighsLogType::kInfo, "Hi %s!", "HiGHS");
  if (dev_run) printf("Log callback yields \"%s\"\n", alt_printed_log);
  REQUIRE(strcmp(alt_printed_log, "Hi HiGHS!") == 0);

  // Check that nothing is printed if the type is VERBOSE when
  // log_dev_level is kHighsLogDevLevelInfo;
  *alt_printed_log = '\0';
  highsLogDev(log_options, HighsLogType::kVerbose, "Hi %s!", "HiGHS");
  REQUIRE(*alt_printed_log == '\0');

  {
    char long_message[sizeof(alt_printed_log)];
    memset(long_message, 'H', sizeof(long_message));
    long_message[sizeof(long_message) - 2] = '\0';
    long_message[sizeof(long_message) - 1] = '\n';
    highsLogDev(log_options, HighsLogType::kInfo, long_message);
    if (dev_run) printf("Log callback yields \"%s\"\n", alt_printed_log);
    REQUIRE(strncmp(alt_printed_log, "HHHH", 4) == 0);
    REQUIRE(strlen(alt_printed_log) <= sizeof(alt_printed_log));
  }

  highsLogUser(log_options, HighsLogType::kInfo, "Hello %s!\n", "HiGHS");
  REQUIRE(strlen(alt_printed_log) > 9);
  REQUIRE(strcmp(alt_printed_log, "Hello HiGHS!\n") == 0);

  {
    char long_message[sizeof(alt_printed_log)];
    memset(long_message, 'H', sizeof(long_message));
    long_message[sizeof(long_message) - 2] = '\0';
    long_message[sizeof(long_message) - 1] = '\n';
    highsLogUser(log_options, HighsLogType::kWarning, long_message);
    if (dev_run) printf("Log callback yields \"%s\"\n", alt_printed_log);
    REQUIRE(strstr(alt_printed_log, "HHHH") != nullptr);
    REQUIRE(strlen(alt_printed_log) <= sizeof(alt_printed_log));
  }
}
