#include <cstdio>
#include <cstring>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

const HighsInt kLogBufferSize = kIoBufferSize;

char printed_log[kLogBufferSize];
void* received_data = NULL;

using std::memset;
using std::strcmp;
using std::strcpy;
using std::strlen;
using std::strncmp;
using std::strstr;

// Callback that saves message for comparison
static void myLogCallback(HighsLogType type, const char* message,
                          void* log_callback_data) {
  strcpy(printed_log, message);
  received_data = log_callback_data;
}

// Callback that provides user logging
static void userLogCallback(HighsLogType type, const char* message,
                            void* log_callback_data) {
  if (dev_run) printf("userLogCallback: %s", message);
}

TEST_CASE("run-callback", "[highs_io]") {
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.setLogCallback(userLogCallback);
  highs.readModel(filename);
  highs.run();
}

TEST_CASE("log-callback", "[highs_io]") {
  int dummy_data = 42;
  bool output_flag = true;
  bool log_to_console = false;
  HighsInt log_dev_level = kHighsLogDevLevelInfo;
  HighsLogOptions log_options;
  log_options.log_file_stream = stdout;
  log_options.output_flag = &output_flag;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  log_options.log_callback = myLogCallback;
  log_options.log_callback_data = (void*)&dummy_data;

  highsLogDev(log_options, HighsLogType::kInfo, "Hi %s!", "HiGHS");
  if (dev_run) printf("Log callback yields \"%s\"\n", printed_log);
  REQUIRE(strcmp(printed_log, "Hi HiGHS!") == 0);
  REQUIRE(received_data == &dummy_data);

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
  REQUIRE(received_data == &dummy_data);

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
