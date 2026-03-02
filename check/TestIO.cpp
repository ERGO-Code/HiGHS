#include <cstdio>
#include <cstring>
#include <fstream>

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
  highs.setOptionValue("output_flag", dev_run);
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

HighsCallbackFunctionType userLoggingCallback =
    [](int callback_type, const std::string& message,
       const HighsCallbackOutput* data_out, HighsCallbackInput* data_in,
       void* user_callback_data) {
      fprintf(static_cast<FILE*>(user_callback_data), "%s", message.c_str());
    };

TEST_CASE("console-file-callback-log", "[highs_io]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  const std::string test_highs_log = test_name + ".log";
  const std::string test_user_log = test_name + "user.log";
  std::string filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.readModel(filename);
  FILE* file;
  // Run through all 8 possibilities of console/file/callback logging
  // on/off
  //
  // Console logging appearing appropriately can only be tested by eye
  // in debug, the others are tested by the log/callback file existing
  h.setOptionValue("output_flag", true);
  for (HighsInt k = 0; k < 8; k++) {
    const bool callback_log = (k & 1) != 0;
    const bool file_log = (k & 2) != 0;
    const bool console_log = (k & 4) != 0 && dev_run;
    if (dev_run)
      printf("\nCase k = %d: %7s %4s %8s\n\n", int(k),
             console_log ? "console" : "       ", file_log ? "file" : "    ",
             callback_log ? "callback" : "        ");
    h.setOptionValue("log_to_console", console_log);
    if (file_log) h.openLogFile(test_highs_log);
    if (callback_log) {
      file = fopen(test_user_log.c_str(), "w");
      //  void* p_user_callback_data = file;
      h.setCallback(userLoggingCallback, file);  // p_user_callback_data);
      h.startCallback(kCallbackLogging);
    }
    h.run();
    if (file_log) {
      h.closeLogFile();
      std::ifstream fin(test_highs_log.c_str());
      const bool exists = fin.is_open();
      REQUIRE(exists);
      fin.close();
      std::remove(test_highs_log.c_str());
      h.openLogFile("");
    } else {
      std::ifstream fin(test_highs_log.c_str());
      const bool exists = fin.is_open();
      REQUIRE(!exists);
    }
    if (callback_log) {
      fclose(file);
      h.stopCallback(kCallbackLogging);
      std::ifstream fin(test_user_log.c_str());
      const bool exists = fin.is_open();
      REQUIRE(exists);
      fin.close();
      std::remove(test_user_log.c_str());
    } else {
      std::ifstream fin(test_user_log.c_str());
      const bool exists = fin.is_open();
      REQUIRE(!exists);
    }
  }

  //  std::remove(test_highs_log.c_str());
  //  std::remove(test_user_log.c_str());
  h.resetGlobalScheduler(true);
}
