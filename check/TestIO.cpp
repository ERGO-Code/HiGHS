#include <cstdio>
#include <cstring>

#include "HighsIO.h"
#include "catch.hpp"

const bool dev_run = false;

char printedmsg[100000];
void* receiveddata = NULL;

// callback that saves message away for comparison
static void myprintmsgcb(int level, const char* msg, void* msgcb_data) {
  strcpy(printedmsg, msg);
  receiveddata = msgcb_data;
}

static void mylogmsgcb(HighsLogType type, const char* msg,
                       void* msgcb_data) {
  strcpy(printedmsg, msg);
  receiveddata = msgcb_data;
}

TEST_CASE("msgcb", "[highs_io]") {
  int dummydata = 42;
  bool output_flag = true;
  bool log_to_console = false;
  int log_dev_level = LOG_DEV_LEVEL_INFO;
  HighsLogOptions log_options;
  log_options.log_file_stream = stdout;
  log_options.output_flag = &output_flag;
  log_options.log_to_console = &log_to_console;
  log_options.log_dev_level = &log_dev_level;
  highsSetLogCallback(myprintmsgcb, mylogmsgcb, (void*)&dummydata);

  highsLogDev(log_options, HighsLogType::INFO, "Hi %s!", "HiGHS");
  REQUIRE(strcmp(printedmsg, "Hi HiGHS!") == 0);
  REQUIRE(receiveddata == &dummydata);

  // Check that nothing is printed if the type is VERBOSE when
  // log_dev_level is LOG_DEV_LEVEL_INFO;
  *printedmsg = '\0';
  highsLogDev(log_options, HighsLogType::VERBOSE, "Hi %s!", "HiGHS");
  REQUIRE(*printedmsg == '\0');

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg) - 2] = '\0';
    longmsg[sizeof(longmsg) - 1] = '\n';
    highsLogDev(log_options, HighsLogType::INFO, longmsg);
    REQUIRE(strncmp(printedmsg, "HHHH", 4) == 0);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }

  highsLogUser(log_options, HighsLogType::INFO, "Hello %s!\n", "HiGHS");
  REQUIRE(strlen(printedmsg) > 9);
  REQUIRE(strcmp(printedmsg, "         Hello HiGHS!\n") == 0);
  REQUIRE(receiveddata == &dummydata);

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg) - 2] = '\0';
    longmsg[sizeof(longmsg) - 1] = '\n';
    highsLogUser(log_options, HighsLogType::WARNING, longmsg);
    REQUIRE(strstr(printedmsg, "HHHH") != NULL);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }
}
