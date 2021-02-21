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

static void mylogmsgcb(HighsMessageType type, const char* msg,
                       void* msgcb_data) {
  strcpy(printedmsg, msg);
  receiveddata = msgcb_data;
}

TEST_CASE("msgcb", "[highs_io]") {
  int dummydata = 42;
  bool output_flag = true;
  bool log_to_console = false;
  int output_dev = OUTPUT_DEV_INFO;
  HighsIoOptions io_options;
  io_options.logging_file = stdout;
  io_options.output_flag = &output_flag;
  io_options.log_to_console = &log_to_console;
  io_options.output_dev = &output_dev;

  highsSetMessageCallback(myprintmsgcb, mylogmsgcb, (void*)&dummydata);

  highsOutputDev(io_options, HighsMessageType::INFO, "Hi %s!", "HiGHS");
  REQUIRE(strcmp(printedmsg, "Hi HiGHS!") == 0);
  REQUIRE(receiveddata == &dummydata);

  // Check that nothing is printed if the type is VERBOSE when
  // output_dev is OUTPUT_DEV_INFO;
  *printedmsg = '\0';
  highsOutputDev(io_options, HighsMessageType::VERBOSE, "Hi %s!", "HiGHS");
  REQUIRE(*printedmsg == '\0');

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg) - 2] = '\0';
    longmsg[sizeof(longmsg) - 1] = '\n';
    highsOutputDev(io_options, HighsMessageType::INFO, longmsg);
    REQUIRE(strncmp(printedmsg, "HHHH", 4) == 0);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }

  highsOutputUser(io_options, HighsMessageType::INFO, "Hello %s!\n", "HiGHS");
  REQUIRE(strlen(printedmsg) > 9);
  REQUIRE(strcmp(printedmsg, "         Hello HiGHS!\n") == 0);
  REQUIRE(receiveddata == &dummydata);

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg) - 2] = '\0';
    longmsg[sizeof(longmsg) - 1] = '\n';
    highsOutputUser(io_options, HighsMessageType::WARNING, longmsg);
    REQUIRE(strstr(printedmsg, "HHHH") != NULL);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }
}
