#include <cstdio>

#include "HighsIO.h"
#include "catch.hpp"

char printedmsg[100000];

// callback that saves message away for comparison
static
void myprintmsgcb(unsigned int level, const char* msg, void* msgcb_data) {
  strcpy(printedmsg, msg);
}

static
void mylogmsgcb(HighsMessageType type, const char* msg, void* msgcb_data) {
  strcpy(printedmsg, msg);
}

TEST_CASE("msgcb", "[highs_io]") {

  HighsSetMessageCallback(myprintmsgcb, mylogmsgcb, NULL);
  
  HighsPrintMessage(4, "Hi %s!", "HiGHS");
  REQUIRE(strcmp(printedmsg, "Hi HiGHS!") == 0);

  /* printed at level 4 when level is 3 should not print */
  *printedmsg = '\0';
  HighsSetMessagelevel(3);
  HighsPrintMessage(4, "Hi %s!", "HiGHS");
  REQUIRE(*printedmsg == '\0');

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg)-1] = '\0';
    HighsPrintMessage(2, longmsg);
    REQUIRE(strncmp(printedmsg, "HHHH", 4) == 0);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }

  HighsLogMessage(HighsMessageType::INFO, "Hello %s!", "HiGHS");
  REQUIRE(strlen(printedmsg) > 8);
  REQUIRE(strcmp(printedmsg+8, " [INFO] Hello HiGHS!\n") == 0);  // begin of printedmsg is a timestamp, which we skip over

  {
    char longmsg[sizeof(printedmsg)];
    memset(longmsg, 'H', sizeof(longmsg));
    longmsg[sizeof(longmsg)-1] = '\0';
    HighsLogMessage(HighsMessageType::WARNING, longmsg);
    REQUIRE(strstr(printedmsg, "HHHH") != NULL);
    REQUIRE(strlen(printedmsg) <= sizeof(printedmsg));
  }


  // TODO test msgcb_data is correctly passed on
}
