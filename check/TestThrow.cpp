//#include "Highs.h"
#include <iostream>
#include <stdexcept>  //Used in HiGHS

#include "catch.hpp"
//#include <exception>
const bool dev_run = false;

void invalidArgument() {
  // Used in .lp file reader
  throw std::invalid_argument("Exception: invalid_argument");
}

void logicError() {
  // Used in IPX
  throw std::logic_error("Exception: logic_error");
}

void badAlloc() {
  // Used in IPX
  throw std::bad_alloc();
}

TEST_CASE("ThrowCatch", "[highs_test_throw]") {
  bool ok;
  // Check that all exceptions usind in HiGHS are caught
  //
  // Catch invalid_argument explicitly
  ok = false;
  try {
    invalidArgument();
  } catch (const std::invalid_argument& exception) {
    if (dev_run) std::cout << exception.what() << std::endl;
    ok = true;
  }
  REQUIRE(ok);
  // Catch logic_error explicitly
  ok = false;
  try {
    logicError();
  } catch (const std::logic_error& exception) {
    if (dev_run) std::cout << exception.what() << std::endl;
    ok = true;
  }
  REQUIRE(ok);
  // Catch bad_alloc explicitly
  ok = false;
  try {
    badAlloc();
  } catch (const std::bad_alloc& exception) {
    if (dev_run) std::cout << exception.what() << std::endl;
    ok = true;
  }
  REQUIRE(ok);
  // Catch invalid_argument via exception
  ok = false;
  try {
    invalidArgument();
  } catch (const std::exception& exception) {
    if (dev_run) std::cout << exception.what() << std::endl;
    ok = true;
  }
  REQUIRE(ok);
  // Catch logic_error via exception
  ok = false;
  try {
    logicError();
  } catch (const std::exception& exception) {
    if (dev_run) std::cout << exception.what() << std::endl;
    ok = true;
  }
  REQUIRE(ok);
  // Catch bad_alloc via exception
  ok = false;
  try {
    badAlloc();
  } catch (const std::exception& exception) {
    if (dev_run) std::cout << exception.what() << std::endl;
    ok = true;
  }
  REQUIRE(ok);
}
