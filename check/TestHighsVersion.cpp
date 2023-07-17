#include <cassert>
#include <sstream>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("HighsVersion", "[highs_version]") {
  std::string version = std::string(highsVersion());
  const HighsInt major = highsVersionMajor();
  const HighsInt minor = highsVersionMinor();
  const HighsInt patch = highsVersionPatch();
  std::stringstream ss;
  ss << major << "." << minor << "." << patch;
  std::string local_version = ss.str();
  if (dev_run) {
    printf("HiGHS version: %s\n", version.c_str());
    printf("HiGHS major version %d\n", int(major));
    printf("HiGHS minor version %d\n", int(minor));
    printf("HiGHS patch version %d\n", int(patch));
    printf("HiGHS githash: %s\n", highsGithash());
    printf("HiGHS compilation date: %s\n", highsCompilationDate());
    printf("HiGHS local version: %s\n", local_version.c_str());
  }
  REQUIRE(major == HIGHS_VERSION_MAJOR);
  REQUIRE(minor == HIGHS_VERSION_MINOR);
  REQUIRE(patch == HIGHS_VERSION_PATCH);
  REQUIRE(local_version == version);
}
