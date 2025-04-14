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
  const std::string compilation = highsCompilationDate();
  const std::string githash = std::string(highsGithash());
  std::stringstream ss;
  ss << major << "." << minor << "." << patch;
  std::string local_version = ss.str();
  if (dev_run) {
    printf("HiGHS version: %s\n", version.c_str());
    printf("HiGHS major version %d\n", int(major));
    printf("HiGHS minor version %d\n", int(minor));
    printf("HiGHS patch version %d\n", int(patch));
    printf("HiGHS githash: %s\n", githash.c_str());
    // Compilation date is deprecated, but make sure that the
    // deprecated method is still tested.
    printf("HiGHS compilation date: %s\n", compilation.c_str());
    printf("HiGHS local version: %s\n", local_version.c_str());
  }
  REQUIRE(major == HIGHS_VERSION_MAJOR);
  REQUIRE(minor == HIGHS_VERSION_MINOR);
  REQUIRE(patch == HIGHS_VERSION_PATCH);
  REQUIRE(githash == std::string(HIGHS_GITHASH));
  REQUIRE(version == local_version);
  // Check that the corresponding methods
  Highs highs;
  const std::string version0 = highs.version();
  REQUIRE(version0 == version);
  const HighsInt major0 = highs.versionMajor();
  REQUIRE(major0 == major);
  const HighsInt minor0 = highs.versionMinor();
  REQUIRE(minor0 == minor);
  const HighsInt patch0 = highs.versionPatch();
  REQUIRE(patch0 == patch);
  const std::string githash0 = highs.githash();
  REQUIRE(githash0 == githash);
  const std::string compilation0 = highs.compilationDate();
  REQUIRE(compilation == compilation);
}

TEST_CASE("sizeof-highs-int", "[highs_version]") {
  Highs highs;
  HighsInt sizeof_highs_int = highs.getSizeofHighsInt();
#ifdef HIGHSINT64
  REQUIRE(sizeof_highs_int == 8);
#else
  REQUIRE(sizeof_highs_int == 4);
#endif

  highs.resetGlobalScheduler(true);
}
