#include "Highs.h"
#include "HighsExternalDeps.h"
#include "catch.hpp"

// simple test for now
// extended tests require workflow changes to include/exclude dependencies
TEST_CASE("HighsExternalDeps", "[highs_external_deps]") {
  // tryLoad without a path — returns true only if extras are available
  bool loaded = HighsExternalDeps::tryLoad();

  std::string status = HighsExternalDeps::getLoadStatus();
  REQUIRE(!status.empty());

  if (loaded) {
    REQUIRE(HighsExternalDeps::isAvailable());
  } else {
    REQUIRE(!HighsExternalDeps::isAvailable());
  }
}

TEST_CASE("HighsExternalDeps-compile-time", "[highs_external_deps]") {
  // isAvailableAtCompile should be consistent with build config
  bool compile = HighsExternalDeps::isAvailableAtCompile();

  // If compile-time available, runtime must also be available
  if (compile) {
    REQUIRE(HighsExternalDeps::isAvailable());
  }
}

TEST_CASE("HighsExternalDeps-getCopyrightInfo", "[highs_external_deps]") {
  HighsExternalDeps::tryLoad();
  std::string info = HighsExternalDeps::getCopyrightInfo();

  if (HighsExternalDeps::isAvailable()) {
    REQUIRE(!info.empty());
  } else {
    REQUIRE(info.empty());
  }
}