#include <cstdio>
#include <unistd.h>
#define GetCurrentDir getcwd

#include "FilereaderEms.h"
#include "HighsIO.h"
#include "HighsLp.h"
#include "Highs.h"
#include "LoadProblem.h"
#include "catch.hpp"

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

// No commas in test case name.
TEST_CASE("read-mps-ems", "[highs_filereader]") {
  HighsOptions options;
  std::string dir = GetCurrentWorkingDir();

  std::cout << dir << std::endl;

  // For debugging use the latter.
  options.filename = dir + "/../../check/instances/adlittle.mps";
  // options.filename = dir + "/check/instances/adlittle.mps";

  // Read mps.
  HighsLp lp_mps;
  HighsStatus read_status = loadLpFromFile(options, lp_mps);
  REQUIRE(read_status == HighsStatus::OK);

  // Write ems.
  FilereaderEms ems;
  ems.writeModelToFile("adlittle.ems", lp_mps);

  // Read ems and compare.
  options.filename = "adlittle.ems"; // todo: check how to specify path

  HighsLp lp_ems;
  HighsStatus ems_read_status = loadLpFromFile(options, lp_ems);
  REQUIRE(ems_read_status == HighsStatus::OK);

  bool are_the_same = lp_mps == lp_ems;
  REQUIRE(are_the_same);

  std::remove(options.filename.c_str());
}
