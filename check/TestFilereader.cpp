#include "catch.hpp"
#include "HighsLp.h"
#include "LoadProblem.h"
#include "HighsSetup.h"

// No commas in test case name.
TEST_CASE("read-mps-ems", "[highs_filereader]") {
  HighsOptions options;
  options.filename_ = "/check/instances/adlittle.mps"; // todo: check how to specify path
  
  // Read mps.
  HighsLp lp_mps;
  HighsInputStatus read_status = loadLpFromFile(options, lp_mps);
  REQUIRE(read_status == HighsInputStatus::OK);

  // Write ems.
  FilereaderEms ems;
  ems.writeModelToFile("filereader-test-ems-file.ems", lp_mps);

  // Read ems and compare.
  options.filename_ = "filereader-test-ems-file.ems" // todo: check how to specify path
  HighsLp lp_ems;
  HighsInputStatus ems_read_status = loadLpFromFile(options, lp_ems);

  REQUIRE(lp_mps == lp_ems); 
}
