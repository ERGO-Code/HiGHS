#include <cstdio>
#include <unistd.h>
#define GetCurrentDir getcwd

#include "FilereaderEms.h"
#include "HighsIO.h"
#include "HMpsFF.h"
#include "HMPSIO.h"
#include "HighsLp.h"
#include "HighsSetup.h"
#include "LoadProblem.h"
#include "catch.hpp"

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  GetCurrentDir(buff, FILENAME_MAX);
  std::string current_working_dir(buff);
  return current_working_dir;
}

TEST_CASE("free-format-parser", "[highs_filereader]") {
  std::string dir = GetCurrentWorkingDir();

  std::cout << dir << std::endl;

  // For debugging use the latter.
  // std::string filename = dir + "/../../check/instances/25fv47.mps";
  std::string filename = dir + "/check/instances/adlittle.mps";

  // Read mps.
  HighsLp lp_free_format, lp_fixed_format;

  HMpsFF parser{};
  int status = parser.loadProblem(filename, lp_free_format);

  std::vector<int> integerColumn;
  status = readMPS(filename.c_str(), -1, -1, lp_fixed_format.numRow_,
                   lp_fixed_format.numCol_, lp_fixed_format.sense_,
                   lp_fixed_format.offset_, lp_fixed_format.Astart_,
                   lp_fixed_format.Aindex_, lp_fixed_format.Avalue_,
                   lp_fixed_format.colCost_, lp_fixed_format.colLower_,
                   lp_fixed_format.colUpper_, lp_fixed_format.rowLower_,
                   lp_fixed_format.rowUpper_, integerColumn,
                   lp_fixed_format.row_names_, lp_fixed_format.col_names_);

  bool are_the_same = lp_free_format == lp_fixed_format;
  REQUIRE(are_the_same);
}

TEST_CASE("read-mps-ems", "[highs_filereader]") {
  HighsOptions options;
  std::string dir = GetCurrentWorkingDir();

  std::cout << dir << std::endl;

  // For debugging use the latter.
  options.filename = dir + "/../../check/instances/adlittle.mps";
  // options.filename = dir + "/check/instances/adlittle.mps";

  // Read mps.
  HighsLp lp_mps;
  HighsInputStatus read_status = loadLpFromFile(options, lp_mps);
  REQUIRE(read_status == HighsInputStatus::OK);

  // Write ems.
  FilereaderEms ems;
  ems.writeModelToFile("adlittle.ems", lp_mps);

  // Read ems and compare.
  options.filename = "adlittle.ems"; // todo: check how to specify path

  HighsLp lp_ems;
  HighsInputStatus ems_read_status = loadLpFromFile(options, lp_ems);
  REQUIRE(ems_read_status == HighsInputStatus::OK);

  bool are_the_same = lp_mps == lp_ems;
  REQUIRE(are_the_same);

  std::remove(options.filename.c_str());
}
