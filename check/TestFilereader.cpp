#include <cstdio>
#include <unistd.h>
#define GetCurrentDir getcwd

#include "FilereaderEms.h"
#include "HMPSIO.h"
#include "HMpsFF.h"
#include "Highs.h"
#include "HighsIO.h"
#include "HighsLp.h"
#include "LoadProblem.h"
#include "catch.hpp"

std::string GetCurrentWorkingDir(void) {
  char buff[FILENAME_MAX];
  auto result = GetCurrentDir(buff, FILENAME_MAX);
  if (result) {
    std::string current_working_dir(buff);
    return current_working_dir;
  }
  return "";
}

TEST_CASE("free-format-parser", "[highs_filereader]") {
  std::string dir = GetCurrentWorkingDir();

  std::cout << dir << std::endl;

  // For debugging use the latter.
  std::string filename = dir + "/../../check/instances/adlittle.mps";
  // std::string filename = dir + "/check/instances/adlittle.mps";

  // Read mps.
  HighsLp lp_free_format, lp_fixed_format;
  bool are_the_same = false;

  std::vector<int> integerColumn;
  int status = readMPS(filename.c_str(), -1, -1, lp_fixed_format.numRow_,
                   lp_fixed_format.numCol_, lp_fixed_format.sense_,
                   lp_fixed_format.offset_, lp_fixed_format.Astart_,
                   lp_fixed_format.Aindex_, lp_fixed_format.Avalue_,
                   lp_fixed_format.colCost_, lp_fixed_format.colLower_,
                   lp_fixed_format.colUpper_, lp_fixed_format.rowLower_,
                   lp_fixed_format.rowUpper_, integerColumn,
                   lp_fixed_format.row_names_, lp_fixed_format.col_names_);
  lp_fixed_format.nnz_ = lp_fixed_format.Avalue_.size();
  if (!status) {
    HMpsFF parser{};
    FreeFormatParserReturnCode result = parser.loadProblem(filename, lp_free_format);
    if (result != FreeFormatParserReturnCode::SUCCESS)
      status = 1;
    if (!status)
      are_the_same = lp_free_format == lp_fixed_format;
  }

  // In case you want to compare.
  // FilereaderEms ems;
  // ems.writeModelToFile("fixed.ems", lp_fixed_format);
  // ems.writeModelToFile("free.ems", lp_free_format);

  REQUIRE(are_the_same);
}

// No commas in test case name.
TEST_CASE("load-options-from-file", "[highs_data]") {
  HighsOptions options;
  std::string dir = GetCurrentWorkingDir();
  
  // For debugging use the latter.
  options.options_file= dir + "/../../check/sample_options_file";
  // options.options_file = dir + "/check/sample_options_file";

  bool success = loadOptionsFromFile(options); 

  REQUIRE(success == true);
  REQUIRE(options.small_matrix_value == 0.001);
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

TEST_CASE("integrality-constraints", "[highs_filereader]") {
  std::string dir = GetCurrentWorkingDir();

  // For debugging use the latter.
  std::string filename = dir + "/../../check/instances/small_mip.mps";
  // std::string filename = dir + "/check/instances/small_mip.mps";

  HighsOptions options;
  options.filename = filename;
  // integer variables are COL03,COL04 so x[2], x[3].
  const std::vector<int> kIntegers{0, 0, 1, 1, 0, 0, 0, 0};

  // Read mps with fixed format parser.
  HighsLp lp_fixed;
  options.parser_type = HighsMpsParserType::fixed;

  HighsStatus read_status = loadLpFromFile(options, lp_fixed);
  REQUIRE(read_status == HighsStatus::OK);
  REQUIRE(lp_fixed.integrality_.size() == lp_fixed.numCol_);
  REQUIRE(lp_fixed.integrality_ == kIntegers);

  // Read mps with free format parser.
  HighsLp lp_free;
  options.parser_type = HighsMpsParserType::free;

  read_status = loadLpFromFile(options, lp_free);
  REQUIRE(read_status == HighsStatus::OK);
  REQUIRE(lp_free.integrality_.size() == lp_free.numCol_);
  REQUIRE(lp_free.integrality_ == kIntegers);
}

TEST_CASE("dualize", "[highs_data]") {

  std::string dir = GetCurrentWorkingDir();

  // For debugging use the latter.
  std::string filename = dir + "/../../check/instances/adlittle.mps";
  // std::string filename = dir + "/check/instances/adlittle.mps";

  // Read mps.
  HighsOptions options;
  options.filename = filename;

  HighsLp lp;
  HMpsFF parser{};
  FreeFormatParserReturnCode result = parser.loadProblem(filename, lp);
  REQUIRE(result == FreeFormatParserReturnCode::SUCCESS);

  HighsLp primal = transformIntoEqualityProblem(lp);
  HighsLp dual = dualizeEqualityProblem(primal);

  Highs highs_lp;
  highs_lp.initializeLp(lp);
  highs_lp.run();

  Highs highs_primal;
  highs_primal.initializeLp(primal);
  highs_primal.run();

  Highs highs_dual;
  highs_dual.initializeLp(dual);
  highs_dual.run();

  double lp_objective = highs_lp.getObjectiveValue();
  double primal_objective = highs_primal.getObjectiveValue();
  double dual_objective = highs_dual.getObjectiveValue();

  double diff_equality = lp_objective - primal_objective;
  REQUIRE(diff_equality < 0.00000001);

  double diff_dual = primal_objective - dual_objective;
  REQUIRE(diff_dual < 0.00000001);
}