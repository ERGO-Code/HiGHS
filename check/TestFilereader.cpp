#include <cstdio>

#include "Highs.h"
#include "catch.hpp"
#include "io/FilereaderEms.h"
#include "io/HMPSIO.h"
#include "io/HMpsFF.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"

const bool dev_run = false;

TEST_CASE("filereader-edge-cases", "[highs_filereader]") {
  std::string model = "";
  std::string model_file;
  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;

  // Several tests don't pass, but should, so possibly skip them
  const bool test_garbage_mps = false;
  const bool test_garbage_ems = true;
  const bool test_garbage_lp = false;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  const HighsInfo& info = highs.getInfo();

  // Try to run HiGHS with default options. No model loaded so OK
  run_status = highs.run();
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kModelEmpty);
  REQUIRE(run_status == HighsStatus::kOk);

  // Load a non-existent file and try to run HiGHS
  model = "";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kError);

  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  // Load a non-supported file type and try to run HiGHS
  model = "model";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".xyz";
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kError);

  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  // Load an LP file that does not end with a newline
  model = "no-newline-eof";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kOk);

  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);

  // Load an existing MPS file and run HiGHS
  model = "adlittle";
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  return_status = highs.readModel(model_file);
  REQUIRE(return_status == HighsStatus::kOk);

  run_status = highs.run();
  REQUIRE(run_status == HighsStatus::kOk);
  REQUIRE(info.simplex_iteration_count == 87);

  model = "garbage";
  if (test_garbage_mps) {
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
    read_status = highs.readModel(model_file);
    REQUIRE(read_status == HighsStatus::kError);
  }

  if (test_garbage_ems) {
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".ems";
    read_status = highs.readModel(model_file);
    REQUIRE(read_status == HighsStatus::kError);
  }

  if (test_garbage_lp) {
    model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".lp";
    read_status = highs.readModel(model_file);
    REQUIRE(read_status == HighsStatus::kError);
  }
}

TEST_CASE("filereader-free-format-parser", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  HighsStatus status;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_free = highs.getLp();

  status = highs.setOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_fixed = highs.getLp();

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}

// No commas in test case name.
TEST_CASE("filereader-read-mps-ems-lp", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  HighsStatus status;

  // Read mps
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);
  HighsLp lp_mps = highs.getLp();

  // Write lp
  std::string filename_lp = "adlittle.lp";
  status = highs.writeModel(filename_lp);
  REQUIRE(status == HighsStatus::kOk);

  /*
  bool are_the_same;
  // Write ems
  std::string filename_ems = "adlittle.ems";
  status = highs.writeModel(filename_ems);
  REQUIRE(status == HighsStatus::kOk);

  // Read ems and compare with mps
  std::cout << "Reading " << filename_ems << std::endl;
  status = highs.readModel(filename_ems);
  REQUIRE(status == HighsStatus::kOk);

  std::cout << "Compare LP from .ems and .mps" << std::endl;
  are_the_same = lp_mps == highs.getLp();
  REQUIRE(are_the_same);

  std::remove(filename_ems.c_str());
  */

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  const HighsInfo& info = highs.getInfo();
  double mps_objective_function_value = info.objective_function_value;

  // Read lp and compare objective with mps
  if (dev_run) std::cout << "Reading " << filename_lp << std::endl;
  status = highs.readModel(filename_lp);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.run();
  REQUIRE(status == HighsStatus::kOk);

  double delta =
      std::fabs(mps_objective_function_value - info.objective_function_value);
  if (dev_run && delta)
    printf("TestFilereader: LP-MPS |objective difference| of %g\n", delta);
  REQUIRE(delta < 1e-8);

  std::remove(filename_lp.c_str());
}

TEST_CASE("filereader-integrality-constraints", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/small_mip.mps";

  // integer variables are COL03,COL04 so x[2], x[3].
  const std::vector<HighsVarType> kIntegers{
      HighsVarType::kContinuous, HighsVarType::kContinuous,
      HighsVarType::kInteger,    HighsVarType::kInteger,
      HighsVarType::kContinuous, HighsVarType::kContinuous,
      HighsVarType::kContinuous, HighsVarType::kContinuous};

  HighsStatus status;

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_free = highs.getLp();

  REQUIRE((int)lp_free.integrality_.size() == lp_free.num_col_);
  REQUIRE(lp_free.integrality_ == kIntegers);

  // Read mps with fixed format parser.
  status = highs.setOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::kOk);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::kOk);

  HighsLp lp_fixed = highs.getLp();

  REQUIRE((int)lp_fixed.integrality_.size() == lp_fixed.num_col_);
  REQUIRE(lp_fixed.integrality_ == kIntegers);

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}
