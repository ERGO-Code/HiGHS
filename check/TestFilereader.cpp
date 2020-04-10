#include <cstdio>

#include "Highs.h"
#include "catch.hpp"
#include "io/FilereaderEms.h"
#include "io/HMPSIO.h"
#include "io/HMpsFF.h"
#include "io/HighsIO.h"
#include "lp_data/HighsLp.h"
#include "lp_data/HighsLpUtils.h"

TEST_CASE("free-format-parser", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  HighsStatus status;

  // Read mps
  HighsOptions options;

  Highs highs(options);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_free = highs.getLp();

  status = highs.setHighsOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::OK);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_fixed = highs.getLp();

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}

// No commas in test case name.
TEST_CASE("read-mps-ems", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  HighsStatus status;

  // Read mps
  HighsOptions options;

  Highs highs(options);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);
  HighsLp lp_mps = highs.getLp();

  // Write ems.
  filename = "adlittle.ems";
  status = highs.writeModel(filename);
  REQUIRE(status == HighsStatus::OK);

  // Read ems and compare.
  options.model_file = "adlittle.ems";  // todo: check how to specify path

  status = highs.setHighsOptionValue("model_file", filename);
  REQUIRE(status == HighsStatus::OK);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_ems = highs.getLp();

  bool are_the_same = lp_mps == lp_ems;
  REQUIRE(are_the_same);

  std::remove("adlittle.ems");
}

TEST_CASE("integrality-constraints", "[highs_filereader]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/small_mip.mps";

  // integer variables are COL03,COL04 so x[2], x[3].
  const std::vector<int> kIntegers{0, 0, 1, 1, 0, 0, 0, 0};

  HighsStatus status;
  HighsOptions options;

  Highs highs(options);
  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_free = highs.getLp();

  REQUIRE(lp_free.integrality_.size() == lp_free.numCol_);
  REQUIRE(lp_free.integrality_ == kIntegers);

  // Read mps with fixed format parser.
  status = highs.setHighsOptionValue("mps_parser_type_free", false);
  REQUIRE(status == HighsStatus::OK);

  status = highs.readModel(filename);
  REQUIRE(status == HighsStatus::OK);

  HighsLp lp_fixed = highs.getLp();

  REQUIRE(lp_fixed.integrality_.size() == lp_fixed.numCol_);
  REQUIRE(lp_fixed.integrality_ == kIntegers);

  bool are_the_same = lp_free == lp_fixed;
  REQUIRE(are_the_same);
}

TEST_CASE("dualize", "[highs_data]") {
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  // Read mps.
  HighsOptions options;
  options.model_file = filename;

  HighsLp lp;
  HMpsFF parser{};
  FreeFormatParserReturnCode result = parser.loadProblem(stdout, filename, lp);
  REQUIRE(result == FreeFormatParserReturnCode::SUCCESS);

  HighsLp primal;
  HighsStatus status;
  status = transformIntoEqualityProblem(lp, primal);
  REQUIRE(status == HighsStatus::OK);

  Highs highs_lp;
  HighsModelStatus model_status;
  status = highs_lp.passModel(lp);
  REQUIRE(status == HighsStatus::OK);
  status = highs_lp.run();
  model_status = highs_lp.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  Highs highs_primal;
  status = highs_primal.passModel(primal);
  REQUIRE(status == HighsStatus::OK);
  status = highs_primal.run();
  model_status = highs_lp.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  double lp_objective;
  highs_lp.getHighsInfoValue("objective_function_value", lp_objective);
  double primal_objective;
  highs_lp.getHighsInfoValue("objective_function_value", primal_objective);

  double diff_equality = lp_objective - primal_objective;
  REQUIRE(diff_equality < 0.00000001);

  HighsLp dual;
  status = dualizeEqualityProblem(primal, dual);
  REQUIRE(status == HighsStatus::OK);
  Highs highs_dual;
  status = assessLp(dual, options);
  REQUIRE(status == HighsStatus::OK);
  status = highs_dual.passModel(dual);
  REQUIRE(status == HighsStatus::OK);
  status = highs_dual.run();
  model_status = highs_lp.getModelStatus();
  REQUIRE(model_status == HighsModelStatus::OPTIMAL);

  double dual_objective;
  highs_dual.getHighsInfoValue("objective_function_value", dual_objective);

  double diff_dual = primal_objective + dual_objective;
  REQUIRE(diff_dual < 0.00000001);
}
