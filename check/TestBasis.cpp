#include "Highs.h"
#include "catch.hpp"

// No commas in test case name.
TEST_CASE("Basis", "[highs_basis]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string model = "";
  std::string model_file;

  HighsOptions options;
  HighsLp lp;
  //  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;

  Highs highs(options);

  // Read mps
  model_file = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  std::string basis_file = "Highs.bas";
  return_status = highs.writeBasis(basis_file);
  REQUIRE(return_status == HighsStatus::OK);

  /*
  return_status = highs.readBasis("Nobasis.bas");
  REQUIRE(return_status == HighsStatus::Error);
  */

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  REQUIRE(highs.getSimplexIterationCount() == 0);
}
