#include "Avgas.h"
#include "Highs.h"
#include "catch.hpp"
#include "io/LoadProblem.h"
//#include "lp_data/HighsLpUtils.h"

// No commas in test case name.
TEST_CASE("LP-simplex", "[highs_simplex]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  HighsOptions options;
  options.model_file = filename;

  // Read mps.
  HighsLp lp;
  HighsStatus read_status = loadLpFromFile(options, lp);
  REQUIRE(read_status == HighsStatus::OK);

  Highs highs(options);
  HighsStatus return_status = highs.passModel(lp);
  REQUIRE(return_status == HighsStatus::OK);

  // Vanilla solve
  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);
  
  int simplex_iteration_count;
  return_status = highs.getHighsInfoValue("simplex_iteration_count", simplex_iteration_count);
  REQUIRE(return_status == HighsStatus::OK);
  REQUIRE(simplex_iteration_count == 86);
  

  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  options.model_file = filename;
  read_status = loadLpFromFile(options, lp);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.passModel(lp);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.setHighsOptionValue("simplex_iteration_limit", 10);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  /*
  return_status = highs.getHighsInfoValue("simplex_iteration_count", simplex_iteration_count);
  REQUIRE(return_status == HighsStatus::OK);
  REQUIRE(simplex_iteration_count == 86);
  */
}
