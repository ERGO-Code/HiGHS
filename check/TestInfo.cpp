#include <cstdio>

#include "FilereaderEms.h"
#include "HMPSIO.h"
#include "Highs.h"
#include "catch.hpp"

TEST_CASE("highs-info", "[highs_info]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  //  filename = std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";
  
  Highs highs;
  HighsStatus return_status = highs.readModel(filename);
  REQUIRE(return_status==HighsStatus::OK);
 
  return_status = highs.writeHighsInfo("");
  REQUIRE(return_status == HighsStatus::Error);

  return_status = highs.writeHighsInfo("Highs.info");
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  double objective_function_value;
  return_status = highs.getHighsInfoValue("objective_value", objective_function_value);
  REQUIRE(return_status == HighsStatus::Error);
  
  return_status = highs.getHighsInfoValue("objective_function_value", objective_function_value);
  REQUIRE(return_status == HighsStatus::OK);

  printf("From getHighsInfoValue: objective_function_value = %g\n", objective_function_value);
  
  int simplex_iteration_count;
  return_status = highs.getHighsInfoValue("iteration_count", simplex_iteration_count);
  REQUIRE(return_status == HighsStatus::Error);
  
  return_status = highs.getHighsInfoValue("simplex_iteration_count", simplex_iteration_count);
  REQUIRE(return_status == HighsStatus::OK);

  printf("From getHighsInfoValue: simplex_iteration_count = %d\n", simplex_iteration_count);
  
}    



