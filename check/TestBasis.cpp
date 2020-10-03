#include <fstream>

#include "Highs.h"
#include "catch.hpp"

// No commas in test case name.
TEST_CASE("Basis", "[highs_basis_io]") {
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
  model_file = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  std::string basis_file = "Highs.bas";
  return_status = highs.writeBasis(basis_file);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.readBasis("Nobasis.bas");
  REQUIRE(return_status == HighsStatus::Error);

  // Check error return for some invalid basis files
  std::string invalid_basis_file = "InvalidBasis.bas";
  std::ofstream f;
  // Write and read a file for unsupported HiGHS version
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS Version 0" << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::Error);

  // Write and read a file for incompatible number of columns
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS Version 1" << std::endl;
  f << highs.getNumCols() - 1 << " " << highs.getNumRows() << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::Error);

  // Write and read a file for incompatible number of rows
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS Version 1" << std::endl;
  f << highs.getNumCols() << " " << 0 << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::Error);

  // Write and read a file for incoplete basis
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS Version 1" << std::endl;
  f << highs.getNumCols() << " " << highs.getNumRows() << std::endl;
  f << "1 1" << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::Error);

  highs.writeSolution("", true);
  // Change some bounds and resolve

  const int changeCol = 3;
  const double old_lower_bound = highs.getLp().colLower_[changeCol];
  const double old_upper_bound = highs.getLp().colUpper_[changeCol];
  const double new_lower_bound = 1;
  highs.changeColBounds(changeCol, new_lower_bound, old_upper_bound);

  return_status = highs.run();

  printf(
      "After modifying lower bound of column %d from %g to %g, solving the LP "
      "requires %d iterations and objective is %g\n",
      changeCol, old_lower_bound, new_lower_bound,
      highs.getSimplexIterationCount(), highs.getObjectiveValue());
  highs.writeSolution("", true);
  /*
  return_status = highs.readBasis(basis_file);
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  REQUIRE(highs.getSimplexIterationCount() == 0);
  */
  //  std::remove(basis_file.c_str());
}
