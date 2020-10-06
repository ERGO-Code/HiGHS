#include <fstream>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const std::string basis_file = "Highs.bas";
HighsBasis basis_data;

void testBasisRestart(Highs& highs, const bool from_file) {
  HighsStatus return_status;
  //  highs.writeSolution("", true);
  // Change a bound and resolve

  const int changeCol = 3;
  const double old_lower_bound = highs.getLp().colLower_[changeCol];
  const double old_upper_bound = highs.getLp().colUpper_[changeCol];
  const double new_lower_bound = 2;
  highs.changeColBounds(changeCol, new_lower_bound, old_upper_bound);

  return_status = highs.run();

  if (dev_run) {
    printf(
	   "After modifying lower bound of column %d from %g to %g, solving the LP "
	   "requires %d iterations and objective is %g\n",
	   changeCol, old_lower_bound, new_lower_bound,
	   highs.getSimplexIterationCount(), highs.getObjectiveValue());
    //  highs.writeSolution("", true);
  }
  // Make sure that the test requires iterations
  assert(highs.getSimplexIterationCount()>0);

  // Recover bound, load optimal basis and resolve

  highs.changeColBounds(changeCol, old_lower_bound, old_upper_bound);

  if (from_file) {
    return_status = highs.readBasis(basis_file);
  } else {
    return_status = highs.setBasis(basis_data);
  }
  REQUIRE(return_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  if (dev_run) {
    printf(
	   "After restoring lower bound of column %d from %g to %g, solving the LP "
	   "requires %d iterations and objective is %g\n",
	   changeCol, new_lower_bound, old_lower_bound,
	   highs.getSimplexIterationCount(), highs.getObjectiveValue());
  }

  REQUIRE(highs.getSimplexIterationCount() == 0);
}

// No commas in test case name.
TEST_CASE("Basis-file", "[highs_basis_file]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string model = "";
  std::string model_file;

  HighsOptions options;
  HighsLp lp;
  //  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;

  Highs highs(options);

  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  // Read mps
  model_file = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

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

  // Write and read a file for incomplete basis
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS Version 1" << std::endl;
  f << highs.getNumCols() << " " << highs.getNumRows() << std::endl;
  f << "1 1" << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::Error);

  testBasisRestart(highs, true);

  highs.clearModel();

  return_status = highs.readBasis(basis_file);
  REQUIRE(return_status == HighsStatus::Error);
  
  std::remove(basis_file.c_str());
}

// No commas in test case name.
TEST_CASE("Basis-save", "[highs_basis_save]") {
  std::cout << std::string(HIGHS_DIR) << std::endl;

  std::string model = "";
  std::string model_file;

  HighsOptions options;
  HighsLp lp;
  //  HighsStatus run_status;
  HighsStatus return_status;
  HighsStatus read_status;

  Highs highs(options);
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }

  // Read mps
  model_file = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";

  return_status = highs.setHighsOptionValue("model_file", model_file);
  REQUIRE(return_status == HighsStatus::OK);

  read_status = highs.readModel(model_file);
  REQUIRE(read_status == HighsStatus::OK);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::OK);

  basis_data = highs.getBasis();
  REQUIRE(return_status == HighsStatus::OK);

  testBasisRestart(highs, false);

  highs.clearModel();

}
