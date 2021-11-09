#include <fstream>

#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
const std::string basis_file = "adlittle.bas";
HighsBasis basis_data;

void testBasisReloadModel(Highs& highs, const bool from_file);
void testBasisRestart(Highs& highs, const bool from_file);

// No commas in test case name.
TEST_CASE("Basis-file", "[highs_basis_file]") {
  HighsStatus return_status;
  std::string model0_file =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  std::string model1_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  Highs highs;
  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
  }
  assert(model0_file != model1_file);

  return_status = highs.readModel(model0_file);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.writeBasis(basis_file);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.readBasis("Nobasis.bas");
  REQUIRE(return_status == HighsStatus::kError);

  // Check error return for some invalid basis files
  std::string invalid_basis_file = "InvalidBasis.bas";
  std::ofstream f;
  // Write and read a file for unsupported HiGHS version
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS v0" << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::kError);

  // Write and read a file for an invalid basis
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS v1" << std::endl;
  f << "None" << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::kOk);

  // Write and read a file for incompatible number of columns
  f.open(invalid_basis_file, std::ios::out);
  f << "HiGHS v1" << std::endl;
  f << "Valid" << std::endl;
  f << "# Columns " << highs.getNumCol() - 1 << std::endl;
  f.close();
  return_status = highs.readBasis(invalid_basis_file);
  REQUIRE(return_status == HighsStatus::kError);

  testBasisRestart(highs, true);
  testBasisReloadModel(highs, true);

  std::remove(basis_file.c_str());
  std::remove(invalid_basis_file.c_str());
}

// No commas in test case name.
TEST_CASE("Basis-data", "[highs_basis_data]") {
  HighsStatus return_status;
  std::string model0_file =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  std::string model1_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  Highs highs;
  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
  }
  assert(model0_file != model1_file);

  return_status = highs.readModel(model0_file);
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  basis_data = highs.getBasis();
  REQUIRE(return_status == HighsStatus::kOk);

  testBasisRestart(highs, false);
  testBasisReloadModel(highs, false);
}

void testBasisReloadModel(Highs& highs, const bool from_file) {
  // Checks that no simplex iterations are required if a saved optimal
  // basis is used for the original LP after solving a different LP
  HighsStatus return_status;
  std::string model0_file =
      std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  std::string model1_file =
      std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  // Clear the current model
  highs.clearModel();

  // Cannot load non-trivial basis for empty model
  if (from_file) {
    return_status = highs.readBasis(basis_file);
  } else {
    return_status = highs.setBasis(basis_data);
  }
  REQUIRE(return_status == HighsStatus::kError);

  // Read and solve a different model
  highs.readModel(model1_file);
  highs.run();

  // Cannot load basis for model of different size
  if (from_file) {
    return_status = highs.readBasis(basis_file);
  } else {
    return_status = highs.setBasis(basis_data);
  }
  REQUIRE(return_status == HighsStatus::kError);

  // Clear and load original model and basis
  highs.clearModel();
  highs.readModel(model0_file);
  if (from_file) {
    return_status = highs.readBasis(basis_file);
  } else {
    return_status = highs.setBasis(basis_data);
  }
  REQUIRE(return_status == HighsStatus::kOk);

  // Ensure that no simplex iterations are required when solved from
  // the optimal basis
  highs.run();
  //  highs.writeSolution("", kSolutionStyleRaw);
  REQUIRE(highs.getInfo().simplex_iteration_count == 0);
}
void testBasisRestart(Highs& highs, const bool from_file) {
  // Checks that no simplex iterations are required if a saved optimal
  // basis is used for the original LP after changing a bound, solving
  // - so that the internal basis changes - and then restoring the
  // original LP
  HighsStatus return_status;
  // highs.writeSolution("", kSolutionStylePretty);
  // Change a bound and resolve

  const HighsLp& lp = highs.getLp();
  const HighsBasis& basis = highs.getBasis();
  const HighsSolution& solution = highs.getSolution();
  const HighsInfo& info = highs.getInfo();
  // Find the first basic variable
  HighsInt iCol;
  for (iCol = 0; iCol < lp.num_col_; iCol++) {
    if (basis.col_status[iCol] == HighsBasisStatus::kBasic) break;
  }
  assert(iCol < lp.num_col_);
  const HighsInt changeCol = iCol;
  const double old_lower_bound = lp.col_lower_[changeCol];
  const double old_upper_bound = lp.col_upper_[changeCol];
  const double new_lower_bound = solution.col_value[changeCol] + 0.1;
  highs.changeColBounds(changeCol, new_lower_bound, old_upper_bound);

  return_status = highs.run();

  if (dev_run) {
    printf("After modifying lower bound of column %" HIGHSINT_FORMAT
           " from %g to %g, solving the "
           "LP "
           "requires %" HIGHSINT_FORMAT " iterations and objective is %g\n",
           changeCol, old_lower_bound, new_lower_bound,
           info.simplex_iteration_count, highs.getObjectiveValue());
    //  highs.writeSolution("", kSolutionStylePretty);
  }
  // Make sure that the test requires iterations
  assert(info.simplex_iteration_count > 0);

  // Recover bound, load optimal basis and resolve

  highs.changeColBounds(changeCol, old_lower_bound, old_upper_bound);

  if (from_file) {
    return_status = highs.readBasis(basis_file);
  } else {
    return_status = highs.setBasis(basis_data);
  }
  REQUIRE(return_status == HighsStatus::kOk);

  return_status = highs.run();
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) {
    printf("After restoring lower bound of column %" HIGHSINT_FORMAT
           " from %g to %g, solving the "
           "LP "
           "requires %" HIGHSINT_FORMAT " iterations and objective is %g\n",
           changeCol, new_lower_bound, old_lower_bound,
           info.simplex_iteration_count, highs.getObjectiveValue());
  }

  REQUIRE(info.simplex_iteration_count == 0);
}
