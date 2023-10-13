#include <fstream>

#include "HCheckConfig.h"
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

TEST_CASE("set-pathological-basis", "[highs_basis_data]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  HighsBasis basis;

  basis.clear();
  highs.clearSolver();
  highs.addCol(1.0, 0, 1, 0, nullptr, nullptr);
  HighsInt index = 0;
  double value = 1.0;
  highs.addRow(0, 1, 1, &index, &value);
  // Set up a basis with everything nonbasic. This will lead to
  // basic_index being empty when passed to
  // HFactor::setupGeneral. Previously this led to the creation of
  // pointer &basic_index[0] that caused Windows faiure referenced in
  // #1129, and reported in #1166. However, now that
  // basic_index.data() is used to create the pointer, there is no
  // Windows failure. Within HFactor::setupGeneral and
  // HFactor::build(), the empty list of basic variables is handled
  // correctly - with a basis of logicals being created
  basis.col_status.push_back(HighsBasisStatus::kLower);
  basis.row_status.push_back(HighsBasisStatus::kLower);
  highs.setBasis(basis);
  highs.run();
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  basis.clear();
  highs.clearModel();
  highs.addCol(1.0, -kHighsInf, kHighsInf, 0, nullptr, nullptr);
  basis.col_status.push_back(HighsBasisStatus::kBasic);
  highs.setBasis(basis);
  highs.run();
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
}

TEST_CASE("Basis-no-basic", "[highs_basis_data]") {
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  highs.addCol(1.0, 0, 1, 0, nullptr, nullptr);
  highs.addCol(-1.0, 0, 1, 0, nullptr, nullptr);
  std::vector<HighsInt> index = {0, 1};
  std::vector<double> value = {1.0, 2.0};
  highs.addRow(0, 1, 2, index.data(), value.data());
  value[0] = -1.0;
  highs.addRow(0, 1, 2, index.data(), value.data());
  // Make all variables basic. This is a 2-row version of
  // set-pathological-basis
  HighsBasis basis;
  basis.col_status.push_back(HighsBasisStatus::kLower);
  basis.col_status.push_back(HighsBasisStatus::kLower);
  basis.row_status.push_back(HighsBasisStatus::kLower);
  basis.row_status.push_back(HighsBasisStatus::kLower);
  REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  REQUIRE(highs.getInfo().objective_function_value == -0.5);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);
}

// No commas in test case name.
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
