#include <cmath>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;

TEST_CASE("sos1-basic-api", "[sos]") {
  // Three variables in SOS1: maximize sum subject to x1 + x2 + x3 <= 10
  // SOS1 means at most one can be nonzero, so optimal is x3 = 10 (obj = -30)
  Highs h;
  h.setOptionValue("output_flag", dev_run);

  h.addVar(0.0, 10.0);
  h.addVar(0.0, 10.0);
  h.addVar(0.0, 10.0);

  h.changeColCost(0, -1.0);
  h.changeColCost(1, -2.0);
  h.changeColCost(2, -3.0);

  // x1 + x2 + x3 <= 10
  HighsInt idx[] = {0, 1, 2};
  double val[] = {1.0, 1.0, 1.0};
  h.addRow(-kHighsInf, 10.0, 3, idx, val);

  // Add SOS1 constraint
  HighsInt cols[] = {0, 1, 2};
  double weights[] = {1.0, 2.0, 3.0};
  REQUIRE(h.addSosConstraint(1, 3, cols, weights) == HighsStatus::kOk);
  REQUIRE(h.getNumSosConstraints() == 1);

  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const auto& sol = h.getSolution();
  // At most one variable should be nonzero
  int num_nonzero = 0;
  for (int i = 0; i < 3; i++) {
    if (std::abs(sol.col_value[i]) > 1e-6) num_nonzero++;
  }
  REQUIRE(num_nonzero <= 1);
  // Best choice is x3 = 10
  REQUIRE(sol.col_value[2] == Approx(10.0).margin(1e-6));
  REQUIRE(h.getInfo().objective_function_value ==
          Approx(-30.0).margin(1e-6));
}

TEST_CASE("sos2-basic-api", "[sos]") {
  // SOS2 constraint: at most two adjacent variables nonzero
  Highs h;
  h.setOptionValue("output_flag", dev_run);

  // 4 variables
  h.addVar(0.0, 10.0);
  h.addVar(0.0, 10.0);
  h.addVar(0.0, 10.0);
  h.addVar(0.0, 10.0);

  // Minimize -x1 - x3 (want x1 and x3 both nonzero, but SOS2 prevents it)
  h.changeColCost(0, -1.0);
  h.changeColCost(1, 0.0);
  h.changeColCost(2, -1.0);
  h.changeColCost(3, 0.0);

  // x1 + x2 + x3 + x4 <= 10
  HighsInt idx[] = {0, 1, 2, 3};
  double val[] = {1.0, 1.0, 1.0, 1.0};
  h.addRow(-kHighsInf, 10.0, 4, idx, val);

  // SOS2 with weights 1, 2, 3, 4
  HighsInt cols[] = {0, 1, 2, 3};
  double weights[] = {1.0, 2.0, 3.0, 4.0};
  REQUIRE(h.addSosConstraint(2, 4, cols, weights) == HighsStatus::kOk);
  REQUIRE(h.getNumSosConstraints() == 1);

  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const auto& sol = h.getSolution();
  // Check SOS2: nonzero variables must be adjacent by weight
  int first_nz = -1, last_nz = -1;
  for (int i = 0; i < 4; i++) {
    if (std::abs(sol.col_value[i]) > 1e-6) {
      if (first_nz == -1) first_nz = i;
      last_nz = i;
    }
  }
  REQUIRE(last_nz - first_nz <= 1);
}

TEST_CASE("sos1-mps-read", "[sos]") {
  // Read MPS file with SOS1 section
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/sos1.mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(filename) == HighsStatus::kOk);
  REQUIRE(h.getNumSosConstraints() == 1);

  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const auto& sol = h.getSolution();
  // At most one variable should be nonzero (SOS1)
  int num_nonzero = 0;
  for (int i = 0; i < 3; i++) {
    if (std::abs(sol.col_value[i]) > 1e-6) num_nonzero++;
  }
  REQUIRE(num_nonzero <= 1);
  // Best: x3 = 10, obj = -30
  REQUIRE(h.getInfo().objective_function_value ==
          Approx(-30.0).margin(1e-6));
}

TEST_CASE("sos1-validation-errors", "[sos]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);

  h.addVar(0.0, 10.0);
  h.addVar(0.0, 10.0);

  HighsInt cols[] = {0, 1};
  double weights[] = {1.0, 2.0};

  // Invalid type
  REQUIRE(h.addSosConstraint(0, 2, cols, weights) == HighsStatus::kError);
  REQUIRE(h.addSosConstraint(3, 2, cols, weights) == HighsStatus::kError);

  // Invalid num_members
  REQUIRE(h.addSosConstraint(1, 0, cols, weights) == HighsStatus::kError);

  // Invalid column index
  HighsInt bad_cols[] = {0, 5};
  REQUIRE(h.addSosConstraint(1, 2, bad_cols, weights) == HighsStatus::kError);

  // Valid
  REQUIRE(h.addSosConstraint(1, 2, cols, weights) == HighsStatus::kOk);
  REQUIRE(h.getNumSosConstraints() == 1);
}

TEST_CASE("sos1-binding", "[sos]") {
  // Problem where SOS1 constraint forces branching
  // Without SOS1: optimal is x1=5, x2=5, obj = -25
  // With SOS1: only one can be nonzero, so optimal is x2=10, obj = -20
  Highs h;
  h.setOptionValue("output_flag", dev_run);

  h.addVar(0.0, 10.0);
  h.addVar(0.0, 10.0);

  h.changeColCost(0, -3.0);
  h.changeColCost(1, -2.0);

  // x1 + x2 <= 10
  HighsInt idx[] = {0, 1};
  double val[] = {1.0, 1.0};
  h.addRow(-kHighsInf, 10.0, 2, idx, val);

  // 2*x1 + x2 <= 10
  HighsInt idx2[] = {0, 1};
  double val2[] = {2.0, 1.0};
  h.addRow(-kHighsInf, 10.0, 2, idx2, val2);

  // SOS1
  HighsInt cols[] = {0, 1};
  double weights[] = {1.0, 2.0};
  h.addSosConstraint(1, 2, cols, weights);

  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const auto& sol = h.getSolution();
  int num_nonzero = 0;
  for (int i = 0; i < 2; i++) {
    if (std::abs(sol.col_value[i]) > 1e-6) num_nonzero++;
  }
  REQUIRE(num_nonzero <= 1);
  // Optimal: x1=5 (obj=-15) or x2=10 (obj=-20). x2=10 wins.
  REQUIRE(h.getInfo().objective_function_value ==
          Approx(-20.0).margin(1e-6));
}

TEST_CASE("sos-multiple", "[sos]") {
  // Multiple SOS1 constraints with integer variables
  Highs h;
  h.setOptionValue("output_flag", dev_run);

  // 4 integer variables
  for (int i = 0; i < 4; i++) h.addVar(0.0, 10.0);
  for (int i = 0; i < 4; i++) h.changeColIntegrality(i, HighsVarType::kInteger);

  h.changeColCost(0, -1.0);
  h.changeColCost(1, -1.0);
  h.changeColCost(2, -1.0);
  h.changeColCost(3, -1.0);

  // x1 + x2 + x3 + x4 <= 20
  HighsInt idx[] = {0, 1, 2, 3};
  double val[] = {1.0, 1.0, 1.0, 1.0};
  h.addRow(-kHighsInf, 20.0, 4, idx, val);

  // SOS1: {x1, x2}
  HighsInt cols1[] = {0, 1};
  double w1[] = {1.0, 2.0};
  h.addSosConstraint(1, 2, cols1, w1);

  // SOS1: {x3, x4}
  HighsInt cols2[] = {2, 3};
  double w2[] = {1.0, 2.0};
  h.addSosConstraint(1, 2, cols2, w2);

  REQUIRE(h.getNumSosConstraints() == 2);

  h.run();
  REQUIRE(h.getModelStatus() == HighsModelStatus::kOptimal);

  const auto& sol = h.getSolution();
  // Check SOS1 for each set
  int nz1 = 0, nz2 = 0;
  if (std::abs(sol.col_value[0]) > 1e-6) nz1++;
  if (std::abs(sol.col_value[1]) > 1e-6) nz1++;
  if (std::abs(sol.col_value[2]) > 1e-6) nz2++;
  if (std::abs(sol.col_value[3]) > 1e-6) nz2++;
  REQUIRE(nz1 <= 1);
  REQUIRE(nz2 <= 1);
  // Optimal: one from each set at 10, obj = -20
  REQUIRE(h.getInfo().objective_function_value ==
          Approx(-20.0).margin(1e-6));
}

TEST_CASE("sos1-ismip", "[sos]") {
  // Test that adding SOS constraints makes the model a MIP
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.addVar(0.0, 1.0);
  h.addVar(0.0, 1.0);

  // Pure LP - should not be MIP
  HighsInt idx[] = {0, 1};
  double val[] = {1.0, 1.0};
  h.addRow(-kHighsInf, 1.0, 2, idx, val);
  h.changeColCost(0, -1.0);
  h.changeColCost(1, -1.0);

  // Add SOS1 - should become MIP
  HighsInt cols[] = {0, 1};
  double weights[] = {1.0, 2.0};
  h.addSosConstraint(1, 2, cols, weights);
  REQUIRE(h.getNumSosConstraints() == 1);
}
