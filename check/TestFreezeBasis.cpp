#include "Highs.h"
#include "catch.hpp"

const double inf = kHighsInf;
const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

TEST_CASE("FreezeBasis", "[highs_test_freeze_basis]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  const HighsLp& lp = highs.getLp();
  const HighsInt num_col = lp.num_col_;
  const HighsInt num_row = lp.num_row_;
  const HighsInt from_col = 0;
  const HighsInt to_col = num_col - 1;
  const HighsInfo& info = highs.getInfo();

  vector<double> original_col_lower = lp.col_lower_;
  vector<double> original_col_upper = lp.col_upper_;

  // Get the continuous solution and iteration count from a logical basis
  highs.run();
  HighsInt continuous_iteration_count = info.simplex_iteration_count;
  double continuous_objective = info.objective_function_value;

  // Get the integer solution to provide bound tightenings
  vector<HighsVarType> integrality;
  integrality.assign(num_col, HighsVarType::kInteger);
  highs.changeColsIntegrality(from_col, to_col, &integrality[0]);

  highs.setOptionValue("output_flag", false);
  highs.run();
  if (dev_run) highs.setOptionValue("output_flag", true);

  vector<double> integer_solution = highs.getSolution().col_value;
  vector<double> local_col_lower;
  vector<double> local_col_upper;

  // Now restore the original integrality and set an explicit logical
  // basis to force reinversion
  integrality.assign(num_col, HighsVarType::kContinuous);
  highs.changeColsIntegrality(from_col, to_col, &integrality[0]);
  HighsBasis basis;
  for (HighsInt iCol = 0; iCol < num_col; iCol++)
    basis.col_status.push_back(HighsBasisStatus::kLower);
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    basis.row_status.push_back(HighsBasisStatus::kBasic);
  highs.setBasis(basis);

  HighsInt frozen_basis_id0;
  // Cannot freeze a basis when there's no INVERT
  REQUIRE(highs.freezeBasis(frozen_basis_id0) == HighsStatus::kError);

  // Get the basic variables to force INVERT with a logical basis
  vector<HighsInt> basic_variables;
  basic_variables.resize(lp.num_row_);
  highs.getBasicVariables(&basic_variables[0]);

  // Can freeze a basis now!
  REQUIRE(highs.freezeBasis(frozen_basis_id0) == HighsStatus::kOk);

  if (dev_run) {
    highs.setOptionValue("output_flag", true);
    highs.setOptionValue("log_dev_level", 2);
    highs.setOptionValue("highs_debug_level", 3);
  }
  highs.run();
  // highs.setOptionValue("output_flag", false);

  // Now freeze the current basis and add the integer solution as lower bounds
  HighsInt frozen_basis_id1;
  REQUIRE(highs.freezeBasis(frozen_basis_id1) == HighsStatus::kOk);
  if (dev_run) printf("\nSolving with bounds (integer solution, upper)\n");
  local_col_lower = integer_solution;
  local_col_upper = original_col_upper;
  highs.changeColsBounds(from_col, to_col, &local_col_lower[0],
                         &local_col_upper[0]);
  if (dev_run) highs.setOptionValue("output_flag", true);
  highs.run();
  // highs.setOptionValue("output_flag", false);
  double semi_continuous_objective = info.objective_function_value;

  // Now freeze the current basis and add the integer solution as upper bounds
  HighsInt frozen_basis_id2;
  REQUIRE(highs.freezeBasis(frozen_basis_id2) == HighsStatus::kOk);
  if (dev_run)
    printf("\nSolving with bounds (integer solution, integer solution)\n");
  local_col_lower = integer_solution;
  local_col_upper = integer_solution;
  highs.changeColsBounds(from_col, to_col, &local_col_lower[0],
                         &local_col_upper[0]);
  if (dev_run) highs.setOptionValue("output_flag", true);
  highs.run();
  // highs.setOptionValue("output_flag", false);

  // Now test unfreeze
  //
  // Restore the upper bounds and unfreeze frozen_basis_id2 (semi-continuous
  // optimal basis)
  if (dev_run)
    printf(
        "\nSolving with bounds (integer solution, upper) after unfreezing "
        "basis %d\n",
        (int)frozen_basis_id2);
  if (dev_run) printf("Change column bounds to (integer solution, upper)\n");
  local_col_lower = integer_solution;
  local_col_upper = original_col_upper;
  highs.changeColsBounds(from_col, to_col, &local_col_lower[0],
                         &local_col_upper[0]);
  if (dev_run) printf("Unfreeze basis %d\n", (int)frozen_basis_id2);
  REQUIRE(highs.unfreezeBasis(frozen_basis_id2) == HighsStatus::kOk);
  // Solving the LP should require no iterations
  if (dev_run) highs.setOptionValue("output_flag", true);
  highs.run();
  // highs.setOptionValue("output_flag", false);
  REQUIRE(info.simplex_iteration_count == 0);
  double dl_objective =
      fabs(semi_continuous_objective - info.objective_function_value);
  REQUIRE(dl_objective < double_equal_tolerance);
  REQUIRE(info.simplex_iteration_count == 0);
  //
  // Restore the lower bounds and unfreeze frozen_basis_id1 (continuous optimal
  // basis)
  if (dev_run)
    printf("\nSolving with bounds (lower, upper) after unfreezing basis %d\n",
           (int)frozen_basis_id1);
  if (dev_run) printf("Change column bounds to (lower, upper)\n");
  local_col_lower = original_col_lower;
  local_col_upper = original_col_upper;
  highs.changeColsBounds(from_col, to_col, &local_col_lower[0],
                         &local_col_upper[0]);
  if (dev_run) printf("Unfreeze basis %d\n", (int)frozen_basis_id1);
  REQUIRE(highs.unfreezeBasis(frozen_basis_id1) == HighsStatus::kOk);
  // Solving the LP should require no iterations
  if (dev_run) highs.setOptionValue("output_flag", true);
  highs.run();
  // highs.setOptionValue("output_flag", false);
  dl_objective = fabs(continuous_objective - info.objective_function_value);
  REQUIRE(dl_objective < double_equal_tolerance);
  REQUIRE(info.simplex_iteration_count == 0);
  //
  // Unfreeze frozen_basis_id0
  if (dev_run)
    printf("\nSolving with bounds (lower, upper) after unfreezing basis %d\n",
           (int)frozen_basis_id0);
  if (dev_run) printf("Unfreeze basis %d\n", (int)frozen_basis_id0);
  REQUIRE(highs.unfreezeBasis(frozen_basis_id0) == HighsStatus::kOk);
  // Solving the LP should require the same number of iterations as before
  if (dev_run) highs.setOptionValue("output_flag", true);
  highs.run();
  // highs.setOptionValue("output_flag", false);
  dl_objective = fabs(continuous_objective - info.objective_function_value);
  REQUIRE(dl_objective < double_equal_tolerance);
  REQUIRE(info.simplex_iteration_count == continuous_iteration_count);

  REQUIRE(highs.frozenBasisAllDataClear() == HighsStatus::kOk);
}
