#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

void quiet_run(Highs& highs) {
  highs.setHighsLogfile();
  highs.setHighsOutput();
  highs.run();
  highs.setHighsLogfile(stdout);
  highs.setHighsOutput(stdout);
}
TEST_CASE("Ranging", "[highs_test_ranging]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }

  std::string model_file;
  HighsModelStatus require_model_status;

  model_file = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  require_model_status = HighsModelStatus::OPTIMAL;
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);

  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  highs.writeSolution("", true);

  double optimal_objective = highs.getObjectiveValue();

  HighsRanging ranging;
  REQUIRE(highs.getRanging(ranging) == HighsStatus::OK);
  HighsBasis basis = highs.getBasis();
  assert(basis.valid_);
  HighsSolution solution = highs.getSolution();

  vector<HighsBasisStatus>& row_status = basis.row_status;
  vector<double>& row_value = solution.row_value;
  vector<double>& row_dual = solution.row_dual;

  HighsLp lp = highs.getLp();
  int numRow = lp.numRow_;

  double total_error = 0;
  const double total_relative_error_tolerance = 1e-15;
  const double realtive_error_tolerance = total_relative_error_tolerance;
  const double relative_error_denominator = max(1.0, fabs(optimal_objective));
  


  // Show all rowwise data
  printf(" --- Row bounds ranging ---\n");
  printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "bound^", "object^",
	 "verify^", "error^", "bound_", "object_", "verify_", "error_");
  for (int i = 0; i < numRow; i++) {
    double solved_up = 0;
    double solved_down = 0;
    double row_bound_up_value = ranging.rowBoundUp.Value_[i];
    double row_bound_up_objective = ranging.rowBoundUp.Objective_[i];
    int row_bound_up_in_col = ranging.rowBoundUp.InCol_[i];
    int row_bound_up_out_col = ranging.rowBoundUp.OutCol_[i];
    double row_bound_down_value = ranging.rowBoundDown.Value_[i];
    double row_bound_down_objective = ranging.rowBoundDown.Objective_[i];
    int row_bound_down_in_col = ranging.rowBoundDown.InCol_[i];
    int row_bound_down_out_col = ranging.rowBoundDown.OutCol_[i];
    double lower = lp.rowLower_[i];
    double upper = lp.rowUpper_[i];
    double new_lower = lower;
    double new_upper = upper;
    if (row_bound_down_value > -HIGHS_CONST_INF) {
      if (row_status[i] == HighsBasisStatus::ZERO) {
	printf("What to do if row is free??\n");
      }
      if (row_status[i] != HighsBasisStatus::BASIC) {
	// Nonbasic
	if (lower == upper) {
	  new_lower = row_bound_down_value;
	  new_upper = row_bound_down_value;
	} else if (row_status[i] == HighsBasisStatus::LOWER) {
	  new_lower = row_bound_down_value;
	} else {
	  new_upper = row_bound_down_value;
	}
      } else {
	new_upper = row_bound_down_value;
      }
      highs.changeRowBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      quiet_run(highs);
      solved_down = highs.getObjectiveValue();
    } else {
      solved_down = row_bound_down_objective;
    }
    double error = abs(solved_down-row_bound_down_objective);
    total_error += error;
    double relative_down_error = error/relative_error_denominator;
    double relative_up_error = 0;
    printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	   i, lower, upper, row_value[i], 0.0, row_dual[i],
	   row_bound_up_value, row_bound_up_objective, solved_up, relative_up_error,
	   row_bound_down_value, row_bound_down_objective, solved_down, relative_down_error);
    highs.changeRowBounds(i, lower, upper);
  }
  printf("\n\n");
  
 
}
