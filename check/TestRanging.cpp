#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

void columnHeader() {
  printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "cost", "dual", "bound^", "object^",
	 "verify^", "error^", "bound_", "object_", "verify_", "error_");
}
void quietRun(Highs& highs) {
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
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  /*
  std::string model_file = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);
  require_model_status = HighsModelStatus::OPTIMAL;
  */
  SpecialLps special_lps;
  special_lps.blendingLp(lp, require_model_status, optimal_objective);
  highs.passModel(lp);

  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);

  REQUIRE(highs.getModelStatus() == require_model_status);
  REQUIRE(highs.getObjectiveValue() == optimal_objective);

  highs.writeSolution("", true);

  HighsRanging ranging;
  REQUIRE(highs.getRanging(ranging) == HighsStatus::OK);
  HighsBasis basis = highs.getBasis();
  assert(basis.valid_);
  HighsSolution solution = highs.getSolution();

  vector<HighsBasisStatus>& row_status = basis.row_status;
  vector<double>& row_value = solution.row_value;
  vector<double>& row_dual = solution.row_dual;

  lp = highs.getLp();
  int numRow = lp.numRow_;

  double total_error = 0;
  const double total_relative_error_tolerance = 1e-15;
  const double realtive_error_tolerance = total_relative_error_tolerance;
  const double relative_error_denominator = max(1.0, fabs(optimal_objective));
  


  // Show all rowwise data
  printf(" --- Row bounds ranging ---\n");
  int to_row = 1;
  for (int i = 0; i < to_row; i++) {
    double solved_up = 0;
    double solved_dn = 0;
    double row_bound_up_value = ranging.row_bound_up.value_[i];
    double row_bound_up_objective = ranging.row_bound_up.objective_[i];
    int row_bound_up_in_col = ranging.row_bound_up.in_col_[i];
    int row_bound_up_ou_col = ranging.row_bound_up.ou_col_[i];
    double row_bound_dn_value = ranging.row_bound_dn.value_[i];
    double row_bound_dn_objective = ranging.row_bound_dn.objective_[i];
    int row_bound_dn_in_col = ranging.row_bound_dn.in_col_[i];
    int row_bound_dn_ou_col = ranging.row_bound_dn.ou_col_[i];
    double lower = lp.rowLower_[i];
    double upper = lp.rowUpper_[i];
    double new_lower = lower;
    double new_upper = upper;
    if (row_bound_dn_value > -HIGHS_CONST_INF) {
      if (row_status[i] == HighsBasisStatus::ZERO) {
	printf("What to do if row is free??\n");
      }
      if (row_status[i] != HighsBasisStatus::BASIC) {
	// Nonbasic
	if (lower == upper) {
	  new_lower = row_bound_dn_value;
	  new_upper = row_bound_dn_value;
	} else if (row_status[i] == HighsBasisStatus::LOWER) {
	  new_lower = row_bound_dn_value;
	} else {
	  new_upper = row_bound_dn_value;
	}
      } else {
	new_upper = row_bound_dn_value;
      }
      highs.changeRowBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      quietRun(highs);
      REQUIRE(highs.getModelStatus() == require_model_status);
      solved_dn = highs.getObjectiveValue();
      highs.writeSolution("", true);

    } else {
      solved_dn = row_bound_dn_objective;
    }
    double error = abs(solved_dn-row_bound_dn_objective);
    total_error += error;
    double relative_dn_error = error/relative_error_denominator;
    double relative_up_error = 0;
    columnHeader();
    printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	   i, lower, upper, row_value[i], 0.0, row_dual[i],
	   row_bound_up_value, row_bound_up_objective, solved_up, relative_up_error,
	   row_bound_dn_value, row_bound_dn_objective, solved_dn, relative_dn_error);
    highs.changeRowBounds(i, lower, upper);
  }
  printf("\n\n");
  
 
}
