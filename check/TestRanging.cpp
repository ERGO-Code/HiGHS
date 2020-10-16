#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;

void quietRun(Highs& highs) {
  highs.setHighsLogfile();
  highs.setHighsOutput();
  highs.run();
  highs.setHighsLogfile(stdout);
  highs.setHighsOutput(stdout);
}

void colCostColumnHeader() {
  if (dev_run) printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "cost", "dual", "cost^", "object^",
	 "verify^", "error^", "cost_", "object_", "verify_", "error_");
}

void colBoundcolumnHeader() {
  if (dev_run) printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "dual", "bound^", "object^",
	 "verify^", "error^", "bound_", "object_", "verify_", "error_");
}

void rowBoundColumnHeader() {
  if (dev_run) printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
	 "lower", "upper", "value", "dual", "bound^", "object^",
	 "verify^", "error^", "bound_", "object_", "verify_", "error_");
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

  const bool from_file = true;
  if (from_file) {
    std::string model_file = std::string(HIGHS_DIR) + "/check/instances/stair.mps";
    REQUIRE(highs.readModel(model_file) == HighsStatus::OK);
    require_model_status = HighsModelStatus::OPTIMAL;
  } else {
    SpecialLps special_lps;
    special_lps.blendingLp(lp, require_model_status, optimal_objective);
    highs.passModel(lp);
  }

  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);

  REQUIRE(highs.getModelStatus() == require_model_status);
  optimal_objective = highs.getObjectiveValue();

  HighsRanging ranging;
  REQUIRE(highs.getRanging(ranging) == HighsStatus::OK);
  HighsBasis basis = highs.getBasis();
  assert(basis.valid_);
  HighsSolution solution = highs.getSolution();

  vector<HighsBasisStatus>& col_status = basis.col_status;
  vector<HighsBasisStatus>& row_status = basis.row_status;
  vector<double>& col_value = solution.col_value;
  vector<double>& col_dual = solution.col_dual;
  vector<double>& row_value = solution.row_value;
  vector<double>& row_dual = solution.row_dual;

  lp = highs.getLp();
  int numRow = lp.numRow_;
  int numCol = lp.numCol_;

  double total_error = 0;
  const double total_relative_error_tolerance = 1e-12;
  const double relative_error_tolerance = total_relative_error_tolerance;
  const double relative_error_denominator = max(1.0, fabs(optimal_objective));
  const double initial_error_report_threshold = 1e-16;//relative_error_tolerance;
  double error_report_threshold;

  double max_col_cost_relative_error = 0;
  int max_col_cost_relative_error_col = -1;
  double max_col_bound_relative_error = 0;
  int max_col_bound_relative_error_col = -1;
  double max_row_bound_relative_error = 0;
  int max_row_bound_relative_error_row = -1;

  if (dev_run) printf(" --- Col cost ranging ---\n");
  bool small_numCol = numCol<10;
  error_report_threshold = initial_error_report_threshold;
  colCostColumnHeader();
  for (int i = 0; i < numCol; i++) {
    double col_cost_up_value = ranging.col_cost_up.value_[i];
    double col_cost_up_objective = ranging.col_cost_up.objective_[i];
    double col_cost_dn_value = ranging.col_cost_dn.value_[i];
    double col_cost_dn_objective = ranging.col_cost_dn.objective_[i];
    double cost = lp.colCost_[i];
    double new_cost;
    double solved_up = 0;
    double solved_dn = 0;
    double error;
    // Col cost up ranging
    if (col_cost_up_value < HIGHS_CONST_INF) {
      new_cost = col_cost_up_value;
      highs.changeColCost(i, new_cost);
      highs.setBasis(basis);
      quietRun(highs);
      REQUIRE(highs.getModelStatus() == require_model_status);
      solved_up = highs.getObjectiveValue();
      highs.changeColCost(i, cost);
      error = fabs(solved_up-col_cost_up_objective);
    } else {
      solved_up = col_cost_up_objective;
      error = 0;
    }
    error = fabs(solved_up-col_cost_up_objective);
    total_error += error;
    double relative_up_error = error/relative_error_denominator;
    REQUIRE(relative_up_error < relative_error_tolerance);
    // Col cost down ranging
    if (col_cost_dn_value > -HIGHS_CONST_INF) {
      new_cost = col_cost_dn_value;
      highs.changeColCost(i, new_cost);
      highs.setBasis(basis);
      quietRun(highs);
      REQUIRE(highs.getModelStatus() == require_model_status);
      solved_dn = highs.getObjectiveValue();
      highs.changeColCost(i, cost);
      error = fabs(solved_dn-col_cost_dn_objective);
    } else {
      solved_dn = col_cost_dn_objective;
      error = 0;
    }
    total_error += error;
    double relative_dn_error = error/relative_error_denominator;
    REQUIRE(relative_dn_error < relative_error_tolerance);
    double relative_error = max(relative_up_error, relative_dn_error);
    if (relative_error > max_col_cost_relative_error) {
      max_col_cost_relative_error = relative_error;
      max_col_cost_relative_error_col = i;
    }
    if (small_numCol || relative_error > error_report_threshold) {
      if (dev_run) printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	     i, cost, col_dual[i],
	     col_cost_up_value, col_cost_up_objective, solved_up, relative_up_error,
	     col_cost_dn_value, col_cost_dn_objective, solved_dn, relative_dn_error);
      error_report_threshold = 10*error_report_threshold;
    }
  }

  if (dev_run) printf(" --- Col bounds ranging ---\n");
  error_report_threshold = initial_error_report_threshold;
  colBoundcolumnHeader();
  for (int i = 0; i < numCol; i++) {
    double col_bound_up_value = ranging.col_bound_up.value_[i];
    double col_bound_up_objective = ranging.col_bound_up.objective_[i];
    double col_bound_dn_value = ranging.col_bound_dn.value_[i];
    double col_bound_dn_objective = ranging.col_bound_dn.objective_[i];
    double lower = lp.colLower_[i];
    double upper = lp.colUpper_[i];
    double new_lower;
    double new_upper;
    double solved_up = 0;
    double solved_dn = 0;
    double error;
    // Col bound up ranging
    if (col_bound_up_value < HIGHS_CONST_INF) {
      // Free cols should not have a finite col_bound_up_value
      assert(col_status[i] != HighsBasisStatus::ZERO);
      new_lower = lower;
      new_upper = upper;
      if (col_status[i] != HighsBasisStatus::BASIC) {
	// Nonbasic
	if (lower == upper) {
	  new_lower = col_bound_up_value;
	  new_upper = col_bound_up_value;
	} else if (col_status[i] == HighsBasisStatus::LOWER) {
	  new_lower = col_bound_up_value;
	} else {
	  new_upper = col_bound_up_value;
	}
      } else {
	new_lower = col_bound_up_value;
      }
      assert(new_lower <= new_upper);
      highs.changeColBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      quietRun(highs);
      REQUIRE(highs.getModelStatus() == require_model_status);
      solved_up = highs.getObjectiveValue();
      highs.changeColBounds(i, lower, upper);
      error = fabs(solved_up-col_bound_up_objective);
    } else {
      solved_up = col_bound_up_objective;
      error = 0;
    }
    error = fabs(solved_up-col_bound_up_objective);
    total_error += error;
    double relative_up_error = error/relative_error_denominator;
    REQUIRE(relative_up_error < relative_error_tolerance);
    // Col bound down ranging
    if (col_bound_dn_value > -HIGHS_CONST_INF) {
      // Free cols should not have a finite col_bound_dn_value
      assert(col_status[i] != HighsBasisStatus::ZERO);
      new_lower = lower;
      new_upper = upper;
      if (col_status[i] != HighsBasisStatus::BASIC) {
	// Nonbasic
	if (lower == upper) {
	  new_lower = col_bound_dn_value;
	  new_upper = col_bound_dn_value;
	} else if (col_status[i] == HighsBasisStatus::LOWER) {
	  new_lower = col_bound_dn_value;
	} else {
	  new_upper = col_bound_dn_value;
	}
      } else {
	new_upper = col_bound_dn_value;
      }
      assert(new_lower <= new_upper);
      highs.changeColBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      quietRun(highs);
      REQUIRE(highs.getModelStatus() == require_model_status);
      solved_dn = highs.getObjectiveValue();
      highs.changeColBounds(i, lower, upper);
      error = fabs(solved_dn-col_bound_dn_objective);
    } else {
      solved_dn = col_bound_dn_objective;
      error = 0;
    }
    total_error += error;
    double relative_dn_error = error/relative_error_denominator;
    REQUIRE(relative_dn_error < relative_error_tolerance);
    double relative_error = max(relative_up_error, relative_dn_error);
    if (relative_error > max_col_bound_relative_error) {
      max_col_bound_relative_error = relative_error;
      max_col_bound_relative_error_col = i;
    }
    if (small_numCol || relative_error > error_report_threshold) {
      if (dev_run) printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	     i, lower, upper, col_value[i], col_dual[i],
	     col_bound_up_value, col_bound_up_objective, solved_up, relative_up_error,
	     col_bound_dn_value, col_bound_dn_objective, solved_dn, relative_dn_error);
      error_report_threshold = 10*error_report_threshold;
    }
  }
  if (dev_run) printf(" --- Row bounds ranging ---\n");
  bool small_numRow = numRow <10;
  error_report_threshold = initial_error_report_threshold;
  rowBoundColumnHeader();
  for (int i = 0; i < numRow; i++) {
    double row_bound_up_value = ranging.row_bound_up.value_[i];
    double row_bound_up_objective = ranging.row_bound_up.objective_[i];
    double row_bound_dn_value = ranging.row_bound_dn.value_[i];
    double row_bound_dn_objective = ranging.row_bound_dn.objective_[i];
    double lower = lp.rowLower_[i];
    double upper = lp.rowUpper_[i];
    double new_lower;
    double new_upper;
    double solved_up = 0;
    double solved_dn = 0;
    double error;
    // Row bound up ranging
    if (row_bound_up_value < HIGHS_CONST_INF) {
      // Free rows should not have a finite row_bound_up_value
      assert(row_status[i] != HighsBasisStatus::ZERO);
      new_lower = lower;
      new_upper = upper;
      if (row_status[i] != HighsBasisStatus::BASIC) {
	// Nonbasic
	if (lower == upper) {
	  new_lower = row_bound_up_value;
	  new_upper = row_bound_up_value;
	} else if (row_status[i] == HighsBasisStatus::LOWER) {
	  new_lower = row_bound_up_value;
	} else {
	  new_upper = row_bound_up_value;
	}
      } else {
	new_lower = row_bound_up_value;
      }
      assert(new_lower <= new_upper);
      highs.changeRowBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      quietRun(highs);
      REQUIRE(highs.getModelStatus() == require_model_status);
      solved_up = highs.getObjectiveValue();
      highs.changeRowBounds(i, lower, upper);
      error = fabs(solved_up-row_bound_up_objective);
    } else {
      solved_up = row_bound_up_objective;
      error = 0;
    }
    error = fabs(solved_up-row_bound_up_objective);
    total_error += error;
    double relative_up_error = error/relative_error_denominator;
    REQUIRE(relative_up_error < relative_error_tolerance);
    // Row bound down ranging
    if (row_bound_dn_value > -HIGHS_CONST_INF) {
      // Free rows should not have a finite row_bound_dn_value
      assert(row_status[i] != HighsBasisStatus::ZERO);
      new_lower = lower;
      new_upper = upper;
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
      assert(new_lower <= new_upper);
      highs.changeRowBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      quietRun(highs);
      REQUIRE(highs.getModelStatus() == require_model_status);
      solved_dn = highs.getObjectiveValue();
      highs.changeRowBounds(i, lower, upper);
      error = fabs(solved_dn-row_bound_dn_objective);
    } else {
      solved_dn = row_bound_dn_objective;
      error = 0;
    }
    total_error += error;
    double relative_dn_error = error/relative_error_denominator;
    REQUIRE(relative_dn_error < relative_error_tolerance);
    double relative_error = max(relative_up_error, relative_dn_error);
    if (relative_error > max_row_bound_relative_error) {
      max_row_bound_relative_error = relative_error;
      max_row_bound_relative_error_row = i;
    }
    if (small_numRow || relative_error > error_report_threshold) {
      if (dev_run) printf("%3d %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
	     i, lower, upper, row_value[i], row_dual[i],
	     row_bound_up_value, row_bound_up_objective, solved_up, relative_up_error,
	     row_bound_dn_value, row_bound_dn_objective, solved_dn, relative_dn_error);
      error_report_threshold = 10*error_report_threshold;
    }
  }
  if (dev_run) {
    if (max_col_cost_relative_error_col >= 0) 
      printf("\nMax col cost  relative objective error = %g in col %d\n", max_col_cost_relative_error,
			  max_col_cost_relative_error_col);
    if (max_col_bound_relative_error_col >= 0) 
      printf("\nMax col bound relative objective error = %g in col %d\n", max_col_bound_relative_error,
			  max_col_bound_relative_error_col);
    if (max_row_bound_relative_error_row >= 0) 
      printf("\nMax row bound relative objective error = %g in row %d\n", max_row_bound_relative_error,
			  max_row_bound_relative_error_row);
  }
 
}
