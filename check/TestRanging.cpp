#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

using std::min;

const bool dev_run = false;

HighsStatus quietRun(Highs& highs);
void colCostColumnHeader();
void colBoundcolumnHeader();
void rowBoundColumnHeader();
void assessNewBounds(double& lower, double& upper);
bool modelStatusOk(Highs& highs);
void testRanging(Highs& highs);

TEST_CASE("Ranging-min", "[highs_test_ranging]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;

  const bool from_file = true;
  if (from_file) {
    std::string model_file =
        std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
    REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  } else {
    SpecialLps special_lps;
    special_lps.blendingLp(lp, require_model_status, optimal_objective);
    highs.passModel(lp);
  }
  testRanging(highs);
}

TEST_CASE("Ranging-max", "[highs_test_ranging]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;

  const bool from_file = true;
  if (from_file) {
    std::string model_file =
        std::string(HIGHS_DIR) + "/check/instances/afiro.mps";
    REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
    REQUIRE(highs.changeObjectiveSense(ObjSense::kMaximize) ==
            HighsStatus::kOk);
  } else {
    SpecialLps special_lps;
    special_lps.blendingMaxLp(lp, require_model_status, optimal_objective);
    highs.passModel(lp);
  }
  testRanging(highs);
}

HighsStatus quietRun(Highs& highs) {
  highs.setOptionValue("output_flag", false);
  HighsStatus call_status = highs.run();
  if (dev_run) highs.setOptionValue("output_flag", true);
  return call_status;
}

void colCostColumnHeader() {
  if (dev_run)
    printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", "cost",
           "dual", "cost^", "object^", "verify^", "error^", "cost_", "object_",
           "verify_", "error_");
}

void colBoundcolumnHeader() {
  if (dev_run)
    printf("Col %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
           "lower", "upper", "value", "dual", "bound^", "object^", "verify^",
           "error^", "bound_", "object_", "verify_", "error_");
}

void rowBoundColumnHeader() {
  if (dev_run)
    printf("Row %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n",
           "lower", "upper", "value", "dual", "bound^", "object^", "verify^",
           "error^", "bound_", "object_", "verify_", "error_");
}

void assessNewBounds(double& lower, double& upper) {
  double difference = lower - upper;
  if (difference > 0) {
    if (dev_run)
      printf("New bounds [%g, %g] inconsistent with difference %g\n", lower,
             upper, difference);
    if (difference > 1e-10) assert(lower < upper);
    double average = (lower + upper) * 0.5;
    lower = average;
    upper = average;
  }
}

bool modelStatusOk(Highs& highs) {
  if (highs.getModelStatus() == HighsModelStatus::kOptimal) return true;
  if (highs.getModelStatus(true) == HighsModelStatus::kOptimal) return true;
  return false;
}

void testRanging(Highs& highs) {
  HighsLp lp;
  double optimal_objective;

  REQUIRE(highs.setBasis() == HighsStatus::kOk);
  highs.run();

  REQUIRE(modelStatusOk(highs));
  optimal_objective = highs.getObjectiveValue();

  HighsRanging ranging;
  REQUIRE(highs.getRanging(ranging) == HighsStatus::kOk);
  HighsBasis basis = highs.getBasis();
  assert(basis.valid);
  HighsSolution solution = highs.getSolution();

  vector<HighsBasisStatus>& col_status = basis.col_status;
  vector<HighsBasisStatus>& row_status = basis.row_status;
  vector<double>& col_value = solution.col_value;
  vector<double>& col_dual = solution.col_dual;
  vector<double>& row_value = solution.row_value;
  vector<double>& row_dual = solution.row_dual;

  lp = highs.getLp();
  HighsInt numRow = lp.num_row_;
  HighsInt numCol = lp.num_col_;

  const double relative_error_tolerance = 1e-10;
  const double relative_error_denominator = max(1.0, fabs(optimal_objective));
  const double initial_error_report_threshold = relative_error_tolerance / 1e6;
  double error_report_threshold;
  HighsInt num_relative_error = 0;
  double sum_relative_error = 0;
  double max_relative_error = 0;
  HighsInt num_lines_printed;

  double max_col_cost_relative_error = 0;
  HighsInt max_col_cost_relative_error_col = -1;
  double max_col_bound_relative_error = 0;
  HighsInt max_col_bound_relative_error_col = -1;
  double max_row_bound_relative_error = 0;
  HighsInt max_row_bound_relative_error_row = -1;
  const HighsInt small_dim = 10;

  if (dev_run) printf(" --- Testing cost ranging ---\n");
  bool small_numCol = numCol < small_dim;
  error_report_threshold = initial_error_report_threshold;
  num_lines_printed = 0;
  const bool test_all = true;
  const bool test_all_col_cost = test_all;
  const bool test_all_col_bound = test_all;
  const bool test_all_row_bound = test_all;
  HighsInt test_col_cost = min(HighsInt{0}, numCol - 1);
  HighsInt test_col_bound = min(HighsInt{0}, numCol - 1);
  HighsInt test_row_bound = min(HighsInt{0}, numRow - 1);
  HighsInt from_i;
  HighsInt to_i;
  if (test_all_col_cost) {
    from_i = 0;
    to_i = numCol;
  } else {
    from_i = test_col_cost;
    to_i = from_i + 1;
  }
  for (HighsInt i = from_i; i < to_i; i++) {
    double col_cost_up_value = ranging.col_cost_up.value_[i];
    double col_cost_up_objective = ranging.col_cost_up.objective_[i];
    double col_cost_dn_value = ranging.col_cost_dn.value_[i];
    double col_cost_dn_objective = ranging.col_cost_dn.objective_[i];
    double cost = lp.col_cost_[i];
    double solved_up = 0;
    double solved_dn = 0;
    double error;
    // Col cost up ranging
    if (col_cost_up_value < kHighsInf) {
      highs.changeColCost(i, col_cost_up_value);
      highs.setBasis(basis);
      if (test_all_col_cost) {
        quietRun(highs);
      } else {
        highs.run();
      }
      REQUIRE(modelStatusOk(highs));
      solved_up = highs.getObjectiveValue();
      highs.changeColCost(i, cost);
      error = fabs(solved_up - col_cost_up_objective);
    } else {
      solved_up = col_cost_up_objective;
      error = 0;
    }
    double relative_up_error = error / relative_error_denominator;
    if (relative_up_error >= relative_error_tolerance) {
      if (dev_run)
        printf("Col %" HIGHSINT_FORMAT
               ": %g = relative_up_error >= relative_error_tolerance = %g | "
               "%g %g %g\n",
               i, relative_up_error, relative_error_tolerance, solved_up,
               col_cost_up_objective, error);
      num_relative_error++;
      if (dev_run) REQUIRE(relative_up_error < relative_error_tolerance);
    }
    max_relative_error = max(max_relative_error, relative_up_error);
    sum_relative_error += relative_up_error;
    // Col cost down ranging
    if (col_cost_dn_value > -kHighsInf) {
      highs.changeColCost(i, col_cost_dn_value);
      highs.setBasis(basis);
      if (test_all_col_cost) {
        quietRun(highs);
      } else {
        highs.run();
      }
      REQUIRE(modelStatusOk(highs));
      solved_dn = highs.getObjectiveValue();
      highs.changeColCost(i, cost);
      error = fabs(solved_dn - col_cost_dn_objective);
    } else {
      solved_dn = col_cost_dn_objective;
      error = 0;
    }
    double relative_dn_error = error / relative_error_denominator;
    if (relative_dn_error >= relative_error_tolerance) {
      if (dev_run)
        printf("Col %" HIGHSINT_FORMAT
               ": %g = relative_dn_error >= relative_error_tolerance = %g | "
               "%g %g %g\n",
               i, relative_dn_error, relative_error_tolerance, solved_dn,
               col_cost_dn_objective, error);
      num_relative_error++;
      if (dev_run) REQUIRE(relative_up_error < relative_error_tolerance);
    }
    max_relative_error = max(max_relative_error, relative_dn_error);
    sum_relative_error += relative_dn_error;
    double relative_error = max(relative_up_error, relative_dn_error);
    if (relative_error > max_col_cost_relative_error) {
      max_col_cost_relative_error = relative_error;
      max_col_cost_relative_error_col = i;
    }
    if (small_numCol || relative_error > error_report_threshold) {
      if (num_lines_printed % 50 == 0) colCostColumnHeader();
      if (dev_run)
        printf("%3" HIGHSINT_FORMAT
               " %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
               i, cost, col_dual[i], col_cost_up_value, col_cost_up_objective,
               solved_up, relative_up_error, col_cost_dn_value,
               col_cost_dn_objective, solved_dn, relative_dn_error);
      error_report_threshold = 10 * error_report_threshold;
      num_lines_printed++;
    }
  }

  if (dev_run) printf(" --- Testing column bounds ranging ---\n");
  error_report_threshold = initial_error_report_threshold;
  num_lines_printed = 0;
  if (test_all_col_bound) {
    from_i = 0;
    to_i = numCol;
  } else {
    from_i = test_col_bound;
    to_i = from_i + 1;
  }
  for (HighsInt i = from_i; i < to_i; i++) {
    double col_bound_up_value = ranging.col_bound_up.value_[i];
    double col_bound_up_objective = ranging.col_bound_up.objective_[i];
    double col_bound_dn_value = ranging.col_bound_dn.value_[i];
    double col_bound_dn_objective = ranging.col_bound_dn.objective_[i];
    double lower = lp.col_lower_[i];
    double upper = lp.col_upper_[i];
    double new_lower;
    double new_upper;
    double solved_up = 0;
    double solved_dn = 0;
    double error;
    // Col bound up ranging
    if (col_bound_up_value < kHighsInf) {
      // Free cols should not have a finite col_bound_up_value
      assert(col_status[i] != HighsBasisStatus::kZero);
      new_lower = lower;
      new_upper = upper;
      if (col_status[i] != HighsBasisStatus::kBasic) {
        // Nonbasic
        if (lower == upper) {
          new_lower = col_bound_up_value;
          new_upper = col_bound_up_value;
        } else if (col_status[i] == HighsBasisStatus::kLower) {
          new_lower = col_bound_up_value;
        } else {
          new_upper = col_bound_up_value;
        }
      } else {
        new_lower = col_bound_up_value;
      }
      assessNewBounds(new_lower, new_upper);
      highs.changeColBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      if (test_all_col_bound) {
        quietRun(highs);
      } else {
        highs.run();
      }
      REQUIRE(modelStatusOk(highs));
      solved_up = highs.getObjectiveValue();
      highs.changeColBounds(i, lower, upper);
      error = fabs(solved_up - col_bound_up_objective);
    } else {
      solved_up = col_bound_up_objective;
      error = 0;
    }
    double relative_up_error = error / relative_error_denominator;
    if (relative_up_error >= relative_error_tolerance) {
      if (dev_run)
        printf("Col %" HIGHSINT_FORMAT
               ": %g = relative_up_error >= relative_error_tolerance = %g | "
               "%g %g %g\n",
               i, relative_up_error, relative_error_tolerance, solved_up,
               col_bound_up_objective, error);
      num_relative_error++;
      if (dev_run) REQUIRE(relative_up_error < relative_error_tolerance);
    }
    max_relative_error = max(max_relative_error, relative_up_error);
    sum_relative_error += relative_up_error;
    // Col bound down ranging
    if (col_bound_dn_value > -kHighsInf) {
      // Free cols should not have a finite col_bound_dn_value
      assert(col_status[i] != HighsBasisStatus::kZero);
      new_lower = lower;
      new_upper = upper;
      if (col_status[i] != HighsBasisStatus::kBasic) {
        // Nonbasic
        if (lower == upper) {
          new_lower = col_bound_dn_value;
          new_upper = col_bound_dn_value;
        } else if (col_status[i] == HighsBasisStatus::kLower) {
          new_lower = col_bound_dn_value;
        } else {
          new_upper = col_bound_dn_value;
        }
      } else {
        new_upper = col_bound_dn_value;
      }
      assessNewBounds(new_lower, new_upper);
      highs.changeColBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      if (test_all_col_bound) {
        quietRun(highs);
      } else {
        highs.run();
      }
      REQUIRE(modelStatusOk(highs));
      solved_dn = highs.getObjectiveValue();
      highs.changeColBounds(i, lower, upper);
      error = fabs(solved_dn - col_bound_dn_objective);
    } else {
      solved_dn = col_bound_dn_objective;
      error = 0;
    }
    double relative_dn_error = error / relative_error_denominator;
    if (relative_dn_error >= relative_error_tolerance) {
      if (dev_run)
        printf("Col %" HIGHSINT_FORMAT
               ": %g = relative_dn_error >= relative_error_tolerance = %g | "
               "%g %g %g\n",
               i, relative_dn_error, relative_error_tolerance, solved_dn,
               col_bound_dn_objective, error);
      num_relative_error++;
      if (dev_run) REQUIRE(relative_up_error < relative_error_tolerance);
    }
    max_relative_error = max(max_relative_error, relative_dn_error);
    sum_relative_error += relative_dn_error;
    double relative_error = max(relative_up_error, relative_dn_error);
    if (relative_error > max_col_bound_relative_error) {
      max_col_bound_relative_error = relative_error;
      max_col_bound_relative_error_col = i;
    }
    if (small_numCol || relative_error > error_report_threshold) {
      if (num_lines_printed % 50 == 0) colBoundcolumnHeader();
      if (dev_run)
        printf("%3" HIGHSINT_FORMAT
               " %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
               i, lower, upper, col_value[i], col_dual[i], col_bound_up_value,
               col_bound_up_objective, solved_up, relative_up_error,
               col_bound_dn_value, col_bound_dn_objective, solved_dn,
               relative_dn_error);
      error_report_threshold = 10 * error_report_threshold;
      num_lines_printed++;
    }
  }
  if (dev_run) printf(" --- Testing row bounds ranging ---\n");
  bool small_numRow = numRow < small_dim;
  error_report_threshold = initial_error_report_threshold;
  num_lines_printed = 0;
  if (test_all_row_bound) {
    from_i = 0;
    to_i = numRow;
  } else {
    from_i = test_row_bound;
    to_i = from_i + 1;
  }
  for (HighsInt i = from_i; i < to_i; i++) {
    double row_bound_up_value = ranging.row_bound_up.value_[i];
    double row_bound_up_objective = ranging.row_bound_up.objective_[i];
    double row_bound_dn_value = ranging.row_bound_dn.value_[i];
    double row_bound_dn_objective = ranging.row_bound_dn.objective_[i];
    double lower = lp.row_lower_[i];
    double upper = lp.row_upper_[i];
    double new_lower = lower;
    double new_upper = upper;
    double solved_up = 0;
    double solved_dn = 0;
    double error;
    // Row bound up ranging
    if (row_bound_up_value < kHighsInf) {
      // Free rows should not have a finite row_bound_up_value
      assert(row_status[i] != HighsBasisStatus::kZero);
      new_lower = lower;
      new_upper = upper;
      if (row_status[i] != HighsBasisStatus::kBasic) {
        // Nonbasic
        if (lower == upper) {
          new_lower = row_bound_up_value;
          new_upper = row_bound_up_value;
        } else if (row_status[i] == HighsBasisStatus::kLower) {
          new_lower = row_bound_up_value;
        } else {
          new_upper = row_bound_up_value;
        }
      } else {
        new_lower = row_bound_up_value;
      }
      assessNewBounds(new_lower, new_upper);
      highs.changeRowBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      if (test_all_row_bound) {
        quietRun(highs);
      } else {
        highs.run();
      }
      REQUIRE(modelStatusOk(highs));
      solved_up = highs.getObjectiveValue();
      highs.changeRowBounds(i, lower, upper);
      error = fabs(solved_up - row_bound_up_objective);
    } else {
      solved_up = row_bound_up_objective;
      error = 0;
    }
    double relative_up_error = error / relative_error_denominator;
    if (relative_up_error >= relative_error_tolerance) {
      if (dev_run)
        printf("Row %" HIGHSINT_FORMAT
               ": %g = relative_up_error >= relative_error_tolerance = %g | "
               "%g %g %g\n",
               i, relative_up_error, relative_error_tolerance, solved_up,
               row_bound_up_objective, error);
      num_relative_error++;
      if (dev_run) REQUIRE(relative_up_error < relative_error_tolerance);
    }
    max_relative_error = max(max_relative_error, relative_up_error);
    sum_relative_error += relative_up_error;
    // Row bound down ranging
    if (row_bound_dn_value > -kHighsInf) {
      // Free rows should not have a finite row_bound_dn_value
      assert(row_status[i] != HighsBasisStatus::kZero);
      new_lower = lower;
      new_upper = upper;
      if (row_status[i] != HighsBasisStatus::kBasic) {
        // Nonbasic
        if (lower == upper) {
          new_lower = row_bound_dn_value;
          new_upper = row_bound_dn_value;
        } else if (row_status[i] == HighsBasisStatus::kLower) {
          new_lower = row_bound_dn_value;
        } else {
          new_upper = row_bound_dn_value;
        }
      } else {
        new_upper = row_bound_dn_value;
      }
      assessNewBounds(new_lower, new_upper);
      highs.changeRowBounds(i, new_lower, new_upper);
      highs.setBasis(basis);
      if (test_all_row_bound) {
        quietRun(highs);
      } else {
        highs.run();
      }
      REQUIRE(modelStatusOk(highs));
      solved_dn = highs.getObjectiveValue();
      highs.changeRowBounds(i, lower, upper);
      error = fabs(solved_dn - row_bound_dn_objective);
    } else {
      solved_dn = row_bound_dn_objective;
      error = 0;
    }
    double relative_dn_error = error / relative_error_denominator;
    if (relative_dn_error >= relative_error_tolerance) {
      if (dev_run)
        printf("Row %" HIGHSINT_FORMAT
               ": %g = relative_dn_error >= relative_error_tolerance = %g | "
               "%g %g %g\n",
               i, relative_dn_error, relative_error_tolerance, solved_dn,
               row_bound_dn_objective, error);
      num_relative_error++;
      if (dev_run) REQUIRE(relative_up_error < relative_error_tolerance);
    }
    max_relative_error = max(max_relative_error, relative_dn_error);
    sum_relative_error += relative_dn_error;
    double relative_error = max(relative_up_error, relative_dn_error);
    if (relative_error > max_row_bound_relative_error) {
      max_row_bound_relative_error = relative_error;
      max_row_bound_relative_error_row = i;
    }
    if (small_numRow || relative_error > error_report_threshold) {
      if (num_lines_printed % 50 == 0) rowBoundColumnHeader();
      if (dev_run)
        printf("%3" HIGHSINT_FORMAT
               " %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g %12g\n",
               i, lower, upper, row_value[i], row_dual[i], row_bound_up_value,
               row_bound_up_objective, solved_up, relative_up_error,
               row_bound_dn_value, row_bound_dn_objective, solved_dn,
               relative_dn_error);
      error_report_threshold = 10 * error_report_threshold;
      num_lines_printed++;
    }
  }
  if (dev_run) {
    if (max_col_cost_relative_error_col >= 0)
      printf(
          "Max col cost  relative objective error = %g in col %" HIGHSINT_FORMAT
          "\n",
          max_col_cost_relative_error, max_col_cost_relative_error_col);
    if (max_col_bound_relative_error_col >= 0)
      printf(
          "Max col bound relative objective error = %g in col %" HIGHSINT_FORMAT
          "\n",
          max_col_bound_relative_error, max_col_bound_relative_error_col);
    if (max_row_bound_relative_error_row >= 0)
      printf(
          "Max row bound relative objective error = %g in row %" HIGHSINT_FORMAT
          "\n",
          max_row_bound_relative_error, max_row_bound_relative_error_row);
    printf("Num / max / sum relative objective error = %" HIGHSINT_FORMAT
           " / %g / %g\n",
           num_relative_error, max_relative_error, sum_relative_error);
  }
  REQUIRE(num_relative_error < 10);
  REQUIRE(max_relative_error < relative_error_tolerance);
}
