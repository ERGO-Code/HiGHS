#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = false;
const double zero_ray_value_tolerance = 1e-8;

void checkRayDirection(const int dim, const vector<double>& ray_value,
                       const vector<double>& expected_ray_value) {
  bool ray_error = false;
  int from_ix = -1;
  for (int ix = 0; ix < dim; ix++) {
    if (fabs(expected_ray_value[ix]) > zero_ray_value_tolerance) {
      // Found a nonzero in the expected ray values
      from_ix = ix;
      break;
    } else {
      // Found a zero in the expected ray values, so make sure that
      // the ray value is also zero
      if (fabs(ray_value[ix]) > zero_ray_value_tolerance) {
        ray_error = true;
        break;
      }
    }
  }
  REQUIRE(!ray_error);
  if (from_ix < 0) return;
  double scale = ray_value[from_ix] / expected_ray_value[from_ix];
  for (int ix = from_ix + 1; ix < dim; ix++) {
    double scaled_expected_ray_value = expected_ray_value[ix] * scale;
    if (fabs(ray_value[ix] - scaled_expected_ray_value) >
        zero_ray_value_tolerance) {
      ray_error = true;
      break;
    }
  }
  REQUIRE(!ray_error);
}

void checkDualRayValue(Highs& highs, const vector<double>& dual_ray_value) {
  const HighsLp& lp = highs.getLp();
  int numCol = lp.numCol_;
  int numRow = lp.numRow_;
  bool ray_error = false;
  const vector<double>& colLower = lp.colLower_;
  const vector<double>& colUpper = lp.colUpper_;
  const vector<double>& rowLower = lp.rowLower_;
  const vector<double>& rowUpper = lp.rowUpper_;
  const vector<HighsBasisStatus>& col_status = highs.getBasis().col_status;
  const vector<HighsBasisStatus>& row_status = highs.getBasis().row_status;
  vector<double> tableau_row;
  tableau_row.assign(numCol, 0.0);
  for (int iCol = 0; iCol < numCol; iCol++) {
    if (col_status[iCol] == HighsBasisStatus::BASIC) continue;
    // Get the tableau row entry for this nonbasic column
    for (int iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1]; iEl++)
      tableau_row[iCol] += dual_ray_value[lp.Aindex_[iEl]] * lp.Avalue_[iEl];
  }

  for (int iCol = 0; iCol < numCol; iCol++) {
    // Nothing to check if basic or fixed (so value can be anything)
    if (col_status[iCol] == HighsBasisStatus::BASIC ||
        colLower[iCol] == colUpper[iCol])
      continue;
    if (col_status[iCol] == HighsBasisStatus::LOWER) {
      // At lower bound so value should be non-positive
      if (tableau_row[iCol] > zero_ray_value_tolerance) {
        if (dev_run)
          printf(
              "Col %3d is at lower bound so dual step should be "
              "non-positive, and is %g\n",
              iCol, tableau_row[iCol]);
        ray_error = true;
      }
    } else if (col_status[iCol] == HighsBasisStatus::UPPER) {
      // At upper bound so value should be non-negative
      if (tableau_row[iCol] < -zero_ray_value_tolerance) {
        if (dev_run)
          printf(
              "Col %3d is at upper bound so dual step should be "
              "non-negative, and is %g\n",
              iCol, tableau_row[iCol]);
        ray_error = true;
      }
    } else {
      // Free so value should be zero
      assert(col_status[iCol] == HighsBasisStatus::ZERO);
      if (fabs(tableau_row[iCol]) > zero_ray_value_tolerance) {
        if (dev_run)
          printf("Col %3d is free so dual step should be zero, and is %g\n",
                 iCol, tableau_row[iCol]);
        ray_error = true;
      }
    }
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    if (row_status[iRow] == HighsBasisStatus::BASIC ||
        rowLower[iRow] == rowUpper[iRow])
      continue;
    if (row_status[iRow] == HighsBasisStatus::LOWER) {
      // At lower bound so value should be non-negative
      if (dual_ray_value[iRow] < -zero_ray_value_tolerance) {
        if (dev_run)
          printf(
              "Row %3d is at lower bound so dual step should be "
              "non-negative, and is %g\n",
              iRow, dual_ray_value[iRow]);
        ray_error = true;
      }
    } else if (row_status[iRow] == HighsBasisStatus::UPPER) {
      // At upper bound so value should be non-positive
      if (dual_ray_value[iRow] > zero_ray_value_tolerance) {
        if (dev_run)
          printf(
              "Row %3d is at upper bound so dual step should be "
              "non-positive, and is %g\n",
              iRow, dual_ray_value[iRow]);
        ray_error = true;
      }
    } else {
      // Free so value should be zero
      assert(row_status[iRow] == HighsBasisStatus::ZERO);
      if (fabs(dual_ray_value[iRow]) > zero_ray_value_tolerance) {
        if (dev_run)
          printf("Row %3d is free so dual step should be zero, and is %g\n",
                 iRow, dual_ray_value[iRow]);
        ray_error = true;
      }
    }
  }
  REQUIRE(!ray_error);
}

void checkPrimalRayValue(Highs& highs, const vector<double>& primal_ray_value) {
  const HighsLp& lp = highs.getLp();
  int numCol = lp.numCol_;
  int numRow = lp.numRow_;
  bool ray_error = false;
  const vector<double>& colLower = lp.colLower_;
  const vector<double>& colUpper = lp.colUpper_;
  const vector<double>& rowLower = lp.rowLower_;
  const vector<double>& rowUpper = lp.rowUpper_;
  double dual_feasibility_tolerance;
  highs.getHighsOptionValue("dual_feasibility_tolerance",
                            dual_feasibility_tolerance);
  vector<double> row_ray_value;
  row_ray_value.assign(numRow, 0.0);
  for (int iCol = 0; iCol < numCol; iCol++) {
    for (int iEl = lp.Astart_[iCol]; iEl < lp.Astart_[iCol + 1]; iEl++)
      row_ray_value[lp.Aindex_[iEl]] +=
          primal_ray_value[iCol] * lp.Avalue_[iEl];
  }

  for (int iCol = 0; iCol < numCol; iCol++) {
    if (primal_ray_value[iCol] > 0) {
      // Upper bound must be infinite
      if (colUpper[iCol] < HIGHS_CONST_INF) {
        ray_error = true;
        if (dev_run)
          printf(
              "Column %d has primal ray value %g and finite upper bound of "
              "%g\n",
              iCol, primal_ray_value[iCol], colUpper[iCol]);
      }
    } else if (primal_ray_value[iCol] < 0) {
      // Lower bound must be infinite
      if (colLower[iCol] > -HIGHS_CONST_INF) {
        ray_error = true;
        if (dev_run)
          printf(
              "Column %d has primal ray value %g and finite lower bound of "
              "%g\n",
              iCol, primal_ray_value[iCol], colLower[iCol]);
      }
    }
  }
  for (int iRow = 0; iRow < numRow; iRow++) {
    if (row_ray_value[iRow] > 0) {
      // Upper bound must be infinite
      if (rowUpper[iRow] > HIGHS_CONST_INF) {
        ray_error = true;
        if (dev_run)
          printf(
              "Row %d has primal ray value %g and finite upper bound of %g\n",
              iRow, row_ray_value[iRow], rowUpper[iRow]);
      }
    } else if (row_ray_value[iRow] < 0) {
      // Lower bound must be infinite
      if (rowLower[iRow] > -HIGHS_CONST_INF) {
        ray_error = true;
        if (dev_run)
          printf(
              "Row %d has primal ray value %g and finite lower bound of %g\n",
              iRow, row_ray_value[iRow], rowLower[iRow]);
      }
    }
  }
  REQUIRE(!ray_error);
}

void testInfeasibleMps(const std::string model) {
  Highs highs;
  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;
  bool has_dual_ray;
  bool has_primal_ray;
  vector<double> dual_ray_value;
  vector<double> primal_ray_value;

  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }

  REQUIRE(highs.setHighsOptionValue("presolve", "off") == HighsStatus::OK);

  // Test dual ray for unbounded LP
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  require_model_status = HighsModelStatus::PRIMAL_INFEASIBLE;
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);
  lp = highs.getLp();
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);
  // Check that there is a dual ray
  dual_ray_value.resize(lp.numRow_);
  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::OK);
  REQUIRE(has_dual_ray == true);
  REQUIRE(highs.getDualRay(has_dual_ray, &dual_ray_value[0]) ==
          HighsStatus::OK);
  checkDualRayValue(highs, dual_ray_value);
  // Check that there is no primal ray
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::OK);
  REQUIRE(has_primal_ray == false);
}

void testUnboundedMps(const std::string model,
                      const ObjSense sense = ObjSense::MINIMIZE) {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }

  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;
  bool has_dual_ray;
  bool has_primal_ray;
  vector<double> dual_ray_value;
  vector<double> primal_ray_value;
  REQUIRE(highs.setHighsOptionValue("presolve", "off") == HighsStatus::OK);

  // Test dual ray for unbounded LP
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  require_model_status = HighsModelStatus::PRIMAL_UNBOUNDED;
  REQUIRE(highs.readModel(model_file) == HighsStatus::OK);
  REQUIRE(highs.changeObjectiveSense(sense));
  lp = highs.getLp();
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);
  // Check that there is no dual ray
  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::OK);
  REQUIRE(has_dual_ray == false);
  // Check that there is a primal ray
  primal_ray_value.resize(lp.numCol_);
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::OK);
  REQUIRE(has_primal_ray == true);
  REQUIRE(highs.getPrimalRay(has_primal_ray, &primal_ray_value[0]) ==
          HighsStatus::OK);
  checkPrimalRayValue(highs, primal_ray_value);
}

TEST_CASE("Rays", "[highs_test_rays]") {
  Highs highs;
  if (!dev_run) {
    highs.setHighsLogfile();
    highs.setHighsOutput();
  }
  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  SpecialLps special_lps;
  bool has_dual_ray;
  bool has_primal_ray;
  vector<double> dual_ray_value;
  vector<double> primal_ray_value;

  //  special_lps.issue285Lp(lp, require_model_status);
  REQUIRE(highs.setHighsOptionValue("presolve", "off") == HighsStatus::OK);

  // Test dual ray for infeasible LP
  special_lps.scipLpi3Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  // Check that there is a dual ray
  dual_ray_value.resize(lp.numRow_);
  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::OK);
  REQUIRE(has_dual_ray == true);
  // Get the dual ray
  REQUIRE(highs.getDualRay(has_dual_ray, &dual_ray_value[0]) ==
          HighsStatus::OK);
  vector<double> expected_dual_ray = {0.5, -1};  // From SCIP
  if (dev_run) {
    printf("Dual ray:\nRow    computed    expected\n");
    for (int iRow = 0; iRow < lp.numRow_; iRow++)
      printf("%3d %11.4g %11.4g\n", iRow, dual_ray_value[iRow],
             expected_dual_ray[iRow]);
  }
  checkDualRayValue(highs, dual_ray_value);
  // Check that there is no primal ray
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::OK);
  REQUIRE(has_primal_ray == false);

  // Check that there are no rays for this LP
  special_lps.issue272Lp(lp, require_model_status, optimal_objective);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);

  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::OK);
  REQUIRE(has_dual_ray == false);
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::OK);
  REQUIRE(has_primal_ray == false);

  // Test primal ray for unbounded LP
  special_lps.scipLpi2Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);
  if (dev_run) highs.writeSolution("", true);

  // Check that there is no dual ray
  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::OK);
  REQUIRE(has_dual_ray == false);

  // Check that there is a primal ray
  primal_ray_value.resize(lp.numCol_);
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::OK);
  REQUIRE(has_primal_ray == true);
  REQUIRE(highs.getPrimalRay(has_primal_ray, &primal_ray_value[0]) ==
          HighsStatus::OK);
  vector<double> expected_primal_ray = {0.5, -1};
  if (dev_run) {
    printf("Primal ray:\nRow    computed    expected\n");
    for (int iRow = 0; iRow < lp.numRow_; iRow++)
      printf("%3d %11.4g %11.4g\n", iRow, primal_ray_value[iRow],
             expected_primal_ray[iRow]);
  }
  checkRayDirection(lp.numRow_, dual_ray_value, expected_dual_ray);
  checkPrimalRayValue(highs, primal_ray_value);

  // Test that there's no primal or dual ray for this LP that is both
  // primal and dual infeasible
  special_lps.primalDualInfeasible1Lp(lp, require_model_status);
  REQUIRE(highs.passModel(lp) == HighsStatus::OK);
  REQUIRE(highs.setBasis() == HighsStatus::OK);
  REQUIRE(highs.run() == HighsStatus::OK);
  REQUIRE(highs.getModelStatus() == require_model_status);
}

TEST_CASE("Rays-gas11", "[highs_test_rays]") { testUnboundedMps("gas11"); }
TEST_CASE("Rays-adlittlemax", "[highs_test_rays]") {
  testUnboundedMps("adlittle", ObjSense::MAXIMIZE);
}

TEST_CASE("Rays-galenet", "[highs_test_rays]") { testInfeasibleMps("galenet"); }

TEST_CASE("Rays-woodinfe", "[highs_test_rays]") {
  testInfeasibleMps("woodinfe");
}

TEST_CASE("Rays-klein1", "[highs_test_rays]") { testInfeasibleMps("klein1"); }

TEST_CASE("Rays-gams10am", "[highs_test_rays]") {
  testInfeasibleMps("gams10am");
}

TEST_CASE("Rays-ex72a", "[highs_test_rays]") { testInfeasibleMps("ex72a"); }

TEST_CASE("Rays-forest6", "[highs_test_rays]") { testInfeasibleMps("forest6"); }

TEST_CASE("Rays-box1", "[highs_test_rays]") { testInfeasibleMps("box1"); }

TEST_CASE("Rays-bgetam", "[highs_test_rays]") { testInfeasibleMps("bgetam"); }
