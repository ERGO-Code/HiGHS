#include "HCheckConfig.h"
#include "Highs.h"
#include "SpecialLps.h"
#include "catch.hpp"
#include "lp_data/HConst.h"

const bool dev_run = true;
const double zero_ray_value_tolerance = 1e-14;

void checkRayDirection(const HighsInt dim, const vector<double>& ray_value,
                       const vector<double>& expected_ray_value) {
  bool ray_error = false;
  HighsInt from_ix = -1;
  for (HighsInt ix = 0; ix < dim; ix++) {
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
  for (HighsInt ix = from_ix + 1; ix < dim; ix++) {
    double scaled_expected_ray_value = expected_ray_value[ix] * scale;
    if (fabs(ray_value[ix] - scaled_expected_ray_value) >
        zero_ray_value_tolerance) {
      ray_error = true;
      break;
    }
  }
  REQUIRE(!ray_error);
}

void checkDualUnboundednessDirection(Highs& highs, const vector<double>& dual_unboundedness_direction_value) {
  const HighsLp& lp = highs.getLp();
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  double dual_unboundedness_direction_error_norm = 0;
  const vector<double>& colLower = lp.col_lower_;
  const vector<double>& colUpper = lp.col_upper_;
  //  const vector<double>& rowLower = lp.row_lower_;
  //  const vector<double>& rowUpper = lp.row_upper_;
  const vector<HighsBasisStatus>& col_status = highs.getBasis().col_status;
  //  const vector<HighsBasisStatus>& row_status = highs.getBasis().row_status;
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    // Nothing to check if basic or fixed (so value can be anything)
    if (col_status[iCol] == HighsBasisStatus::kBasic ||
        colLower[iCol] == colUpper[iCol])
      continue;
    if (col_status[iCol] == HighsBasisStatus::kLower) {
      // At lower bound so value should be non-positive
      if (dual_unboundedness_direction_value[iCol] > 0) {
        dual_unboundedness_direction_error_norm += fabs(dual_unboundedness_direction_value[iCol]);
        if (dual_unboundedness_direction_value[iCol] > zero_ray_value_tolerance && dev_run)
          printf("Col %3" HIGHSINT_FORMAT
                 " is at lower bound so dual step should be "
                 "non-positive, and is %g\n",
                 iCol, dual_unboundedness_direction_value[iCol]);
      }
    } else if (col_status[iCol] == HighsBasisStatus::kUpper) {
      // At upper bound so value should be non-negative
      if (dual_unboundedness_direction_value[iCol] < 0) {
        dual_unboundedness_direction_error_norm += fabs(dual_unboundedness_direction_value[iCol]);
        if (dual_unboundedness_direction_value[iCol] < -zero_ray_value_tolerance && dev_run)
          printf("Col %3" HIGHSINT_FORMAT
                 " is at upper bound so dual step should be "
                 "non-negative, and is %g\n",
                 iCol, dual_unboundedness_direction_value[iCol]);
      }
    } else {
      // Free so value should be zero
      assert(col_status[iCol] == HighsBasisStatus::kZero);
      if (fabs(dual_unboundedness_direction_value[iCol]) > 0) {
        dual_unboundedness_direction_error_norm += fabs(dual_unboundedness_direction_value[iCol]);
        if (fabs(dual_unboundedness_direction_value[iCol]) > zero_ray_value_tolerance && dev_run)
          printf("Col %3" HIGHSINT_FORMAT
                 " is free so dual step should be zero, and is %g\n",
                 iCol, dual_unboundedness_direction_value[iCol]);
      }
    }
  }
  if (dev_run)
    printf("checkDualUnboundednessDirection: dual_unboundedness_direction_error_norm = %g\n", dual_unboundedness_direction_error_norm);
  REQUIRE(dual_unboundedness_direction_error_norm < 1e-6);
}

void checkDualRayValue(Highs& highs, const vector<double>& dual_ray_value) {
  const HighsLp& lp = highs.getLp();
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  double ray_error_norm = 0;
  //  const vector<double>& colLower = lp.col_lower_;
  //  const vector<double>& colUpper = lp.col_upper_;
    const vector<double>& rowLower = lp.row_lower_;
    const vector<double>& rowUpper = lp.row_upper_;
  const vector<HighsBasisStatus>& col_status = highs.getBasis().col_status;
    const vector<HighsBasisStatus>& row_status = highs.getBasis().row_status;
  vector<double> dual_unboundedness_direction_value;
  dual_unboundedness_direction_value.assign(numCol, 0.0);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    if (col_status[iCol] == HighsBasisStatus::kBasic) continue;
    // Get the tableau row entry for this nonbasic column
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
         iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
      dual_unboundedness_direction_value[iCol] +=
          dual_ray_value[lp.a_matrix_.index_[iEl]] * lp.a_matrix_.value_[iEl];
  }
  checkDualUnboundednessDirection(highs, dual_unboundedness_direction_value);
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    if (row_status[iRow] == HighsBasisStatus::kBasic ||
        rowLower[iRow] == rowUpper[iRow])
      continue;
    if (row_status[iRow] == HighsBasisStatus::kLower) {
      // At lower bound so value should be non-negative
      if (dual_ray_value[iRow] < 0) {
        ray_error_norm += fabs(dual_ray_value[iRow]);
        if (dual_ray_value[iRow] < -zero_ray_value_tolerance && dev_run)
          printf("Row %3" HIGHSINT_FORMAT
                 " is at lower bound so dual step should be "
                 "non-negative, and is %g\n",
                 iRow, dual_ray_value[iRow]);
      }
    } else if (row_status[iRow] == HighsBasisStatus::kUpper) {
      // At upper bound so value should be non-positive
      if (dual_ray_value[iRow] > 0) {
        ray_error_norm += fabs(dual_ray_value[iRow]);
        if (dual_ray_value[iRow] > zero_ray_value_tolerance && dev_run)
          printf("Row %3" HIGHSINT_FORMAT
                 " is at upper bound so dual step should be "
                 "non-positive, and is %g\n",
                 iRow, dual_ray_value[iRow]);
      }
    } else {
      // Free so value should be zero
      assert(row_status[iRow] == HighsBasisStatus::kZero);
      if (fabs(dual_ray_value[iRow]) > 0) {
        ray_error_norm += fabs(dual_ray_value[iRow]);
        if (fabs(dual_ray_value[iRow]) > zero_ray_value_tolerance && dev_run)
          printf("Row %3" HIGHSINT_FORMAT
                 " is free so dual step should be zero, and is %g\n",
                 iRow, dual_ray_value[iRow]);
      }
    }
  }
  if (dev_run)
    printf("checkDualRayValue: ray_error_norm = %g\n", ray_error_norm);
  REQUIRE(ray_error_norm < 1e-6);
}

void checkPrimalRayValue(Highs& highs, const vector<double>& primal_ray_value) {
  const HighsLp& lp = highs.getLp();
  HighsInt numCol = lp.num_col_;
  HighsInt numRow = lp.num_row_;
  double ray_error_norm = 0;
  const vector<double>& colLower = lp.col_lower_;
  const vector<double>& colUpper = lp.col_upper_;
  const vector<double>& rowLower = lp.row_lower_;
  const vector<double>& rowUpper = lp.row_upper_;
  double dual_feasibility_tolerance;
  highs.getOptionValue("dual_feasibility_tolerance",
                       dual_feasibility_tolerance);
  vector<double> row_ray_value;
  row_ray_value.assign(numRow, 0.0);
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    for (HighsInt iEl = lp.a_matrix_.start_[iCol];
         iEl < lp.a_matrix_.start_[iCol + 1]; iEl++)
      row_ray_value[lp.a_matrix_.index_[iEl]] +=
          primal_ray_value[iCol] * lp.a_matrix_.value_[iEl];
  }
  for (HighsInt iCol = 0; iCol < numCol; iCol++) {
    if (primal_ray_value[iCol] > 0) {
      // Upper bound must be infinite
      if (colUpper[iCol] < kHighsInf) {
        ray_error_norm += fabs(primal_ray_value[iCol]);
        if (primal_ray_value[iCol] > zero_ray_value_tolerance && dev_run)
          printf("Column %" HIGHSINT_FORMAT
                 " has primal ray value %g and finite upper bound of "
                 "%g\n",
                 iCol, primal_ray_value[iCol], colUpper[iCol]);
      }
    } else if (primal_ray_value[iCol] < 0) {
      // Lower bound must be infinite
      if (colLower[iCol] > -kHighsInf) {
        ray_error_norm += fabs(primal_ray_value[iCol]);
        if (primal_ray_value[iCol] < -zero_ray_value_tolerance && dev_run)
          printf("Column %" HIGHSINT_FORMAT
                 " has primal ray value %g and finite lower bound of "
                 "%g\n",
                 iCol, primal_ray_value[iCol], colLower[iCol]);
      }
    }
  }
  for (HighsInt iRow = 0; iRow < numRow; iRow++) {
    if (row_ray_value[iRow] > 0) {
      // Upper bound must be infinite
      if (rowUpper[iRow] > kHighsInf) {
        ray_error_norm += fabs(row_ray_value[iRow]);
        if (row_ray_value[iRow] > zero_ray_value_tolerance && dev_run)
          printf("Row %" HIGHSINT_FORMAT
                 " has primal ray value %g and finite upper bound of %g\n",
                 iRow, row_ray_value[iRow], rowUpper[iRow]);
      }
    } else if (row_ray_value[iRow] < -0) {
      // Lower bound must be infinite
      if (rowLower[iRow] > -kHighsInf) {
        ray_error_norm += fabs(row_ray_value[iRow]);
        if (row_ray_value[iRow] < -zero_ray_value_tolerance && dev_run)
          printf("Row %" HIGHSINT_FORMAT
                 " has primal ray value %g and finite lower bound of %g\n",
                 iRow, row_ray_value[iRow], rowLower[iRow]);
      }
    }
  }
  if (dev_run)
    printf("checkPrimalRayValue: ray_error_norm = %g\n", ray_error_norm);
  REQUIRE(ray_error_norm < 1e-6);
}

void testInfeasibleMpsLp(const std::string model,
                         const bool has_dual_ray_ = true) {
  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;
  bool has_dual_ray;
  bool has_primal_ray;
  vector<double> dual_ray_value;
  vector<double> primal_ray_value;

  Highs highs;
  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
  } else {
    highs.setOptionValue("log_dev_level", 2);
  }

  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);

  // Test dual ray for infeasible LP
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  require_model_status = HighsModelStatus::kInfeasible;
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  lp = highs.getLp();
  REQUIRE(highs.setBasis() == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == require_model_status);
  // Check that there is a dual ray
  dual_ray_value.resize(lp.num_row_);
  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::kOk);
  REQUIRE(has_dual_ray == has_dual_ray_);
  REQUIRE(highs.getDualRay(has_dual_ray, dual_ray_value.data()) ==
          HighsStatus::kOk);
  checkDualRayValue(highs, dual_ray_value);
  // Check that there is no primal ray
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::kOk);
  REQUIRE(has_primal_ray == false);
  // Now test with presolve on, and forcing ray calculation if possible
}

void testInfeasibleMpsMip(const std::string model) {
  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;
  bool has_dual_ray;
  bool has_primal_ray;

  Highs highs;
  if (!dev_run) {
    highs.setOptionValue("output_flag", false);
  } else {
    highs.setOptionValue("log_dev_level", 2);
  }

  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);

  // Test dual ray for MIP of infeasible LP
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  require_model_status = HighsModelStatus::kInfeasible;
  REQUIRE(highs.readModel(model_file) == HighsStatus::kOk);
  lp = highs.getLp();
  lp.integrality_.clear();
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    lp.integrality_.push_back(HighsVarType::kInteger);
  highs.passModel(lp);
  REQUIRE(highs.setBasis() == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == require_model_status);
  // Check that there is no dual ray
  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::kOk);
  REQUIRE(!has_dual_ray);
  // Check that there is no primal ray
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::kOk);
  REQUIRE(has_primal_ray == false);
}

void testUnboundedMpsLp(const std::string model,
                        const ObjSense sense = ObjSense::kMinimize) {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  if (dev_run) highs.setOptionValue("log_dev_level", 1);

  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;
  HighsStatus require_status;
  bool has_dual_ray;
  bool has_primal_ray;
  vector<double> dual_ray_value;
  vector<double> primal_ray_value;
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);

  // Test dual ray for unbounded LP
  model_file = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  require_model_status = HighsModelStatus::kUnbounded;
  // gas11 contains small nonzero matrix entries, so readModel yield
  // HighsStatus::kWarning
  require_status = model == "gas11" ? HighsStatus::kWarning : HighsStatus::kOk;
  REQUIRE(highs.readModel(model_file) == require_status);
  REQUIRE(highs.changeObjectiveSense(sense) == HighsStatus::kOk);
  lp = highs.getLp();
  lp.model_name_ = model;
  REQUIRE(highs.setBasis() == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);

  if (dev_run)
    printf("Solved %s with presolve: status = %s\n", lp.model_name_.c_str(),
           highs.modelStatusToString(highs.getModelStatus()).c_str());
  REQUIRE(highs.getModelStatus() == require_model_status);
  // Check that there is no dual ray
  REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::kOk);
  REQUIRE(has_dual_ray == false);
  // Check that there is a primal ray
  primal_ray_value.resize(lp.num_col_);
  REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::kOk);
  REQUIRE(has_primal_ray == true);
  REQUIRE(highs.getPrimalRay(has_primal_ray, primal_ray_value.data()) ==
          HighsStatus::kOk);
  checkPrimalRayValue(highs, primal_ray_value);
}

TEST_CASE("Rays", "[highs_test_rays]") {
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  std::string model_file;
  HighsLp lp;
  HighsModelStatus require_model_status;
  double optimal_objective;
  SpecialLps special_lps;
  bool has_dual_ray;
  bool has_primal_ray;
  vector<double> dual_ray_value;
  vector<double> dual_unboundedness_direction_value;
  vector<double> primal_ray_value;
  const bool test_scipLpi3Lp = true;
  const bool test_other = true;
  const bool test_scipLpi2Lp = test_other;
  const bool test_issue272Lp = test_other;
  const bool test_primalDualInfeasible1Lp = test_other;
  const HighsInt num_pass = 2;

  //  special_lps.issue285Lp(lp, require_model_status);
  REQUIRE(highs.setOptionValue("presolve", kHighsOffString) == HighsStatus::kOk);
  std::string presolve_status = "off";
  bool presolve_off = true;
  
  if (test_scipLpi3Lp) {
    // Test dual ray for infeasible LP
    special_lps.scipLpi3Lp(lp, require_model_status);
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    for (HighsInt k = 0; k < num_pass; k++) {
      // Loop twice, without and with presolve
      if (dev_run)
	printf("\nSolved %s with presolve %s\n", lp.model_name_.c_str(),
	       presolve_status.c_str());
      
      REQUIRE(highs.setBasis() == HighsStatus::kOk);
      REQUIRE(highs.run() == HighsStatus::kOk);
      if (dev_run) printf("Model status = %s\n", 
	       highs.modelStatusToString(highs.getModelStatus()).c_str());
   
      REQUIRE(highs.getModelStatus() == require_model_status);
      
      // Get the dual ray twice, to check that, second time, it's
      // copied from ekk_instance_.dual_ray_
      for (HighsInt p = 0; p < num_pass; p++) {
	printf("\nPass k = %1d; p = %1d\n", int(k), int(p));
	// Check that there is a dual ray
	dual_ray_value.resize(lp.num_row_);
	REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::kOk);
	// There is a dual ray iff this is the second pass or presolve is off 
	bool require_dual_ray = p == 1 || presolve_off;
	REQUIRE(has_dual_ray == require_dual_ray);

	// Get the dual ray
	REQUIRE(highs.getDualRay(has_dual_ray, dual_ray_value.data()) ==
		HighsStatus::kOk);
	vector<double> expected_dual_ray = {0.5, -1};  // From SCIP
	if (dev_run) {
	  printf("Dual ray:\nRow    computed    expected\n");
	  for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
	    printf("%3" HIGHSINT_FORMAT " %11.4g %11.4g\n", iRow,
		   dual_ray_value[iRow], expected_dual_ray[iRow]);
	}
	checkDualRayValue(highs, dual_ray_value);
      }

      // Check that there is no primal ray
      REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::kOk);
      REQUIRE(has_primal_ray == false);
      highs.clearSolver();

      // Now check dual unboundedness direction
      dual_unboundedness_direction_value.resize(lp.num_col_);
      bool has_dual_unboundedness_direction;
      REQUIRE(highs.getDualUnboundednessDirection(has_dual_unboundedness_direction,
						  dual_unboundedness_direction_value.data()) ==
	      HighsStatus::kOk);
      REQUIRE(has_dual_unboundedness_direction);

      checkDualUnboundednessDirection(highs, dual_unboundedness_direction_value);
      presolve_status = "on";
      presolve_off = false;
      REQUIRE(highs.setOptionValue("presolve", kHighsOnString) == HighsStatus::kOk);
      highs.clearSolver();
    }
  }
  
  if (test_issue272Lp) {
    // Check that there are no rays for this LP
    REQUIRE(highs.setOptionValue("presolve", kHighsOffString) == HighsStatus::kOk);
    special_lps.issue272Lp(lp, require_model_status, optimal_objective);
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.setBasis() == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    REQUIRE(highs.getModelStatus() == require_model_status);
    
    REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::kOk);
    REQUIRE(has_dual_ray == false);
    REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::kOk);
    REQUIRE(has_primal_ray == false);
  }
  
  if (test_scipLpi2Lp) {
    // Test primal ray for unbounded LP
    REQUIRE(highs.setOptionValue("presolve", kHighsOffString) == HighsStatus::kOk);
    special_lps.scipLpi2Lp(lp, require_model_status);
    vector<double> expected_dual_ray = {0.5, -1};  // From SCIP
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.setBasis() == HighsStatus::kOk);
    if (dev_run) highs.setOptionValue("log_dev_level", 1);
    REQUIRE(highs.run() == HighsStatus::kOk);
    if (dev_run)
      printf("Solved %s with presolve: status = %s\n", lp.model_name_.c_str(),
	     highs.modelStatusToString(highs.getModelStatus()).c_str());
    
    if (dev_run) highs.writeSolution("", kSolutionStylePretty);
    REQUIRE(highs.getModelStatus() == require_model_status);
    
    // Check that there is no dual ray
    REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::kOk);
    REQUIRE(has_dual_ray == false);
    
    // Check that a primal ray can be obtained
    primal_ray_value.resize(lp.num_col_);
    REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::kOk);
    REQUIRE(has_primal_ray == true);
    REQUIRE(highs.getPrimalRay(has_primal_ray, primal_ray_value.data()) ==
	    HighsStatus::kOk);
    vector<double> expected_primal_ray = {0.5, -1};
    if (dev_run) {
      printf("Primal ray:\nRow    computed    expected\n");
      for (HighsInt iRow = 0; iRow < lp.num_row_; iRow++)
	printf("%3" HIGHSINT_FORMAT " %11.4g %11.4g\n", iRow,
	       primal_ray_value[iRow], expected_primal_ray[iRow]);
    }
    checkRayDirection(lp.num_row_, dual_ray_value, expected_dual_ray);
    checkPrimalRayValue(highs, primal_ray_value);
  }
  if (test_primalDualInfeasible1Lp) {
    // Test that there's no primal or dual ray for this LP that is both
    // primal and dual infeasible
    REQUIRE(highs.setOptionValue("presolve", kHighsOffString) == HighsStatus::kOk);
    special_lps.primalDualInfeasible1Lp(lp, require_model_status);
    REQUIRE(highs.passModel(lp) == HighsStatus::kOk);
    REQUIRE(highs.setBasis() == HighsStatus::kOk);
    REQUIRE(highs.run() == HighsStatus::kOk);
    if (dev_run)
      printf("Solved %s with presolve: status = %s\n", lp.model_name_.c_str(),
	     highs.modelStatusToString(highs.getModelStatus()).c_str());
    REQUIRE(highs.getModelStatus() == require_model_status);
    // Check that there is no dual ray
    REQUIRE(highs.getDualRay(has_dual_ray) == HighsStatus::kOk);
    REQUIRE(has_dual_ray == false);
    // Check that there is no primal ray
    REQUIRE(highs.getPrimalRay(has_primal_ray) == HighsStatus::kOk);
    REQUIRE(has_primal_ray == false);
  }
}

TEST_CASE("Rays-gas11", "[highs_test_rays]") { testUnboundedMpsLp("gas11"); }
TEST_CASE("Rays-adlittlemax", "[highs_test_rays]") {
  testUnboundedMpsLp("adlittle", ObjSense::kMaximize);
}

TEST_CASE("Rays-galenet", "[highs_test_rays]") {
  testInfeasibleMpsLp("galenet");
}

TEST_CASE("Rays-woodinfe", "[highs_test_rays]") {
  testInfeasibleMpsLp("woodinfe");
}

TEST_CASE("Rays-klein1", "[highs_test_rays]") {
  testInfeasibleMpsLp("klein1");
}

TEST_CASE("Rays-gams10am", "[highs_test_rays]") {
  testInfeasibleMpsLp("gams10am");
  testInfeasibleMpsMip("gams10am");
}

TEST_CASE("Rays-ex72a", "[highs_test_rays]") { testInfeasibleMpsLp("ex72a"); }

TEST_CASE("Rays-forest6", "[highs_test_rays]") {
  testInfeasibleMpsLp("forest6");
}

TEST_CASE("Rays-box1", "[highs_test_rays]") { testInfeasibleMpsLp("box1"); }

TEST_CASE("Rays-bgetam", "[highs_test_rays]") { testInfeasibleMpsLp("bgetam"); }

TEST_CASE("Rays-464a", "[highs_test_rays]") {
  // The model is:
  //    min -x - y
  //         x - y == 0
  //
  // which has a primal ray: [d, d], for all d > 0.
  Highs highs;
  highs.setOptionValue("output_flag", dev_run);
  double inf = highs.getInfinity();
  highs.addCol(-1.0, -inf, inf, 0, NULL, NULL);
  highs.addCol(-1.0, -inf, inf, 0, NULL, NULL);
  HighsInt aindex[2] = {0, 1};
  double avalue[2] = {1.0, -1.0};
  highs.addRow(0.0, 0.0, 2, aindex, avalue);
  highs.setOptionValue("presolve", kHighsOffString);
  highs.run();
  if (dev_run)
    printf("Solved 464a without presolve: status = %s\n",
           highs.modelStatusToString(highs.getModelStatus()).c_str());
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
  bool has_ray = false;
  REQUIRE(highs.getPrimalRay(has_ray) == HighsStatus::kOk);
  REQUIRE(has_ray == true);
  vector<double> ray_value;
  ray_value.assign(2, NAN);
  highs.getPrimalRay(has_ray, ray_value.data());
  checkPrimalRayValue(highs, ray_value);
  REQUIRE(has_ray);
  REQUIRE(ray_value[0] == ray_value[1]);
  REQUIRE(ray_value[0] > 0);
  // Now repeat with presolve
  highs.setOptionValue("presolve", kHighsOnString);
  highs.clearSolver();
  highs.run();
  if (dev_run)
    printf("Solved 464a with presolve: status = %s\n",
           highs.modelStatusToString(highs.getModelStatus()).c_str());
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
  has_ray = false;
  REQUIRE(highs.getPrimalRay(has_ray) == HighsStatus::kOk);
  REQUIRE(has_ray == true);
  ray_value.assign(2, NAN);
  highs.getPrimalRay(has_ray, ray_value.data());
  checkPrimalRayValue(highs, ray_value);
  REQUIRE(has_ray);
  REQUIRE(ray_value[0] == ray_value[1]);
  REQUIRE(ray_value[0] > 0);
}

TEST_CASE("Rays-464b", "[highs_test_rays]") {
  // The model is:
  //    min -x - y
  //         x - y == 0
  //         x,  y >= 0
  //
  // which has a primal ray: [d, d], for all d > 0.
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  double inf = highs.getInfinity();
  highs.addCol(-1.0, 0.0, inf, 0, NULL, NULL);
  highs.addCol(-1.0, 0.0, inf, 0, NULL, NULL);
  HighsInt aindex[2] = {0, 1};
  double avalue[2] = {1.0, -1.0};
  highs.addRow(0.0, 0.0, 2, aindex, avalue);
  //  highs.setOptionValue("presolve", "off");
  highs.run();
  if (dev_run)
    printf("Solved 464b without presolve: status = %s\n",
           highs.modelStatusToString(highs.getModelStatus()).c_str());
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kUnbounded);
  bool has_ray = false;
  REQUIRE(highs.getPrimalRay(has_ray) == HighsStatus::kOk);
  REQUIRE(has_ray == true);
  vector<double> ray_value;
  ray_value.assign(2, NAN);
  highs.getPrimalRay(has_ray, ray_value.data());
  checkPrimalRayValue(highs, ray_value);
  REQUIRE(has_ray);
  REQUIRE(ray_value[0] == ray_value[1]);
  REQUIRE(ray_value[0] > 0);
}

/*
  TEST_CASE("Rays-infeasible-qp", "[highs_test_rays]") {
  HighsModel model;
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;
  lp.num_col_ = 2;
  lp.num_row_ = 1;
  lp.col_cost_ = {0, 0};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-1};
  lp.row_upper_ = {-1};
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, 2};
  lp.a_matrix_.index_ = {0, 1};
  lp.a_matrix_.value_ = {1, 1};
  hessian.dim_ = 2;
  hessian.start_ = {0, 1, 2};
  hessian.index_ = {0, 1};
  hessian.value_ = {1, 1};
  Highs highs;
  //highs.setOptionValue("output_flag", dev_run);
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  highs.run();
  //  if (dev_run)
  printf("Solved infeasible QP: status = %s\n",
  highs.modelStatusToString(highs.getModelStatus()).c_str());
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kInfeasible);
  bool has_ray = false;
  REQUIRE(highs.getDualRay(has_ray) == HighsStatus::kOk);
  REQUIRE(has_ray == true);
  std::vector<double> ray_value;
  ray_value.assign(2, NAN);
  highs.getDualRay(has_ray, ray_value.data());
  }
*/
