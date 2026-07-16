#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
// #include "io/FilereaderLp.h"

const bool dev_run = true;  // false;
const double inf = kHighsInf;
const double double_equal_tolerance = 1e-5;

HighsHessian getHessianDiagonal4();
HighsHessian getHessianDiagonalWithZero4();
HighsHessian getHessian4();
HighsHessian getHessian5();
HighsModel getQpQjh();
bool valuesRelEqual(const double v0, const double v1);
bool vectorsRelEqual(const HighsInt dim, const double* v0, const double* v1);
void addVars(HighsModel& model);
void testUnconConOracleSolve(const HighsHessian& model);
void testOracleSolve(const HighsModel& model);

// On entry:
//
// The values of x are in x_value
//
// If x_index is a null pointer, then it is assumed that the
// values of x are scattered in x_value
//
// If x_index is not a null pointer, it is assumed that there
// are x_num_entries values of x, packed in x_value, with the
// corresponding indices in x_index.
//
// If q_x_index is not a null pointer, and q_x_num_entries is
// non-negative, then it is assumed that only the q_x_num_entries
// indices in q_x_index of the result are needed. Typical use
// case: getting an individual Hessian entry - particularly the
// diagonal
//
// On exit:
//
// If q_x_index is a null pointer, then it is assumed that the
// values of Qx are scattered in q_x_value. Typical use case:
// forming the full vector Qx
//
// If q_x_index is not a null pointer, and q_x_num_entries was
// non-negative on entry, then the values of Qx corresponding to
// q_x_index are packed in q_x_value. Typical use case: getting
// an individual Hessian entry - particularly the diagonal
//
// If q_x_index is not a null pointer, and q_x_num_entries was
// negative on entry, then the values of Qx are packed in
// q_x_value, with corresponding indices in q_x_index. Typical
// use case: getting a column of the Hessian

HighsHessianFunctionType oracleCallSquareHessian =
  [](const HighsInt call_type,
     const HighsInt x_num_entries, const HighsInt* x_index, const double* x_value, 
     HighsInt& q_x_num_entries, HighsInt* q_x_index, double* q_x_value,
     void* hessian_p) {
      assert(kHessianOracleCallTypeMin <= call_type && call_type <= kHessianOracleCallTypeMax);
      assert(x_value != nullptr);
      assert(q_x_value != nullptr);

      // Lambda for zeroing q_x_value
      auto zeroQx = [&](const HighsInt dim) {
        for (HighsInt iCol = 0; iCol < dim; iCol++) {
	  assert(q_x_value[iCol] == 0);
	  q_x_value[iCol] = 0;
	}
      };

      HighsHessian hessian = *(static_cast<HighsHessian*>(hessian_p));
      assert(hessian.format_ == HessianFormat::kSquare);

      // Lambda for adding multiple of Hessian column into q_x_value
      auto addScaledQcol = [&](const HighsInt iCol, const double x_value) {
        for (HighsInt iEl = hessian.start_[iCol];
             iEl < hessian.start_[iCol + 1]; iEl++) {
          HighsInt iRow = hessian.index_[iEl];
          q_x_value[iRow] += hessian.value_[iEl] * x_value;
        }
      };

      if (call_type == kHessianOracleCallTypeEntry) {
	HighsInt iCol = x_index[0];
	HighsInt iRow = q_x_index[0];
	// Zero Qx value in case the Hessian entry requested is zero
	q_x_value[0] = 0;
	for (HighsInt iEl = hessian.start_[iCol];
	     iEl < hessian.start_[iCol + 1]; iEl++) {
	  if (hessian.index_[iEl] == iRow) {
	    q_x_value[0] = hessian.value_[iEl];
	    return;
	  }
	}
	// Hessian entry is zero
	return;
      } else if (call_type == kHessianOracleCallTypeColumn) {
	// Get the entries in column iCol
	q_x_num_entries = 0;
	HighsInt iCol = x_index[0];
	for (HighsInt iEl = hessian.start_[iCol];
	     iEl < hessian.start_[iCol + 1]; iEl++) {
	  q_x_index[q_x_num_entries] = hessian.index_[iEl];
	  q_x_value[q_x_num_entries] = hessian.value_[iEl] * x_value[0];
	  q_x_num_entries++;
	}
	return;
      } else {
	assert(call_type == kHessianOracleCallTypeProduct);
	assert(x_index == nullptr || x_num_entries >= 1);
	assert(q_x_index == nullptr);
	if (x_index == nullptr) {
	  // Simple product with full vector x, full vector q_x
	  zeroQx(hessian.dim_);
	  for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++)
	    addScaledQcol(iCol, x_value[iCol]);
	  return;
	} else if (x_num_entries > 1) {
	  // x is sparse with x_num_entries entries in rows x_index
	  zeroQx(hessian.dim_);
	  for (HighsInt iX = 0; iX < x_num_entries; iX++)
	    addScaledQcol(x_index[iX], x_value[iX]);
	  return;
	} else if (x_num_entries == 1) {
	  // x is sparse with one entry in row x_index
	  q_x_num_entries = 0;
	  // Get the entries in column iCol
	  HighsInt iCol = x_index[0];
	  for (HighsInt iEl = hessian.start_[iCol];
	       iEl < hessian.start_[iCol + 1]; iEl++) {
	    q_x_value[hessian.index_[iEl]] = hessian.value_[iEl] * x_value[0];
	    q_x_num_entries++;
	  }
	  return;
	}
      }

      if (x_index == nullptr) {
        // Simple product with full vector x, full vector q_x, and no
        // Qx indices required
        assert(q_x_index == nullptr);
        zeroQx(hessian.dim_);
        for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++)
          addScaledQcol(iCol, x_value[iCol]);
        return;
      } else if (x_num_entries > 1) {
        // x is sparse with x_num_entries entries in rows x_index, and
        // no Qx indices required
        assert(q_x_index == nullptr);
        zeroQx(hessian.dim_);
        for (HighsInt iX = 0; iX < x_num_entries; iX++)
          addScaledQcol(x_index[iX], x_value[iX]);
        return;
      } else if (x_num_entries == 1) {
        if (q_x_index == nullptr) {
	  // x is sparse with one entry in row x_index, and no Qx
	  // index required
	  q_x_num_entries = 0;
	  // Get the entries in column iCol
	  HighsInt iCol = x_index[0];
	  for (HighsInt iEl = hessian.start_[iCol];
	       iEl < hessian.start_[iCol + 1]; iEl++) {
	    q_x_value[hessian.index_[iEl]] = hessian.value_[iEl] * x_value[0];
	    q_x_num_entries++;
	  }
	  return;
	} else if (q_x_num_entries < 0) {
	  // x is sparse with one entry in row x_index, and all Qx index
	  // required
	  q_x_num_entries = 0;
	  // Get the entries in column iCol
	  HighsInt iCol = x_index[0];
	  for (HighsInt iEl = hessian.start_[iCol];
	       iEl < hessian.start_[iCol + 1]; iEl++) {
	    q_x_index[q_x_num_entries] = hessian.index_[iEl];
	    q_x_value[q_x_num_entries] = hessian.value_[iEl] * x_value[0];
	    q_x_num_entries++;
	  }
	  return;
	} else if (q_x_num_entries == 1) {
	  // x is sparse with one entry in row x_index, and one Qx index
	  // required
	  HighsInt iCol = x_index[0];
	  HighsInt iRow = q_x_index[0];
	  // Zero Qx value in case the Hessian entry requested is zero
	  q_x_value[0] = 0;
	  for (HighsInt iEl = hessian.start_[iCol];
	       iEl < hessian.start_[iCol + 1]; iEl++) {
	    if (hessian.index_[iEl] == iRow) {
	      q_x_value[0] = hessian.value_[iEl] * x_value[0];
	      return;
	    }
	  }
	  // Hessian entry is zero
	  return;
	}
      }
      // Case not coded, since it may be unnecessary
      assert(1234 == 5678);
    };

TEST_CASE("hessian-oracle-check", "[qp-oracle]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);

  const bool qp4 = false;
  HighsModel model;
  model.hessian_ = qp4 ? getHessian4(): getHessian5(); //
  // Add an LP so that the Hessian can be passed
  addVars(model);
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;

  HighsInt zero_diagonal_el;
  if (!qp4) {
     HighsInt zero_diagonal_col = 2;
     zero_diagonal_el = 6;
     REQUIRE(hessian.index_[zero_diagonal_el] == zero_diagonal_col);
     double zero_diagonal_original_value = hessian.value_[zero_diagonal_el];
     REQUIRE(zero_diagonal_original_value != 0);
  }
  h.passModel(model);
  //  h.writeModel("");
  h.writeModel("qp5.mps");
  h.setOptionValue("qp_regularization_value", 0);
  h.run();
  h.writeSolution("", 1);

  if (!qp4) hessian.value_[zero_diagonal_el] = 0.0;

  HighsHessian square_hessian = hessian.toSquare();
  if (dev_run) square_hessian.print();

  // Now define the triangular instance without the explicit diagonal
  // zero or diagonal entry first in each packed column (cf col
  // 1). Can't do this from the outset, because HighsHessian::toSquare
  // assumes that the Hessian is triangular with the first entry in
  // each packed column being the diagonal entry
  hessian.start_ = {0, 4, 6, 7, 9, 10};
  hessian.index_ = {0, 1, 3, 4, 4, 1, 3, 3, 4, 4};
  hessian.value_ = {5, 1, -1, 2, 1, 4, -1, 3, -2, 2};
  if (dev_run) hessian.print();

  HighsStatus return_status = h.passModel(lp);
  REQUIRE(return_status == HighsStatus::kOk);

  // No oracle defined returns error
  REQUIRE(h.checkHessianOracle() == HighsStatus::kError);

  REQUIRE(h.passHessian(lp.num_col_, oracleCallSquareHessian, &square_hessian)  == HighsStatus::kOk);


  // Check the method to extract Hessian entries from the oracle
  //
  // This is used to form the Hessian as a matrix in
  // Highs::checkHessianOracle(), so has to be tested externally
  const HessianOracle& oracle = h.getModel().hessian_.oracle_;
  std::vector<double> column;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    column.assign(lp.num_col_, 0);
    for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol + 1];
	 iEl++)
      column[hessian.index_[iEl]] = hessian.value_[iEl];
    REQUIRE(oracle.diagonal(iCol) == column[iCol]);
    for (HighsInt iRow = iCol + 1; iRow < lp.num_col_; iRow++) {
      REQUIRE(oracle.entry(iRow, iCol) == column[iRow]);
      REQUIRE(oracle.entry(iCol, iRow) == column[iRow]);
    }
  }
  // Test the asymmetry check:
  //
  // By replacing the zero (2, 4) entry with a nonzero
  HighsInt zero_el = square_hessian.numNz();
  HighsInt zero_entry_row = 2;
  HighsInt zero_entry_col = 4;
  REQUIRE(oracle.entry(zero_entry_row, zero_entry_col) == 0);

  // Replacing the zero (2, 4) entry with a nonzero
  square_hessian.start_[zero_entry_col + 1]++;
  square_hessian.index_.push_back(zero_entry_row);
  square_hessian.value_.push_back(1.0);
  if (dev_run) square_hessian.print("Added nonzero");
  // Asymmetry yields error
  REQUIRE(h.checkHessianOracle(true) == HighsStatus::kError);

  // Remove the nonzero from the (2, 4) entry
  square_hessian.start_[zero_entry_col + 1]--;
  square_hessian.index_.resize(zero_el);
  square_hessian.value_.resize(zero_el);
  if (dev_run) square_hessian.print("Reversion 0");

  // By perturbing the (0, 3) entry by less than the asymmetry
  // tolerance
  HighsInt nonzero_entry_row = 0;
  HighsInt nonzero_entry_col = 3;
  HighsInt nonzero_el = -1;
  for (HighsInt iEl = square_hessian.start_[nonzero_entry_col];
       iEl < square_hessian.start_[nonzero_entry_col + 1]; iEl++) {
    if (square_hessian.index_[iEl] == nonzero_entry_row) {
      nonzero_el = iEl;
      break;
    }
  }
  REQUIRE(nonzero_el >= 0);
  double nonzero_value = oracle.entry(nonzero_entry_row, nonzero_entry_col);
  REQUIRE(nonzero_value != 0);
  REQUIRE(square_hessian.index_[nonzero_el] == nonzero_entry_row);
  REQUIRE(square_hessian.value_[nonzero_el] == nonzero_value);

  square_hessian.value_[nonzero_el] +=
    0.5 * kSquareHessianAsymmetryTolerance;
  if (dev_run) square_hessian.print("Perturbed nonzero");
  // Perturbed nonzero yields warning
  REQUIRE(h.checkHessianOracle() == HighsStatus::kWarning);
  // Restore the nonzero value
  square_hessian.value_[nonzero_el] = nonzero_value;

  if (dev_run) square_hessian.print("Reversion 1");
  
  // Oracle should be OK
  const HighsStatus status = h.checkHessianOracle();
  printf("Check for square Hessian oracle has status %s\n",
	 h.highsStatusToString(status).c_str());
  REQUIRE(status == HighsStatus::kOk);
}

TEST_CASE("hessian-oracle-solve", "[qp-oracle]") {
  HighsHessian hessian;
  auto header = [&] (const HighsInt k, const std::string& message) {
    if (!dev_run) return;
    printf("\n==============================\n");
    printf("Case %d: %s\n", int(k), message.c_str());
    printf("==============================\n\n");
  };
  const HighsInt from_k = 0;
  for (HighsInt k = from_k; k < 4; k++) {
    if (k == 0) {
      header(k, "Hessian - diagonal of dimension 4");
      hessian = getHessianDiagonal4();
      testUnconConOracleSolve(hessian);
    } else if (k == 1) {
      header(k, "Hessian - diagonal of dimension 4 with zero entry");
      hessian = getHessianDiagonalWithZero4();
      testUnconConOracleSolve(hessian);
    } else if (k == 2) {
      header(k, "Hessian - of dimension 5");
      hessian = getHessian5();
      testUnconConOracleSolve(hessian);
    } else {
      header(k, "QpQjh");
      HighsModel model = getQpQjh();
      testOracleSolve(model);
    }
  }
}

HighsHessian getHessianDiagonal4() {
  HighsHessian hessian;
  hessian.dim_ = 4;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 1, 2, 3, 4};
  hessian.index_ = {0, 1, 2, 3};
  hessian.value_ = {1, 2, 3, 4};
  return hessian;
}

HighsHessian getHessianDiagonalWithZero4() {
  HighsHessian hessian;
  hessian.dim_ = 4;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 1, 2, 3, 4};
  hessian.index_ = {0, 1, 2, 3};
  hessian.value_ = {1, 2, 0, 4};
  return hessian;
}

HighsHessian getHessian4() {
  // Row|    0    1    2    3
  //-------------------------
  //   0|    4   -2         2
  //   1|   -2    2    1   -2
  //   2|         1    5    1
  //   3|    2         1    4
  HighsHessian hessian;
  hessian.dim_ = 4;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 3, 6, 8, 9};
  hessian.index_ = {0,  1, 3,  1, 2,  3,  2, 3,  3};
  hessian.value_ = {4, -2, 2,  2, 1, -2,  5, 1,  4};
  return hessian;
}

HighsHessian getHessian5() {
  // Row|    0    1    2    3    4
  //------------------------------
  //   0|    5    1        -1    2
  //   1|    1    4              1
  //   2|             10   -1
  //   3|   -1        -1    3   -2
  //   4|    2    1        -2    2
  HighsHessian hessian;
  hessian.dim_ = 5;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 4, 6, 8, 10, 11};
  hessian.index_ = {0, 1, 3, 4, 1, 4, 2, 3, 3, 4, 4};
  hessian.value_ = {5, 1, -1, 2, 4, 1, 10, -1, 3, -2, 2};
  return hessian;
}

HighsModel getQpQjh() {
  HighsModel model;
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;

  lp.model_name_ = "qjh";
  lp.num_col_ = 3;
  lp.num_row_ = 1;
  lp.col_cost_ = {0.0, -1.0, -3.0};
  lp.col_lower_ = {0.0, 0.0, 0.0};
  lp.col_upper_ = {inf, inf, inf};
  lp.sense_ = ObjSense::kMinimize;
  lp.offset_ = 0;
  lp.row_lower_ = {-inf};
  lp.row_upper_ = {2};
  lp.a_matrix_.start_ = {0, 1, 1, 2};
  lp.a_matrix_.index_ = {0, 0};
  lp.a_matrix_.value_ = {1.0, 1.0};
  lp.a_matrix_.format_ = MatrixFormat::kColwise;
  hessian.dim_ = lp.num_col_;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 2, 3, 4};
  hessian.index_ = {0, 2, 1, 2};
  hessian.value_ = {2.0, -1.0, 0.2, 2.0};

  return model;
}

void addVars(HighsModel& model) {
  HighsLp& lp = model.lp_;
  lp.num_col_ = model.hessian_.dim_;
  lp.num_row_ = 0;
  lp.col_cost_.assign(lp.num_col_, 0);
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, kHighsInf);
}

bool valuesRelEqual(const double v0, const double v1) {
  return std::fabs(v0-v1)/(1.0 + std::fabs(v0) + std::fabs(v1)) <
double_equal_tolerance;
}

bool vectorsRelEqual(const HighsInt dim, const double* v0, const double* v1) {
  for (HighsInt iCol = 0; iCol < dim; iCol++)
    if (!valuesRelEqual(v0[iCol], v1[iCol])) return false;
  return true;
}

void testUnconConOracleSolve(const HighsHessian& hessian) {

  HighsModel model;
  HighsLp& lp = model.lp_;
  model.hessian_ = hessian;
  addVars(model);
  // Set up the solution
  std::vector<double>solution(lp.num_col_);
  // First use a solution with a zero component
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    solution[iCol] = double(iCol);
  // Set up the costs
  lp.col_cost_.resize(lp.num_col_);
  hessian.product(solution, lp.col_cost_);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    lp.col_cost_[iCol] *= -1;

  testOracleSolve(model);

  // Now set up a solution without a zero component, and a constraint
  // that lies on it but cuts off the unconstrained minimizer
  std::vector<double>a_vector(lp.num_col_, 1);
  double rhs = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    solution[iCol] = double(iCol+1);
    rhs += solution[iCol] * a_vector[iCol];
  }
  lp.num_row_ = 1;
  lp.a_matrix_.num_col_ = lp.num_col_;
  lp.a_matrix_.num_row_ = 1;
  lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  lp.a_matrix_.start_ = {0, lp.num_col_};
  lp.a_matrix_.index_.clear();
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    lp.a_matrix_.index_.push_back(iCol);
  lp.a_matrix_.value_ = a_vector;
  lp.row_lower_ = {rhs};
  lp.row_upper_ = {kHighsInf};
  double row_dual = 1;
  hessian.product(solution, lp.col_cost_);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    lp.col_cost_[iCol] = -lp.col_cost_[iCol] + row_dual * a_vector[iCol];

  testOracleSolve(model);

}

void testOracleSolve(const HighsModel& model) {
  const HighsLp& lp = model.lp_;
  HighsHessian local_hessian = model.hessian_;
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("solver", kQpAsmString);
  const HighsInfo& info = h.getInfo();
  const double& objective_function_value = info.objective_function_value;
  double required_objective_function_value = 0;
  std::vector<double> required_solution(lp.num_col_, 0);
  // Four passes:
  //
  // 0: with triangular Hessian
  //
  // 1: with triangular Hessian and internal test_oracle
  //
  // 2: with triangular Hessian oracle
  //
  // 3: with square Hessian oracle
  
  for (HighsInt k = 0; k < 3; k++) {
    if (dev_run) {
      printf("\n========================================================\n");
      if (k == 0) {
	printf("Vanilla QPASM with triangular Hessian\n");
      }	else if (k == 1) {
	printf("QPASM with triangular Hessian and internal test_oracle\n");
      }	else if (k == 2) {
	printf("QPASM with square Hessian oracle\n");
      }	
      printf("========================================================\n\n");
    }
    const bool internal_test_oracle = k == 1;
    if (k == 0) {
      REQUIRE(h.passModel(model) == HighsStatus::kOk);
      //      h.writeModel("qp.mps");
    } else {
      REQUIRE(h.passModel(lp) == HighsStatus::kOk);
      if (internal_test_oracle) {
	REQUIRE(h.passHessian(local_hessian) == HighsStatus::kOk);
	REQUIRE(h.setOptionValue("test_qp_oracle", true) == HighsStatus::kOk);
      } else {
	void* oracle_data = &local_hessian;
	HighsStatus return_status;
	REQUIRE(local_hessian.format_ == HessianFormat::kSquare);
	return_status =
	  h.passHessian(lp.num_col_, oracleCallSquareHessian, oracle_data);
      }
    }
    //    if (dev_run) h.writeModel("");
    REQUIRE(h.run() == HighsStatus::kOk);

    if (dev_run) h.writeSolution("", 1);
    if (k == 0) {
      required_objective_function_value = objective_function_value;
      required_solution = h.getSolution().col_value;
    } else {
      REQUIRE(valuesRelEqual(objective_function_value,
			     required_objective_function_value));
      REQUIRE(vectorsRelEqual(lp.num_col_, h.getSolution().col_value.data(),
			      required_solution.data()));
    }
    if (internal_test_oracle)
      REQUIRE(h.setOptionValue("test_qp_oracle", false) == HighsStatus::kOk);
    if (local_hessian.format_ == HessianFormat::kSquare) break;
    if (k == 1) local_hessian = local_hessian.toSquare();
  }

  h.resetGlobalScheduler(true);
}
