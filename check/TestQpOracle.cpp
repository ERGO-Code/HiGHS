#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"
//#include "io/FilereaderLp.h"

const bool dev_run = true;//false;
const double inf = kHighsInf;
const double double_equal_tolerance = 1e-5;

HighsModel getQp();
bool vectorsEqual(const HighsInt dim, double* v0, double* v1);

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

HighsHessianFunctionType oracleCallTriangularHessian =
    [](const HighsInt x_num_entries, const HighsInt* x_index, const double* x_value,
       HighsInt& q_x_num_entries, HighsInt* q_x_index, double* q_x_value,
       void* hessian_p) {

      assert(x_value != nullptr);
      assert(q_x_value != nullptr);

      // Lambda for zeroing q_x_value
      auto zeroQx = [&] (const HighsInt dim) {
	for (HighsInt iCol = 0; iCol < dim; iCol++) q_x_value[iCol] = 0;
      };

      HighsHessian hessian = *(static_cast<HighsHessian*>(hessian_p));
      assert(hessian.format_ == HessianFormat::kTriangular);
      // With a triangular Hessian, have to scatter any packed values
      // of x unless only one Qx index is required
      const bool scatter = x_index != nullptr && (q_x_index == nullptr || q_x_num_entries != 1);
      std::vector<double>scattered_x;
      if (scatter) {
	scattered_x.assign(hessian.dim_, 0);
	for (HighsInt iX = 0; iX < x_num_entries; iX++)
	  scattered_x[x_index[iX]] = x_value[iX];
      }
      const double* use_x_value = scatter ? scattered_x.data() : x_value;
      
      // Lambda for adding multiple of Hessian column into q_x_value
      auto addScaledQcol = [&] (const HighsInt iCol) {
	HighsInt iEl = hessian.start_[iCol];
	for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol+1]; iEl++) {
	  HighsInt iRow = hessian.index_[iEl];
	  q_x_value[iRow] += hessian.value_[iEl] * use_x_value[iCol];
	  if (iRow != iCol) q_x_value[iCol] += hessian.value_[iEl] * use_x_value[iRow];
	}
      };

      if (x_index == nullptr) {
	// Simple product with full vector x, full vector q_x, and no
	// Qx indices required
	assert(q_x_index == nullptr);
	zeroQx(hessian.dim_);
	for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++) 
	  addScaledQcol(iCol);
	return;
      } else if (x_num_entries > 1) {
	// x is sparse with x_num_entries entries in rows x_index, and
	// no Qx indices required
	assert(q_x_index == nullptr);
	zeroQx(hessian.dim_);
	for (HighsInt iX = 0; iX < x_num_entries; iX++) 
	  addScaledQcol(x_index[iX]);
	return;
      } else if (x_num_entries == 1) {
	assert(q_x_index != nullptr);
	if (q_x_num_entries < 0) {
	  // x is sparse with one entry in row x_index, and all Qx index
	  // required
	  q_x_num_entries = 0;
	  // Get the entries below the diagonal in column x_index[0]
	  HighsInt iCol = x_index[0];
	  for (HighsInt iEl = hessian.start_[iCol]; iEl <hessian.start_[iCol+1]; iEl++) {
	    q_x_index[q_x_num_entries] = hessian.index_[iEl];
	    q_x_value[q_x_num_entries] = hessian.value_[iEl] * use_x_value[iCol];
	    q_x_num_entries++;
	  }
	  // Get the entries in row x_index[0] in previous columns
	  for (HighsInt iCol = 0; iCol < x_index[0]; iCol++) {
	    for (HighsInt iEl = hessian.start_[iCol]; iEl <hessian.start_[iCol+1]; iEl++) {
	      HighsInt iRow = hessian.index_[iEl];
	      if (iRow == x_index[0]) {
		q_x_index[q_x_num_entries] = iCol;
		q_x_value[q_x_num_entries] = hessian.value_[iEl] * use_x_value[iRow];
		q_x_num_entries++;
		break;
	      }
	    }
	  }
	  return;
	} else if (q_x_num_entries == 1) {
	  // x is sparse with one entry in row x_index, and one Qx index
	  // required
	  // With a triangular Hessian, need to identify which column to
	  // search down, and which row to look for
	  HighsInt iCol = std::min(x_index[0], q_x_index[0]);
	  HighsInt iRow = std::max(x_index[0], q_x_index[0]);
	  // Zero Qx value in case the Hessian entry requested is zero
	  q_x_value[0] = 0;
	  for (HighsInt iEl = hessian.start_[iCol]; iEl <hessian.start_[iCol+1]; iEl++) {
	    if (hessian.index_[iEl] == iRow) {
	      q_x_value[0] = hessian.value_[iEl] * use_x_value[0];
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

HighsHessianFunctionType oracleCallSquareHessian =
    [](const HighsInt x_num_entries, const HighsInt* x_index, const double* x_value, 
       HighsInt& q_x_num_entries, HighsInt* q_x_index, double* q_x_value, 
       void* hessian_p) {

      assert(x_value != nullptr);
      assert(q_x_value != nullptr);

      // Lambda for zeroing q_x_value
      auto zeroQx = [&] (const HighsInt dim) {
	for (HighsInt iCol = 0; iCol < dim; iCol++) q_x_value[iCol] = 0;
      };

      HighsHessian hessian = *(static_cast<HighsHessian*>(hessian_p));
      assert(hessian.format_ == HessianFormat::kSquare);
      
      // Lambda for adding multiple of Hessian column into q_x_value
      auto addScaledQcol = [&] (const HighsInt iCol) {
	for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol+1]; iEl++) {
	  HighsInt iRow = hessian.index_[iEl];
	  q_x_value[iRow] += hessian.value_[iEl] * x_value[iCol];
	}
      };

      if (x_index == nullptr) {
	// Simple product with full vector x, full vector q_x, and no
	// Qx indices required
	assert(q_x_index == nullptr);
	zeroQx(hessian.dim_);
	for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++) 
	  addScaledQcol(iCol);
	return;
      } else if (x_num_entries > 1) {
	// x is sparse with x_num_entries entries in rows x_index, and
	// no Qx indices required
	assert(q_x_index == nullptr);
	zeroQx(hessian.dim_);
	for (HighsInt iX = 0; iX < x_num_entries; iX++) 
	  addScaledQcol(x_index[iX]);
	return;
      } else if (x_num_entries == 1) {
	assert(q_x_index != nullptr);
	if (q_x_num_entries < 0) {
	  // x is sparse with one entry in row x_index, and all Qx index
	  // required
	  q_x_num_entries = 0;
	  // Get the entries below the diagonal in column iCol
	  HighsInt iCol = x_index[0];
	  for (HighsInt iEl = hessian.start_[iCol]; iEl <hessian.start_[iCol+1]; iEl++) {
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
	  for (HighsInt iEl = hessian.start_[iCol]; iEl <hessian.start_[iCol+1]; iEl++) {
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

TEST_CASE("hessian-oracle-check", "[qpsolver]") {
  HighsLp lp;
  lp.num_col_ = 5;
  lp.num_row_ = 0;
  lp.col_cost_.assign(lp.num_col_, 0);
  lp.col_lower_.assign(lp.num_col_, 0);
  lp.col_upper_.assign(lp.num_col_, 1);
  // Test using both triangular and square instances of this Hessian
  //
  // .  0  1  2  3  4
  // 0  5
  // 1  1  4
  // 2        0
  // 3 -1    -1  3
  // 4  2  1    -2  2
  HighsHessian hessian;
  hessian.dim_ = lp.num_col_;
  hessian.format_ = HessianFormat::kTriangular;
  hessian.start_ = {0, 4, 6, 8, 10, 11};
  hessian.index_ = {0, 1,  3, 4,  1, 4,  2,  3,  3,  4,  4};
  hessian.value_ = {5, 1, -1, 2,  4, 1,  0, -1,  3, -2,  2};
  HighsHessian square_hessian = hessian.toSquare();
  if (dev_run) square_hessian.print();

  // Now define the triangular instance without the explicit diagonal
  // zero or diagonal entry first in each packed column (cf col
  // 1). Can't do this from the outset, because HighsHessian::toSquare
  // assumes that the Hessian is triangular with the first entry in
  // each packed column being the diagonal entry
  hessian.start_ = {0, 4, 6, 7, 9, 10};
  hessian.index_ = {0, 1,  3, 4,  4, 1,   3,  3,  4,  4};
  hessian.value_ = {5, 1, -1, 2,  1, 4,  -1,  3, -2,  2};
  if (dev_run) hessian.print();

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  HighsStatus return_status = h.passModel(lp);
  REQUIRE(return_status == HighsStatus::kOk);

  // No oracle defined returns error
  REQUIRE(h.checkHessianOracle() == HighsStatus::kError);
  
  for (HighsInt k = 0; k < 2; k++) {
    // First pass is with square Hessian
    bool square = k == 0;
    if (square) {
      return_status = h.passHessian(lp.num_col_, oracleCallSquareHessian, &square_hessian);
    } else {
      return_status = h.passHessian(lp.num_col_, oracleCallTriangularHessian, &hessian);
    }
    REQUIRE(return_status == HighsStatus::kOk);

    // Check the method to extract Hessian entries from the oracle
    //
    // This is used to form the Hessian as a matrix in
    // Highs::checkHessianOracle(), so has to be tested externally
    const HessianOracle& oracle = h.getModel().hessian_.oracle_;
    std::vector<double> column;
    for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
      column.assign(lp.num_col_, 0);
      for (HighsInt iEl = hessian.start_[iCol]; iEl < hessian.start_[iCol+1]; iEl++) 
	column[hessian.index_[iEl]] = hessian.value_[iEl];
      REQUIRE(oracle.diag(iCol) == column[iCol]);
      for (HighsInt iRow = iCol+1; iRow < lp.num_col_; iRow++) {
	REQUIRE(oracle.entry(iRow, iCol) == column[iRow]);
	REQUIRE(oracle.entry(iCol, iRow) == column[iRow]);
      }
    }
    if (square) {
      // Test the asymmetry check:
      //
      // By replacing the zero (2, 4) entry with a nonzero
      HighsInt zero_el = square_hessian.numNz();
      HighsInt zero_entry_row = 2;
      HighsInt zero_entry_col = 4;
      REQUIRE(oracle.entry(zero_entry_row, zero_entry_col) == 0);

      // Replacing the zero (2, 4) entry with a nonzero
      square_hessian.start_[zero_entry_col+1]++;
      square_hessian.index_.push_back(zero_entry_row);
      square_hessian.value_.push_back(1.0);
      if (dev_run) square_hessian.print("Added nonzero");
      // Asymmetry yields error
      REQUIRE(h.checkHessianOracle(true) == HighsStatus::kError);

      // Remove the nonzero from the (2, 4) entry
      square_hessian.start_[zero_entry_col+1]--;
      square_hessian.index_.resize(zero_el);
      square_hessian.value_.resize(zero_el);
      if (dev_run) square_hessian.print("Reversion 0");

      // By perturbing the (0, 3) entry by less than the asymmetry
      // tolerance
      HighsInt nonzero_entry_row = 0;
      HighsInt nonzero_entry_col = 3;
      HighsInt nonzero_el = -1;
      for (HighsInt iEl = square_hessian.start_[nonzero_entry_col];
	   iEl < square_hessian.start_[nonzero_entry_col+1]; iEl++) {
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

      square_hessian.value_[nonzero_el] += 0.5 * kSquareHessianAsymmetryTolerance;
      if (dev_run) square_hessian.print("Perturbed nonzero");
      // Perturbed nonzero yields warning
      REQUIRE(h.checkHessianOracle() == HighsStatus::kWarning);
      // Restore the nonzero value
      square_hessian.value_[nonzero_el] = nonzero_value;

      if (dev_run) square_hessian.print("Reversion 1");
    
      
    }
    // Oracle should be OK
    const HighsStatus status = h.checkHessianOracle();
    printf("Check for %s Hessian oracle has status %s\n", square ? "square" : "triangular",
	   h.highsStatusToString(status).c_str());
    REQUIRE(status == HighsStatus::kOk);
  }
}

TEST_CASE("hessian-oracle-solve", "[qpsolver]") {
  HighsModel model = getQp();
  HighsLp lp = model.lp_;
  HighsHessian hessian = model.hessian_;

  Highs h;
  h.setOptionValue("output_flag", dev_run);
  h.setOptionValue("solver", kQpAsmString);
  const HighsInfo& info = h.getInfo();
  const double& objective_function_value = info.objective_function_value;

  HighsStatus return_status = h.passModel(lp);
  REQUIRE(return_status == HighsStatus::kOk);

  void* oracle_data = &hessian;

  return_status = h.passHessian(hessian.dim_, oracleCallTriangularHessian, oracle_data);
  REQUIRE(return_status == HighsStatus::kOk);

  if (dev_run) h.writeModel("");
  REQUIRE(h.writeModel("Null.mps") == HighsStatus::kError);
  return_status = h.run();
  REQUIRE(return_status == HighsStatus::kOk);

  double required_objective_function_value = -5.50;
  
  REQUIRE(fabs(required_objective_function_value-objective_function_value) < double_equal_tolerance);
  
  h.resetGlobalScheduler(true);
}

HighsModel getQp() {
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

bool vectorsEqual(const HighsInt dim, double* v0, double* v1) {
  for (HighsInt iCol = 0; iCol < dim; iCol++) 
    if (v0[iCol] != v1[iCol]) return false;
  return true;
}
