#include <algorithm>
#include <cstdio>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;  // true;  //
const double inf = kHighsInf;
const double double_equal_tolerance = 1e-5;

const std::vector<std::string> solvers{kQpAsmString
#ifdef HIPO
                                       ,
                                       kHipoString
#endif
};

HighsHessian getHessianDiagonal4();
HighsHessian getHessianDiagonalWithZero4();
HighsHessian getHessian5();
HighsModel getQpQjh();
bool valuesRelEqual(const double v0, const double v1);
bool vectorsRelEqual(const HighsInt dim, const double* v0, const double* v1);
void addVars(HighsModel& model);
void testUnconConOracleSolve(const HighsHessian& model);
void testOracleSolve(const HighsModel& model);

// User-supplied Hessian oracle is called with one of three value for call_type:
//
// kHessianOracleCallTypeEntry
//
// Set *hessian_x_value as Hessian entry (*x_index, *hessian_x_index)
//
// Return 0 if the Hessian entry is available, otherwise, return a nonzero value
//
// kHessianOracleCallTypeColumn
//
// Set *hessian_x_num_entries, *hessian_x_index and *hessian_x_value as Hessian
// column *x_index, where *hessian_x_index contains the *hessian_x_num_entries
// indices of the column nonzeros, and their values are assumed to be scattered
// in *hessian_x_value
//
// Return 0 if the Hessian column is available, otherwise, return a nonzero
// value
//
// kHessianOracleCallTypeProduct
//
// Set *hessian_x_value as the values of the product between the Hessian and
// vector which is either
//
// * a full vector *x_value (if x_index is nullptr)
//
// * a scattered vector *x_value with *x_num_entries of nonzeros in *x_index
//
// Return 0 if the Hessian product is available, otherwise, return a
// nonzero value - in which case the oracle is not valid
//
HighsHessianFunctionType oracleCallSquareHessian =
    [](const HighsInt call_type, const HighsInt* x_num_entries,
       const HighsInt* x_index, const double* x_value,
       HighsInt* hessian_x_num_entries, HighsInt* hessian_x_index,
       double* hessian_x_value, void* hessian_p) {
      assert(kHessianOracleCallTypeMin <= call_type &&
             call_type <= kHessianOracleCallTypeMax);

      HighsHessian hessian = *(static_cast<HighsHessian*>(hessian_p));
      assert(hessian.format_ == HessianFormat::kSquare);

      // Lambda for adding multiple of Hessian column into hessian_x_value
      auto addScaledQcol = [&](const HighsInt iCol, const double x_value) {
        for (HighsInt iEl = hessian.start_[iCol];
             iEl < hessian.start_[iCol + 1]; iEl++) {
          HighsInt iRow = hessian.index_[iEl];
          hessian_x_value[iRow] += hessian.value_[iEl] * x_value;
        }
      };

      if (call_type == kHessianOracleCallTypeEntry) {
        assert(x_num_entries == nullptr);
        assert(x_value == nullptr);
        assert(x_index != nullptr);
        assert(hessian_x_num_entries == nullptr);
        assert(hessian_x_index != nullptr);
        assert(hessian_x_value != nullptr);
        HighsInt iCol = x_index[0];
        HighsInt iRow = hessian_x_index[0];
        // Zero Qx value in case the Hessian entry requested is zero
        hessian_x_value[0] = 0;
        for (HighsInt iEl = hessian.start_[iCol];
             iEl < hessian.start_[iCol + 1]; iEl++) {
          if (hessian.index_[iEl] == iRow) {
            hessian_x_value[0] = hessian.value_[iEl];
            return 0;
          }
        }
      } else if (call_type == kHessianOracleCallTypeColumn) {
        // Get the entries in column iCol
        assert(x_num_entries == nullptr);
        assert(x_value == nullptr);
        assert(x_index != nullptr);
        assert(hessian_x_num_entries != nullptr);
        assert(hessian_x_index != nullptr);
        assert(hessian_x_value != nullptr);
        (*hessian_x_num_entries) = 0;
        HighsInt iCol = x_index[0];
        for (HighsInt iEl = hessian.start_[iCol];
             iEl < hessian.start_[iCol + 1]; iEl++) {
          hessian_x_index[*hessian_x_num_entries] = hessian.index_[iEl];
          hessian_x_value[*hessian_x_num_entries] = hessian.value_[iEl];
          (*hessian_x_num_entries)++;
        }
      } else {
        //	printf("x_index == nullptr = %s\n", x_index == nullptr ? "T" :
        //"F"); 	printf("*x_num_entries = %d: *x_num_entries >= 0 = %s \n",
        //int(*x_num_entries), *x_num_entries >= 0 ? "T" : "F");
        assert(x_index == nullptr || *x_num_entries >= 0);
        assert(hessian_x_index == nullptr);
        assert(hessian_x_value != nullptr);
        if (x_index == nullptr) {
          // Simple product with full vector x, full vector q_x
          for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++)
            addScaledQcol(iCol, x_value[iCol]);
        } else {
          // x is scattered with x_num_entries entries in rows x_index
          for (HighsInt iX = 0; iX < *x_num_entries; iX++) {
            HighsInt iCol = x_index[iX];
            addScaledQcol(iCol, x_value[iCol]);
          }
        }
      }
      return 0;
    };

TEST_CASE("hessian-oracle-check", "[qp-oracle]") {
  Highs h;
  h.setOptionValue("output_flag", dev_run);

  HighsModel model;
  model.hessian_ = getHessian5();
  // Add an LP so that the Hessian can be passed
  addVars(model);
  HighsLp& lp = model.lp_;
  HighsHessian& hessian = model.hessian_;

  HighsInt zero_diagonal_el;
  HighsInt zero_diagonal_col = 2;
  zero_diagonal_el = 6;
  REQUIRE(hessian.index_[zero_diagonal_el] == zero_diagonal_col);
  double zero_diagonal_original_value = hessian.value_[zero_diagonal_el];
  REQUIRE(zero_diagonal_original_value != 0);

  h.passModel(model);
  h.setOptionValue("qp_regularization_value", 0);
  h.run();
  h.writeSolution("", 1);

  hessian.value_[zero_diagonal_el] = 0.0;

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

  REQUIRE(h.passHessian(lp.num_col_, oracleCallSquareHessian,
                        &square_hessian) == HighsStatus::kOk);

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

  square_hessian.value_[nonzero_el] += 0.5 * kSquareHessianAsymmetryTolerance;
  if (dev_run) square_hessian.print("Perturbed nonzero");
  // Perturbed nonzero yields warning
  REQUIRE(h.checkHessianOracle() == HighsStatus::kWarning);
  // Restore the nonzero value
  square_hessian.value_[nonzero_el] = nonzero_value;

  if (dev_run) square_hessian.print("Reversion 1");

  // Oracle should be OK
  HighsStatus status = h.checkHessianOracle();
  if (dev_run)
    printf("Check for square Hessian oracle has status %s\n\n",
           h.highsStatusToString(status).c_str());
  REQUIRE(status == HighsStatus::kOk);

  // Define a customisable oracle to test recovery from absence of
  // entry or column calls, error handling and code coveage

  bool no_entry_call = false;
  bool no_column_call = false;
  bool no_product_call = false;
  bool column_error = false;
  bool product_error = false;

  HighsHessianFunctionType oracleCallSquareHessianCustomised =
      [&](const HighsInt call_type, const HighsInt* x_num_entries,
          const HighsInt* x_index, const double* x_value,
          HighsInt* hessian_x_num_entries, HighsInt* hessian_x_index,
          double* hessian_x_value, void* hessian_p) {
        assert(kHessianOracleCallTypeMin <= call_type &&
               call_type <= kHessianOracleCallTypeMax);

        HighsHessian hessian = *(static_cast<HighsHessian*>(hessian_p));
        assert(hessian.format_ == HessianFormat::kSquare);

        // Lambda for adding multiple of Hessian column into hessian_x_value
        auto addScaledQcol = [&](const HighsInt iCol, const double x_value) {
          for (HighsInt iEl = hessian.start_[iCol];
               iEl < hessian.start_[iCol + 1]; iEl++) {
            HighsInt iRow = hessian.index_[iEl];
            hessian_x_value[iRow] += hessian.value_[iEl] * x_value;
          }
        };

        if (call_type == kHessianOracleCallTypeEntry) {
          assert(x_num_entries == nullptr);
          assert(x_value == nullptr);
          assert(x_index != nullptr);
          assert(hessian_x_num_entries == nullptr);
          assert(hessian_x_index != nullptr);
          assert(hessian_x_value != nullptr);
          if (no_entry_call) {
            hessian_x_value[0] = 0;
            return -1;
          }
          HighsInt iCol = x_index[0];
          HighsInt iRow = hessian_x_index[0];
          // Zero Qx value in case the Hessian entry requested is zero
          hessian_x_value[0] = 0;
          for (HighsInt iEl = hessian.start_[iCol];
               iEl < hessian.start_[iCol + 1]; iEl++) {
            if (hessian.index_[iEl] == iRow) {
              hessian_x_value[0] = hessian.value_[iEl];
              return 0;
            }
          }
        } else if (call_type == kHessianOracleCallTypeColumn) {
          // Get the entries in column iCol
          assert(x_num_entries == nullptr);
          assert(x_value == nullptr);
          assert(x_index != nullptr);
          assert(hessian_x_num_entries != nullptr);
          assert(hessian_x_index != nullptr);
          assert(hessian_x_value != nullptr);
          if (no_column_call) return -1;
          (*hessian_x_num_entries) = 0;
          HighsInt iCol = x_index[0];
          if (column_error && iCol == 1) return 0;
          for (HighsInt iEl = hessian.start_[iCol];
               iEl < hessian.start_[iCol + 1]; iEl++) {
            hessian_x_index[*hessian_x_num_entries] = hessian.index_[iEl];
            hessian_x_value[*hessian_x_num_entries] = hessian.value_[iEl];
            (*hessian_x_num_entries)++;
          }
        } else {
          assert(x_index == nullptr || *x_num_entries >= 0);
          assert(hessian_x_index == nullptr);
          assert(hessian_x_value != nullptr);
          if (no_product_call) return -1;
          if (product_error && *x_num_entries == 4) return 0;
          if (x_index == nullptr) {
            // Simple product with full vector x, full vector q_x
            for (HighsInt iCol = 0; iCol < hessian.dim_; iCol++)
              addScaledQcol(iCol, x_value[iCol]);
          } else {
            // x is scattered with x_num_entries entries in rows x_index
            for (HighsInt iX = 0; iX < *x_num_entries; iX++) {
              HighsInt iCol = x_index[iX];
              addScaledQcol(iCol, x_value[iCol]);
            }
          }
        }
        return 0;
      };

  REQUIRE(h.passHessian(lp.num_col_, oracleCallSquareHessianCustomised,
                        &square_hessian) == HighsStatus::kOk);

  auto checkLog = [&](const std::string& insert) {
    if (dev_run)
      printf("Check for customised square Hessian oracle %shas status %s\n\n",
             insert.c_str(), h.highsStatusToString(status).c_str());
  };
  // Default oracle should be OK
  status = h.checkHessianOracle();
  checkLog("");
  REQUIRE(status == HighsStatus::kOk);

  // With a column error
  column_error = true;
  status = h.checkHessianOracle();
  checkLog("with column error ");
  REQUIRE(status == HighsStatus::kError);
  column_error = false;

  // With a product error
  product_error = true;
  status = h.checkHessianOracle();
  checkLog("with product error ");
  REQUIRE(status == HighsStatus::kError);
  product_error = false;

  // With no entry call or column call, the Hessian oracle is valid

  // With no entry call
  no_entry_call = true;
  status = h.checkHessianOracle();
  checkLog("with no entry call ");
  REQUIRE(status == HighsStatus::kOk);
  no_entry_call = false;

  // With no column call
  no_column_call = true;
  status = h.checkHessianOracle();
  checkLog("with no column call ");
  REQUIRE(status == HighsStatus::kOk);
  no_column_call = false;

  // With no entry call and column call
  no_entry_call = true;
  no_column_call = true;
  status = h.checkHessianOracle();
  checkLog("with no entry or column call ");
  REQUIRE(status == HighsStatus::kOk);
  no_entry_call = false;
  no_column_call = false;

  // With no product call, the Hessian oracle is not valid
  no_product_call = true;
  status = h.checkHessianOracle();
  checkLog("with no product call ");
  REQUIRE(status == HighsStatus::kError);

  // Code coverage of code testing for passing an invalid Hessian oracle
  status = h.passHessian(0, oracleCallSquareHessianCustomised, &square_hessian);
  REQUIRE(status == HighsStatus::kWarning);
  if (dev_run)
    printf(
        "Check for passing Hessian oracle with zero dimension has "
        "status %s\n\n",
        h.highsStatusToString(status).c_str());

  status = h.passHessian(lp.num_col_, oracleCallSquareHessianCustomised,
                         &square_hessian);
  REQUIRE(status == HighsStatus::kError);
  if (dev_run)
    printf(
        "Check for passing invalid Hessian oracle has "
        "status %s\n\n",
        h.highsStatusToString(status).c_str());
}

TEST_CASE("hessian-oracle-solve", "[qp-oracle]") {
  HighsHessian hessian;
  auto header = [&](const HighsInt k, const std::string& message) {
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

TEST_CASE("hessian-oracle-primal1", "[qp-oracle]") {
  const std::string test_name = Catch::getResultCapture().getCurrentTestName();
  std::string write_model_filename = test_name + ".mps";
  std::string filename =
      std::string(HIGHS_DIR) + "/check/instances/primal1.mps";
  Highs h;
  h.setOptionValue("output_flag", dev_run);
  REQUIRE(h.readModel(filename) == HighsStatus::kOk);
  HighsLp lp = h.getModel().lp_;
  HighsHessian square_hessian = h.getModel().hessian_.toSquare();
  REQUIRE(h.passModel(lp) == HighsStatus::kOk);
  REQUIRE(h.passHessian(lp.num_col_, oracleCallSquareHessian,
                        &square_hessian) == HighsStatus::kOk);
  double optimal_obective_value;
  std::vector<double> solution;
  for (auto& solver : solvers) {
    h.setOptionValue("solver", solver);
    h.run();
    if (solver == kQpAsmString) {
      REQUIRE(h.writeModel(write_model_filename) == HighsStatus::kError);
      optimal_obective_value = h.getObjectiveValue();
      solution = h.getSolution().col_value;
    } else {
      REQUIRE(valuesRelEqual(h.getObjectiveValue(), optimal_obective_value));
      REQUIRE(vectorsRelEqual(lp.num_col_, h.getSolution().col_value.data(),
                              solution.data()));
    }
  }
  //  std::remove(write_model_filename.c_str());

  h.resetGlobalScheduler(true);
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
  return std::fabs(v0 - v1) / (1.0 + std::fabs(v0) + std::fabs(v1)) <
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
  std::vector<double> solution(lp.num_col_);
  // First use a solution with a zero component
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++)
    solution[iCol] = double(iCol);
  // Set up the costs
  lp.col_cost_.resize(lp.num_col_);
  hessian.product(solution, lp.col_cost_);
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) lp.col_cost_[iCol] *= -1;

  testOracleSolve(model);

  // Now set up a solution without a zero component, and a constraint
  // that lies on it but cuts off the unconstrained minimizer
  std::vector<double> a_vector(lp.num_col_, 1);
  double rhs = 0;
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    solution[iCol] = double(iCol + 1);
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
      } else if (k == 1) {
        printf("QPASM with triangular Hessian and internal test_oracle\n");
      } else if (k == 2) {
        printf("QPASM with square Hessian oracle\n");
      }
      printf("========================================================\n\n");
    }
    const bool internal_test_oracle = k == 1;
    if (k == 0) {
      REQUIRE(h.passModel(model) == HighsStatus::kOk);
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
