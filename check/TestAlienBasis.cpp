#include "Highs.h"
#include "catch.hpp"
#include "util/HighsRandom.h"

const double inf = kHighsInf;
const bool dev_run = true;
const double double_equal_tolerance = 1e-5;

void testAlienBasis(const bool avgas);
TEST_CASE("AlienBasis-avgas", "[highs_test_alien_basis]") {
  bool avgas = true;
  testAlienBasis(avgas);
  avgas = false;
  testAlienBasis(avgas);
}

void testAlienBasis(const bool avgas) {
  std::string filename;
  std::string model;
  if (avgas) {
    model = "avgas";
  } else {
    model = "israel";
  }

  filename = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  HighsLp lp = highs.getLp();
  HighsInt num_col = lp.num_col_;
  HighsInt num_row = lp.num_row_;
  // Assumes that the test LP has fewer columns than rows
  // (portrait). Lansdcape test is performed on its dual.
  assert(num_col < num_row);
  const HighsInt num_var = num_col + num_row;
  HighsBasis basis;
  basis.col_status.resize(num_col);
  basis.row_status.resize(num_row);
  const bool run_square_test = true;
  if (run_square_test) {
    basis.debug_origin_name = "AlienBasis: " + model + " square";
    // Create a full-dimension basis using struturals and then enough logicals
    HighsBasisStatus status = HighsBasisStatus::kBasic;
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      if (iCol >= num_row) status = HighsBasisStatus::kNonbasic;
      basis.col_status[iCol] = status;
    }
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (num_col + iRow >= num_row) status = HighsBasisStatus::kNonbasic;
      basis.row_status[iRow] = status;
    }
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }

  std::string profile = num_col < num_row ? "portrait" : "landscape";
  const bool run_primal_test = true;
  if (run_primal_test) {
    // Create a rectangular basis using just struturals
    basis.debug_origin_name = "AlienBasis: " + model + " primal " + profile;
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      basis.col_status[iCol] = HighsBasisStatus::kBasic;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      basis.row_status[iRow] = HighsBasisStatus::kNonbasic;
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }
  const bool run_primal_random_test = false;
  if (run_primal_random_test) {
    // Create a rectangular basis using random selection of num_col variables
    basis.col_status.assign(num_col, HighsBasisStatus::kNonbasic);
    basis.row_status.assign(num_row, HighsBasisStatus::kNonbasic);
    basis.debug_origin_name =
        "AlienBasis: " + model + " primal random " + profile;
    HighsRandom random;
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      HighsInt iVar = random.integer(num_var);
      if (iVar < num_col) {
        basis.col_status[iVar] = HighsBasisStatus::kBasic;
      } else {
        basis.row_status[iVar - num_col] = HighsBasisStatus::kBasic;
      }
    }
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }

  // Test opposite profile using dual LP.
  //
  // Primal must be either
  //
  // min c^Tx s.t. Ax >= b; x >= 0 (avgas-primal)
  //
  // min -b^Ty s.t. A^Ty <= c; y >= 0 (avgas-dual)
  //
  // Or
  //
  // min c^Tx s.t. Ax <= b; x >= 0 (israel-primal)
  //
  // min b^Ty s.t. A^Ty >= -c; y >= 0 (israel-dual)
  //

  HighsLp dual_lp;
  dual_lp.num_col_ = lp.num_row_;
  dual_lp.num_row_ = lp.num_col_;
  dual_lp.sense_ = ObjSense::kMinimize;
  dual_lp.col_lower_.assign(dual_lp.num_col_, 0);
  dual_lp.col_upper_.assign(dual_lp.num_col_, inf);
  if (lp.row_lower_[0] > -inf) {
    // avgas
    for (HighsInt iCol = 0; iCol < dual_lp.num_col_; iCol++)
      dual_lp.col_cost_.push_back(-lp.row_lower_[iCol]);
    dual_lp.row_lower_.assign(dual_lp.num_row_, -inf);
    dual_lp.row_upper_ = lp.col_cost_;
  } else {
    // israel
    dual_lp.col_cost_ = lp.row_upper_;
    for (HighsInt iRow = 0; iRow < dual_lp.num_row_; iRow++)
      dual_lp.row_lower_.push_back(-lp.col_cost_[iRow]);
    dual_lp.row_upper_.assign(dual_lp.num_row_, inf);
  }
  dual_lp.a_matrix_ = lp.a_matrix_;
  dual_lp.a_matrix_.num_col_ = dual_lp.num_col_;
  dual_lp.a_matrix_.num_row_ = dual_lp.num_row_;
  dual_lp.a_matrix_.format_ = MatrixFormat::kRowwise;
  dual_lp.a_matrix_.ensureColwise();
  highs.passModel(dual_lp);

  num_col = dual_lp.num_col_;
  num_row = dual_lp.num_row_;
  basis.col_status.resize(num_col);
  basis.row_status.resize(num_row);
  profile = num_col < num_row ? "portrait" : "landscape";
  const bool run_dual_test = true;
  if (run_dual_test) {
    // Create a rectangular basis using just struturals
    basis.debug_origin_name = "AlienBasis: " + model + " dual " + profile;
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      basis.col_status[iCol] = HighsBasisStatus::kBasic;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      basis.row_status[iRow] = HighsBasisStatus::kNonbasic;
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }
  const bool run_dual_random_test = false;
  if (run_dual_random_test) {
    // Create a rectangular basis using random selection of num_col variables
    basis.col_status.assign(num_col, HighsBasisStatus::kNonbasic);
    basis.row_status.assign(num_row, HighsBasisStatus::kNonbasic);
    basis.debug_origin_name =
        "AlienBasis: " + model + " dual random " + profile;
    HighsRandom random;
    for (HighsInt iCol = 0; iCol < num_col; iCol++) {
      HighsInt iVar = random.integer(num_var);
      if (iVar < num_col) {
        basis.col_status[iVar] = HighsBasisStatus::kBasic;
      } else {
        basis.row_status[iVar - num_col] = HighsBasisStatus::kBasic;
      }
    }
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }
}
