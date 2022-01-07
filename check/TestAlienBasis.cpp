#include <sstream>

#include "Highs.h"
#include "catch.hpp"
#include "util/HFactor.h"
#include "util/HighsRandom.h"

const double inf = kHighsInf;
const bool dev_run = false;
const double double_equal_tolerance = 1e-5;

void testAlienBasis(const bool avgas, const HighsInt seed);
void getDependentCols(const HighsSparseMatrix& matrix,
                      std::vector<HighsInt>& col_set,
                      std::vector<HighsInt>& dependent_col_set,
                      const HighsInt required_rank_deficiency);
void reportColSet(const std::string message,
                  const std::vector<HighsInt>& col_set);
void reportDependentCols(const std::vector<HighsInt>& dependent_col_set);

TEST_CASE("AlienBasis-rank-detection", "[highs_test_alien_basis]") {
  // To find the dependent rows in
  //
  // [1 1  ]
  // [2 2  ]
  // [1   1]
  // [2   2]
  // [6 3 3]
  //
  // Define the transpose of the matrix column-wise
  HighsSparseMatrix matrix;
  matrix.num_col_ = 5;
  matrix.num_row_ = 3;
  matrix.start_ = {0, 2, 4, 6, 8, 11};
  matrix.index_ = {0, 1, 0, 1, 0, 2, 0, 2, 0, 1, 2};
  matrix.value_ = {1, 1, 2, 2, 1, 1, 2, 2, 6, 3, 3};
  std::vector<HighsInt> col_set;
  std::vector<HighsInt> dependent_col_set;
  HighsInt required_rank_deficiency;
  // getDependentCols uses HFactor::build() to determine col_set as
  // the maximal linearly independent subset of columns defined by the
  // original col_set, together with the indices of logical columns so
  // that the returned col_set is of full rank.
  //
  // For the transpose of the matrix above, here are the 5 columns and logical
  // columns
  //
  // 0 1 2 3 4 | 5 6 7
  // ----------+------
  // 1 2 1 2 6 | 1
  // 1 2     3 |   1
  //     1 2 3 |     1
  //
  // Data about the linearly dependent columns defined by the original
  // col_set is available in the following data members of HFactor,
  // that have size rank_deficiency
  //
  //   row_with_no_pivot: Rows in which no pivot was found when
  //   factorizing the matrix
  //
  //   col_with_no_pivot: Positions in col_set corresponding to
  //   columns in which no pivot was found when factorizing the matrix
  //
  //   var_with_no_pivot: Entries in col_set for which no pivot was
  //   found when factorizing the matrix
  //
  // The terms "var" and "col" relate to the set of basic variables
  // and columns of the basis matrix.
  //
  // Case 1
  // ======
  //
  // With all columns in the set, the distinction between "var" and
  // "col" is not demonstrated, but its still a case worth documenting
  //
  // The matrix defined by col_set has rank 2, so the rank deficiency is 3
  //
  col_set = {0, 1, 2, 3, 4};
  required_rank_deficiency = 3;
  if (dev_run) reportColSet("\nOriginal", col_set);
  getDependentCols(matrix, col_set, dependent_col_set,
                   required_rank_deficiency);
  if (dev_run) reportColSet("Returned", col_set);
  if (dev_run) reportDependentCols(dependent_col_set);
  //
  // The entries in the returned column set correspond to the matrix
  //
  // 4 6 3 | 8 9
  // ------+----
  // 6   2 |
  // 3 1   |
  // 3   2 |
  //       | 1
  //       |   1
  //
  // The last two entries can be ignored. They are fictitious logical
  // columns 8 and 9 so that the whole set of size 5 is a non-singular
  // 5x5 matrix
  //
  // Case 2
  // ======
  //
  // With a subset of columns in the set, particularly if the indices
  // are not ordered, the distinction between "var" and "col" in the
  // data produced by HFactor::build() is demonstrated
  //
  // The matrix defined by col_set has rank 2, so the rank deficiency is 2
  //
  col_set = {2, 0, 1, 3};
  required_rank_deficiency = 2;
  if (dev_run) reportColSet("\nOriginal", col_set);
  getDependentCols(matrix, col_set, dependent_col_set,
                   required_rank_deficiency);
  if (dev_run) reportColSet("Returned", col_set);
  if (dev_run) reportDependentCols(dependent_col_set);
  //
  // The entries in the returned column set correspond to the
  // matrix
  //
  // 1 6 3 | 8
  // ------+--
  // 1   2 |
  // 1 1   |
  //     2 |
  //       | 1
  //
  // The last entry can be ignored. It is a fictitious logical column
  // 8 (num_col is still 5!), so that the whole set of size 4 is a
  // non-singular 4x4 matrix
  //
}

TEST_CASE("AlienBasis-delay-singularity", "[highs_test_alien_basis]") {
  // Test the use of HFactor to complete a rectangular matrix when
  // (near) cancellation yields a (near-)zero row or column, to form a
  // nonsingular square matrix
  //
  // Set up a matrix with 6 columns, 5 rows and a column rank deficiency of 2
  HighsSparseMatrix matrix;
  matrix.num_col_ = 6;
  matrix.num_row_ = 5;
  matrix.start_ = {0, 2, 4, 9, 12, 15, 18};
  matrix.index_ = {0, 1, 0, 1, 0, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4, 2, 3, 4};
  matrix.value_ = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 4, 1, 2, 3, 1, -1, -1};
  const HighsInt perturbed_entry = 3;
  const double perturbation = -1e-12;
  matrix.value_[perturbed_entry] += perturbation;
  const HighsInt num_row = matrix.num_row_;
  const HighsInt num_col = matrix.num_col_;
  const HighsInt num_basic_col = num_col;
  HighsInt rank_deficiency;
  HighsInt required_rank_deficiency;
  required_rank_deficiency = 2;
  // The column set is all matrix columns except 2
  std::vector<HighsInt> col_set = {0, 1, 2, 3, 4, 5};

  if (dev_run) reportColSet("\nOriginal", col_set);
  HFactor factor;
  factor.setup(matrix, col_set);
  rank_deficiency = factor.build();
  REQUIRE(rank_deficiency == required_rank_deficiency);
}

TEST_CASE("AlienBasis-LP", "[highs_test_alien_basis]") {
  const HighsInt num_seed = 10;
  bool avgas = true;
  for (HighsInt seed = 0; seed < num_seed; seed++) testAlienBasis(avgas, seed);
  avgas = false;
  for (HighsInt seed = 0; seed < num_seed; seed++) testAlienBasis(avgas, seed);
}

void getDependentCols(const HighsSparseMatrix& matrix,
                      std::vector<HighsInt>& col_set,
                      std::vector<HighsInt>& dependent_col_set,
                      const HighsInt required_rank_deficiency) {
  HFactor factor;
  factor.setup(matrix, col_set);
  HighsInt rank_deficiency = factor.build();
  REQUIRE(rank_deficiency == required_rank_deficiency);
  if (dev_run) {
  }
  if (dev_run)
    printf("Returned rank_deficiency = %d:\n  No pivot in\nk Row Col Var\n",
           (int)rank_deficiency);
  dependent_col_set.clear();
  for (HighsInt k = 0; k < rank_deficiency; k++) {
    if (dev_run)
      printf("%1d %3d %3d %3d\n", (int)k, (int)factor.row_with_no_pivot[k],
             (int)factor.col_with_no_pivot[k],
             (int)factor.var_with_no_pivot[k]);
    dependent_col_set.push_back(factor.var_with_no_pivot[k]);
  }
}

void reportDependentCols(const std::vector<HighsInt>& dependent_col_set) {
  printf("Dependent column(s) in col_set:");
  for (HighsInt k = 0; k < (HighsInt)dependent_col_set.size(); k++)
    printf(" %1d", (int)dependent_col_set[k]);
  printf("\n");
}

void reportColSet(const std::string message,
                  const std::vector<HighsInt>& col_set) {
  printf("%s col_set:\n", message.c_str());
  for (HighsInt k = 0; k < (HighsInt)col_set.size(); k++)
    printf(" %1d", (int)col_set[k]);
  printf("\n");
}

void testAlienBasis(const bool avgas, const HighsInt seed) {
  std::string filename;
  std::string model;
  if (avgas) {
    model = "avgas";
  } else {
    model = "israel";
  }

  filename = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  std::stringstream ss;

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
  if (run_square_test && !seed) {
    ss.str(std::string());
    ss << "AlienBasis: " << model << " square";
    basis.debug_origin_name = ss.str();
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
  const bool run_square_random_test = true;
  if (run_square_random_test) {
    ss.str(std::string());
    ss << "AlienBasis: " << model << " random-" << seed << " square";
    basis.debug_origin_name = ss.str();
    // Create a full-dimension basis using random selection of num_col variables
    basis.col_status.assign(num_col, HighsBasisStatus::kNonbasic);
    basis.row_status.assign(num_row, HighsBasisStatus::kNonbasic);
    HighsRandom random(seed);
    HighsInt num_basic = 0;
    for (;;) {
      HighsInt iVar = random.integer(num_var);
      if (iVar < num_col) {
        if (basis.col_status[iVar] == HighsBasisStatus::kNonbasic) {
          basis.col_status[iVar] = HighsBasisStatus::kBasic;
          num_basic++;
        }
      } else {
        if (basis.row_status[iVar - num_col] == HighsBasisStatus::kNonbasic) {
          basis.row_status[iVar - num_col] = HighsBasisStatus::kBasic;
          num_basic++;
        }
      }
      if (num_basic == num_row) break;
    }
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }

  std::string profile = num_col < num_row ? "portrait" : "landscape";
  const bool run_primal_test = true;
  if (run_primal_test && !seed) {
    // Create a rectangular basis using just struturals
    ss.str(std::string());
    ss << "AlienBasis: " << model << " primal " << profile;
    basis.debug_origin_name = ss.str();
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      basis.col_status[iCol] = HighsBasisStatus::kBasic;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      basis.row_status[iRow] = HighsBasisStatus::kNonbasic;
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }
  const bool run_primal_random_test = true;
  if (run_primal_random_test) {
    // Create a rectangular basis using random selection of num_col variables
    basis.col_status.assign(num_col, HighsBasisStatus::kNonbasic);
    basis.row_status.assign(num_row, HighsBasisStatus::kNonbasic);
    ss.str(std::string());
    ss << "AlienBasis: " << model << " primal random-" << seed << " "
       << profile;
    basis.debug_origin_name = ss.str();
    HighsRandom random(seed);
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
  if (run_dual_test && !seed) {
    // Create a rectangular basis using just struturals
    ss.str(std::string());
    ss << "AlienBasis: " << model << " dual " << profile;
    basis.debug_origin_name = ss.str();
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      basis.col_status[iCol] = HighsBasisStatus::kBasic;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      basis.row_status[iRow] = HighsBasisStatus::kNonbasic;
    REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
    highs.run();
  }
  const bool run_dual_random_test = true;
  if (run_dual_random_test) {
    // Create a rectangular basis using random selection of num_col variables
    basis.col_status.assign(num_col, HighsBasisStatus::kNonbasic);
    basis.row_status.assign(num_row, HighsBasisStatus::kNonbasic);
    basis.debug_origin_name =
        "AlienBasis: " + model + " dual random " + profile;
    ss.str(std::string());
    ss << "AlienBasis: " << model << " dual random-" << seed << " " << profile;
    basis.debug_origin_name = ss.str();
    HighsRandom random(seed);
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
