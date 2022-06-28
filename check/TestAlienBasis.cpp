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

// Generally HFactor is initialised with a rectangular constraint
// matrix, and a subset of basic columns defined by
// basic_index.
//
// In these examples, an instance of HighsSparseMatrix is used to
// define the coefficient matrix, and col_set is used to generalise
// the name of the identifier basic_index, where col_set.size()
// defines the number of entries in the set.
//
// The examples use
//
// factor.setup(matrix, col_set);
//
// HighsInt rank_deficiency = factor.build();
//
// However, the traditional pointer-based HFactor::setup is still
// available for a rectangular constraint matrix, and a subset of
// basic columns defined by basic_index.
//
//  void setup(const HighsInt num_col,
//             const HighsInt num_row,
//             const HighsInt* a_start,
//             const HighsInt* a_index,
//             const double* a_value,
//             HighsInt* basic_index);
//
// Also, for rank-detection of rectangular matrices, and completion of
// tall rectangular matrices into an LU factorization, there is
//
//  void setupGeneral(const HighsInt num_col,
//                    const HighsInt num_row,
//                    const HighsInt num_basic,
//                    const HighsInt* a_start,
//                    const HighsInt* a_index,
//                    const double* a_value,
//                    HighsInt* basic_index);
//

TEST_CASE("Hessian-rank-detection", "[highs_test_alien_basis]") {
  // To find the rank and corresponding rows/columns in a square
  // symmetric matrix
  //
  // [1     2      ]
  // [  0          ]
  // [    1       1]
  // [2     4      ]
  // [        1   2]
  // [          0  ]
  // [    1   2   8]
  //
  // NB HFactor has no facility to handle symmetric matrices
  // efficiently
  //
  // The matrix must not contain explicit zero values, and there is no
  // check in HFactor::setup
  //
  HighsSparseMatrix matrix;
  matrix.num_col_ = 7;
  matrix.num_row_ = 7;
  matrix.start_ = {0, 2, 2, 4, 6, 8, 8, 11};
  matrix.index_ = {0, 3, 2, 6, 0, 3, 4, 6, 2, 4, 6};
  matrix.value_ = {1, 2, 1, 1, 2, 4, 1, 2, 1, 2, 8};
  std::vector<HighsInt> col_set = {0, 1, 2, 3, 4, 5, 6};
  std::vector<HighsInt> dependent_col_set;
  const HighsInt required_rank_deficiency = 3;
  if (dev_run) reportColSet("\nOriginal", col_set);
  getDependentCols(matrix, col_set, dependent_col_set,
                   required_rank_deficiency);
  if (dev_run) reportColSet("Returned", col_set);
  if (dev_run) reportDependentCols(dependent_col_set);
}

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

TEST_CASE("AlienBasis-rectangular-completion", "[highs_test_alien_basis]") {
  // Test the use of HFactor to complete a rectangular matrix in order
  // to form a nonsingular square matrix
  //
  // Set up a matrix with 5 columns, 6 rows and a rank deficiency of 1
  HighsSparseMatrix matrix;
  matrix.num_col_ = 5;
  matrix.num_row_ = 6;
  matrix.start_ = {0, 4, 6, 12, 17, 20};
  matrix.index_ = {0, 2, 4, 5, 0, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 4, 5, 1, 4, 5};
  matrix.value_ = {1, 1, 2,  1, 2,  1, 1, 1, 1, 1,
                   1, 1, -1, 2, -1, 2, 1, 1, 2, 1};
  const HighsInt num_row = matrix.num_row_;
  const HighsInt num_col = matrix.num_col_;
  const HighsInt num_basic_col = num_col - 1;
  HighsInt rank_deficiency;
  HighsInt required_rank_deficiency;
  required_rank_deficiency = 1;
  // The column set is all matrix columns except 2
  std::vector<HighsInt> col_set = {0, 1, 3, 4};
  HFactor factor;
  factor.setup(matrix, col_set);
  rank_deficiency = factor.build();
  REQUIRE(rank_deficiency == required_rank_deficiency);
  if (dev_run) reportColSet("Returned", col_set);
  // Note that col_set already has the index 4 = num_col+0 to replace
  // the deficient column 2 with the logical 0.
  if (dev_run)
    printf("Returned rank_deficiency = %d:\n  No pivot in\nk Row Col Var\n",
           (int)rank_deficiency);
  if (dev_run) {
    // Report on the row with no pivot, index of the deficient entry
    // in col_set, and entry itself.
    for (HighsInt k = 0; k < rank_deficiency; k++)
      printf("%1d %3d %3d %3d\n", (int)k, (int)factor.row_with_no_pivot[k],
             (int)factor.col_with_no_pivot[k],
             (int)factor.var_with_no_pivot[k]);
  }
  // Now illustrate how col_set can be extended with two more logicals
  for (HighsInt k = rank_deficiency; k < num_row - num_basic_col + 1; k++) {
    if (dev_run) printf("%1d %3d\n", (int)k, (int)factor.row_with_no_pivot[k]);
    // Identify the index of the logical that is required
    const HighsInt introduce_logical = factor.row_with_no_pivot[k];
    const HighsInt introduce_column = num_col + introduce_logical;
    col_set.push_back(introduce_column);
  }
  // Need to call HFactor::setup again as the dimension of col_set has
  // changed
  factor.setup(matrix, col_set);
  required_rank_deficiency = 0;
  rank_deficiency = factor.build();
  REQUIRE(rank_deficiency == required_rank_deficiency);
  // Demonstrate the existence of a factorizaion
  HighsRandom random;
  vector<double> solution;
  HVector rhs;
  rhs.setup(num_row);
  rhs.clear();
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    HighsInt iVar = col_set[iRow];
    double solution_value = random.fraction();
    solution.push_back(solution_value);
    if (iVar < num_col) {
      for (HighsInt iEl = matrix.start_[iVar]; iEl < matrix.start_[iVar + 1];
           iEl++)
        rhs.array[matrix.index_[iEl]] += solution_value * matrix.value_[iEl];
    } else {
      rhs.array[iVar - num_col] += solution_value;
    }
  }
  std::iota(rhs.index.begin(), rhs.index.end(), 0);
  rhs.count++;
  factor.ftranCall(rhs, 1);
  double solution_error = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    solution_error += std::abs(rhs.array[iRow] - solution[iRow]);
  if (dev_run)
    printf("AlienBasis-rectangular-completion: solution_error = %g\n",
           solution_error);
  fflush(stdout);

  REQUIRE(solution_error < 1e-8);
}

TEST_CASE("AlienBasis-delay-singularity0", "[highs_test_alien_basis]") {
  // Test the use of HFactor to complete a rectangular matrix when
  // (near) cancellation yields a (near-)zero row or column, to form a
  // nonsingular square matrix.
  //
  // Set up a matrix with 6 columns, 5 rows and a column rank deficiency of 2
  //
  // Generates a singleton column with pivot 1e-12, then a singleton
  // row with pivot 1e-12
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

TEST_CASE("AlienBasis-delay-singularity1", "[highs_test_alien_basis]") {
  // Test the use of HFactor to complete a rectangular matrix when
  // (near) cancellation yields a (near-)zero row or column, to form a
  // nonsingular square matrix
  //
  // Set up a matrix with 5 columns, 5 rows and a column rank deficiency of 1
  HighsSparseMatrix matrix;
  matrix.num_col_ = 5;
  matrix.num_row_ = 5;
  matrix.start_ = {0, 5, 10, 14, 17, 20};
  matrix.index_ = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 1, 2, 3, 4, 2, 3, 4, 2, 3, 4};
  matrix.value_ = {1, 1, 2, -1, 3, 1, 1, 2, -1, 3,
                   1, 1, 2, 4,  1, 1, 1, 1, 2,  3};
  const HighsInt perturbed_from_entry = 6;
  const HighsInt perturbed_to_entry = 10;
  const double perturbation_multiplier = 1 + 1e-12;
  for (HighsInt iEl = perturbed_from_entry; iEl < perturbed_to_entry; iEl++)
    matrix.value_[iEl] *= perturbation_multiplier;
  const HighsInt num_row = matrix.num_row_;
  const HighsInt num_col = matrix.num_col_;
  const HighsInt num_basic_col = num_col;
  HighsInt rank_deficiency;
  HighsInt required_rank_deficiency;
  required_rank_deficiency = 1;
  // The column set is all matrix columns except 2
  std::vector<HighsInt> col_set = {0, 1, 2, 3, 4};

  if (dev_run) reportColSet("\nOriginal", col_set);
  HFactor factor;
  factor.setup(matrix, col_set);
  rank_deficiency = factor.build();
  REQUIRE(rank_deficiency == required_rank_deficiency);
}

TEST_CASE("AlienBasis-qap10", "[highs_test_alien_basis]") {
  // Alien basis example derived from sub-MIP when solving qap10
  //
  // Has logical in col_set that duplicates structural column
  HighsSparseMatrix matrix;
  matrix.num_col_ = 15;
  matrix.num_row_ = 3;
  matrix.start_ = {0, 3, 5, 6, 8, 9, 10, 12, 14, 17, 20, 23, 26, 29, 31, 32};
  matrix.index_ = {2, 0, 1, 2, 0, 0, 2, 0, 0, 0, 0, 2, 0, 2, 2, 0,
                   1, 1, 0, 2, 0, 1, 2, 2, 1, 0, 2, 1, 0, 2, 0, 1};
  matrix.value_ = {1,     0.5, 0.5, 0.75, 1,    1, 0.5, 1,    1,  1,   2,
                   0.5,   1,   -1,  1,    1,    1, 1,   1,    -1, 0.5, 0.5,
                   1.125, 1.5, 1,   1,    -0.5, 2, 2,   -0.5, 1,  1};

  std::vector<HighsInt> col_set = {0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 15};
  const HighsInt required_rank_deficiency = 10;

  if (dev_run) reportColSet("\nOriginal", col_set);
  HFactor factor;
  factor.setup(matrix, col_set);
  HighsInt rank_deficiency = factor.build();
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

TEST_CASE("AlienBasis-reuse-basis", "[highs_test_alien_basis]") {
  HighsLp lp;
  lp.num_col_ = 2;
  lp.num_row_ = 3;
  lp.col_cost_ = {400, 650};
  lp.col_lower_ = {0, 0};
  lp.col_upper_ = {inf, inf};
  lp.row_lower_ = {-inf, -inf, -inf};
  lp.row_upper_ = {140, 14, 85};
  lp.a_matrix_.start_ = {0, 3, 6};
  lp.a_matrix_.index_ = {0, 1, 2, 0, 1, 2};
  lp.a_matrix_.value_ = {15, 2, 10, 25, 3, 20};
  lp.sense_ = ObjSense::kMaximize;
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.passModel(lp);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
  HighsBasis basis = highs.getBasis();
  // Add another variable
  vector<HighsInt> new_index = {0, 1, 2};
  vector<double> new_value = {50, 4, 30};
  highs.addCol(850, 0, inf, 3, &new_index[0], &new_value[0]);
  // Add a new constraint
  new_value[0] = 15;
  new_value[1] = 24;
  new_value[2] = 30;
  highs.addRow(-inf, 108, 3, &new_index[0], &new_value[0]);
  const bool singlar_also = true;
  if (singlar_also) {
    const HighsInt from_col = 0;
    const HighsInt to_col = 0;
    HighsInt get_num_col;
    double get_cost;
    double get_lower;
    double get_upper;
    HighsInt get_num_nz;
    highs.getCols(from_col, to_col, get_num_col, &get_cost, &get_lower,
                  &get_upper, get_num_nz, NULL, NULL, NULL);
    vector<HighsInt> get_start(get_num_col + 1);
    vector<HighsInt> get_index(get_num_nz);
    vector<double> get_value(get_num_nz);
    highs.getCols(from_col, to_col, get_num_col, &get_cost, &get_lower,
                  &get_upper, get_num_nz, &get_start[0], &get_index[0],
                  &get_value[0]);

    // Make the first two columns parallel, so that the saved basis is
    // singular, as well as having too few basic variables
    REQUIRE(highs.changeCoeff(0, 1, 30) == HighsStatus::kOk);
    REQUIRE(highs.changeCoeff(1, 1, 4) == HighsStatus::kOk);
    REQUIRE(highs.changeCoeff(2, 1, 20) == HighsStatus::kOk);
    REQUIRE(highs.changeCoeff(3, 1, 30) == HighsStatus::kOk);
  }

  if (dev_run) highs.setOptionValue("log_dev_level", 3);
  highs.setOptionValue("simplex_scale_strategy", 0);
  // Make the basis status for new row and column nonbasic
  basis.col_status.push_back(HighsBasisStatus::kNonbasic);
  basis.row_status.push_back(HighsBasisStatus::kNonbasic);
  basis.alien = true;
  REQUIRE(highs.setBasis(basis) == HighsStatus::kOk);
  highs.run();
  if (dev_run) highs.writeSolution("", 1);
}
