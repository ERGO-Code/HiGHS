#include "HCheckConfig.h"
#include "catch.hpp"
#include "lp_data/HighsOptions.h"
#include "model/HighsHessian.h"
#include "model/HighsHessianUtils.h"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("HighsHessian", "[highs_hessian]") {
  HighsOptions options;
  if (!dev_run) options.output_flag = false;

  HighsHessian square_hessian;
  // .  0  1  2  3  4
  // 0  5  1    -1  2
  // 1  1  4        1
  // 2        3 -1
  // 3 -1    -1  4
  // 4  2  1        5
  square_hessian.dim_ = 5;
  square_hessian.format_ = HessianFormat::kSquare;
  square_hessian.start_ = {0, 4, 7, 9, 12, 15};
  square_hessian.index_ = {0, 1, 3, 4, 0, 1, 4, 2, 3, 0, 2, 3, 0, 1, 4};
  square_hessian.value_ = {5, 1, -1, 2, 1, 4, 1, 3, -1, -1, -1, 4, 2, 1, 5};
  HighsHessian square_hessian0 = square_hessian;

  HighsHessian triangular_hessian;
  // .  0  1  2  3  4
  // 0  5
  // 1  1  4
  // 2        3
  // 3 -1    -1  4
  // 4  2  1        5
  triangular_hessian.dim_ = 5;
  triangular_hessian.format_ = HessianFormat::kTriangular;
  triangular_hessian.start_ = {0, 4, 6, 8, 9, 10};
  triangular_hessian.index_ = {0, 1, 3, 4, 1, 4, 2, 3, 3, 4};
  triangular_hessian.value_ = {5, 1, -1, 2, 4, 1, 3, -1, 4, 5};
  HighsHessian triangular_hessian0 = triangular_hessian;

  // Check that the positive diagonal entries are recognised as being
  // OK for default assessHessian (minimization) but not for
  // maximization.
  REQUIRE(assessHessian(square_hessian, options) == HighsStatus::kOk);
  REQUIRE(okHessianDiagonal(options, square_hessian, ObjSense::kMinimize));
  REQUIRE(!okHessianDiagonal(options, square_hessian, ObjSense::kMaximize));
  if (dev_run) {
    printf("\nReturned square Hessian\n");
    square_hessian.print();
  }

  // Check that the positive diagonal entries are recognised as being
  // OK for default assessHessian (minimization) but not for
  // maximization.
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  REQUIRE(okHessianDiagonal(options, square_hessian, ObjSense::kMinimize));
  REQUIRE(!okHessianDiagonal(options, square_hessian, ObjSense::kMaximize));
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }

  REQUIRE((triangular_hessian == triangular_hessian0));

  // Extract the triangluar Hessian from the square Hessian
  REQUIRE(extractTriangularHessian(options, square_hessian) ==
          HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangularised square Hessian\n");
    square_hessian.print();
  }

  // Extract the triangluar Hessian from the triangular Hessian
  REQUIRE(extractTriangularHessian(options, triangular_hessian) ==
          HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangularised triangular Hessian\n");
    triangular_hessian.print();
  }

  HighsHessian negative_diagonal_hessian = triangular_hessian;
  negative_diagonal_hessian.value_[0] = -negative_diagonal_hessian.value_[0];
  negative_diagonal_hessian.value_[4] = -negative_diagonal_hessian.value_[4];
  negative_diagonal_hessian.value_[6] = -negative_diagonal_hessian.value_[6];
  negative_diagonal_hessian.value_[8] = -negative_diagonal_hessian.value_[8];
  negative_diagonal_hessian.value_[9] = -negative_diagonal_hessian.value_[9];
  REQUIRE(assessHessian(negative_diagonal_hessian, options) ==
          HighsStatus::kOk);
  REQUIRE(!okHessianDiagonal(options, negative_diagonal_hessian,
                             ObjSense::kMinimize));
  REQUIRE(okHessianDiagonal(options, negative_diagonal_hessian,
                            ObjSense::kMaximize));

  // Square Hessian with only triangular entries - doubled strictly triangular
  // entries.
  HighsHessian hessian0;
  hessian0.dim_ = 5;
  hessian0.format_ = HessianFormat::kSquare;
  hessian0.start_ = {0, 1, 3, 4, 7, 10};
  hessian0.index_ = {0, 0, 1, 2, 0, 2, 3, 0, 1, 4};
  hessian0.value_ = {5, 2, 4, 3, -2, -2, 4, 4, 2, 5};

  //  REQUIRE(assessHessian(hessian0, options) == HighsStatus::kOk);
  REQUIRE(assessHessian(hessian0, options) == HighsStatus::kError);
  if (dev_run) {
    printf("\nReturned\n");
    hessian0.print();
  }
  //  REQUIRE((hessian0 == triangular_hessian0));

  // Nonsymmetric Hessian - with entries resulting in cancellation
  HighsHessian hessian1;
  hessian1.format_ = HessianFormat::kSquare;
  hessian1.dim_ = 5;
  hessian1.start_ = {0, 3, 5, 7, 10, 14};
  hessian1.index_ = {0, 3, 4, 0, 1, 2, 4, 0, 2, 3, 0, 1, 2, 4};
  hessian1.value_ = {5, -5, 1, 2, 4, 3, 1, 3, -2, 4, 3, 2, -1, 5};
  //  REQUIRE(assessHessian(hessian1, options) == HighsStatus::kOk);
  REQUIRE(assessHessian(hessian1, options) == HighsStatus::kError);
  if (dev_run) {
    printf("\nReturned\n");
    hessian1.print();
  }
  //  REQUIRE((hessian1 == triangular_hessian0));

  HighsHessian indefinite;
  indefinite.dim_ = 3;
  indefinite.format_ = HessianFormat::kTriangular;
  indefinite.start_ = {0, 2, 2, 3};
  indefinite.index_ = {0, 2, 2};
  indefinite.value_ = {2, 1, 4};
  HighsInt indefinite_num_nz0 = indefinite.numNz();
  REQUIRE(assessHessian(indefinite, options) == HighsStatus::kOk);
  // Check that there is one more entry due to the explicit zero
  REQUIRE(indefinite.numNz() == indefinite_num_nz0 + 1);
  // Check that all first indices are for the diagonal
  for (HighsInt iCol = 0; iCol < indefinite.dim_; iCol++) {
    const HighsInt iEl = indefinite.start_[iCol];
    REQUIRE(indefinite.index_[iEl] == iCol);
  }
  // Negate the diagonal entries
  for (HighsInt iCol = 0; iCol < indefinite.dim_; iCol++) {
    const HighsInt iEl = indefinite.start_[iCol];
    indefinite.value_[iEl] = -indefinite.value_[iEl];
  }
  if (dev_run) {
    printf("\nIndefinite\n");
    indefinite.print();
  }
  REQUIRE(assessHessian(indefinite, options) == HighsStatus::kOk);
  // Can tell that the indefinite Hessian is not OK for minimization
  REQUIRE(!okHessianDiagonal(options, indefinite, ObjSense::kMinimize));
  // Cannot tell that the indefinite Hessian is not OK for
  // maximization since its diagonal entries are in [-4 0)
  REQUIRE(okHessianDiagonal(options, indefinite, ObjSense::kMaximize));

  // Now add tests for the new normaliseHessian method for Hessians
  // supplied by users or read from MPS files

  HighsInt square_hessian_04_el = 12;
  square_hessian = square_hessian0;
  assert(square_hessian.index_[square_hessian_04_el] == 0);
  assert(square_hessian.value_[square_hessian_04_el] == 2);

  // As defined, the triangular-format matrix is equivalent to
  // triangular_hessian0 (above), but it has three explicit zeros
  // which together with the upper triangular nonzero, can be accessed
  // directly
  //
  // .  0  1  2  3  4
  // 0  5
  // 1  1  4  0     1
  // 2        3
  // 3 -1    -1  4  0
  // 4  2  0        5
  triangular_hessian.dim_ = 5;
  triangular_hessian.format_ = HessianFormat::kTriangular;
  triangular_hessian.start_ = {0, 4, 6, 9, 10, 13};
  triangular_hessian.index_ = {0, 1, 3, 4, 1, 4, 1, 2, 3, 3, 1, 3, 4};
  triangular_hessian.value_ = {5, 1, -1, 2, 4, 0, 0, 3, -1, 4, 1, 0, 5};
  HighsInt triangular_hessian_41_el = 5;
  HighsInt triangular_hessian_12_el = 6;
  HighsInt triangular_hessian_14_el = 10;
  HighsInt triangular_hessian_34_el = 11;
  assert(triangular_hessian.index_[triangular_hessian_41_el] == 4);
  assert(triangular_hessian.value_[triangular_hessian_41_el] == 0);

  assert(triangular_hessian.index_[triangular_hessian_12_el] == 1);
  assert(triangular_hessian.value_[triangular_hessian_12_el] == 0);

  assert(triangular_hessian.index_[triangular_hessian_14_el] == 1);
  assert(triangular_hessian.value_[triangular_hessian_14_el] == 1);

  assert(triangular_hessian.index_[triangular_hessian_34_el] == 3);
  assert(triangular_hessian.value_[triangular_hessian_34_el] == 0);

  HighsHessian triangular_hessian1 = triangular_hessian;

  if (dev_run) {
    printf("\nOriginal\n");
    triangular_hessian.print();
  }
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }
  REQUIRE(triangular_hessian == triangular_hessian0);

  // Now replace the explicit zero in (4, 1) by -1 and the 1 in (1, 4)
  // by 2 so they are summed to give the original 1 in (1, 4)
  triangular_hessian = triangular_hessian1;
  triangular_hessian.value_[triangular_hessian_41_el] = -1;
  triangular_hessian.value_[triangular_hessian_14_el] = 2;

  if (dev_run) {
    printf("\nOriginal\n");
    triangular_hessian.print();
  }
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }
  REQUIRE(triangular_hessian == triangular_hessian0);

  // Now replace the explicit zeros in (1, 2) and (3, 4) by 1
  triangular_hessian = triangular_hessian1;
  triangular_hessian.value_[triangular_hessian_12_el] = -1;
  triangular_hessian.value_[triangular_hessian_34_el] = -2;

  if (dev_run) {
    printf("\nOriginal\n");
    triangular_hessian.print();
  }
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }

  triangular_hessian = square_hessian0;
  triangular_hessian.format_ = HessianFormat::kTriangular;

  if (dev_run) {
    printf("\nSquare Hessian as triangular original\n");
    triangular_hessian.print();
  }
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }

  // Define the triangular matrix as its transpose
  //
  // .  0  1  2  3  4
  // 0  5
  // 1  1  4  0     1
  // 2        3
  // 3 -1    -1  4  0
  // 4  2  0        5
  triangular_hessian.dim_ = 5;
  triangular_hessian.format_ = HessianFormat::kTriangular;
  triangular_hessian.start_ = {0, 1, 5, 6, 10, 13};
  triangular_hessian.index_ = {0, 0, 1, 2, 4, 2, 0, 2, 3, 4, 0, 1, 4};
  triangular_hessian.value_ = {5, 1, 4, 0, 1, 3, -1, -1, 4, 0, 2, 0, 5};
  triangular_hessian_41_el = 11;
  triangular_hessian_12_el = 3;
  triangular_hessian_14_el = 4;
  triangular_hessian_34_el = 9;
  assert(triangular_hessian.index_[triangular_hessian_41_el] == 1);
  assert(triangular_hessian.value_[triangular_hessian_41_el] == 0);

  assert(triangular_hessian.index_[triangular_hessian_12_el] == 2);
  assert(triangular_hessian.value_[triangular_hessian_12_el] == 0);

  assert(triangular_hessian.index_[triangular_hessian_14_el] == 4);
  assert(triangular_hessian.value_[triangular_hessian_14_el] == 1);

  assert(triangular_hessian.index_[triangular_hessian_34_el] == 4);
  assert(triangular_hessian.value_[triangular_hessian_34_el] == 0);

  HighsHessian triangular_hessian2 = triangular_hessian;

  if (dev_run) {
    printf("\nOriginal\n");
    triangular_hessian.print();
  }
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }
  REQUIRE(triangular_hessian == triangular_hessian0);

  // Now replace the explicit zero in (4, 1) by -1 and the 1 in (1, 4)
  // by 2 so they are summed to give the original 1 in (1, 4)
  triangular_hessian = triangular_hessian2;
  triangular_hessian.value_[triangular_hessian_41_el] = -1;
  triangular_hessian.value_[triangular_hessian_14_el] = 2;

  if (dev_run) {
    printf("\nOriginal\n");
    triangular_hessian.print();
  }
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }
  REQUIRE(triangular_hessian == triangular_hessian0);

  // Now replace the explicit zeros in (1, 2) and (3, 4) by 1
  triangular_hessian = triangular_hessian2;
  triangular_hessian.value_[triangular_hessian_12_el] = -1;
  triangular_hessian.value_[triangular_hessian_34_el] = -2;

  if (dev_run) {
    printf("\nOriginal\n");
    triangular_hessian.print();
  }
  REQUIRE(assessHessian(triangular_hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned triangular Hessian\n");
    triangular_hessian.print();
  }
}
