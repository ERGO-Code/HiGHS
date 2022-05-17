#include "catch.hpp"
#include "lp_data/HighsOptions.h"
#include "model/HighsHessian.h"
#include "model/HighsHessianUtils.h"
//#include "<cstdio>"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("HighsHessian", "[highs_hessian]") {
  HighsOptions options;
  if (!dev_run) options.output_flag = false;

  HighsHessian square_hessian;
  square_hessian.dim_ = 5;
  square_hessian.format_ = HessianFormat::kSquare;
  square_hessian.start_ = {0, 4, 7, 9, 12, 15};
  square_hessian.index_ = {0, 1, 3, 4, 0, 1, 4, 2, 3, 0, 2, 3, 0, 1, 4};
  square_hessian.value_ = {5, 1, -1, 2, 1, 4, 1, 3, -1, -1, -1, 4, 2, 1, 5};

  HighsHessian triangular_hessian;
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

  REQUIRE(assessHessian(hessian0, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned\n");
    hessian0.print();
  }
  REQUIRE((hessian0 == triangular_hessian0));

  // Nonsymmetric Hessian - with entries resulting in cancellation
  HighsHessian hessian1;
  hessian1.format_ = HessianFormat::kSquare;
  hessian1.dim_ = 5;
  hessian1.start_ = {0, 3, 5, 7, 10, 14};
  hessian1.index_ = {0, 3, 4, 0, 1, 2, 4, 0, 2, 3, 0, 1, 2, 4};
  hessian1.value_ = {5, -5, 1, 2, 4, 3, 1, 3, -2, 4, 3, 2, -1, 5};

  REQUIRE(assessHessian(hessian1, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned\n");
    hessian1.print();
  }
  REQUIRE((hessian1 == triangular_hessian0));

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
}
