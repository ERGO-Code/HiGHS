#include "catch.hpp"
#include "lp_data/HighsOptions.h"
#include "model/HighsHessian.h"
#include "model/HighsHessianUtils.h"
//#include "<cstdio>"

const bool dev_run = false;

// No commas in test case name.
TEST_CASE("HighsHessian", "[highs_hessian]") {
  HighsHessian hessian;
  hessian.dim_ = 5;
  hessian.q_start_ = {0, 4, 7, 9, 12, 15};
  hessian.q_index_ = {0, 1, 3, 4, 0, 1, 4, 2, 3, 0, 2, 3, 0, 1, 4};
  hessian.q_value_ = {5, 1, -1, 2, 1, 4, 1, 3, -1, -1, -1, 4, 2, 1, 5};
  HighsHessian hessian_copy = hessian;
  HighsOptions options;

  if (!dev_run) options.output_flag = false;
  REQUIRE(assessHessian(hessian, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned\n");
    hessian.print();
  }
  REQUIRE((hessian == hessian_copy));

  // Pure upper triangular Hessian - with doubled strictly upper
  // triangular entries.
  HighsHessian hessian0;
  hessian0.dim_ = 5;
  hessian0.q_start_ = {0, 1, 3, 4, 7, 10};
  hessian0.q_index_ = {0, 0, 1, 2, 0, 2, 3, 0, 1, 4};
  hessian0.q_value_ = {5, 2, 4, 3, -2, -2, 4, 4, 2, 5};

  REQUIRE(assessHessian(hessian0, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned\n");
    hessian0.print();
  }
  REQUIRE((hessian0 == hessian));

  // Nonsymmetric Hessian - with entries resulting in cancellation
  HighsHessian hessian1;
  hessian1.dim_ = 5;
  hessian1.q_start_ = {0, 3, 5, 7, 10, 14};
  hessian1.q_index_ = {0, 3, 4, 0, 1, 2, 4, 0, 2, 3, 0, 1, 2, 4};
  hessian1.q_value_ = {5, -5, 1, 2, 4, 3, 1, 3, -2, 4, 3, 2, -1, 5};

  REQUIRE(assessHessian(hessian1, options) == HighsStatus::kOk);
  if (dev_run) {
    printf("\nReturned\n");
    hessian1.print();
  }
  REQUIRE((hessian1 == hessian));
}
