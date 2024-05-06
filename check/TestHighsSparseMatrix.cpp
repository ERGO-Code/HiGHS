#include <algorithm>
#include <cassert>

#include "HCheckConfig.h"
#include "Highs.h"
#include "catch.hpp"

const bool dev_run = false;
bool infNormDiffOk(const std::vector<double> x0, const std::vector<double> x1) {
  assert(x1.size() >= x0.size());
  double norm_diff = 0;
  for (HighsInt ix = 0; ix < HighsInt(x0.size()); ix++)
    norm_diff = std::max(std::abs(x0[ix] - x1[ix]), norm_diff);
  return norm_diff < 1e-12;
}

// No commas in test case name.
TEST_CASE("Sparse-matrix-products", "[highs_sparse_matrix]") {
  HighsStatus status;
  Highs highs;
  HighsRandom random;
  for (int k = 0; k < 2; k++) {
    std::string filename;
    if (k == 0) {
      filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
    } else {
      filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
    }
    if (!dev_run) highs.setOptionValue("output_flag", false);
    highs.readModel(filename);
    HighsSparseMatrix matrix = highs.getLp().a_matrix_;
    REQUIRE(matrix.isColwise());
    if (k == 0) {
      REQUIRE(matrix.num_row_ > matrix.num_col_);
    } else {
      REQUIRE(matrix.num_col_ > matrix.num_row_);
    }
    // Set up x and y for computing Ax and A^Ty
    std::vector<double> x;
    std::vector<double> y;
    for (HighsInt ix = 0; ix < matrix.num_col_; ix++)
      x.push_back(random.fraction());
    for (HighsInt ix = 0; ix < matrix.num_row_; ix++)
      y.push_back(random.fraction());
    std::vector<double> at_y;
    std::vector<double> a_x;
    // Form Ax
    std::vector<double> exact_a_x;
    matrix.product(exact_a_x, x);
    // Form A^Ty
    std::vector<double> exact_at_y;
    matrix.productTranspose(exact_at_y, y);
    // Now repeat with matrix rowwise
    matrix.ensureRowwise();

    matrix.product(a_x, x);
    REQUIRE(infNormDiffOk(a_x, exact_a_x));

    matrix.productTranspose(at_y, y);
    REQUIRE(infNormDiffOk(at_y, exact_at_y));

    // Now test the quad precision products
    matrix.ensureColwise();

    matrix.productQuad(a_x, x);
    REQUIRE(infNormDiffOk(a_x, exact_a_x));

    std::vector<HighsInt> at_y_index;
    std::vector<double> at_y_value;
    HVector y_HVector;
    y_HVector.setup(matrix.num_row_);
    y_HVector.array = y;
    matrix.productTransposeQuad(at_y_value, at_y_index, y_HVector);
    at_y.assign(matrix.num_col_, 0);
    for (HighsInt ix = 0; ix < int(at_y_index.size()); ix++)
      at_y[at_y_index[ix]] = at_y_value[ix];
    REQUIRE(infNormDiffOk(at_y, exact_at_y));

    // Now repeat with matrix rowwise
    matrix.ensureRowwise();
    matrix.productQuad(a_x, x);
    REQUIRE(infNormDiffOk(a_x, exact_a_x));

    at_y_index.clear();
    at_y_value.clear();
    matrix.productTransposeQuad(at_y_value, at_y_index, y_HVector);
    at_y.assign(matrix.num_col_, 0);
    for (HighsInt ix = 0; ix < int(at_y_index.size()); ix++)
      at_y[at_y_index[ix]] = at_y_value[ix];
    REQUIRE(infNormDiffOk(at_y, exact_at_y));

    // Now test y := y + alpha Ax and x := x + alphaA^Ty
    //
    // Test y := y + alpha Ax
    const double alpha = random.fraction();

    // First column-wise
    matrix.ensureColwise();
    for (HighsInt orientation = 0; orientation < 2; orientation++) {
      std::vector<double> initial_y;
      for (HighsInt ix = 0; ix < matrix.num_row_; ix++) initial_y.push_back(1);
      std::vector<double> exact_result;
      for (HighsInt ix = 0; ix < matrix.num_row_; ix++)
        exact_result.push_back(initial_y[ix] + alpha * exact_a_x[ix]);

      std::vector<double> result = initial_y;
      matrix.alphaProductPlusY(alpha, x, result);
      REQUIRE(infNormDiffOk(result, exact_result));

      // Test x := x + alphaA^Ty
      std::vector<double> initial_x;
      for (HighsInt ix = 0; ix < matrix.num_col_; ix++) initial_x.push_back(1);
      exact_result.clear();
      for (HighsInt ix = 0; ix < matrix.num_col_; ix++)
        exact_result.push_back(initial_x[ix] + alpha * exact_at_y[ix]);

      result = initial_x;
      matrix.alphaProductPlusY(alpha, y, result, true);
      REQUIRE(infNormDiffOk(result, exact_result));

      // Switch to row-wise for second pass
      matrix.ensureRowwise();
    }

    highs.clear();
  }
}
