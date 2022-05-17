#include "Highs.h"
#include "catch.hpp"
#include "util/HFactor.h"

const bool dev_run = false;

HVector rhs;
HVector col_aq;
HVector row_ep;
std::vector<HighsInt> basic_set;
std::vector<double> solution;
HighsLp lp;
HighsInt num_col;
HighsInt num_row;
HighsInt basis_change;
HFactor factor;
HighsInt rowOut(const HighsInt variable_out);
bool iterate(const HighsInt variable_out, const HighsInt variable_in);
bool testSolve();
bool testSolveDense();

TEST_CASE("Factor-dense-tran", "[highs_test_factor]") {
  std::string filename;
  const bool avgas = false;  // true;//
  std::string model = avgas ? "avgas" : "adlittle";
  filename = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  lp = highs.getLp();
  num_row = lp.num_row_;
  highs.run();
  basic_set.resize(num_row);
  // Get the optimal set of basic variables
  highs.getBasicVariables(&basic_set[0]);
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    basic_set[iRow] =
        basic_set[iRow] < 0 ? num_col - basic_set[iRow] + 1 : basic_set[iRow];
  factor.setup(lp.a_matrix_, basic_set);
  factor.build();
  solution.resize(num_row);
  HighsRandom random;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    solution[iRow] = random.fraction();
  rhs.setup(num_row);
  REQUIRE(testSolveDense());
}

TEST_CASE("Factor-put-get-iterate", "[highs_test_factor]") {
  std::string filename;
  const bool avgas = false;  // true;//
  std::string model = avgas ? "avgas" : "adlittle";
  filename = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";
  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  // Switch off presolve so that recovery of DSE weights can be
  // assessed
  highs.setOptionValue("presolve", kHighsOffString);
  highs.run();
  HighsSolution solution = highs.getSolution();
  HighsBasis basis = highs.getBasis();
  const HighsLp& lp = highs.getLp();
  double value;
  double integer_value;
  double lower;
  double upper;
  double save_lower;
  double save_upper;
  HighsInt num_test = 0;
  const HighsInt max_num_test = 10;
  bool put_iterate = false;
  if (dev_run) {
    highs.setOptionValue("highs_debug_level", 2);
    highs.setOptionValue("log_dev_level", 2);
  }
  for (HighsInt iCol = 0; iCol < lp.num_col_; iCol++) {
    value = solution.col_value[iCol];
    integer_value = round(value);
    if (std::fabs(value - integer_value) < 1e-1) continue;
    assert(basis.col_status[iCol] == HighsBasisStatus::kBasic);
    save_lower = lp.col_lower_[iCol];
    save_upper = lp.col_upper_[iCol];
    if (value < integer_value) {
      lower = integer_value;
      upper = lp.col_upper_[iCol];
    } else {
      lower = integer_value + 1;
      upper = lp.col_upper_[iCol];
    }
    if (lower > upper) continue;
    if (dev_run)
      printf(
          "\nChanging bounds on column %d (value %g) from [%g, %g] to [%g, "
          "%g]\n",
          (int)iCol, value, save_lower, save_upper, lower, upper);
    highs.changeColBounds(iCol, lower, upper);
    if (!put_iterate) {
      REQUIRE(highs.putIterate() == HighsStatus::kOk);
    } else {
      REQUIRE(highs.getIterate() == HighsStatus::kOk);
    }
    num_test++;
    highs.run();
    if (num_test == max_num_test) break;
    highs.changeColBounds(iCol, save_lower, save_upper);
  }
}
TEST_CASE("Factor-get-set-invert", "[highs_test_factor]") {
  std::string filename;
  const bool avgas = false;  // true;//
  std::string model = avgas ? "avgas" : "adlittle";
  filename = std::string(HIGHS_DIR) + "/check/instances/" + model + ".mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);
  highs.readModel(filename);
  lp = highs.getLp();
  num_col = lp.num_col_;
  num_row = lp.num_row_;
  std::vector<HighsInt> variable_out;
  std::vector<HighsInt> variable_in;
  if (avgas) {
    variable_out = {16, 9, 15, 12, 8, 14};
    variable_in = {5, 2, 0, 4, 3, 6};
  } else {
    variable_out = {97,  151, 124, 101, 138, 130, 102, 143, 146, 140, 142,
                    116, 48,  1,   126, 134, 144, 117, 69,  3,   110, 101,
                    31,  30,  56,  100, 139, 129, 128, 127, 53,  150, 114,
                    131, 113, 111, 108, 136, 63,  120, 106, 112, 123, 59,
                    45,  71,  141, 132, 135, 121, 4,   144, 26,  145, 137,
                    98,  52,  6,   125, 134, 105, 27,  55,  147, 90,  103,
                    64,  3,   99,  50,  58,  80,  117, 119};
    variable_in = {1,  69, 76, 95,  75, 71, 48,  56, 3,   77, 80, 6,  50,
                   55, 30, 31, 64,  53, 72, 101, 3,  134, 90, 51, 0,  2,
                   61, 60, 59, 117, 52, 47, 63,  35, 38,  26, 41, 4,  144,
                   25, 44, 29, 45,  32, 24, 5,   68, 66,  56, 94, 67, 91,
                   27, 7,  58, 18,  69, 92, 31,  63, 12,  14, 6,  74, 30,
                   11, 49, 79, 53,  81, 42, 82,  58, 1};
  }
  HighsRandom random;
  solution.resize(num_row);
  basic_set.clear();
  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    solution[iRow] = random.fraction();
    basic_set.push_back(num_col + iRow);
  }
  rhs.setup(num_row);
  col_aq.setup(num_row);
  row_ep.setup(num_row);
  factor.setup(lp.a_matrix_, basic_set);
  factor.build();
  HighsInt from_basis_change = 0;
  HighsInt to_basis_change = avgas ? 3 : 65;
  for (basis_change = from_basis_change; basis_change < to_basis_change;
       basis_change++) {
    if (basis_change == 50) factor.build();
    REQUIRE(iterate(variable_out[basis_change], variable_in[basis_change]));
  }

  std::vector<HighsInt> get_basic_set = basic_set;
  InvertibleRepresentation invert = factor.getInvert();
  std::vector<InvertibleRepresentation> invert_set;
  from_basis_change = to_basis_change;
  to_basis_change = variable_out.size();
  for (basis_change = from_basis_change; basis_change < to_basis_change;
       basis_change++) {
    REQUIRE(iterate(variable_out[basis_change], variable_in[basis_change]));
    invert_set.push_back(factor.getInvert());
  }
  basis_change = -1;
  basic_set = get_basic_set;
  factor.setInvert(invert);
  REQUIRE(testSolve());

  for (basis_change = from_basis_change; basis_change < to_basis_change;
       basis_change++)
    REQUIRE(iterate(variable_out[basis_change], variable_in[basis_change]));
}

HighsInt rowOut(const HighsInt variable_out) {
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    if (basic_set[iRow] == variable_out) return iRow;
  return -1;
}

bool iterate(const HighsInt variable_out, const HighsInt variable_in) {
  const HighsInt row_out = rowOut(variable_out);
  if (row_out < 0) return false;
  row_ep.clear();
  row_ep.count = 1;
  row_ep.index[0] = row_out;
  row_ep.array[row_out] = 1;
  row_ep.packFlag = true;
  factor.btranCall(row_ep, 1);

  col_aq.clear();
  col_aq.packFlag = true;
  lp.a_matrix_.collectAj(col_aq, variable_in, 1);
  factor.ftranCall(col_aq, 1);

  basic_set[row_out] = variable_in;
  HighsInt rebuild_reason = 0;
  HighsInt lc_row_out = row_out;
  factor.update(&col_aq, &row_ep, &lc_row_out, &rebuild_reason);
  if (rebuild_reason) return false;

  return testSolve();
}

bool testSolve() {
  const bool test_all = true;
  bool test_sparse_ftran = test_all;
  bool test_sparse_btran = test_all;
  bool test_dense_ftran = test_all;
  bool test_dense_btran = test_all;
  const HighsInt check_basis_change = -66;
  const bool debug_report = basis_change == check_basis_change;
  if (debug_report) {
    printf("basis_change = %d\n", (int)basis_change);
  }
  factor.setDebugReport(debug_report);
  double error_norm;
  HighsInt iCol;
  HighsRandom random(basis_change);
  iCol = random.integer(num_row);
  if (test_sparse_ftran) {
    // Sparse FTRAN
    rhs.clear();
    lp.a_matrix_.collectAj(rhs, basic_set[iCol], 1);
    factor.ftranCall(rhs, 1);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow == iCol) {
        error_norm = std::max(std::fabs(1 - rhs.array[iRow]), error_norm);
      } else {
        error_norm = std::max(std::fabs(rhs.array[iRow]), error_norm);
      }
    }
    if (debug_report)
      printf("Sparse FTRAN %2d: %g\n", (int)basis_change, error_norm);
    if (error_norm > 1e-4) return false;
  }

  iCol = random.integer(num_row);
  if (test_sparse_btran) {
    // Sparse BTRAN
    std::vector<double> unit;
    unit.assign(num_row, 0);
    unit[iCol] = 1.0;
    rhs.clear();
    for (HighsInt iCol = 0; iCol < num_row; iCol++) {
      rhs.array[iCol] = lp.a_matrix_.computeDot(unit, basic_set[iCol]);
      if (rhs.array[iCol]) rhs.index[rhs.count++] = iCol;
    }
    factor.btranCall(rhs, 1);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow == iCol) {
        error_norm = std::max(std::fabs(1 - rhs.array[iRow]), error_norm);
      } else {
        error_norm = std::max(std::fabs(rhs.array[iRow]), error_norm);
      }
    }
    if (debug_report)
      printf("Sparse BTRAN %2d: %g\n", (int)basis_change, error_norm);
    if (error_norm > 1e-4) return false;
  }

  if (test_dense_ftran) {
    // FTRAN
    rhs.clear();
    for (HighsInt iCol = 0; iCol < num_row; iCol++)
      lp.a_matrix_.collectAj(rhs, basic_set[iCol], solution[iCol]);
    factor.ftranCall(rhs, 1);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      error_norm =
          std::max(std::fabs(solution[iRow] - rhs.array[iRow]), error_norm);
    if (debug_report)
      printf("Dense  FTRAN %2d: %g\n", (int)basis_change, error_norm);
    if (error_norm > 1e-4) return false;
  }

  if (test_dense_btran) {
    // BTRAN
    rhs.clear();
    for (HighsInt iCol = 0; iCol < num_row; iCol++) {
      rhs.array[iCol] = lp.a_matrix_.computeDot(solution, basic_set[iCol]);
      if (rhs.array[iCol]) rhs.index[rhs.count++] = iCol;
    }
    factor.btranCall(rhs, 1);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      error_norm =
          std::max(std::fabs(solution[iRow] - rhs.array[iRow]), error_norm);
    if (debug_report)
      printf("Dense  BTRAN %2d: %g\n", (int)basis_change, error_norm);
    if (error_norm > 1e-4) return false;
  }
  factor.setDebugReport(false);
  return true;
}

bool testSolveDense() {
  const bool test_all = true;
  bool test_sparse_ftran = test_all;
  bool test_sparse_btran = test_all;
  bool test_dense_ftran = test_all;
  bool test_dense_btran = test_all;
  std::vector<double> rhs_dense;
  double error_norm = 0;
  HighsInt iCol;
  if (test_sparse_ftran) {
    HighsRandom random;
    iCol = random.integer(num_row);
    rhs.clear();
    lp.a_matrix_.collectAj(rhs, basic_set[iCol], 1);
    rhs_dense = rhs.array;
    factor.ftranCall(rhs_dense);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow == iCol) {
        error_norm = std::max(std::fabs(1 - rhs_dense[iRow]), error_norm);
      } else {
        error_norm = std::max(std::fabs(rhs_dense[iRow]), error_norm);
      }
    }
    if (dev_run) printf("Sparse  FTRAN: %g\n", error_norm);
    if (error_norm > 1e-4) return false;
  }

  if (test_dense_ftran) {
    rhs.clear();
    for (HighsInt iCol = 0; iCol < num_row; iCol++)
      lp.a_matrix_.collectAj(rhs, basic_set[iCol], solution[iCol]);
    rhs_dense = rhs.array;
    factor.ftranCall(rhs_dense);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      error_norm =
          std::max(std::fabs(solution[iRow] - rhs_dense[iRow]), error_norm);
    if (dev_run) printf("Dense   FTRAN: %g\n", error_norm);
    if (error_norm > 1e-4) return false;
  }

  if (test_sparse_btran) {
    // Sparse BTRAN
    std::vector<double> unit;
    unit.assign(num_row, 0);
    unit[iCol] = 1.0;
    rhs_dense.clear();
    for (HighsInt iCol = 0; iCol < num_row; iCol++)
      rhs_dense[iCol] = lp.a_matrix_.computeDot(unit, basic_set[iCol]);
    factor.btranCall(rhs_dense);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++) {
      if (iRow == iCol) {
        error_norm = std::max(std::fabs(1 - rhs_dense[iRow]), error_norm);
      } else {
        error_norm = std::max(std::fabs(rhs_dense[iRow]), error_norm);
      }
    }
    if (dev_run) printf("Sparse  BTRAN: %g\n", error_norm);
    if (error_norm > 1e-4) return false;
  }

  if (test_dense_btran) {
    // Dense BTRAN
    rhs_dense.clear();
    for (HighsInt iCol = 0; iCol < num_row; iCol++)
      rhs_dense[iCol] = lp.a_matrix_.computeDot(solution, basic_set[iCol]);
    factor.btranCall(rhs_dense);
    error_norm = 0;
    for (HighsInt iRow = 0; iRow < num_row; iRow++)
      error_norm =
          std::max(std::fabs(solution[iRow] - rhs_dense[iRow]), error_norm);
    if (dev_run) printf("Sparse  BTRAN: %g\n", error_norm);
    if (error_norm > 1e-4) return false;
  }

  return true;
}
