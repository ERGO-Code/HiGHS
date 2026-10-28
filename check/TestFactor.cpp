#include "HCheckConfig.h"
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
  highs.getBasicVariables(basic_set.data());
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

  highs.resetThreadScheduler(true);
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

  highs.resetThreadScheduler(true);
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

TEST_CASE("Factor-zero-col-pivot-tolerance", "[highs_test_factor]") {
  Highs highs;
  REQUIRE(highs.setOptionValue("output_flag", dev_run) == HighsStatus::kOk);
  REQUIRE(highs.setOptionValue("presolve", "off") == HighsStatus::kOk);

  // Near-dependent columns used to select a pivot below the absolute
  // tolerance, then discard a larger entry in that column.
  HighsLp model;
  model.num_col_ = 4;
  model.num_row_ = 5;
  model.col_cost_ = {0, 0, 0, 0};
  model.col_lower_ = {0, 0, 0, -kHighsInf};
  model.col_upper_ = {kHighsInf, kHighsInf, kHighsInf, kHighsInf};
  model.row_lower_ = {-kHighsInf, -kHighsInf, -kHighsInf, -kHighsInf,
                      -kHighsInf};
  model.row_upper_ = {0, 0, 0, 0, -1};
  model.a_matrix_.format_ = MatrixFormat::kColwise;
  model.a_matrix_.start_ = {0, 4, 8, 11, 14};
  model.a_matrix_.index_ = {0, 2, 3, 4, 0, 2, 3, 4, 0, 1, 4, 1, 2, 4};
  model.a_matrix_.value_ = {2, -1, 3e6, -1, -1e13, 1,    -3000001,
                            1, 1,  1,   -2, 1e7,   -1e7, 1};
  REQUIRE(highs.passModel(model) == HighsStatus::kOk);
  REQUIRE(highs.run() == HighsStatus::kOk);
  REQUIRE(highs.getModelStatus() == HighsModelStatus::kOptimal);

  highs.resetThreadScheduler(true);
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
      rhs_dense.push_back(lp.a_matrix_.computeDot(unit, basic_set[iCol]));
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
      rhs_dense.push_back(lp.a_matrix_.computeDot(solution, basic_set[iCol]));
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

TEST_CASE("study-factor", "[highs_test_factor]") {
  HighsInt dim = 3;
  HighsRandom random;
  double epsilon = 1e-12;
  double kPivotTolerance = 1e-8;
  const double kNumericalTroubleTolerance = 1e-7;

  HighsSparseMatrix matrix;
  matrix.num_col_ = 6;
  matrix.num_row_ = dim;
  matrix.start_ = {0, 2, 4, 6, 8, 10, 12};
  matrix.index_ = {0, 1, 1, 2, 0, 2, 0, 2, 0, 2, 0, 2};
  matrix.value_ = {1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1,
                   1 + epsilon * random.fraction(),
                   1 + epsilon * random.fraction(),
                   1,
                   -1};
  const HighsInt kParallelTo2 = dim;
  const HighsInt kAlmostParallelTo2 = kParallelTo2 + 1;
  const HighsInt kOrthogonalTo2 = kAlmostParallelTo2 + 1;
  HVector hat_a_q, pi_p, rhs;
  hat_a_q.setup(dim);
  pi_p.setup(dim);
  rhs.setup(dim);
  std::vector<HighsInt> basis = {0, 1, 2};
  HFactor factor;
  factor.setup(matrix, basis);

  HighsInt update_count;
  double update_solve_synthetic_tick;
  double invert_synthetic_tick;
  HighsInt illegal_hint = -kHighsIInf;
  HighsInt hint = illegal_hint;

  const bool logging = dev_run;  // true; //

  auto testInvert = [&]() {
    std::vector<double> solution;
    rhs.clear();
    for (HighsInt iCol = 0; iCol < dim; iCol++) {
      double value = random.fraction();
      solution.push_back(value);
      HighsInt iVar = basis[iCol];
      if (iVar < matrix.num_col_) {
        for (HighsInt iEl = matrix.start_[iVar]; iEl < matrix.start_[iVar + 1];
             iEl++)
          rhs.array[matrix.index_[iEl]] += matrix.value_[iEl] * value;
      } else {
        rhs.array[iVar - matrix.num_col_] += value;
      }
    }
    rhs.count = dim;
    factor.ftranCall(rhs, 1);
    double error = 0;
    for (HighsInt iCol = 0; iCol < dim; iCol++)
      error += std::fabs(rhs.array[iCol] - solution[iCol]);
    return error;
  };

  auto analyseFactor = [&]() {
    double error = testInvert();
    if (logging)
      printf(
          "After update %3d: INVERT error = %1.1e; Ticks: INVERT = %9d; UPDATE "
          "solve = %9d\n",
          int(update_count), error, int(invert_synthetic_tick),
          int(update_solve_synthetic_tick));
    return error < 1e-10;
  };

  auto printBasis = [&](const std::string& when) {
    if (!logging) return;
    printf("%s, basis is [", when.c_str());
    for (HighsInt iRow = 0; iRow < dim; iRow++)
      printf("%d%s", int(basis[iRow]), iRow < dim - 1 ? ", " : "]\n");
  };

  auto invert = [&]() {
    printBasis("\nBefore INVERT");
    HighsInt rank_deficiency = factor.build();
    update_count = 0;
    update_solve_synthetic_tick = 0;
    invert_synthetic_tick = factor.build_synthetic_tick;
    // If the basis matrix is singular, then the simplex solver
    // back-tracks to the basis half way through the previous sequence
    // of updates, and refactors. For now the QP ASM solver should log
    // an error and fail cleanly
    if (rank_deficiency && logging)
      printf("Basis matrix has rank deficiency of %d\n", int(rank_deficiency));
    printBasis("After  INVERT");
    if (logging) printf("INVERT error = %g\n\n", testInvert());
  };

  auto updateOk = [&](const HighsInt entering_column,
                      const HighsInt leaving_row) {
    hat_a_q.clear();
    matrix.collectAj(hat_a_q, entering_column, 1);
    hat_a_q.packFlag = true;
    factor.ftranCall(hat_a_q, 1);
    update_solve_synthetic_tick += hat_a_q.synthetic_tick;

    pi_p.clear();
    pi_p.count = 1;
    pi_p.index[0] = leaving_row;
    pi_p.array[leaving_row] = 1;
    pi_p.packFlag = true;
    factor.btranCall(pi_p, 1);
    update_solve_synthetic_tick += pi_p.synthetic_tick;

    // The pivot is e_p^T.B^{-1}.a_q, which can be computed as either
    //
    // pivot0 = e_p^T.(B^{-1}.a_q) = e_p^T.hat_a_q
    //
    // pivot1 = (e_p^T.B^{-1}).a_q = pi_p^T.a_q
    double pivot0 = hat_a_q.array[leaving_row];
    double pivot1;
    if (entering_column < matrix.num_col_) {
      pivot1 = 0;
      for (HighsInt iEl = matrix.start_[entering_column];
           iEl < matrix.start_[entering_column + 1]; iEl++)
        pivot1 += matrix.value_[iEl] * pi_p.array[matrix.index_[iEl]];
    } else {
      pivot1 = pi_p.array[entering_column - matrix.num_col_];
    }

    double numerical_trouble = std::fabs(pivot0 - pivot1) /
                               std::min(std::fabs(pivot0), std::fabs(pivot1));
    const bool bad_pivot = fabs(pivot0) < kPivotTolerance;
    const bool bad_numerics =
        bad_pivot || numerical_trouble > kNumericalTroubleTolerance;
    if ((bad_pivot || bad_numerics) && logging)
      printf(
          "%4s pivot; %4s numerics: pivot0 = %g; pivot1 = %g => "
          "numerical_trouble = %g\n",
          bad_pivot ? "Bad " : "Good", bad_numerics ? "Bad " : "Good", pivot0,
          pivot1, numerical_trouble);

    // Bad pivots should be avoided by tolerances in the ratio test,
    // but bad numerics can evolve - maybe as a consequence of a
    // sequence of poor (but not bad) pivots
    //
    // If bad numerics are detected after at least one update has
    // been performed, the simplex solver rejects the basis change
    // and reinverts the basis in the hope that the numerical
    // clean-up leads to better behaviour
    update_count++;
    // Because leaving_row is not const in HFactor::update
    HighsInt pass_leaving_row = leaving_row;
    factor.update(&hat_a_q, &pi_p, &pass_leaving_row, &hint);
    // hint is no longer used in HiGHS, since it no longer uses Qi's
    // experimental simplex update schemes that forced inversion
    // from within HFactor::update
    //
    // In the simplex solvers, decisions to refactor on numerical or
    // efficiency grounds are taken in the simplex solver based on
    // tests that are illustrated above in updateOk
    assert(hint == illegal_hint);
    basis[leaving_row] = entering_column;
    return analyseFactor();
  };

  auto replaceStructurals = [&]() {
    // Have to introduce logicals in this order to avoid singularity
    for (HighsInt iRow = dim - 1; iRow >= 0; iRow--)
      REQUIRE(updateOk(matrix.num_col_ + iRow, iRow));
  };

  auto replaceLogicals = [&]() {
    for (HighsInt iRow = 0; iRow < dim; iRow++) REQUIRE(updateOk(iRow, iRow));
  };

  invert();

  replaceStructurals();
  replaceLogicals();
  replaceStructurals();
  replaceLogicals();

  if (update_solve_synthetic_tick > invert_synthetic_tick) invert();

  std::vector<HighsInt> ok_basis = basis;
  // Column kParallelTo2 is the same as column 2, so replacing column
  // 0 (or 1) - basis entry 2 (or 1) by kParallelTo2 yields a singular
  // matrix
  HighsInt pivot_row = 2;
  assert(basis[pivot_row] != 2);
  REQUIRE(!updateOk(kParallelTo2, pivot_row));

  // Note that the zero pivot doesn't cause a NaN since the it
  // corresponds to a RHS entry of 0 (due to the basis matrix having
  // no entries in row 2) hence the zero pivot is skipped. There's a
  // meaningful INVERT error because the singular system has a
  // non-unique solution.

  // Factor the singular basis

  invert();
  // After INVERT, basis is [3, 1, 8], where the 8 is the identity
  // column variable replacing column 2 to ensure non-singularity
  //
  // If

  // Recover the non-singular basis [2, 1, 0]
  basis = ok_basis;
  invert();
  // After INVERT the basis is [0, 1, 2]

  // Column kAlmostParallelTo2 is almost the same as column 2, so
  // replacing column 0 (or 1) - basis entry 0 (or 1) by kParallelTo2
  // yields an ill-conditioned matrix
  pivot_row = 0;
  assert(basis[pivot_row] != 2);
  REQUIRE(!updateOk(kAlmostParallelTo2, pivot_row));

  REQUIRE(updateOk(kOrthogonalTo2, pivot_row));
}
