#include <algorithm>

//#include "HConfig.h"
#include "Highs.h"
#include "HighsRandom.h"
#include "catch.hpp"

using std::min;

const bool dev_run = false;

bool GetBasisSolvesSolutionNzOk(HighsInt numRow,
                                const vector<double>& pass_solution_vector,
                                HighsInt* solution_num_nz,
                                vector<HighsInt>& solution_indices) {
  if (solution_num_nz == NULL) return true;
  vector<double> solution_vector;
  solution_vector.resize(numRow);
  bool solution_nz_ok = true;
  for (HighsInt row = 0; row < numRow; row++)
    solution_vector[row] = pass_solution_vector[row];
  // Check that the indexed entries are nonzero
  for (HighsInt ix = 0; ix < *solution_num_nz; ix++) {
    HighsInt row = solution_indices[ix];
    if (!solution_vector[row]) {
      if (dev_run)
        printf("SolutionNzOk: Indexed entry solution_vector[%2" HIGHSINT_FORMAT
               "] = %11.4g\n",
               row, solution_vector[row]);
      solution_nz_ok = false;
    } else {
      solution_vector[row] = 0;
    }
  }
  // Solution should now be zero
  for (HighsInt row = 0; row < numRow; row++) {
    if (solution_vector[row]) {
      if (dev_run)
        printf(
            "SolutionNzOk: Non-indexed entry solution_vector[%2" HIGHSINT_FORMAT
            "] = %11.4g\n",
            row, solution_vector[row]);
      solution_nz_ok = false;
    }
  }
  return solution_nz_ok;
}

double GetBasisSolvesCheckSolution(const HighsLp& lp,
                                   const vector<HighsInt>& basic_variables,
                                   const vector<double>& rhs,
                                   const vector<double>& solution,
                                   const bool transpose = false) {
  const double residual_tolerance = 1e-8;
  double residual_norm = 0;
  if (transpose) {
    for (HighsInt k = 0; k < lp.num_row_; k++) {
      double residual = 0;
      HighsInt var = basic_variables[k];
      if (var < 0) {
        HighsInt row = -(1 + var);
        residual = fabs(rhs[k] - solution[row]);
        if (residual > residual_tolerance) {
          if (dev_run)
            printf("Row |[B^Tx-b]_{%2" HIGHSINT_FORMAT "}| = %11.4g\n", k,
                   residual);
        }
      } else {
        HighsInt col = var;
        for (HighsInt el = lp.a_matrix_.start_[col];
             el < lp.a_matrix_.start_[col + 1]; el++) {
          HighsInt row = lp.a_matrix_.index_[el];
          residual += lp.a_matrix_.value_[el] * solution[row];
        }
        residual = fabs(rhs[k] - residual);
        if (residual > residual_tolerance) {
          if (dev_run)
            printf("Col |[B^Tx-b]_{%2" HIGHSINT_FORMAT "}| = %11.4g\n", k,
                   residual);
        }
      }
      residual_norm += residual;
    }
  } else {
    vector<double> basis_matrix_times_solution;
    basis_matrix_times_solution.assign(lp.num_row_, 0);
    for (HighsInt k = 0; k < lp.num_row_; k++) {
      HighsInt var = basic_variables[k];
      if (var < 0) {
        HighsInt row = -(1 + var);
        basis_matrix_times_solution[row] += solution[k];
      } else {
        HighsInt col = var;
        for (HighsInt el = lp.a_matrix_.start_[col];
             el < lp.a_matrix_.start_[col + 1]; el++) {
          HighsInt row = lp.a_matrix_.index_[el];
          basis_matrix_times_solution[row] +=
              lp.a_matrix_.value_[el] * solution[k];
        }
      }
    }
    for (HighsInt k = 0; k < lp.num_row_; k++) {
      double residual = fabs(rhs[k] - basis_matrix_times_solution[k]);
      if (residual > residual_tolerance) {
        if (dev_run)
          printf("|[B^Tx-b]_{%2" HIGHSINT_FORMAT "}| = %11.4g\n", k, residual);
      }
      residual_norm += residual;
    }
  }
  return residual_norm;
}

void GetBasisSolvesFormRHS(HighsLp& lp, vector<HighsInt>& basic_variables,
                           vector<double>& solution, vector<double>& rhs,
                           const bool transpose = false) {
  if (transpose) {
    for (HighsInt k = 0; k < lp.num_row_; k++) {
      rhs[k] = 0;
      HighsInt var = basic_variables[k];
      if (var < 0) {
        HighsInt row = -(1 + var);
        rhs[k] = solution[row];
      } else {
        HighsInt col = var;
        for (HighsInt el = lp.a_matrix_.start_[col];
             el < lp.a_matrix_.start_[col + 1]; el++) {
          HighsInt row = lp.a_matrix_.index_[el];
          rhs[k] += lp.a_matrix_.value_[el] * solution[row];
        }
      }
    }
  } else {
    for (HighsInt k = 0; k < lp.num_row_; k++) rhs[k] = 0;
    for (HighsInt k = 0; k < lp.num_row_; k++) {
      HighsInt var = basic_variables[k];
      if (var < 0) {
        HighsInt row = -(1 + var);
        rhs[row] += solution[k];
      } else {
        HighsInt col = var;
        for (HighsInt el = lp.a_matrix_.start_[col];
             el < lp.a_matrix_.start_[col + 1]; el++) {
          HighsInt row = lp.a_matrix_.index_[el];
          rhs[row] += lp.a_matrix_.value_[el] * solution[k];
        }
      }
    }
  }
}

void testBasisSolve(Highs& highs) {
  HighsStatus highs_status;

  vector<HighsInt> basic_variables, solution_row_indices, solution_col_indices;
  vector<double> rhs, known_solution, solution_row, solution_col;

  HighsLp lp = highs.getLp();
  HighsInt numRow = lp.num_row_;
  HighsInt numCol = lp.num_col_;
  basic_variables.resize(numRow);
  known_solution.resize(numRow);
  rhs.resize(numRow);
  solution_row.resize(numCol);
  solution_col.resize(numRow);
  solution_row_indices.resize(numCol);
  solution_col_indices.resize(numRow);

  HighsInt solution_num_nz;

  HighsInt check_row = 0;
  HighsInt check_col = 0;

  double residual_norm;
  const double residual_norm_tolerance = 1e-8;
  const double solution_error_tolerance = 1e-8;
  HighsRandom random;

  HighsInt basic_col;

  highs_status = highs.getBasicVariables(&basic_variables[0]);
  REQUIRE(highs_status == HighsStatus::kOk);

  for (HighsInt ix = 0; ix < numRow; ix++) known_solution[ix] = 0;
  bool transpose = true;
  HighsInt num_ix = 3;
  HighsInt col;
  col = 6;
  basic_col = basic_variables[col];
  known_solution[col] = random.fraction();

  if (num_ix > 1) {
    col = 15;
    basic_col = basic_variables[col];
    known_solution[col] = random.fraction();
  }

  if (num_ix > 2) {
    col = 12;
    basic_col = basic_variables[col];
    known_solution[col] = random.fraction();
  }

  GetBasisSolvesFormRHS(lp, basic_variables, known_solution, rhs, transpose);
  if (transpose) {
    highs_status = highs.getBasisTransposeSolve(&rhs[0], &solution_col[0]);
  } else {
    highs_status = highs.getBasisSolve(&rhs[0], &solution_col[0]);
  }
  REQUIRE(highs_status == HighsStatus::kOk);
  residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                              solution_col, transpose);
  REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
  double solution_error_norm = 0;
  for (HighsInt ix = 0; ix < numRow; ix++) {
    double solution_error = fabs(known_solution[ix] - solution_col[ix]);
    if (solution_error > solution_error_tolerance) {
      if (dev_run)
        printf("Row %2" HIGHSINT_FORMAT ": |x-x^|_i = %11.4g\n", ix,
               solution_error);
      solution_error_norm += solution_error;
    }
  }
  if (dev_run)
    printf(
        "Test 0:     residual_norm = %11.4g\n      solution_error_norm = "
        "%11.4g "
        "(Known solution)\n",
        residual_norm, solution_error_norm);

  double max_residual_norm;
  HighsInt max_k = min(numRow, HighsInt{9});
  HighsInt k;

  k = 0;
  max_residual_norm = 0;
  for (HighsInt row = 0; row < numRow; row++) {
    HighsInt var = basic_variables[row];
    if (var >= 0) {
      basic_col = var;
      for (HighsInt ix = 0; ix < numRow; ix++) rhs[ix] = 0;
      for (HighsInt el = lp.a_matrix_.start_[basic_col];
           el < lp.a_matrix_.start_[basic_col + 1]; el++)
        rhs[lp.a_matrix_.index_[el]] = lp.a_matrix_.value_[el];

      highs_status =
          highs.getBasisSolve(&rhs[0], &solution_col[0], &solution_num_nz,
                              &solution_col_indices[0]);
      REQUIRE(highs_status == HighsStatus::kOk);
      bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
          numRow, solution_col, &solution_num_nz, solution_col_indices);
      REQUIRE(solution_nz_ok == true);
      residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                  solution_col, false);
      max_residual_norm = std::max(residual_norm, max_residual_norm);
      if (residual_norm > residual_norm_tolerance && dev_run)
        printf("getBasisSolve(%" HIGHSINT_FORMAT "): residual_norm = %g\n", k,
               residual_norm);
      REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
      if (k < max_k)
        k++;
      else
        k *= 2;
    }
    if (k >= numRow) break;
  }
  if (dev_run)
    printf("Test 1: max_residual_norm = %11.4g (Basic column)\n",
           max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    check_row = k;
    // Determine row check_row of B^{-1}
    highs_status =
        highs.getBasisInverseRow(check_row, &solution_col[0], &solution_num_nz,
                                 &solution_col_indices[0]);
    REQUIRE(highs_status == HighsStatus::kOk);
    bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
        numRow, solution_col, &solution_num_nz, solution_col_indices);
    REQUIRE(solution_nz_ok == true);
    // Check solution
    // Set up RHS as e_{check_row}
    for (HighsInt row = 0; row < numRow; row++) rhs[row] = 0;
    rhs[check_row] = 1;
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_col, true);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance && dev_run)
      printf("getBasisInverseRow(%" HIGHSINT_FORMAT "): residual_norm = %g\n",
             k, residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  if (dev_run)
    printf("Test 2: max_residual_norm = %11.4g (getBasisInverseRow)\n",
           max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    check_col = k;
    // Determine col check_col of B^{-1}
    highs_status =
        highs.getBasisInverseCol(check_col, &solution_col[0], &solution_num_nz,
                                 &solution_col_indices[0]);
    REQUIRE(highs_status == HighsStatus::kOk);
    bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
        numRow, solution_col, &solution_num_nz, solution_col_indices);
    REQUIRE(solution_nz_ok == true);
    // Check solution
    // Set up RHS as e_{check_col}
    for (HighsInt row = 0; row < numRow; row++) rhs[row] = 0;
    rhs[check_col] = 1;
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_col, false);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance && dev_run)
      printf("getBasisInverseCol(%" HIGHSINT_FORMAT "): residual_norm = %g\n",
             k, residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  if (dev_run)
    printf("Test 3: max_residual_norm = %11.4g (getBasisInverseCol)\n",
           max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    for (HighsInt row = 0; row < numRow; row++) rhs[row] = random.fraction();
    highs_status = highs.getBasisSolve(&rhs[0], &solution_col[0]);
    REQUIRE(highs_status == HighsStatus::kOk);
    // Check solution
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_col, false);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance && dev_run)
      printf("getBasisSolve(%" HIGHSINT_FORMAT "): residual_norm = %g\n", k,
             residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  if (dev_run)
    printf("Test 4: max_residual_norm = %11.4g (getBasisSolve)\n",
           max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  for (;;) {
    for (HighsInt row = 0; row < numRow; row++) rhs[row] = random.fraction();
    highs_status = highs.getBasisTransposeSolve(&rhs[0], &solution_col[0]);
    REQUIRE(highs_status == HighsStatus::kOk);
    // Check solution
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_col, true);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance && dev_run)
      printf("getBasisTransposeSolve(%" HIGHSINT_FORMAT
             "): residual_norm = %g\n",
             k, residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  if (dev_run)
    printf("Test 5: max_residual_norm = %11.4g (getBasisTransposeSolve)\n",
           max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  max_k = min(numRow, HighsInt{9});
  for (;;) {
    check_row = k;
    highs_status =
        highs.getReducedRow(check_row, &solution_row[0], &solution_num_nz,
                            &solution_row_indices[0]);
    REQUIRE(highs_status == HighsStatus::kOk);
    bool solution_nz_ok = GetBasisSolvesSolutionNzOk(
        numCol, solution_row, &solution_num_nz, solution_row_indices);
    REQUIRE(solution_nz_ok == true);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numRow) break;
  }
  if (dev_run)
    printf("Test 6: max_residual_norm = %11.4g (getReducedRow)\n",
           max_residual_norm);

  k = 0;
  max_residual_norm = 0;
  max_k = min(numCol, HighsInt{9});
  for (;;) {
    check_col = k;
    highs_status =
        highs.getReducedColumn(check_col, &solution_col[0], &solution_num_nz,
                               &solution_col_indices[0]);
    REQUIRE(highs_status == HighsStatus::kOk);
    // Check solution
    for (HighsInt row = 0; row < numRow; row++) rhs[row] = 0;
    for (HighsInt el = lp.a_matrix_.start_[check_col];
         el < lp.a_matrix_.start_[check_col + 1]; el++)
      rhs[lp.a_matrix_.index_[el]] = lp.a_matrix_.value_[el];
    residual_norm = GetBasisSolvesCheckSolution(lp, basic_variables, rhs,
                                                solution_col, false);
    max_residual_norm = std::max(residual_norm, max_residual_norm);
    if (residual_norm > residual_norm_tolerance && dev_run)
      printf("getBasisTransposeSolve(%" HIGHSINT_FORMAT
             "): residual_norm = %g\n",
             k, residual_norm);
    REQUIRE(fabs(residual_norm) < residual_norm_tolerance);
    if (k < max_k)
      k++;
    else
      k *= 2;
    if (k >= numCol) break;
  }
  if (dev_run)
    printf("Test 7: max_residual_norm = %11.4g (getReducedColumn)\n",
           max_residual_norm);
}

// No commas in test case name.
TEST_CASE("Basis-solves", "[highs_basis_solves]") {
  std::string filename;
  filename = std::string(HIGHS_DIR) + "/check/instances/chip.mps";
  filename = std::string(HIGHS_DIR) + "/check/instances/avgas.mps";
  filename = std::string(HIGHS_DIR) + "/check/instances/adlittle.mps";
  //  filename = std::string(HIGHS_DIR) + "/check/instances/25fv47.mps";

  Highs highs;
  if (!dev_run) highs.setOptionValue("output_flag", false);

  vector<HighsInt> basic_variables;
  vector<double> rhs, solution_row, solution_col;
  basic_variables.resize(1);
  rhs.resize(1);
  solution_row.resize(1);
  solution_col.resize(1);

  HighsStatus highs_status;

  // Read the LP given by filename
  highs_status = highs.readModel(filename);
  REQUIRE(highs_status == HighsStatus::kOk);

  highs_status = highs.writeModel("");
  REQUIRE(highs_status == HighsStatus::kOk);

  HighsInt numRow = highs.getNumRow();
  HighsInt numCol = highs.getNumCol();
  basic_variables.resize(numRow);
  rhs.resize(numRow);
  solution_row.resize(numCol);
  solution_col.resize(numRow);

  // Check the NULL pointer trap for RHS
  highs_status = highs.getBasisSolve(NULL, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisSolve(NULL, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  // Check the NULL pointer trap for the basic variables
  highs_status = highs.getBasicVariables(NULL);
  REQUIRE(highs_status == HighsStatus::kError);

  // Check the NULL pointer trap for the solution vector
  highs_status = highs.getBasisInverseRow(0, NULL);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisInverseCol(0, NULL);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisSolve(&rhs[0], NULL);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisTransposeSolve(&rhs[0], NULL);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getReducedRow(0, NULL);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getReducedColumn(0, NULL);
  REQUIRE(highs_status == HighsStatus::kError);

  // Check the indexing traps
  highs_status = highs.getBasisInverseRow(-1, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisInverseRow(numRow, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisInverseCol(-1, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisInverseCol(numCol, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  // Check the no INVERSE traps - these should all work, as the first should
  // force inversion!!!
  highs_status = highs.getBasicVariables(&basic_variables[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisInverseRow(0, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisInverseCol(0, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisSolve(&rhs[0], &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getBasisTransposeSolve(&rhs[0], &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getReducedRow(0, &solution_row[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  highs_status = highs.getReducedColumn(0, &solution_col[0]);
  REQUIRE(highs_status == HighsStatus::kError);

  // Solve and perform numerical tests
  highs_status = highs.run();
  REQUIRE(highs_status == HighsStatus::kOk);

  testBasisSolve(highs);

  // Save the optimal basis and clear the model
  HighsBasis basis = highs.getBasis();

  highs.clearModel();

  // Read the LP given by filename
  highs.readModel(filename);

  // Load the optimal basis
  highs_status = highs.setBasis(basis);
  REQUIRE(highs_status == HighsStatus::kOk);

  testBasisSolve(highs);

  // Solve
  highs.run();
  REQUIRE(highs.getInfo().simplex_iteration_count == 0);
}
