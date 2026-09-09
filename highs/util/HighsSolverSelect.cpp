/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSolverSelect.cpp
 * @brief
 */
#include "util/HighsSolverSelect.h"

#include <algorithm>
#include <cmath>

#include "lp_data/HConst.h"
#include "lp_data/HighsLp.h"

namespace {

bool isFinite(const double value) {
  return value > -kHighsInf && value < kHighsInf;
}

// Fraction of entries of `values` that share their value (exact equality) with
// at least one other entry.  `values` is sorted in place.
double relativeNumberOfDuplicates(std::vector<double>& values) {
  const size_t n = values.size();
  if (n < 2) return 0.0;
  std::sort(values.begin(), values.end());
  size_t num_in_group = 0;
  size_t run = 1;
  for (size_t i = 1; i <= n; i++) {
    if (i < n && values[i] == values[i - 1]) {
      run++;
    } else {
      if (run >= 2) num_in_group += run;
      run = 1;
    }
  }
  return static_cast<double>(num_in_group) / static_cast<double>(n);
}

// Fraction of `abs_values` that lie in a cluster of size >= 2, where two values
// join the same cluster when they agree to the relative tolerance `tol`.
// `abs_values` is sorted in place.
double relativeNumberOfAlmostIdentical(std::vector<double>& abs_values,
                                      const double tol) {
  const size_t n = abs_values.size();
  if (n < 2) return 0.0;
  std::sort(abs_values.begin(), abs_values.end());
  size_t num_in_cluster = 0;
  size_t run = 1;
  for (size_t i = 1; i <= n; i++) {
    const bool close =
        i < n && (abs_values[i] - abs_values[i - 1]) <=
                     tol * std::max(1.0, std::fabs(abs_values[i]));
    if (close) {
      run++;
    } else {
      if (run >= 2) num_in_cluster += run;
      run = 1;
    }
  }
  return static_cast<double>(num_in_cluster) / static_cast<double>(n);
}

double safeRatio(const double numerator, const double denominator) {
  if (denominator <= 0.0) return 0.0;
  return numerator / denominator;
}

double maxOverMin(const double max_abs, const double min_abs) {
  if (min_abs <= 0.0 || !isFinite(min_abs)) return 0.0;
  return max_abs / min_abs;
}

}  // namespace

HighsLpFeatures computeLpFeatures(const HighsLp& lp,
                                 const HighsLpFeatureParams& params) {
  HighsLpFeatures f;
  const HighsInt num_col = lp.num_col_;
  const HighsInt num_row = lp.num_row_;
  const HighsSparseMatrix& matrix = lp.a_matrix_;
  const HighsInt num_nz = matrix.numNz();

  f.num_row = num_row;
  f.num_col = num_col;
  f.num_nz = num_nz;

  const double d_num_col = static_cast<double>(num_col);
  const double d_num_row = static_cast<double>(num_row);
  const double d_num_nz = static_cast<double>(num_nz);

  // --------------------------------------------------------------------------
  // Column bounds and objective coefficients.
  // --------------------------------------------------------------------------
  HighsInt num_cols_without_upper = 0;
  HighsInt num_cols_without_lower = 0;
  HighsInt num_free_cols = 0;
  HighsInt num_boxed_cols = 0;
  HighsInt num_singly_bounded_cols = 0;
  HighsInt num_fixed_cols = 0;

  double max_abs_obj = 0.0;
  double min_abs_obj = kHighsInf;
  std::vector<double> obj_values;
  obj_values.reserve(num_col);

  for (HighsInt iCol = 0; iCol < num_col; iCol++) {
    const bool has_lower = lp.col_lower_[iCol] > -kHighsInf;
    const bool has_upper = lp.col_upper_[iCol] < kHighsInf;
    if (!has_upper) num_cols_without_upper++;
    if (!has_lower) num_cols_without_lower++;
    if (!has_lower && !has_upper) {
      num_free_cols++;
    } else if (has_lower && has_upper) {
      if (lp.col_lower_[iCol] == lp.col_upper_[iCol])
        num_fixed_cols++;
      else
        num_boxed_cols++;
    } else {
      num_singly_bounded_cols++;
    }

    const double cost = lp.col_cost_[iCol];
    obj_values.push_back(cost);
    const double abs_cost = std::fabs(cost);
    if (abs_cost > 0.0) {
      max_abs_obj = std::max(max_abs_obj, abs_cost);
      min_abs_obj = std::min(min_abs_obj, abs_cost);
    }
  }
  if (!lp.integrality_.empty()) {
    for (HighsInt iCol = 0; iCol < num_col; iCol++)
      if (lp.integrality_[iCol] == HighsVarType::kInteger ||
          lp.integrality_[iCol] == HighsVarType::kSemiContinuous ||
          lp.integrality_[iCol] == HighsVarType::kSemiInteger)
        f.num_integer_col++;
  }
  if (min_abs_obj == kHighsInf) min_abs_obj = 0.0;

  // --------------------------------------------------------------------------
  // Row activity bounds ("right-hand sides").
  // --------------------------------------------------------------------------
  HighsInt num_equalities = 0;
  HighsInt num_inequalities = 0;
  HighsInt num_ranged_rows = 0;
  HighsInt num_free_rows = 0;

  double max_abs_rhs = 0.0;
  double min_abs_rhs = kHighsInf;
  std::vector<double> rhs_values;        // one representative per row
  rhs_values.reserve(num_row);
  std::vector<double> all_finite_rhs;    // every finite bound, for the range
  all_finite_rhs.reserve(num_row);

  for (HighsInt iRow = 0; iRow < num_row; iRow++) {
    const double lower = lp.row_lower_[iRow];
    const double upper = lp.row_upper_[iRow];
    const bool has_lower = lower > -kHighsInf;
    const bool has_upper = upper < kHighsInf;
    if (!has_lower && !has_upper) {
      num_free_rows++;
    } else if (has_lower && has_upper) {
      if (lower == upper)
        num_equalities++;
      else
        num_ranged_rows++;
    } else {
      num_inequalities++;
    }

    // Representative rhs: finite upper if present, else finite lower.
    if (has_upper)
      rhs_values.push_back(upper);
    else if (has_lower)
      rhs_values.push_back(lower);

    for (int k = 0; k < 2; k++) {
      const double bound = k == 0 ? lower : upper;
      if (!isFinite(bound)) continue;
      all_finite_rhs.push_back(bound);
      const double abs_bound = std::fabs(bound);
      if (abs_bound > 0.0) {
        max_abs_rhs = std::max(max_abs_rhs, abs_bound);
        min_abs_rhs = std::min(min_abs_rhs, abs_bound);
      }
    }
  }
  if (min_abs_rhs == kHighsInf) min_abs_rhs = 0.0;

  // --------------------------------------------------------------------------
  // Single pass over the constraint matrix: per-row / per-column nonzero
  // counts and coefficient magnitude statistics.
  // --------------------------------------------------------------------------
  std::vector<HighsInt> row_count(num_row, 0);
  std::vector<HighsInt> col_count(num_col, 0);
  double max_abs_matrix = 0.0;
  double min_abs_matrix = kHighsInf;
  std::vector<double> abs_matrix_values;
  abs_matrix_values.reserve(num_nz);

  const bool colwise = matrix.isColwise();
  const HighsInt num_major = colwise ? num_col : num_row;
  const bool have_matrix =
      num_nz > 0 &&
      static_cast<HighsInt>(matrix.start_.size()) >= num_major + 1;
  if (have_matrix) {
    for (HighsInt major = 0; major < num_major; major++) {
      const HighsInt from = matrix.start_[major];
      const HighsInt to = matrix.start_[major + 1];
      if (colwise)
        col_count[major] = to - from;
      else
        row_count[major] = to - from;
      for (HighsInt el = from; el < to; el++) {
        const HighsInt minor = matrix.index_[el];
        if (colwise)
          row_count[minor]++;
        else
          col_count[minor]++;
        const double abs_value = std::fabs(matrix.value_[el]);
        abs_matrix_values.push_back(abs_value);
        if (abs_value > 0.0) {
          max_abs_matrix = std::max(max_abs_matrix, abs_value);
          min_abs_matrix = std::min(min_abs_matrix, abs_value);
        }
      }
    }
  }
  if (min_abs_matrix == kHighsInf) min_abs_matrix = 0.0;

  HighsInt max_row_count = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    max_row_count = std::max(max_row_count, row_count[iRow]);
  HighsInt max_col_count = 0;
  for (HighsInt iCol = 0; iCol < num_col; iCol++)
    max_col_count = std::max(max_col_count, col_count[iCol]);

  // "Dense" rows: longer than `dense_row_factor` times the mean row length and
  // at least `dense_row_min_count` nonzeros.
  const double mean_row_count = safeRatio(d_num_nz, d_num_row);
  const double dense_row_threshold = std::max(
      params.dense_row_factor * mean_row_count,
      static_cast<double>(params.dense_row_min_count));
  HighsInt num_dense_rows = 0;
  for (HighsInt iRow = 0; iRow < num_row; iRow++)
    if (static_cast<double>(row_count[iRow]) > dense_row_threshold)
      num_dense_rows++;

  // --------------------------------------------------------------------------
  // Assemble the feature record.
  // --------------------------------------------------------------------------
  // How constrained is the problem?
  f.relative_num_equalities = safeRatio(num_equalities, d_num_row);
  f.relative_num_inequalities = safeRatio(num_inequalities, d_num_row);
  f.relative_num_ranged_rows = safeRatio(num_ranged_rows, d_num_row);
  f.relative_num_free_rows = safeRatio(num_free_rows, d_num_row);
  f.relative_num_cols_without_upper = safeRatio(num_cols_without_upper, d_num_col);
  f.relative_num_cols_without_lower = safeRatio(num_cols_without_lower, d_num_col);
  f.relative_num_free_cols = safeRatio(num_free_cols, d_num_col);
  f.relative_num_boxed_cols = safeRatio(num_boxed_cols, d_num_col);
  f.relative_num_singly_bounded_cols =
      safeRatio(num_singly_bounded_cols, d_num_col);
  f.relative_num_fixed_cols = safeRatio(num_fixed_cols, d_num_col);

  // How degenerate is the problem?
  f.relative_num_almost_identical_nonzeros =
      relativeNumberOfAlmostIdentical(abs_matrix_values, params.almost_identical_tol);
  f.relative_num_identical_obj_values = relativeNumberOfDuplicates(obj_values);
  f.relative_num_identical_rhs_values = relativeNumberOfDuplicates(rhs_values);

  // How stable is the problem?
  f.max_abs_matrix_coefficient = max_abs_matrix;
  f.min_abs_matrix_coefficient = min_abs_matrix;
  f.max_abs_obj_coefficient = max_abs_obj;
  f.min_abs_obj_coefficient = min_abs_obj;
  f.max_abs_rhs_value = max_abs_rhs;
  f.min_abs_rhs_value = min_abs_rhs;
  f.relative_max_matrix_coefficient = maxOverMin(max_abs_matrix, min_abs_matrix);
  f.relative_max_obj_coefficient = maxOverMin(max_abs_obj, min_abs_obj);
  f.relative_max_rhs_value = maxOverMin(max_abs_rhs, min_abs_rhs);

  // How dense is the problem?
  f.matrix_nonzero_density = safeRatio(d_num_nz, d_num_row * d_num_col);
  f.avg_nonzeros_per_col = safeRatio(d_num_nz, d_num_col);
  f.avg_nonzeros_per_row = safeRatio(d_num_nz, d_num_row);
  f.relative_num_dense_rows = safeRatio(num_dense_rows, d_num_row);
  f.max_nonzeros_in_row = static_cast<double>(max_row_count);
  f.max_nonzeros_in_col = static_cast<double>(max_col_count);

  return f;
}

std::vector<std::pair<std::string, double>> highsLpFeatureVector(
    const HighsLpFeatures& f) {
  return {
      {"num_row", static_cast<double>(f.num_row)},
      {"num_col", static_cast<double>(f.num_col)},
      {"num_nz", static_cast<double>(f.num_nz)},
      {"num_integer_col", static_cast<double>(f.num_integer_col)},

      {"relative_num_equalities", f.relative_num_equalities},
      {"relative_num_inequalities", f.relative_num_inequalities},
      {"relative_num_ranged_rows", f.relative_num_ranged_rows},
      {"relative_num_free_rows", f.relative_num_free_rows},
      {"relative_num_cols_without_upper", f.relative_num_cols_without_upper},
      {"relative_num_cols_without_lower", f.relative_num_cols_without_lower},
      {"relative_num_free_cols", f.relative_num_free_cols},
      {"relative_num_boxed_cols", f.relative_num_boxed_cols},
      {"relative_num_singly_bounded_cols", f.relative_num_singly_bounded_cols},
      {"relative_num_fixed_cols", f.relative_num_fixed_cols},

      {"relative_num_almost_identical_nonzeros",
       f.relative_num_almost_identical_nonzeros},
      {"relative_num_identical_obj_values", f.relative_num_identical_obj_values},
      {"relative_num_identical_rhs_values", f.relative_num_identical_rhs_values},

      {"relative_max_matrix_coefficient", f.relative_max_matrix_coefficient},
      {"relative_max_obj_coefficient", f.relative_max_obj_coefficient},
      {"relative_max_rhs_value", f.relative_max_rhs_value},
      {"max_abs_matrix_coefficient", f.max_abs_matrix_coefficient},
      {"min_abs_matrix_coefficient", f.min_abs_matrix_coefficient},
      {"max_abs_obj_coefficient", f.max_abs_obj_coefficient},
      {"min_abs_obj_coefficient", f.min_abs_obj_coefficient},
      {"max_abs_rhs_value", f.max_abs_rhs_value},
      {"min_abs_rhs_value", f.min_abs_rhs_value},

      {"matrix_nonzero_density", f.matrix_nonzero_density},
      {"avg_nonzeros_per_col", f.avg_nonzeros_per_col},
      {"avg_nonzeros_per_row", f.avg_nonzeros_per_row},
      {"relative_num_dense_rows", f.relative_num_dense_rows},
      {"max_nonzeros_in_row", f.max_nonzeros_in_row},
      {"max_nonzeros_in_col", f.max_nonzeros_in_col},
  };
}

std::vector<std::string> highsLpFeatureNames() {
  std::vector<std::string> names;
  for (const auto& name_value : highsLpFeatureVector(HighsLpFeatures()))
    names.push_back(name_value.first);
  return names;
}

HighsSolverSelect selectSolverByFeatures(const HighsLpFeatures& f) {
  // NOTE: provisional heuristic.  Replace the body with the classifier fitted
  // offline on `highsLpFeatureVector` once the PCA / decision-tree study is
  // done; the feature definitions above are the contract with that study.
  if (f.num_nz == 0 || f.num_row == 0 || f.num_col == 0)
    return HighsSolverSelect::kDualSimplex;

  const bool large = f.num_nz > 500000;
  const bool very_large = f.num_nz > 5000000;
  const bool dense =
      f.matrix_nonzero_density > 1e-2 || f.avg_nonzeros_per_row > 25.0 ||
      f.relative_num_dense_rows > 1e-2;
  const bool ill_conditioned = f.relative_max_matrix_coefficient > 1e8 ||
                               f.relative_max_rhs_value > 1e8;

  // Simplex is the safe default for small, sparse, well-scaled models: it
  // thrives on sparsity and warm-starts cheaply.
  if (!large && !dense) return HighsSolverSelect::kDualSimplex;

  // Large and/or dense: an interior-point method amortises the factorisation
  // cost better.  Keep ill-conditioned models on IPX, whose linear algebra and
  // crossover are the more battle-tested; send the very largest well-scaled
  // models to HiPO.
  if (very_large && !ill_conditioned) return HighsSolverSelect::kHipo;
  return HighsSolverSelect::kIpx;
}

HighsSolverSelect selectSolverByFeatures(const HighsLp& lp,
                                         const HighsLpFeatureParams& params) {
  return selectSolverByFeatures(computeLpFeatures(lp, params));
}

HighsSolverSelect selectSolver(const HighsLp& lp) {
  return HighsSolverSelect::kDualSimplex;
}
