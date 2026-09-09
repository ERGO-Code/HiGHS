/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsSolverSelect.h
 * @brief Choose an LP algorithm (simplex / IPX / HiPO) from cheap
 *        structural features of the (presolved) model.
 *
 * The feature set mirrors the one described by Zonghao Gu in the Gurobi 11.0
 * "New Performance Techniques" talk (slide 6): a handful of scale-free
 * ("relative") numbers computed from the constraint matrix, the objective and
 * the bounds, grouped by the four questions that matter for algorithm choice:
 *
 *   - How constrained is the problem?  Barrier likes having an interior.
 *   - How degenerate is the problem?   Simplex can get stuck.
 *   - How stable is the problem?       Simplex and barrier suffer differently.
 *   - How dense is the problem?        Simplex thrives on sparsity.
 *
 * `computeLpFeatures` fills a `HighsLpFeatures` from a `HighsLp`; it is meant to
 * be run over a corpus of presolved instances and dumped to CSV
 * (`highsLpFeatureVector`) so that a PCA / decision model can be fitted
 * offline.  `selectSolverByFeatures` applies a (currently heuristic)
 * classifier to those features for a single new instance.
 */
#ifndef UTIL_HIGHS_SOLVER_SELECT_H_
#define UTIL_HIGHS_SOLVER_SELECT_H_

#include <string>
#include <utility>
#include <vector>

#include "util/HighsInt.h"

class HighsLp;

enum class HighsSolverSelect {
    kDualSimplex,
    kPrimalSimplex,
    kIpx,
    kHipo,
    kCupdlp,
    kHipdlp
};

// Tuning constants for the feature computation.  Collected here so that the
// offline study and the online selector use exactly the same definitions.
struct HighsLpFeatureParams {
  // Two coefficients are "almost identical" when their absolute values agree to
  // this relative tolerance.  Feeds `relative_num_almost_identical_nonzeros`.
  double almost_identical_tol = 1e-9;
  // A row is "dense" when its nonzero count exceeds this multiple of the mean
  // row length, and is also at least `dense_row_min_count` long.  Feeds
  // `relative_num_dense_rows`.
  double dense_row_factor = 10.0;
  HighsInt dense_row_min_count = 10;
};

// Structural features of a single (presolved) LP.  Every "relative_*" member is
// scale free: a count divided by #rows, #cols or #nonzeros, or a max-over-min
// magnitude ratio.  The plain counts / extrema are kept too as raw material for
// deriving further ratios offline.
//
// The numbered comments refer to the slide-6 "features by sorted predictive
// power" list.
struct HighsLpFeatures {
  // --- Raw size --------------------------------------------------------------
  HighsInt num_row = 0;
  HighsInt num_col = 0;
  HighsInt num_nz = 0;
  HighsInt num_integer_col = 0;

  // --- How constrained is the problem? (barrier likes an interior) ----------
  double relative_num_equalities = 0.0;            //  1. equality rows / #rows
  double relative_num_inequalities = 0.0;          //     one-sided rows / #rows
  double relative_num_ranged_rows = 0.0;           //     l < u, both finite
  double relative_num_free_rows = 0.0;             //     -inf..inf rows / #rows
  double relative_num_cols_without_upper = 0.0;    // 12. u_j = +inf / #cols
  double relative_num_cols_without_lower = 0.0;    //     l_j = -inf / #cols
  double relative_num_free_cols = 0.0;             //     -inf..inf cols / #cols
  double relative_num_boxed_cols = 0.0;            //     both bounds finite
  double relative_num_singly_bounded_cols = 0.0;   //     exactly one finite
  double relative_num_fixed_cols = 0.0;            //     l_j == u_j

  // --- How degenerate is the problem? (simplex can get stuck) ---------------
  double relative_num_almost_identical_nonzeros = 0.0;  //  2. see params
  double relative_num_identical_obj_values = 0.0;       // 10. repeated c_j
  double relative_num_identical_rhs_values = 0.0;       // 11. repeated rhs

  // --- How stable is the problem? (conditioning) ---------------------------
  double relative_max_matrix_coefficient = 0.0;   //  3. max|a| / min|a|
  double relative_max_obj_coefficient = 0.0;      //  4. max|c| / min|c|
  double relative_max_rhs_value = 0.0;            //  6. max|rhs| / min|rhs|
  double max_abs_matrix_coefficient = 0.0;
  double min_abs_matrix_coefficient = 0.0;
  double max_abs_obj_coefficient = 0.0;
  double min_abs_obj_coefficient = 0.0;
  double max_abs_rhs_value = 0.0;
  double min_abs_rhs_value = 0.0;

  // --- How dense is the problem? (simplex thrives on sparsity) --------------
  double matrix_nonzero_density = 0.0;   //  5. #nz / (#rows * #cols)
  double avg_nonzeros_per_col = 0.0;     //  7. #nz / #cols
  double avg_nonzeros_per_row = 0.0;     //  9. #nz / #rows
  double relative_num_dense_rows = 0.0;  //  8. "dense" rows / #rows
  double max_nonzeros_in_row = 0.0;
  double max_nonzeros_in_col = 0.0;

  void clear() { *this = HighsLpFeatures(); }
};

// Compute the structural features of `lp`.  `lp` is not modified; the constraint
// matrix is read in whatever orientation it is already stored in.
HighsLpFeatures computeLpFeatures(
    const HighsLp& lp, const HighsLpFeatureParams& params = HighsLpFeatureParams());

// Flatten `features` into (name, value) pairs, in a stable order, for writing a
// feature matrix to CSV.  `highsLpFeatureNames()` returns just the names (same
// order) so a header row can be emitted without a features instance.
std::vector<std::pair<std::string, double>> highsLpFeatureVector(
    const HighsLpFeatures& features);
std::vector<std::string> highsLpFeatureNames();

// Original entry point: pick a solver for `lp` (expected to be the presolved
// LP).  Currently a stub returning dual simplex.
HighsSolverSelect selectSolver(const HighsLp& lp);

// Pick a solver from the slide-6 features.  Until the PCA / decision model from
// the offline study is wired in, this applies a transparent heuristic over
// `computeLpFeatures(lp)` and only ever returns kDualSimplex, kIpx or kHipo.
HighsSolverSelect selectSolverByFeatures(
    const HighsLp& lp, const HighsLpFeatureParams& params = HighsLpFeatureParams());
HighsSolverSelect selectSolverByFeatures(const HighsLpFeatures& features);

#endif /* UTIL_HIGHS_SOLVER_SELECT_H_ */
