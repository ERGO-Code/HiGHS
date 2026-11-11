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
 * Regenerating the models in HighsSolverSelectModels.h
 *
 * The models are fitted offline by the scripts in select/ of the HiGHS
 * test repo:
 *
 * 1. presolve_set.sh writes the presolved model of each instance.
 * 2. run_solvers.sh <version> LP|MIP times dual simplex, IPX and HiPO on
 *    the presolved LPs and on the LP relaxations of presolved MIPs.
 * 3. extract_features/ applies LP presolve, as solveLp sees the LP, and
 *    writes computeLpFeatures() of each model to a CSV.
 * 4. collect_results.py --version <version> writes select.csv.
 * 5. visualize.py --cpp-out <HiGHS>/highs/util/HighsSolverSelectModels.h
 *    fits the gradient boosted trees and writes the header.
 *
 * The models index highsLpFeatureVector() by position: changing the order
 * of the features, or what computeLpFeatures() computes, needs a refit.
 * The test solver-select-feature-order fixes the order.
 */
#ifndef UTIL_HIGHS_SOLVER_SELECT_H_
#define UTIL_HIGHS_SOLVER_SELECT_H_

#include <string>
#include <utility>
#include <vector>

#include "lp_data/HConst.h"
#include "util/HighsInt.h"

class HighsLp;

enum class HighsSolverSelect { kDualSimplex, kIpx, kHipo };

enum class SolverSelectStrategy {
  kSimplexOnly = 0,
  kColsRowsNnz = 1,
  kRankedFeatures = 2,
  kHeuristic = 3,
  // Always choose one solver, to compare them
  kDualSimplex = 4,
  kIpx = 5,
  kHipo = 6,
};

// Whether LPs go through solver selection with this strategy
inline bool useSolverSelection(const HighsInt strategy) {
  return strategy > static_cast<HighsInt>(SolverSelectStrategy::kSimplexOnly);
}

// Whether solver selection with this strategy may choose IPM. Presolve
// must then allow for dual postsolve
inline bool solverSelectionMayChooseIpm(const HighsInt strategy) {
  return useSolverSelection(strategy) &&
         strategy != static_cast<HighsInt>(SolverSelectStrategy::kDualSimplex);
}

struct HighsLpFeatures {
  // --- Raw size ---
  HighsInt num_row = 0;
  HighsInt num_col = 0;
  HighsInt num_nz = 0;
  HighsInt num_integer_col = 0;  // Set only by offline feature extraction

  // --- How constrained is the problem? (barrier likes an interior) ---
  double relative_num_equalities = 0.0;           // equality rows / #rows
  double relative_num_inequalities = 0.0;         // one-sided rows / #rows
  double relative_num_ranged_rows = 0.0;          // l < u, both finite
  double relative_num_free_rows = 0.0;            // -inf..inf rows / #rows
  double relative_num_cols_without_upper = 0.0;   // u_j = +inf / #cols
  double relative_num_cols_without_lower = 0.0;   // l_j = -inf / #cols
  double relative_num_free_cols = 0.0;            // -inf..inf cols / #cols
  double relative_num_boxed_cols = 0.0;           // both bounds finite
  double relative_num_singly_bounded_cols = 0.0;  // exactly one finite
  double relative_num_fixed_cols = 0.0;           // l_j == u_j

  // --- How degenerate is the problem? (simplex can get stuck) ---
  double relative_num_almost_identical_nonzeros = 0.0;
  double relative_num_identical_obj_values = 0.0;  // repeated c_j
  double relative_num_identical_rhs_values = 0.0;  // repeated rhs

  // --- How stable is the problem? (conditioning) ---
  double relative_max_matrix_coefficient = 0.0;  // max|a| / min|a|
  double relative_max_obj_coefficient = 0.0;     // max|c| / min|c|
  double relative_max_rhs_value = 0.0;           // max|rhs| / min|rhs|
  double max_abs_matrix_coefficient = 0.0;
  double min_abs_matrix_coefficient = 0.0;
  double max_abs_obj_coefficient = 0.0;
  double min_abs_obj_coefficient = 0.0;
  double max_abs_rhs_value = 0.0;
  double min_abs_rhs_value = 0.0;

  // --- How dense is the problem? (simplex thrives on sparsity) ---
  double matrix_nonzero_density = 0.0;   // #nz / (#rows * #cols)
  double avg_nonzeros_per_col = 0.0;     // #nz / #cols
  double avg_nonzeros_per_row = 0.0;     // #nz / #rows
  double relative_num_dense_rows = 0.0;  // "dense" rows / #rows
  double max_nonzeros_in_row = 0.0;
  double max_nonzeros_in_col = 0.0;
};

HighsLpFeatures computeLpFeatures(const HighsLp& lp);

// Flatten `features` into (name, value) pairs, in the stable order that
// the fitted models index.  `highsLpFeatureNames()` returns just the names
// (same order) so a header row can be emitted without a features instance.
std::vector<std::pair<std::string, double>> highsLpFeatureVector(
    const HighsLpFeatures& features);
std::vector<std::string> highsLpFeatureNames();

// kSimplexOnly => dual simplex; kColsRowsNnz => boosted trees
// on num_col, num_row and num_nz; kRankedFeatures => boosted trees
// on the ranked features chosen when fitting (see the top of this file);
// kHeuristic => provisional hand-written heuristic; kDualSimplex, kIpx
// and kHipo => that solver
HighsSolverSelect selectSolver(const HighsLp& lp,
                               const SolverSelectStrategy strategy);

// kColsRowsNnz and kRankedFeatures: evaluate the gradient boosted
// trees in HighsSolverSelectModels.h, fitted offline as described at the
// top of this file
HighsSolverSelect selectSolverByModel(const HighsLpFeatures& features,
                                      const SolverSelectStrategy strategy);

// kHeuristic: the provisional hand-written heuristic
HighsSolverSelect selectSolverByHeuristic(const HighsLpFeatures& features);

#endif /* UTIL_HIGHS_SOLVER_SELECT_H_ */
