# Solver selection by LP features

Pick an LP algorithm — dual simplex, IPX or HiPO — from cheap structural
features of the presolved model, and extract those features in bulk so a
selection model can be fitted.

## Idea

Features of the presolved model, fitted a decision tree to choose between simplex and IPM. Features are scale free: nrows, ncols, nonzeros, or a max/min ratio, so they
compare across instances of any size. Four questions:


| Question | Reason |
|----------|--------|
| How constrained? | Barrier needs a non-empty interior |
| How degenerate?  | Simplex can stall |
| How stable?      | Wide value ranges hurt simplex and barrier differently |
| How dense?       | Simplex thrives on sparsity |

## Workflow

1. Extract — for each presolved instance, write `highsLpFeatureVector` as a
   CSV row (`tests/select/extract_features`, in the `tests` repo —
   see [Batch tool](#batch-tool) below).
2. Label — `tests/select/collect_results.py` collects the 1.15.1 dual
   simplex / primal simplex / IPX / HiPO benchmark runs into `select.csv`,
   checking from the logs that the requested solver ran and how crossover
   ended, so that each solve is flagged as valid with and without a basis
   being required.
3. Fit — `tests/select/visualize.py` standardises the (log-scaled)
   features, projects them onto the leading principal components, and fits
   a shallow "regret tree" on the PC scores: each leaf picks the solver with
   the smallest mean log-slowdown relative to the fastest. PCs kept and tree
   depth are chosen by cross validation. With `--cpp-out` it writes the four
   fitted models (strategy 1 / 2 × basis required or not) to
   `HighsSolverSelectModels.h`, which `selectSolverByModel` evaluates.
   Without a basis requirement, PDLP is a proxy class for instances that no
   other solver solved within the time limit (there are no PDLP timings).

## Strategies — `solver_select_strategy`

`solveLp` dispatches to `selectSolver(lp, strategy, require_basis,
num_integer_col)` whenever `solver = choose` (or, inside the MIP solver,
`mip_lp_solver = choose`). `strategy` is a `SolverSelectStrategy` (enum
class) in C++; the **`solver_select_strategy`** *option* that drives it is
the plain `HighsInt` the enum's values match — range `[0, 3]`, default `0` —
settable via `Highs::setOptionValue` or the CLI flag
`--solver_select_strategy`. The one place the raw option int becomes the
typed enum is the `solveLp` call site, via `static_cast<SolverSelectStrategy>`.

| Value | `SolverSelectStrategy` | Behaviour |
|:-----:|-------------------------|-----------|
| 0 | `kSimplexOnly` | Always `kDualSimplex`. |
| 1 | `kNColsRowsNnzIntegrality` | `selectSolverByModel`: the fitted PCA + decision-tree model on 4 size features (`num_col`, `num_row`, `num_nz`, `num_integer_col`) — cheap, skips the rest of `computeLpFeatures`. |
| 2 | `kFull` | `selectSolverByModel` on the full feature set (below). |
| 3 | `kHeuristic` | `selectSolverByFeatures`, the hand-written heuristic (see [below](#strategy-3--hand-written-heuristic)). |

`require_basis` is the **`solver_select_require_basis`** option — bool,
default `false`, CLI flag `--solver_select_require_basis` — which restricts
strategies 1–3 to a solver yielding a basis: simplex, or IPX/HiPO followed by
crossover (`solverYieldsBasis`). If the chosen solver can't yield one (e.g.
IPM with `run_crossover = off`), `solveLp` returns an error rather than
silently substituting another solver. With `require_basis = false`, PDLP may
also be chosen (by strategies 1 and 2; strategy 3 uses the same heuristic for
both).

A chosen `kHipo` is never actually returned unless HiPO is available in this
build (`HighsExternalApi::isAvailable<HighsExtras::hipo>()` — compiled in and
its deps, AMD/BLAS/METIS/RCM, loaded); otherwise it's substituted for `kIpx`,
the other interior-point method with identical basis/crossover semantics.
Callers never need to check availability themselves.

### MIP solves

A MIP node's LP relaxation always needs a basis (warm-starting, branching,
pseudocosts), so `HighsLpRelaxation::run()` hardcodes `require_basis = true`
on its `lpsolver` before dispatching — every call, every node, every
(sub-)MIP — rather than relying on an inherited option value. `strategy` *is*
read fresh from `mipsolver.options_mip_->solver_select_strategy` each time,
so it does propagate via the plain `HighsOptions` copy `solveSubMip` makes
for each sub-MIP it spawns.

Because a node LP's `integrality_` is cleared before it's passed to the LP
solver (it has to be a pure LP to solve as one), `num_integer_col` would
otherwise always read 0 there. `HighsLpRelaxation::loadModel()` counts it
from the original model before clearing and stashes it via
`Highs::setSolverSelectNumIntegerCol`, which `selectSolver`'s
`num_integer_col` parameter uses in place of recomputing it from (the now
cleared) `lp.integrality_`.

Re-selection only happens when the relaxation's `lpsolver` has no valid
warm-start basis (the first LP of a (sub-)MIP, or after a basis-invalidating
event) — most nodes just reuse simplex directly without calling
`selectSolver` again. Different sub-MIPs (RINS/RENS/diving spawn their own
`HighsMipSolver`, each with its own relaxation) routinely select differently
from the parent tree and from each other, since they're usually much smaller
after fixing. `lpsolver` is a plain default-constructed `Highs` — only
`output_flag`, `random_seed` and the feasibility tolerances are copied from
the MIP options — so `run_crossover` stays at the library default `"on"`
there regardless of what the user set at the top level: the
`require_basis`-with-crossover-off error case can't arise inside a MIP solve.

## API

```cpp
#include "util/HighsSolverSelect.h"

HighsLpFeatures f = computeLpFeatures(lp);         // one pass, does not modify lp
auto row  = highsLpFeatureVector(f);               // vector<pair<string,double>> — a CSV row
auto head = highsLpFeatureNames();                 // vector<string> — the CSV header

HighsInt n_int = countIntegerCols(lp.integrality_);  // kInteger/kSemiContinuous/kSemiInteger

HighsSolverSelect s = selectSolverByFeatures(lp);  // strategy 3 -> kDualSimplex | kIpx | kHipo
HighsSolverSelect b = selectSolverByFeatures(lp, /*require_basis=*/true);
HighsSolverSelect m = selectSolverByModel(f, SolverSelectStrategy::kFull,
                                          /*require_basis=*/true);
bool ok = solverYieldsBasis(b, /*run_crossover=*/true);

// full dispatch: strategy is SolverSelectStrategy, not a raw int
HighsSolverSelect e = selectSolver(lp, SolverSelectStrategy::kHeuristic,
                                   /*require_basis=*/false);
```

`lp` is expected to be the **presolved** LP. Helpers are safe on empty /
degenerate input (features fall back to 0).

### Tuning — `HighsLpFeatureParams`

| Field | Default | Meaning |
|-------|---------|---------|
| `almost_identical_tol` | `1e-9` | relative tolerance for "almost identical" coefficients |
| `dense_row_factor` | `10.0` | a row is "dense" past this multiple of the mean row length … |
| `dense_row_min_count` | `10` | … and at least this many nonzeros |

## Features

`m` = #rows, `n` = #cols, `nnz` = #matrix nonzeros. A bound is *finite* when
strictly inside `(-kHighsInf, kHighsInf)`. Ratios are guarded (0 if the
denominator is 0). Magnitude minima are taken over strictly nonzero values.

The **Rank** column is the predictive-power order (1 = strongest).

### Size

| Member | Value |
|--------|-------|
| `num_row`, `num_col`, `num_nz` | `m`, `n`, `nnz` |
| `num_integer_col` | integer / semi-continuous / semi-integer columns (0 for a pure LP) |

### How constrained?

| Member | Rank | Definition |
|--------|:----:|------------|
| `relative_num_equalities` | 1 | equality rows (`l^r_i = u^r_i`, finite) / `m` |
| `relative_num_cols_without_upper` | 12 | columns with `u_j = +inf` / `n` |
| `relative_num_inequalities` | | one-sided rows / `m` |
| `relative_num_ranged_rows` | | rows with both bounds finite and `l^r_i < u^r_i` / `m` |
| `relative_num_free_rows` | | rows with no finite bound / `m` |
| `relative_num_cols_without_lower` | | columns with `l_j = -inf` / `n` |
| `relative_num_free_cols` | | free columns / `n` |
| `relative_num_boxed_cols` | | columns with both bounds finite, `l_j < u_j` / `n` |
| `relative_num_singly_bounded_cols` | | columns with exactly one finite bound / `n` |
| `relative_num_fixed_cols` | | columns with `l_j = u_j` / `n` |

The free / boxed / singly-bounded / fixed ratios partition the columns and sum
to 1.

### How degenerate?

| Member | Rank | Definition |
|--------|:----:|------------|
| `relative_num_almost_identical_nonzeros` | 2 | fraction of `\|a_ij\|` lying in a cluster of ≥ 2 values that agree to `almost_identical_tol` |
| `relative_num_identical_obj_values` | 10 | fraction of `c_j` exactly equal to another `c_j` (zeros included) |
| `relative_num_identical_rhs_values` | 11 | same, over one rhs per row (finite `u^r_i`, else finite `l^r_i`); free rows excluded |

### How stable?

`max/min` magnitude ratios, with the raw extrema exposed for deriving other
proxies (e.g. `log10`).

| Member | Rank | Definition |
|--------|:----:|------------|
| `relative_max_matrix_coefficient` | 3 | `max\|a_ij\| / min\|a_ij\|` |
| `relative_max_obj_coefficient` | 4 | `max\|c_j\| / min\|c_j\|` |
| `relative_max_rhs_value` | 6 | `max\|b\| / min\|b\|` over finite row bounds (both bounds of a ranged row count) |
| `max_abs_*`, `min_abs_*` | | numerator / denominator of the three ratios above (`matrix_coefficient`, `obj_coefficient`, `rhs_value`) |

### How dense?

| Member | Rank | Definition |
|--------|:----:|------------|
| `matrix_nonzero_density` | 5 | `nnz / (m·n)` |
| `avg_nonzeros_per_col` | 7 | `nnz / n` |
| `relative_num_dense_rows` | 8 | rows longer than `max(dense_row_factor·(nnz/m), dense_row_min_count)` / `m` |
| `avg_nonzeros_per_row` | 9 | `nnz / m` |
| `max_nonzeros_in_row`, `max_nonzeros_in_col` | | longest row / column |

## Strategy 3 — hand-written heuristic

The provisional, hand-written fallback (`selectSolverByFeatures` /
`selectSolverRequiringBasis`/`selectSolverNotRequiringBasis`), predating the
fitted models above and still what strategy 3 evaluates:

- empty input → `kDualSimplex`
- small **and** sparse → `kDualSimplex`
- large / dense, very large and well-conditioned → `kHipo`
- large / dense, otherwise → `kIpx`

where *large* = `nnz > 5·10⁵`, *very large* = `nnz > 5·10⁶`,
*dense* = `density > 10⁻²` or `avg_nonzeros_per_row > 25` or
`relative_num_dense_rows > 10⁻²`, *ill-conditioned* =
`relative_max_matrix_coefficient > 10⁸` or `relative_max_rhs_value > 10⁸`.
Thresholds are provisional.

## Example — dump a feature matrix

```cpp
#include <fstream>
#include "util/HighsSolverSelect.h"

std::ofstream csv("features.csv");
csv << "instance";
for (const std::string& name : highsLpFeatureNames()) csv << ',' << name;
csv << '\n';

for (const auto& [name, presolved_lp] : corpus) {
  csv << name;
  for (const auto& [_, value] : highsLpFeatureVector(computeLpFeatures(presolved_lp)))
    csv << ',' << value;
  csv << '\n';
}
```

## Batch tool

A CLI that walks a directory of saved (presolved) model files and writes one
CSV row per file — the loop above, wrapped up — lives in the sibling `tests`
repo at `select/extract_features`, next to where the PCA on this feature set
happens.
