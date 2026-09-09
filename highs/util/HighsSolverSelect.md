# Solver selection by structural features

Pick an LP algorithm — dual simplex, IPX or HiPO — from cheap structural
features of the presolved model, and extract those features in bulk so a
selection model can be fitted offline.

Files: [`HighsSolverSelect.h`](HighsSolverSelect.h) ·
[`HighsSolverSelect.cpp`](HighsSolverSelect.cpp)

## Idea

From Zonghao Gu's Gurobi 11.0 talk *"New Performance Techniques"* (2023),
slide 6. Gurobi computed ~90 features of the presolved model over ~2200 hard
instances and fitted a decision tree to choose between a simplex-style vertex
solution and a barrier-style interior solution.

Features are **scale free** — a count over `#rows`/`#cols`/`#nonzeros`, or a
`max/min` magnitude ratio — so they compare across instances of any size. They
answer four questions:

| Question | Why it matters |
|----------|----------------|
| How **constrained**? | Barrier wants a non-empty interior |
| How **degenerate**?  | Simplex can stall on ties / alternative optima |
| How **stable**?      | Wide value ranges hurt simplex and barrier differently |
| How **dense**?       | Simplex thrives on sparsity |

## Workflow

1. **Extract** — for each presolved instance, write `highsLpFeatureVector` as a
   CSV row (header from `highsLpFeatureNames`); add a label for which solver was
   fastest.
2. **Fit** — run PCA / train a classifier on that matrix offline.
3. **Deploy** — replace the body of `selectSolverByFeatures(const
   HighsLpFeatures&)` with the fitted classifier. Until then it uses a
   [placeholder heuristic](#placeholder-heuristic).

The feature definitions here are the contract between steps 1 and 3 — keep them
in sync.

## API

```cpp
#include "util/HighsSolverSelect.h"

HighsLpFeatures f = computeLpFeatures(lp);         // one pass, does not modify lp
auto row  = highsLpFeatureVector(f);               // vector<pair<string,double>> — a CSV row
auto head = highsLpFeatureNames();                 // vector<string> — the CSV header

HighsSolverSelect s = selectSolverByFeatures(lp);  // -> kDualSimplex | kIpx | kHipo
```

`lp` is expected to be the **presolved** LP. All helpers are safe on empty /
degenerate input (features fall back to 0). `selectSolver(lp)` is the original
entry point (still a stub returning `kDualSimplex`); point it at
`selectSolverByFeatures` once the classifier lands.

### Tuning — `HighsLpFeatureParams`

| Field | Default | Meaning |
|-------|---------|---------|
| `almost_identical_tol` | `1e-9` | relative tolerance for "almost identical" coefficients |
| `dense_row_factor` | `10.0` | a row is "dense" past this multiple of the mean row length … |
| `dense_row_min_count` | `10` | … and at least this many nonzeros |

## Features

`m` = #rows, `n` = #cols, `nnz` = #matrix nonzeros. A bound is *finite* when
strictly inside `(-kHighsInf, kHighsInf)`. Every ratio is guarded (0 if the
denominator is 0). Magnitude minima are taken over strictly nonzero values.

The **Rank** column is the slide-6 predictive-power order (1 = strongest).

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

## Placeholder heuristic

Stand-in until the fitted model is wired in:

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
