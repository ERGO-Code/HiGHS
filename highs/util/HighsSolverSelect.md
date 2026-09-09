# HiGHS solver selection by structural features

Facility for choosing an LP algorithm — dual simplex, IPX or HiPO — from cheap
structural features of the (presolved) model, and for extracting those features
in bulk so that a selection model can be fitted offline.

Source files: [`HighsSolverSelect.h`](HighsSolverSelect.h),
[`HighsSolverSelect.cpp`](HighsSolverSelect.cpp).

## Background

The feature set mirrors the one described by Zonghao Gu in the Gurobi 11.0
_"New Performance Techniques"_ talk (16 October 2023), slide 6 ("SolutionTarget —
Yes we can!"). Gurobi collected ~69 features of the presolved model over a set
of ~2200 instances that do not solve in presolve, added ~20 derived
("relative") features, and cross-fit a decision tree to decide between a vertex
solution (simplex-style) and an interior solution (barrier), weighting
instances by their performance ratio and iteratively dropping unimportant
features and replacing absolute features with matching relative ones.

The features that survived, sorted by predictive power:

| #  | Slide-6 feature                                | Intuition group |
|----|------------------------------------------------|-----------------|
| 1  | relative number of equalities                  | how constrained |
| 2  | relative number of almost identical nonzeros   | how degenerate  |
| 3  | relative maximum nonzero coefficient           | how stable      |
| 4  | relative maximum objective coefficient         | how stable      |
| 5  | nonzero density of constraint matrix           | how dense       |
| 6  | relative maximum of rhs values                 | how stable      |
| 7  | number of nonzeros per column                  | how dense       |
| 8  | relative number of "dense" rows                | how dense       |
| 9  | number of nonzeros per row                     | how dense       |
| 10 | relative number of identical obj values        | how degenerate  |
| 11 | relative number of identical rhs values        | how degenerate  |
| 12 | relative number of variables without upper bound | how constrained |

The four intuition groups from the slide:

- **How constrained is the problem?** Barrier likes having an interior.
- **How degenerate is the problem?** Simplex could get stuck.
- **How stable is the problem?** Simplex and barrier suffer in different ways.
- **How dense is the problem?** Simplex thrives with sparsity.

Every "relative" feature is **scale free**: a count divided by the number of
rows, columns or nonzeros, or a max-over-min magnitude ratio. This is what
makes the features comparable across instances of very different size, and
suitable for PCA / a decision tree.

## Intended workflow

1. **Offline, over a corpus of presolved instances.** For each instance call
   `computeLpFeatures` and append `highsLpFeatureVector` as a row of a CSV whose
   header is `highsLpFeatureNames`. Join with a label per instance (which of
   simplex / IPX / HiPO was fastest, and by how much).
2. **Offline, model fitting.** Run PCA / fit a decision tree or small
   classifier on that feature matrix, weighting instances by performance ratio
   (large win/loss ⇒ more weight), as in the Gurobi procedure.
3. **Online, per new instance.** Replace the body of
   `selectSolverByFeatures(const HighsLpFeatures&)` with the fitted classifier.
   The feature definitions in this file are the contract between the offline
   study and the online selector — keep them in lock-step.

Until step 3 is done, `selectSolverByFeatures` applies a transparent
placeholder heuristic (see [below](#placeholder-heuristic)).

## API

```cpp
#include "util/HighsSolverSelect.h"
```

### `HighsLpFeatures computeLpFeatures(const HighsLp& lp, const HighsLpFeatureParams& params = {})`

Computes the feature record for `lp`. `lp` is not modified; the constraint
matrix is read in whatever orientation (row-wise or column-wise) it is already
stored in, in a single pass. Safe on degenerate inputs (no rows, no columns, no
nonzeros, an unsized matrix): the corresponding features are left at 0.

`lp` is expected to be the **presolved** LP — that is what the Gurobi study used
and what makes the features discriminating.

### `std::vector<std::pair<std::string, double>> highsLpFeatureVector(const HighsLpFeatures&)`

Flattens the record into `(name, value)` pairs in a fixed order — one CSV row.

### `std::vector<std::string> highsLpFeatureNames()`

The names in the same order — the CSV header. Derived from
`highsLpFeatureVector` so the two cannot drift apart.

### `HighsSolverSelect selectSolverByFeatures(const HighsLp& lp, const HighsLpFeatureParams& params = {})`
### `HighsSolverSelect selectSolverByFeatures(const HighsLpFeatures& features)`

Returns one of `HighsSolverSelect::kDualSimplex`, `::kIpx`, `::kHipo`. The first
overload is `computeLpFeatures` followed by the second.

### `HighsSolverSelect selectSolver(const HighsLp& lp)`

The original entry point, called from `HighsSolve.cpp` when
`options.solver == "choose"`. Currently an unchanged stub returning
`kDualSimplex`; wire it to `selectSolverByFeatures` once the classifier is in.

## `HighsLpFeatureParams`

Tuning constants for the feature computation, collected in one struct so the
offline study and the online selector use identical definitions.

| Field                  | Default | Meaning |
|------------------------|---------|---------|
| `almost_identical_tol` | `1e-9`  | Two coefficients are "almost identical" when their absolute values agree to this relative tolerance. Feeds feature 2. |
| `dense_row_factor`     | `10.0`  | A row is "dense" when its nonzero count exceeds this multiple of the mean row length … |
| `dense_row_min_count`  | `10`    | … and is also at least this many nonzeros. Feeds feature 8. |

## Feature reference

Notation: `m` = number of rows, `n` = number of columns, `nnz` = number of
matrix nonzeros. `l^c_j, u^c_j` are column bounds; `l^r_i, u^r_i` are row
(activity) bounds. A bound is "finite" when strictly inside
`(-kHighsInf, kHighsInf)`. All ratios are guarded: a zero denominator yields 0.

### Raw size

| Member            | Definition |
|-------------------|------------|
| `num_row`         | `m` |
| `num_col`         | `n` |
| `num_nz`          | `nnz` (`lp.a_matrix_.numNz()`) |
| `num_integer_col` | columns whose `integrality_` is `kInteger`, `kSemiContinuous` or `kSemiInteger` (0 for a pure LP) |

### How constrained is the problem?

Barrier methods want a non-empty interior; lots of equalities or fixed
variables squeeze it.

| Member | Slide # | Definition |
|--------|:------:|------------|
| `relative_num_equalities`          | 1  | `#{i : l^r_i = u^r_i, both finite} / m` |
| `relative_num_inequalities`        | –  | `#{i : exactly one of l^r_i, u^r_i finite} / m` |
| `relative_num_ranged_rows`         | –  | `#{i : both finite, l^r_i < u^r_i} / m` |
| `relative_num_free_rows`           | –  | `#{i : neither bound finite} / m` |
| `relative_num_cols_without_upper`  | 12 | `#{j : u^c_j = +inf} / n` |
| `relative_num_cols_without_lower`  | –  | `#{j : l^c_j = -inf} / n` |
| `relative_num_free_cols`           | –  | `#{j : l^c_j = -inf and u^c_j = +inf} / n` |
| `relative_num_boxed_cols`          | –  | `#{j : both bounds finite, l^c_j < u^c_j} / n` |
| `relative_num_singly_bounded_cols` | –  | `#{j : exactly one bound finite} / n` |
| `relative_num_fixed_cols`          | –  | `#{j : l^c_j = u^c_j, both finite} / n` |

The four column classes (free / boxed / singly-bounded / fixed) partition the
columns, so those four ratios sum to 1.

### How degenerate is the problem?

Repeated coefficients / right-hand sides / objective values create ties and
alternative optima where simplex can stall.

| Member | Slide # | Definition |
|--------|:------:|------------|
| `relative_num_almost_identical_nonzeros` | 2  | Sort the `\|a_ij\|`. Walk them, joining consecutive values into a cluster while `next - cur <= almost_identical_tol * max(1, \|next\|)`. Result = (number of nonzeros in clusters of size ≥ 2) / `nnz`. |
| `relative_num_identical_obj_values`      | 10 | Fraction of the `n` objective coefficients `c_j` that are **exactly** equal to at least one other `c_j` (sort, count runs of length ≥ 2). Zeros are included. |
| `relative_num_identical_rhs_values`      | 11 | Same duplicate-fraction measure over one representative right-hand side per row: the finite `u^r_i` if present, else the finite `l^r_i`. Free rows contribute nothing, so the denominator is the number of non-free rows. |

### How stable is the problem?

Wide coefficient / bound ranges mean poor conditioning; simplex and barrier
degrade differently. Each "relative maximum" is a max-over-min magnitude ratio;
minima are taken over strictly nonzero magnitudes. The raw extrema are exposed
too so other conditioning proxies (e.g. `log10` of the ratio) can be derived
offline.

| Member | Slide # | Definition |
|--------|:------:|------------|
| `relative_max_matrix_coefficient` | 3 | `max\|a_ij\| / min\|a_ij\|` over nonzeros (0 if `nnz = 0`) |
| `relative_max_obj_coefficient`    | 4 | `max\|c_j\| / min\|c_j\|` over `c_j != 0` (0 if objective is all zero) |
| `relative_max_rhs_value`          | 6 | `max\|b\| / min\|b\|` over every finite row bound `b` with `\|b\| > 0` (both bounds of a ranged row count) |
| `max_abs_matrix_coefficient`, `min_abs_matrix_coefficient` | – | numerator / denominator of feature 3 |
| `max_abs_obj_coefficient`, `min_abs_obj_coefficient`       | – | numerator / denominator of feature 4 |
| `max_abs_rhs_value`, `min_abs_rhs_value`                   | – | numerator / denominator of feature 6 |

### How dense is the problem?

Simplex thrives on sparsity; dense rows in particular are expensive.

| Member | Slide # | Definition |
|--------|:------:|------------|
| `matrix_nonzero_density` | 5 | `nnz / (m * n)` |
| `avg_nonzeros_per_col`   | 7 | `nnz / n` |
| `avg_nonzeros_per_row`   | 9 | `nnz / m` |
| `relative_num_dense_rows`| 8 | `#{i : rowlen_i > max(dense_row_factor * (nnz / m), dense_row_min_count)} / m` |
| `max_nonzeros_in_row`    | – | longest row |
| `max_nonzeros_in_col`    | – | longest column |

## Placeholder heuristic

`selectSolverByFeatures(const HighsLpFeatures&)` currently implements, as a
stand-in for the fitted model:

- degenerate input (`nnz`, `m` or `n` zero) ⇒ `kDualSimplex`;
- otherwise classify size (`nnz > 5e5` "large", `> 5e6` "very large"),
  density (`matrix_nonzero_density > 1e-2` or `avg_nonzeros_per_row > 25` or
  `relative_num_dense_rows > 1e-2`) and conditioning
  (`relative_max_matrix_coefficient > 1e8` or `relative_max_rhs_value > 1e8`
  ⇒ "ill-conditioned");
- small **and** sparse ⇒ `kDualSimplex` (thrives on sparsity, cheap warm start);
- otherwise an interior-point method: `kHipo` for very large and
  well-conditioned models, `kIpx` otherwise (its linear algebra / crossover are
  the more battle-tested on hard cases).

The thresholds are provisional and exist only so the entry point returns
something sensible before the study is done. Replace the whole body with the
classifier from step 2 of the workflow.

## Example: dumping a feature matrix

```cpp
#include <fstream>
#include "util/HighsSolverSelect.h"

std::ofstream csv("features.csv");

// header
csv << "instance";
for (const std::string& name : highsLpFeatureNames()) csv << ',' << name;
csv << '\n';

// one row per presolved instance
for (const auto& [instance_name, presolved_lp] : corpus) {
  csv << instance_name;
  for (const auto& [name, value] : highsLpFeatureVector(computeLpFeatures(presolved_lp)))
    csv << ',' << value;
  csv << '\n';
}
```
