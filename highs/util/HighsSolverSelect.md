# Solver selection by LP features

Pick an LP solver (dual simplex, IPX or HiPO) from cheap structural features of
the presolved model.

## Workflow

The fitting tools are in the sibling `tests` repo, under `select/`.

1. Presolve: `presolve_set.sh` writes the presolved model of each instance.
2. Time: `run_solvers.sh <version> LP|MIP` runs the three solvers on the
   presolved LP sets and on the LP relaxations of presolved MIPs.
3. Extract: `extract_features/` applies LP presolve to each presolved model and
   writes `computeLpFeatures` to a CSV, so that fitting sees what `solveLp` sees.
4. Collect: `collect_results.py --version <version>` writes `select.csv`,
   checking from the logs that the requested solver ran and how crossover ended.
   Each solve is flagged by whether it is solved and whether it yields a basis.
5. Fit: `visualize.py --cpp-out .../HighsSolverSelectModels.h` fits gradient
   boosted trees: for each solver, 100 small trees that predict how much slower
   it is than the fastest. The solver with the smallest prediction is selected.
   The settings are fixed rather than tuned, and scores come from repeated cross
   validation.

## Strategies: `solver_select_strategy`

Applies when `solver = choose`, or in the MIP solver when
`mip_lp_solver = choose`. The range is `[0, 3]` and the default is 0.

| Value | `SolverSelectStrategy` | Behaviour |
|:-----:|------------------------|-----------|
| 0 | `kSimplexOnly` | Dual simplex, with no selection. |
| 1 | `kColsRowsNnz` | Fitted model on `num_col`, `num_row` and `num_nz`. |
| 2 | `kRankedFeatures` | Fitted model on the 12 ranked features below. |
| 3 | `kHeuristic` | The provisional heuristic below. |

Strategies 1 and 2 send LPs smaller than the smallest training instance to dual
simplex. Selected HiPO becomes IPX when HiPO is not available, in `solveLp`.

`solver_select_require_basis` (default `false`) makes `solveLp` end with a basis.
If IPM is selected and `run_crossover` is off, it returns an error rather than
substituting another solver.

IPM does not always yield a basis, even with `run_crossover = on`. A solve with
`Crossover not run` in its log, `Optimal` status and no simplex iterations has no
basis, so it does not count when fitting. With simplex iterations, IPM did not
finish and simplex found the basis. Selecting that solver costs the recorded time,
so it counts. On the 1.15.1 data every solved instance has a basis, so one model
per strategy is enough.

If a selected IPM fails, `solveLp` runs simplex. It also does so if a basis is
required and IPM claims optimality, infeasibility or unboundedness without one.

HiPO stops with `Memory limit reached` if its predicted memory exceeds what is
available, including any container limit. If HiPO was selected, `solveLp` then runs
IPX, which needs far less memory. An explicit `solver = hipo` still ends with the
memory error. This check applies to HiPO with any strategy, unlike selection.

If the incumbent has a basis, dual simplex is used without evaluating a
strategy. `Highs::run` has skipped presolve by then, on the assumption that the
basis will be used.

### MIP solves

`HighsLpRelaxation::run` uses selection only when `lpsolver` has no basis, and
always requires a basis. This covers the root LP of each (sub-)MIP and the
LP-based heuristics that start without a basis, so a MIP can select many times.
Selection happens in `solveLp`, so the features are those of the LP after any
presolve by `lpsolver`. If the selected solver fails, simplex is tried.

Selection is skipped when the MIP solver has set primal simplex as a last
resort. Iteration limits in the MIP solver are multiples of the root LP's count,
so an IPM solve counts its IPM and crossover iterations. An LP with a simplex
iteration limit, such as in a heuristic, limits selected IPM to 200 iterations.
Crossover has no iteration limit.

With `log_dev_level` at least 1, each selection in the main MIP logs its IPM and
simplex iteration counts.

## Strategy 3: hand-written heuristic

- empty input: `kDualSimplex`
- small and sparse: `kDualSimplex`
- large or dense, very large and well-conditioned: `kHipo`
- large or dense, otherwise: `kIpx`

*large* is `nnz > 5e5`, *very large* is `nnz > 5e6`. *dense* is `density > 1e-2`,
`avg_nonzeros_per_row > 25` or `relative_num_dense_rows > 1e-2`.
*ill-conditioned* is `relative_max_matrix_coefficient > 1e8` or
`relative_max_rhs_value > 1e8`. The thresholds are provisional.

## Features

`m` = #rows, `n` = #cols, `nnz` = #matrix nonzeros. A bound is *finite* when
strictly inside `(-kHighsInf, kHighsInf)`. Ratios are guarded (0 if the
denominator is 0). Magnitude minima are taken over strictly nonzero values.

The **Rank** column is the predictive-power order (1 = strongest). Strategy 2
uses the ranked features. The rest are extracted for analysis.

### Size

| Member | Value |
|--------|-------|
| `num_row`, `num_col`, `num_nz` | `m`, `n`, `nnz` |
| `num_integer_col` | integer / semi-continuous / semi-integer columns. Extracted only: `select/test_integer_count.py` checks that it adds little |

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
| `relative_num_almost_identical_nonzeros` | 2 | fraction of `\|a_ij\|` lying in a cluster of ≥ 2 values that agree to `1e-9 * max(1, |a_ij|)` |
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
| `relative_num_dense_rows` | 8 | rows longer than `max(10·(nnz/m), 10)` / `m` |
| `avg_nonzeros_per_row` | 9 | `nnz / m` |
| `max_nonzeros_in_row`, `max_nonzeros_in_col` | | longest row / column |
