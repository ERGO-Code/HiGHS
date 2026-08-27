# Legacy HiPDLP Characterization

## Scope and provenance

This record characterizes legacy HiPDLP before the rewrite changes its
algorithm, transformations, or scaling. The frozen Milestone 0 artifacts are
under:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/
  m0-legacy-serial-e29a6c0-20260825
```

They were produced from commit
`e29a6c04bc28a0b5ebcfa7b42b53f221d4f3d73d` with the serial CPU build,
presolve disabled, and `kkt_tolerance = 1e-4`. Exact commands, options,
checksums, stdout, stderr, and reference-solver results remain in that
immutable directory.

The executable characterization tests are in `check/TestPdlpHi.cpp`. They use
the independent verifier in
`highs/pdlp/hipdlp/verification/kkt_verifier.{hpp,cc}` and separately call
`getLpKktFailures` as a cross-check. The verifier does not call `PDLPSolver`,
legacy convergence checks, or transformed/scaled model data.

## Independent acceptance gates

The verifier consumes the original `HighsLp` and final unscaled
`HighsSolution`. It independently computes row activities, reduced costs,
bound and dual-sign violations, row-activity and stationarity residuals,
primal and dual objectives, objective gap, complementarity, and worst
offending original row/column indices. It validates sparse storage before any
matrix product and accepts only valid HiGHS infinity bounds; all other numeric
inputs and all computed metrics must be finite. A solution is accepted only
when:

```text
relative primal residual <= 1e-4
relative dual residual   <= 1e-4
relative objective gap   <= 1e-4
```

The primal bound residual is normalized by one plus the infinity norm of
original-LP bounds active within `sqrt(tolerance)`. An infeasible value outside
that neighborhood also includes its violated bound in the local measure, so a
large inactive bound cannot hide an unrelated violation. Row-activity
residuals use the same active-bound norm.

The dual cost measure is the infinity norm of original objective coefficients
whose independently recomputed reduced costs are within `sqrt(tolerance)` of
zero. A larger column sign violation also includes its own objective
coefficient in the local measure; stationarity residuals use the near-zero
cost norm. This excludes large inactive costs while preserving model-relative
scaling. Dual-objective bound contributions are selected from the
sense-adjusted sign of each independently recomputed reduced cost or row dual,
not from the reported primal value. Objective gap uses
`abs(primal - dual) / (1 + abs(primal) + abs(dual))`.
Complementarity is retained as a raw diagnostic with its own worst indices;
the acceptance decision follows the three primary gates above. For a
primal-feasible, dual-feasible LP solution, the objective-gap gate also
captures aggregate complementarity error.

## Frozen traces

| Model and policy | Legacy trace | Public result | Independent conclusion |
| --- | --- | --- | --- |
| `afiro`, default option translated to PID, presolve off | Iteration 0: primal infeasibility `4.40e+01`, dual infeasibility `9.23e-01`; iteration 280: primal objective `-464.759739`, dual objective `-464.782381`, gap `2.26e-02` | `Optimal`, objective `-464.75973945` | Passes the strict `1e-4` gates |
| `scrs8`, default option translated to PID, presolve off | Iteration 865,080: primal objective `1.20e+06`, dual objective `-1.16e+11`; iteration 9,516,960: gap `2.42e+10`; iteration 10,381,800 at 60 seconds | `Time limit reached` | Divergent; not acceptable |
| `scrs8`, fixed step, presolve off | Iteration 21,040: primal objective `907.714122`, dual objective `907.834353`, solver gap `1.20e-01` | `Optimal`, objective `907.71412156` | Reject: relative primal infeasibility is `1.76e-4`; trusted simplex objective is `904.29695380` |

The serial-CPU focused test bounds the default `scrs8` replay at 40,000
iterations and asserts `HighsStatus::kWarning`,
`HighsModelStatus::kIterationLimit`, and independent rejection. This preserves
the defect without adding a 60-second unit test. Its final primal and row-dual
vectors must exactly match a second run configured explicitly with PID, and
both the independent verifier and `getLpKktFailures` reject the result. The
fixed replay remains fast and deterministically reaches the frozen false
optimum at iteration 21,040. Exact iteration characterizations are excluded
from CUDA builds because these artifacts and traces were frozen on the serial
backend.

## Option inventory

This table records legacy HiPDLP behavior in `PDLPSolver::setup`. It does not
describe in-tree cuPDLP behavior.

| Option | Legacy HiPDLP behavior |
| --- | --- |
| `pdlp_iteration_limit` | Honored as the solve-loop iteration limit. |
| `time_limit` | Honored using the shared HiGHS timer. |
| `pdlp_optimality_tolerance` | Honored unless a nondefault `kkt_tolerance` overrides it. |
| `kkt_tolerance` | When nondefault, overrides the internal HiPDLP optimality tolerance. |
| `pdlp_scaling_mode` | Ruiz, Pock-Chambolle, and L2 bits are honored by the existing legacy scaling sequence. |
| `pdlp_ruiz_iterations` | Honored by the legacy Ruiz stage. |
| `pdlp_features_off` scaling bit | Honored and disables all legacy scaling stages. |
| `pdlp_restart_strategy` | Ignored. The solve loop always uses adaptive Halpern restart and logs that the option was ignored. |
| `pdlp_features_off` restart bit | Ignored. Halpern restart remains enabled and a message is logged. |
| `pdlp_cupdlpc_restart_method` | Not read by HiPDLP; it is specific to cuPDLP-C. |
| `pdlp_step_size_strategy = 0` | Honored as fixed step size. |
| `pdlp_step_size_strategy = 3` | Honored as PID. |
| `pdlp_step_size_strategy = 1` (global default adaptive) | Unsupported and silently forced to PID after an informational log message. |
| `pdlp_step_size_strategy = 2` (Malitsky-Pock) | Unsupported and forced to PID after an informational log message. |
| `presolve` | Applied by the surrounding HiGHS solve path, not by HiPDLP option translation. |

Legacy HiPDLP has no Curtis-Reid scaling option or stage. Curtis-Reid belongs
to the later transformation/scaling milestone and is intentionally absent from
this characterization work.

## Known defects retained for rewrite arbitration

- Default `scrs8` adaptive-to-PID translation diverges and reaches the time
  limit with a relative gap above `1e10`.
- Fixed-step `scrs8` reports `Optimal` even though independent and core KKT
  checks reject it at the strict `1e-4` gate.
- In the 13-instance frozen batch, `25fv47`, `adlittle`, `avgas`, `standata`,
  and `standgub` also report `Optimal` while failing at least one strict final
  KKT metric.
- Legacy restart and unsupported step-size options are reinterpreted instead
  of being rejected or reported as an effective policy through the public
  result.

These defects are characterization expectations, not behavior to preserve in
the rewritten convergence policy.