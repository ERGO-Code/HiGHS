# HiPDLP Phase 7 Cross-Solver Validation Report

Date: 2026-09-03

Status: **INCOMPLETE. Do not proceed to Phase 8.**

## Findings

1. Phase 7 is not complete. The full Netlib matrix is missing for rewrite CUDA
   and all comparison solvers, and no energy regression or full-energy run
   exists.
2. The result root contained 26 physical directories when audited, not 10.
   Four cannot support a completion claim: one empty cuPDLPx attempt, the old
   rewrite CUDA regression that crossed a source-change boundary, and the old
   smoke and reference-triage directories that have no run-level provenance.
3. Fresh, provenance-stable smoke runs pass rewrite CPU and CUDA 13/13. Their
   statuses and iteration counts agree on all 13 instances; the largest paired
   objective relative difference is `1.97e-12` on `scrs8`.
4. Fresh Netlib regression runs accept rewrite CPU on 97/103 and CUDA on
   96/103. Every rewrite `OPTIMAL` result passes the original-LP `1e-4` primal,
   dual, and gap checks. All rejected rewrite rows are time limits.
5. `pilot4` is the only fresh CPU/CUDA terminal disagreement. Equal iteration
   budgets produce KKT-passing `ITERATION_LIMIT` results on both backends. The
   first restart split is a strict comparison at the `1e-13` scale, and both
   backends take the same controlled PID update one 40-iteration cadence apart.
   This is `TIME_OR_ITERATION_LIMIT`, not a CUDA backend defect.
6. The earlier false infeasible claims on `fffff800`, `gfrd-pnc`, `greenbea`,
   and `vtp.base` were real rewrite defects and are resolved by the certificate
   fix. The old pre-fix directories remain characterization evidence only.
7. cuOpt 26.08 reports finite but incorrect all-zero reduced costs on AFIRO.
   The raw values are retained and marked inconsistent; KKT uses reduced costs
   independently derived in C++ from the original LP and cuOpt row duals.
8. cuPDLPx output naming truncates dotted instance names at the first dot. The
   adapter previously misclassified `vtp.base`, `pilot.ja`, and `pilot.we` as
   `SOLVE_ERROR`. A focused test and one-line adapter fix preserve the real
   outcomes: two accepted `OPTIMAL` results and one `TIME_LIMIT`.

## Repository Identity

The audited HiPDLP checkout is branch `feature/hipdlp-rewrite` at
`db17f8cfbea6e7ce79b0c76239570c3e32cca635`. The working tree is dirty and
contains the rewrite solver, certificate verifier, adapters, tests, and plans.
Results must therefore be identified by the recorded executable, manifest, and
configuration hashes plus the before/after provenance snapshots, not by `HEAD`
alone.

The public default selector in this uncommitted working tree already points to
the rewrite. That pre-existing overlap with Phase 8 does not authorize legacy
deletion or any additional cutover work.

Machine-readable evidence is under
`/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase7`:

- `evidence-audit.csv`: one row for each of the original 26 directories;
- `summary.csv`: 977 canonical solver/instance rows for executed scopes;
- `disagreements.csv`: 270 classified regression or historical disagreements;
- `generate_canonical_summary.py`: aggregation logic with row-count,
  uniqueness, and zero-`UNRESOLVED` assertions.

## Original Directory Audit

| Directory | Rows / accepted | Decision | Reason |
| --- | ---: | --- | --- |
| `certificate-fix-cpu-20260902-01` | 4 / 3 | defect evidence | Stable four-case CPU rerun; `greenbea` remains a time limit. |
| `certificate-fix-cuda-20260902-01` | 4 / 3 | defect evidence | Stable four-case CUDA rerun; `greenbea` remains a time limit. |
| `netlib-full-hipdlp-rewrite-cpu-20260903-01` | 114 / 111 | partial Phase 7 evidence | Complete rewrite CPU full Netlib run only. |
| `netlib-full-hipdlp-rewrite-cpu-resume-b01-20260903` | 10 / 10 | duplicate supplement | Superseded by the complete 114-row run. |
| `netlib-full-hipdlp-rewrite-cpu-resume-b02-20260903` | 10 / 9 | duplicate supplement | Superseded by the complete 114-row run. |
| `netlib-regression-cuopt-26.08-20260903-01` | 103 / 45 | baseline evidence | Pinned adapter and cuOpt; reported reduced costs remain defective. |
| `netlib-regression-cupdlp-in-tree-20260903-01` | 103 / 40 | baseline evidence | Binary, manifest, and config identities are stable. |
| `netlib-regression-cupdlpx-20260903-01` | 0 / 0 | not evidence | Empty initial attempt with the wrong configuration. |
| `netlib-regression-cupdlpx-20260903-02` | 103 / 47 | partial baseline evidence | Three dotted names were misparsed by the adapter; corrected separately. |
| `netlib-regression-highs-ipm-20260903-01` | 103 / 103 | reference evidence | Stable high-accuracy IPM results. |
| `netlib-regression-highs-simplex-20260903-01` | 103 / 102 | reference evidence | `qap15` reaches the 60-second limit. |
| `netlib-regression-hipdlp-legacy-20260902-01` | 101 / 28 | merged legacy shard | Incomplete first pass; deduplicate with later shards. |
| `netlib-regression-hipdlp-legacy-20260902-02` | 30 / 11 | merged legacy shard | Overlapping retry shard. |
| `netlib-regression-hipdlp-legacy-resume-b01-20260902` | 10 / 4 | merged legacy shard | Stable continuation shard. |
| `netlib-regression-hipdlp-legacy-resume-b02-20260902` | 10 / 2 | merged legacy shard | Stable continuation shard. |
| `netlib-regression-hipdlp-legacy-resume-b03-20260903` | 10 / 2 | merged legacy shard | Stable continuation shard. |
| `netlib-regression-hipdlp-legacy-resume-b04-20260903` | 10 / 2 | merged legacy shard | Stable continuation shard. |
| `netlib-regression-hipdlp-legacy-resume-b05-20260903` | 10 / 1 | merged legacy shard | Stable continuation shard. |
| `netlib-regression-hipdlp-legacy-resume-b06-20260903` | 10 / 3 | merged legacy shard | Stable continuation shard. |
| `netlib-regression-hipdlp-legacy-resume-b07-20260903` | 10 / 0 | merged legacy shard | Stable continuation shard. |
| `netlib-regression-hipdlp-legacy-resume-b08-20260903` | 3 / 3 | merged legacy shard | Final continuation shard. |
| `netlib-regression-hipdlp-rewrite-cpu-20260902-01` | 103 / 93 | historical defect evidence | Pre-certificate-fix and pre-reference-correction. |
| `netlib-regression-hipdlp-rewrite-cpu-20260902-02` | 103 / 97 | superseded evidence | Stable old binary; replaced by the fresh provenance run. |
| `netlib-regression-hipdlp-rewrite-cuda-20260902-01` | 103 / 95 | not evidence | Run spanned a source/build edit at 12:30 local time. |
| `netlib-regression-reference-triage-20260902-01` | 20 / n/a | not evidence | No metadata, command identity, or hashes. |
| `smoke-20260902-01` | 39 / 31 | not evidence | Numerically useful but no run-level provenance. |

The canonical summary replaces invalid or superseded rewrite evidence with
`smoke-recheck-20260903-02` and
`netlib-regression-recheck-20260903-01`. It merges the legacy shards by latest
timestamp and replaces the three dotted cuPDLPx rows with
`cupdlpx-exit-fix-20260903-02`.

## Task Status

| Task | Status | Evidence |
| --- | --- | --- |
| 1. Freeze execution point | Complete, 3/3 | CPU/CUDA builds and focused routing checks pass; dirty tree recorded. |
| 2. HiGHS-family adapter | Partial, 5/6 | Fresh AFIRO simplex/IPM/cuPDLP results pass; historical missing-target RED was not retained. |
| 3. External verifier | Partial, 4/5 | Valid AFIRO passes and one-value perturbation fails; historical missing-target RED was not retained. |
| 4. cuOpt adapter | Complete, 5/5 | Pinned test passes; reduced-cost baseline defect is explicit. |
| 5. Harness configs | Complete, 6/6 | Six runner tests and seven dry runs pass; dry runs create no result directory. |
| 6. Thirteen-instance smoke | Complete, 6/6 | Rewrite CPU and CUDA both 13/13 with stable provenance. |
| 7. Cross-solver suites | Partial, 3/6 | Regression executed and classified; full matrix and both energy scopes missing. |
| 8. Triage disagreements | Partial, 4/5 | Regression has no `UNRESOLVED`; later scopes do not yet exist. |
| 9. Report | Complete, 5/5 | This report and the three machine-readable audit files agree. |

Overall checklist state after this report: 41 checked and 6 unchecked out of
47. An unchecked historical RED is not retroactively treated as proof.

## Result Counts

### Smoke

| Solver | Accepted | Statuses |
| --- | ---: | --- |
| rewrite CPU | 13/13 | 13 `OPTIMAL` |
| rewrite CUDA | 13/13 | 13 `OPTIMAL` |
| legacy | 5/13 | 10 `OPTIMAL`, 3 process timeouts; 5 optimal claims fail KKT |

### Netlib Regression

| Solver | Accepted | Statuses |
| --- | ---: | --- |
| rewrite CPU | 97/103 | 97 `OPTIMAL`, 6 `TIME_LIMIT` |
| rewrite CUDA | 96/103 | 96 `OPTIMAL`, 7 `TIME_LIMIT` |
| legacy | 28/103 | 81 `OPTIMAL`, 22 missing after process timeout |
| in-tree cuPDLP | 40/103 | 95 `OPTIMAL`, 4 `UNKNOWN`, 4 `TIME_LIMIT` |
| cuPDLPx | 49/103 | 99 `OPTIMAL`, 3 `TIME_LIMIT`, 1 parser `SOLVE_ERROR` |
| cuOpt 26.08 | 45/103 | 99 `OPTIMAL`, 3 `TIME_LIMIT`, 1 parser error |
| HiGHS simplex | 102/103 | 102 `OPTIMAL`, 1 `TIME_LIMIT` |
| HiGHS IPM | 103/103 | 103 `OPTIMAL` |

Low baseline acceptance counts reflect the strict independent `1e-4` KKT and
`1e-3` reference-objective gates. They are not substituted for solver-reported
status.

### Full Netlib

Only rewrite CPU has run: 111/114 accepted. `pilot`, `fit2d`, and `greenbea`
reach the 300-second limit. All 111 `OPTIMAL` rows pass the original-LP KKT
gate. The other seven solver/backend paths are not run.

### Energy

Energy regression: not run for any solver. Full energy: not run for any solver.
The resource probe found sufficient host/disk capacity and a CUDA-visible
unified memory pool, so no instance is currently classified
`NOT_RUN_RESOURCE_LIMIT`; it remains `NOT_RUN`.

## CPU, CUDA, And Legacy

Smoke CPU/CUDA decisions and iteration counts match on all 13 instances.
Regression has one terminal difference, `pilot4`, explained above. The six
common rewrite time limits are `fit1d`, `d2q06c`, `greenbeb`, `greenbea`,
`pilot`, and `sierra`.

Against legacy on the same 103-instance regression set:

- rewrite CPU passes while legacy fails on 71 instances;
- both pass on 26;
- both fail on 4;
- legacy passes while rewrite times out on `fit1d` and `greenbeb`.

The objective differences above `1e-3` are `modszk1`, `bnl2`, `lotfi`, and
the time-limited rewrite result for `pilot`. Every corresponding legacy row is
rejected by the independent gate. No legacy CUDA comparison is claimed because
there is no distinct legacy build with `CUPDLP_GPU=ON`.

## Disagreement Classes

The current `disagreements.csv` contains:

| Classification | Rows |
| --- | ---: |
| `TIME_OR_ITERATION_LIMIT` | 46 |
| `BASELINE_DEFECT` | 215 |
| `UNSUPPORTED_STATUS_OR_CERTIFICATE` | 4 |
| resolved `REWRITE_DEFECT` | 4 |
| resolved `REFERENCE_OR_MODEL_ISSUE` | 1 |
| `UNRESOLVED` | 0 |

This zero applies only to the executed smoke and Netlib regression scopes plus
the historical fixes. It is not a Phase 7 completion claim for the missing
full-Netlib and energy scopes.

## Checks Not Run

- Full Netlib: rewrite CUDA, legacy, in-tree cuPDLP, cuPDLPx, cuOpt, simplex,
  and IPM.
- Energy regression: all eight requested solvers.
- Full energy: all available solvers.
- A distinct legacy CUDA comparison, because no such build is configured.
- The two historical pre-implementation target-absence checks; no retained
  command or session artifact can prove them retroactively.

The missing full-Netlib and energy matrices were not started as partial
background work: completing the required sequential eight-solver, 300-second
matrix is a separate multi-hour benchmark campaign. No relevant benchmark or
debugger process remains active. This is a scheduling gap, not a resource
refusal or a solver pass; every missing row remains `NOT_RUN`.

The exact `unittest` command printed in Task 4 rejects an absolute file path as
a module name. The equivalent pinned direct-file invocation passed and is
recorded in the implementation plan.

A repository-wide `git diff --check` in `lp-benchmarks` reports pre-existing
trailing whitespace in `manifests/netlib-feasible.csv` and
`suites/netlib-feasible-reference.csv`. The scoped check for the touched
cuPDLPx adapter and test is clean; the unrelated reference rows were not
modified.

All recorded timings are descriptive single runs. They are not the formal
five-repetition performance gate and must not be used for the Phase 8 cutover
performance decision.

## Phase 8 Decision

**Do not continue to Phase 8.** The Phase 7 exit criteria lack evidence for the
full Netlib comparison and both energy scopes, and those absent scopes have not
been triaged. Do not delete legacy code. The pre-existing uncommitted default
rewrite routing should remain reversible until these gates and the separate
performance gate are complete.