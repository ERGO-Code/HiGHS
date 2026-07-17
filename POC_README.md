# HiPO-for-PyPSA proof of concept — build & run guide

This branch (`claude/hipo-parallel-solve-sketch`, based on HiGHS master v1.15.1)
carries an experimental set of changes aimed at PyPSA-class energy-system LPs:

1. **Runtime tuning options** for HiPO's factorisation (amalgamation, pivoting,
   regularisation, Metis seed — `hipo_sn_*`, `hipo_alpha_bk`, ... in
   `HighsOptions.h`), so tuning experiments need no recompile.
2. **Serial solve-phase fixes** (workspace reuse, swap-skip flags) that combine
   with upstream's new small-supernode fast path.
3. **An experimental parallel triangular solve** (`hipo_parallel_solve`):
   mode 1 = level schedule, mode 2 = block-level DAG. Default 0 (off).

Background and measurements: `HIPO_PYPSA_NOTES.md` and
`highs/ipm/hipo/PHASE3_DESIGN.md`.

## Build

Dependencies (Debian/Ubuntu): `sudo apt install cmake g++ libopenblas-dev`
(macOS: `brew install cmake openblas`). Metis is vendored.

```sh
cmake -S . -B build-hipo -DCMAKE_BUILD_TYPE=Release -DFAST_BUILD=ON -DHIPO=ON
cmake --build build-hipo --parallel
# solver binary: build-hipo/bin/highs
```

## Run on your own models

CLI, on an `.lp`/`.mps` exported by linopy/PyPSA:

```sh
build-hipo/bin/highs --solver hipo --options_file hipo.opts  mymodel.lp
```

with `hipo.opts`:

```
threads = 8              # your core count; helps the factorisation
parallel = on
user_bound_scale = -14   # if HiGHS warns about excessively large bounds
run_crossover = off      # if you do not need a basic solution (duals are
                         # still available from the interior point)
```

From PyPSA/linopy:

```python
n.optimize(solver_name="highs",
           solver_options={"solver": "hipo", "threads": 8, "parallel": "on",
                           "run_crossover": "off"})
# point linopy at this binary via PATH, or highspy built from this tree
```

The experimental parallel solve is opt-in: add `hipo_parallel_solve = 2`
(and optionally tune `hipo_dag_min_task_ops`, `hipo_dag_min_seg_rows`).

## Measured results (honest numbers)

Environment: 4-vCPU cloud VM (Intel Xeon, shared memory bandwidth), OpenBLAS.
Instances from the public benchmark bucket:

```sh
curl -O https://storage.googleapis.com/solver-benchmarks/instances/pypsa-de-elec-10-1h.mps.gz
curl -O https://storage.googleapis.com/solver-benchmarks/archive/pypsa-eur-elec-op-2-3h.lp.gz
gunzip *.gz
```

`pypsa-eur-elec-op-2-3h.lp` (310k rows, full solve to optimality, objective
8.8605241567e9 in every configuration):

| configuration | solve-phase time |
|---|---|
| HiGHS 1.14-dev, before this work | 10.2s |
| + serial fixes (this branch, pre-upstream-merge) | 9.5s |
| **v1.15.1 upstream fast path + serial fixes (current, default)** | **5.5s** |
| parallel solve mode 2, 4 threads | 20.5s |

`pypsa-de-elec-10-1h.mps` (4.5M rows raw — the hourly family where HiPO trails
Gurobi), iterations reached in a 600s budget:

| configuration | iterations | solve-phase time |
|---|---|---|
| HiGHS 1.14-dev, before this work | 41 | 239.6s |
| + serial fixes | 48 | 236.9s |
| **current default (serial)** | 45 | **159.7s** |
| parallel solve mode 2, 4 threads | 26 | 340.0s |

Takeaways:

- **The default serial path got ~33% faster on the solve phase** on the hourly
  family (upstream's small-supernode fast path and this branch's fixes compound).
- **The parallel solve (mode 2) is correct but LOSES on this VM**: triangular
  solves are memory-bandwidth-bound, and 4 shared vCPUs cannot feed 4 workers —
  parallel execution just bounces cache lines. It is bit-deterministic across
  thread counts and repeats and reaches the same optimum (verified), so it is
  safe to *evaluate*; whether it ever wins is an open question for real
  hardware with more memory channels (8+ core desktop/server). **That
  evaluation is exactly what this PoC is for.** If it loses on your machine
  too, the honest conclusion is that solve-phase parallelism needs the
  factorisation's data layout to change (blocked multi-RHS, or GPU offload) —
  valuable input for the upstream conversation either way.
- The biggest *convergence* lever found remains data scaling:
  `user_bound_scale=-14` reached a 2.5× smaller duality gap in equal time on
  hourly instances (their 1e10 "big-M" bounds are numerically hostile; the
  root fix belongs in PyPSA/linopy emitting true infinities).

## Status & caveats

- Experimental fork branch; not upstream-reviewed. HiPO itself is under very
  active development by its authors (this branch is merged up to v1.15.1).
- Validated: identical optimal objectives and clean KKT residuals on the
  instances above; parallel-solve determinism across runs and thread counts.
  Not validated: broad corpora (Netlib/Mittelmann), QP paths under mode 2.
- Next step: share findings with the HiPO authors (arXiv:2508.04370) before
  building further on the parallel-solve design.
