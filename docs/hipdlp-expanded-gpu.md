# Expanded GPU comparison after preprocessing optimization

Date: 2026-10-04 (Europe/London). GPU: RTX 4070 Ti SUPER, 16 GB.
CPU: Intel Core i7-14700. CUDA 12.8.61; driver 580.173.02.
`OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1` for every solver.

The rewrite source is `f36962d34e14f9ae10ac2c440cc653dcf0a919d1`. HPR-LP uses
the official C implementation at `79dbaa3fb082a1b826785d13dcd30d9ba350e056`;
cuPDLPx is `083b90e482f9f0e8229eac965d189a1ba386ca20`.

## Measurement and acceptance

All runs use GPU device 0 with presolve disabled and run serially.
Rewrite uses auto scaling and PID. HPR-C requests its default CR
(with its upstream guard). cuPDLPx is measured both with its default
12-sweep geometric stage and with `--geo_mean_iter 0`.

Every exported primal, row-dual and reduced-cost vector has 17 digits.
A persistent, read-only host verifier loads each LP once and checks
`hipdlp-original-kkt/2` at `1e-4` between solver runs. Its work is outside
all benchmark timers. Rewrite also retains its production KKT screen
and host acceptance check inside the solver.

For external solvers the predeclared native tolerance grid is `1e-4`,
`1e-6`, `1e-8`. The first run on this grid passing common KKT is selected;
a time limit or missing export ends the search. This does not identify
the earliest possible KKT-passing iterate or the fastest intermediate
tolerance. It can require extra accuracy from an external solver.
Each attempt has a 300-second native budget and an outer safety timeout.
Native budgets retain each implementation's scope; cuPDLPx uses CPU
`clock()` internally. The reported comparison timers all use wall time.

A table entry is one selected attempt. Tolerance calibration and failed
attempts are retained in the evidence and are excluded from that entry.
These are single runs with a fixed order; small timing differences do
not establish a statistically reliable ranking.

Six added models come from `/srv/mps_da/*.mps.gz`. Canonicalization
removed 0 entries and added 0 bound-only columns on each of these six.
Exact HiGHS input audits pass for matrices, costs, bounds, names, sense
and offset. NESO uses the previously audited canonical inputs described
in [the original comparison](hipdlp-scaling-study/NESO_GPU_COMPARISON.md).
The failed default-geometric NESO runs are retained from that study;
the current NESO runs retest rewrite, HPR-C and cuPDLPx without the
geometric stage using the new common timers.

## End-to-end process wall time

Seconds from process launch to exit, including program startup, input
reading, GPU preparation, solve, full-precision export and cleanup.
`FAIL` means the selected run did not pass common KKT.

| Model | Rewrite | HPR-LP-C | cuPDLPx default | cuPDLPx no geometric |
|---|---:|---:|---:|---:|
| 25fv47 | 0.868 | 1.520 | 1.972 | 1.570 |
| pds-20 | 1.871 | 1.571 | 1.571 | 1.420 |
| pypsa-eur-240t5n | 1.871 | 2.172 | 2.775 | 25.104 |
| nug20 | 1.621 | 1.420 | 1.570 | 1.470 |
| mod2 | 5.083 | 61.680 | 185.031 | 21.139 |
| a2864 | 51.741 | 29.816 | 16.926 | 17.877 |
| neso-2245 | 29.214 | 18.376 | prior FAIL | 14.013 |
| neso-2005 | 123.762 | 89.276 | prior FAIL | 48.413 |

## Model sizes

| Model | Rows | Columns | Nonzeros |
|---|---:|---:|---:|
| 25fv47 | 821 | 1,571 | 10,400 |
| pds-20 | 33,874 | 105,728 | 230,200 |
| pypsa-eur-240t5n | 119,299 | 56,899 | 237,928 |
| nug20 | 15,240 | 72,600 | 304,800 |
| mod2 | 34,774 | 31,728 | 165,129 |
| a2864 | 22,117 | 200,787 | 20,078,717 |
| neso-2245 | 45,715 | 93,151 | 12,980,419 |
| neso-2005 | 60,083 | 108,047 | 36,876,017 |

## Solver-call wall time

Seconds around the call receiving the already parsed LP and returning
the solution. Includes each solver's model preparation, scaling, GPU
initialization, iterations and postsolve; excludes text parsing and
solution export. HPR-C and cuPDLPx have passive CLI timer/export patches.

| Model | Rewrite | HPR-LP-C | cuPDLPx default | cuPDLPx no geometric |
|---|---:|---:|---:|---:|
| 25fv47 | 0.518 | 1.224 | 1.671 | 1.254 |
| pds-20 | 0.887 | 0.792 | 0.887 | 0.841 |
| pypsa-eur-240t5n | 0.732 | 1.313 | 2.197 | 24.461 |
| nug20 | 0.677 | 0.780 | 0.976 | 0.940 |
| mod2 | 4.435 | 61.071 | 184.533 | 20.664 |
| a2864 | 21.714 | 7.514 | 9.063 | 9.032 |
| neso-2245 | 11.806 | 4.470 | prior FAIL | 6.029 |
| neso-2005 | 74.859 | 50.834 | prior FAIL | 31.484 |

## Input reading and rewrite preprocessing

All values are seconds. Preprocessing is already included in the
rewrite solver-call time above. Parsing is outside that call.

| Model | Rewrite read | HPR-C read | cuPDLPx default read | cuPDLPx no geometric read | Rewrite preprocessing |
|---|---:|---:|---:|---:|---:|
| 25fv47 | 0.017 | 0.009 | 0.006 | 0.007 | 0.009 |
| pds-20 | 0.536 | 0.333 | 0.200 | 0.158 | 0.298 |
| pypsa-eur-240t5n | 0.679 | 0.410 | 0.201 | 0.199 | 0.309 |
| nug20 | 0.554 | 0.291 | 0.201 | 0.162 | 0.227 |
| mod2 | 0.316 | 0.250 | 0.122 | 0.122 | 0.135 |
| a2864 | 29.564 | 21.752 | 7.385 | 8.369 | 13.023 |
| neso-2245 | 17.011 | 13.487 | — | 7.530 | 6.801 |
| neso-2005 | 48.395 | 37.996 | — | 16.483 | 19.364 |

## Selected native tolerances and iterations

| Model | Solver | Native tolerance | Iterations | Native status | Common KKT |
|---|---|---:|---:|---|---|
| 25fv47 | rewrite | 1e-4 | 9560 | OPTIMAL | PASS |
| 25fv47 | hpr | 1e-6 | 56400 | OPTIMAL | PASS |
| 25fv47 | cupdlpx | 1e-6 | 66600 | OPTIMAL | PASS |
| 25fv47 | cupdlpx-no-gm | 1e-6 | 45000 | OPTIMAL | PASS |
| pds-20 | rewrite | 1e-4 | 6360 | OPTIMAL | PASS |
| pds-20 | hpr | 1e-6 | 10350 | OPTIMAL | PASS |
| pds-20 | cupdlpx | 1e-6 | 10600 | OPTIMAL | PASS |
| pds-20 | cupdlpx-no-gm | 1e-6 | 10600 | OPTIMAL | PASS |
| pypsa-eur-240t5n | rewrite | 1e-4 | 2440 | OPTIMAL | PASS |
| pypsa-eur-240t5n | hpr | 1e-6 | 28200 | OPTIMAL | PASS |
| pypsa-eur-240t5n | cupdlpx | 1e-6 | 19400 | OPTIMAL | PASS |
| pypsa-eur-240t5n | cupdlpx-no-gm | 1e-6 | 501600 | OPTIMAL | PASS |
| nug20 | rewrite | 1e-4 | 3840 | OPTIMAL | PASS |
| nug20 | hpr | 1e-6 | 15000 | OPTIMAL | PASS |
| nug20 | cupdlpx | 1e-6 | 17600 | OPTIMAL | PASS |
| nug20 | cupdlpx-no-gm | 1e-6 | 17600 | OPTIMAL | PASS |
| mod2 | rewrite | 1e-4 | 140520 | OPTIMAL | PASS |
| mod2 | hpr | 1e-8 | 3785850 | OPTIMAL | PASS |
| mod2 | cupdlpx | 1e-8 | 5718800 | OPTIMAL | PASS |
| mod2 | cupdlpx-no-gm | 1e-6 | 626600 | OPTIMAL | PASS |
| a2864 | rewrite | 1e-4 | 8440 | OPTIMAL | PASS |
| a2864 | hpr | 1e-6 | 8700 | OPTIMAL | PASS |
| a2864 | cupdlpx | 1e-6 | 10200 | OPTIMAL | PASS |
| a2864 | cupdlpx-no-gm | 1e-6 | 10200 | OPTIMAL | PASS |
| neso-2245 | rewrite | 1e-4 | 6680 | OPTIMAL | PASS |
| neso-2245 | hpr | 1e-6 | 10500 | OPTIMAL | PASS |
| neso-2245 | cupdlpx-no-gm | 1e-6 | 9800 | OPTIMAL | PASS |
| neso-2005 | rewrite | 1e-4 | 34360 | OPTIMAL | PASS |
| neso-2005 | hpr | 1e-8 | 39750 | OPTIMAL | PASS |
| neso-2005 | cupdlpx-no-gm | 1e-6 | 20000 | OPTIMAL | PASS |

## Interpretation

The rewrite has the shortest solver-call time on `25fv47`,
`pypsa-eur-240t5n`, `nug20`, and `mod2` in this sample. Including input
reading changes the ordering on `nug20`. `pds-20` is close at the
solver-call level. `a2864` still exposes a substantial rewrite
preprocessing/solve cost and slower input reading.

The NESO preprocessing optimization improves rewrite time but does
not close its gap to HPR-C on these two models. See the controlled
[before/after measurements](hipdlp-preprocessing.md), where iteration
counts and exported solution vectors are unchanged.

The current cuPDLPx no-geometric NESO runs already pass common KKT
at native `1e-6`; the previous study used `1e-8` for its tight
ablation. In particular, the shorter 2005 time here must not be
attributed to a cuPDLPx code optimization. Both native tolerance
and the reported timing scope differ from that older table.

cuPDLPx geometric scaling helps `pypsa-eur-240t5n` and hurts `mod2`
in this experiment; disabling it globally is not supported by these
results. Its min/max geometric scaling is distinct from CR.

These rankings are for reaching the common KKT acceptance contract
under the stated tolerance grid. They do not isolate kernel speed
or measure all solvers at identical internal stopping conditions.

## Reproduction and evidence

The benchmark tools are under `docs/hipdlp-scaling-study/tools/`:
`gpu_comparison_rewrite.cpp`, `gpu_comparison_verify.cpp`,
`gpu_comparison_run.py`, and `gpu_comparison_report.py`.
The input normalization and audit reuse `neso_canonicalize.py` and
`neso_input_audit.cpp`. Numerical solver sources are unchanged by
this expanded comparison. External CLI patches, binary hashes, input
and solution hashes, exact argv, all calibration attempts, KKT metrics and compressed
logs are retained in [the evidence directory](hipdlp-scaling-study/data/2026-10-04-expanded-gpu/).

Build the two adapters against the same Release CUDA library. For example:

```sh
c++ -O3 -std=c++11 -Ihighs -Ibuild-scaling-cuda \
  docs/hipdlp-scaling-study/tools/gpu_comparison_rewrite.cpp \
  -Lbuild-scaling-cuda/lib -Wl,-rpath,"$PWD/build-scaling-cuda/lib" \
  -lhighs -o /tmp/hipdlp-expanded-gpu/rewrite
c++ -O3 -std=c++11 -Ihighs -Ibuild-scaling-cuda \
  docs/hipdlp-scaling-study/tools/gpu_comparison_verify.cpp \
  -Lbuild-scaling-cuda/lib -Wl,-rpath,"$PWD/build-scaling-cuda/lib" \
  -lhighs -o /tmp/hipdlp-expanded-gpu/verify-server
```

External build options are recorded in the original comparison.
Apply the archived passive patches to the pinned external sources.
Place audited canonical inputs in the `models/` subdirectory using
the filenames recorded in each command JSON. Keep all solver binaries
and the shared library fixed for the whole serial run.

```sh
python3 docs/hipdlp-scaling-study/tools/gpu_comparison_run.py \
  --root /tmp/hipdlp-expanded-gpu --external /tmp/neso-gpu-comparison \
  --names 25fv47,pds-20,pypsa-eur-240t5n,nug20,mod2,a2864,neso-2245,neso-2005 \
  --limit 300
```
