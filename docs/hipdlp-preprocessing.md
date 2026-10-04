# HiPDLP preprocessing profiling and sparse construction

Measured on 2026-10-04 with an RTX 4070 Ti SUPER. The starting commit is
`8703b5cde`, pulled from `origin/feature/hipdlp-rewrite`. It already includes
overlapped host verification and the CUDA original-KKT screen. The baseline
adds optional preprocessing timers only; the optimized build adds the sparse
construction and working-storage changes described below.

## Results

Both runs use CUDA, PID, auto scaling, presolve off, and original KKT and solver
tolerances of `1e-4`. Source files are `/srv/mps_da/neso-{2245,2005}.mps.gz`.
They use the same audited canonical inputs as the
[external solver comparison](hipdlp-scaling-study/NESO_GPU_COMPARISON.md).
The table records wall time within the solve call, excluding MPS reading and
the additional independent post-run KKT check. GPU runs were serial.

| Instance | Preprocessing before / after, s | Reduction | Solve before / after, s | Iterations, both |
|---|---:|---:|---:|---:|
| neso-2245 | 16.096 / 6.968 | 56.7% | 21.210 / 11.975 | 6,680 |
| neso-2005 | 52.714 / 20.596 | 60.9% | 108.209 / 76.106 | 34,360 |

The original-coordinate primal, row-dual, and reduced-cost vectors are
**byte-for-byte identical** between builds on each instance. Both solutions
pass the independent original-model KKT verifier at `1e-4`. SHA-256 comparisons
are retained with the measurements. A separate preprocessing-only run also
matches the exact-field FNV-1a fingerprints of the matrices, transpose, bounds,
scales, norms and row mapping, both before and after scaling, on both models.
See [equivalence.json](hipdlp-scaling-study/data/2026-10-04-preprocessing/equivalence.json).
These are individual diagnostic GPU runs,
not statistical speedup estimates across hardware or model families.

The phase breakdown identifies where the time was saved:

| Phase | 2245 before / after, s | 2005 before / after, s |
|---|---:|---:|
| Model construction, total | 8.436 / 1.534 | 29.506 / 4.631 |
| Input matrix canonicalization and construction | 3.055 / 0.436 | 8.869 / 1.039 |
| Transformed matrix construction | 1.670 / 0.187 | 9.140 / 0.483 |
| Transpose construction | 3.057 / 0.293 | 9.282 / 0.901 |
| Unscaled norms | 0.607 / 0.611 | 2.084 / 2.199 |
| Scaling, total | 7.660 / 5.434 | 23.208 / 15.964 |
| Ruiz stages, including scale composition | 4.687 / 2.833 | 14.056 / 7.821 |
| Iteration loop | 4.307 / 4.236 | 53.751 / 53.752 |

Child timings are nested; they must not be added to their parents. Both auto
pipelines contribute to stage totals because CR is rejected on these inputs.
The working-copy timer covers the initial matrix/value copies, not every
allocation within a stage. Model construction includes validation, row metadata,
and cleanup not itemized in this table.

The GPU screen behaves identically before and after: 2245 makes 21 screens
(20 rejected), and 2005 makes 166 screens (165 rejected). Each solve invokes
the host KKT verifier once. The previous 6,880/34,520 iteration results used
the older `35b83fd9c` code and are not the baseline for this optimization.

## Implementation

`SolveDiagnostics::preprocessing` exposes optional model/scaling subphase wall
timings. Model and scaling APIs accept an optional diagnostics pointer; without
one the timers do not read the clock. The existing asynchronous acceptance path
preserves the subphase timings when restoring the accepted candidate's counters.

Sparse construction now proceeds as follows:

1. Build CSR directly from column-wise input by counting row entries, taking
   prefix sums, and scattering columns in increasing order. Row-wise input is
   copied directly into the owned CSR representation.
2. Skip sorting strictly ordered, unique rows. For a row with duplicates or
   unordered columns, use the existing column/value/sign ordering and
   long-double duplicate summation, then compact in place and remove exact
   zeros. This preserves deterministic handling of noncanonical inputs.
3. Construct transformed rows in final order and append their synthetic slack
   columns. No second canonicalization is necessary because row mapping is
   bijective and slack columns are distinct and follow the original columns.
   With no row permutation or slacks, move the original CSR storage instead.
4. Construct the transpose using counts and a row-ordered scatter. Every
   transpose row is already canonical, eliminating its global sort and triplet
   buffers. Guard the transformed column/nonzero counts against index overflow.

Ruiz reuses the immutable sparse pattern and copies only its working values.
Its row/column scratch vectors are reused between sweeps. Scale composition
updates private working values in place, avoiding allocation and zero-filling
of another full value array on each sweep. Failed working stages are discarded;
the public `Problem::applyDiagonalScaling` retains its transactional behavior.
The arithmetic, scale guards, stopping rules, and GPU screen are unchanged.

## Validation and reproduction

The Release CUDA build succeeds. `unit_tests '[hipdlp]'` reports:

| Build | Cases passed / total | Failed assertions |
|---|---:|---:|
| Instrumented baseline | 346 / 355 | 16 |
| Optimized | 348 / 357 | 16 |

The same nine cases fail before and after (four CPU/CUDA overflow oracle cases,
shared step ownership, two fixed-restart trace cases, and two legacy solver
characterizations). The two new cases add 83 passing assertions covering
row/column layouts, permutations, duplicates, cancelling entries, signed zeros,
transpose consistency, and invalid/nonfinite/overflowing entries. Existing
scaling and CUDA KKT screen tests also run in this suite.

The study tools are
`hipdlp-scaling-study/tools/neso_rewrite_profile.cpp` (full CUDA solve, timings,
KKT check and optional binary solution export) and
`hipdlp-scaling-study/tools/neso_preprocess_profile.cpp` (preprocessing-only
timings and exact-field fingerprints). The latter excludes fingerprint work
from the phase timers.

```sh
cmake --build build-scaling-cuda -j 8
build-scaling-cuda/bin/unit_tests '[hipdlp]'
c++ -O3 -std=c++11 -Ihighs -Ibuild-scaling-cuda \
  docs/hipdlp-scaling-study/tools/neso_rewrite_profile.cpp \
  -Lbuild-scaling-cuda/lib \
  -Wl,-rpath,/home/yzhou/Documents/Github/hipdlp2/build-scaling-cuda/lib \
  -lhighs -o /tmp/neso-profile
/tmp/neso-profile MODEL.mps auto /tmp/neso.solution
```

Raw measurements, argv, test logs, solution hashes, build hashes and source
patches are retained in
[the evidence directory](hipdlp-scaling-study/data/2026-10-04-preprocessing/).
Scratch binaries and large solution vectors remain in
`/tmp/neso-preprocessing-study/{baseline,optimized}`. Each scratch binary uses
its own copied shared library, so rebuilding the workspace does not change
the baseline being measured.
