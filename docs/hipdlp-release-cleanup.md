# HiPDLP release cleanup validation

Date: 2026-10-04. Starting source revision: `6d0b6cda6`.

This is a behavior-preserving cleanup, not a new numerical algorithm or a
release certification. The pre-existing test failures below remain open.

## Changes

- Moved the fixed-iteration replay implementation and header (about 1,100 lines)
  from `highs/pdlp/hipdlp/algorithm/` to `check/hipdlp/`. CMake and Meson now
  compile it into the trace tests, not the solver library. Source inspection
  tests follow the new paths.
- Made the migration adapter and validation executables opt-in CMake targets;
  the validation executables are also excluded from Meson's default build.
  These tools are not installed. The production entry point is `bin/highs`.
- Removed the migration adapter's ineffective `--presolve` flag. It directly
  calls the solver core and now consistently reports `presolve = off`.
  Presolved runs use the normal HiGHS driver or the four-solver study workflow.
- Kept context and diagnostic defaults beside their declarations, defaulted
  redundant move operations and destructors, simplified fixed-array loops,
  removed a redundant options alias, and formatted the solver control loop,
  option translation and backend factory with the repository's style.
- Documented the rewrite build, options, accuracy boundaries and source layout
  in the root README; added `hipdlp` to the main CLI's solver help.

The CR guard, scaling formulas, iteration scheduling, GPU KKT screening, host
verification, certificates, limits and optional experiments retain their
existing behavior. In particular, the default-off scaling pilot and early-KKT
paths have not been deleted as part of this cleanup. OpenMP/ROCm remain
explicitly unavailable reservations. There is no claimed solve-time speedup.

## Checks

| Check | Result |
| --- | --- |
| CUDA Release build: solver, unit tests, migration adapter | Passed |
| CPU-only Release build with `BUILD_TESTING=OFF`, then install | Passed; installed executable is `highs` only |
| CPU-only driver, `adlittle`, presolve on, backend `cpu` / `auto` | Both optimal; effective CPU |
| CPU-only driver, explicit backend `cuda` | Correct unavailable-backend error, nonzero exit |
| Migration adapter, removed `--presolve on` / `off` | Both rejected as unknown options |
| Migration adapter, GPU `afiro` without removed flag | Optimal; reports actual `presolve = off` |
| Shared-library symbol inspection, CPU and CUDA | No `hipdlp::algorithm::runFixedTrace` symbol |
| Installed headers | No `fixed_trace` header |
| Default CMake build dependencies | No migration, validation, verification or audit executables |
| HiPDLP unit suite before / after | Same 9 failing cases and 16 failing assertions; identical failure locations |
| Whitespace/error check | `git diff --check` passed |

The CUDA build uses the existing `build-scaling-cuda` configuration with tests
and `ALL_TESTS` enabled. The baseline executable and shared library were copied
before edits, and the baseline test suite ran against the unmodified sources.
The after-cleanup JUnit report records 77,779 assertions, including the 16
failures below. Meson source lists were updated but a Meson build was not run
because Meson is not installed in this environment.

## Actual-model regression

GPU: NVIDIA GeForce RTX 4070 Ti SUPER. Both versions use presolve on,
`pdlp_optimality_tolerance = kkt_tolerance = 1e-4`, fixed step size, Halpern
restart, 10 Ruiz sweeps and PC. `auto` additionally enables guarded CR.
The same HiGHS driver was run with the saved baseline library and the updated
library. All six pairs returned `Optimal`, the same iteration count, and
byte-identical complete solution files.

| Model | Scaling | Iterations before and after |
| --- | --- | ---: |
| afiro | Ruiz + PC (`legacy`) | 200 |
| adlittle | Ruiz + PC (`legacy`) | 7,160 |
| e226 | Ruiz + PC (`legacy`) | 29,600 |
| afiro | Guarded CR + Ruiz + PC (`auto`) | 280 |
| adlittle | Guarded CR + Ruiz + PC (`auto`) | 960 |
| e226 | Guarded CR + Ruiz + PC (`auto`) | 5,480 |

These runs check unchanged results, not speed: builds were running concurrently
and timings were not used as performance evidence. The unchanged numerical
implementation also means the original-model postsolve KKT limitations in the
[20-model study](hipdlp-four-solvers-presolve.md) still apply.

## Baseline failures still requiring investigation

No failing assertion was removed or relaxed to make this cleanup pass.

| Failing test | Failing assertions | Observed mismatch |
| --- | ---: | --- |
| `hi-pdlp-cpu-check-packet-oracle-fused-product-overflow` | 1 | Host oracle returns non-finite convergence error |
| `hi-pdlp-cpu-check-packet-oracle-dual-product-overflow` | 1 | Host convergence metric computation returns an error |
| `hi-pdlp-cuda-check-packet-oracle-fused-product-overflow` | 1 | Same oracle failure in CUDA comparison |
| `hi-pdlp-cuda-check-packet-oracle-dual-product-overflow` | 1 | Same oracle failure in CUDA comparison |
| `hi-pdlp-phase5h-shared-step-ownership-boundary` | 4 | Source inspection forbids symbols present in the current implementation |
| `hi-pdlp-fixed-trace-fixed-point-restart-cadence-is-deterministic` | 2 | Device-to-device transfer counts exceed the frozen expectation |
| `hi-pdlp-fixed-trace-fixed-point-restart-cadence-backend-agreement` | 2 | Same transfer-count mismatch |
| `hi-pdlp-legacy-optimal-kkt-characterization` | 2 | `e226` verification differs from expected pass/core-KKT agreement |
| `hi-pdlp-scrs8-legacy-characterization` | 2 | Fixed-step iteration count and objective differ from the frozen result |

Full release qualification must resolve these failures and repeat the relevant
checks. Their presence predates this cleanup; they are not evidence that the
current branch is fully release-ready.
