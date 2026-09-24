# CUDA Check-Path Optimization Plan

**Goal:** Complete the approved remaining CUDA reduction, check reuse, and kernel fusion work without relaxing numerical acceptance.

**Architecture:** Keep algorithm ordering in the shared solver. Backends provide checked vector differences and batched norm/dot primitives. Resident results carry independent validity for raw stopping metrics, prepared fixed-point products, optional reductions, and certificate screens. CUDA Graph uses the same shared step.

**Scope:** User-approved follow-up to the CUDA audit on 2026-09-22. No legacy, scaling, PID formula, certified norm bound, or KKT v2 tolerance changes. Preserve existing user edits and archived experiments; no commit or branch.

## Stages

- [x] Freeze before sources/library in `20260922-cuda-check-optimizations`.
- [x] Add RED long difference-norm precision test; use checked parallel differences and cuBLAS with setup-preallocated scratch. Verify zero runtime allocations and sticky errors; 20 cases/2720 assertions pass.
- [x] Add `Backend::reduceFixedPointMetrics(delta_primal, delta_dual, aty_delta, primal_norm, dual_norm, cross_term)`. CUDA uses two copies/one wait; RED3 counter assertions, GREEN83 consumer cases/8245 assertions.
- [x] Route FPE through the batch and reuse valid prepared deltas/products. Fuse copy-plus-subtract in check preparation; prepared-path audit has no extra SpMV, no D2D copies, and one wait. Optional-error/terminal precedence tests pass.
- [x] Separate raw, fixed-point, anchor-distance, and prepared-product validity. Keep the 112-byte packet using explicit flag bits; independent supported-environment certificate screen. RED11 group assertions; GREEN143 packet/solver/certificate cases/39182 assertions.
- [x] Fuse primal/dual blend into one launch. Fuse forward SpMV validation only into its immediate dual projection consumer, checking the product before clamping; standalone and transpose SpMV keep their scans. Direct/Graph and Inf-clamping regression pass. New grid-stride loops use overflow-checked wide indices.
- [x] Use CUDA 12+ cuBLAS wide-index entry points; retain guarded older-toolkit paths. Keep already fused checked scale/axpy. Also batch compatibility convergence norms: RED3 transfer/wait assertions, GREEN8 cases/305 assertions.
- [x] Run focused tests after every stage. CPU-only regression: 170 cases/22893 assertions. CUDA regression: 247 cases/57722 assertions. After the reviewed grid-index fix: 18 focused cases/1584 assertions. Compute Sanitizer: 18 cases/1855 assertions, zero errors/leaks. No gigantic near-INT_MAX allocation or cuBLAS vendor-status injection was run.
- [x] Freeze production before/after with matching sources/drivers and 2118 pins. Both cases use 1 warmup + 5 measurements, CPU5, 120s/100000-step budgets; all 24 records accepted independently. Energy50 median 21.771453159 -> 9.372434466s (57.0% less), afiro 0.017405984 -> 0.016233597s (6.7% less). Energy50 solution bits differ but KKT/reference pass; steps unchanged.
- [x] Record results, retained guards, resource costs, ABI changes, untested limits, and protocol caveats in `docs/hipdlp-original-stopping.md`. No full legacy release-gate claim.

## Verification Commands

The new isolated CUDA build is `/home/yanyu/builds/hipdlp2-cuda-checks-cuda-20260922`.
Focused tests use `[check-optimization]`, `[parallel-reduction]`, `[fixed-point-error]`, `[phase7p]`, and `[allocation-audit]` as appropriate. Production measurements never use a test-hook library. Runs are serial on CPU5; no simultaneous GPU profiling/build workload and no changes to power or clocks.

An optimization is retained only after its correctness checks and scoped measurement or operation-count evidence pass. Unsafe finite-check removal is not an acceptable shortcut; explain and retain a guard if no equivalent immediate consumer exists.