# CR scale-spread guard for rewrite auto scaling

Implemented on 2026-10-04, based on the policy in
[HPR-LP-C scaling.cu](https://github.com/PolyU-IOR/HPR-LP-C/blob/79dbaa3fb082a1b826785d13dcd30d9ba350e056/src/solver/scaling.cu#L237).
The numerical source before this change is the preprocessing optimization in
`f36962d34`; the complete starting checkout is `d3e4aca59`.

## Policy and computational work

`hipdlp_scaling=auto` now evaluates the CR candidate's **raw** row and column
log-scale ranges. Both ranges must be finite, nonnegative and at most `log(1e8)`.
The check happens before exponent clamping and before modifying the working
matrix. Empty dimensions have range zero; empty rows/columns retain zero log
scales and participate in the range, as in HPR-C.

If the candidate is rejected, the following enabled Ruiz/PC/L2 stages run on
the original working matrix. If it is accepted, they run after CR. Each stage
runs once. `curtis_reid` explicitly requests unguarded CR; `legacy` still means
the Ruiz/PC/L2 scaling pipeline, not the old solver implementation.

This replaces the former auto comparison against a second complete Ruiz+PC
pipeline. The old `scaledNormInflation` helper remains available for historical
diagnostic tools, but no longer determines auto selection. The historical
reference-pipeline timer remains zero. Scaling metadata records the two raw
log ranges and the rejection flag for later trigger experiments.

The guard is a **post-trial safety heuristic**, not a predictor of convergence
speed: the CR logarithms and configured sweeps must still be computed. A large
uniform scale can pass even when the former cost/RHS/bound-inflation guard would
reject it. The explicit regression test records this intentional policy change.
The solver's independent original-model KKT acceptance remains in effect.

The threshold policy matches HPR-C; it is evaluated on rewrite's transformed
matrix using rewrite's CR arithmetic. This does not claim identical CR factors
or decisions to HPR-C on every formulation.

## GPU scaling

This change removes duplicate CPU work. Scaling itself remains on the CPU.
Moving it to CUDA is a reasonable next implementation step: upload the sparse
pattern and values once, compute CR and its guard on the device, then run the
selected Ruiz/PC stages while the matrix remains resident. Avoid round-trip
matrix transfers on each sweep.

The current `Problem`, postsolve and original-KKT screen use consistent host
matrix/scaling metadata, so a CUDA implementation needs an explicit handoff
contract rather than merely calling kernels inside the current CPU loop.
Keep a CPU reference path and compare unscaled KKT results, extreme and empty
matrix behavior, and total preparation wall time including transfers. These
device changes are not part of this guard implementation.

## Validation

Release CUDA build succeeds. The scaling suite passes all 16 cases / 3,765
assertions. The full HiPDLP suite passes 350 of 359 cases; the same 16 failed
assertions in nine pre-existing cases match the archived baseline exactly.

New coverage checks the inclusive spread threshold, either dimension exceeding
it, nonfinite ranges, exponent clamping that would hide a bad spread, empty
dimensions, unchanged Ruiz+PC factors/matrix after rejection, forced CR,
rejection metadata when following stages are disabled, and the distinction
between spread and absolute data inflation.

Measurements below use the same eight audited inputs, GPU, options and public
KKT `1e-4` contract as [the earlier comparison](hipdlp-expanded-gpu.md).
The baseline is the earlier selected rewrite run; these are single diagnostic
runs and small timing differences are not statistically conclusive. No MPS
reading, extra post-run KKT check, fingerprint or binary export is included in
the reported solve-call time.

| Model | CR | Preprocessing before / after, s | Solve before / after, s | Iterations, both |
|---|---|---:|---:|---:|
| 25fv47 | kept | 0.009 / 0.007 | 0.518 / 0.511 | 9,560 |
| pds-20 | kept | 0.298 / 0.207 | 0.887 / 0.791 | 6,360 |
| pypsa-eur-240t5n | kept | 0.309 / 0.212 | 0.732 / 0.656 | 2,440 |
| nug20 | kept | 0.227 / 0.197 | 0.677 / 0.683 | 3,840 |
| mod2 | kept | 0.135 / 0.116 | 4.435 / 4.534 | 140,520 |
| a2864 | kept | 13.023 / 10.991 | 21.714 / 19.784 | 8,440 |
| neso-2245 | rejected | 6.801 / 6.315 | 11.806 / 11.448 | 6,680 |
| neso-2005 | rejected | 19.364 / 16.175 | 74.859 / 71.800 | 34,360 |

All eight runs pass common KKT and return byte-identical primal, row-dual and
reduced-cost vectors to the prior auto runs. Original/scaled preprocessing
fingerprints on both NESO inputs also match the earlier optimized study.
The six non-NESO models retain CR; both NESO models reject it. This is a
local validation set, not evidence that every former auto decision is unchanged.

Evidence, exact commands, input references, solution hashes, binary hashes,
tests and the source patch are in
[the study directory](hipdlp-scaling-study/data/2026-10-04-cr-spread-guard/).
The updated `neso_preprocess_profile.cpp` prints the raw CR ranges.
`neso_rewrite_profile.cpp MODEL.mps auto OUTPUT.bin 120` runs the timed
solve with independent post-run KKT and a 120-second native budget.
