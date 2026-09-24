# HiPDLP Large-Instance CPU/CUDA Screening

Date: 2026-09-22

Subsequent correctness work: [original-model KKT v2](hipdlp-original-stopping.md)
rechecks saved solutions under an explicit revised contract. The measurements,
acceptance definitions, and tables in this screening report remain historical
and unchanged; they are not silently reclassified with the new verifier.

## Outcome

**Completed all 80 planned jobs: 20 real continuous LPs, true legacy and
current rewrite, on both CPU and CUDA.** This is a bounded, single-cold-run
screening, not a repeated performance acceptance campaign. The current
rewrite includes CPU redundant-clear removal, dedicated CUDA
projection/reflection fusion, and production command-buffer replay.

Fifteen instances have at least one million matrix nonzeros, including three
with more than ten million. Five are explicitly labeled medium, not large.
No case was selected, replaced, or removed based on this campaign's outcomes.

| Role | Planned / completed | Reports OPTIMAL | Original KKT passes | KKT plus available reference | Full acceptance | Solver time limits | Process watchdogs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Legacy CPU | 20 / 20 | 8 | 1 | 1 | 0 | 12 | 0 |
| Rewrite CPU | 20 / 20 | 6 | 6 | 5 | 5 | 14 | 0 |
| Legacy CUDA | 20 / 20 | 14 | 2 | 2 | 0 | 6 | 0 |
| Rewrite CUDA | 20 / 20 | 6 | 6 | 5 | 5 | 13 | 1 |

Full acceptance additionally requires the reported reduced costs to match
those recomputed from the original objective and row duals within an absolute
1e-4 tolerance. This is an extra, scale-dependent condition, not the same as
the normalized KKT test. It rejects the three otherwise KKT-passing legacy
records. Thus **zero full legacy accepts does not mean every legacy result
failed KKT**, and does not establish that rewrite is faster.

There are zero fully accepted legacy/rewrite pairs on either backend. Both
full-cohort and joint-success diagnostic time geometric means are unavailable.
The earlier 23-case acceptance gate is unchanged; this different cohort is
not a replacement denominator.

## Protocol

- The user requested approximately 20 large instances on both CPU and GPU.
  With the user unavailable to choose a longer budget, the recommended
  first-pass option was used: zero warmups, one measured fresh process per
  role/instance, 60 seconds configured solver time, and a 180-second process
  watchdog including model input and output. No eight-hour repeated campaign
  was silently started.
- All four roles use FP64, one thread, presolve off, PID, Halpern restart,
  seed 20260922, and all four solver tolerances at 1e-4. Legacy keeps its
  native scaling; rewrite keeps auto scaling and the default-off early-KKT
  option. Identical numeric tolerances do not imply identical stopping rules.
- Solver processes run serially on CPU5 of NVIDIA GB10 / AArch64. Instance
  order is frozen with seed 20260922; adjacent role/backend order alternates.
  No clock, power, privilege, other-user process, or driver settings changed.
- True legacy is the frozen pre-rewrite source, not rewrite with legacy
  scaling. Each backend uses matched compiler families and Release settings:
  CPU Clang 18; CUDA GCC 13 / NVCC 13, sm_121; no LTO or testing hooks.
  Current CPU and CUDA production libraries were rebuilt in separate new
  directories, with matching driver headers and explicit runtime paths.
- The same wrapper timer includes preprocessing, scaling, allocation/upload,
  norm initialization, graph construction, iteration/check work, postsolve,
  in-solver verification where performed, and local teardown. Model reading,
  solution serialization, and the independent verifier are outside that timer
  and retained separately. Returned host solutions imply GPU completion.
- Each of the 79 returned solution records was checked by the separately
  frozen original-model verifier, without deriving replacement reduced costs.
  One watchdog job returned no complete solution record. Process failures,
  missing outputs, nonfinite values, KKT failures, and time overruns are not
  accepted even if the solver printed OPTIMAL.
- Trusted objective references exist for only `osa-60` and
  `pypsa-eur-elec_50_168h`; relative objective error must be at most 1e-3.
  The other 18 instances are explicitly KKT-only with respect to objective
  reference coverage. No reference was inferred from either compared solver.

No confidence interval or stable speedup is claimed from one cold solve.
The 60-second threshold cannot distinguish methods that both need longer.
Thresholds, cohort, and repetitions were not changed after seeing results.

## Instance Sizes

Sixteen files were downloaded from the public
[Mittelmann LP collection](https://plato.asu.edu/ftp/lptestset/), with four
existing Netlib/energy cases. The public source describes many of these as
LP relaxations already converted from MIPLIB models; no integer constraints
were silently removed here. The existing MPS inventory tool verified the
downloaded files contain no integer/semi-continuous or nonlinear sections.

Sizes below are original-MPS inventory dimensions, not rewrite's transformed
or presolved matrix. Medium means 100000-999999 nonzeros, large means
1000000-9999999, and huge means at least 10000000.

| Instance | Rows | Columns | Nonzeros | Scale |
| --- | ---: | ---: | ---: | --- |
| a2864 | 22117 | 200787 | 20078717 | huge |
| bharat | 1299953 | 590519 | 3565462 | large |
| chromaticindex1024-7 | 67583 | 73728 | 270324 | medium |
| datt256_lp | 11077 | 262144 | 1503732 | large |
| ex10 | 69608 | 17680 | 1162000 | large |
| fhnw-binschedule1 | 772872 | 1141653 | 8611326 | large |
| irish-electricity | 104259 | 61728 | 523257 | medium |
| neos-3025225 | 91572 | 69846 | 9357951 | large |
| neos-5052403-cygnet | 38268 | 32868 | 4898304 | large |
| neos-5251015 | 486531 | 136971 | 1955388 | large |
| physiciansched3-3 | 266227 | 79555 | 1062479 | large |
| rmine15 | 358395 | 42438 | 879732 | medium |
| s100 | 14733 | 364417 | 1777917 | large |
| s250r10 | 10962 | 273142 | 1318607 | large |
| supportcase10 | 165684 | 14770 | 555082 | medium |
| woodlands09 | 194599 | 382147 | 2646003 | large |
| osa-60 | 10280 | 232966 | 1397793 | large |
| pypsa-eur-elec_50_168h | 134766 | 64877 | 258209 | medium |
| pypsa-eur-elec_16_1h_2050 | 6754193 | 3171294 | 12684994 | huge |
| pypsa-de-elec-trex_vopt_50_1h | 18370653 | 8638003 | 37770662 | huge |

## All Results

Times are **seconds**, not medians. S = fully accepted; T = solver time
limit; K = reports OPTIMAL but original-model KKT fails; C = KKT passes but
the extra absolute reduced-cost consistency check fails; R = KKT passes but
the trusted reference objective check fails; W = process watchdog with no
completed solver-time record. K may also fail reduced-cost consistency.
Limited and rejected rows are not valid times-to-an-accepted-solution.

| Instance | Legacy CPU | Rewrite CPU | Legacy CUDA | Rewrite CUDA |
| --- | ---: | ---: | ---: | ---: |
| a2864 | 60.371 T | 61.146 T | 13.512 K | 60.999 T |
| bharat | 60.229 T | 60.601 T | 60.121 T | 61.064 T |
| chromaticindex1024-7 | 0.218 C | 0.497 S | 0.344 C | 1.052 S |
| datt256_lp | 1.155 K | 3.583 S | 0.773 K | 4.659 S |
| ex10 | 0.682 K | 1.830 S | 0.684 K | 2.356 S |
| fhnw-binschedule1 | 60.761 T | 61.180 T | 60.188 T | 61.390 T |
| irish-electricity | 60.019 T | 60.067 T | 25.863 K | 60.065 T |
| neos-3025225 | 60.363 T | 60.655 T | 23.962 K | 60.635 T |
| neos-5052403-cygnet | 60.164 T | 60.302 T | 6.617 K | 60.225 T |
| neos-5251015 | 6.452 K | 60.274 T | 3.329 K | 60.318 T |
| physiciansched3-3 | 60.069 T | 60.122 T | 60.020 T | 60.118 T |
| rmine15 | 23.781 K | 60.135 T | 3.973 K | 60.196 T |
| s100 | 60.111 T | 60.203 T | 60.027 T | 60.260 T |
| s250r10 | 60.026 T | 60.151 T | 11.363 K | 60.177 T |
| supportcase10 | 13.152 K | 60.068 T | 2.743 K | 60.027 T |
| woodlands09 | 1.819 K | 6.530 S | 1.151 K | 11.601 S |
| osa-60 | 6.876 K | 29.795 S | 1.436 K | 36.690 S |
| pypsa-eur-elec_50_168h | 60.007 T | 10.247 R | 12.037 C | 21.288 R |
| pypsa-eur-elec_16_1h_2050 | 61.603 T | 63.824 T | 60.447 T | 80.369 T |
| pypsa-de-elec-trex_vopt_50_1h | 64.530 T | 69.437 T | 61.741 T | W (180.974 wall) |

Full precision, iteration counts, each KKT component, reference error, process
wall time, and rejection reasons are in
[the result CSV](hipdlp-phase8-large20-screen.csv). All 40 backend pairs remain
present; ratio fields are empty because no pair passes full acceptance.

## Interpretation

1. **Legacy's shorter reported times frequently stop short of common
   accuracy.** Nineteen of its 22 OPTIMAL records fail original KKT. For
   example, `woodlands09` legacy CUDA has gap 1.6787e-3, and `osa-60` legacy
   CUDA has primal 2.3090e-4 and dual 3.4370e-3. These cannot be treated as
   valid speedups against 1e-4-accurate rewrite solutions.
2. **The extra reduced-cost gate matters and must remain visible.** Both
   legacy `chromaticindex1024-7` results pass KKT; their reduced-cost absolute
   errors are about 1.09585e-4, just above the separate 1e-4 bound. Legacy
   CUDA on the 50-node energy case passes KKT and reference objective but has
   reduced-cost absolute error 1.7689. Ignoring only this extra gate would
   yield legacy counts 1/20 CPU and 2/20 CUDA, versus rewrite 5/20 on each.
   The sole common KKT-plus-reference case would be `chromaticindex1024-7`,
   with rewrite/legacy time ratios 2.2831 CPU and 3.0561 CUDA. This is an
   explicitly conditional sensitivity observation, not a changed gate or
   a cohort-wide performance claim.
3. **KKT at 1e-4 does not guarantee the requested reference-objective
   agreement.** Both rewrite backends finish the 50-node energy model at
   9000 iterations with objective about 31756489973.75, versus the frozen
   IPM reference 31872370698.00. Relative reference error is 0.00363577,
   above 0.001, although KKT passes. This reproduces the previously documented
   tolerance interaction; no criterion was relaxed in this screen.
4. **More nonzeros alone do not guarantee a GPU whole-solve win.** The five
   fully accepted rewrite cases finish at the same CPU/CUDA iteration counts
   (200, 480, 480, 640, 9280). Their GPU/CPU cold whole-solve ratios range
   from 1.2314 to 2.1179. Setup and checking remain inside these times; this
   does not isolate SpMV or kernel throughput, and one run is not a stable
   performance estimate.
5. **Huge-model initialization needs separate attention.** Rewrite CUDA on
   the 12.68-million-nonzero energy model reports TIME_LIMIT after 80.369
   wrapper seconds with zero iterations. Work before the iteration loop and
   returning its result already exceeds the configured budget; the evidence
   does not identify a particular setup subphase. The 37.77-million-nonzero
   model reaches the 180-second process watchdog, leaving only a begin event
   and no completed solve record. Its exact phase and cause are unproven.
   Neither number should be attributed to Graph replay without profiling.

## Next Work

Prioritize a common original-model stopping/accuracy contract and the huge
GPU initialization/timeout path before broadening the timing campaign. Any
diagnostic change to legacy stopping must remain a separately labeled variant,
not overwrite the frozen native legacy results.

Then rerun this fixed cohort with a predeclared longer budget, such as 300
seconds, and repeated measurements for meaningful completed comparisons.
Do not silently relax KKT/reference criteria or replace hard cases. Repeating
all current timeouts five times would mainly refine the timeout count, not
establish solve-time speedups.

## Evidence

All source snapshots, drivers, runtime bindings, downloads, input dimensions,
plans, raw logs, solutions, verifier outputs, process statuses, and seals are
under:

```text
/home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-large20
```

The immutable execution plan is `plan-v2.json`, digest
`01ecbc102cd867a32a760db826c8709a73ee5fe920fd1c02735370115ed555c7`.
The first preflight harness/plan/smokes were retained separately; no formal
job used that weaker preflight protocol. The hardened harness has 39
assessment self-tests and four successful readiness smokes. It binds verifier
identity and process results, seals the plan before execution, rejects
optimized Python, and records malformed or missing results as rejections.

Finalization accounts for all 80 jobs. An independent Node audit re-derived
the acceptance criteria from raw driver/verifier files and performed 2375
hash comparisons, including 1400 source/build/driver/runtime pins and all
20 model identities. The watchdog is a completed rejected attempt, not a
fabricated solver record. Verification can be repeated without solving:

```sh
node /home/yanyu/benchmark-results/lp/hipdlp-rewrite/phase8/20260922-large20/audit.mjs
```

No solver code, default options, prior evidence, branches, or commits were
changed by this screening task. Full release coverage and the legacy
replacement performance gate remain unestablished.