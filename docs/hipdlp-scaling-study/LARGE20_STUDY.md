# Twenty large LP scaling timings (2026-10-02)

On this sample and an RTX 4070 Ti SUPER, **PC-only is the fastest of the four tested policies under the 180-second per-instance budget**. Ruiz+PC is close. The selector does not improve on either fixed policy. No policy solved all twenty: eleven instances timed out for every policy, so these measurements do not establish the minimum time to converge on the entire twenty-instance workload.

All twenty meet at least one size criterion after `readModel`: rows+columns >= 100,000, or nonzeros >= 1,000,000. The sample combines size-stratified, family-deduplicated choices with neso-2245/2005, pds-100 and four independent larger heldouts. Fourteen names are in the original 232-instance list. Five smaller draft entries were removed before any of their timings were run.

## Results

| Policy | Verified optimal | Observed total, including timeouts (s) | Common eight solved total (s) |
|---|---:|---:|---:|
| PC-only | 9/20 | 2178.684 | 166.886 |
| Ruiz+PC | 9/20 | 2210.761 | 197.134 |
| two-pilot selector | 9/20 | 2232.881 | 214.896 |
| auto: guarded CR+Ruiz+PC | 8/20 | 2367.954 | 203.697 |

A/L/S/P all solve the same eight cases. L/S/P also solve watson_1; auto times out there. Their common nine verified solves cost PC-only **194.405 s**, Ruiz+PC **226.325 s**, and selector **248.689 s**. The eleven all-policy timeouts are shown below, rather than treated as successful fast solves.

PAR2 (each failed solve charged 360 s): PC-only 4154.405 s; Ruiz+PC 4186.325 s; selector 4208.689 s; auto 4523.697 s.

A second serial pass on the nine L/P common verified solves confirms the order:

- PC-only: 194.602 s, 9/9 verified.
- Ruiz+PC: 220.811 s, 9/9 verified.
- PC-only uses 11.9% less time in this repeat, compared with 14.1% less in the primary run.
- Most of that advantage is neso-2245: PC-only 88.250 s versus Ruiz+PC 114.008 s in the repeat. Excluding neso-2245, repeat totals are 106.352 s and 106.803 s, only 0.4% apart.

For the measured PC-only pipeline, set `hipdlp_scaling=legacy`, `pdlp_scaling_mode=4`, and `hipdlp_scaling_selector=false`. All configurations here use the rewrite implementation; “legacy” in that setting names the scaling pipeline.

## Selector findings

The selector continues CR+Ruiz+PC on 11 instances, Ruiz+PC on 7, and skips both pilots after the committed Curtis–Reid guard rejects CR on the two neso instances. Eight instances have exactly tied candidate scores: L1_sixm250, graph40-40, neos1, pds-100, savsched1, square41, steiner_t, support19.

On the eight cases solved by every policy, selector total time is **5.5% greater than auto**. Across the nine L/S/P solved cases it is **9.9% greater than Ruiz+PC**, and **27.9% greater than PC-only**. Its better overall count than auto comes from watson_1, where choosing Ruiz+PC succeeds; it does not give the shortest total time.

Worst selector runtime ratios versus auto among common verified solves:

- savsched1: 6.804 / 2.816 s = 2.42x.
- graph40-40: 4.353 / 2.266 s = 1.92x.
- pypsa-eur-240t5n: 6.163 / 4.858 s = 1.27x.

PC-only is slower than auto on pypsa-eur-240t5n: 8.103 / 4.858 s = 1.67x. This is why the recommendation is scoped to aggregate time on this sample. Auto loses a verified solve on watson_1: its 180-second timeout compares with 27.519 s for PC-only and 29.191 s for Ruiz+PC.

The independent larger heldouts L1_sixm250, istefanov_med, support19 and square41 all time out for all four policies (0/4 each). They therefore provide no confirmation of the estimated 20% selector convergence gain. The prototype remains opt-in; the default is unchanged.

## Method and integrity

- Intel i7-14700, RTX 4070 Ti SUPER 16 GB; CUDA 12.8, Release build with `ALL_TESTS=ON` and `HIGHS_HIPDLP_CUDA=ON`. Exact hardware and source/binary hashes are in `data/2026-10-02-large20/`.
- Both `kkt_tolerance` and `pdlp_optimality_tolerance` are 1e-4. Fixed step size, presolve off, 180 s and 2,147,483,647 total iterations per CUDA solve. No change to `verifyKkt`.
- GPU solves run sequentially. Each LP is read once, the backend is warmed on a separate LP, and every timed solve has fresh state. Configuration order rotates across models. Frozen library and driver prevent rebuilds affecting a running experiment.
- Observed solver times include model preparation, scaling, pilot work and the solver’s final KKT check. Shared LP input reading is excluded; adding the same recorded 174.785 s to each policy does not change the ranking. The additional independent external KKT check is separately timed.
- Shared `/srv/mps_da` files are copied into scratch. All twenty compressed source hashes were unchanged after GPU timing, and every final scratch payload was validated. p28 contains a valid gzip member followed by 1144 bytes of SCIP text; `prepare_large.py` decodes its complete MPS payload while ignoring that non-gzip suffix. An old truncated p28 cache failed reading before any solve; it was corrected and the other nineteen payloads matched byte-for-byte.
- watson_1 was remeasured in a complete four-policy batch after a selection-refinement interruption. That complete batch supplies its primary records; the initial partial records are retained separately.

## CPU iteration check

All 80 CPU probes completed with identical 2,000-total-iteration and 180-second budgets. A/L/S/P each solve **2/20**: savsched1 and graph40-40. Their CPU convergence counts are:

| Instance | auto | Ruiz+PC | selector | PC-only |
|---|---:|---:|---:|---:|
| savsched1 | 1400 | 1400 | 1600 | 1400 |
| graph40-40 | 240 | 240 | 440 | 240 |

On just these two common solved cases, CPU iteration geometric means are 579.655 for A/L/P and 839.047 for S (ratio 1.447). The selector includes its 200 discarded pilot iterations. The other eighteen instances are censored by the iteration or time budget; their convergence iteration counts and ratios are unknown. These probes check pilot behavior and budget accounting rather than establish a CPU winner on the whole large-instance set. CUDA iteration counts remain raw trace data and are not used for policy iteration comparisons. CPU probes ran after GPU timings, in four workers; their wall times are not compared with the serial GPU measurements.

## Per-instance primary CUDA times

T marks a timeout with failed KKT acceptance; all other entries are verified optimal. Times include scaling and pilots.

| Instance | auto | Ruiz+PC | selector | PC-only |
|---|---:|---:|---:|---:|
| neso-2245 | 121.764 | 117.098 | 116.170 | 87.856 |
| neso-2005 | 181.090 T | 180.925 T | 181.128 T | 181.122 T |
| steiner_t | 13.698 | 12.955 | 16.539 | 12.601 |
| savsched1 | 2.816 | 1.993 | 6.804 | 1.658 |
| 798kRows | 180.375 T | 180.378 T | 180.486 T | 180.405 T |
| watson_1 | 180.107 T | 29.191 | 33.793 | 27.519 |
| ts-palko | 180.075 T | 180.062 T | 180.065 T | 180.060 T |
| graph40-40 | 2.266 | 1.705 | 4.353 | 1.541 |
| psched3-3 | 180.073 T | 180.114 T | 180.127 T | 180.144 T |
| osa-60 | 180.068 T | 180.095 T | 180.083 T | 180.086 T |
| cont11 | 180.050 T | 180.029 T | 180.063 T | 180.051 T |
| fome13 | 25.283 | 26.169 | 27.036 | 23.601 |
| p28 | 180.190 T | 180.168 T | 180.135 T | 180.184 T |
| neos1 | 8.638 | 8.196 | 8.855 | 8.029 |
| pypsa-eur-240t5n | 4.858 | 5.332 | 6.163 | 8.103 |
| pds-100 | 24.373 | 23.686 | 28.976 | 23.497 |
| L1_sixm250 | 180.803 T | 181.010 T | 180.647 T | 180.706 T |
| istefanov_med | 180.085 T | 180.089 T | 180.084 T | 180.089 T |
| support19 | 180.742 T | 181.010 T | 180.817 T | 180.866 T |
| square41 | 180.598 T | 180.558 T | 180.558 T | 180.565 T |

Raw primary results: [CSV](data/2026-10-02-large20/cuda180/results.csv), [JSON](data/2026-10-02-large20/cuda180/results.json). Repeat, CPU probes, compressed complete logs, input hashes and commands are retained alongside them.

Reproduce preparation with `tools/prepare_large.py`, time with `tools/large_compare.cpp` and `tools/large_sweep.py`, and aggregate one experiment directory with `tools/large_report.py`. Inputs and output paths are explicit command arguments; logs record every driver invocation.
