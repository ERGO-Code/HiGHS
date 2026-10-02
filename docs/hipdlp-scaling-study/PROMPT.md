Continue the HiPDLP scaling study in this repo (branch `feature/hipdlp-rewrite`).
First read `docs/hipdlp-scaling-study/HANDOUT.md` fully; tools are in `docs/hipdlp-scaling-study/tools`
and prior results in `docs/hipdlp-scaling-study/data`. Benchmark LPs (neso-2245/2005 and the 232-instance
list in `data/instances.txt`) must be treated as read-only; copy them to a scratch directory.

Tasks, in order:
1. Build the rewrite (Release, `-DALL_TESTS=ON`, CUDA if available) and confirm the committed Curtis-Reid
   guard reproduces: neso-2245-filtered `--scaling auto` logs "Curtis-Reid rejected" and converges in
   ~106k iterations; run `unit_tests "[hipdlp]"` and compare failures against HEAD~1.
2. Prototype the two-pilot scaling selector described in the handout (200 iterations each with
   CR+Ruiz+PC and Ruiz+PC, geometric mean of the Kw estimates, continue the winner) behind an option,
   keeping the guard in front. Validate on neso (CUDA) and on instances NOT in `data/instances.txt`
   before trusting the ~20% gain; report solved counts, iteration geometric means and worst regressions.
3. Investigate why `stair` regresses 18x under the selector.
Use fair A/B: same tolerance for `kkt_tolerance` and `pdlp_optimality_tolerance`, iteration counts on the
CPU backend for policy comparisons. Do not change `verifyKkt` without asking me.
